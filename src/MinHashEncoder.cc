#include "MinHashEncoder.h"

MinHashEncoder::MinHashEncoder() :
mpParameters(NULL), mpData(NULL), mInverseIndexPub(mInverseIndex) {

}

MinHashEncoder::MinHashEncoder(Parameters* apParameters, Data* apData, INDEXType apIndexType):mInverseIndexPub(mInverseIndex) {
	Init(apParameters, apData, apIndexType);
}

void MinHashEncoder::Init(Parameters* apParameters, Data* apData, INDEXType apIndexType) {
	mpParameters = apParameters;
	mpData = apData;
	numKeys=0;
	numFullBins=0;
	//init inverse index data structure
	mInverseIndex.clear();
	for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k){
		mInverseIndex.push_back(umap_uint_vec_uint());
	}
	mHashBitMask = numeric_limits<unsigned>::max() >> 1;
	mHashBitMask = (2 << (mpParameters->mHashBitSize - 1)) - 1;
	indexType = apIndexType;
	cout << "indexType " << indexType << endl;
}

inline vector<unsigned> MinHashEncoder::HashFuncNSPDK(const string& aString, unsigned aStart, unsigned aMaxRadius, unsigned aBitMask) {
	unsigned int hash = 0xAAAAAAAA;
	unsigned effective_end = min((unsigned) aString.size() - 1, aStart + aMaxRadius);
	unsigned radius = 0;
	vector<unsigned> code_list(aMaxRadius + 1, 0);
	for (std::size_t i = aStart; i <= effective_end; i++) {
		hash ^= ((i & 1) == 0) ? ((hash << 7) ^ aString[i] * (hash >> 3)) : (~(((hash << 11) + aString[i]) ^ (hash >> 5)));
		code_list[radius] = hash & aBitMask;
		radius++;
	}
	return code_list;
}

inline unsigned MinHashEncoder::HashFuncNSPDK(const vector<unsigned>& aList, unsigned aBitMask) {
	unsigned int hash = 0xAAAAAAAA;
	for (std::size_t i = 0; i < aList.size(); i++) {
		hash ^= ((i & 1) == 0) ? ((hash << 7) ^ aList[i] * (hash >> 3)) : (~(((hash << 11) + aList[i]) ^ (hash >> 5)));
	}
	return hash & aBitMask;
}

void MinHashEncoder::generate_feature_vector(const GraphClass& aG, SVector& x) {
	//assume there is a mMinRadius and a mMinDistance
	unsigned mRadius = mpParameters->mRadius;
	unsigned mDistance = mpParameters->mDistance;
	unsigned mMinRadius = mpParameters->mMinRadius;
	unsigned mMinDistance = mpParameters->mMinDistance;

	//assume 1 vertex with all info on the label
	string seq = aG.GetVertexLabel(0);
	unsigned size = seq.size();

	vector<vector<unsigned> > mFeatureCache;
	for (unsigned i = 0; i < size; ++i) {
		mFeatureCache.push_back(vector<unsigned>(mRadius, 0));
	}

	//create neighborhood features
	for (unsigned start = 0; start < size; ++start)
		mFeatureCache[start] = HashFuncNSPDK(seq, start, mRadius, mHashBitMask);

	vector<unsigned> endpoint_list(4);
	for (unsigned r = mMinRadius; r <= mRadius; r++) {
		endpoint_list[0] = r;
		for (unsigned d = mMinDistance; d <= mDistance; d++) {
			endpoint_list[1] = d;
			SVector z(pow(2, mpParameters->mHashBitSize));
			for (unsigned start = 0; start < size; ++start) {
				unsigned src_code = mFeatureCache[start][r];
				unsigned effective_dest = min(start + d, size - 1);
				unsigned dest_code = mFeatureCache[effective_dest][r];
				if (src_code > dest_code) {
					endpoint_list[2] = src_code;
					endpoint_list[3] = dest_code;
				} else {
					endpoint_list[3] = src_code;
					endpoint_list[2] = dest_code;
				}
				unsigned code = HashFunc(endpoint_list, mHashBitMask);
				//				z.coeffRef(code) += 1;
				//
				//				//add feature that ignore src endpoint
				//				//NOTE: this is important when the label of src vertex is considered noisy, in this way, it is only the context that will define the features
				//				endpoint_list[2] = 0;
				//				endpoint_list[3] = dest_code;
				//				unsigned nosrc_code = HashFunc(endpoint_list, mHashBitMask);
				//				z.coeffRef(nosrc_code) += 1;
				x.coeffRef(code) = 1;
			}
			z /= z.norm();
			x += z;
		}
	}
	x /= x.norm();
}

void MinHashEncoder::worker_readFiles(int numWorkers){

	while (!done){

		SeqDataSetP myData;
		unique_lock<mutex> lk(mut1);
		cv1.wait(lk,[&]{if ( (done) || (readFile_queue.try_pop(myData))) return true; else return false;});
		lk.unlock();
		if (!done && myData->filename != ""){
			cout << endl << "read next file " << myData->filename << " index " << myData->idx << endl;
			igzstream fin;
			fin.open(myData->filename.c_str());

			if (!fin)
				throw range_error("ERROR Data::LoadData: Cannot open file: " + myData->filename);

			bool valid_input = true;
			int file_instances = 0;
			string currSeq;
			while (!fin.eof() && valid_input ) {

				unsigned maxB = max(100,(int)log2((double)signature_counter)*50);
				unsigned currBuff = rand()%(maxB*2 - maxB + 1) + maxB;

				workQueueS myDataChunk;
				myDataChunk.gr.resize(currBuff);
				myDataChunk.offset= instance_counter; // offset goes over all files
				myDataChunk.dataSet = &(*myData);

				unsigned i = 0;
				while (i < currBuff && !fin.eof() && valid_input) {

					switch (myData->filetype) {
					case FASTA:
						mpData->SetGraphFromFASTAFile(fin, myDataChunk.gr.at(i),currSeq);
						break;
					case STRINGSEQ:
						mpData->SetGraphFromFile(fin, myDataChunk.gr.at(i));
						break;
					default:
						throw range_error("ERROR Data::LoadData: file type not recognized: " + myData->filetype);
					}

					if (myDataChunk.gr.at(i).IsEmpty()) {
						valid_input = false;
					} else {
						i++;
						instance_counter++;
						file_instances++;
					}
				}
				myDataChunk.gr.resize(i);

				//cout << " instances read " << instance_counter << " " << myDataChunk.gr.size() << " buffer " << currBuff << " full..." << graph_queue.size() << endl;
				if (graph_queue.size()>=numWorkers*5){
					unique_lock<mutex> lk(mut2);
					cv2.wait(lk,[&]{if ((done) || (graph_queue.size()<=numWorkers*2)) return true; else return false;});
					lk.unlock();
				}

				graph_queue.push(std::make_shared<workQueueS> (myDataChunk));
				cv2.notify_all();
			}
			fin.close();
			files_done++;
			myData->numSequences = file_instances;
			cout << " instances produced from file " << file_instances << endl;
		}
	}
}

void MinHashEncoder::worker_Graph2Signature(){

	while (!done){

		workQueueT myData;
		unique_lock<mutex> lk(mut2);
		cv2.wait(lk,[&]{if ( (done) || (graph_queue.try_pop( (myData) ))) return true; else return false;});
		lk.unlock();

		if (!done && myData->gr.size()>0) {
			cv2.notify_all();
			myData->sigs.resize(myData->gr.size());

			//cout << "  graph2sig thread got chunk " << myData->gr.size() << " offset " << myData->offset << " idx " << myData->dataSet->idx <<  " " << mpParameters->mHashBitSize << endl;

			for (unsigned j = 0; j < myData->gr.size(); j++) {
				SVector x(pow(2, mpParameters->mHashBitSize));
				generate_feature_vector(myData->gr[j], x);
				myData->sigs[j] = ComputeHashSignature(x);
			}

			sig_queue.push(myData);
			cv3.notify_all();
		}
	}
}

void MinHashEncoder::finisher(vector<vector<unsigned> >* myCache){
	ProgressBar progress_bar(1000);
	while (!done){

		workQueueT myData;
		unique_lock<mutex> lk(mut3);
		cv3.wait(lk,[&]{if ( (done) || (sig_queue.try_pop( (myData) ))) return true; else return false;});
		lk.unlock();

		if (!done && myData->sigs.size()>0) {
			if (myData->dataSet->updateIndex){
				switch (indexType) {
				case CLUSTER:
					for (unsigned j = 0; j < myData->sigs.size(); j++) {
						UpdateInverseIndex(myData->sigs[j], myData->offset+j);
					}
					break;
				case CLASSIFY:
					for (unsigned j = 0; j < myData->sigs.size(); j++) {
						UpdateInverseIndexHist(myData->sigs[j], myData->dataSet->idx-1);
					}
					break;
				default:
					throw range_error("Cannot update index! Index Type not recognized.");
				}
			}

			if (myData->dataSet->updateSigCache){
				if (myCache->size()<myData->offset+myData->sigs.size())
					myCache->resize(myData->offset+myData->sigs.size());
				for (unsigned j = 0; j < myData->sigs.size(); j++) {
					myCache->at(myData->offset+j) = myData->sigs[j];
				}
			}
			signature_counter += myData->sigs.size();

			cvm.notify_all();
			//cout << "    finisher updated index with " << mySigSet->sigs.size() << " signatures all_sigs=" <<  signature_counter << " inst=" << instance_counter << endl;
		}
		progress_bar.Count(signature_counter);
	}
}

void MinHashEncoder::LoadDataIntoIndexThreaded(vector<SeqDataSet>& myFiles, vector<vector<unsigned> >* myCache){

	if (mpParameters->mNumRepeatsHashFunction == 0 || mpParameters->mNumRepeatsHashFunction > mpParameters->mNumHashShingles * mpParameters->mNumHashFunctions){
		mpParameters->mNumRepeatsHashFunction = mpParameters->mNumHashShingles * mpParameters->mNumHashFunctions;
	}

	cout << "Using " << mpParameters->mNumHashFunctions << " hash functions (with factor " << mpParameters->mNumRepeatsHashFunction << " for single minhash)" << endl;
	cout << "Using " << mpParameters->mNumHashShingles << " as hash shingle factor" << endl;
	cout << "Using Feature radius   " << mpParameters->mMinRadius<<".."<<mpParameters->mRadius << endl;
	cout << "Using Feature distance " << mpParameters->mMinDistance<<".."<<mpParameters->mDistance << endl;
	cout << "Using FASTA win " << mpParameters->mSeqWindow<<" shift"<<mpParameters->mSeqShift << endl;

	done = false;
	files_done=0;
	signature_counter = 0;
	instance_counter = 0;

	int graphWorkers = std::thread::hardware_concurrency();
	if (mpParameters->mNumThreads>0)
		graphWorkers = mpParameters->mNumThreads;
	cout << "Using " << graphWorkers << " worker threads and 2 helper threads..." << endl;

	cout << "Computing MinHash signatures on the fly while reading " << myFiles.size() << " file(s)..." << endl;

	vector<vector<unsigned> >* myMinHashCache;
	{
		for (unsigned i=0;i<myFiles.size(); i++){
			if (myFiles[i].updateIndex)
				mIndexDataSets.push_back(myFiles[i]);
		}

		if (indexType == CLASSIFY)
				cout << "INDEX histogram length " << mIndexDataSets.size() << endl;

		for (unsigned i=0;i<myFiles.size(); i++){
			readFile_queue.push(std::make_shared<SeqDataSet>(myFiles[i]));
		}

		// use external MinHash cache
		if (myCache != NULL) {
			myMinHashCache = myCache;
		} else {
			// use standard cache
			myMinHashCache = &mMinHashCache;
		}
	}

	//create all threads
	threads.clear();
	threads.push_back( std::thread(&MinHashEncoder::finisher,this,myMinHashCache));
	for (int i=0;i<graphWorkers;i++){
		threads.push_back( std::thread(&MinHashEncoder::worker_Graph2Signature,this));
	}
	threads.push_back( std::thread(&MinHashEncoder::worker_readFiles,this,graphWorkers));

	{
		join_threads joiner(threads);

		unique_lock<mutex> lk(mutm);
		cv1.notify_all();
		while(!done){
			cvm.wait(lk,[&]{if ( (files_done<(int)myFiles.size()) || (signature_counter < instance_counter)) return false; else return true;});
			done=true;
			cv1.notify_all();
			cv2.notify_all();
			cv3.notify_all();
		};
	}

	cout << endl << "Inverse index ratio of overfull bins (maxSizeBin): " << ((double)numFullBins)/((double)numKeys) << " "<< numFullBins << "/" << numKeys << " instances " << instance_counter << endl;

	if (instance_counter > 0){
		mpData->SetDataSize(instance_counter);
	} else
		throw range_error("ERROR Data::LoadData: something went wrong: data file was expected but no data is available");
}

void MinHashEncoder::UpdateInverseIndex(vector<unsigned>& aSignature, unsigned aIndex) {
	for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k) { //for every hash value
		unsigned key = aSignature[k];
		if (key != MAXUNSIGNED && key != 0) { //if key is equal to markers for empty bins then skip insertion instance in data structure
			if (mInverseIndex[k].count(key) == 0) { //if this is the first time that an instance exhibits that specific value for that hash function, then store for the first time the reference to that instance
				vector<unsigned> tmp;
				tmp.push_back(aIndex);
				mInverseIndex[k].insert(make_pair(key, tmp));
				numKeys++; // just for bin statistics
			} else if (mInverseIndex[k][key].size() < mpParameters->mMaxSizeBin && mInverseIndex[k][key].front() != MAXUNSIGNED) {
				// add key to bin if we not have a full bin
				mInverseIndex[k][key].push_back(aIndex);
			} else if (mInverseIndex[k][key].size() == mpParameters->mMaxSizeBin){
				// if a bin is full we clear it and add a key with MAXUNSIGNED to indicate an overfull bin
				numFullBins++; // just for bin statistics
				mInverseIndex[k][key].clear();
				mInverseIndex[k][key].push_back(MAXUNSIGNED);
			}
		}
	}
}

void MinHashEncoder::UpdateInverseIndexHist(vector<unsigned>& aSignature, unsigned aIndex) {
	for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k) { //for every hash value
		unsigned key = aSignature[k];
		if (key != MAXUNSIGNED && key != 0) { //if key is equal to markers for empty bins then skip insertion instance in data structure
			if (mInverseIndex[k].count(key) == 0) { //if this is the first time that an instance exhibits that specific value for that hash function, then store for the first time the reference to that instance
				vector<unsigned> tmp(mIndexDataSets.size(),0);
				tmp[aIndex]++;
				mInverseIndex[k].insert(make_pair(key, tmp));
				numKeys++; // just for bin statistics
			} else {
				// add key to bin if we not have a full bin
				mInverseIndex[k][key][aIndex]++;
			}
		}
	}
}

void MinHashEncoder::CleanUpInverseIndex() {
	for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k) {
		for (umap_uint_vec_uint::const_iterator jt = mInverseIndex[k].begin(); jt != mInverseIndex[k].end(); ++jt) {
			unsigned hash_id = jt->first;
			if (hash_id != 0 && hash_id != MAXUNSIGNED) { //do not consider buckets corresponding to null bins
				unsigned collision_size = mInverseIndex[k][hash_id].size();

				if (collision_size < mpParameters->mMaxSizeBin) {
				} else {//remove bins that are too full from inverse index
					mInverseIndex[k].erase(hash_id);
				}
			}
		}
	}
}

vector<unsigned> MinHashEncoder::ComputeHashSignatureSize(vector<unsigned>& aSignature) {
	vector<unsigned> signature_size(mpParameters->mNumHashFunctions);
	assert(aSignature.size()==mpParameters->mNumHashFunctions);
	for (unsigned i = 0; i < aSignature.size(); ++i) {
		unsigned key = aSignature[i];
		signature_size[i] = mInverseIndex[i][key].size();
	}
	return signature_size;
}

vector<unsigned> MinHashEncoder::ComputeHashSignature(unsigned aID) {
	if (mMinHashCache.size()>0 && mMinHashCache[aID].size() > 0)
		return mMinHashCache[aID];
	else {
		throw range_error("ERROR internal MinHashCache is not filled!");
	}
}

vector<unsigned> MinHashEncoder::ComputeHashSignature(SVector& aX) {

	unsigned numHashFunctionsFull = mpParameters->mNumHashFunctions * mpParameters->mNumHashShingles;
	unsigned sub_hash_range = numHashFunctionsFull / mpParameters->mNumRepeatsHashFunction;
	vector<unsigned> signature(numHashFunctionsFull);

	//init with MAXUNSIGNED
	for (unsigned k = 0; k < numHashFunctionsFull; ++k)
		signature[k] = MAXUNSIGNED;

	//prepare a vector containing the signature as the k min values
	//for each element of the sparse vector
	for (SVector::InnerIterator it(aX); it; ++it) {
		unsigned feature_id = it.index();
		//for each sub_hash
		for (unsigned l = 1; l <= mpParameters->mNumRepeatsHashFunction; ++l) {
			unsigned key = IntHash(feature_id, MAXUNSIGNED, l);
			for (unsigned kk = 0; kk < sub_hash_range; ++kk) { //for all k values
				unsigned lower_bound = MAXUNSIGNED / sub_hash_range * kk;
				unsigned upper_bound = MAXUNSIGNED / sub_hash_range * (kk + 1);
				// upper bound can be different from MAXUNSIGNED due to rounding effects, correct this
				if (kk+1==sub_hash_range) upper_bound=MAXUNSIGNED;
				if (key >= lower_bound && key < upper_bound) { //if we are in the k-th slot
					unsigned signature_feature = kk + (l - 1) * sub_hash_range;
					if (key < signature[signature_feature]) //keep the min hash within the slot
						signature[signature_feature] = key;
				}
			}
		}
	}

	// compute shingles, i.e. rehash mNumHashShingles hash values into one hash value

	if (mpParameters->mNumHashShingles == 1 ) {
		return signature;
	} else {
		vector<unsigned> signatureFinal(mpParameters->mNumHashFunctions);
		for (unsigned i=0;i<mpParameters->mNumHashFunctions;i++){
			vector<unsigned> sigShinglet;
			for (unsigned j=i*mpParameters->mNumHashShingles;j<=i*mpParameters->mNumHashShingles+mpParameters->mNumHashShingles-1; j++){
				sigShinglet.push_back(signature[j]);
			}
			signatureFinal[i] = HashFunc(sigShinglet);
		}
		return signatureFinal;
	}
}

void  MinHashEncoder::ComputeApproximateNeighborhoodCore(const vector<unsigned>& aSignature, umap_uint_int& neighborhood, unsigned& collisions) {

	collisions = 0;
	for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k) {

		unsigned hash_id = aSignature[k];
		//cout << "hash func " << k << " id " << hash_id << endl;
		if (hash_id != 0 && hash_id != MAXUNSIGNED && mInverseIndex[k][hash_id].front() != MAXUNSIGNED) {

			//fill neighborhood set counting number of occurrences
			for (vector<unsigned>::iterator it = mInverseIndex[k][hash_id].begin(); it != mInverseIndex[k][hash_id].end(); ++it) {
				unsigned instance_id = *it;
				//cout << "inst "<< instance_id << endl;
				if (neighborhood.count(instance_id) > 0)
					neighborhood[instance_id]++;
				else
					neighborhood[instance_id] = 1;
			}
		} else {
			collisions++;
		}
	}
	//cout << " coll " << collisions << endl;
}

vector<unsigned> MinHashEncoder::ComputeApproximateNeighborhood(const vector<unsigned>& aSignature, unsigned& collisions, double& density) {

	umap_uint_int neighborhood;
	collisions = 0;
	density = 0;
	ComputeApproximateNeighborhoodCore(aSignature,neighborhood,collisions);
	//	cout <<"here "<< aSignature.size() << " " << collisions << " " << mpParameters->mMaxSizeBin << endl;
	return TrimNeighborhood(neighborhood, collisions, density);
}


umap_uint_int MinHashEncoder::ComputeApproximateNeighborhoodExt(const vector<unsigned>& aSignature, unsigned& collisions, double& density) {

	umap_uint_int neighborhood;
	collisions = 0;

	ComputeApproximateNeighborhoodCore(aSignature,neighborhood,collisions);
	vector<unsigned> myN = TrimNeighborhood(neighborhood, collisions, density);

	return neighborhood;
}

vector<unsigned> MinHashEncoder::TrimNeighborhood(umap_uint_int& aNeighborhood, unsigned collisions, double& density) {

	vector<unsigned> neighborhood_list;
	density = 0;
	double myC = (mpParameters->mPureApproximateSim * (mpParameters->mNumHashFunctions - collisions));
	//cout << "here 2 myC "<< myC << " nsize " << aNeighborhood.size() << " coll " << collisions << endl;
	for (umap_uint_int::const_iterator it = aNeighborhood.begin(); it != aNeighborhood.end();) {
		if ( (double)it->second >= myC ) {
			neighborhood_list.push_back(it->first);
			density += it->second;
			it++;
		} else {
			it=aNeighborhood.erase(it);
		}
	}
	if (collisions < mpParameters->mNumHashFunctions && neighborhood_list.size()!=0){
		density = ( density / neighborhood_list.size() ) / (mpParameters->mNumHashFunctions - collisions);
	} else
		density = 0;

	return neighborhood_list;
}

double MinHashEncoder::ComputeApproximateSim(const unsigned& aID, const unsigned& bID) {

	vector<unsigned> signatureA = ComputeHashSignature(aID);
	vector<unsigned> signatureB = ComputeHashSignature(bID);

	pair<unsigned,unsigned> simAB = ComputeApproximateSim(aID,signatureB);
	pair<unsigned,unsigned> simBA = ComputeApproximateSim(bID,signatureA);

	return (double)(simAB.first + simBA.first)/(2*mpParameters->mNumHashFunctions-simAB.second-simBA.second);
}


pair<unsigned,unsigned> MinHashEncoder::ComputeApproximateSim(const unsigned& aID, const vector<unsigned>& bSignature) {

	//umap_uint_int neighborhood;
	unsigned collisions = 0;
	unsigned counts_aID = 0;
	for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k) {
		unsigned hash_id = bSignature[k];
		if (hash_id != 0 && hash_id != MAXUNSIGNED && mInverseIndex[k][hash_id].front() != MAXUNSIGNED) {

			//fill neighborhood set counting number of occurrences
			for (vector<unsigned>::iterator it = mInverseIndex[k][hash_id].begin(); it != mInverseIndex[k][hash_id].end(); ++it) {
				if (*it == aID) counts_aID++;
			}

		} else
			collisions++;
	}
	pair<unsigned,unsigned> tmp=make_pair(counts_aID,collisions);
	return tmp;
}

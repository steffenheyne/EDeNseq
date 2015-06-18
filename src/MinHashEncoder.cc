#include "MinHashEncoder.h"

MinHashEncoder::~MinHashEncoder(){

}

MinHashEncoder::MinHashEncoder(Parameters* apParameters, Data* apData, INDEXType apIndexType) {
	Init(apParameters, apData, apIndexType);
}

void MinHashEncoder::Init(Parameters* apParameters, Data* apData, INDEXType apIndexType) {
	mpParameters = apParameters;
	mpData = apData;
	numKeys=0;
	numFullBins=0;
	mHashBitMask = numeric_limits<unsigned>::max() >> 1;
	mHashBitMask = (2 << (mpParameters->mHashBitSize - 1)) - 1;
	indexType = apIndexType;

	if (mpParameters->mNumRepeatsHashFunction == 0 || mpParameters->mNumRepeatsHashFunction > mpParameters->mNumHashShingles * mpParameters->mNumHashFunctions){
		mpParameters->mNumRepeatsHashFunction = mpParameters->mNumHashShingles * mpParameters->mNumHashFunctions;
	}

	cout << "MinHashEncoder object created of indexType " << indexType << endl;
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

	int file_instances = 0;
	while (!done){

		SeqDataSetP myData;
		unique_lock<mutex> lk(mut1);
		cv1.wait(lk,[&]{if ( ( (done) ||  ((mInstanceCounter <= mSignatureCounter) && (readFile_queue.try_pop(myData))))) return true; else return false;});
		lk.unlock();
		if (file_instances != 0) cout << " instances produced from file " << file_instances << endl;

		if (!done && myData->filename != ""){

			cout << endl << "read next file " << myData->filename << " index " << myData->idx << " " << " sig_all_counter " << mSignatureCounter << " inst_counter "<< mInstanceCounter << " thread " << std::this_thread::get_id() << endl;
			igzstream fin;
			fin.open(myData->filename.c_str(),std::ios::in);

			if (!fin)
				throw range_error("ERROR Data::LoadData: Cannot open file: " + myData->filename);

			bool valid_input = true;
			file_instances = 0;
			string currSeq;
			unsigned pos;

			while (  valid_input ) {

				unsigned maxB = max(100,(int)log2((double)mSignatureCounter)*1);
				unsigned currBuff = rand()%(maxB*2 - maxB + 1) + maxB;

				workQueueS myDataChunk;
				myDataChunk.gr.resize(currBuff);
				myDataChunk.names.resize(currBuff);
				myDataChunk.offset= mInstanceCounter; // offset goes over all files
				myDataChunk.dataSet = &(*myData);

				unsigned i = 0;
				while (i < currBuff && !fin.eof() && valid_input) {

					switch (myData->filetype) {
					case FASTA:
						mpData->SetGraphFromFASTAFile(fin, myDataChunk.gr.at(i),currSeq, pos, myDataChunk.names.at(i));
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
						mInstanceCounter++;
						file_instances++;
						if (myDataChunk.names.at(i) == "") {
							myDataChunk.names.at(i) = std::to_string(mInstanceCounter);
						}
						i++;
					}
				}
				myDataChunk.gr.resize(i);
				myDataChunk.names.resize(i);

				//cout << " instances read " << mInstanceCounter << " " << myDataChunk.gr.size() << " buffer " << currBuff << " full..." << graph_queue.size() << " " << std::this_thread::get_id() << endl;
				if (graph_queue.size()>=numWorkers*100){
					unique_lock<mutex> lk(mut2);
					cv2.wait(lk,[&]{if ((done) || (graph_queue.size()<=numWorkers*20)) return true; else return false;});
					lk.unlock();
				}

				graph_queue.push(std::make_shared<workQueueS> (myDataChunk));
				cv2.notify_all();
			}
			fin.close();
			files_done++;
			myData->numSequences = file_instances;
			cvm.notify_all();
		}
	}
}

void MinHashEncoder::worker_Graph2Signature(){

	while (!done){

		workQueueT myData;
		unique_lock<mutex> lk(mut2);
		cv2.wait(lk,[&]{if ( (done) ||  (graph_queue.try_pop( (myData) )) ) return true; else return false;});
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
			cv2.notify_all();
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
						UpdateInverseIndex(myData->sigs[j], myData->dataSet->idx);
					}
					break;
				default:
					throw range_error("Cannot update index! Index Type not recognized.");
				}
			}

			uint chunkSize = myData->sigs.size();
			if (myData->dataSet->updateSigCache){
				if (myCache->size()<myData->offset + chunkSize) {
					myCache->resize(myData->offset + chunkSize);
					idx2nameMap.resize(myData->offset + chunkSize);
				}
				for (unsigned j = 0; j < chunkSize; j++) {
					myCache->at(myData->offset+j) = myData->sigs[j];
					name2idxMap.insert(make_pair(myData->names[j],myData->offset+j));
					idx2nameMap.at(myData->offset+j) = myData->names[j];
				}
			}
			mSignatureCounter += chunkSize;
			//cout << "    finisher updated index with " << chunkSize << " signatures all_sigs=" <<  mSignatureCounter << " inst=" << mInstanceCounter << endl;
		}
		progress_bar.Count(mSignatureCounter);
		cv1.notify_all();
		cv2.notify_all();
		cvm.notify_all();
	}
}

void MinHashEncoder::LoadDataIntoIndexThreaded(vector<SeqDataSet>& myFiles, vector<vector<unsigned> >* myCache){


	uint numIndexUpdates = 0;
	for (unsigned i=0;i<myFiles.size(); i++){
		readFile_queue.push(std::make_shared<SeqDataSet>(myFiles[i]));
		if (myFiles[i].updateIndex)
			numIndexUpdates++;
	}

	// check which cache to use for signatures
	// use either external MinHash cache or class member
	vector<vector<unsigned> >* myMinHashCache;
	if (myCache != NULL) {
		myMinHashCache = myCache;
	} else {
		// use standard cache
		myMinHashCache = &mMinHashCache;
	}
	myMinHashCache->clear();

	cout << "Using " << mpParameters->mNumHashFunctions << " hash functions (with factor " << mpParameters->mNumRepeatsHashFunction << " for single minhash)" << endl;
	cout << "Using " << mpParameters->mNumHashShingles << " as hash shingle factor" << endl;
	cout << "Using feature radius   " << mpParameters->mMinRadius<<".."<<mpParameters->mRadius << endl;
	cout << "Using feature distance " << mpParameters->mMinDistance<<".."<<mpParameters->mDistance << endl;
	cout << "Using sequence window  " << mpParameters->mSeqWindow<<" shift"<<mpParameters->mSeqShift << endl;
	cout << "Computing MinHash signatures on the fly while reading " << myFiles.size() << " file(s)..." << endl;
	cout << "Update inverse index with " << numIndexUpdates << " file(s)..." << endl;

	// threaded producer-consumer model for signature creation and index update
	// created threads:
	// 	1 finisher that updates the index and signature cache,
	// 	1 to read files and produces sequence instances
	//		n worker threads that create the signatures

	int graphWorkers = std::thread::hardware_concurrency();
	if (mpParameters->mNumThreads>0)
		graphWorkers = mpParameters->mNumThreads;

	cout << "Using " << graphWorkers << " worker threads and 2 helper threads..." << endl;

	done = false;
	files_done=0;
	mSignatureCounter = 0;
	mInstanceCounter = 0;

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
			cvm.wait(lk,[&]{if ( (files_done<(int)myFiles.size()) || (mSignatureCounter < mInstanceCounter)) return false; else return true;});
			lk.unlock();
			done=true;
			cv3.notify_all();
			cv2.notify_all();
			cv1.notify_all();
		}

	} // leave block only if threads are finished and joined, destroys joiner

	// threads finished

	if (numIndexUpdates>0)
		cout << endl << "Inverse index ratio of overfull bins (maxSizeBin): " << ((double)numFullBins)/((double)numKeys) << " "<< numFullBins << "/" << numKeys << " instances " << mInstanceCounter << endl;

	if (mInstanceCounter == 0) {
		throw range_error("ERROR in MinHashEncoder::LoadData: something went wrong; no instances/signatures produced");
	} else
		cout << "Instances/signatures produced " << mInstanceCounter << endl;

	mpData->SetDataSize(mInstanceCounter);
}

vector<unsigned>& MinHashEncoder::ComputeHashSignature(unsigned aID) {
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
				//if (kk+1==sub_hash_range) upper_bound=MAXUNSIGNED;
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


void NeighborhoodIndex::UpdateInverseIndex(vector<unsigned>& aSignature, unsigned aIndex) {
	for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k) { //for every hash value
		unsigned key = aSignature[k];
		if (key != MAXUNSIGNED && key != 0) { //if key is equal to markers for empty bins then skip insertion instance in data structure
			if (mInverseIndex[k].count(key) == 0) { //if this is the first time that an instance exhibits that specific value for that hash function, then store for the first time the reference to that instance
				indexBinTy tmp;
				tmp.push_back(aIndex);
				mInverseIndex[k].insert(make_pair(key, tmp));
				numKeys++; // just for bin statistics
			} else if (mInverseIndex[k][key].size() < mpParameters->mMaxSizeBin && mInverseIndex[k][key].front() != MAXBINKEY) {
				// add key to bin if we not have a full bin
				mInverseIndex[k][key].push_back(aIndex);
			} else if (mInverseIndex[k][key].size() == mpParameters->mMaxSizeBin){
				// if a bin is full we clear it and add a key with MAXUNSIGNED to indicate an overfull bin
				numFullBins++; // just for bin statistics
				mInverseIndex[k][key].clear();
				mInverseIndex[k][key].push_back(MAXBINKEY);
			}
		}
	}
}


void NeighborhoodIndex::ComputeApproximateNeighborhoodCore(const vector<unsigned>& aSignature, umap_uint_int& neighborhood, unsigned& collisions) {

	collisions = 0;
	for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k) {

		unsigned hash_id = aSignature[k];
		if (hash_id != 0 && hash_id != MAXUNSIGNED && mInverseIndex[k][hash_id].front() != MAXBINKEY) {

			//fill neighborhood set counting number of occurrences
			for (indexBinTy::iterator it = mInverseIndex[k][hash_id].begin(); it != mInverseIndex[k][hash_id].end(); ++it) {
				binKeyTy instance_id = *it;
				if (neighborhood.count(instance_id) > 0)
					neighborhood[instance_id]++;
				else
					neighborhood[instance_id] = 1;
			}
		} else {
			collisions++;
		}
	}
}

vector<unsigned> NeighborhoodIndex::ComputeApproximateNeighborhood(const vector<unsigned>& aSignature, unsigned& collisions, double& density) {

	umap_uint_int neighborhood;
	collisions = 0;
	density = 0;
	ComputeApproximateNeighborhoodCore(aSignature,neighborhood,collisions);
	//	cout <<"here "<< aSignature.size() << " " << collisions << " " << mpParameters->mMaxSizeBin << endl;
	return TrimNeighborhood(neighborhood, collisions, density);
}


umap_uint_int NeighborhoodIndex::ComputeApproximateNeighborhoodExt(const vector<unsigned>& aSignature, unsigned& collisions, double& density) {

	umap_uint_int neighborhood;
	collisions = 0;

	ComputeApproximateNeighborhoodCore(aSignature,neighborhood,collisions);
	vector<unsigned> myN = TrimNeighborhood(neighborhood, collisions, density);

	return neighborhood;
}

vector<unsigned> NeighborhoodIndex::TrimNeighborhood(umap_uint_int& aNeighborhood, unsigned collisions, double& density) {

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

double NeighborhoodIndex::ComputeApproximateSim(const unsigned& aID, const unsigned& bID) {

	vector<unsigned> signatureA = ComputeHashSignature(aID);
	vector<unsigned> signatureB = ComputeHashSignature(bID);

	pair<unsigned,unsigned> simAB = ComputeApproximateSim(aID,signatureB);
	pair<unsigned,unsigned> simBA = ComputeApproximateSim(bID,signatureA);

	return (double)(simAB.first + simBA.first)/(2*mpParameters->mNumHashFunctions-simAB.second-simBA.second);
}


pair<unsigned,unsigned> NeighborhoodIndex::ComputeApproximateSim(const unsigned& aID, const vector<unsigned>& bSignature) {

	//umap_uint_int neighborhood;
	unsigned collisions = 0;
	unsigned counts_aID = 0;
	for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k) {
		unsigned hash_id = bSignature[k];
		if (hash_id != 0 && hash_id != MAXUNSIGNED && mInverseIndex[k][hash_id].front() != MAXBINKEY) {

			//fill neighborhood set counting number of occurrences
			for (indexBinTy::iterator it = mInverseIndex[k][hash_id].begin(); it != mInverseIndex[k][hash_id].end(); ++it) {
				if (*it == aID) counts_aID++;
			}

		} else
			collisions++;
	}
	pair<unsigned,unsigned> tmp=make_pair(counts_aID,collisions);
	return tmp;
}

void HistogramIndex::SetHistogramSize(unsigned size){
	mHistogramSize = size;
}


unsigned HistogramIndex::GetHistogramSize(){
	return mHistogramSize;
}


bool compare(const SeqDataSet& first, const SeqDataSet& second) {

	return (first.uIdx<second.uIdx);
}

void  HistogramIndex::PrepareIndexDataSets(vector<SeqDataSet>& myFileList){

	if (!myFileList.size())
		throw range_error("ERROR no datasets to prepare ...");

	sort(myFileList.begin(),myFileList.end(),compare);

	// check if we have datasets to update the index
	uint bin = 0;
	uint userIdx = myFileList[0].uIdx;
	for (unsigned i=0;i<myFileList.size(); i++){
		if (myFileList[i].updateIndex){
			if ( myFileList[i].uIdx != userIdx )
				bin++;
			mHistBin2DatasetIdx.insert(make_pair(bin,i));
			userIdx = myFileList[i].uIdx;
			myFileList[i].idx = bin;
			mIndexDataSets.push_back(myFileList[i]);
			//cout << "bin" << bin << " " << userIdx << " "<<myFileList[i].idx << " " << myFileList[i].uIdx << endl;
		}
	}
	SetHistogramSize(bin+1);
}

//void HistogramIndex::UpdateInverseIndex(vector<unsigned>& aSignature, unsigned aIndex) {
//	for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k) { //for every hash value
//		unsigned key = aSignature[k];
//		if (key != MAXUNSIGNED && key != 0) { //if key is equal to markers for empty bins then skip insertion instance in data structure
//			if (mInverseIndex[k].count(key) == 0) { //if this is the first time that an instance exhibits that specific value for that hash function, then store for the first time the reference to that instance
//			/*	indexBinTy t(1,aIndex);
//				t.reserve(1);
//				mInverseIndex[k][key]= t;*/
//				std::unique_ptr<bin> myP(new bin);
//				myP->val=aIndex;
//				mInverseIndex[k][key] = std::move(myP);
//				//mInverseIndex[k][key].reserve(1);
//				numKeys++; // just for bin statistics
//			} else if (mInverseIndex[k][key].back() != aIndex){
//				mInverseIndex[k][key].push_back(aIndex);
//			}
//		}
//	}
//}


void HistogramIndex::UpdateInverseIndex(vector<unsigned>& aSignature, unsigned aIndex) {
	for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k) { //for every hash value
		unsigned key = aSignature[k];
		if (key != MAXUNSIGNED && key != 0) { //if key is equal to markers for empty bins then skip insertion instance in data structure
			if (mInverseIndex[k].count(key) == 0) { //if this is the first time that an instance exhibits that specific value for that hash function, then store for the first time the reference to that instance

				binKeyTy * foo;
				foo = new binKeyTy[2];
				foo[0]= (binKeyTy)aIndex;
				foo[1]= MAXBINKEY;

				mInverseIndex[k][key] = foo;
				numKeys++; // just for bin statistics
			} else if (mInverseIndex[k][key][0] != (binKeyTy)aIndex){
				binKeyTy* foo = mInverseIndex[k][key];
				int i = 0;
				while (foo[i]!=MAXBINKEY){
					i++;
				}

				binKeyTy * fooNew;
				fooNew = new binKeyTy[i+2];
				fooNew[0]=(binKeyTy)aIndex;

				for (int j=0;j<i;j++){
					fooNew[j+1]=foo[j];
				}

				fooNew[i+1] = MAXBINKEY;
				mInverseIndex[k][key] = fooNew;
				delete[] foo;
			}
		}
	}
}


//void HistogramIndex::ComputeHistogram(const vector<unsigned>& aSignature, std::valarray<double>& hist, unsigned& emptyBins) {
//
//	hist.resize(GetHistogramSize());
//	hist *= 0;
//	emptyBins = 0;
//
//	for (unsigned k = 0; k < aSignature.size(); ++k) {
//		if (mInverseIndex[k].count(aSignature[k]) > 0) {
//
//			std::valarray<double> t(0.0, hist.size());
//
//			for (typename indexBinTy::const_iterator it = mInverseIndex[k][aSignature[k]].begin(); it != mInverseIndex[k][aSignature[k]].end();it++){
//				t[*it]=1;
//			}
//
//			hist += t;
//
//		} else {
//			emptyBins++;
//		}
//	}
//}

void HistogramIndex::ComputeHistogram(const vector<unsigned>& aSignature, std::valarray<double>& hist, unsigned& emptyBins) {

	hist.resize(GetHistogramSize());
	hist *= 0;
	emptyBins = 0;

	for (unsigned k = 0; k < aSignature.size(); ++k) {
		if (mInverseIndex[k].count(aSignature[k]) > 0) {

			std::valarray<double> t(0.0, hist.size());
			int i = 0;
			while (mInverseIndex[k][aSignature[k]][i] != MAXBINKEY){
				t[mInverseIndex[k][aSignature[k]][i]]=1;
				i++;
			}

			hist += t;

		} else {
			emptyBins++;
		}
	}
}

void HistogramIndex::writeBinaryIndex2(ostream &out, const indexTy& index) {
	// create binary reverse index representation
	// format:
	out.write((const char*) &mpParameters->mHashBitSize, sizeof(unsigned));
	out.write((const char*) &mpParameters->mRandomSeed, sizeof(unsigned));
	out.write((const char*) &mpParameters->mRadius, sizeof(unsigned));
	out.write((const char*) &mpParameters->mMinRadius, sizeof(unsigned));
	out.write((const char*) &mpParameters->mDistance, sizeof(unsigned));
	out.write((const char*) &mpParameters->mMinDistance, sizeof(unsigned));
	out.write((const char*) &mpParameters->mNumHashShingles, sizeof(unsigned));
	out.write((const char*) &mpParameters->mNumRepeatsHashFunction, sizeof(unsigned));
	out.write((const char*) &mpParameters->mSeqWindow, sizeof(unsigned));
	out.write(reinterpret_cast<char *>(&mpParameters->mSeqShift), sizeof(mpParameters->mSeqShift));

	unsigned numHashFunc = index.size();
	out.write((const char*) &numHashFunc, sizeof(unsigned));
	for (typename indexTy::const_iterator it = index.begin(); it!= index.end(); it++){
		unsigned numBins = it->size();
		out.write((const char*) &numBins, sizeof(unsigned));
		for (typename indexSingleTy::const_iterator itBin = it->begin(); itBin!=it->end(); itBin++){
			unsigned binId = itBin->first;

			unsigned i = 0;
			while (itBin->second[i]!=MAXBINKEY){
				i++;
			}

			unsigned numBinEntries = i;
			out.write((const char*) &binId, sizeof(unsigned));
			out.write((const char*) &numBinEntries, sizeof(unsigned));
			for (i=0;i<numBinEntries;i++){
				binKeyTy s = itBin->second[i];
				out.write((const char*) &(s), sizeof(binKeyTy));
			}

		}
	}
}

bool HistogramIndex::readBinaryIndex2(string filename, indexTy &index){
	igzstream fin;
	fin.open(filename.c_str());

	fin.read((char*) &mpParameters->mHashBitSize, sizeof(unsigned));
	fin.read((char*) &mpParameters->mRandomSeed, sizeof(unsigned));
	fin.read((char*) &mpParameters->mRadius, sizeof(unsigned));
	fin.read((char*) &mpParameters->mMinRadius, sizeof(unsigned));
	fin.read((char*) &mpParameters->mDistance, sizeof(unsigned));
	fin.read((char*) &mpParameters->mMinDistance, sizeof(unsigned));
	fin.read((char*) &mpParameters->mNumHashShingles, sizeof(unsigned));
	fin.read((char*) &mpParameters->mNumRepeatsHashFunction, sizeof(unsigned));
	fin.read((char*) &mpParameters->mSeqWindow, sizeof(unsigned));
	fin.read( reinterpret_cast<char*>( &mpParameters->mSeqShift ), sizeof mpParameters->mSeqShift);

	cout << " shift "<< mpParameters->mSeqShift << endl;

	unsigned numHashFunc = 0;
	fin.read((char*) &numHashFunc, sizeof(unsigned));
	if (numHashFunc <= 0)
		fin.setstate(std::ios::badbit);
	if (!fin.good())
		return false;
	index.resize(numHashFunc);
	for (unsigned  hashFunc = 0; hashFunc < numHashFunc; hashFunc++){

		unsigned numBins = 0;
		fin.read((char*) &numBins, sizeof(unsigned));
		if (numBins < 0)
			fin.setstate(std::ios::badbit);
		if (!fin.good())
			return false;
		index[hashFunc].resize(numBins);
		for (unsigned  bin = 0; bin < numBins; bin++){

			unsigned binId = 0;
			unsigned numBinEntries = 0;
			fin.read((char*) &binId, sizeof(unsigned));
			fin.read((char*) &numBinEntries, sizeof(unsigned));
			if (numBinEntries < 0 || binId <= 0)
				fin.setstate(std::ios::badbit);
			if (!fin.good())
				return false;
			indexBinTy tmp = new binKeyTy[numBinEntries+1];

			for (unsigned entry = 0; entry < numBinEntries; entry++ ){

				binKeyTy s;
				fin.read((char*) &s, sizeof(binKeyTy));
				if (s < 0)
					fin.setstate(std::ios::badbit);
				if (!fin.good())
					return false;
				tmp[entry] = s;
			}
			tmp[numBinEntries] = MAXBINKEY;
			index[hashFunc][binId] = tmp;
		}
	}
	fin.close();
	return true;
}

//double indicator(double i){if (i>0) return 1.0; else return 0.0;}

//void  HistogramIndex::ComputeHistogram(const vector<unsigned>& aSignature, std::valarray<double>& hist, unsigned& emptyBins) {
//
//	hist.resize(mIndexDataSets.size());
//	hist *= 0;
//	emptyBins = 0;
//
//	for (unsigned k = 0; k < aSignature.size(); ++k) {
//		if (mInverseIndex[k].count(aSignature[k]) > 0) {
//			vector<double> md(mInverseIndex[k][aSignature[k]].begin(),mInverseIndex[k][aSignature[k]].end());
//			std::valarray<double> t(md.data(),hist.size() );
//			std::valarray<double> ti = t.apply(indicator);
//			hist += ti;
//		} else {
//			emptyBins++;
//		}
//	}
//}


//bool HistogramIndex::readBinaryIndex(string filename, indexTy &index){
//	igzstream fin;
//	fin.open(filename.c_str());
//	unsigned numHashFunc = 0;
//	fin.read((char*) &numHashFunc, sizeof(unsigned));
//	if (numHashFunc <= 0)
//		fin.setstate(std::ios::badbit);
//	if (!fin.good())
//		return false;
//	index.resize(numHashFunc);
//	for (unsigned  hashFunc = 0; hashFunc < numHashFunc; hashFunc++){
//
//		unsigned numBins = 0;
//		fin.read((char*) &numBins, sizeof(unsigned));
//		if (numBins < 0)
//			fin.setstate(std::ios::badbit);
//		if (!fin.good())
//			return false;
//		for (unsigned  bin = 0; bin < numBins; bin++){
//
//			unsigned binId = 0;
//			unsigned numBinEntries = 0;
//			fin.read((char*) &binId, sizeof(unsigned));
//			fin.read((char*) &numBinEntries, sizeof(unsigned));
//			if (numBinEntries < 0 || binId <= 0)
//				fin.setstate(std::ios::badbit);
//			if (!fin.good())
//				return false;
//			indexBinTy tmp(numBinEntries);
//			index[hashFunc].insert(make_pair(binId,tmp));
//
//			for (unsigned entry = 0; entry < numBinEntries; entry++ ){
//
//				binTy t = 0;
//				fin.read((char*) &t, sizeof(binTy));
//				if (t < 0)
//					fin.setstate(std::ios::badbit);
//				if (!fin.good())
//					return false;
//				index[hashFunc][binId][entry]=t;
//			}
//		}
//	}
//	fin.close();
//	return true;
//}



//void HistogramIndex::UpdateInverseIndexHist(vector<unsigned>& aSignature, unsigned aIndex) {
//	for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k) { //for every hash value
//		unsigned key = aSignature[k];
//		if (key != MAXUNSIGNED && key != 0) { //if key is equal to markers for empty bins then skip insertion instance in data structure
//			if (mInverseIndex[k].count(key) == 0) { //if this is the first time that an instance exhibits that specific value for that hash function, then store for the first time the reference to that instance
//				indexBinTy tmp(mIndexDataSets.size(),0);
//				tmp[aIndex]++;
//				mInverseIndex[k].insert(make_pair(key, tmp));
//				numKeys++; // just for bin statistics
//			} else if (mInverseIndex[k][key][aIndex]<MAXBIN){
//				mInverseIndex[k][key][aIndex]++;
//			}
//		}
//	}
//}


//void MinHashEncoder::CleanUpInverseIndex() {
//	for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k) {
//		for (indexSingleTy::const_iterator jt = mInverseIndex[k].begin(); jt != mInverseIndex[k].end(); ++jt) {
//			unsigned hash_id = jt->first;
//			if (hash_id != 0 && hash_id != MAXUNSIGNED) { //do not consider buckets corresponding to null bins
//				unsigned collision_size = mInverseIndex[k][hash_id].size();
//
//				if (collision_size < mpParameters->mMaxSizeBin) {
//				} else {//remove bins that are too full from inverse index
//					mInverseIndex[k].erase(hash_id);
//				}
//			}
//		}
//	}
//}

//vector<unsigned> MinHashEncoder::ComputeHashSignatureSize(vector<unsigned>& aSignature) {
//	vector<unsigned> signature_size(mpParameters->mNumHashFunctions);
//	assert(aSignature.size()==mpParameters->mNumHashFunctions);
//	for (unsigned i = 0; i < aSignature.size(); ++i) {
//		unsigned key = aSignature[i];
//		signature_size[i] = mInverseIndex[i][key].size();
//	}
//	return signature_size;
//}

//void HistogramIndex::writeBinaryIndex(ostream &out, const indexTy& index) {
//	// create binary reverse index representation
//	// format:
//	unsigned numHashFunc = index.size();
//	out.write((const char*) &numHashFunc, sizeof(unsigned));
//	for (indexTy::const_iterator it = index.begin(); it!= index.end(); it++){
//		unsigned numBins = it->size();
//		out.write((const char*) &numBins, sizeof(unsigned));
//		for (indexSingleTy::const_iterator itBin = it->begin(); itBin!=it->end(); itBin++){
//			unsigned binId = itBin->first;
//			unsigned numBinEntries = itBin->second.size();
//			out.write((const char*) &binId, sizeof(unsigned));
//			out.write((const char*) &numBinEntries, sizeof(unsigned));
//			for (indexBinTy::const_iterator binEntry = itBin->second.begin(); binEntry != itBin->second.end(); binEntry++){
//				binTy t= *binEntry;
//				out.write((const char*) &(t), sizeof(binTy));
//			}
//		}
//	}
//}

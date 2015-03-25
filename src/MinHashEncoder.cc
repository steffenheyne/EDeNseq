#include "MinHashEncoder.h"

MinHashEncoder::MinHashEncoder() :
mpParameters(NULL), mpData(NULL) {
}

MinHashEncoder::MinHashEncoder(Parameters* apParameters, Data* apData) {
	Init(apParameters, apData);
}

void MinHashEncoder::Init(Parameters* apParameters, Data* apData) {
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

bool MinHashEncoder::SetGraphFromFASTAFile(istream& in, GraphClass& oG, string& currSeq) {
	vector<bool> vertex_status(5, false);
	vertex_status[0] = true; //kernel point
	vertex_status[1] = true; //kind
	vertex_status[2] = true; //viewpoint
	vertex_status[3] = false; //dead
	vertex_status[4] = false; //abstraction

	bool success_status = false;

	unsigned win=mpParameters->mSeqWindow;
	unsigned shift = (unsigned)((double)win*mpParameters->mSeqShift);

	if (currSeq.size() == 0){
		in >> std::ws;
	}

	char c = in.peek();
	string header;
	bool newSeq=false;
	//cout << "here "<< " " << currSeq.size() << " " << in.eof() << " c "  << c << endl;
	if (!in.eof() && c != EOF && c=='>' && currSeq.size() == 0 ){
		getline(in, header,'>');
		getline(in, header);
		getline(in, currSeq,'>');

		currSeq.erase(std::remove(currSeq.begin(), currSeq.end(), '\n'),currSeq.end());
		currSeq.erase(std::remove(currSeq.begin(), currSeq.end(), ' '),currSeq.end());
		std::transform(currSeq.begin(), currSeq.end(), currSeq.begin(), ::toupper);
		in.unget();
		if (currSeq.size()==0 || header.size()==0)
			throw range_error("ERROR FASTA reader - empty Sequence or header found! Header:"+header);
		//cout << " found seq " << header << " " << currSeq.size() << " length" << " EOF? "<< in.eof() << endl;
		newSeq=true;
	} else if (c != '>' && c != EOF && c!= '\n') {
		throw range_error("ERROR FASTA format error!");
	}

	if (currSeq.size() > 0 ) {

		// default case for window/shift
		unsigned currSize = win;
		// case now window/shift
		if (win==0){
			currSize = currSeq.size();
		} else if (win>currSeq.size()) {
			// case seq left is smaller than win
			currSize=currSeq.size();
		}

		if (currSize>=win){

			string graphSeq = currSeq.substr(0,currSize);
			//cout << graphSeq << " " << currSeq.size() << " " << currSize << " " << in.eof() << endl;
			unsigned real_vertex_index = oG.InsertVertex();
			vector<string> vertex_symbolic_attribute_list(1);
			vertex_symbolic_attribute_list[0] = graphSeq;
			oG.SetVertexSymbolicAttributeList(real_vertex_index, vertex_symbolic_attribute_list);
			oG.SetVertexStatusAttributeList(real_vertex_index, vertex_status);
		} else if (newSeq){
			throw range_error("ERROR FASTA reader! Too short sequence found. Either use win=0 or increase window size!");
		}

		if (win>0 && currSeq.size()-shift>=win){
			currSeq.erase(0,shift);
		} else if (currSeq.size()-shift<win || win == 0)
			currSeq="";
		success_status = true;
	}
	return success_status;
}

void MinHashEncoder::worker_readFiles(int numWorkers){

	while (!done){

		SeqDataSet myData;
		myData.filename = "";
		unique_lock<mutex> lk(mut1);
		cv1.wait(lk,[&]{if ( (done) || (readFile_queue.try_pop(myData))) return true; else return false;});
		lk.unlock();
		if (!done && myData.filename != ""){
			cout << endl << "read next file " << myData.filename << " index " << myData.idx << endl;
			igzstream fin;
			fin.open(myData.filename.c_str());

			if (!fin)
				throw range_error("ERROR Data::LoadData: Cannot open file: " + myData.filename);

			bool valid_input = true;
			int file_instances = 0;
			string currSeq;
			while (!fin.eof() && valid_input ) {

				unsigned maxB = max(100,(int)log2((double)signature_counter)*50);
				unsigned currBuff = rand()%(maxB*2 - maxB + 1) + maxB;

				graphQueueS grSet;
				grSet.gr.resize(currBuff);
				grSet.id = myData.idx;
				grSet.offset= file_instances;

				unsigned i = 0;
				while (i < currBuff && !fin.eof() && valid_input) {

					switch (myData.filetype) {
					case FASTA:
						SetGraphFromFASTAFile(fin, grSet.gr.at(i),currSeq);
						break;
					case STRINGSEQ:
						mpData->SetGraphFromFile(fin, grSet.gr.at(i));
						break;
					default:
						throw range_error("ERROR Data::LoadData: file type not recognized: " + myData.filetype);
					}

					if (grSet.gr.at(i).IsEmpty()) {
						valid_input = false;
					} else {
						i++;
						lk.lock();
						instance_counter++;
						lk.unlock();
						file_instances++;
					}
				}
				grSet.gr.resize(i);

				//cout << " instances read " << instance_counter << " buffer " << currBuff << " full..." << graph_queue.size() << endl;
				if (graph_queue.size()>=numWorkers*5){
					unique_lock<mutex> lk(mut2);
					cv2.wait(lk,[&]{if ((done) || (graph_queue.size()<=numWorkers*2)) return true; else return false;});
					lk.unlock();
				}
				graph_queue.push(std::make_shared<graphQueueS>(grSet));
				cv2.notify_all();
			}
			fin.close();
			files_done++;
			cout << " instances produced from file " << file_instances << endl;
		}
	}
}

void MinHashEncoder::worker_Graph2Signature(int id){

	while (!done){

		graphQueueT myGraphSet;
		unique_lock<mutex> lk(mut2);
		cv2.wait(lk,[&]{if ( (done) || (graph_queue.try_pop( (myGraphSet) ))) return true; else return false;});
		lk.unlock();

		if (!done && myGraphSet->gr.size()>0) {
			cv2.notify_all();
			sigQueueS sigSet;
			sigSet.sigs.resize(myGraphSet->gr.size());
			sigSet.offset = myGraphSet->offset;
			sigSet.id = myGraphSet->id;

			//cout << "  graph2sig thread got chunk " << myGraphSet->gr.size() << " offset " << myGraphSet->offset << " threadID " << " id " << id << endl;

			for (unsigned j = 0; j < myGraphSet->gr.size(); j++) {
				SVector x(pow(2, mpParameters->mHashBitSize));
				generate_feature_vector(myGraphSet->gr.at(j), x);
				//mpData->mKernel.GenerateFeatureVector(myGraphSet->gr.at(j), x, id);
				sigSet.sigs[j] = ComputeHashSignature(x);
			}

			sig_queue.push(std::make_shared<sigQueueS>(sigSet));
			cv3.notify_all();
		}
	}
}

void MinHashEncoder::finisher(bool idxUpdate, vector<vector<unsigned> >* myC){
	ProgressBar progress_bar(1000);
	while (!done){

		sigQueueT mySigSet;
		unique_lock<mutex> lk(mut3);
		cv3.wait(lk,[&]{if ( (done) || (sig_queue.try_pop( (mySigSet) ))) return true; else return false;});
		lk.unlock();

		if (!done && mySigSet->sigs.size()>0) {

			if (idxUpdate){
				for (unsigned j = 0; j < mySigSet->sigs.size(); j++) {
					UpdateInverseIndex(mySigSet->sigs[j], mySigSet->offset+j);
				}
			}

			if (myC){
				if (myC->size()<mySigSet->offset+mySigSet->sigs.size())
					myC->resize(mySigSet->offset+mySigSet->sigs.size());
				for (unsigned j = 0; j < mySigSet->sigs.size(); j++) {
					myC->at(mySigSet->offset+j) = mySigSet->sigs[j];
				}
			}
			signature_counter += mySigSet->sigs.size();
			cvm.notify_all();
			//cout << "    finisher updated index with " << mySigSet->sigs.size() << " signatures all_sigs=" <<  signature_counter << " inst=" << instance_counter << endl;
		}
		progress_bar.Count(signature_counter);
	}
}

void MinHashEncoder::LoadDataIntoIndexThreaded(vector<SeqDataSet> myFiles, bool useMinHashCache, vector<vector<unsigned> >* myCache){

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

	cout << "Computing MinHash signatures on the fly while reading " << myFiles.size() << " files..." << endl;
	bool updateIndex = true;
	vector<vector<unsigned> >* myMinHashCache;
	{
		unique_lock<mutex> lk(mutm);

		for (unsigned i=0;i<myFiles.size(); i++){
			readFile_queue.push(myFiles[i]);
		}
		// use external MinHash cache
		if ( (useMinHashCache) && (myCache != NULL)) {
			myMinHashCache = myCache;
		} else if (useMinHashCache){
			// use standard cache
			myMinHashCache = &mMinHashCache;
		} else
			// no cache at all, eg. when only index building
			myMinHashCache = NULL;
	}

	//create all threads
	threads.push_back( std::thread(&MinHashEncoder::finisher,this,updateIndex,myMinHashCache));
	for (int i=0;i<graphWorkers;i++){
		threads.push_back( std::thread(&MinHashEncoder::worker_Graph2Signature,this,i));
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

	cout << endl << "Inverse index ratio of overfull bins (maxSizeBin): " << ((double)numFullBins)/((double)numKeys) << endl;

	if (instance_counter > 0){
		mpData->SetDataSize(instance_counter);
		mpData->mDataIsLoaded = true;
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
	if (mMinHashCache[aID].size() > 0)
		return mMinHashCache[aID];
	else {
		throw range_error("ERROR cache should be filled!");
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
	const int MIN_BINS_IN_COMMON = 1; //Minimum number of bins that two instances have to have in common in order to be considered similar
	//given a list of neighbors with an associated occurrences count, return only a fraction of the highest count ones

	density = 0;
	double myC = (mpParameters->mPureApproximateSim * (mpParameters->mNumHashFunctions - collisions));
	//if (myC<MIN_BINS_IN_COMMON) myC=MIN_BINS_IN_COMMON;
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

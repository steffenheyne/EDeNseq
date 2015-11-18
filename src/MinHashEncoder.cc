#include "MinHashEncoder.h"
#include "Utility.h"



MinHashEncoder::~MinHashEncoder(){

}

MinHashEncoder::MinHashEncoder(Parameters* apParameters, Data* apData)
{
	Init(apParameters, apData);
}


void MinHashEncoder::Init(Parameters* apParameters, Data* apData) {
	mpParameters = apParameters;
	mpData = apData;
	numKeys = 0;
	numFullBins = 0;

	cout << "MAXUNSIGNED "<< MAXUNSIGNED << endl;

	if (mpParameters->mMinRadius>mpParameters->mRadius  || mpParameters->mMinDistance>mpParameters->mDistance) {
		throw range_error("Radius and Distance cannot be smaller than minRadius/minDistance!");
	}

	if (mpParameters->mNumRepeatsHashFunction == 0 || mpParameters->mNumRepeatsHashFunction > mpParameters->mNumHashShingles * mpParameters->mNumHashFunctions){
		mpParameters->mNumRepeatsHashFunction = mpParameters->mNumHashShingles * mpParameters->mNumHashFunctions;
	}

	while ( (mpParameters->mNumHashFunctions * mpParameters->mNumHashShingles) % mpParameters->mNumRepeatsHashFunction != 0){
		--mpParameters->mNumRepeatsHashFunction;
	}

	cout << "Parameter num_repeat_hash_functions adjusted to " <<  mpParameters->mNumRepeatsHashFunction << endl;

	if (mpParameters->mSeqWindow != 0 && mpParameters->mSeqShift == 0){
		throw range_error("Please provide seq_shift > 0 if seq_window is > 0!");
	}
}

inline vector<unsigned> MinHashEncoder::HashFuncNSPDK(const string& aString, unsigned& aStart, unsigned& aMinRadius, unsigned& aMaxRadius, unsigned aBitMask) {
	unsigned int hash = 0xAAAAAAAA;
	unsigned effective_end = min((unsigned) aString.size() - 1, aStart + aMaxRadius);
	unsigned radius = 0;
	vector<unsigned> code_list(aMaxRadius + 1, 0);
	for (std::size_t i = aStart; i <= effective_end; i++) {
		hash ^= ((radius & 1) == 0) ? ((hash << 7) ^ aString[i] * (hash >> 3)) : (~(((hash << 11) + aString[i]) ^ (hash >> 5)));
		code_list[radius] = hash & aBitMask;
		radius++;
	}
	return code_list;
}

inline void  MinHashEncoder::generate_feature_vector(const string& seq, SVector& x) {

	//x.resize(MAXUNSIGNED);
	//x.reserve(500);
	//assume there is a mMinRadius and a mMinDistance
	unsigned& mRadius = mpParameters->mRadius;
	unsigned& mDistance = mpParameters->mDistance;
	unsigned& mMinRadius = mpParameters->mMinRadius;
	unsigned& mMinDistance = mpParameters->mMinDistance;

	unsigned size = seq.size();

	vector<vector<unsigned> > mFeatureCache;
	for (unsigned i = 0; i < size; ++i) {
		mFeatureCache.push_back(vector<unsigned>(mRadius, 0));
	}

	//create neighborhood features
	for (unsigned start = 0; start < size; ++start)
		mFeatureCache[start] = HashFuncNSPDK(seq, start, mMinRadius, mRadius, MAXUNSIGNED);

	vector<unsigned> endpoint_list(3);
	for (unsigned r = mMinRadius; r <= mRadius; r++) {
		//	endpoint_list[0] = r;
		for (unsigned d = mMinDistance; d <= mDistance; d++) {
			endpoint_list[0] = d;
			for (unsigned start = 0; start < size-r-d; ++start) {
				endpoint_list[1] = mFeatureCache[start][r];
				endpoint_list[2] = mFeatureCache[start + d][r];
				//cout << start << " " << start+d << "   " << endpoint_list[2] << "  " << endpoint_list[3]<< endl;
				//  0 1 2 3 4   5 6 7 8 9   r=4 d=5
				//unsigned code = HashFunc(endpoint_list, MAXUNSIGNED);
				//x.coeffRef(code) = 1;
				x.push_back(HashFunc(endpoint_list, MAXUNSIGNED));
				//x.insert(code);
			}
		}
	}
}

void MinHashEncoder::worker_readFiles(unsigned numWorkers, unsigned chunkSizeFactor){

	while (!done){

		SeqFileP myData;

		bool succ = readFile_queue.try_pop(myData);

		if (!done && succ && myData->filename != ""){

			cout << endl << "read next file " << myData->filename << " sig_all_counter " << mSignatureCounter << " inst_counter "<< mInstanceCounter  << endl <<endl;
			igzstream fin;
			fin.open(myData->filename.c_str(),std::ios::in);

			if (!fin)
				throw range_error("ERROR Data::LoadData: Cannot open file: " + myData->filename);

			std::tr1::unordered_map<string, uint8_t> seq_names_seen;

			unsigned pos = 0; // tracks the current seq start pos (window/shift)
			unsigned end = 0; // tracks the current seq end pos, set from BED entry or to full seq end
			unsigned idx = 0; // holds the instance id for the inverse index

			bool valid_input = false; // set to false so that we get new seq in while further down directly
			string currSeq;
			string currFullSeq;
			string currSeqName;

			std::pair<Data::BEDdataIt,Data::BEDdataIt> annoEntries;
			Data::BEDdataIt it; // iteratore over all bed entries of current seq

			uint curr_q = 0;

			while (!fin.eof()) {

				unsigned maxB = max((uint)1000,(uint)log2((double)mSignatureCounter)^2*chunkSizeFactor);
				unsigned currBuff = maxB*3; //rand()%(maxB*4	 - maxB*2 + 1) + maxB; // curr chunk size
				unsigned i = 0;			// current fragment in currBuff
				bool lastSeqGr = false; // indicates that we have the last fragment from current seq, used to get all fragments from current seq into current chunk
				// necessary to have all fragments for one seq/feature if we want to combine signatures in finisher

				ChunkP 		myChunkP = std::make_shared<ChunkT>();
				myChunkP->reserve(currBuff);
				while ( ((i<currBuff) && !fin.eof()) || (myData->signatureAction==CLASSIFY && i>=currBuff && lastSeqGr == false) ) {

					//					cout << "valid? " << valid_input << " name :" << currSeqName << ": pos " << pos << " end " << end <<  endl;
					if (!valid_input) {
						if  ( it == annoEntries.second ) {
							// last seq and all bed entries for it are finished, get next seq from file

							switch (myData->filetype) {
							case FASTA:
								mpData->GetNextFastaSeq(fin, currFullSeq, currSeqName);
								if (fin.eof() )
									continue;
								mSequenceCounter++;
								if (myData->checkUniqueSeqNames && seq_names_seen.count(currSeqName) > 0) {
									throw range_error("Sequence names are not unique in FASTA file! "+currSeqName);
								} else if (myData->checkUniqueSeqNames) {
									seq_names_seen.insert(make_pair(currSeqName,1));
								}
								break;
							case STRINGSEQ:
								mpData->GetNextStringSeq(fin, currFullSeq);
								if (fin.eof() )
									continue;
								mSequenceCounter++;
								currSeqName =  std::to_string(mSequenceCounter);
								break;
							default:
								throw range_error("ERROR Data::LoadData: file type not recognized: " + myData->filetype);
							}

							// log output
							if (myData->signatureAction==INDEX){
								//		cout << endl << " next found Seq #" <<  seq_names_seen.size() << " length " << currFullSeq.size() << ":" << currSeqName << ": " << endl;
							}

							// if we have bed entries for a seq, find them and set iterator to first bed entry
							if (myData->dataBED && myData->dataBED->find(currSeqName) != myData->dataBED->end()){
								annoEntries = myData->dataBED->equal_range(currSeqName);
								it = annoEntries.first;
							} else if (myData->dataBED){
								//	cout << "no bed entry found! "<< seq_names_seen.size()<< endl;
								// bed is present, but no entry for current seq found -> we take next seq
								valid_input = false;
								continue;
							}
						} // if no bed entries left for current seq get new seq

						// check if we use the same idx-group for the whole seq, either by seq name or feature id from BED
						// idx also defines the value under which we insert features into the index
						switch (myData->groupGraphsBy){
						// use seq name as value for inverse index
						case SEQ_NAME:
							if (mFeature2IndexValue.find(currSeqName) != mFeature2IndexValue.end()){
								idx = mFeature2IndexValue[currSeqName];
							} else {
								myData->lastMetaIdx++;
								idx=myData->lastMetaIdx;
								mFeature2IndexValue.insert(make_pair(currSeqName,idx));
							}
							break;
							// use given value/name in BED file col4 as  value for inverse index
						case SEQ_FEATURE:
							if (mFeature2IndexValue.find(it->second->NAME) != mFeature2IndexValue.end()){
								idx = mFeature2IndexValue[it->second->NAME];
							} else {
								myData->lastMetaIdx++;
								idx=myData->lastMetaIdx;
								mFeature2IndexValue.insert(make_pair(it->second->NAME,idx));
							}
							break;
						default:
							break;
						}

						// only true if we have a found a BED entry for current seq
						// set region according to BED entry
						if ( it != annoEntries.second ) {
							pos = it->second->START;
							end = it->second->END;
							//mIndexValue2Feature.insert(make_pair(idx,it->second));
							//cout << endl << "BED entry found for seq name " << currSeqName << " " << it->second->NAME << " MetaIdx "<< idx << " " << pos << "-"<< end << endl;
							it++;
						} else {
							// no bed is present, then we set start/end to full seq, eg. in case for clustering
							pos=0;
							end=currFullSeq.size();
							//cout << "no BED data present! "<< currSeqName << " " << pos << "-" << end << endl;
						}

						// check if start/end is within bounds of found seq
						if (pos>currSeq.size() || end > currFullSeq.size())
							throw range_error(" BED entry start/end is outside current seq ");

						// apply current seq start/end
						currSeq = currFullSeq.substr(pos,end-pos);

					} // valid_input?

					// new instance for this chunk
					InstanceT	myInstance;
					mpData->GetNextWinFromSeq(currSeq, pos, lastSeqGr,myInstance.seq);

					if (myInstance.seq.size() == 0 && !lastSeqGr) {
						valid_input = false;
					} else {
						// fill current Instance with all data
						// make graph from seq

						valid_input = true;

						switch (myData->groupGraphsBy){
						case NONE:
						case SEQ_WINDOW:
							idx = mInstanceCounter+1;
							break;
						default:
							break;
						}

						if (myData->strandType != REV){

							myInstance.seqFile = myData;
							myInstance.name = currSeqName;
							myInstance.idx = idx;
							myInstance.pos = pos;
							myInstance.rc = false;

							myChunkP->push_back(myInstance);
							i++;
							mInstanceCounter++;
						}

						if (myData->strandType != FWD){
							InstanceT	myInstanceRC;
							myInstanceRC.seqFile = myData;
							myInstanceRC.name = currSeqName;
							myInstanceRC.idx = idx;
							myInstanceRC.pos = pos;
							myInstanceRC.rc = true;
							mpData->GetRevComplSeq(myInstance.seq,myInstanceRC.seq);

							myChunkP->push_back(myInstanceRC);
							mInstanceCounter++;
							i++;
						}

					}
					//cout << "Gr: " << myChunkP->size() << " "<< i << " " << currBuff<< " "<< pos << " " << currSeqName<<  " " << currSeq.size() << " " << lastSeqGr << endl;
				} // while buffer not full or eof

				//cout << "Gr: " << myChunkP->size() << " "<< i << " " << currBuff<< " "<< pos << " " << currSeqName<<  " " << currSeq.size() << " " << lastSeqGr << endl;
				if (i==0)
					continue;

				graph_queue[curr_q].push(myChunkP);

				unsigned fillstatus = MAXUNSIGNED;
				for (unsigned i=0; i < graph_queue.size(); ++i){
					if ((uint)graph_queue[i].size() < fillstatus){
						fillstatus = graph_queue[i].size();
						curr_q = i;
					}
				}

				if (fillstatus>min((uint)20,(uint)(graph_queue.size()*2))){
					unique_lock<mutex> lk(mut1);
					cv1.wait(lk,[&]{fillstatus = MAXUNSIGNED;for (uint i=0; i<graph_queue.size(); ++i){ if ((uint)graph_queue[i].size()<fillstatus) { fillstatus = graph_queue[i].size();}}; if ((done) || (fillstatus<(uint)max((uint)10,(uint)graph_queue.size()))) return true; else return false;});
					lk.unlock();
				}

			} // while eof
			fin.close();
			files_done++;
			//cout << endl << "file " << files_done << " seqs " << file_seqs << " " << mInstanceCounter << " " << mSignatureCounter << " instances produced from file " << file_instances << endl;
		}
	}
}

void MinHashEncoder::worker_Graph2Signature(int numWorkers, unsigned id){

	Signature* tmpSig = new Signature(numHashFunctionsFull);

	while (!done){

		ChunkP myData;
		vector<ChunkP> myQ;

		// take multiple chunks at once, gives better throughput
		while (graph_queue[id].size()>=1){
			graph_queue[id].wait_and_pop(myData);
			myQ.push_back(myData);
		}

		while (!done && myQ.size()){
			myData = myQ.back();
			myQ.pop_back();
			//cout << "  graph2sig thread got chunk " << myData->size() << endl;
			for (unsigned j = 0; j < myData->size(); j++) {
				generate_feature_vector((*myData)[j].seq, (*myData)[j].svec);
				ComputeHashSignature((*myData)[j].svec,(*myData)[j].sig,tmpSig);
			}
			{
				lock_guard<mutex> lk(mut2);

			for (unsigned i=0;i<index_queue.size();i++){
				index_queue[i].push(myData);
			}
			}
			uint fillstatus=0;
			for (uint i=0; i<index_queue.size(); ++i){ fillstatus += index_queue[i].size();}

			if (fillstatus>index_queue.size()*40){
				unique_lock<mutex> lk(mut2);
				cv2.wait(lk,[&]{fillstatus = 0;for (uint i=0; i<index_queue.size(); ++i){ fillstatus += index_queue[i].size();} if ((done) || (fillstatus<index_queue.size()*40)) return true; else return false;});
				lk.unlock();
			}
		}
	}
	delete tmpSig;
}


void MinHashEncoder::finisher_IndexUpdate(unsigned id, unsigned min, unsigned max){
	ProgressBar progress_bar(1000);
	while (!done){

		ChunkP myData;
		bool succ;
		{
			lock_guard<mutex> lk(mut2);
			succ = index_queue[id].try_pop(myData);
		}
		if (!done && succ && myData->size()>0) {

			finishUpdate(myData,min,max);
			mSignatureUpdateCounter += myData->size()*( (max-min+1));
			//	cout << "finishUpdate "<< myData->size() << " " << id << " " << min << " " << max << " " << index_queue[id].size() <<  " " << myData->size()*( (max-min+1)*mpParameters->mNumHashFunctions) << " " << mSignatureUpdateCounter << endl;

			if (id==0){
				mSignatureCounter += myData->size();
				double elap = progress_bar.getElapsed()/1000;
				cout.setf(ios::fixed);
				cout << "\r" <<  std::setprecision(1) << elap << " sec elapsed  numSeqs=" << std::setprecision(0) << setw(10);
				cout << mSequenceCounter  << "("<<mSequenceCounter/(elap) <<" seq/s)  sign=" << setw(10);
				cout << mSignatureCounter << "("<<(double)mSignatureCounter/((elap)) <<" sig/s - "<< (double)(mSignatureCounter/graph_queue.size())/((elap)) << " per thread) inst=";
				cout << mInstanceCounter << " graphQueue=";
				uint avg = 0;
				for (uint i=0; i<graph_queue.size(); ++i){
					//cout << graph_queue[i].size() << " ";
					avg += graph_queue[i].size();
				};
				cout << avg/graph_queue.size();
				cout << " indexQueue=";
				avg = 0;
				for (uint i=0; i<index_queue.size(); ++i){
					//cout << index_queue[i].size() << " ";
					avg += index_queue[i].size();
				}
				cout << avg;
			}
		}
	}
}


void MinHashEncoder::LoadData_Threaded(SeqFilesT& myFiles){

	for (unsigned i=0;i<myFiles.size(); i++){
		readFile_queue.push(myFiles[i]);
	}

	HashSignatureHelper();

	cout << "Using " << mHashBitMask << " \t as bitmask"<< endl;
	cout << "Using " << mpParameters->mHashBitSize << " \t bits to encode features" << endl;
	cout << "Using " << mpParameters->mRandomSeed << " \t as random seed" << endl;
	cout << "Using " << mpParameters->mNumHashFunctions << " \t final hash VALUES (param 'num_hash_functions' is used as FINAL signature size!), " << endl;
	cout << "Using " << mpParameters->mNumHashFunctions*mpParameters->mNumHashShingles << " \t hash values internally (num_hash_functions * num_hash_shingles) from " << mpParameters->mNumRepeatsHashFunction << " independent hash functions (param 'num_repeat_hash_functions')" << endl;
	cout << "Using " << mpParameters->mNumHashShingles << " \t as hash shingle factor (param 'num_hash_shingles')" << endl;
	cout << "Using " << mpParameters->mPureApproximateSim << " \t as required approximate similarity (param 'pure_approximate_sim')" << endl;
	cout << "Using feature radius   " << mpParameters->mMinRadius<<".."<<mpParameters->mRadius << endl;
	cout << "Using feature distance " << mpParameters->mMinDistance<<".."<<mpParameters->mDistance << endl;
	cout << "Using sequence window  " << mpParameters->mSeqWindow<<" shift "<<mpParameters->mSeqShift << " nt - clip " << mpParameters->mSeqClip << endl;

	cout << endl << "Computing MinHash signatures on the fly while reading " << myFiles.size() << " file(s)..." << endl;

	// threaded producer-consumer model for signature creation and index update
	// created threads:
	// 	 1 worker_readFiles thread that produces chunks of sequences/windows and put these into the n graph_queues
	//	 n worker_Graph2Signature threads that create the signatures and put them into the m index_queues
	//   m worker_IndexUpdate threads, each updates a range of mNumHashFunctions so that we can update the index in "parallel"

	int graphWorkers = std::thread::hardware_concurrency();
	if (mpParameters->mNumThreads>0)
		graphWorkers = mpParameters->mNumThreads;

	cout << "Using " << graphWorkers << " worker threads and " << 1+min(max((unsigned)1,mpParameters->mNumIndexThreads),mpParameters->mNumHashFunctions) << " helper threads..." << endl;

	done 					= false;
	files_done				= 0;
	mSignatureCounter		= 0;
	mInstanceCounter 		= 0;
	mSequenceCounter 		= 0;
	mSignatureUpdateCounter = 0;

	vector<std::thread> threads;
	graph_queue.resize(graphWorkers);

	// create m worker_IndexUpdate threads
	// distribute mpParameters->mNumHashFunctions as equal as possible between the
	// requested number of mpParameters->mNumIndexThreads, adjust numIndexThreads to have maximal 1 thread per signature value
	unsigned numIndexThreads = max((unsigned)1,mpParameters->mNumIndexThreads);
	index_queue.resize(min(numIndexThreads,mpParameters->mNumHashFunctions));
	unsigned hf_left = mpParameters->mNumHashFunctions;
	for (unsigned m=min(numIndexThreads,mpParameters->mNumHashFunctions);m>=1;m--){
		unsigned range = std::ceil((double)hf_left / (double)(m));
		// launch thread with range over mNumHashFunctions
		threads.push_back( std::thread(&MinHashEncoder::finisher_IndexUpdate,this,m-1,hf_left-range,hf_left-1));
		cout << "index update thread " << m << " range " << hf_left-range+1 << ".." << hf_left << endl;
		hf_left -= range;
	}

	// create worker_readFiles thread
	threads.push_back( std::thread(&MinHashEncoder::worker_readFiles,this,graphWorkers, 500));
	std::this_thread::sleep_for(std::chrono::milliseconds(500));

	// create n worker_Graph2Signature threads
	for (int n=0;n<graphWorkers;n++){
		threads.push_back( std::thread(&MinHashEncoder::worker_Graph2Signature,this,graphWorkers,n));
	}

	{
		join_threads joiner(threads);

		while(!done){

			std::this_thread::sleep_for(std::chrono::milliseconds(20));
			if ( (files_done<myFiles.size()) || (mSignatureUpdateCounter < mInstanceCounter*mpParameters->mNumHashFunctions))
				done=false;
			else done=true;
			cv2.notify_all();
			cv1.notify_all();
		}

	} // by leaving this block threads get joined by destruction of joiner

	cout << "SignatureUpdateCounter="<< mSignatureUpdateCounter << endl;

	// threads finished
	if (numKeys>0)
		cout << endl << "Inverse index ratio of overfull bins (maxSizeBin): " << ((double)numFullBins)/((double)numKeys) << " "<< numFullBins << "/" << numKeys << " instances " << mInstanceCounter << endl;

	if (mInstanceCounter == 0) {
		throw range_error("ERROR in MinHashEncoder::LoadData: something went wrong; no instances/signatures produced");
	} else
		cout << "Instances/signatures produced " << mInstanceCounter << endl;
}


unsigned MinHashEncoder::GetLoadedInstances() {
	return mInstanceCounter;
}

void MinHashEncoder::HashSignatureHelper() {

	// set parameters for ComputeHashSignature
	mHashBitMask = (2 << (mpParameters->mHashBitSize - 1)) - 1;

	// sub_hash_range is always floor(numHashFunctionsFull / mpParameters->mNumRepeatsHashFunction)
	numHashFunctionsFull = mpParameters->mNumHashFunctions * mpParameters->mNumHashShingles;
	sub_hash_range = numHashFunctionsFull / mpParameters->mNumRepeatsHashFunction;

	if (mpParameters->mNumHashShingles == 1){
		mHashBitMask_feature = mHashBitMask;
	} else {
		mHashBitMask_feature = MAXUNSIGNED;
		mHashBitMask_shingle = mHashBitMask;
	}

	mBounds.resize(sub_hash_range+1);
	for (unsigned kk = 0; kk < sub_hash_range; ++kk) { //for all k values
		mBounds[kk] = mHashBitMask_feature / sub_hash_range * kk;
	}
	mBounds[sub_hash_range] = mHashBitMask_feature;
}

void MinHashEncoder::ComputeHashSignature(const SVector& aX, Signature& signature, Signature* tmpSig) {

	// sub_hash_range is always floor(numHashFunctionsFull / mpParameters->mNumRepeatsHashFunction)
	// if we use shingles we need a temp signature of length numHashFunctionsFull
	// otherwise we directly put values in the provided signature object
	Signature *signatureP;
	if (mpParameters->mNumHashShingles>1){
		signatureP = tmpSig;
	} else {
		signature.resize(numHashFunctionsFull);
		signatureP = &signature;
	}

	//init with MAXUNSIGNED
	for (unsigned k = 0; k < numHashFunctionsFull; ++k)
		(*signatureP)[k] = MAXUNSIGNED;
	//prepare a vector containing the signature as the k min values
	//for each element of the sparse vector
	//for (SVector::InnerIterator it(aX); it; ++it) {
	for (SVector::const_iterator it=aX.begin(); it!=aX.end(); ++it) {
		//for each sub_hash
		for (unsigned l = 1; l <= mpParameters->mNumRepeatsHashFunction; ++l) {
			//unsigned key = IntHash(it.index(), mHashBitMask_feature, l);
			unsigned key = IntHash(*it, mHashBitMask_feature, l);
			for (unsigned kk = 0; kk < sub_hash_range; ++kk) { //for all k values
				if (key >= mBounds[kk] && key < mBounds[kk+1]) { //if we are in the k-th slot
					unsigned signature_feature = kk + (l - 1) * sub_hash_range;
					if (key < (*signatureP)[signature_feature]) {//keep the min hash within the slot
						(*signatureP)[signature_feature] = key;
					}
					//cout << MAXUNSIGNED << " "<< mpParameters->mNumRepeatsHashFunction << " " << sub_hash_range << " " << lower_bound <<" " << upper_bound << " " << signature_feature << " " << l << " " << kk << " " <<  key<< endl;
				}
			}
		}
	}

	// compute shingles, i.e. rehash mNumHashShingles hash values into one hash value
	if (mpParameters->mNumHashShingles > 1 ) {
		vector<unsigned> signatureFinal(mpParameters->mNumHashFunctions);
		for (unsigned i=0;i<mpParameters->mNumHashFunctions;i++){
			signatureFinal[i] = HashFunc(signatureP->begin()+(i*mpParameters->mNumHashShingles),signatureP->begin()+(i*mpParameters->mNumHashShingles+mpParameters->mNumHashShingles),mHashBitMask_shingle);
		}
		signature.swap(signatureFinal);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
//
//	CLASS NEIGHBORHOODINDEX
//
///////////////////////////////////////////////////////////////////////////////////////////

void NeighborhoodIndex::NeighborhoodCacheReset() {
	cout << "... nearest neighbor cache reset ..." << endl;
	mNeighborhoodCache.clear();
	mNeighborhoodCacheExt.clear();
	mNeighborhoodCacheInfo.clear();

	mNeighborhoodCache.resize(GetLoadedInstances());
	mNeighborhoodCacheExt.resize(GetLoadedInstances());
	mNeighborhoodCacheInfo.resize(GetLoadedInstances());
}


vector<unsigned>& NeighborhoodIndex::ComputeHashSignature(unsigned aID) {

	if (mMinHashCache && mMinHashCache->size()>0 && mMinHashCache->at(aID).size() > 0)
		return mMinHashCache->at(aID);
	else {
		throw range_error("ERROR: MinHashCache is not filled!");
	}
}

void NeighborhoodIndex::UpdateInverseIndex(vector<unsigned>& aSignature, unsigned& aIndex) {
	const binKeyTy& aIndexT = static_cast<binKeyTy>(aIndex);
	for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k) { //for every hash value
		unsigned key = aSignature[k];
		if (key != MAXUNSIGNED && key != 0) { //if key is equal to markers for empty bins then skip insertion instance in data structure
			if (!mInverseIndex[k][key]) { //if this is the first time that an instance exhibits that specific value for that hash function, then store for the first time the reference to that instance

				indexBinTy foo;
				foo = new binKeyTy[2];
				foo[1]= aIndexT;
				foo[0]= 1; //index of last element is stored at idx[0]

				mInverseIndex[k][key] = foo;
				numKeys++; // just for bin statistics

			} else if (mInverseIndex[k][key][0] < mpParameters->mMaxSizeBin && mInverseIndex[k][key][0] != MAXBINKEY) {

				binKeyTy*& myValue = mInverseIndex[k][key];
				binKeyTy newSize = (myValue[0])+1;
				indexBinTy fooNew = new binKeyTy[newSize+1];
				memcpy(fooNew, myValue, (newSize)*sizeof(binKeyTy));
				fooNew[newSize] = aIndexT;
				fooNew[0] = newSize;
				delete[] myValue;
				myValue = fooNew;

			} else if (mInverseIndex[k][key][0] == mpParameters->mMaxSizeBin){
				// if a bin is full we clear it and add a key with MAXUNSIGNED to indicate an overfull/not useful bin
				numFullBins++; // just for bin statistics
				indexBinTy fooNew = new binKeyTy[1];
				fooNew[0] = MAXBINKEY;
				delete[] mInverseIndex[k][key];
				mInverseIndex[k][key] = fooNew;
			}
		}
	}
}


void NeighborhoodIndex::ComputeApproximateNeighborhoodCore(const vector<unsigned>& aSignature, umap_uint_int& neighborhood, unsigned& collisions) {

	collisions = 0;
	for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k) {

		unsigned hash_id = aSignature[k];
		if (hash_id != 0 && hash_id != MAXUNSIGNED && mInverseIndex[k][hash_id][0] != MAXBINKEY) {

			//fill neighborhood set counting number of occurrences
			for (binKeyTy i = 1; i<=mInverseIndex[k][hash_id][0];i++){
				//	for (indexBinTy::iterator it = mInverseIndex[k][hash_id].begin(); it != mInverseIndex[k][hash_id].end(); ++it) {
				binKeyTy instance_id = (mInverseIndex[k][hash_id][i])-1;
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

vector<unsigned> NeighborhoodIndex::ComputeNeighborhood(const unsigned aID, unsigned& collisions, double& density) {

	vector<unsigned> neighborhood_list;

	if (mNeighborhoodCache[aID].size() != 0) {
		neighborhood_list = mNeighborhoodCache[aID];
		pair<unsigned, double> tmp = mNeighborhoodCacheInfo[aID];
		collisions = tmp.first;
		density = tmp.second;
	} else {
		vector<unsigned> signature = ComputeHashSignature(aID);
		vector<unsigned> neighborhood_list = ComputeApproximateNeighborhood(signature, collisions, density);
		mNeighborhoodCacheInfo[aID] = make_pair(collisions,density);
		mNeighborhoodCache[aID] = neighborhood_list;
	}
	return neighborhood_list;
}

umap_uint_int NeighborhoodIndex::ComputeNeighborhoodExt(unsigned aID, unsigned& collisions, double& density) {
	//cache neighborhoods (if opted for)
	umap_uint_int neighborhood_list;
	if (mNeighborhoodCacheExt[aID].size() != 0) {
		neighborhood_list = mNeighborhoodCacheExt[aID];
		pair<unsigned, double> tmp = mNeighborhoodCacheInfo[aID];
		collisions = tmp.first;
		density = tmp.second;
	} else {
		vector<unsigned> signature = ComputeHashSignature(aID);
		umap_uint_int neighborhood_list = ComputeApproximateNeighborhoodExt(signature, collisions, density);
		mNeighborhoodCacheExt[aID] = neighborhood_list;
		mNeighborhoodCacheInfo[aID] = make_pair(collisions,density);
	}
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
		if (hash_id != 0 && hash_id != MAXUNSIGNED && mInverseIndex[k][hash_id][0] != MAXBINKEY) {

			//fill neighborhood set counting number of occurrences
			for (binKeyTy i = 1; i<=mInverseIndex[k][hash_id][0];i++){
				//for (indexBinTy::iterator it = mInverseIndex[k][hash_id].begin(); it != mInverseIndex[k][hash_id].end(); ++it) {
				if (mInverseIndex[k][hash_id][i] == aID) counts_aID++;
			}

		} else
			collisions++;
	}
	pair<unsigned,unsigned> tmp=make_pair(counts_aID,collisions);
	return tmp;
}

///////////////////////////////////////////////////////////////////////////////////////////
//
//	CLASS HISTOGRAMINDEX
//
///////////////////////////////////////////////////////////////////////////////////////////


void HistogramIndex::SetHistogramSize(binKeyTy size){
	mHistogramSize = size;
}


HistogramIndex::binKeyTy HistogramIndex::GetHistogramSize(){
	return mHistogramSize;
}


void HistogramIndex::InitInverseIndex() {

	mInverseIndex.resize(mpParameters->mNumHashFunctions);

	mMemPool_2.resize(mpParameters->mNumHashFunctions);
	mMemPool_3.resize(mpParameters->mNumHashFunctions);
	mMemPool_4.resize(mpParameters->mNumHashFunctions);
	mMemPool_5.resize(mpParameters->mNumHashFunctions);
	mMemPool_6.resize(mpParameters->mNumHashFunctions);
	mMemPool_7.resize(mpParameters->mNumHashFunctions);
	mMemPool_8.resize(mpParameters->mNumHashFunctions);
	mMemPool_9.resize(mpParameters->mNumHashFunctions);
	mMemPool_10.resize(mpParameters->mNumHashFunctions);

	for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k){
		mInverseIndex[k].max_load_factor(0.5);
		mInverseIndex[k].set_resizing_parameters(0.0,0.5);
		mInverseIndex[k].rehash(134217728);
		//mInverseIndex[k].set_deleted_key(0);
		//mInverseIndex[k].set_empty_key(0);
		mMemPool_2[k] = new MemoryPool<newIndexBin_2,mMemPool_BlockSize>();
		mMemPool_3[k] = new MemoryPool<newIndexBin_3,mMemPool_BlockSize>();
		mMemPool_4[k] = new MemoryPool<newIndexBin_4,mMemPool_BlockSize>();
		mMemPool_5[k] = new MemoryPool<newIndexBin_5,mMemPool_BlockSize>();
		mMemPool_6[k] = new MemoryPool<newIndexBin_6,mMemPool_BlockSize>();
		mMemPool_7[k] = new MemoryPool<newIndexBin_7,mMemPool_BlockSize>();
		mMemPool_8[k] = new MemoryPool<newIndexBin_8,mMemPool_BlockSize>();
		mMemPool_9[k] = new MemoryPool<newIndexBin_9,mMemPool_BlockSize>();
		mMemPool_10[k] = new MemoryPool<newIndexBin_10,mMemPool_BlockSize>();
	}
}

void HistogramIndex::UpdateInverseIndex(const vector<unsigned>& aSignature, const unsigned& aIndex) {
	unsigned min = 0;
	unsigned max = mpParameters->mNumHashFunctions-1;
	UpdateInverseIndex(aSignature, aIndex, min, max);
}


inline void* HistogramIndex::memcpy2(void* dest,const void* src, size_t count){
	binKeyTy* dst8 = (binKeyTy*)dest;
	binKeyTy* src8 = (binKeyTy*)src;

	while (count--){
		*dst8++ = *src8++;
	}
	return dest;
}

void HistogramIndex::UpdateInverseIndex(const vector<unsigned>& aSignature, const unsigned& aIndex, unsigned& min, unsigned& max) {
	//cout << "UpdateInverseIndex "<< aSignature.size() << " " << aIndex << " " << min << " " << max << endl;
	const binKeyTy& aIndexT =(binKeyTy)aIndex;
	for (unsigned k = min; k <= max; ++k) { //for every hash value
		const unsigned& key = aSignature[k];

		if (key != MAXUNSIGNED && key != 0) { //if key is equal to markers for empty bins then skip insertion instance in data structure
			if (mInverseIndex[k].count(key)==0) { //if this is the first time that an instance exhibits that specific value for that hash function, then store for the first time the reference to that instance

				binKeyTy* foo;
				foo = reinterpret_cast<indexBinTy>(mMemPool_2[k]->newElement());
				foo[1] = aIndexT;
				foo[0] = 1; //index of last element is stored at idx[0]

				mInverseIndex[k].rehash((numKeys/mpParameters->mNumHashFunctions)+5000000);
				mInverseIndex[k][key] = foo;
				numKeys++; // just for bin statistics
			} else {

				binKeyTy*& myValue = mInverseIndex[k][key];

				if ( myValue[myValue[0]] != aIndexT){

					// find pos for insert, assume sorted array
					binKeyTy i = myValue[0];
					while ((myValue[i]> aIndexT) && (i>1)){
						i--;
					}

					// only insert if element is not there
					if (myValue[i]<aIndexT){
						binKeyTy newSize = (myValue[0])+1;
						binKeyTy * fooNew;

						switch (newSize){
						case 2:
							fooNew = reinterpret_cast<binKeyTy(*)>(mMemPool_3[k]->newElement());
							break;
						case 3:
							fooNew = reinterpret_cast<binKeyTy(*)>(mMemPool_4[k]->newElement());
							break;
						case 4:
							fooNew = reinterpret_cast<binKeyTy(*)>(mMemPool_5[k]->newElement());
							break;
//						case 5:
//							fooNew = reinterpret_cast<binKeyTy(*)>(mMemPool_6[k]->newElement());
//							break;
//						case 6:
//							fooNew = reinterpret_cast<binKeyTy(*)>(mMemPool_7[k]->newElement());
//							break;
//						case 7:
//							fooNew = reinterpret_cast<binKeyTy(*)>(mMemPool_8[k]->newElement());
//							break;
//						case 8:
//							fooNew = reinterpret_cast<binKeyTy(*)>(mMemPool_9[k]->newElement());
//							break;
//						case 9:
//							fooNew = reinterpret_cast<binKeyTy(*)>(mMemPool_10[k]->newElement());
//							break;
						default:
							fooNew = new binKeyTy[newSize+1];
							break;
						}

						memcpy(fooNew,myValue,(i+1)*sizeof(binKeyTy));
						//memcpy2(fooNew,myValue,(i+1));
						fooNew[i+1] = aIndexT;
						memcpy(&fooNew[i+2],&myValue[i+1],(myValue[0]-i)*sizeof(binKeyTy));
						//memcpy2(&fooNew[i+2],&myValue[i+1],(myValue[0]-i));
						fooNew[0] = newSize;

						switch (myValue[0]) {
						case 1:
							mMemPool_2[k]->deleteElement(reinterpret_cast<newIndexBin_2(*)>(myValue));
							break;
						case 2:
							mMemPool_3[k]->deleteElement(reinterpret_cast<newIndexBin_3(*)>(myValue));
							break;
						case 3:
							mMemPool_4[k]->deleteElement(reinterpret_cast<newIndexBin_4(*)>(myValue));
							break;
						case 4:
							mMemPool_5[k]->deleteElement(reinterpret_cast<newIndexBin_5(*)>(myValue));
							break;
//						case 5:
//							mMemPool_6[k]->deleteElement(reinterpret_cast<newIndexBin_6(*)>(myValue));
//							break;
//						case 6:
//							mMemPool_7[k]->deleteElement(reinterpret_cast<newIndexBin_7(*)>(myValue));
//							break;
//						case 7:
//							mMemPool_8[k]->deleteElement(reinterpret_cast<newIndexBin_8(*)>(myValue));
//							break;
//						case 8:
//							mMemPool_9[k]->deleteElement(reinterpret_cast<newIndexBin_9(*)>(myValue));
//							break;
//						case 9:
//							mMemPool_10[k]->deleteElement(reinterpret_cast<newIndexBin_10(*)>(myValue));
//							break;
						default:
							delete[] myValue;
							break;
						}
						mInverseIndex[k][key] = fooNew;
					}
					//				cout << "bin " << key << " k " <<  k << " aIdx "<< aIndex << "\t";
					//				for (unsigned j=0; j<=mInverseIndex[k][key][0];j++){
					//					cout << mInverseIndex[k][key][j] <<"\t";
					//				}
					//				cout << endl;
				}

			}
		}
	}
}


void HistogramIndex::ComputeHistogram(const vector<unsigned>& aSignature, std::valarray<double>& hist, unsigned& emptyBins) {

	hist.resize(GetHistogramSize());
	hist *= 0;
	emptyBins = 0;
	for (unsigned k = 0; k < aSignature.size(); ++k) {
		if (mInverseIndex[k].count(aSignature[k])!=0) {
		//if (ahmV[k]->count(aSignature[k])!=0) {
			std::valarray<double> t(0.0, hist.size());

			indexBinTy& myValue = mInverseIndex[k][aSignature[k]];
			//indexBinTy& myValue = ahmV[k]->find(aSignature[k])->second;

			//		if (myValue[0]<=1){
			for (unsigned i=1;i<=myValue[0];i++){
				t[myValue[i]-1]=1;
			}
			//		}
			hist += t;

		} else {
			emptyBins++;
		}
	}
}

void HistogramIndex::writeBinaryIndex2(ostream &out, const indexTy& index) {
	// create binary reverse index representation
	// format:
	out.write((const char*) &INDEX_FORMAT_VERSION, sizeof(unsigned));
	out.write((const char*) &mpParameters->mHashBitSize, sizeof(unsigned));
	out.write((const char*) &mpParameters->mRandomSeed, sizeof(unsigned));
	out.write((const char*) &mpParameters->mRadius, sizeof(unsigned));
	out.write((const char*) &mpParameters->mMinRadius, sizeof(unsigned));
	out.write((const char*) &mpParameters->mDistance, sizeof(unsigned));
	out.write((const char*) &mpParameters->mMinDistance, sizeof(unsigned));
	out.write((const char*) &mpParameters->mNumHashShingles, sizeof(unsigned));
	out.write((const char*) &mpParameters->mNumRepeatsHashFunction, sizeof(unsigned));
	out.write((const char*) &mpParameters->mSeqWindow, sizeof(unsigned));
	out.write((const char*) &mpParameters->mIndexSeqShift, sizeof(unsigned));
	unsigned tmp = GetHistogramSize();
	out.write((const char*) &tmp, sizeof(unsigned));

	if (mFeature2IndexValue.size() != GetHistogramSize()){
		throw range_error("Ups! Histogramsize is different to mFeature2IndexValue.size()");
	}

	std::pair<string,uint> highest = *mFeature2IndexValue.rbegin();          // last element
	std::map<string,uint>::iterator it = mFeature2IndexValue.begin();
	do {
		out.write((const char*) &(it->second), sizeof(unsigned));
		tmp = it->first.size();
		out.write((const char*) &(tmp), sizeof(unsigned));
		out.write(const_cast<char*>(it->first.c_str()), it->first.size());
	} while ( mFeature2IndexValue.value_comp()(*it++, highest) );

	unsigned numHashFunc = index.size();
	out.write((const char*) &numHashFunc, sizeof(unsigned));
	for (typename indexTy::const_iterator it = index.begin(); it!= index.end(); it++){
		unsigned numBins = it->size();
		out.write((const char*) &numBins, sizeof(unsigned));
		for (typename indexSingleTy::const_iterator itBin = it->begin(); itBin!=it->end(); itBin++){
			unsigned binId = itBin->first;

			unsigned numBinEntries = itBin->second[0];
			out.write((const char*) &binId, sizeof(unsigned));
			out.write((const char*) &numBinEntries, sizeof(unsigned));

			for (binKeyTy i=1;i<=numBinEntries;i++){
				binKeyTy s = itBin->second[i];
				out.write((const char*) &(s), sizeof(binKeyTy));
			}
		}
	}
}

bool HistogramIndex::readBinaryIndex2(string filename, indexTy &index){
	igzstream fin;
	fin.open(filename.c_str());
	unsigned tmp;
	fin.read((char*) &tmp, sizeof(unsigned));
	if (tmp != INDEX_FORMAT_VERSION) {
		cout << endl << endl << "Incompatible Index file format - Index file has version " << tmp  << ", but this program uses version " << INDEX_FORMAT_VERSION << "!" << endl;
		cout << "Please re-create index with this program!" << endl;
		return false;
	}

	fin.read((char*) &mpParameters->mHashBitSize, sizeof(unsigned));
	fin.read((char*) &mpParameters->mRandomSeed, sizeof(unsigned));
	fin.read((char*) &mpParameters->mRadius, sizeof(unsigned));
	fin.read((char*) &mpParameters->mMinRadius, sizeof(unsigned));
	fin.read((char*) &mpParameters->mDistance, sizeof(unsigned));
	fin.read((char*) &mpParameters->mMinDistance, sizeof(unsigned));
	fin.read((char*) &mpParameters->mNumHashShingles, sizeof(unsigned));
	fin.read((char*) &mpParameters->mNumRepeatsHashFunction, sizeof(unsigned));
	fin.read((char*) &mpParameters->mSeqWindow, sizeof(unsigned));
	fin.read((char*) &mpParameters->mIndexSeqShift, sizeof(unsigned));
	fin.read((char*) &tmp, sizeof(unsigned));
	SetHistogramSize(tmp);

	mFeature2IndexValue.clear();
	for (unsigned idx=1;idx<=GetHistogramSize();idx++){
		unsigned hist_idx;
		unsigned size;

		fin.read((char*) &hist_idx, sizeof(unsigned));
		fin.read((char*) &size, sizeof(unsigned));
		string feature;
		feature.resize(size);
		fin.read(const_cast<char*>(feature.c_str()), size);
		mFeature2IndexValue.insert(make_pair(feature,hist_idx));
	}

	//unsigned numHashFunc = 0;
	fin.read((char*) &mpParameters->mNumHashFunctions, sizeof(unsigned));
	if (mpParameters->mNumHashFunctions <= 0)
		fin.setstate(std::ios::badbit);
	if (!fin.good())
		return false;

	InitInverseIndex();

	cout << endl << "read "<< mpParameters->mNumHashFunctions << " sub indices ..." << endl;
	for (unsigned  hashFunc = 0; hashFunc < mpParameters->mNumHashFunctions; hashFunc++){

		unsigned numBins = 0;
		fin.read((char*) &numBins, sizeof(unsigned));
		if (numBins < 0)
			fin.setstate(std::ios::badbit);
		if (!fin.good())
			return false;

		index[hashFunc].rehash(numBins);
		cout << "sub index "<< hashFunc+1 << " (keys="<< numBins << ") : "<< flush;
		for (unsigned  bin = 0; bin < numBins; bin++){
			if (bin%(unsigned)(std::ceil((double)numBins/100.0))==0){
				cout << "." << flush;
			}
			unsigned binId = 0;
			unsigned numBinEntries = 0;
			fin.read((char*) &binId, sizeof(unsigned));
			fin.read((char*) &numBinEntries, sizeof(unsigned));
			if (numBinEntries < 0 || binId <= 0)
				fin.setstate(std::ios::badbit);
			if (!fin.good())
				return false;
			//		indexBinTy tmp = indexBinTy(numBinEntries);
			binKeyTy* tmp;
			switch (numBinEntries){
			case 1:
				tmp = reinterpret_cast<binKeyTy(*)>(mMemPool_2[hashFunc]->newElement());
				break;
			case 2:
				tmp = reinterpret_cast<binKeyTy(*)>(mMemPool_3[hashFunc]->newElement());
				break;
			case 3:
				tmp = reinterpret_cast<binKeyTy(*)>(mMemPool_4[hashFunc]->newElement());
				break;
			case 4:
				tmp = reinterpret_cast<binKeyTy(*)>(mMemPool_5[hashFunc]->newElement());
				break;
			case 5:
				tmp = reinterpret_cast<binKeyTy(*)>(mMemPool_6[hashFunc]->newElement());
				break;
			case 6:
				tmp = reinterpret_cast<binKeyTy(*)>(mMemPool_7[hashFunc]->newElement());
				break;
			case 7:
				tmp = reinterpret_cast<binKeyTy(*)>(mMemPool_8[hashFunc]->newElement());
				break;
			case 8:
				tmp = reinterpret_cast<binKeyTy(*)>(mMemPool_9[hashFunc]->newElement());
				break;
			case 9:
				tmp = reinterpret_cast<binKeyTy(*)>(mMemPool_10[hashFunc]->newElement());
				break;
			default:
				tmp = new binKeyTy[numBinEntries+1];
				//		if ((numBinEntries+1)<=100){
				//			tmp = (binKeyTy*)mpool_alloc(mp[hashFunc], (numBinEntries+1) * sizeof(binKeyTy));
				//		} else
				//			tmp = new binKeyTy[numBinEntries+1];
				break;
			}
			//cout << "new bin " << binId << " " << numBinEntries << " ";
			for (unsigned entry = 1; entry <= numBinEntries; entry++ ){

				binKeyTy s;
				fin.read((char*) &s, sizeof(binKeyTy));
				if (s < 0)
					fin.setstate(std::ios::badbit);
				if (!fin.good())
					return false;
				tmp[entry] = s;
				//cout << s << " ";
			}
			//cout << endl;

			tmp[0] = numBinEntries;
			//MyMapV[hashFunc][binId] = tmp;
			index[hashFunc][binId] = tmp;
		}
		cout << endl;
	}
	fin.close();
	return true;
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

//bool compare(const Data::SeqFileT& first, const Data::SeqFileT& second) {
//
//	return (first.uIdx<second.uIdx);
//}

//void  HistogramIndex::PrepareIndexDataSets(vector<SeqDataSet>& myFileList){
//
//	if (!myFileList.size())
//		throw range_error("ERROR no datasets to prepare ...");
//
//	sort(myFileList.begin(),myFileList.end(),compare);
//
//	// check if we have datasets to update the index
//	uint bin = 0;
//	uint userIdx = myFileList[0].uIdx;
//	for (unsigned i=0;i<myFileList.size(); i++){
//		if (myFileList[i].updateIndex){
//			if ( myFileList[i].uIdx != userIdx )
//				bin++;
//			mHistBin2DatasetIdx.insert(make_pair(bin,i));
//			userIdx = myFileList[i].uIdx;
//			myFileList[i].idx = bin;
//			mIndexDataSets.push_back(myFileList[i]);
//			//cout << "bin" << bin << " " << userIdx << " "<<myFileList[i].idx << " " << myFileList[i].uIdx << endl;
//		}
//	}
//	SetHistogramSize(bin+1);
//}

//void NeighborhoodIndex::UpdateInverseIndex(vector<unsigned>& aSignature, unsigned aIndex) {
//	for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k) { //for every hash value
//		unsigned key = aSignature[k];
//		if (key != MAXUNSIGNED && key != 0) { //if key is equal to markers for empty bins then skip insertion instance in data structure
//			if (mInverseIndex[k].count(key) == 0) { //if this is the first time that an instance exhibits that specific value for that hash function, then store for the first time the reference to that instance
//				indexBinTy tmp;
//				tmp.push_back(aIndex);
//				mInverseIndex[k].insert(make_pair(key, tmp));
//				numKeys++; // just for bin statistics
//			} else if (mInverseIndex[k][key].size() < mpParameters->mMaxSizeBin && mInverseIndex[k][key].front() != MAXBINKEY) {
//				// add key to bin if we not have a full bin
//				mInverseIndex[k][key].push_back(aIndex);
//			} else if (mInverseIndex[k][key].size() == mpParameters->mMaxSizeBin){
//				// if a bin is full we clear it and add a key with MAXUNSIGNED to indicate an overfull bin
//				numFullBins++; // just for bin statistics
//				mInverseIndex[k][key].clear();
//				mInverseIndex[k][key].push_back(MAXBINKEY);
//			}
//		}
//	}
//}


//if (myValue[0] == 0){
//				//binKeyTy* myValue = Hash[k].Find(key)->m_data;
////					binKeyTy* myValue = mInverseIndex[k][key];
//				binKeyTy (*foo)[2] = reinterpret_cast<newIndexBin(*)>(mInverseIndex[k][key]);
//				if ((*foo)[(*foo)[0]] != aIndexT){
//					//mInverseIndex[k][key][0] < 2 &&
//					//				binKeyTy*& myValue = mInverseIndex[k][key];
//					//binKeyTy* myValue = Hash[k].Find(key)->m_data;
//
//					// find pos for insert, assume sorted array
//					binKeyTy i = (*foo)[0];
//					while (((*foo)[i]> aIndexT) && (i>1)){
//						i--;
//					}
//
//					// only insert if element is not there
//					if ((*foo)[i]<aIndexT){
//						binKeyTy newSize = ((*foo)[0])+1;
//						binKeyTy * fooNew;
//						fooNew = new binKeyTy[newSize+1];
//
//						memcpy(fooNew,foo,(i+1)*sizeof(binKeyTy));
//						fooNew[i+1] = aIndexT;
//						memcpy(&fooNew[i+2],&(*foo)[i+1],((*foo)[0]-i)*sizeof(binKeyTy));
//						fooNew[0] = newSize;
//
//						mMemPool[k]->deleteElement(foo);
//						mInverseIndex[k][key] = fooNew;
//					}
//				}
//			} else {

//binKeyTy (*foo)[2];
//		foo = mMemPool[k]->newElement();
//		//foo = new binKeyTy[2];
//		(*foo)[1] = (binKeyTy)aIndex;
//		(*foo)[0]= 1; //index of last element is stored at idx[0]

//void HistogramIndex::ComputeHistogram(const vector<unsigned>& aSignature, std::valarray<double>& hist, unsigned& emptyBins) {
//
//	hist.resize(GetHistogramSize());
//	hist *= 0;
//	emptyBins = 0;
//	for (unsigned k = 0; k < aSignature.size(); ++k) {
//
//		if (mInverseIndex[k].count(aSignature[k]) != 0) {
//
//			indexBinTy& myValue = mInverseIndex[k][aSignature[k]];
//
//			std::valarray<double> t(0.0, hist.size());
//
//			for (uint i=1;i<=myValue.size();i++){
//				t[myValue[i]-1]=1;
//			}
//
//			hist += t;
//
//		} else {
//			emptyBins++;
//		}
//	}
//}

//void HistogramIndex::UpdateInverseIndex(const vector<unsigned>& aSignature,const unsigned& aIndex) {
//	const binKeyTy aIndexT =(binKeyTy)aIndex;
//	for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k) { //for every hash value
//		unsigned key = aSignature[k];
//		if (key != MAXUNSIGNED && key != 0) { //if key is equal to markers for empty bins then skip insertion instance in data structure
//			if (mInverseIndex[k].count(key) == 0) { //if this is the first time that an instance exhibits that specific value for that hash function, then store for the first time the reference to that instance
//				indexBinTy tmp;
//				tmp.push_back(aIndex);
//				mInverseIndex[k].insert(make_pair(key, tmp));
//				numKeys++; // just for bin statistics
//			} else if (mInverseIndex[k][key].back() != aIndexT) {
//				// add key to bin if we not have a full bin
//				indexBinTy& myValue = mInverseIndex[k][key];
//				indexBinTy::iterator it = myValue.begin();
//				while (*it<aIndex && it !=myValue.end() ){
//					++it;
//				}
//				mInverseIndex[k][key].insert(it,aIndexT);
//			}
//		}
//	}
//}

//pair<unsigned,unsigned> NeighborhoodIndex::ComputeApproximateSim(const unsigned& aID, const vector<unsigned>& bSignature) {
//
//	//umap_uint_int neighborhood;
//	unsigned collisions = 0;
//	unsigned counts_aID = 0;
//	for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k) {
//		unsigned hash_id = bSignature[k];
//		if (hash_id != 0 && hash_id != MAXUNSIGNED && mInverseIndex[k][hash_id].front() != MAXBINKEY) {
//
//			//fill neighborhood set counting number of occurrences
//			for (indexBinTy::iterator it = mInverseIndex[k][hash_id].begin(); it != mInverseIndex[k][hash_id].end(); ++it) {
//				if (*it == aID) counts_aID++;
//			}
//
//		} else
//			collisions++;
//	}
//	pair<unsigned,unsigned> tmp=make_pair(counts_aID,collisions);
//	return tmp;
//}

//void NeighborhoodIndex::ComputeApproximateNeighborhoodCore(const vector<unsigned>& aSignature, umap_uint_int& neighborhood, unsigned& collisions) {
//
//	collisions = 0;
//	for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k) {
//
//		unsigned hash_id = aSignature[k];
//		if (hash_id != 0 && hash_id != MAXUNSIGNED && mInverseIndex[k][hash_id].front() != MAXBINKEY) {
//
//			//fill neighborhood set counting number of occurrences
//			for (indexBinTy::iterator it = mInverseIndex[k][hash_id].begin(); it != mInverseIndex[k][hash_id].end(); ++it) {
//				binKeyTy instance_id = *it;
//				if (neighborhood.count(instance_id) > 0)
//					neighborhood[instance_id]++;
//				else
//					neighborhood[instance_id] = 1;
//			}
//		} else {
//			collisions++;
//		}
//	}
//}

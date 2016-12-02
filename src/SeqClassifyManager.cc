/*
 * SeqClassifyManager.cc
 *
 *  Created on: 20.03.2015
 *      Author: heyne
 */

#include "SeqClassifyManager.h"
#include "MinHashEncoder.h"

SeqClassifyManager::SeqClassifyManager(Parameters* apParameters, Data* apData):
HistogramIndex(apParameters,apData)
{

}


void SeqClassifyManager::Exec() {

	ProgressBar pb(1000);
	cout << endl << SEP << endl << "INVERSE INDEX"<< endl << SEP << endl;

	CheckParameters();

	SeqFileT mySet;
	mySet.filename            	= mpParameters->mIndexSeqFile;
	mySet.filename_BED		  	= mpParameters->mIndexBedFile;
	mySet.filetype            	= FASTA;
	mySet.checkUniqueSeqNames 	= true;
	mySet.signatureAction	  	= INDEX;
	mySet.groupGraphsBy		  	= SEQ_FEATURE;
	Data::BEDdataP indexBED   	= mpData->LoadBEDfile(mpParameters->mIndexBedFile.c_str());
	mySet.dataBED             	= indexBED;
	mySet.lastMetaIdx			   = 0;
	mySet.strandType			   = FWD;

	mIndexDataSet = std::make_shared<SeqFileT>(mySet);

	string indexName = mpParameters->mIndexBedFile;
	const unsigned pos = mpParameters->mIndexBedFile.find_last_of("/");
	if (std::string::npos != pos)
		indexName = mpParameters->mIndexBedFile.substr(pos+1);

	// create/load new inverse MinHash index against that we can classify other sequences
	if (mpParameters->mNoIndexCacheFile || !std::ifstream(mpParameters->mIndexBedFile+".bhi").good()){
		cout << endl << " *** Creating inverse index *** "<< endl << endl;

		// use desired shift value for index, "LoadData_Threaded" only uses variable mpParameters->mSeqShift
		// assume 0 as not set
		unsigned tmp_shift = mpParameters->mSeqShift;
		if (mpParameters->mIndexSeqShift > 0){
			mpParameters->mSeqShift = mpParameters->mIndexSeqShift;
		} else mpParameters->mIndexSeqShift = mpParameters->mSeqShift;

		SeqFilesT myList;
		myList.push_back(mIndexDataSet);

		wobbleDist = 1;
		InitInverseIndex();
		LoadData_Threaded(myList);

		SetHistogramSize(mIndexDataSet->lastMetaIdx);
		mpParameters->mSeqShift = tmp_shift;

		// write index to file
		if (!mpParameters->mNoIndexCacheFile){
			cout << "inverse index file : " << mpParameters->mIndexBedFile+".bhi" << endl;
			cout << " write index file ... ";
			OutputManager om((indexName + ".bhi").c_str(), mpParameters->mDirectoryPath);
			writeBinaryIndex2(om.mOut,mInverseIndex);
			om.mOut.close();
			mIndexDataSet->filename_index = mpParameters->mDirectoryPath+"/"+indexName+".bhi";
			cout << endl;
		} else {
			cout << "Index is NOT saved to file!"<< endl;
			mIndexDataSet->filename_index = "IN_MEMORY_INDEX_ONLY";
		}

	} else {
		// read existing index file (*.bhi)
		mIndexDataSet->filename_index = mpParameters->mIndexBedFile+".bhi";

		cout << endl << " *** Read inverse index *** "<< endl << endl;
		cout << "inverse index file : " << mpParameters->mIndexBedFile+".bhi" << endl << "read index ...";

		bool indexState = readBinaryIndex2(mpParameters->mIndexBedFile+".bhi",mInverseIndex);

		if (indexState == false)
			throw range_error("\nCannot read index from file " + mpParameters->mIndexBedFile+".bhi\n");

		cout << "finished! ";
		cout << " Index OK! Format Version "<< INDEX_FORMAT_VERSION << endl << endl << "Read index parameters:"<< endl << endl;
		cout << setw(30) << std::right << " hist size  " << GetHistogramSize() << endl;
		cout << setw(30) << std::right << " hash_bit_size  " << mpParameters->mHashBitSize << endl;
		cout << setw(30) << std::right << " random_seed  " << mpParameters->mRandomSeed << endl;
		cout << setw(30) << std::right << " num_hash_functions  " << mpParameters->mNumHashFunctions << endl;
		cout << setw(30) << std::right << " num_repeat_hash_function  " << mpParameters->mNumRepeatsHashFunction << endl;
		cout << setw(30) << std::right << " num_hash_shingles  " << mpParameters->mNumHashShingles << endl;
		cout << setw(30) << std::right << " radius  " << mpParameters->mMinRadius<<".."<<mpParameters->mRadius << endl;
		cout << setw(30) << std::right << " distance  " << mpParameters->mMinDistance<<".."<<mpParameters->mDistance << endl;
		cout << setw(30) << std::right << " seq_window  " << mpParameters->mSeqWindow << endl;
		cout << setw(30) << std::right << " index_seq_shift  " << mpParameters->mIndexSeqShift << " nt" << endl;

		CheckParameters();
	}

	// update IndexValue2Feature map from provided Index BED file
	for (Data::BEDdataIt it=mIndexDataSet->dataBED->begin(); it!=mIndexDataSet->dataBED->end(); ++it ) {
		map<string, uint>::iterator It2 = mFeature2IndexValue.find(it->second->NAME);
		if (It2 != mFeature2IndexValue.end()){
			mIndexValue2Feature.insert(make_pair(It2->second,it->second));
		}
	}

	// just a sanity check if provided index BED matches features present in index
	for (map<string,uint>::iterator it=mFeature2IndexValue.begin();it!=mFeature2IndexValue.end();++it){
		if (mIndexValue2Feature.count(it->second)==0){
			cout << "Provided index BED file " << indexName << " does not contain feature " << it->first << endl;
			cout << "Create dummy feature for it! "<<endl;
			Data::BEDentryP myBED = std::make_shared<Data::BEDentryT>();
			myBED->SEQ = "UNKNOWN_FEAT_"+ it->first;
			myBED->START=0;
			myBED->END=1;
			myBED->NAME=it->first;
			myBED->SCORE=0;
			myBED->STRAND='.';
			myBED->COLS.push_back("DUMMY_BED_ENTRY_FOR_FEATURE_"+it->first);
			myBED->COLS.push_back("DUMMY_BED_ENTRY_FOR_FEATURE_"+it->first);
			mIndexValue2Feature.insert(make_pair(it->second,myBED));
		}
	}

	// do the classification
	ClassifySeqs();

	pb.PrintElapsed();
}

void SeqClassifyManager::worker_Classify(int numWorkers, unsigned id){

	Signature* tmpSig = new Signature(numHashFunctionsFull);
	while (!done){

		ChunkP myData;
		vector<ChunkP> myQ;

		while (graph_queue[id].size()>=1){
			graph_queue[id].wait_and_pop(myData);
			myQ.push_back(myData);
		}
		//if (myData != NULL) {
		while (!done && myQ.size()){
			myData = myQ.back();
			myQ.pop_back();

			//cout << "  graph2sig thread got chunk " << myData->size() << " offset " << (*myData)[0].idx << " " << mpParameters->mHashBitSize << endl;

			ResultChunkP myResultChunk = std::make_shared<ResultChunkT>();

			//			for (unsigned j = 0; j < myData->size(); j++) {
			//
			//				generate_feature_vector((*myData)[j].seq, (*myData)[j].svec);
			//				ComputeHashSignature((*myData)[j].svec,(*myData)[j].sig,tmpSig);
			//
			//			}

			for (ChunkT::iterator j=myData->begin(); j!=myData->end();j++){
				sliding_window_minhash(j->minHashes,j->seq,mpParameters->mMinRadius,mpParameters->mRadius, mpParameters->mMinDistance, mpParameters->mDistance, mpParameters->mSeqWindow, mpParameters->mSeqShift);
				//cout << "final " << j->name << " len=" << j->seq.size() << "idx=" << j->idx << " minH_hf " << j->minHashes.size() <<  " minh_len "<< j->minHashes[0].size() << " win=" << mpParameters->mSeqWindow << " step=" << mpParameters->mSeqShift << endl;
			}

			finishUpdate(myData,myResultChunk);
			res_queue.push(myResultChunk);
			if (res_queue.size()>=numWorkers*25){
				unique_lock<mutex> lk(mut2);
				cv2.wait(lk,[&]{if ((done) || (res_queue.size()<=numWorkers*10)) return true; else return false;});
				lk.unlock();
			}

		}
	}
	delete tmpSig;
}

void SeqClassifyManager::finisher_Results(ogzstream* fout_res){
	ProgressBar progress_bar(1000);
	unsigned sigCounter = 0;

	while (!done){

		ResultChunkP myResults;
		bool succ = res_queue.try_pop(myResults);

		if (!done && succ && myResults->size()>0) {

			for (unsigned i=0; i<myResults->size(); i++){
				*fout_res << (*myResults)[i].output_line;
				//mResultCounter += (*myResults)[i].numInstances;
				sigCounter += (*myResults)[i].numInstances;
			}
			mResultCounter += myResults->size()*2;
			double elap = progress_bar.getElapsed()/1000;
			cout.setf(ios::fixed);
			cout << "\r" <<  std::setprecision(1) << elap << " sec elapsed   Finished numSeqs=" << std::setprecision(0) << setw(10);
			cout << mNumSequences  << "("<<mNumSequences/(elap) <<" seq/s)  signatures=" << setw(10);
			cout << sigCounter << "("<<(double)sigCounter/((elap)) <<" sig/s - "<<(double)sigCounter/elap/mpParameters->mNumThreads <<" per thread)  inst=";
			cout << mInstanceCounter << " resQueue=" << res_queue.size() << " graphQueue=";
			uint avg = 0;
			for (uint i=0; i<graph_queue.size(); ++i){
				avg += graph_queue[i].size();
			};
			cout << avg/graph_queue.size() << "   ";
		}
	}
	progress_bar.PrintElapsed();
}

void SeqClassifyManager::Classify_Signatures(SeqFilesT& myFiles){

	for (unsigned i=0;i<myFiles.size(); i++){
		readFile_queue.push(myFiles[i]);
	}

	HashSignatureHelper();

	cout << "hist size: " << GetHistogramSize() << endl;
	cout << "Using " << mHashBitMask << " as bitmask"<< endl;
	cout << "Using " << mpParameters->mHashBitSize << " bits to encode features" << endl;
	cout << "Using " << mpParameters->mRandomSeed << " as random seed" << endl;
	cout << "Using " << mpParameters->mNumHashFunctions << " final hash VALUES (param 'num_hash_functions' is used as FINAL signature size!), " << endl;
	cout << "Using " << mpParameters->mNumHashFunctions*mpParameters->mNumHashShingles << " hash values internally (num_hash_functions * num_hash_shingles) from " << mpParameters->mNumRepeatsHashFunction << " independent hash functions (param 'num_repeat_hash_function')" << endl;
	cout << "Using " << mpParameters->mNumHashShingles << " as hash shingle factor" << endl;
	cout << "Using " << mpParameters->mPureApproximateSim << " as minimal approximate similarity (pure_approximate_sim)" << endl;
	cout << "Using feature radius   " << mpParameters->mMinRadius<<".."<<mpParameters->mRadius << endl;
	cout << "Using feature distance " << mpParameters->mMinDistance<<".."<<mpParameters->mDistance << endl;
	cout << "Using sequence window  " << mpParameters->mSeqWindow<<" shift "<<mpParameters->mSeqShift << " nt - clip " << mpParameters->mSeqClip << endl;
	cout << endl << "Computing MinHash signatures on the fly while reading " << myFiles.size() << " file(s)..." << endl;

	int graphWorkers = std::thread::hardware_concurrency();
	if (mpParameters->mNumThreads>0)
		graphWorkers = mpParameters->mNumThreads;

	cout << "Using " << graphWorkers << " worker threads and 2 helper threads..." << endl;

	done					= false;
	files_done			= 0;
	mSignatureCounter = 0;
	mInstanceCounter 	= 0;
	mSequenceCounter  = 0;
	mResultCounter 	= 0;

	vector<std::thread> threads;
	graph_queue.resize(graphWorkers);

	// launch all threads
	threads.push_back( std::thread(&SeqClassifyManager::finisher_Results,this,myFiles[0]->out_results_fh));
	for (int i=0;i<graphWorkers;i++){
		threads.push_back( std::thread(&SeqClassifyManager::worker_Classify,this,graphWorkers,i));
	}
	threads.push_back( std::thread(&MinHashEncoder::worker_readFiles,this,graphWorkers,2));

	{
		join_threads joiner(threads);

		while(!done){
			std::this_thread::sleep_for(std::chrono::milliseconds(20));
			if ( (files_done>=myFiles.size()) && (mResultCounter >= mInstanceCounter))
				done = true;
			else done = false;
			cv2.notify_all();
			cv1.notify_all();
		}

	} // by leaving this block threads get joined by destruction of joiner

	if (mInstanceCounter == 0) {
		throw range_error("ERROR in MinHashEncoder::LoadData: something went wrong; no instances/signatures produced");
	} else
		cout << "Instances/signatures produced " << mInstanceCounter << " " << mResultCounter << endl;

	//cout << endl << "Classification finished - signatures=" << mSignatureCounter << " instances=" << mNumSequences << " classified=" << mClassifiedInstances<< endl;

	cout << SEP << endl;
	cout << " CLASSIFICATION FINISHED" << endl << SEP << endl;
}


void SeqClassifyManager::finishUpdate(ChunkP& myData, unsigned& min, unsigned& max) {

	// this is the overloaded virtual function from MinHashEncoder
	// we assume signatureAction==INDEX always here
	if (mUseSlidingWindowMinHash){
		for (ChunkT::iterator j=myData->begin(); j!= myData->end();j++) {
			for (uint hf=min;hf<=max;hf++){
				unsigned last = MAXUNSIGNED;
				for (auto i : j->minHashes[hf] ){
					if (i != last)
						UpdateInverseIndex(i, j->idx, hf);
					last = i;
					if (hf==0) {
						mSignatureUpdateCounter++;
					}
				}
			}
		}
	} else {

		for (ChunkT::iterator j=myData->begin(); j!= myData->end();j++) {
			UpdateInverseIndex(j->sig,j->idx,min,max);
			if (min==0){
				mSignatureUpdateCounter++;
			}
		}
	}
}

inline double changeNAN(double i) {if (std::isnan(i)) return 0.0; else return i;}

inline double indicator(double i) {if (i>0) return 1; else return 0;}


void SeqClassifyManager::finishUpdate(ChunkP& myData, ResultChunkP& myResultChunk) {

	ChunkT::iterator j = myData->begin();
	unsigned hist_size = GetHistogramSize();

	while (j != myData->end()) {
		valarray<double> hist(0.0,hist_size);
		valarray<double> histRC(0.0,hist_size);

		unsigned emptyBins = 0;
		unsigned emptyBinsRC = 0;

		ChunkT::iterator k=j;
		unsigned matchingSigs 	= 0;
		unsigned matchingSigsRC = 0;
		unsigned numSigs = 0;
		unsigned numSigsRC = 0;

		do {
			valarray<double> hist_tmp;
			vector<unsigned> emptyBins_tmp;
			// compute hist for all sliding windows at once
			ComputeHistogram(k->minHashes,hist_tmp,emptyBins_tmp);

			switch (k->rc){
			case true:
				histRC += hist_tmp;
				for (uint i = 0; i < emptyBins_tmp.size(); ++i){
					if ( emptyBins_tmp[i] < mpParameters->mNumHashFunctions) matchingSigsRC++;
					emptyBinsRC += emptyBins_tmp[i];
				}
				numSigsRC += emptyBins_tmp.size();
				break;
			case false:
				hist += hist_tmp;
				for (uint i = 0; i< emptyBins_tmp.size(); ++i){
					if ( emptyBins_tmp[i] < mpParameters->mNumHashFunctions) matchingSigs++;
					emptyBins += emptyBins_tmp[i];
				}
				numSigs += emptyBins_tmp.size();
				break;
			}
			k++;
		} while (k != myData->end() && k->idx == j->idx);

		ResultT myResult;
		unsigned max = hist.max();
		unsigned maxRC = histRC.max();

		switch (j->seqFile->strandType){
		case FWD:
			myResult.numInstances = numSigs;
			getResultString(myResult.output_line,hist,emptyBins,matchingSigs,numSigs,j->name,FWD);
			myResultChunk->push_back(myResult);
			break;
		case REV:
			myResult.numInstances = numSigsRC;
			getResultString(myResult.output_line,histRC,emptyBinsRC,matchingSigsRC,numSigsRC,j->name,REV);
			myResultChunk->push_back(myResult);
			break;
		case FR:
			switch (mpParameters->mOutputTypeCode){
			case ALL:
			case MAX:
				myResult.numInstances = numSigs+numSigsRC;
				getResultString(myResult.output_line,hist+histRC,emptyBins+emptyBinsRC,matchingSigs+matchingSigsRC,numSigs+numSigsRC,j->name,FR);
				myResultChunk->push_back(myResult);
				break;
			case ALL_STRAND:
			case MAX_STRAND:

				if (max > maxRC){
					myResult.numInstances = numSigs;
					getResultString(myResult.output_line,hist,emptyBins,matchingSigs,numSigs,j->name,FWD);
					myResultChunk->push_back(myResult);
				} else if (maxRC>max){
					myResult.numInstances = numSigsRC;
					getResultString(myResult.output_line,histRC,emptyBinsRC,matchingSigsRC,numSigsRC,j->name,REV);
					myResultChunk->push_back(myResult);
				} else {
					myResult.numInstances = numSigs+numSigsRC;
					getResultString(myResult.output_line,hist+histRC,emptyBins+emptyBinsRC,matchingSigs+matchingSigsRC,numSigs+numSigsRC,j->name,FR);
					myResultChunk->push_back(myResult);
				}
				break;
			default:
				break;
			}
			break;
		case FR_sep:
				myResult.numInstances = numSigs;
				getResultString(myResult.output_line,hist,emptyBins,matchingSigs,numSigs,j->name,FWD);
				myResultChunk->push_back(myResult);
				myResult.numInstances = numSigsRC;
				getResultString(myResult.output_line,histRC,emptyBinsRC,matchingSigsRC,numSigsRC,j->name,REV);
				myResultChunk->push_back(myResult);
				break;
		default:
			break;
		}

		j = k;
		++mNumSequences;

		// meta analysis, only for screen output summary
		if (mpParameters->mVerbose){

			hist += histRC;
			unsigned sum = hist.sum();
			unsigned max = hist.max();

			if (sum>0)
				mClassifiedInstances++;

			valarray<double> hist_t = hist;
			hist_t /= ( (numSigs+numSigsRC)*mpParameters->mNumHashFunctions);

			hist_t /= sum;
			hist_t = hist.apply(changeNAN);

			for (unsigned i = 0; i<hist_size; i++){
				if (hist[i] >= max  && max != 0) {hist[i] = 1;} else { hist[i]=0;};
			}

			{
				std::lock_guard<std::mutex> lk(mut_meta);
				metaHist += hist_t;
				metaHistNum += hist;
			}
		}
	}

}


/*void SeqClassifyManager::finishUpdate(ChunkP& myData, ResultChunkP& myResultChunk) {

	unsigned j = 0;
	unsigned hist_size = GetHistogramSize();

	while (j < (*myData).size()) {

		switch ((*myData)[j].seqFile->signatureAction){
		case CLASSIFY: {

			valarray<double> hist(0.0,hist_size);
			valarray<double> histRC(0.0,hist_size);

			unsigned emptyBins = 0;
			unsigned emptyBinsRC = 0;

			unsigned k=0;
			unsigned matchingSigs 	= 0;
			unsigned matchingSigsRC = 0;

			do {
				valarray<double> hist_tmp;
				unsigned emptyBins_tmp;
				ComputeHistogram((*myData)[j+k].sig,hist_tmp,emptyBins_tmp);

				switch ((*myData)[j+k].rc){
				case true:
					histRC += hist_tmp;
					if ( emptyBins_tmp < hist_size ) matchingSigsRC++;
					emptyBinsRC += emptyBins_tmp;
					break;
				case false:
					hist += hist_tmp;
					if ( emptyBins_tmp < hist_size ) matchingSigs++;
					emptyBins += emptyBins_tmp;
					break;
				}
				k++;
			} while (j+k<myData->size() && (*myData)[j].name==(*myData)[j+k].name);

			ResultT myResult;

			switch ((*myData)[j].seqFile->strandType){
			case FWD:
				myResult.numInstances = k;
				getResultString(myResult.output_line,hist,emptyBins,matchingSigs,k,(*myData)[j].name,FWD);
				myResultChunk->push_back(myResult);
				break;
			case REV:
				myResult.numInstances = k;
				getResultString(myResult.output_line,histRC,emptyBinsRC,matchingSigsRC,k,(*myData)[j].name,REV);
				myResultChunk->push_back(myResult);
				break;
			case FR:
				myResult.numInstances = k;
				getResultString(myResult.output_line,hist+histRC,emptyBins+emptyBinsRC,matchingSigs+matchingSigsRC,k,(*myData)[j].name,FR);
				myResultChunk->push_back(myResult);
				break;
			case FR_sep:
				myResult.numInstances = k;
				getResultString(myResult.output_line,hist,emptyBins,matchingSigs,k/2,(*myData)[j].name,FWD);
				myResultChunk->push_back(myResult);
				myResult.numInstances = k;
				getResultString(myResult.output_line,histRC,emptyBinsRC,matchingSigsRC,k/2,(*myData)[j].name,REV);
				myResultChunk->push_back(myResult);
				break;
			default:
				break;
			}

			j += k;
			++mNumSequences;

			// meta analysis, only for screen output summary
			if (mpParameters->mVerbose){

				hist += histRC;
				unsigned sum = hist.sum();
				unsigned max = hist.max();

				if (sum>0)
					mClassifiedInstances++;

				valarray<double> hist_t = hist;
				hist_t /= (k*mpParameters->mNumHashFunctions);

				hist_t /= sum;
				hist_t = hist.apply(changeNAN);

				for (unsigned i = 0; i<hist_size; i++){
					if (hist[i] >= max  && max != 0) {hist[i] = 1;} else { hist[i]=0;};
				}

				{
					std::lock_guard<std::mutex> lk(mut_meta);
					metaHist += hist_t;
					metaHistNum += hist;
				}
			}
			break;
		}
		default:
			break;
		}
	} // while
}*/


void SeqClassifyManager::getResultString(string& resT,histogramT hist,unsigned emptyBins, unsigned matchingSigs, unsigned numSigs, string& name, strandTypeT strand){

	stringstream res;

	uint sum = hist.sum();
	uint max = hist.max();

	if (mpParameters->mPureApproximateSim != 0) {
		for (uint i = 0; i<hist.size(); i++){
			if (hist[i] / (numSigs*mpParameters->mNumHashFunctions) < mpParameters->mPureApproximateSim ) {
				hist[i] = 0.0;
			};
		}
	}

	string str;
	switch (strand) {
	case FWD:
		str="+";
		break;
	case REV:
		str="-";
		break;
	case FR:
		str=".";
		break;
	default:
		break;
	}

	res << name << "\t" << str << "\t"<< numSigs << "\t" << matchingSigs <<  "\t" << (numSigs*mpParameters->mNumHashFunctions)-emptyBins << "\t" << sum << "\t" << max << "\t";

	string maxIndices;
	string indices;
	string values;

	switch (mpParameters->mOutputTypeCode){
	case ALL_STRAND:
	case ALL:
		for (unsigned i=0; i<hist.size();i++){
			if (hist[i]>0.0) {
				char buf[30];
				sprintf( buf,"%.0f",hist[i] );
				values += std::string(buf) +",";
				//values += std::to_string((int)hist[i])+",";
				//sprintf( buf,"%i",i+1 );
				indices += std::to_string(i+1)+ ",";
				//indices += std::string(buf) +",";
				if (hist[i]==max && max!=0){
					//maxIndices += std::string(buf) +",";
					maxIndices+=std::to_string(i+1)+ ","; // << i+1 << ",";
				}
			}
		}
		break;
	case MAX:
	case MAX_STRAND:
		for (unsigned i=0; i<hist.size();i++){
			if (hist[i]==max && max!=0){
				char buf[30];
				sprintf( buf,"%.0f",hist[i] );
				values += std::string(buf) +",";
				indices += std::to_string(i+1)+ ",";
				maxIndices+=std::to_string(i+1)+ ",";
			}
		}
		break;
	default:
		break;
	}

	if (max!=0) {
		res << indices << "\t" <<values << "\t" << maxIndices;
	}

	res << "\n";
	resT = res.str();
}

void SeqClassifyManager::ClassifySeqs(){

	cout << endl << SEP << endl << "CLASSIFICATION" << endl << SEP << endl;

	cout << endl << " *** Read sequences for classification and create their MinHash signatures *** " << endl << endl;

	// prepare sequence set for classification
	SeqFileP mySet = std::make_shared<SeqFileT>();
	mySet->filename            = mpParameters->mInputDataFileName;
	mySet->filetype            = mpParameters->mFileTypeCode;
	mySet->groupGraphsBy       = SEQ_NUM; // actually we check by InstanceT.name field for graphs from one seq
	mySet->checkUniqueSeqNames = false;
	mySet->signatureAction	   = CLASSIFY;
	mySet->strandType          = FR;

	// write results file header and get results file handle
	string resultsName = mpParameters->mInputDataFileName;
	const unsigned pos = mpParameters->mInputDataFileName.find_last_of("/");
	if (std::string::npos != pos)
		resultsName = mpParameters->mInputDataFileName.substr(pos+1);

	mySet->out_results_fh = PrepareResultsFile(mpParameters->mDirectoryPath+resultsName+".classified.tab.gz");

	metaHist.resize(GetHistogramSize());
	metaHist *= 0;

	metaHistNum.resize(GetHistogramSize());
	metaHistNum *= 0;

	mClassifiedInstances = 0;
	mNumSequences = 0;

	wobbleDist = 1;

	SeqFilesT myList;
	myList.push_back(mySet);

	//do the real work
	Classify_Signatures(myList);

	mySet->out_results_fh->close();

	/////////////////////////////////////////////////////////////////////////////
	// classification finished
	/////////////////////////////////////////////////////////////////////////////

	// metahistogram

	if (mpParameters->mVerbose){

		cout << endl << endl << "META histogram - classified seqs: " << setprecision(3) << (double)mClassifiedInstances/((double)mNumSequences) << " (" << mClassifiedInstances << ") - TOP 20" << endl;

		metaHist = metaHist/metaHist.sum();
		vector<pair<double,uint> > sortedHist;
		for (unsigned j=0; j<metaHist.size();j++){
			sortedHist.push_back(make_pair(-metaHistNum[j],j));
		}
		sort(sortedHist.begin(), sortedHist.end());
		for (unsigned j=0; j<std::min((unsigned)20,(unsigned)sortedHist.size());j++){
			multimap<uint, Data::BEDentryP>::iterator it = mIndexValue2Feature.find(sortedHist[j].second+1);
			uint num = mIndexValue2Feature.count(sortedHist[j].second+1);
			cout << setprecision(2) << j+1 << "\t" << metaHist[sortedHist[j].second] << "\t" << (uint)metaHistNum[sortedHist[j].second] << "\t" << sortedHist[j].second+1 << "\t";
			if (it != mIndexValue2Feature.end()) cout << "feature\t"<< it->second->NAME << "\t#features=" << num;
			if (it != mIndexValue2Feature.end() && it->second->COLS.size()>=2) cout << "\t" << it->second->COLS[1];
			cout << endl;
		}
		cout << "SUM\t"<< metaHist.sum() << endl << endl;

		// bin size statistics

		valarray<double> indexHist;
		indexHist.resize(GetHistogramSize());
		indexHist *= 0;

		for (typename HistogramIndex::indexTy::const_iterator it = mInverseIndex.begin(); it!= mInverseIndex.end(); it++){
			for (typename HistogramIndex::indexSingleTy::const_iterator itBin = it->begin(); itBin!=it->end(); itBin++){
				indexHist[itBin->second[0]-1] += 1;
			}
		}


		const char * indexBinSizePath;
		if (mpParameters->mDirectoryPath != "") indexBinSizePath = (mpParameters->mDirectoryPath + "/IndexBinSize.tab").c_str();
		else indexBinSizePath = (mpParameters->mDirectoryPath + "IndexBinSize.tab").c_str();

		cout << indexBinSizePath << endl;
		FILE * indexBinSizeFile = fopen(indexBinSizePath, "w"); 	// opening index bin sizes file


		for (unsigned i=0; i<indexHist.size();i++){
			string outputString = " index bin size " + to_string(i+1) + "\t" + to_string(indexHist[i]);
			cout << outputString << "\t" << setprecision(5) << indexHist[i]/indexHist.sum() << endl; 	// writing to stdout
			const char * outputLine = (outputString + "\n").c_str();
			if (indexBinSizeFile != nullptr) {
				fputs(outputLine, indexBinSizeFile); 	// writing to file (if opening was successful)
			}
		}
		fclose(indexBinSizeFile); 	// closing file again
	}
}

ogzstream* SeqClassifyManager::PrepareResultsFile(string filename){

	ogzstream* fout = new ogzstream(filename.c_str(),std::ios::out);
	// write header to output results file
	// parameters
	*fout << "##INDEX PARAMETERS" << endl;
	*fout << "##" << endl;
	*fout << "#PARAM\tINDEX\t" << mIndexDataSet->filename_index << endl;
	*fout << "#PARAM\tSEQFILE\t" << mIndexDataSet->filename << endl;
	*fout << "#PARAM\tBEDFILE\t" << mpParameters->mIndexBedFile << endl;
	*fout << "#PARAM\tHASHBITSIZE\t" << mpParameters->mHashBitSize << endl;
	*fout << "#PARAM\tRANDOMSEED\t" << mpParameters->mRandomSeed << endl;
	*fout << "#PARAM\tNUMHASHFUNC\t" << mpParameters->mNumHashFunctions << endl;
	*fout << "#PARAM\tNUMREPEATHASHFUNC\t" << mpParameters->mNumRepeatsHashFunction << endl;
	*fout << "#PARAM\tNUMHASHSHINGLES\t" << mpParameters->mNumHashShingles << endl;
	*fout << "#PARAM\tMINRADIUS\t" << mpParameters->mMinRadius << endl;
	*fout << "#PARAM\tRADIUS\t" << mpParameters->mRadius << endl;
	*fout << "#PARAM\tMINDISTANCE\t" << mpParameters->mMinDistance<< endl;
	*fout << "#PARAM\tDISTANCE\t" << mpParameters->mDistance << endl;
	*fout << "#PARAM\tSEQWINDOW\t" << mpParameters->mSeqWindow<< endl;
	*fout << "#PARAM\tINDEXSEQSHIFT\t" << mpParameters->mIndexSeqShift << endl;
	*fout << "#PARAM\tHISTOGRAMSIZE\t" << GetHistogramSize() << endl;
	*fout << "##" << endl;
	*fout << "##CLASSIFY PARAMETERS" << endl;
	*fout << "##" << endl;
	*fout << "#PARAM\tINPUTFILE\t" << mpParameters->mInputDataFileName << endl;
	*fout << "#PARAM\tSEQSHIFT\t" <<mpParameters->mSeqShift << endl;
	*fout << "#PARAM\tSEQCLIP\t" <<mpParameters->mSeqClip << endl;
	*fout << "#PARAM\tAPPROXSIM\t" <<mpParameters->mPureApproximateSim << endl;
	*fout << "##" << endl;
	*fout << "##INDEX MAPPING TABLE" << endl;
	*fout << "##" << endl;

	// mapping table histogram idx -> feature
	vector<pair<uint,string> > items;
	for (std::map<string,uint>::iterator it = mFeature2IndexValue.begin(); it != mFeature2IndexValue.end();++it) {
		items.push_back(make_pair(it->second,it->first));
	}

	sort(items.begin(),items.end());

	for (vector<pair<uint,string> >::iterator item=items.begin(); item != items.end();++item){
		*fout << "#HIST_IDX\t"<< item->first << "\t" << "feature\t"<< item->second;
		multimap<uint,Data::BEDentryP>::iterator it2 = mIndexValue2Feature.find(item->first);
		if (it2 != mIndexValue2Feature.end() && it2->second->COLS.size()>=2) *fout << "\t" << it2->second->COLS[1];
		*fout << endl;
	}

	*fout << "##" << endl;
	*fout << "##SEQUENCE CLASSIFICATION RESULTS" << endl;
	*fout << "##" << endl;
	*fout << "#SEQ\tSTR\tSIGS\tSIG_HITS\tHF_HITS\tSUM\tMAX\tIDX\tVALS\tMAX_IDX"<< endl;
	return fout;
}

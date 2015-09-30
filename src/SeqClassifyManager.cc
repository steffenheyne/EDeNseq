/*
 * SeqClassifyManager.cc
 *
 *  Created on: 20.03.2015
 *      Author: heyne
 */

#include "SeqClassifyManager.h"
#include "MinHashEncoder.h"


SeqClassifyManager::SeqClassifyManager(Parameters* apParameters, Data* apData):
HistogramIndex(apParameters,apData),pb(1000)
{

}


void SeqClassifyManager::Exec() {

	// load and prepare index list file to get files for indexing
	//	vector<SeqDataSet> fileList = mpData->LoadIndexDataList(mpParameters->mIndexDataList.c_str());

	SeqFileT mySet;
	mySet.filename            	= mpParameters->mIndexSeqFile;
	mySet.filename_BED		  	= mpParameters->mIndexBedFile;
	mySet.filetype            	= FASTA;
	mySet.checkUniqueSeqNames 	= true;
	mySet.signatureAction	  	= INDEX;
	mySet.groupGraphsBy		  	= SEQ_FEATURE;
	Data::BEDdataP indexBED   	= mpData->LoadBEDfile(mpParameters->mIndexBedFile.c_str());
	mySet.dataBED             	= indexBED;
	mySet.lastMetaIdx				= 0;
	mySet.strandType				= FWD;

	mIndexDataSet = std::make_shared<SeqFileT>(mySet);

	string indexName = mpParameters->mIndexBedFile;
	const unsigned pos = mpParameters->mIndexBedFile.find_last_of("/");
	if (std::string::npos != pos)
		indexName = mpParameters->mIndexBedFile.substr(pos+1);

	// error if the shift is 0 and we have a specific window size
	if (mpParameters->mSeqWindow != 0 && mpParameters->mSeqShift==0)
		throw range_error("\nERROR! 'seq_shift' cannot be 0 for a specific window size!");

	// create/load new inverse MinHash index against that we can classify other sequences
	if (!std::ifstream(mpParameters->mIndexBedFile+".bhi").good()){
		cout << endl << " *** Creating inverse index *** "<< endl << endl;

		// use desired shift value for index, "LoadData_Threaded" only uses variable mpParameters->mSeqShift
		// assume 0 as not set
		double tmp_shift = mpParameters->mSeqShift;
		if (mpParameters->mIndexSeqShift > 0){
			mpParameters->mSeqShift = mpParameters->mIndexSeqShift;
		} else mpParameters->mIndexSeqShift = mpParameters->mSeqShift;

		SeqFilesT myList;
		myList.push_back(mIndexDataSet);
		LoadData_Threaded(myList);
		SetHistogramSize(mIndexDataSet->lastMetaIdx);
		mpParameters->mSeqShift = tmp_shift;

		// write index to file
		if (!mpParameters->mNoIndexCacheFile){
			cout << "inverse index file : " << mpParameters->mIndexBedFile+".bhi" << endl;
			cout << " write index file ... ";
			OutputManager om((indexName+".bhi").c_str(), mpParameters->mDirectoryPath);
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
		cout << endl << " *** Read inverse index *** "<< endl << endl;
		cout << "inverse index file : " << mpParameters->mIndexBedFile+".bhi" << endl << "read index ..." << endl;
		bool indexState = readBinaryIndex2(mpParameters->mIndexBedFile+".bhi",mInverseIndex);
		cout << "finished! ";
		if (indexState == true){
			cout << " Index OK!" << endl;

		} else
			throw range_error("Cannot read index from file " + mpParameters->mIndexBedFile+".bhi");
		mIndexDataSet->filename_index = mpParameters->mIndexBedFile+".bhi";
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

	ClassifySeqs();
}

void SeqClassifyManager::worker_Classify(int numWorkers){

	while (!done){

		ChunkP myData;
		unique_lock<mutex> lk(mut2);
		cv2.wait(lk,[&]{if ( (done) ||  (graph_queue.try_pop( (myData) )) ) return true; else return false;});
		lk.unlock();

		if (!done && myData->size()>0) {
			//cout << "  graph2sig thread got chunk " << myData->size() << " offset " << (*myData)[0].idx << " " << mpParameters->mHashBitSize << endl;

			ResultChunkP myResultChunk = std::make_shared<ResultChunkT>();

			for (unsigned j = 0; j < myData->size(); j++) {
				//SVector x(pow(2, mpParameters->mHashBitSize));
				//(*myData)[j].svec.resize(pow(2, mpParameters->mHashBitSize));

				generate_feature_vector((*myData)[j].seq, (*myData)[j].svec);

				(*myData)[j].sig = MinHashEncoder::ComputeHashSignature((*myData)[j].svec);

			}
			finishUpdate(myData,myResultChunk);
			//			if (res_queue.size()>=numWorkers*25){
			//				unique_lock<mutex> lk(mut_res);
			//				cv_res.wait(lk,[&]{if ((done) || (res_queue.size()<=numWorkers*10)) return true; else return false;});
			//				lk.unlock();
			//			}
			res_queue.push(myResultChunk);
		}
		cv2.notify_all();
		cv_res.notify_all();
	}
}

void SeqClassifyManager::finisher_Results(ogzstream* fout_res){
	ProgressBar progress_bar(1000);
	while (!done){

		ResultChunkP myResults;
		unique_lock<mutex> lk(mut_res);
		cv_res.wait(lk,[&]{if ( (done) || (res_queue.try_pop( (myResults) ))) return true; else return false;});
		lk.unlock();

		if (!done && myResults->size()>0) {

			uint chunkSize = myResults->size();

			for (unsigned i=0; i<myResults->size(); i++){
				*fout_res << (*myResults)[i].output_line;
				mResultCounter += (*myResults)[i].numInstances;
			}

			if (mResultCounter%1000000 <= 10*chunkSize){
				cout << endl << "    Resultfinisher updated index with " << chunkSize << " signatures all_sigs=" <<  mSignatureCounter << " inst=" << mInstanceCounter << " resQueue=" << res_queue.size() << endl;
			}
			progress_bar.Count(mResultCounter);
		}
		cv1.notify_all();
		cv2.notify_all();
		cvm.notify_all();
		cv_res.notify_all();
	}
}

void SeqClassifyManager::Classify_Signatures(SeqFilesT& myFiles){

	for (unsigned i=0;i<myFiles.size(); i++){
		readFile_queue.push(myFiles[i]);
	}
	cout << "Using " << mpParameters->mHashBitSize << " bits to encode features" << endl;
	cout << "Using " << mpParameters->mRandomSeed << " as random seed" << endl;
	cout << "Using " << mpParameters->mNumHashFunctions << " hash functions (with factor " << mpParameters->mNumRepeatsHashFunction << " for single minhash)" << endl;
	cout << "Using " << mpParameters->mNumHashShingles << " as hash shingle factor" << endl;
	cout << "Using " << mpParameters->mPureApproximateSim << " as approximate similarity " << endl;
	cout << "Using feature radius   " << mpParameters->mMinRadius<<".."<<mpParameters->mRadius << endl;
	cout << "Using feature distance " << mpParameters->mMinDistance<<".."<<mpParameters->mDistance << endl;
	cout << "Using sequence window  " << mpParameters->mSeqWindow<<" shift "<<mpParameters->mSeqShift << " ("<< (unsigned)std::max((double)1,(double)mpParameters->mSeqWindow*mpParameters->mSeqShift) << ") clip " << mpParameters->mSeqClip << endl;

	cout << "Computing MinHash signatures on the fly while reading " << myFiles.size() << " file(s)..." << endl;

	int graphWorkers = std::thread::hardware_concurrency();
	if (mpParameters->mNumThreads>0)
		graphWorkers = mpParameters->mNumThreads;

	cout << "Using " << graphWorkers << " worker threads and 2 helper threads..." << endl;

	done = false;
	files_done=0;
	mSignatureCounter = 0;
	mInstanceCounter = 0;
	mResultCounter = 0;
	vector<std::thread> threads;

	threads.push_back( std::thread(&SeqClassifyManager::finisher_Results,this,myFiles[0]->out_results_fh));
	for (int i=0;i<graphWorkers;i++){
		threads.push_back( std::thread(&SeqClassifyManager::worker_Classify,this,graphWorkers));
	}
	threads.push_back( std::thread(&SeqClassifyManager::worker_readFiles,this,graphWorkers));

	{
		join_threads joiner(threads);

		unique_lock<mutex> lk(mutm);
		cv1.notify_all();
		while(!done){
			cvm.wait(lk,[&]{if ( (files_done<myFiles.size()) || (mInstanceCounter > mResultCounter) ) return false; else return true;});
			lk.unlock();
			done=true;
			cv3.notify_all();
			cv2.notify_all();
			cv1.notify_all();
		}

	} // by leaving this block threads get joined by destruction of joiner

	cout << " sig classifier finished" << endl;

	if (mInstanceCounter == 0) {
		throw range_error("ERROR in MinHashEncoder::LoadData: something went wrong; no instances/signatures produced");
	} else
		cout << "Instances/signatures produced " << mInstanceCounter << " " << mResultCounter << endl;
}



inline double changeNAN(double i) {if (std::isnan(i)) return 0.0; else return i;}

inline double indicator(double i) {if (i>0) return 1; else return 0;}

//inline double minSim(double i) {if (i<0.1) return 0; else return i;}

void SeqClassifyManager::finishUpdate(ChunkP& myData) {

	// as this is the overloaded virtual function from MinHshEncoder
	//we assume signatureAction==INDEX always here
	for (unsigned j = 0; j < myData->size(); j++) {
		UpdateInverseIndex((*myData)[j].sig, (*myData)[j].idx);
	}
}

void SeqClassifyManager::finishUpdate(ChunkP& myData, ResultChunkP& myResultChunk) {

	unsigned j = 0;
	while (j < (*myData).size()) {

		switch ((*myData)[j].seqFile->signatureAction){
		case INDEX:{
			UpdateInverseIndex((*myData)[j].sig, (*myData)[j].idx);
			j++;
			break;
		}
		case CLASSIFY: {

			valarray<double> hist(0.0,GetHistogramSize());
			valarray<double> histRC(0.0,GetHistogramSize());

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
					if ( hist_tmp.sum() != 0 ) matchingSigsRC++;
					emptyBinsRC += emptyBins_tmp;
					break;
				case false:
					hist += hist_tmp;
					if ( hist_tmp.sum() != 0 ) matchingSigs++;
					emptyBins += emptyBins_tmp;
					break;
				}
				//cout << (*myData)[j+k].name << " " << (*myData)[j+k].rc << endl;
				k++;
			} while (j+k<myData->size() && (*myData)[j].name==(*myData)[j+k].name);

			ResultT myResult;
			//ogzstream *fout = (*myData)[j].seqFile->out_results_fh;

			switch ((*myData)[j].seqFile->strandType){
			case FWD:
				myResult.numInstances = k;
				myResult.output_line = getResultString(hist,emptyBins,matchingSigs,k,(*myData)[j].name,FWD);
				myResultChunk->push_back(myResult);
				//*fout << myRes;
				break;
			case REV:
				myResult.numInstances = k;
				myResult.output_line = getResultString(histRC,emptyBinsRC,matchingSigsRC,k,(*myData)[j].name,REV);
				myResultChunk->push_back(myResult);
				//*fout << myRes;
				break;
			case FR:
				myResult.numInstances = k;
				myResult.output_line = getResultString(hist+histRC,emptyBins+emptyBinsRC,matchingSigs+matchingSigsRC,k,(*myData)[j].name,FR);
				myResultChunk->push_back(myResult);
				//*fout << myRes;
				break;
			case FR_sep:
				myResult.numInstances = k;
				myResult.output_line = getResultString(hist,emptyBins,matchingSigs,k/2,(*myData)[j].name,FWD);
				myResultChunk->push_back(myResult);
				//*fout << myRes;
				myResult.numInstances = k;
				myResult.output_line = getResultString(histRC,emptyBinsRC,matchingSigsRC,k/2,(*myData)[j].name,REV);
				myResultChunk->push_back(myResult);
				//*fout << myRes;
				break;
			default:
				break;
			}

			// meta analysis, only for screen output summary
			hist += histRC;
			unsigned sum = hist.sum();
			unsigned max = hist.max();

			++mNumSequences;

			if (sum!=0)
				mClassifiedInstances++;

			valarray<double> hist_t = hist;
			hist_t /= (k*mpParameters->mNumHashFunctions);

			hist_t /= sum;
			hist_t = hist.apply(changeNAN);

			metaHist += hist_t;

			for (unsigned i = 0; i<hist.size(); i++){
				if (hist[i] >= max ) {hist[i] = 1;} else { hist[i]=0.0;};
			}
			metaHistNum += hist; //.apply(indicator);

			j += k;
			break;
		}
		case INDEX_SIGCACHE:{

		}
		break;
		default:
			break;
		}
	} // while
}


string SeqClassifyManager::getResultString(histogramT hist,unsigned emptyBins, unsigned matchingSigs, unsigned numSigs, string name, strandTypeT strand){

	stringstream res;

	uint sum = hist.sum();
	uint max = hist.max();

	for (uint i = 0; i<hist.size(); i++){
		if (hist[i] / (numSigs*mpParameters->mNumHashFunctions) < mpParameters->mPureApproximateSim ) {
			hist[i] = 0.0;
		};
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

	string values;
	string indices;
	string maxIndices;
	for (unsigned i=0; i<hist.size();i++){
		if (hist[i]>0) {
			indices += std::to_string(i+1)+ ",";
			values += std::to_string((int)hist[i])+",";
		}
		if (hist[i]==max && max!=0)
			maxIndices += std::to_string(i+1)+",";
	}

	if (max!=0) {
		res << indices << "\t" << values << "\t";
		res << maxIndices;
	}

	res << "\n";
	return res.str();
}

void SeqClassifyManager::ClassifySeqs(){

	cout << endl << " *** Read sequences for classification and create their MinHash signatures *** " << endl << endl;
	cout << "hist size: " << GetHistogramSize() << endl;

	// prepare sequence set for classification
	SeqFileP mySet = std::make_shared<SeqFileT>();
	mySet->filename = mpParameters->mInputDataFileName;
	mySet->filetype = mpParameters->mFileTypeCode;
	mySet->groupGraphsBy=SEQ_NAME; // actually we check by InstanceT.name field for graphs from one seq
	mySet->checkUniqueSeqNames = true;
	mySet->signatureAction	= CLASSIFY;
	mySet->strandType			= FR;

	// get results file handle
	mySet->out_results_fh = PrepareResultsFile();

	metaHist.resize(GetHistogramSize());
	metaHist *= 0;
	metaHistNum.resize(GetHistogramSize());
	metaHistNum *= 0;

	mClassifiedInstances = 0;
	mNumSequences = 0;
	pb.Begin();
	SeqFilesT myList;
	myList.push_back(mySet);
	mHashBitMask = (2 << (mpParameters->mHashBitSize - 1)) - 1;
	cout << "bitmask " << mHashBitMask << endl;
	Classify_Signatures(myList);
	//LoadData_Threaded(myList);

	cout << "Classification finished - signatures=" << mSignatureCounter << " instances=" << mNumSequences << " classified=" << mClassifiedInstances<< endl;
	mySet->out_results_fh->close();

	/////////////////////////////////////////////////////////////////////////////
	// classification finished
	/////////////////////////////////////////////////////////////////////////////

	// metahistogram
	cout << endl << endl << "META histogram - classified seqs: " << setprecision(3) << (double)mClassifiedInstances/((double)mNumSequences) << " (" << mClassifiedInstances << ")" << endl;
	metaHist = metaHist/metaHist.sum();
	vector<pair<double,uint> > sortedHist;
	for (unsigned j=0; j<metaHist.size();j++){
		sortedHist.push_back(make_pair(-metaHist[j],j));
	}
	sort(sortedHist.begin(), sortedHist.end());
	for (unsigned j=0; j<std::min((unsigned)20,(unsigned)sortedHist.size());j++){
		multimap<uint, Data::BEDentryP>::iterator it = mIndexValue2Feature.find(sortedHist[j].second+1);
		uint num = mIndexValue2Feature.count(sortedHist[j].second+1);
		cout << setprecision(2) << j+1 << "\t" << -sortedHist[j].first << "\t" << setprecision(10) << metaHistNum[sortedHist[j].second] << "\t" << sortedHist[j].second+1 << "\t";
		if (it != mIndexValue2Feature.end()) cout << "feature\t"<< it->second->NAME << "\t#features=" << num;
		if (it != mIndexValue2Feature.end() && it->second->COLS.size()>=2) cout << "\t" << it->second->COLS[1];
		cout << endl;
	}
	cout << "SUM\t"<< metaHist.sum() << endl << endl;

	// bin size statistics
	if (mpParameters->mVerbose){
		valarray<double> indexHist;
		indexHist.resize(GetHistogramSize());
		indexHist *= 0;

		for (typename HistogramIndex::indexTy::const_iterator it = mInverseIndex.begin(); it!= mInverseIndex.end(); it++){
			for (typename HistogramIndex::indexSingleTy::const_iterator itBin = it->begin(); itBin!=it->end(); itBin++){
				indexHist[itBin->second[0]-1] += 1;
			}
		}

		for (unsigned i=0; i<indexHist.size();i++){
			cout << " index bin size " << i+1 << "\t" << indexHist[i] << "\t" << setprecision(5) << indexHist[i]/indexHist.sum() << endl;
		}
	}
}

ogzstream* SeqClassifyManager::PrepareResultsFile(){

	string resultsName = mpParameters->mInputDataFileName;
	const unsigned pos = mpParameters->mInputDataFileName.find_last_of("/");
	if (std::string::npos != pos)
		resultsName = mpParameters->mInputDataFileName.substr(pos+1);

	ogzstream* fout = new ogzstream((mpParameters->mDirectoryPath+resultsName+".classified.tab.gz").c_str(),std::ios::out);
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
	for (std::map<string,uint>::iterator it = mFeature2IndexValue.begin(); it != mFeature2IndexValue.end();++it) {
		*fout << "#HIST_IDX\t"<< it->second << "\t" << "feature\t"<< it->first;
		multimap<uint,Data::BEDentryP>::iterator it2 = mIndexValue2Feature.find(it->second);
		if (it2 != mIndexValue2Feature.end() && it2->second->COLS.size()>=2) *fout << "\t" << it2->second->COLS[1];
		*fout << endl;
	}
	*fout << "##" << endl;
	*fout << "#SEQ\tSTR\tSIGS\tSIG_HITS\tHF_HITS\tSUM\tMAX\tIDX\tVALS\tMAX_IDX"<< endl;
	return fout;
}

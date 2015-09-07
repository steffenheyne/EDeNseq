/*
 * SeqClassifyManager.cc
 *
 *  Created on: 20.03.2015
 *      Author: heyne
 */

#include "SeqClassifyManager.h"
#include "MinHashEncoder.h"


SeqClassifyManager::SeqClassifyManager(Parameters* apParameters, Data* apData):
HistogramIndex(apParameters,apData),pb(1000) //,mHistogramIndex(apPa
{
}


void SeqClassifyManager::Exec() {

	// load and prepare index list file to get files for indexing
	//	vector<SeqDataSet> fileList = mpData->LoadIndexDataList(mpParameters->mIndexDataList.c_str());

	SeqFileT mySet;
	mySet.filename          = mpParameters->mIndexSeqFile;
	mySet.filename_BED		= mpParameters->mIndexBedFile;
	mySet.filetype          = FASTA;
	mySet.updateIndex       = SEQ_FEATURE;
	mySet.updateSigCache    = false;
	Data::BEDdataP indexBED = mpData->LoadBEDfile(mpParameters->mIndexBedFile.c_str());
	mySet.dataBED           = indexBED;
	mySet.lastMetaIdx       = 0;

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

inline double changeNAN(double i) {if (std::isnan(i)) return 0.0; else return i;}

inline double indicator(double i) {if (i>0) return 1; else return 0;}

//inline double minSim(double i) {if (i<0.1) return 0; else return i;}

void SeqClassifyManager::finishUpdate(workQueueP& myData) {

	// return if we build the inverse index, otherwise we classify
	if (myData->seqFile->updateIndex != NONE) return;

	ogzstream *fout = myData->seqFile->out_results_fh;
	unsigned j = 0;

	while (j < myData->sigs.size()) {

		valarray<double> hist;
		hist.resize(GetHistogramSize());
		unsigned emptyBins = 0;
		unsigned k=0;
		unsigned matchingSigs = 0;
		do {
			valarray<double> hist_tmp;
			unsigned emptyBins_tmp;
			ComputeHistogram(myData->sigs[j+k],hist_tmp,emptyBins_tmp);
			hist += hist_tmp;
			if ( hist_tmp.sum() != 0 ) matchingSigs++;
			emptyBins += emptyBins_tmp;
			k++;
		} while (j+k<myData->sigs.size() && myData->names[j]==myData->names[j+k]);

		uint sum = hist.sum();
		uint max = hist.max();

		for (uint i = 0; i<hist.size(); i++){
			if (hist[i] / (k*mpParameters->mNumHashFunctions) < mpParameters->mPureApproximateSim ) {hist[i] = 0.0;};
		}

		mNumSequences++;
		if (sum!=0)
			mClassifiedInstances++;
		*fout << myData->names[j] << "\t"<< k << "\t" << matchingSigs <<  "\t" << (k*mpParameters->mNumHashFunctions)-emptyBins << "\t" << sum << "\t" << max << "\t";

		string values;
		string indices;
		string maxIndices;
		for (unsigned i=0; i<hist.size();i++){
			if (hist[i]>0) {
				indices += std::to_string(i+1)+ ",";
				values += std::to_string((int)hist[i])+",";
			}
			if (hist[i]==max && max!=0) maxIndices += std::to_string(i+1)+",";
		}

		if (max!=0) {
			*fout << indices << "\t" << values << "\t";
			*fout << maxIndices;
		}
		*fout << endl;

		hist /= (k*mpParameters->mNumHashFunctions);

		//		hist = hist.apply([tr] {if (i<0.1) return 0; else return i;});

		hist /= sum;
		//hist /= mpParameters->mNumHashFunctions-emptyBins;
		hist = hist.apply(changeNAN);
		//hist = hist.apply(minSim);

		metaHist += hist;
		metaHistNum += hist.apply(indicator);
		j += k;
	}
}


void SeqClassifyManager::ClassifySeqs(){

	cout << endl << " *** Read sequences for classification and create their MinHash signatures *** " << endl << endl;
	cout << "hist size: " << GetHistogramSize() << endl;

	// prepare sequence set for classification
	SeqFileP mySet = std::make_shared<SeqFileT>();
	mySet->filename = mpParameters->mInputDataFileName;
	mySet->filetype = mpParameters->mFileTypeCode;
	mySet->updateIndex=NONE;
	mySet->updateSigCache=false;

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

	LoadData_Threaded(myList);

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
	*fout << "#SEQ\tSIGS\tSIG_HITS\tHF_HITS\tSUM\tMAX\tIDX\tVALS\tMAX_IDX"<< endl;
	return fout;
}

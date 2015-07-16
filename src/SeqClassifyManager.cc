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
	SeqFilesT myList;
	SeqFileT mySet;

	mySet.filename          = mpParameters->mIndexSeqFile;
	mySet.filetype          = FASTA;
	mySet.updateIndex       = SEQ_FEATURE;
	mySet.updateSigCache    = false;
	Data::BEDdataP indexBED = mpData->LoadBEDfile(mpParameters->mIndexBedFile.c_str());
	mySet.dataBED           = indexBED;
	mySet.lastMetaIdx       = 0;

	mIndexDataSet = std::make_shared<SeqFileT>(mySet);

	myList.push_back(mIndexDataSet);
	string indexName = mpParameters->mIndexBedFile;
	const unsigned pos = mpParameters->mIndexBedFile.find_last_of("/");
	if (std::string::npos != pos)
		indexName = mpParameters->mIndexBedFile.substr(pos+1);

	// create/load new inverse MinHash index against that we can classify other sequences
	if (!std::ifstream(mpParameters->mIndexBedFile+".bhi").good()){
		cout << endl << " *** Creating inverse index *** "<< endl << endl;

		LoadData_Threaded(myList);
		SetHistogramSize(mIndexDataSet->lastMetaIdx);

		// write index
		if (!mpParameters->mNoIndexCacheFile){
			cout << "inverse index file : " << mpParameters->mIndexBedFile+".bhi" << endl;
			cout << " write index file ... ";
			OutputManager om((indexName+".bhi").c_str(), mpParameters->mDirectoryPath);
			writeBinaryIndex2(om.mOut,mInverseIndex);
			om.mOut.close();
			cout << endl;
		} else
			cout << "Index is NOT saved to file!"<< endl;
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
	}

	ClassifySeqs();
}

inline double changeNAN(double i) {if (std::isnan(i)) return 0.0; else return i;}

inline double indicator(double i) {if (i>0) return 1; else return 0;}

inline double minSim(double i) {if (i<0.1) return 0; else return i;}

void SeqClassifyManager::finishUpdate(workQueueP& myData) {

	ogzstream *fout = myData->seqFile->out_results_fh;

	for (unsigned j = 0; j < myData->sigs.size(); j++) {

		valarray<double> hist;
		unsigned emptyBins = 0;

		ComputeHistogram(myData->sigs[j],hist,emptyBins);

		unsigned sum = hist.sum();
		unsigned max = hist.max();

		if (sum!=0)
			mClassifiedInstances++;
		*fout << myData->names[j] << "\t"<<mpParameters->mNumHashFunctions-emptyBins << "\t" << sum << "\t";
		string values;
		string maxVals;
		for (unsigned i=0; i<hist.size();i++){
			if (hist[i]>0) values += std::to_string(i+1)+ ":" +std::to_string((int)hist[i])+"\t";
			if (hist[i]==max && max!=0) maxVals += std::to_string(i+1)+",";
		}
		if (max!=0) {
			*fout << "MAX:"<<max<<":"<<maxVals<< "\t";
			*fout << values;
		}
		*fout << endl;

		hist /= sum;
		//hist /= mpParameters->mNumHashFunctions-emptyBins;
		hist = hist.apply(changeNAN);
		//hist = hist.apply(minSim);

		metaHist += hist;
		metaHistNum += hist.apply(indicator);
	}
}

//void SeqClassifyManager::OutputResults() {
//	while (!done_output){
//
//		workQueueP myData;
//		unique_lock<mutex> lk(mut4);
//		cv4.wait(lk,[&]{if ( (done_output) || (finish_queue.try_pop( (myData) ))) return true; else return false;});
//		lk.unlock();
//
//		if (!done_output && myData->sigs.size()>0) {
//			cout << "results "<< myData->offset << endl;
//			mClassifiedInstances += myData->sigs.size();
//
//		}
//		cvr.notify_all();
//	}
//}

void SeqClassifyManager::ClassifySeqs(){

	// prepare sequence set for classification
	SeqFileP mySet = std::make_shared<SeqFileT>();
	mySet->filename = mpParameters->mInputDataFileName;
	mySet->filetype = mpParameters->mFileTypeCode;
	mySet->updateIndex=NONE;
	mySet->updateSigCache=false;

	string resultsName = mpParameters->mInputDataFileName;
	const unsigned pos = mpParameters->mInputDataFileName.find_last_of("/");
	if (std::string::npos != pos)
		resultsName = mpParameters->mInputDataFileName.substr(pos+1);

	ogzstream fout((mpParameters->mDirectoryPath+resultsName+".classified.tab.gz").c_str(),std::ios::out);
	//OutputManager out_results((resultsName+".classified.tab").c_str(), mpParameters->mDirectoryPath);
	mySet->out_results_fh = &fout;

	for (unsigned j=1; j<=GetHistogramSize();j++){
		multimap<uint, Data::BEDentryP>::iterator it = mIndexValue2Feature.find(j);
		uint num = mIndexValue2Feature.count(j);
		fout << "#HIST_IDX\t"<< j;
		if (it != mIndexValue2Feature.end()) fout << "\t" << "feature\t"<< it->second->NAME << "\t#BED_features\t" << num;
		if (it != mIndexValue2Feature.end() && it->second->COLS.size()>=2) fout << "\t" << it->second->COLS[1];
		fout << endl;
	}

	fout << "#SEQNAME\tHASH_SIG_HITS\tHIST_SUM\tHIST_VALUES"<< endl;

	SeqFilesT myList;
	myList.push_back(mySet);

	cout << endl << " *** Read sequences for classification and create their MinHash signatures *** " << endl << endl;

	cout << "hist size: " << GetHistogramSize() << endl;

	metaHist.resize(GetHistogramSize());
	metaHist *= 0;
	metaHistNum.resize(GetHistogramSize());
	metaHistNum *= 0;

	// currently we classify the whole seq always
	mpParameters->mSeqWindow = 0;

	mClassifiedInstances = 0;

	pb.Begin();
	LoadData_Threaded(myList);
	cout << "Load finished " << mSignatureCounter << " " << mClassifiedInstances<< endl;
	fout.close();

	// metahistogram
	cout << endl << endl << "META histogram - classified seqs: " << setprecision(3) << (double)mClassifiedInstances/((double)GetLoadedInstances()) << " (" << mClassifiedInstances << ")" << endl;
	metaHist = metaHist/metaHist.sum();
	vector<pair<double,uint> > sortedHist;
	for (unsigned j=0; j<metaHist.size();j++){
		sortedHist.push_back(make_pair(-metaHist[j],j));
	}
	sort(sortedHist.begin(), sortedHist.end());

	for (unsigned j=0; j<sortedHist.size();j++){
		multimap<uint, Data::BEDentryP>::iterator it = mIndexValue2Feature.find(sortedHist[j].second+1);
		uint num = mIndexValue2Feature.count(sortedHist[j].second+1);
		cout << setprecision(2) << j+1 << "\t" << -sortedHist[j].first << "\t" << setprecision(10) << metaHistNum[sortedHist[j].second] << "\t" << sortedHist[j].second+1 << "\t";
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

	for (unsigned i=0; i<indexHist.size();i++){
		cout << " index bin size " << i+1 << "\t" << indexHist[i] << "\t" << setprecision(5) << indexHist[i]/indexHist.sum() << endl;
	}
}

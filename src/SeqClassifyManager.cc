/*
 * SeqClassifyManager.cc
 *
 *  Created on: 20.03.2015
 *      Author: heyne
 */

#include "SeqClassifyManager.h"
#include "MinHashEncoder.h"


SeqClassifyManager::SeqClassifyManager(Parameters* apParameters, Data* apData):
HistogramIndex(apParameters, apData),pb(1000) //,mHistogramIndex(apParameters, apData)
{
	classifiedInstances = 0;
	//	mpParameters = apParameters;
	//	mpData=apData;
}


void SeqClassifyManager::Exec() {

	// load and prepare index list file to get files for indexing
	vector<SeqDataSet> fileList = mpData->LoadIndexDataList(mpParameters->mIndexDataList.c_str());
	string indexName = mpParameters->mIndexDataList;
	const unsigned pos = mpParameters->mIndexDataList.find_last_of("/");
	if (std::string::npos != pos)
		indexName = mpParameters->mIndexDataList.substr(pos+1);
	PrepareIndexDataSets(fileList);

	// create/load new inverse MinHash index against that we can classify other sequences
	if (!std::ifstream(mpParameters->mIndexDataList+".bhi").good()){
		cout << endl << " *** Creating inverse index *** "<< endl << endl;

		LoadData_Threaded(fileList,NULL);

		// write index
		if (!mpParameters->mNoIndexCacheFile){
			cout << "inverse index file : " << mpParameters->mIndexDataList+".bhi" << endl;
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
		cout << "inverse index file : " << mpParameters->mIndexDataList+".bhi" << endl << "read index ..." << endl;
		bool indexState = readBinaryIndex2(mpParameters->mIndexDataList+".bhi",mInverseIndex);
		cout << "finished! ";
		if (indexState == true){
			cout << " Index OK!" << endl;
		} else
			throw range_error("Cannot read index from file " + mpParameters->mIndexDataList+".bhi");
	}

	ClassifySeqs();
}

inline double changeNAN(double i) {if (std::isnan(i)) return 0.0; else return i;}

inline double indicator(double i) {if (i>0) return 1; else return 0;}

inline double minSim(double i) {if (i<0.1) return 0; else return i;}


void SeqClassifyManager::finishUpdate(workQueueT& myData,vector<vector<unsigned> >* myCache) {

	for (unsigned j = 0; j < myData->sigs.size(); j++) {

		valarray<double> hist;
		unsigned emptyBins = 0;
		ComputeHistogram(myData->sigs[j],hist,emptyBins);

		unsigned sum = hist.sum();
		hist /= sum;
		//hist /= mpParameters->mNumHashFunctions-emptyBins;
		hist = hist.apply(changeNAN);
		//hist = hist.apply(minSim);

		metaHist += hist;
		metaHistNum += hist.apply(indicator);

		if (sum!=0)
			classifiedInstances++;
//		pb.Count();
		//#ifdef USEMULTITHREAD
		//#pragma omp critical
		//#endif
		/*	if (emptyBins<mpParameters->mNumHashFunctions){
			cout << i << ":"<< emptyBins << "  \t";
			for (unsigned j=0; j<hist.size();j++){
				cout << setprecision(2) << hist[j] << "\t";
			}
			cout << endl;
		}*/
	}

}

void SeqClassifyManager::ClassifySeqs(){

	// prepare sequence set for classification
	SeqDataSet mySet;
	mySet.filename = mpParameters->mInputDataFileName.c_str();
	mySet.uIdx = 1;
	mySet.idx = 0;
	mySet.filetype = mpParameters->mFileTypeCode;
	mySet.desc = "seqs_to_classify";
	mySet.updateIndex=false;
	mySet.updateSigCache=false;
	vector<SeqDataSet> myList;
	myList.push_back(mySet);
	mpParameters->mSeqWindow = 0;
	cout << endl << " *** Read sequences for classification and create their MinHash signatures *** " << endl << endl;

	metaHist.resize(GetHistogramSize());
	metaHist *= 0;
	metaHistNum.resize(GetHistogramSize());
	metaHistNum *= 0;

	pb.Begin();
	LoadData_Threaded(myList,NULL);

	// metahistogram
	cout << endl << endl << "META histogram - classified seqs: " << setprecision(3) << (double)classifiedInstances/((double)GetLoadedInstances()) << " (" << classifiedInstances << ")" << endl;
	metaHist = metaHist/metaHist.sum();
	vector<pair<double,uint> > sortedHist;
	for (unsigned j=0; j<metaHist.size();j++){
		sortedHist.push_back(make_pair(-metaHist[j],j));
	}
	sort(sortedHist.begin(), sortedHist.end());

	for (unsigned j=0; j<sortedHist.size();j++){
		std::pair< std::multimap<uint,uint>::iterator, std::multimap<uint,uint>::iterator > ret;
		ret = mHistBin2DatasetIdx.equal_range(sortedHist[j].second);
		cout << setprecision(2) << j+1 << "\t" << -sortedHist[j].first << "\t" << setprecision(10) << metaHistNum[sortedHist[j].second] << "\t" << sortedHist[j].second << "\t";
		for (std::multimap<uint,uint>::iterator it = ret.first; it != ret.second; ++it ){
			cout << mIndexDataSets[it->second].desc <<";";
		}
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

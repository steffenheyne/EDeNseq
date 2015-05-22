/*
 * SeqClassifyManager.cc
 *
 *  Created on: 20.03.2015
 *      Author: heyne
 */

#include "SeqClassifyManager.h"
#include "MinHashEncoder.h"


SeqClassifyManager::SeqClassifyManager(Parameters* apParameters, Data* apData):
BaseManager(apParameters, apData), mHistogramIndex(apParameters, apData)
{
	mpParameters = apParameters;
	mpData=apData;
}


void SeqClassifyManager::Exec() {

	// load data source for index (genomes)
	vector<SeqDataSet> fileList = mpData->LoadIndexDataList(mpParameters->mIndexDataList.c_str());

	// sequence set for classification
	SeqDataSet mySet;
	mySet.filename = mpParameters->mInputDataFileName.c_str();
	mySet.idx = 1;
	mySet.filetype = mpParameters->mFileTypeCode;
	mySet.desc = "approx_cluster_set";
	mySet.updateIndex=false;
	mySet.updateSigCache=true;
	vector<SeqDataSet> myList;
	myList.push_back(mySet);

	// create/load inverse MinHash index against that we can classify other sequences
	const unsigned pos = mpParameters->mIndexDataList.find_last_of("/");
	string indexName = mpParameters->mIndexDataList;
	if (std::string::npos != pos)
		indexName = mpParameters->mIndexDataList.substr(pos+1);

	// create new index
	if (!std::ifstream(mpParameters->mIndexDataList+".bhi").good()){
		cout << endl << " *** Creating inverse index *** "<< endl << endl;
		mHistogramIndex.LoadDataIntoIndexThreaded(fileList,NULL);

		// write index
		if (mpParameters->mNoIndexCacheFile){
			cout << "inverse index file : " << mpParameters->mIndexDataList+".bhi" << endl;
			cout << " write index file ... ";
			OutputManager om((indexName+".bhi").c_str(), mpParameters->mDirectoryPath);
			mHistogramIndex.writeBinaryIndex2(om.mOut,mHistogramIndex.mInverseIndex);
			om.mOut.close();
			cout << endl;
		} else
			cout << "Index is NOT saved to file!"<< endl;
	} else {

		mHistogramIndex.mInverseIndex.clear();
		for (unsigned i=0;i<fileList.size(); i++){
			mHistogramIndex.mIndexDataSets.push_back(fileList[i]);
		}
		// read index
		cout << endl << " *** Read inverse index *** "<< endl << endl;
		cout << "inverse index file : " << mpParameters->mIndexDataList+".bhi" << endl << "read index ..." << endl;
		bool indexState = mHistogramIndex.readBinaryIndex2(mpParameters->mIndexDataList+".bhi",mHistogramIndex.mInverseIndex);
		cout << "finished! ";
		if (indexState == true){
			cout << " Index OK!" << endl;
		} else
			throw range_error("Cannot read index from file " + mpParameters->mIndexDataList+".bhi");
	}

	cout << endl << " *** Read sequences for classification and create their MinHash signatures *** " << endl << endl;
	mHistogramIndex.LoadDataIntoIndexThreaded(myList,NULL);

	ClassifySeqs();
}

double changeNAN(double i) {if (std::isnan(i)) return 0.0; else return i;}

void SeqClassifyManager::ClassifySeqs(){

	valarray<double> metaHist;
	metaHist.resize(mHistogramIndex.mIndexDataSets.size());
	metaHist *= 0;
	unsigned classifiedInstances = 0;
	ProgressBar pb(1000);
//#ifdef USEMULTITHREAD
//#pragma omp parallel for schedule(dynamic,100)
//#endif
	for (unsigned i = 0; i < mpData->Size(); ++i) {

		valarray<double> hist;
		unsigned emptyBins = 0;
		mHistogramIndex.ComputeHistogram(mHistogramIndex.ComputeHashSignature(i),hist,emptyBins);

		hist /= mpParameters->mNumHashFunctions-emptyBins;
		hist = hist.apply(changeNAN);
		metaHist += hist;

		if (hist.sum()!=0)
			classifiedInstances++;
		pb.Count();
	//#ifdef USEMULTITHREAD
	//#pragma omp critical
	//#endif
			{
			/*cout << i << ":"<< emptyBins << "  \t";
			/for (unsigned j=0; j<hist.size();j++){
					cout << setprecision(2) << hist[j] << "\t";
			}
			cout << endl;*/
		}
	}

	// metahistogram
	cout << endl << endl << "META histogram - classified seqs: " << setprecision(3) << (double)classifiedInstances/((double)mpData->Size()) << " (" << classifiedInstances << ")" << endl;
	metaHist = metaHist/metaHist.sum();
	for (unsigned j=0; j<metaHist.size();j++){
		cout << setprecision(2) << metaHist[j] << "\t";
	}
	cout << "  SUM "<< metaHist.sum() << endl << endl;
}

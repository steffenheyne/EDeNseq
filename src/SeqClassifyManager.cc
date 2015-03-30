/*
 * SeqClassifyManager.cc
 *
 *  Created on: 20.03.2015
 *      Author: heyne
 */

#include "SeqClassifyManager.h"


SeqClassifyManager::SeqClassifyManager(Parameters* apParameters, Data* apData):
BaseManager(apParameters, apData), mMinHashEncoder(apParameters, apData, MinHashEncoder::CLASSIFY) {
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
		mMinHashEncoder.LoadDataIntoIndexThreaded(fileList,NULL);

		// write index
		if (mpParameters->mNoIndexCacheFile){
			cout << "inverse index file : " << mpParameters->mIndexDataList+".bhi" << endl;
			cout << " write index file ... ";
			OutputManager om((indexName+".bhi").c_str(), mpParameters->mDirectoryPath);
			mpData->writeBinaryIndex(om.mOut,mMinHashEncoder.mInverseIndexPub);
			om.mOut.close();
			cout << endl;
		} else
			cout << "Index is NOT saved to file!"<< endl;
	} else {
		// read index
		cout << endl << " *** Read inverse index *** "<< endl << endl;
		cout << "inverse index file : " << mpParameters->mIndexDataList+".bhi" << endl << "read index ..." << endl;
		vector<umap_uint_vec_uint> tIndex;
		bool indexState = mpData->readBinaryIndex(mpParameters->mIndexDataList+".bhi",tIndex);
		cout << "finished! ";
		if (indexState == true){
			cout << " Index OK!" << endl;
		} else
			throw range_error("Cannot read index from file " + mpParameters->mIndexDataList+".bhi");
	}

	cout << endl << " *** Read sequences for classification and create their MinHash signatures *** " << endl << endl;
	mMinHashEncoder.LoadDataIntoIndexThreaded(myList,NULL);

	ClassifySeqs();

}

void SeqClassifyManager::ClassifySeqs(){



}

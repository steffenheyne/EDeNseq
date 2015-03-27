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

	vector<SeqDataSet> fileList = mpData->LoadIndexDataList(mpParameters->mIndexDataList.c_str());
	// load data (genomes) and create inverse MinHash index against that we can classify other sequences
	bool cacheSignatures = false;
	mMinHashEncoder.LoadDataIntoIndexThreaded(fileList,cacheSignatures,NULL);

	// save Index to disk
}

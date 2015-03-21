/*
 * SeqClassifyManager.cc
 *
 *  Created on: 20.03.2015
 *      Author: heyne
 */

#include "SeqClassifyManager.h"


SeqClassifyManager::SeqClassifyManager(Parameters* apParameters, Data* apData):
BaseManager(apParameters, apData), mMinHashEncoder(apParameters, apData) {
}

vector<SeqDataSet> SeqClassifyManager::LoadIndexDataList(string filename){

	vector<SeqDataSet> myList;
	bool valid_input = true;
	igzstream fin;
	string line;

	fin.open(filename.c_str());
	if (!fin)
		throw range_error("ERROR LoadData: Cannot open index data file: " + filename);
	while (!fin.eof() && valid_input) {
		SeqDataSet mySet;
		mySet.filetype=FASTA;
		if (fin >> mySet.idx >> mySet.filename >> mySet.desc){
			cout << "found file idx " << mySet.idx << "\t" << mySet.filename << "\t" << mySet.desc << endl;
			myList.push_back(mySet);
		}
		getline(fin, line);
	}
	if (!myList.size())
		throw range_error("ERROR LoadIndexData: No data found in " + filename + "!");
	return myList;
}

void SeqClassifyManager::Exec() {

	vector<SeqDataSet> fileList = LoadIndexDataList(mpParameters->mIndexDataList.c_str());
	mMinHashEncoder.LoadDataIntoIndexThreaded(fileList,false,NULL);
}

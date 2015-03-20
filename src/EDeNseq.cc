/*
 * EDeN_threadpool.cc
 *
 *  Created on: 26.02.2015
 *      Author: heyne
 */

#include <iostream>


#include "BaseManager.h"
#include "Data.h"
#include "Parameters.h"
//#include "NearestNeighbor.h"
//#include "Utility.h"
//#include "Kernel.h"
#include "MinHashEncoder.h"
#include "ClusterManager.h"

using namespace std;

//FlagsService& The_FlagsService = FlagsService::get_instance();


class MetaGenomeManager: public BaseManager {

public:

	MinHashEncoder mMinHashEncoder;

	MetaGenomeManager(Parameters* apParameters, Data* apData):
		BaseManager(apParameters, apData), mMinHashEncoder() {
		Init(apParameters, apData);
	}

	void Init(Parameters* apParameters, Data* apData){
		BaseManager::Init(apParameters, apData);
		mMinHashEncoder.Init(apParameters, apData);
	}

	vector<SeqDataSet> LoadIndexDataList(string filename){

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

	void Exec() {

		vector<SeqDataSet> fileList = LoadIndexDataList(mpParameters->mIndexDataList.c_str());
		mMinHashEncoder.LoadDataIntoIndexThreaded(fileList,false,NULL);
	}

};

int main(int argc, const char **argv) {

	Parameters mParameters;
	Data mData;
	mParameters.Init(argc, argv);
	srand(mParameters.mRandomSeed);
	mData.Init(&mParameters);

	MetaGenomeManager meta_genome_manager(&mParameters,&mData);
	meta_genome_manager.Exec();

	return 0;
}

/*
 * EDeNseq.cc
 *
 *  Created on: 26.02.2015
 *      Author: heyne
 */

#include <iostream>

#include "BaseManager.h"
#include "Data.h"
#include "Parameters.h"
#include "SeqClusterManager.h"
#include "SeqClassifyManager.h"
#include "MinHashEncoder.h"

using namespace std;

class Dispatcher {
protected:
	Parameters mParameters;
	Data mData;

public:
	Dispatcher() {
	}

	void Init(int argc, const char **argv) {
		mParameters.Init(argc, argv);
		srand(mParameters.mRandomSeed);
		mData.Init(&mParameters);
	}

	void Exec() {
		ProgressBar pb(100);

		cout << SEP << endl << PROG_NAME << endl << "Version: " << PROG_VERSION << endl << "Last Update: " << PROG_DATE << endl << PROG_CREDIT << endl << SEP << endl;

		switch (mParameters.mActionCode) {
		case CLASSIFY:{
			SeqClassifyManager seq_classify_manager(&mParameters,&mData);
			seq_classify_manager.Exec();
		}
		break;
		case CLUSTER:{
			SeqClusterManager cluster_manager(&mParameters, &mData);
			cluster_manager.Exec();
		}
		break;
		default:
			throw range_error("ERROR: Unknown action parameter: " + mParameters.mAction);
		}
		pb.Count(1);
		cout << "Total run-time:" << endl;
	}
};



int main(int argc, const char **argv) {
	try {
		Eigen::initParallel();
		Dispatcher myDispatcher;
		myDispatcher.Init(argc,argv);
		myDispatcher.Exec();
	} catch (exception& e) {
		cerr << e.what() << endl;
	}

	return 0;
}

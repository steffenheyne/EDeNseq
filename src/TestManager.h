/*
 * TestManager.h
 *
 *  Created on: 30.07.2015
 *      Author: Steffen Heyne
 */

#ifndef TESTMANAGER_H_
#define TESTMANAGER_H_

#include "MinHashEncoder.h"
#include "BaseManager.h"
#include "Data.h"
#include "Parameters.h"

#include <math.h>

class TestManager: public HistogramIndex {

public:
	TestManager(Parameters* apParameters, Data* apData);

	unsigned mNumSequences;
	valarray<double> metaHist;
	valarray<double> metaHistNum;
	unsigned mClassifiedInstances;
	ProgressBar pb;

	threadsafe_queue<workQueueP> finish_queue;
	std::atomic_bool done_output;

	void finishUpdate(workQueueP& myData);
	void Exec();
	void ClassifySeqs();
	ogzstream* PrepareResultsFile();
};

#endif /* TESTMANAGER_H_ */

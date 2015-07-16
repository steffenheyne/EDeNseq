/*
 * SeqClassifyManager.h
 *
 *  Created on: 20.03.2015
 *      Author: Steffen Heyne
 */

#ifndef SEQCLASSIFYMANAGER_H_
#define SEQCLASSIFYMANAGER_H_

#include "MinHashEncoder.h"
#include "BaseManager.h"
#include "Data.h"
#include "Parameters.h"

#include <math.h>

class SeqClassifyManager: public HistogramIndex {

public:
	SeqClassifyManager(Parameters* apParameters, Data* apData);

	unsigned mNumSequences;
	valarray<double> metaHist;
	valarray<double> metaHistNum;
	unsigned mClassifiedInstances;
	ProgressBar pb;

	threadsafe_queue<workQueueP> finish_queue;
	std::atomic_bool done_output;

//	std::condition_variable cv4;
//	std::condition_variable cvr;
//	mutable std::mutex mut4;
//	mutable std::mutex mutr;

	//void OutputResults();
	void finishUpdate(workQueueP& myData);
	void Exec();
	void ClassifySeqs();

};

#endif /* SEQCLASSIFYMANAGER_H_ */

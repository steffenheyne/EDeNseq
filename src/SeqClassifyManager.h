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

	std::atomic_bool done_output;

	void finishUpdate(ChunkP& myData);
	void Exec();
	void ClassifySeqs();
	ogzstream* PrepareResultsFile();

	//inline double minSim(double i) { if (i<mpParameters->mPureApproximateSim) return 0; else return i; };
};

#endif /* SEQCLASSIFYMANAGER_H_ */

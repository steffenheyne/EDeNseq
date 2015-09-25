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

	void 			Exec();
	void 			finishUpdate(ChunkP& myData);
	void 			ClassifySeqs();
	string		getResultString(histogramT hist,unsigned emptyBins, unsigned matchingSigs, unsigned numSigs, string name, strandTypeT strand);
	ogzstream* 	PrepareResultsFile();

	//inline double minSim(double i) { if (i<mpParameters->mPureApproximateSim) return 0; else return i; };
};

#endif /* SEQCLASSIFYMANAGER_H_ */

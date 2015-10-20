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

	struct resultS{
		string output_line;
		unsigned numInstances;
	};

	typedef resultS ResultT;
	typedef vector<ResultT> ResultChunkT;
	typedef std::shared_ptr<ResultChunkT> ResultChunkP;

	std::atomic_uint mNumSequences;
	histogramT metaHist;
	histogramT metaHistNum;
	unsigned mClassifiedInstances;

	std::atomic_bool done_output;

	threadsafe_queue<ResultChunkP> res_queue;
	mutable std::mutex mut_res;
   std::condition_variable cv_res;
	std::atomic_uint mResultCounter;


	void 			Exec();
	void 			finishUpdate(ChunkP& myData);
	void 			finishUpdate(ChunkP& myData, ResultChunkP& myResult);
	void 			finishUpdate(ChunkP& myData, unsigned& min, unsigned& max);

	void 			ClassifySeqs();
	void 			Classify_Signatures(SeqFilesT& myFiles);
	void 			worker_Classify(int numWorkers);
	void 			finisher_Results(ogzstream* fout_res);
	void 			getResultString(string& resT, histogramT hist, unsigned emptyBins, unsigned matchingSigs, unsigned numSigs, string& name, strandTypeT strand);
	ogzstream* 	PrepareResultsFile();

	//inline double minSim(double i) { if (i<mpParameters->mPureApproximateSim) return 0; else return i; };
};

#endif /* SEQCLASSIFYMANAGER_H_ */

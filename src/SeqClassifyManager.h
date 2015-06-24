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

class SeqClassifyManager: public BaseManager {

public:
	SeqClassifyManager(Parameters* apParameters, Data* apData);

	HistogramIndex mHistogramIndex;
	unsigned mNumSequences;

	void Exec();
	void ClassifySeqs();

};

#endif /* SEQCLASSIFYMANAGER_H_ */

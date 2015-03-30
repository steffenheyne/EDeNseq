/*
 * SeqClassifyManager.h
 *
 *  Created on: 20.03.2015
 *      Author: Steffen Heyne
 */

#ifndef SEQCLASSIFYMANAGER_H_
#define SEQCLASSIFYMANAGER_H_

#include "BaseManager.h"
#include "Data.h"
#include "Parameters.h"
#include "MinHashEncoder.h"

class SeqClassifyManager: public BaseManager {
private:
	MinHashEncoder mMinHashEncoder;

public:
	SeqClassifyManager(Parameters* apParameters, Data* apData);
	void Exec();
	void SeqClassifyManager::ClassifySeqs();
};

#endif /* SEQCLASSIFYMANAGER_H_ */

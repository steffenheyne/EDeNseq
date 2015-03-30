/* -*- mode:c++ -*- */
#ifndef CLUSTER_MANAGER_H
#define CLUSTER_MANAGER_H

#include "NearestNeighbor.h"

using namespace std;

class SeqClusterManager: public BaseManager {
public:
	SeqClusterManager(Parameters* apParameters, Data* apData);
   void Exec();
protected:
   SeqDataSet mDataSet;
   NearestNeighbor mNearestNeighbor;
   void DenseCluster();
};

#endif /* CLUSTER_MANAGER_H */

/* -*- mode:c++ -*- */
#ifndef CLUSTER_MANAGER_H
#define CLUSTER_MANAGER_H

#include "MinHashEncoder.h"


using namespace std;

class SeqClusterManager : public NeighborhoodIndex {
public:
	SeqClusterManager(Parameters* apParameters, Data* apData);
   void Exec();
protected:
   //SeqDataSet mDataSet;
   //NearestNeighbor mNearestNeighbor;
   void DenseCluster();
};

#endif /* CLUSTER_MANAGER_H */

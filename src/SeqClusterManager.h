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
   NearestNeighbor mNearestNeighbor;
   void DenseCluster(ostream& out_c, ostream& out_n);
};

#endif /* CLUSTER_MANAGER_H */

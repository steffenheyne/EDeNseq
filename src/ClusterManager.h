/* -*- mode:c++ -*- */
#ifndef CLUSTER_MANAGER_H
#define CLUSTER_MANAGER_H

#include "NearestNeighbor.h"

using namespace std;

class ClusterManager: public BaseManager {
protected:
	NearestNeighbor mNearestNeighbor;
public:
	ClusterManager(Parameters* apParameters, Data* apData);
	void Exec();
};

class ApproxDensityClusterManager: public ClusterManager {
public:
	ApproxDensityClusterManager(Parameters* apParameters, Data* apData);
   void Exec();
protected:
	void DenseCluster(ostream& out_c, ostream& out_n);
};

#endif /* CLUSTER_MANAGER_H */

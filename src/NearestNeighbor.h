/* -*- mode:c++ -*- */
#ifndef NEAREST_NEIGHBOR_H
#define NEAREST_NEIGHBOR_H

#include "BaseManager.h"
#include "MinHashEncoder.h"

using namespace std;

class NearestNeighbor: public BaseManager {
public:
	NeighborhoodIndex mMinHashEncoder;
protected:
	vector<vector<unsigned> > 			mNeighborhoodCache;
	vector<umap_uint_int> 				mNeighborhoodCacheExt;
	vector<pair<unsigned,double> >	mNeighborhoodCacheInfo;

public:
	NearestNeighbor(Parameters* apParameters, Data* apData);
	void CacheReset();

	vector<unsigned> ComputeNeighborhood(unsigned aID, unsigned& collisions, double& density);
	umap_uint_int    ComputeNeighborhoodExt(unsigned aID, unsigned& collisions, double& density);
	vector<unsigned> ComputeApproximateNeighborhood(unsigned aID, unsigned& collisions, double& density);
	umap_uint_int 	  ComputeApproximateNeighborhoodExt(unsigned aID, unsigned& collisions, double& density);

	//	vector<unsigned> ComputeNeighborhood(unsigned aID);
	//	vector<unsigned> ComputeApproximateNeighborhood(SVector& aX);
	//	vector<unsigned> ComputeSharedNeighborhood(unsigned aID);
	//	double 			  ComputeSharedNeighborhoodSimilarity(unsigned aI, unsigned aJ);
	//	unsigned 		  ComputeNeighborhoodIntersection(unsigned aI, unsigned aJ);
	//	void             ComputeSharedNeighborhoods();
	//	vector<unsigned> ComputeNeighborhood(SVector& aX);
};

#endif /* NEAREST_NEIGHBOR_H */

#include "NearestNeighbor.h"

NearestNeighbor::NearestNeighbor(Parameters* apParameters, Data* apData) :
BaseManager(apParameters, apData) {
	Init(apParameters, apData);
}

void NearestNeighbor::Init(Parameters* apParameters, Data* apData) {
	BaseManager::Init(apParameters, apData);
	mMinHashEncoder.Init(apParameters, apData);
}

void NearestNeighbor::CacheReset() {
	cout << "... nearest neighbor cache reset ..." << endl;
	mNeighborhoodCache.clear();
	mNeighborhoodCacheExt.clear();
	mNeighborhoodCacheInfo.clear();
	if (mpData->IsDataLoaded() == false)
		throw range_error("ERROR: Cannot clean cache if data is not loaded");
	mNeighborhoodCache.resize(mpData->Size());
	mNeighborhoodCacheExt.resize(mpData->Size());
	mNeighborhoodCacheInfo.resize(mpData->Size());
}

vector<unsigned> NearestNeighbor::ComputeNeighborhood(unsigned aID) {

	unsigned collisions;
	double density;

	vector<unsigned> neighborhood_list = ComputeNeighborhood(aID,collisions,density);

	return neighborhood_list;
}


vector<unsigned> NearestNeighbor::ComputeNeighborhood(unsigned aID, unsigned& collisions, double& density) {
	//cache neighborhoods (if opted for)
	vector<unsigned> neighborhood_list;
	if (mNeighborhoodCache[aID].size() != 0) {
			neighborhood_list = mNeighborhoodCache[aID];
			pair<unsigned, double> tmp = mNeighborhoodCacheInfo[aID];
			collisions = tmp.first;
			density = tmp.second;
		} else {
				neighborhood_list = ComputeApproximateNeighborhood(aID, collisions, density);
				mNeighborhoodCacheInfo[aID] = make_pair(collisions,density);
			mNeighborhoodCache[aID] = neighborhood_list;
		}
	return neighborhood_list;
}

umap_uint_int NearestNeighbor::ComputeNeighborhoodExt(unsigned aID, unsigned& collisions, double& density) {
	//cache neighborhoods (if opted for)
	umap_uint_int neighborhood_list;
		if (mNeighborhoodCacheExt[aID].size() != 0) {
			neighborhood_list = mNeighborhoodCacheExt[aID];
			pair<unsigned, double> tmp = mNeighborhoodCacheInfo[aID];
			collisions = tmp.first;
			density = tmp.second;
		} else {
			neighborhood_list = ComputeApproximateNeighborhoodExt(aID,collisions,density);
			mNeighborhoodCacheExt[aID] = neighborhood_list;
			mNeighborhoodCacheInfo[aID] = make_pair(collisions,density);
		}
	return neighborhood_list;
}


vector<unsigned> NearestNeighbor::ComputeApproximateNeighborhood(unsigned aID, unsigned& collisions, double& density) {
	vector<unsigned> signature = mMinHashEncoder.ComputeHashSignature(aID);
	vector<unsigned> approximate_neighborhood = mMinHashEncoder.ComputeApproximateNeighborhood(signature, collisions, density);
//	if (mpParameters->mPureApproximateSim==0){
//		vector<unsigned> approximate_true_neighborhood = ComputeTrueSubNeighborhood(aID, approximate_neighborhood);
//		return approximate_true_neighborhood;
//	} else {
		return approximate_neighborhood;
//	}
}

umap_uint_int NearestNeighbor::ComputeApproximateNeighborhoodExt(unsigned aID, unsigned& collisions, double& density) {
	vector<unsigned> signature = mMinHashEncoder.ComputeHashSignature(aID);
	umap_uint_int approximate_neighborhood = mMinHashEncoder.ComputeApproximateNeighborhoodExt(signature, collisions, density);
	return approximate_neighborhood;
}


/**
 Computes the fraction of neighbors that are common between instance I and J
 */
double NearestNeighbor::ComputeSharedNeighborhoodSimilarity(unsigned aI, unsigned aJ) {
	vector<unsigned> neighborhood_i = ComputeNeighborhood(aI);
	vector<unsigned> neighborhood_j = ComputeNeighborhood(aJ);
	unsigned intersection_size = ComputeNeighborhoodIntersection(aI, aJ);
	double shared_neighborhood_value = (double) intersection_size / sqrt((double) neighborhood_i.size() * (double) neighborhood_j.size());
	return shared_neighborhood_value;
}

vector<unsigned> NearestNeighbor::ComputeSharedNeighborhood(unsigned aID) {
	vector<unsigned> shared_neighborhood;
	vector<unsigned> neighborhood = ComputeNeighborhood(aID);
	//for each element in the neighborhood consider their neighborhood and check if there is aID, if yes then the element can stay in the neighborhood otherwise not
	for (unsigned i = 0; i < neighborhood.size(); ++i) {
		unsigned nn_id = neighborhood[i];
		vector<unsigned> nn_neighborhood = ComputeNeighborhood(nn_id);
		bool is_present = false;
		for (unsigned j = 0; j < nn_neighborhood.size() && is_present == false; j++) {
			if (nn_neighborhood[j] == aID)
				is_present = true;
		}
		if (is_present == true)
			shared_neighborhood.push_back(nn_id);
	}
	return shared_neighborhood;
}

unsigned NearestNeighbor::ComputeNeighborhoodIntersection(unsigned aI, unsigned aJ) {
	vector<unsigned> neighborhood_i = ComputeNeighborhood(aI);
	set<unsigned> neighborhood_i_set;
	neighborhood_i_set.insert(neighborhood_i.begin(), neighborhood_i.end());
	vector<unsigned> neighborhood_j = ComputeNeighborhood(aJ);
	set<unsigned> neighborhood_j_set;
	neighborhood_j_set.insert(neighborhood_j.begin(), neighborhood_j.end());

	set<unsigned> intersection;
	set_intersection(neighborhood_i.begin(), neighborhood_i.end(), neighborhood_j_set.begin(), neighborhood_j_set.end(), inserter(intersection, intersection.begin()));
	return intersection.size();
}

void NearestNeighbor::ComputeSharedNeighborhoods() {
	cout << "Computing shared neighborhoods" << endl;
	ProgressBar pb;
	vector<vector<unsigned> > new_neighborhood_cache(mpData->Size());
#ifdef USEMULTITHREAD
#pragma omp parallel for schedule(dynamic)
#endif
	for (unsigned i = 0; i < mpData->Size(); ++i) {
		vector<unsigned> neighborhood = ComputeSharedNeighborhood(i);
		new_neighborhood_cache[i] = neighborhood;
		pb.Count();
	}
	mNeighborhoodCache = new_neighborhood_cache;
}

vector<unsigned> NearestNeighbor::ComputeNeighborhood(SVector& aX) {
	vector<unsigned> neighborhood_list;
		neighborhood_list = ComputeApproximateNeighborhood(aX);
	return neighborhood_list;
}

vector<unsigned> NearestNeighbor::ComputeApproximateNeighborhood(SVector& aX) {
	unsigned collisions;
	double density;
	vector<unsigned> signature = mMinHashEncoder.ComputeHashSignature(aX);
	vector<unsigned> approximate_neighborhood = mMinHashEncoder.ComputeApproximateNeighborhood(signature,collisions,density);
	return approximate_neighborhood;
}

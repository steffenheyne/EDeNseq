/* -*- mode:c++ -*- */
#ifndef MIN_HASH_ENCODER_H
#define MIN_HASH_ENCODER_H

#include "Utility.h"
#include "Parameters.h"
//#include "Kernel.h"
#include "Data.h"
#include "BaseManager.h"
#include <streambuf>

using namespace std;


class MinHashEncoder {
protected:
	Parameters* mpParameters;
	Data* mpData;
	vector<umap_uint_vec_uint> mInverseIndex;
	vector<vector<unsigned> > mMinHashCache;
	unsigned numKeys;
	unsigned numFullBins;

	struct graphQueueS {
		vector<GraphClass> gr;
		unsigned id;
		unsigned offset;
	};

	struct sigQueueS {
		vector<vector<unsigned int> > sigs;
		unsigned id;
		unsigned offset;
	};

	typedef std::shared_ptr<graphQueueS> graphQueueT;
	typedef std::shared_ptr<sigQueueS> sigQueueT;

	mutable std::mutex mut;
	std::condition_variable cv;
	vector<std::thread> threads;
	threadsafe_queue<SeqDataSet> readFile_queue;
	threadsafe_queue<graphQueueT > graph_queue;
	threadsafe_queue<sigQueueT > sig_queue;

	std::atomic_bool done;
	std::atomic_int files_done;
	std::atomic_int instance_counter;
	std::atomic_int signature_counter;

	void worker_readFiles(int numWorkers);
	void worker_Graph2Signature(int id);
	void finisher(bool idxUpdate, vector<vector<unsigned> >* myCache);

	void generate_feature_vector(const GraphClass& aG, SVector& x);
	void InitFeatureCache(unsigned aSize, unsigned aRadius);
	vector<unsigned> HashFuncNSPDK(const string& aString, unsigned aStart, unsigned aMaxRadius, unsigned aBitMask);
	unsigned HashFuncNSPDK(const vector<unsigned>& aList, unsigned aBitMask);
	bool SetGraphFromFASTAFile(istream& in, GraphClass& oG, string& currSeq);
	unsigned mHashBitMask;

public:
	MinHashEncoder();
	MinHashEncoder(Parameters* apParameters, Data* apData);
	void Init(Parameters* apParameters, Data* apData);
	void CacheReset();
	void ComputeInverseIndex();
	void UpdateInverseIndex(vector<unsigned>& aSignature, unsigned aIndex);
	void CleanUpInverseIndex();
	void LoadDataIntoIndex();
	void LoadDataIntoIndexThreaded(vector<SeqDataSet> myFiles,bool useMinHashCache, vector<vector<unsigned> >* myCache);
	vector<unsigned> ComputeHashSignature(unsigned aID);
	vector<unsigned> ComputeHashSignature(SVector& aX);
	vector<unsigned> ComputeHashSignatureSize(vector<unsigned>& aSignature);
	void             ComputeApproximateNeighborhoodCore(const vector<unsigned>& aSignature, umap_uint_int& neighborhood, unsigned& collisions);
	umap_uint_int    ComputeApproximateNeighborhoodExt(const vector<unsigned>& aSignature, unsigned& collisions, double& density);

	vector<unsigned> ComputeApproximateNeighborhood(const vector<unsigned>& aSignature, unsigned& collisions, double& density);

	vector<unsigned> TrimNeighborhood(umap_uint_int& aNeighborhood, unsigned aNeighborhoodSize, unsigned collisions, double& density);

	double           ComputeApproximateSim(const unsigned& aID, const unsigned& bID);
	pair<unsigned,unsigned> ComputeApproximateSim(const unsigned& aID, const vector<unsigned>& bSignature);
};

#endif /* MIN_HASH_ENCODER_H */

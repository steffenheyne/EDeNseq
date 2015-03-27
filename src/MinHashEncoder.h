/* -*- mode:c++ -*- */
#ifndef MIN_HASH_ENCODER_H
#define MIN_HASH_ENCODER_H

#include "Utility.h"
#include "Parameters.h"
#include "Data.h"
#include "BaseManager.h"
#include <streambuf>

using namespace std;


class MinHashEncoder {

public:
	enum INDEXTypeE {
			CLUSTER, CLASSIFY
		};

		typedef INDEXTypeE INDEXType;

protected:

	INDEXType	indexType;

	Parameters* mpParameters;
	Data* mpData;
	vector<umap_uint_vec_uint> mInverseIndex;
	vector<vector<unsigned> > mMinHashCache;
	unsigned numKeys;
	unsigned numFullBins;

	vector<SeqDataSet>	mIndexDataSets;
	vector<SeqDataSet>	mClassifyDataSets;

	struct workQueueS {
		vector<GraphClass> gr;
		vector<vector<unsigned> > sigs;
		unsigned offset;
		SeqDataSet* dataSet;
	};

	typedef std::shared_ptr<SeqDataSet> SeqDataSetP;
	typedef std::shared_ptr<workQueueS> workQueueT;

	mutable std::mutex mutm;
	mutable std::mutex mut1;
	mutable std::mutex mut2;
	mutable std::mutex mut3;


	std::condition_variable cvm;
	std::condition_variable cv1;
	std::condition_variable cv2;
	std::condition_variable cv3;

	vector<std::thread> threads;
	threadsafe_queue<SeqDataSetP> readFile_queue;
	threadsafe_queue<workQueueT> graph_queue;
	threadsafe_queue<workQueueT> sig_queue;

	std::atomic_bool done;
	std::atomic_int files_done;
	std::atomic_int instance_counter;
	std::atomic_int signature_counter;

	unsigned mHashBitMask;

	void worker_readFiles(int numWorkers);
	void worker_Graph2Signature();
	void finisher(vector<vector<unsigned> >* myCache);
	void generate_feature_vector(const GraphClass& aG, SVector& x);
	void InitFeatureCache(unsigned aSize, unsigned aRadius);
	vector<unsigned> HashFuncNSPDK(const string& aString, unsigned aStart, unsigned aMaxRadius, unsigned aBitMask);
	unsigned HashFuncNSPDK(const vector<unsigned>& aList, unsigned aBitMask);
	bool SetGraphFromFASTAFile(istream& in, GraphClass& oG, string& currSeq);

public:
	MinHashEncoder();
	MinHashEncoder(Parameters* apParameters, Data* apData, INDEXType apIndexType=CLUSTER);
	void Init(Parameters* apParameters, Data* apData, INDEXType apIndexType=CLUSTER);
	void UpdateInverseIndex(vector<unsigned>& aSignature, unsigned aIndex);
	void UpdateInverseIndexHist(vector<unsigned>& aSignature, unsigned aIndex);
	void CleanUpInverseIndex();
	void LoadDataIntoIndexThreaded(vector<SeqDataSet>& myFiles,bool useMinHashCache, vector<vector<unsigned> >* myCache);
	vector<unsigned> ComputeHashSignature(unsigned aID);
	vector<unsigned> ComputeHashSignature(SVector& aX);
	vector<unsigned> ComputeHashSignatureSize(vector<unsigned>& aSignature);
	void             ComputeApproximateNeighborhoodCore(const vector<unsigned>& aSignature, umap_uint_int& neighborhood, unsigned& collisions);
	umap_uint_int    ComputeApproximateNeighborhoodExt(const vector<unsigned>& aSignature, unsigned& collisions, double& density);
	vector<unsigned> ComputeApproximateNeighborhood(const vector<unsigned>& aSignature, unsigned& collisions, double& density);
	vector<unsigned> TrimNeighborhood(umap_uint_int& aNeighborhood, unsigned collisions, double& density);
	double           ComputeApproximateSim(const unsigned& aID, const unsigned& bID);
	pair<unsigned,unsigned> ComputeApproximateSim(const unsigned& aID, const vector<unsigned>& bSignature);
};

#endif /* MIN_HASH_ENCODER_H */

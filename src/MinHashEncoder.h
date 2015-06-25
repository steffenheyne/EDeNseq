/* -*- mode:c++ -*- */
#ifndef MIN_HASH_ENCODER_H
#define MIN_HASH_ENCODER_H

#include <streambuf>
#include <valarray>
#include <list>
#include "Utility.h"
#include "Parameters.h"
#include "Data.h"
#include "BaseManager.h"
#include "sparsehash-2.0.2/sparsehash/sparse_hash_map"

using namespace std;

class MinHashEncoder {

public:

	enum INDEXTypeE {
		CLUSTER, CLASSIFY
	};

	struct workQueueS {
		vector<GraphClass> gr;
		vector<vector<unsigned> > sigs;
		vector<string> names;
		unsigned offset;
		SeqDataSet* dataSet;
	};

	typedef INDEXTypeE INDEXType;
	typedef std::shared_ptr<SeqDataSet> SeqDataSetP;
	typedef std::shared_ptr<workQueueS> workQueueT;

	vector<SeqDataSet>	mIndexDataSets;
	std::tr1::unordered_map<string, unsigned> name2idxMap;
	vector<string>	idx2nameMap;

	Parameters* mpParameters;
	Data* mpData;

protected:

	INDEXType	indexType;

	unsigned numKeys;
	unsigned numFullBins;

private:
	vector<vector<unsigned> > mMinHashCache;

	mutable std::mutex mutm;
	mutable std::mutex mut1;
	mutable std::mutex mut2;
	mutable std::mutex mut3;

	std::condition_variable cvm;
	std::condition_variable cv1;
	std::condition_variable cv2;
	std::condition_variable cv3;

	threadsafe_queue<SeqDataSetP> readFile_queue;
	threadsafe_queue<workQueueT> graph_queue;
	threadsafe_queue<workQueueT> sig_queue;

	std::atomic_bool done;
	std::atomic_int files_done;
	std::atomic_int mInstanceCounter;
	std::atomic_int mSignatureCounter;

	unsigned mHashBitMask;

	void 					worker_readFiles(int numWorkers);
	void					worker_Graph2Signature();
	void 					finisher(vector<vector<unsigned> >* myCache);
	void 					generate_feature_vector(const GraphClass& aG, SVector& x);
	vector<unsigned>	HashFuncNSPDK(const string& aString, unsigned aStart, unsigned aMaxRadius, unsigned aBitMask);

public:

	MinHashEncoder(Parameters* apParameters, Data* apData, INDEXType apIndexType=CLUSTER);
	virtual	~MinHashEncoder();
	void		Init(Parameters* apParameters, Data* apData, INDEXType apIndexType=CLUSTER);

	virtual void 		UpdateInverseIndex(vector<unsigned>& aSignature, unsigned aIndex) {};
	virtual void 		finishUpdate(workQueueT& myData, vector<vector<unsigned> >* myCache) {};

	void 					LoadData_Threaded(vector<SeqDataSet>& myFiles, vector<vector<unsigned> >* myCache);
	unsigned				GetLoadedInstances();
	vector<unsigned>& ComputeHashSignature(unsigned aID);
	vector<unsigned>	ComputeHashSignature(SVector& aX);
};

class NeighborhoodIndex : public MinHashEncoder
{
public:
	NeighborhoodIndex(Parameters* apParameters, Data* apData)
		:MinHashEncoder(apParameters,apData,CLUSTER)
	{
		mInverseIndex.clear();
		for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k){
			mInverseIndex.push_back(indexSingleTy());
		}
	}

	// inverse index
	typedef unsigned binKeyTy;
	typedef vector<binKeyTy>	indexBinTy;
	typedef std::tr1::unordered_map<unsigned,indexBinTy> indexSingleTy;
	typedef vector<indexSingleTy> indexTy;

	const binKeyTy MAXBINKEY = std::numeric_limits<binKeyTy>::max();

	indexTy mInverseIndex;

	void 				  UpdateInverseIndex(vector<unsigned>& aSignature, unsigned aIndex);
	void             ComputeApproximateNeighborhoodCore(const vector<unsigned>& aSignature, umap_uint_int& neighborhood, unsigned& collisions);
	umap_uint_int    ComputeApproximateNeighborhoodExt(const vector<unsigned>& aSignature, unsigned& collisions, double& density);
	vector<unsigned> ComputeApproximateNeighborhood(const vector<unsigned>& aSignature, unsigned& collisions, double& density);

	vector<unsigned> TrimNeighborhood(umap_uint_int& aNeighborhood, unsigned collisions, double& density);
	double           ComputeApproximateSim(const unsigned& aID, const unsigned& bID);
	pair<unsigned,unsigned> ComputeApproximateSim(const unsigned& aID, const vector<unsigned>& bSignature);

};

class HistogramIndex : public MinHashEncoder
{
public:

	unsigned mHistogramSize;

	typedef uint16_t binKeyTy;
	typedef binKeyTy* indexBinTy;

//	typedef std::tr1::unordered_map<unsigned, indexBinTy> indexSingleTy;
//	typedef google::dense_hash_map<unsigned, indexBinTy> indexSingleTy;
	typedef google::sparse_hash_map<unsigned, indexBinTy> indexSingleTy;
	typedef vector<indexSingleTy> indexTy;

	const binKeyTy MAXBINKEY = std::numeric_limits<binKeyTy>::max();

	HistogramIndex(Parameters* apParameters, Data* apData)
		:MinHashEncoder(apParameters,apData,CLASSIFY)
	{
		mInverseIndex.clear();
		for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k){
			mInverseIndex.push_back(indexSingleTy());
			//mInverseIndex.back().set_empty_key(0);
		}

	}

	multimap<uint, uint> mHistBin2DatasetIdx;
	indexTy mInverseIndex;

	unsigned	GetHistogramSize();
	void  SetHistogramSize(unsigned size);
	void	PrepareIndexDataSets(vector<SeqDataSet>& myFileList);
	void	UpdateInverseIndex(vector<unsigned>& aSignature, unsigned aIndex);
	void  ComputeHistogram(const vector<unsigned>& aSignature, std::valarray<double>& hist, unsigned& emptyBins);
	void	writeBinaryIndex2(ostream &out, const indexTy& index);
	bool	readBinaryIndex2(string filename, indexTy& index);
};

#endif /* MIN_HASH_ENCODER_H */

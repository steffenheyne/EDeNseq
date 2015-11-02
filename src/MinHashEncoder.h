/* -*- mode:c++ -*- */
#ifndef MIN_HASH_ENCODER_H
#define MIN_HASH_ENCODER_H

#include <streambuf>
#include <valarray>
#include <list>
#include <chrono>
#include "Utility.h"
#include "Parameters.h"
#include "Data.h"
#include "sparsehash-2.0.2/sparsehash/sparse_hash_map"
#include "sparsehash-2.0.2/sparsehash/dense_hash_map"
#include "eigen-eigen-3.20/Eigen/Sparse"

#include "MemoryPool.h"

using namespace std;

class MinHashEncoder {


public:

	enum signatureActionE {
		INDEX, INDEX_SIGCACHE, CLASSIFY
	};

	enum groupGraphsByE {
		SEQ_WINDOW, SEQ_NAME, SEQ_FEATURE, NONE
	};

	enum strandTypeE {
		FWD, REV, FR, FR_sep
	};

	typedef strandTypeE			strandTypeT;
	typedef groupGraphsByE		groupGraphsByT;
	typedef signatureActionE	signatureActionT;

	typedef vector<unsigned>				Signature;
	typedef vector<Signature>				SigCacheT;
	typedef std::shared_ptr<SigCacheT>	SigCacheP;

	struct SeqFileS {
		string filename;
		string filename_BED;
		string filename_index;
		InputFileType filetype;
		groupGraphsByT groupGraphsBy;
		bool checkUniqueSeqNames;
		strandTypeT strandType;
		signatureActionT	signatureAction;
		SigCacheP sigCache;
		Data::BEDdataP	dataBED;
		unsigned lastMetaIdx;
		ogzstream* out_results_fh;
	};

	typedef SeqFileS 							SeqFileT;
	typedef std::shared_ptr<SeqFileT> 	SeqFileP;
	typedef vector<SeqFileP> 				SeqFilesT;

	typedef Eigen::SparseVector<unsigned> SVector;

	struct instanceS {
		Signature 	sig;
		string 		name;
		unsigned 	idx;
		unsigned 	pos;
		string 		seq;
		SVector 		svec;
		bool			rc;
		SeqFileP 	seqFile;
	};

	typedef instanceS InstanceT;
	typedef vector<InstanceT> ChunkT;
	typedef std::shared_ptr<ChunkT> ChunkP;

	Parameters* mpParameters;
	Data* 		mpData;

	SeqFileP 	mIndexDataSet;

	unsigned numHashFunctionsFull;
	unsigned sub_hash_range;
	vector<unsigned>  bounds;

protected:

	unsigned numKeys;
	unsigned numFullBins;

	multimap<uint, Data::BEDentryP> mIndexValue2Feature;
	map<string, uint> mFeature2IndexValue;

	mutable std::mutex mutm;
	mutable std::mutex mut1;
	mutable std::mutex mut2;
	mutable std::mutex mut3;

	std::condition_variable cvm;
	std::condition_variable cv1;
	std::condition_variable cv2;
	std::condition_variable cv3;

	threadsafe_queue<SeqFileP> readFile_queue;
	threadsafe_queue<ChunkP> graph_queue;
	threadsafe_queue<ChunkP> sig_queue;
	vector<threadsafe_queue<ChunkP>> index_queue;

	std::atomic_bool done;
	std::atomic_uint files_done;
	std::atomic_uint mSequenceCounter;
	std::atomic_uint mInstanceCounter;
	std::atomic_uint mSignatureCounter;
	std::atomic_uint mSignatureUpdateCounter;


	void					worker_Graph2Signature(int numWorkers);
	void 					finisher();
	void 					generate_feature_vector(const string& seq, SVector& x);
	vector<unsigned>	HashFuncNSPDK(const string& aString, unsigned& aStart, unsigned& aMinRadius, unsigned& aMaxRadius, unsigned& aBitMask);
	void 					worker_readFiles(int numWorkers);
	void 					HashSignatureHelper();

public:


	unsigned 			mHashBitMask;

	MinHashEncoder(Parameters* apParameters, Data* apData);
	virtual	~MinHashEncoder();
	void		Init(Parameters* apParameters, Data* apData);

	virtual void 		UpdateInverseIndex(vector<unsigned>& aSignature, unsigned& aIndex) {};
//	virtual void 		finishUpdate(ChunkP& myData) {};
	virtual void 		finishUpdate(ChunkP& myData, unsigned& min, unsigned& max) {};
	void 					LoadData_Threaded(SeqFilesT& myFiles);
	unsigned				GetLoadedInstances();
	void 					worker_IndexUpdate(unsigned id, unsigned min, unsigned max);

	void					ComputeHashSignature(const SVector& aX, Signature& signaure, Signature* tmpSig);
};

class NeighborhoodIndex : public MinHashEncoder
{
private:

	// inverse index
	typedef unsigned binKeyTy;
	//typedef vector<binKeyTy>	indexBinTy;
	typedef binKeyTy* indexBinTy;
	//	typedef std::tr1::unordered_map<unsigned,indexBinTy> indexSingleTy;
	//	typedef google::sparse_hash_map<unsigned, indexBinTy> indexSingleTy;
	typedef google::dense_hash_map<unsigned, indexBinTy> indexSingleTy;
	typedef vector<indexSingleTy> indexTy;

	const binKeyTy MAXBINKEY = std::numeric_limits<binKeyTy>::max();

public:
	NeighborhoodIndex(Parameters* apParameters, Data* apData)
	:MinHashEncoder(apParameters,apData)
	 {
		mInverseIndex.clear();
		for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k){
			mInverseIndex.push_back(indexSingleTy());
			mInverseIndex.back().set_empty_key(0);
		}
	 }

	indexTy 									mInverseIndex;
	SigCacheP								mMinHashCache;
	vector<vector<unsigned> > 			mNeighborhoodCache;
	vector<umap_uint_int> 				mNeighborhoodCacheExt;
	vector<pair<unsigned,double> >	mNeighborhoodCacheInfo;

	std::tr1::unordered_map<string, unsigned> name2idxMap;
	vector<string>	idx2nameMap;

	void 				  NeighborhoodCacheReset();
	vector<unsigned>& ComputeHashSignature(unsigned aID);
	void 				  UpdateInverseIndex(vector<unsigned>& aSignature, unsigned& aIndex);
	void             ComputeApproximateNeighborhoodCore(const vector<unsigned>& aSignature, umap_uint_int& neighborhood, unsigned& collisions);
	umap_uint_int    ComputeApproximateNeighborhoodExt(const vector<unsigned>& aSignature, unsigned& collisions, double& density);
	vector<unsigned> ComputeApproximateNeighborhood(const vector<unsigned>& aSignature, unsigned& collisions, double& density);

	vector<unsigned> TrimNeighborhood(umap_uint_int& aNeighborhood, unsigned collisions, double& density);
	vector<unsigned> ComputeNeighborhood(const unsigned aID, unsigned& collisions, double& density);
	umap_uint_int 	  ComputeNeighborhoodExt(unsigned aID, unsigned& collisions, double& density);
	double           ComputeApproximateSim(const unsigned& aID, const unsigned& bID);
	pair<unsigned,unsigned> ComputeApproximateSim(const unsigned& aID, const vector<unsigned>& bSignature);

	virtual ~NeighborhoodIndex(){
		for (typename indexTy::const_iterator it = mInverseIndex.begin(); it!= mInverseIndex.end(); it++){
			for (typename indexSingleTy::const_iterator itBin = it->begin(); itBin!=it->end(); itBin++){
				delete[] itBin->second;
			}
		}
	};

};


class HistogramIndex : public MinHashEncoder
{

public:

	const unsigned INDEX_FORMAT_VERSION = 2;

	typedef uint16_t binKeyTy;
	typedef binKeyTy* indexBinTy;
	const binKeyTy MAXBINKEY = std::numeric_limits<binKeyTy>::max();


	struct hashFunc {
		size_t operator()(unsigned a) const {
			return static_cast<size_t>(a);
		}

		bool operator()(unsigned a, unsigned b) const {
			return a == b;
		}
	};

	// we use one of these hash maps, the interface is quite similar to each other
	// should test space an speed
	//typedef std::tr1::unordered_map<unsigned, indexBinTy> indexSingleTy;
	typedef google::dense_hash_map<unsigned, indexBinTy, hashFunc,hashFunc> indexSingleTy;
	//typedef google::sparse_hash_map<unsigned, indexBinTy, hashFunc, hashFunc> indexSingleTy;
	//typedef google::sparse_hash_map<unsigned, indexBinTy> indexSingleTy;

	// the index
	typedef vector<indexSingleTy> indexTy;
	// used for memory pool
	typedef binKeyTy newIndexBin[2];

	typedef valarray<double> histogramT;

	const static unsigned mMemPool_BlockSize = 524288; // (2^19)

	// a vector of memory pools, we use for each single Index a different memory pool
	vector<MemoryPool<newIndexBin,mMemPool_BlockSize>*>  mMemPool;

	binKeyTy mHistogramSize;
	indexTy mInverseIndex;

	/////////////////////////////
	// member functions
	////////////////////////////

	// constructor
	HistogramIndex(Parameters* apParameters, Data* apData)
	:MinHashEncoder(apParameters,apData)
	{
		mInverseIndex.resize(mpParameters->mNumHashFunctions, indexSingleTy(2^22));
		for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k){
			mInverseIndex[k].max_load_factor(0.999);
			mInverseIndex[k].set_empty_key(0);
			mMemPool.push_back(new MemoryPool<newIndexBin,mMemPool_BlockSize>());
		}
	}

	binKeyTy	GetHistogramSize();
	void		SetHistogramSize(binKeyTy size);
	void		UpdateInverseIndex(const vector<unsigned>& aSignature, const unsigned& aIndex);
	void		UpdateInverseIndex(const vector<unsigned>& aSignature, const unsigned& aIndex, unsigned& min, unsigned& max);
	void  	ComputeHistogram(const vector<unsigned>& aSignature, std::valarray<double>& hist, unsigned& emptyBins);
	void		writeBinaryIndex2(ostream &out, const indexTy& index);
	bool		readBinaryIndex2(string filename, indexTy& index);

	// destructor
	virtual ~HistogramIndex(){
		uint k=0;
		for (typename indexTy::const_iterator it = mInverseIndex.begin(); it!= mInverseIndex.end(); it++){
			for (typename indexSingleTy::const_iterator itBin = it->begin(); itBin!=it->end(); itBin++){
				if ( itBin->second[0]>1){
					delete[] itBin->second;
				} else {
					mMemPool[k]->deleteElement(reinterpret_cast<newIndexBin(*)>(itBin->second));
				}
			}
			k++;
		}
	};
};


//###########################################################################
//
// NOT USED YET -- UNDER DEVELOPMENT ----------------------------------------
//
template<typename T>
class MinHashIndex {

public:
	typedef	T 				ValueBinTy;
	typedef	ValueBinTy*	ValueTy;
	typedef	unsigned 	KeyTy;

	struct hashFunc {
		size_t operator()(KeyTy a) const {
			return static_cast<size_t>(a);
		}

		bool operator()(KeyTy a, KeyTy b) const {
			return a == b;
		}
	};
	//	typedef std::tr1::unordered_map<unsigned, indexBinTy,hashFunc> indexSingleTy;
	//	typedef google::dense_hash_map<unsigned, indexBinTy, hashFunc> indexSingleTy;
	//	typedef google::sparse_hash_map<unsigned, indexBinTy> indexSingleTy;

	typedef google::sparse_hash_map<KeyTy, ValueTy> indexSingleTy;

	typedef vector<indexSingleTy> indexTy;
	const ValueTy MAXBINKEY = std::numeric_limits<ValueTy>::max();

	int 		mIndexSize;
	indexTy mInverseIndex;

	MinHashIndex(){};
	MinHashIndex(int l):mIndexSize(l){
		mInverseIndex.clear();
				for (unsigned i = 0; i < mIndexSize; ++i){
					mInverseIndex.push_back(indexSingleTy());
					mInverseIndex.back().set_empty_key(0);
				}
	};

};

#endif /* MIN_HASH_ENCODER_H */

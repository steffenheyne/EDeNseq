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
#include "sparsehash-2.0.2/sparsehash/dense_hash_map"
//#include "sparsehash-2.0.2/sparsehash/dense_hash_set"
//#include "sparsehash-2.0.2/sparsehash/sparse_hash_set"
#include <Eigen/Sparse>
#include <chrono>

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
//	typedef google::dense_hash_set<unsigned> SVectorMap;


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

	struct resultS{
		string output_line;
		unsigned numInstances;
	};

	typedef resultS ResultT;
	typedef vector<ResultT> ResultChunkT;
	typedef std::shared_ptr<ResultChunkT> ResultChunkP;

	Parameters* mpParameters;
	Data* 		mpData;

	SeqFileP 	mIndexDataSet;



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

	std::atomic_bool done;
	std::atomic_uint files_done;
	std::atomic_uint mSequenceCounter;
	std::atomic_uint mInstanceCounter;
	std::atomic_uint mSignatureCounter;



	void					worker_Graph2Signature(int numWorkers);
	void 					finisher();
	void 					generate_feature_vector(const string& seq, SVector& x);
	vector<unsigned>	HashFuncNSPDK(const string& aString, unsigned aStart, unsigned aMaxRadius, unsigned aBitMask);

public:


	unsigned 			mHashBitMask;
	void 					worker_readFiles(int numWorkers);
	MinHashEncoder(Parameters* apParameters, Data* apData);
	virtual	~MinHashEncoder();
	void		Init(Parameters* apParameters, Data* apData);

	virtual void 		UpdateInverseIndex(vector<unsigned>& aSignature, unsigned aIndex) {};
	virtual void 		finishUpdate(ChunkP& myData) {};
	void 					LoadData_Threaded(SeqFilesT& myFiles);
	unsigned				GetLoadedInstances();

	void					ComputeHashSignature(const SVector& aX, Signature& signaure);
};

class NeighborhoodIndex : public MinHashEncoder
{
private:

	// inverse index
	typedef unsigned binKeyTy;
	typedef vector<binKeyTy>	indexBinTy;
	typedef std::tr1::unordered_map<unsigned,indexBinTy> indexSingleTy;
	typedef vector<indexSingleTy> indexTy;

	const binKeyTy MAXBINKEY = std::numeric_limits<binKeyTy>::max();

public:
	NeighborhoodIndex(Parameters* apParameters, Data* apData)
		:MinHashEncoder(apParameters,apData)
	{
		mInverseIndex.clear();
		for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k){
			mInverseIndex.push_back(indexSingleTy());
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
	void 				  UpdateInverseIndex(vector<unsigned>& aSignature, unsigned aIndex);
	void             ComputeApproximateNeighborhoodCore(const vector<unsigned>& aSignature, umap_uint_int& neighborhood, unsigned& collisions);
	umap_uint_int    ComputeApproximateNeighborhoodExt(const vector<unsigned>& aSignature, unsigned& collisions, double& density);
	vector<unsigned> ComputeApproximateNeighborhood(const vector<unsigned>& aSignature, unsigned& collisions, double& density);

	vector<unsigned> TrimNeighborhood(umap_uint_int& aNeighborhood, unsigned collisions, double& density);
	vector<unsigned> ComputeNeighborhood(const unsigned aID, unsigned& collisions, double& density);
	umap_uint_int 	  ComputeNeighborhoodExt(unsigned aID, unsigned& collisions, double& density);
	double           ComputeApproximateSim(const unsigned& aID, const unsigned& bID);
	pair<unsigned,unsigned> ComputeApproximateSim(const unsigned& aID, const vector<unsigned>& bSignature);

};

class HistogramIndex : public MinHashEncoder
{

public:
	typedef uint16_t binKeyTy;
	typedef binKeyTy* indexBinTy;
	//typedef vector<binKeyTy> indexBinTy;

	struct hashFunc {
	 size_t operator()(unsigned a) const {
	    return static_cast<size_t>(a);
	  }

	 bool operator()(unsigned a, unsigned b) const {
	     return a == b;
	   }
	};
	//	typedef std::tr1::unordered_map<unsigned, indexBinTy,hashFunc> indexSingleTy;
	//	typedef google::dense_hash_map<unsigned, indexBinTy, hashFunc> indexSingleTy;
		typedef google::sparse_hash_map<unsigned, indexBinTy, hashFunc, hashFunc> indexSingleTy;
	//typedef google::sparse_hash_map<unsigned, indexBinTy> indexSingleTy;

	typedef vector<indexSingleTy> indexTy;
	const binKeyTy MAXBINKEY = std::numeric_limits<binKeyTy>::max();

	typedef valarray<double> histogramT;

	binKeyTy mHistogramSize;
	indexTy mInverseIndex;

	HistogramIndex(Parameters* apParameters, Data* apData)
		:MinHashEncoder(apParameters,apData)
	{
		mInverseIndex.resize(mpParameters->mNumHashFunctions, indexSingleTy(500000000));
		for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k){
			//mInverseIndex.push_back(indexSingleTy(500000000));
//			mInverseIndex[k].min_load_factor(0.1);
			mInverseIndex[k].max_load_factor(0.9);
	//		mInverseIndex.back().resize(500000000);
	//		mInverseIndex.back().min_load_factor(1000.0);
	//		mInverseIndex.back().set_empty_key(0);
		}

	}

	binKeyTy	GetHistogramSize();
	void  SetHistogramSize(binKeyTy size);
	void	UpdateInverseIndex(const vector<unsigned>& aSignature, const unsigned& aIndex);
	void  ComputeHistogram(const vector<unsigned>& aSignature, std::valarray<double>& hist, unsigned& emptyBins);
	void	writeBinaryIndex2(ostream &out, const indexTy& index);
	bool	readBinaryIndex2(string filename, indexTy& index);
};

#endif /* MIN_HASH_ENCODER_H */

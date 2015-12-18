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
	const unsigned MAXUNSIGNED = (2 << 31)-1 ;

	enum signatureActionE {
		INDEX, INDEX_SIGCACHE, CLASSIFY
	};

	enum groupGraphsByE {
		SEQ_WINDOW, SEQ_NAME, SEQ_FEATURE, SEQ_NUM, NONE
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

	//typedef Eigen::SparseVector<unsigned> SVector;
	typedef vector<unsigned> SVector;


	struct instanceS {
		Signature 	sig;
		string 		name;
		unsigned 	idx;
		unsigned 	pos;
		string 		seq;
		SVector 		svec;
		vector<vector<unsigned> > minHashes;
		bool			rc;
		SeqFileP 	seqFile;
	};

	typedef instanceS InstanceT;
	typedef vector<InstanceT> ChunkT;
	typedef std::shared_ptr<ChunkT> ChunkP;

	Parameters* mpParameters;
	Data* 		mpData;

	SeqFileP 	mIndexDataSet;

	unsigned 			numHashFunctionsFull;
	unsigned 			sub_hash_range;
	vector<unsigned>  mBounds;


	unsigned numKeys;
	unsigned numFullBins;
	bool		mUseSlidingWindowMinHash;

	multimap<uint, Data::BEDentryP> mIndexValue2Feature;
	map<string, uint> mFeature2IndexValue;

	mutable std::mutex mutm;
	mutable std::mutex mut1;
	mutable std::mutex mut2;
	mutable std::mutex mut3;
	mutable std::mutex mut4;
	mutable std::mutex mut5;

	std::condition_variable cvm;
	std::condition_variable cv1;
	std::condition_variable cv2;
	std::condition_variable cv3;
	std::condition_variable cv4;

	threadsafe_queue<SeqFileP> readFile_queue;
	vector<threadsafe_queue<ChunkP>> graph_queue;
	vector<threadsafe_queue<ChunkP>> index_queue;

	std::atomic_bool done;
	std::atomic_uint files_done;
	std::atomic_uint mSequenceCounter;
	std::atomic_uint mInstanceCounter;
	std::atomic_uint mInstanceProcCounter;
	std::atomic_uint mSignatureCounter;
	std::atomic_uint mSignatureUpdateCounter;

	unsigned 			mHashBitMask;
	unsigned 			mHashBitMask_feature;
	unsigned 			mHashBitMask_shingle;

				MinHashEncoder(Parameters* apParameters, Data* apData);
	virtual		~MinHashEncoder();

	void					Init(Parameters* apParameters, Data* apData);

	void 					LoadData_Threaded(SeqFilesT& myFiles);
	void 					worker_readFiles(unsigned numWorkers, unsigned chunkSizeFactor);
	void					worker_Seq2Signature_SlidingWin(int numWorkers,unsigned id);
	void              worker_Seq2Signature_SingleWin(int numWorkers, unsigned id);
	void 					finisher_IndexUpdate(unsigned id, unsigned min, unsigned max);
	virtual void			finishUpdate(ChunkP& myData, unsigned& min, unsigned& max) {};

	vector<unsigned>		HashFuncNSPDK(const string& aString, unsigned& aStart, unsigned& aMinRadius, unsigned& aMaxRadius, unsigned aBitMask);
	void 					HashSignatureHelper();
	void 					generate_feature_vector(const string& seq, SVector& x);
	void					ComputeHashSignature(const SVector& aX, Signature& signaure, Signature* tmpSig);

	vector<unsigned> 				iterated_hash(string& kmer, unsigned minRadius);
	void								running_hash(vector<vector<unsigned>>& paired_kmer_hashes_array, string& seq, unsigned& minRadius, unsigned& maxRadius, unsigned& minDist, unsigned& maxDist);
	void								sliding_window_minimum(vector<vector<unsigned>>& array, unsigned winsize, unsigned& maxFeatLen);
	void								sliding_window_minhash(vector<vector<unsigned>>& res, string& seq, unsigned& minRadius, unsigned maxRadius, unsigned minDistance, unsigned maxDistance, unsigned winsize, unsigned step);

	virtual void 			UpdateInverseIndex(vector<unsigned>& aSignature, unsigned& aIndex) {};

	unsigned				GetLoadedInstances();
};

class NeighborhoodIndex : public MinHashEncoder
{
private:

	// inverse index
	typedef unsigned binKeyTy;
	//typedef vector<binKeyTy>	indexBinTy;
	typedef binKeyTy* indexBinTy;
	//	typedef std::tr1::unordered_map<unsigned,indexBinTy> indexSingleTy;
	typedef google::sparse_hash_map<unsigned, indexBinTy> indexSingleTy;
	//typedef google::dense_hash_map<unsigned, indexBinTy> indexSingleTy;
	typedef vector<indexSingleTy> indexTy;

	const binKeyTy MAXBINKEY = std::numeric_limits<binKeyTy>::max();

public:
	NeighborhoodIndex(Parameters* apParameters, Data* apData)
	:MinHashEncoder(apParameters,apData)
	{
		mInverseIndex.clear();
		for (unsigned k = 0; k < mpParameters->mNumHashFunctions; ++k){
			mInverseIndex.push_back(indexSingleTy());
			//		mInverseIndex.back().set_empty_key(0);
		}
	}

	indexTy 									mInverseIndex;
	SigCacheP								mMinHashCache;
	vector<vector<unsigned> > 			mNeighborhoodCache;
	vector<umap_uint_int> 				mNeighborhoodCacheExt;
	vector<pair<unsigned,double> >	mNeighborhoodCacheInfo;

	std::tr1::unordered_map<string, unsigned> name2idxMap;
	vector<string>	idx2nameMap;

	void 			  		NeighborhoodCacheReset();
	vector<unsigned>& 		ComputeHashSignature(unsigned aID);
	void 			 		UpdateInverseIndex(vector<unsigned>& aSignature, unsigned& aIndex);
	void             		ComputeApproximateNeighborhoodCore(const vector<unsigned>& aSignature, umap_uint_int& neighborhood, unsigned& collisions);
	umap_uint_int    		ComputeApproximateNeighborhoodExt(const vector<unsigned>& aSignature, unsigned& collisions, double& density);
	vector<unsigned> 		ComputeApproximateNeighborhood(const vector<unsigned>& aSignature, unsigned& collisions, double& density);

	vector<unsigned> 		TrimNeighborhood(umap_uint_int& aNeighborhood, unsigned collisions, double& density);
	vector<unsigned> 		ComputeNeighborhood(const unsigned aID, unsigned& collisions, double& density);
	umap_uint_int 	  		ComputeNeighborhoodExt(unsigned aID, unsigned& collisions, double& density);
	double           		ComputeApproximateSim(const unsigned& aID, const unsigned& bID);
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
	//typedef google::dense_hash_map<unsigned, indexBinTy, hashFunc, hashFunc> indexSingleTy;
	//typedef google::sparse_hash_map<unsigned, indexBinTy, hashFunc, hashFunc,std::allocator<indexBinTy>> indexSingleTy;
	typedef google::sparse_hash_map<unsigned, indexBinTy,hashFunc, hashFunc> indexSingleTy;

	// the index
	typedef vector<indexSingleTy> indexTy;

	// used for memory pool
	typedef binKeyTy newIndexBin_2[2];
	typedef binKeyTy newIndexBin_3[3];
	typedef binKeyTy newIndexBin_4[4];
	typedef binKeyTy newIndexBin_5[5];
	typedef binKeyTy newIndexBin_6[6];
	typedef binKeyTy newIndexBin_7[7];
	typedef binKeyTy newIndexBin_8[8];
	typedef binKeyTy newIndexBin_9[9];
	typedef binKeyTy newIndexBin_10[10];

	typedef valarray<double> histogramT;

	//const static unsigned mMemPool_BlockSize = 720*4096; // (2*3*4*5*6*4096)
	//const static unsigned mMemPool_BlockSize = 491520; // kgv(2*3*4*5*6)=60*4096*2
	//const static unsigned mMemPool_BlockSize = 1966080; // kgv(2*3*4*5*6)=60*4096*2
	//const static unsigned mMemPool_BlockSize = 1290240 ; // kgv(2*3*4*5*6)=60*7*8*9*10=2520,4096=1290240
	const static unsigned mMemPool_BlockSize = 409600; // (2*3*4*5*6*4096)
	//const static unsigned mMemPool_maxEnt = 2;

	// a vector of memory pools,
	// in case of building a new index, we use for each index updater thread a different memory pool,
	// in case of reading index from file, we use only a single Memory pool: mMemPool[0]
	//vector<MemoryPool<newIndexBin,mMemPool_BlockSize>*>  mMemPool;
	vector<MemoryPool<newIndexBin_2,mMemPool_BlockSize>*>  mMemPool_2;
	vector<MemoryPool<newIndexBin_3,mMemPool_BlockSize>*>  mMemPool_3;
	vector<MemoryPool<newIndexBin_4,mMemPool_BlockSize>*>  mMemPool_4;
	vector<MemoryPool<newIndexBin_5,mMemPool_BlockSize>*>  mMemPool_5;
	vector<MemoryPool<newIndexBin_6,mMemPool_BlockSize>*>  mMemPool_6;
	vector<MemoryPool<newIndexBin_7,mMemPool_BlockSize>*>  mMemPool_7;
	vector<MemoryPool<newIndexBin_8,mMemPool_BlockSize>*>  mMemPool_8;
	vector<MemoryPool<newIndexBin_9,mMemPool_BlockSize>*>  mMemPool_9;
	vector<MemoryPool<newIndexBin_10,mMemPool_BlockSize>*>  mMemPool_10;

	binKeyTy mHistogramSize;
	indexTy mInverseIndex;

	/////////////////////////////
	// member functions
	////////////////////////////

	// constructor
	HistogramIndex(Parameters* apParameters, Data* apData)
	:MinHashEncoder(apParameters,apData) { };

	void		InitInverseIndex();
	binKeyTy	GetHistogramSize();
	void		SetHistogramSize(binKeyTy size);
	void		UpdateInverseIndex(const vector<unsigned>& aSignature, const unsigned& aIndex);
	void		UpdateInverseIndex(const vector<unsigned>& aSignature, const unsigned& aIndex, unsigned& min, unsigned& max);
	void 		UpdateInverseIndex(const unsigned& key, const unsigned& aIndex, unsigned& k);
	//void		ComputeHistogram(const vector<unsigned>& aSignature, std::valarray<double>& hist, unsigned& emptyBins);
	void 		ComputeHistogram(const vector<vector<unsigned>>& aSigArray, std::valarray<double>& hist, vector<unsigned>& emptyBins);
	void		writeBinaryIndex2(ostream &out, const indexTy& index);
	bool		readBinaryIndex2(string filename, indexTy& index);

	// destructor
	virtual ~HistogramIndex(){

		uint k=0;
		for (typename indexTy::const_iterator it = mInverseIndex.begin(); it!= mInverseIndex.end(); it++){
			for (typename indexSingleTy::const_iterator itBin = it->begin(); itBin!=it->end(); itBin++){
				switch (itBin->second[0]){
				case 1:
					mMemPool_2[k]->deleteElement(reinterpret_cast<newIndexBin_2(*)>(itBin->second));
					break;
				case 2:
					mMemPool_3[k]->deleteElement(reinterpret_cast<newIndexBin_3(*)>(itBin->second));
					break;
				case 3:
					mMemPool_4[k]->deleteElement(reinterpret_cast<newIndexBin_4(*)>(itBin->second));
					break;
				case 4:
					mMemPool_5[k]->deleteElement(reinterpret_cast<newIndexBin_5(*)>(itBin->second));
					break;
				case 5:
					mMemPool_6[k]->deleteElement(reinterpret_cast<newIndexBin_6(*)>(itBin->second));
					break;
				case 6:
					mMemPool_7[k]->deleteElement(reinterpret_cast<newIndexBin_7(*)>(itBin->second));
					break;
				case 7:
					mMemPool_8[k]->deleteElement(reinterpret_cast<newIndexBin_8(*)>(itBin->second));
					break;
				case 8:
					mMemPool_9[k]->deleteElement(reinterpret_cast<newIndexBin_9(*)>(itBin->second));
					break;
				case 9:
					mMemPool_10[k]->deleteElement(reinterpret_cast<newIndexBin_10(*)>(itBin->second));
					break;
				default:
					delete[] itBin->second;
					break;
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

template <typename T>
class IndexMemoryPool {

	typedef T newIndexBin_2[2];
	typedef T newIndexBin_3[3];
	typedef T newIndexBin_4[4];
	typedef T newIndexBin_5[5];
	typedef T newIndexBin_6[6];
	typedef T newIndexBin_7[7];
	typedef T newIndexBin_8[8];
	typedef T newIndexBin_9[9];
	typedef T newIndexBin_10[10];

	const static unsigned mMemPool_BlockSize = 1512000;

	vector<MemoryPool<newIndexBin_2,mMemPool_BlockSize>*>  mMemPool_2;
	vector<MemoryPool<newIndexBin_3,mMemPool_BlockSize>*>  mMemPool_3;
	vector<MemoryPool<newIndexBin_4,mMemPool_BlockSize>*>  mMemPool_4;
	vector<MemoryPool<newIndexBin_5,mMemPool_BlockSize>*>  mMemPool_5;
	vector<MemoryPool<newIndexBin_6,mMemPool_BlockSize>*>  mMemPool_6;
	vector<MemoryPool<newIndexBin_7,mMemPool_BlockSize>*>  mMemPool_7;
	vector<MemoryPool<newIndexBin_8,mMemPool_BlockSize>*>  mMemPool_8;
	vector<MemoryPool<newIndexBin_9,mMemPool_BlockSize>*>  mMemPool_9;
	vector<MemoryPool<newIndexBin_10,mMemPool_BlockSize>*>  mMemPool_10;


}
;

class ArrayMemoryPool {

public:

	const int length;
	int BlockSize ;


	typedef uint16_t	element_type;

	typedef element_type    value_type;
	typedef element_type*   pointer;
	typedef element_type&        reference;
	typedef const element_type*  const_pointer;
	typedef const element_type&  const_reference;
	typedef size_t          size_type;
	typedef ptrdiff_t       difference_type;

	pointer arrayP;

private:
	union Slot_ {
		value_type element;
		Slot_* next;
	};

	typedef char* data_pointer_;
	typedef Slot_ slot_type_;
	typedef Slot_* slot_pointer_;

	slot_pointer_ currentBlock_;
	slot_pointer_ currentSlot_;
	slot_pointer_ lastSlot_;
	slot_pointer_ freeSlots_;

	/*  void test2(){
	  void* my;
  	 void* my2;
	  my = *reinterpret_cast<value_type(*)[length]>(my2);
  }*/

public:

	ArrayMemoryPool(int size, unsigned blocksize): length(size),BlockSize(blocksize),arrayP(new element_type[size])
	{
		currentBlock_ = nullptr;
		currentSlot_ = nullptr;
		lastSlot_ = nullptr;
		freeSlots_ = nullptr;
	};


	inline size_type padPointer(data_pointer_ p, size_type align) const
	{
		uintptr_t result = reinterpret_cast<uintptr_t>(p);
		return ((align - result) % align);
	};


	void allocateBlock()
	{
		// Allocate space for the new block and store a pointer to the previous one
		data_pointer_ newBlock = reinterpret_cast<data_pointer_>(operator new(BlockSize));
		reinterpret_cast<slot_pointer_>(newBlock)->next = currentBlock_;
		currentBlock_ = reinterpret_cast<slot_pointer_>(newBlock);

		// Pad block body to staisfy the alignment requirements for elements
		data_pointer_ body = newBlock + sizeof(slot_pointer_);
		size_type bodyPadding = padPointer(body, alignof(slot_type_));
		currentSlot_ = reinterpret_cast<slot_pointer_>(body + bodyPadding);
		lastSlot_ = reinterpret_cast<slot_pointer_>
		(newBlock + BlockSize - sizeof(slot_type_) + 1);
	};


	inline pointer allocate()
	{
		if (freeSlots_ != nullptr) {
			pointer result = reinterpret_cast<pointer>(freeSlots_);
			freeSlots_ = freeSlots_->next;
			return result;
		}
		else {
			if (currentSlot_ >= lastSlot_)
				allocateBlock();
			return reinterpret_cast<pointer>(currentSlot_++);
		}
	};



	inline void deallocate(pointer p)
	{
		if (p != nullptr) {
			reinterpret_cast<slot_pointer_>(p)->next = freeSlots_;
			freeSlots_ = reinterpret_cast<slot_pointer_>(p);
		}
	};

	inline void construct(value_type* p)
	{
		new (p) value_type();
	};

	inline pointer newElement()
	{
		pointer result = allocate();
		construct(result);
		return result;
	};



	inline void deleteElement(pointer p)
	{
		if (p != nullptr) {
			//p->~value_type();
			deallocate(p);
		}
	};



};

#endif /* MIN_HASH_ENCODER_H */

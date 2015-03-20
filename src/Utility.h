#ifndef UTILITY_H
#define UTILITY_H

#include <sys/types.h>
#include <sys/stat.h>
#include <cstring>
#include <errno.h>
#include <cstdio>
#include <cstdlib>
#include <ctype.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cassert>
#include <vector>
#include <list>
#include <stdexcept>
#include <map>
#include <set>
#include <memory>
#include <tr1/unordered_map>
#include <queue>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <ctime>
#include <iomanip>
#include <limits>
#include <thread>
#include <condition_variable>
#include <mutex>
#include <atomic>

#ifdef USEMULTITHREAD
#include <omp.h>
#endif

#include "vectors.h"

using namespace std;

typedef std::tr1::unordered_map<unsigned, int> umap_uint_int;
typedef std::tr1::unordered_map<unsigned, vector<unsigned> > umap_uint_vec_uint;

//typedef std::tr1::unordered_map<pair<unsigned, unsigned>, int> umap_pair_uint_uint_int;
//typedef std::tr1::unordered_map<pair<unsigned, int>, vector<unsigned> > umap_pair_uint_int_vec_uint;

//------------------------------------------------------------------------------------------------------------------------
const string SEP = "-----------------------------------------------------------------";
const string TAB = "    ";

//------------------------------------------------------------------------------------------------------------------------
const unsigned MAXUNSIGNED = 2 << 30;

///Returns a random number uniformly distributed between 0 and 1
inline double random01() {
	return (double) rand() / (double) RAND_MAX;
}

///Returns a random number approximately normally distributed with mean 0 and std 1
inline double randomn() {
	return  -log(1 / random01() - 1) / 1.702; //Rao, K. R., Boiroju, N. K., & Reddy, M. K. (2011). GENERATION OF STANDARD NORMAL RANDOM VARIABLES. Indian Journal of Scientific Research, 2(4).
}

///Returns a random integer uniformly distributed between 0 and the argument aMax-1
inline unsigned randomUnsigned(unsigned aMax) {
	return (unsigned) (rand() % aMax);
}

template<class InType>
void PermuteVector(vector<InType> & t) {
	for (unsigned i = 0; i < t.size(); ++i) {
		unsigned j = randomUnsigned(t.size());
		swap(t[i], t[j]);
	}
}

///Implements a type converter via streams
template<class OutType, class InType>
OutType stream_cast(const InType & t) {
	stringstream ss;
	ss << t; // first insert value to stream
	OutType result; // value will be converted to OutType
	ss >> result; // write value to result
	return result;
}

void MakeShuffledDataIndicesList(vector<unsigned>& oDataIdList, unsigned aSize);

///Return an integer hash value for a given input integer in a given domain range
inline int IntHashSimple(int key, int aModulo) {
	key = ~key + (key << 15); // key = (key << 15) - key - 1;
	key = key ^ (key >> 12);
	key = key + (key << 2);
	key = key ^ (key >> 4);
	key = key * 2057; // key = (key + (key << 3)) + (key << 11);
	key = key ^ (key >> 16);
	return key % aModulo;
}

inline int IntHashPair(int key1, int key2, int aModulo = 2147483647) {
	const double A = sqrt(2) - 1;
	int key = key1 * (key2 + 1) * A;
	return IntHashSimple(key, aModulo);
}

///Return an integer hash value for a given input integer in a given domain range given an additional seed to select the random hash function
inline int IntHash(int key, int aModulo, unsigned aSeed) {
	const double A = sqrt(2) - 1;
	return IntHashSimple(key * (aSeed + 1) * A, aModulo);
}

unsigned RSHash(const string& str);
unsigned RSHash(const vector<unsigned>& aV);
unsigned APHash(const string& str);
unsigned APHash(const vector<unsigned>& aV);
unsigned HashFunc(const string& str, unsigned aBitMask = 2147483647);
unsigned HashFunc(const vector<unsigned>& aList, unsigned aBitMask = 2147483647);

inline vector<unsigned> ComputeMinHashSignature(SVector& aX, unsigned aNumHashFunctions) {
	const unsigned MAXUNSIGNED = 2147483647;
	//prepare a vector containing the k min values
	vector<unsigned> signature(aNumHashFunctions, MAXUNSIGNED);
	for (SparseVector<double>::InnerIterator it(aX); it; ++it) {
		//extract only the feature id (i.e. ignore the actual feature value)
		unsigned hash_id = it.index();
		if (hash_id == 0) {
			hash_id = 1; //force collision between feature 0 and 1 to avoid features with ID=0
		}
		for (unsigned k = 0; k < aNumHashFunctions; ++k) { //for all k hashes
			unsigned new_hash = IntHash(hash_id, MAXUNSIGNED, k); //rehash the feature id with a procedure that is aware of the index k
			if (new_hash < signature[k])
				signature[k] = new_hash; //keep the minimum value only
		}
	}
	return signature;
}

//------------------------------------------------------------------------------------------------------------------------
///Returns the time between the creation and destruction of an object of TimerClass
class TimerClass {
public:
	TimerClass();
	~TimerClass() {
	}
	void Output();
private:
	std::time_t mStartSec;
	std::clock_t mStart;
};

//------------------------------------------------------------------------------------------------------------------------
///Plots the increase of an internal counter (which has to be
///explicitly increased by calling the Count member). Returns also
///the time elapsed between the creation and the destruction of an
///object of ProgressBar.
class ProgressBar {
public:
	ProgressBar(unsigned aStep = 100);
	~ProgressBar();
	void Begin();
	void Count();
	void Count(int co);
	unsigned End();
private:
	unsigned mStep;
	atomic_int mCounter;
	TimerClass mTimer;
};

//------------------------------------------------------------------------------------------------------------------------
///Implements safe access policies to a vector container and offers
///members to compute various statistical estimators.
class VectorClass {
	friend ostream& operator<<(ostream& out, const VectorClass& aV);
public:
	VectorClass();
	VectorClass(unsigned aSize);
	void operator=(const VectorClass& aVector);
	VectorClass(const VectorClass& aVector);
	VectorClass(const vector<double>& aVector);
	VectorClass(const vector<int>& aVector);
	void Init(unsigned aSize);
	void Import(const string& aFileName);
	void Clear();
	unsigned Size() const;
	ostream& Output(ostream& out) const;
	void PushBack(double aValue);
	double& operator[](unsigned i);
	double operator[](unsigned i) const;
	double Prod() const;
	double Sum() const;
	double Mean() const;
	double StandardDeviation() const;
	double Order(double aOrder) const;
	double Median() const;
	double MedianAbsoluteDifference() const;
	double Min() const;
	double Max() const;
	VectorClass RemoveNulls();
	ostream& OutputStatistics(ostream& out);
protected:
	vector<double> mV;
};

//------------------------------------------------------------------------------------------------------------------------
template<class ForwardIterator, class T>
void iota(ForwardIterator first, ForwardIterator last, T val) {
	while (first != last) {
		*first = val;
		++first;
		++val;
	}
}

//------------------------------------------------------------------------------------------------------------------------
class OutputManager{
public:
	ofstream mOut;
protected:
	string mFileName;
	string mDirectoryPath;
public:
	OutputManager(string aFileName,string aDirectoryPath);
	string GetFullPathFileName();
};


//-------------------------------------------------------------------------------------------------------------------------

template<typename T>
class threadsafe_queue
{
private:
	mutable std::mutex mut;
	std::queue<T> data_queue;
	std::condition_variable data_cond;
public:
	threadsafe_queue()
	{}
	threadsafe_queue(threadsafe_queue const& other)
	{
		std::lock_guard<std::mutex> lk(other.mut);
		data_queue=other.data_queue;
	}
	void push(T new_value)
	{
		std::lock_guard<std::mutex> lk(mut);
		data_queue.push(new_value);
		data_cond.notify_one();
	}
	void wait_and_pop(T& value)
	{
		std::unique_lock<std::mutex> lk(mut);
		data_cond.wait(lk,[this]{return !data_queue.empty();});
		value=data_queue.front();
		data_queue.pop();
	}
	std::shared_ptr<T> wait_and_pop()
																																										{
		std::unique_lock<std::mutex> lk(mut);
		data_cond.wait(lk,[this]{return !data_queue.empty();});
		std::shared_ptr<T> res(std::make_shared<T>(data_queue.front()));
		data_queue.pop();
		return res;
																																										}

	bool try_pop(T& value)
	{
		std::lock_guard<std::mutex> lk(mut);
		if(data_queue.empty())
			return false;
		value=data_queue.front();
		data_queue.pop();
		return true;
	}
	std::shared_ptr<T> try_pop(){
		std::lock_guard<std::mutex> lk(mut);
		if(data_queue.empty())
			return std::shared_ptr<T>();
		std::shared_ptr<T> res(std::make_shared<T>(data_queue.front()));
		data_queue.pop();
		return res;
}

	bool empty() const
	{
		std::lock_guard<std::mutex> lk(mut);
		return data_queue.empty();
	}

	int size() const
	{
		std::lock_guard<std::mutex> lk(mut);
		return data_queue.size();
	}

};

class join_threads
{
	vector<thread>& threads;
public:
	explicit join_threads(vector<thread>& threads_):
	threads(threads_)
	{}
	~join_threads()
	{
		for(unsigned long i=0;i<threads.size();++i)
		{
			if(threads[i].joinable())
				threads[i].join();
		}
		cout << threads.size() << " threads joined" << endl;
	}
};



#endif

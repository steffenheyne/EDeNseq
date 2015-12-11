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
#include <chrono>
#include <deque>

#ifdef USEMULTITHREAD
#include <omp.h>
#endif

//#include "vectors.h"

using namespace std;

typedef std::tr1::unordered_map<unsigned, int> umap_uint_int;
typedef std::tr1::unordered_map<unsigned, vector<unsigned> > umap_uint_vec_uint;
typedef std::tr1::unordered_map<unsigned, vector<int> > umap_uint_vec_int;

//typedef std::unordered_map<unsigned, int> umap_uint_int;
//typedef std::unordered_map<unsigned, vector<unsigned> > umap_uint_vec_uint;
//typedef std::unordered_map<unsigned, vector<int> > umap_uint_vec_int;


//typedef std::tr1::unordered_map<pair<unsigned, unsigned>, int> umap_pair_uint_uint_int;
//typedef std::tr1::unordered_map<pair<unsigned, int>, vector<unsigned> > umap_pair_uint_int_vec_uint;

//------------------------------------------------------------------------------------------------------------------------
const string SEP = "-----------------------------------------------------------------";
const string TAB = "    ";

//------------------------------------------------------------------------------------------------------------------------

///Returns a random number uniformly distributed between 0 and 1
//inline double random01() {
//	return (double) rand() / (double) RAND_MAX;
//}
//
/////Returns a random number approximately normally distributed with mean 0 and std 1
//inline double randomn() {
//	return  -log(1 / random01() - 1) / 1.702; //Rao, K. R., Boiroju, N. K., & Reddy, M. K. (2011). GENERATION OF STANDARD NORMAL RANDOM VARIABLES. Indian Journal of Scientific Research, 2(4).
//}
//
/////Returns a random integer uniformly distributed between 0 and the argument aMax-1
//inline unsigned randomUnsigned(unsigned aMax) {
//	return (unsigned) (rand() % aMax);
//}
//
//template<class InType>
//void PermuteVector(vector<InType> & t) {
//	for (unsigned i = 0; i < t.size(); ++i) {
//		unsigned j = randomUnsigned(t.size());
//		swap(t[i], t[j]);
//	}
//}

///Implements a type converter via streams
template<class OutType, class InType>
OutType stream_cast(const InType & t) {
	stringstream ss;
	ss << t; // first insert value to stream
	OutType result; // value will be converted to OutType
	ss >> result; // write value to result
	return result;
}

//void MakeShuffledDataIndicesList(vector<unsigned>& oDataIdList, unsigned aSize);

//Return an integer hash value for a given input integer in a given domain range
inline uint IntHashSimple(uint key, uint aModulo) {
	key = ~key + (key << 15); // key = (key << 15) - key - 1;
	key = key ^ (key >> 12);
	key = key + (key << 2);
	key = key ^ (key >> 4);
	key = key * 2057; // key = (key + (key << 3)) + (key << 11);
	key = key ^ (key >> 16);
	return key % aModulo;
}

//inline int IntHashPair(int key1, int key2, int aModulo = 2147483647) {
//	const double A = sqrt(2) - 1;
//	int key = key1 * (key2 + 1) * A;
//	return IntHashSimple(key, aModulo);
//}

//Return an integer hash value for a given input integer in a given domain range given an additional seed to select the random hash function
inline uint IntHash(const uint& key, uint& aModulo, unsigned aSeed) {
	const double A = sqrt(2) - 1;
	return IntHashSimple(key * (aSeed + 1) * A, aModulo);
}

//unsigned RSHash(const string& str);
//unsigned RSHash(const vector<unsigned>& aV);
//unsigned APHash(const string& str);
//unsigned APHash(const vector<unsigned>& aV);
//unsigned APHash(const vector<unsigned>::const_iterator& aV_begin, vector<unsigned>::const_iterator& aV_end);
//unsigned HashFunc(const string& str, unsigned aBitMask = 2147483647);
//unsigned HashFunc(const vector<unsigned>& aList, unsigned aBitMask = 2147483647)

inline unsigned APHash(const vector<unsigned>::const_iterator& aV_begin,const vector<unsigned>::const_iterator& aV_end){
	unsigned int hash = 0xAAAAAAAA;
	size_t i=0;
	//for (;aV_begin!=aV_end;aV_begin++) {
	for (auto it = aV_begin; it !=aV_end;it++) {
		hash ^= ((i & 1) == 0) ? ((hash << 7) ^ *it * (hash >> 3)) : (~(((hash << 11) + *it) ^ (hash >> 5)));
		i++;
	}
	return hash;

}

inline unsigned HashFunc(const vector<vector<unsigned>>& array, unsigned x1_min, unsigned x1_max, unsigned& x2, unsigned& bitmask){
	unsigned int hash = 0xAAAAAAAA;

	for (unsigned i = x1_min; i<= x1_max; i++) {
		hash ^= (( ((i-x1_min)) & 1) == 0) ? ((hash << 7) ^ array[i][x2] * (hash >> 3)) : (~(((hash << 11) + array[i][x2]) ^ (hash >> 5)));
	}
	return hash & bitmask;
}


inline unsigned APHashSpec(unsigned aV, const unsigned& bitmask, const unsigned& seed) {
	unsigned int hash = aV;
	hash ^=  ((hash << 7) ^ seed * (hash >> 3));
	//hash ^=   ~((((hash << 11) + seed) ^ (hash >> 5)));
	return hash & bitmask;
}


inline unsigned APHash(const vector<unsigned>& aV) {
	unsigned int hash = 0xAAAAAAAA;
	for (std::size_t i = 0; i < aV.size(); i++) {
		hash ^= ((i & 1) == 0) ? ((hash << 7) ^ aV[i] * (hash >> 3)) : (~(((hash << 11) + aV[i]) ^ (hash >> 5)));
	}
	return hash;
}

inline unsigned HashFunc3(const unsigned& v1, const unsigned& v2, const unsigned& v3, unsigned aBitMask) {
	unsigned int hash = 0xAAAAAAAA;
	hash ^=  ((hash << 7) ^ v1 * (hash >> 3));
	hash ^=  (~(((hash << 11) + v2) ^ (hash >> 5)));
	hash ^= ((hash << 7) ^ v3 * (hash >> 3));
	return hash & aBitMask;
}

inline unsigned HashFunc4(const unsigned& v1, const unsigned& v2, const unsigned& v3, const unsigned& v4, unsigned aBitMask) {
	unsigned int hash = 0xAAAAAAAA;
	hash ^=  ((hash << 7) ^ v1 * (hash >> 3));
	hash ^=  (~(((hash << 11) + v2) ^ (hash >> 5)));
	hash ^= ((hash << 7) ^ v3 * (hash >> 3));
	hash ^=  (~(((hash << 11) + v4) ^ (hash >> 5)));
	return hash & aBitMask;
}


inline unsigned HashFunc(const vector<unsigned>& aList, unsigned aBitMask) {
	return APHash(aList) & aBitMask;
}

//unsigned HashFunc(vector<unsigned>::iterator aList_begin,vector<unsigned>::iterator aList_end,unsigned aBitMask = 2147483647);
inline unsigned HashFunc(vector<unsigned>::iterator aList_begin,vector<unsigned>::iterator aList_end,unsigned aBitMask) {
	return APHash(aList_begin,aList_end) & aBitMask;
}

//------------------------------------------------------------------------------------------------------------------------
///Returns the time between the creation and destruction of an object of TimerClass
class TimerClass {
public:
	TimerClass();
	~TimerClass() {
	}
	void Output();
	double getElapsed();
private:
	std::time_t mStartSec;
	std::clock_t mStart;
	std::chrono::time_point<std::chrono::steady_clock> mStart2;
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

	void 		Begin();
	void 		Count();
	void 		Count(int co);
	double 	getElapsed();
	void 		PrintElapsed();
	unsigned End();
private:
	mutable std::mutex mut;
	unsigned mStep;
	unsigned mCounter;
	TimerClass mTimer;
};

//------------------------------------------------------------------------------------------------------------------------
///Implements safe access policies to a vector container and offers
///members to compute various statistical estimators.
//class VectorClass {
//	friend ostream& operator<<(ostream& out, const VectorClass& aV);
//public:
//	VectorClass();
//	VectorClass(unsigned aSize);
//	void operator=(const VectorClass& aVector);
//	VectorClass(const VectorClass& aVector);
//	VectorClass(const vector<double>& aVector);
//	VectorClass(const vector<int>& aVector);
//	void Init(unsigned aSize);
//	void Import(const string& aFileName);
//	void Clear();
//	unsigned Size() const;
//	ostream& Output(ostream& out) const;
//	void PushBack(double aValue);
//	double& operator[](unsigned i);
//	double operator[](unsigned i) const;
//	double Prod() const;
//	double Sum() const;
//	double Mean() const;
//	double StandardDeviation() const;
//	double Order(double aOrder) const;
//	double Median() const;
//	double MedianAbsoluteDifference() const;
//	double Min() const;
//	double Max() const;
//	VectorClass RemoveNulls();
//	ostream& OutputStatistics(ostream& out);
//protected:
//	vector<double> mV;
//};
//
////------------------------------------------------------------------------------------------------------------------------
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
	std::atomic_uint sizeA;
public:
	threadsafe_queue()
	{ sizeA=0; }
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
		sizeA++;
	}

	void push_unsafe(T new_value)
	{
		//	std::lock_guard<std::mutex> lk(mut);
		data_queue.push(new_value);
		data_cond.notify_one();
		sizeA++;
	}


	void wait_and_pop(T& value)
	{
		std::unique_lock<std::mutex> lk(mut);
		data_cond.wait(lk,[this]{return !data_queue.empty();});
		value=data_queue.front();
		data_queue.pop();
		sizeA--;
	}
	std::shared_ptr<T> wait_and_pop()
																																												{
		std::unique_lock<std::mutex> lk(mut);
		data_cond.wait(lk,[this]{return !data_queue.empty();});
		std::shared_ptr<T> res(std::make_shared<T>(data_queue.front()));
		data_queue.pop();
		sizeA--;
		return res;
																																												}

	bool try_pop(T& value)
	{
		std::lock_guard<std::mutex> lk(mut);
		if(data_queue.empty())
			return false;
		value=data_queue.front();
		data_queue.pop();
		sizeA--;
		return true;
	}
	bool try_pop_unsafe(T& value)
	{
		//	std::lock_guard<std::mutex> lk(mut);
		if(data_queue.empty())
			return false;
		value=data_queue.front();
		data_queue.pop();
		sizeA--;
		return true;
	}

	std::shared_ptr<T> try_pop(){
		std::lock_guard<std::mutex> lk(mut);
		if(data_queue.empty())
			return std::shared_ptr<T>();
		std::shared_ptr<T> res(std::make_shared<T>(data_queue.front()));
		data_queue.pop();
		sizeA--;
		return res;
	}

	bool empty() const
	{
		std::lock_guard<std::mutex> lk(mut);
		return data_queue.empty();
	}

	int size() const
	{
		//		std::lock_guard<std::mutex> lk(mut);
		return sizeA;
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
			if (i==0) cout << endl<<flush;
			cout << i+1 << " ";
		}
		cout << " threads joined" << endl;
	}
};



#endif

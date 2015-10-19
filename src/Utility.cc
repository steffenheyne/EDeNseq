#include "Utility.h"


//void MakeShuffledDataIndicesList(vector<unsigned>& oDataIdList, unsigned aSize) {
//	for (unsigned i = 0; i < aSize; ++i)
//		oDataIdList.push_back(i);
//	std::random_shuffle ( oDataIdList.begin(), oDataIdList.end() );
//	//old manual version
////	for (unsigned i = 0; i < aSize; ++i) {
////		unsigned j = rand() * aSize / RAND_MAX;
////		swap(oDataIdList[i], oDataIdList[j]);
////	}
//}


//----------------------------------------------------------------------------------------------------------------------------
//hash functions

//unsigned RSHash(const string& aString) {
//	unsigned int b = 378551;
//	unsigned int a = 63689;
//	unsigned int hash = 0;
//	for (std::size_t i = 0; i < aString.length(); i++) {
//		hash = hash * a + aString[i];
//		a = a * b;
//	}
//	return hash;
//}
//
//unsigned RSHash(const vector<unsigned>& aV) {
//	unsigned int b = 378551;
//	unsigned int a = 63689;
//	unsigned int hash = 0;
//	for (std::size_t i = 0; i < aV.size(); i++) {
//		hash = hash * a + aV[i];
//		a = a * b;
//	}
//	return hash;
//}

//unsigned APHash(const string& aString) {
//	unsigned int hash = 0xAAAAAAAA;
//	for (std::size_t i = 0; i < aString.length(); i++) {
//		hash ^= ((i & 1) == 0) ? ((hash << 7) ^ aString[i] * (hash >> 3)) : (~(((hash << 11) + aString[i]) ^ (hash >> 5)));
//	}
//	return hash;
//}

//unsigned APHash(const vector<unsigned>& aV) {
//	unsigned int hash = 0xAAAAAAAA;
//	for (std::size_t i = 0; i < aV.size(); i++) {
//		hash ^= ((i & 1) == 0) ? ((hash << 7) ^ aV[i] * (hash >> 3)) : (~(((hash << 11) + aV[i]) ^ (hash >> 5)));
//	}
//	return hash;
//}

//inline unsigned APHash(const vector<unsigned>::const_iterator& aV_begin,const vector<unsigned>::const_iterator& aV_end){
//	unsigned int hash = 0xAAAAAAAA;
//	size_t i=0;
//	//for (;aV_begin!=aV_end;aV_begin++) {
//	for (auto it = aV_begin; it !=aV_end;it++) {
//		hash ^= ((i & 1) == 0) ? ((hash << 7) ^ *it * (hash >> 3)) : (~(((hash << 11) + *it) ^ (hash >> 5)));
//		i++;
//	}
//	return hash;
//
//}

//unsigned HashFunc(const string& aString, unsigned aBitMask) { //NOTE: extract the least significant bits from the hash
//	return APHash(aString) & aBitMask;
//}

//inline unsigned HashFunc(const vector<unsigned>& aList, unsigned aBitMask) {
//	return APHash(aList) & aBitMask;
//}


//unsigned HashFunc(vector<unsigned>::iterator aList_begin,vector<unsigned>::iterator aList_end,unsigned aBitMask) {
//	return APHash(aList_begin,aList_end) & aBitMask;
//}



//------------------------------------------------------------------------------------------------------------------------
TimerClass::TimerClass() :
		mStartSec(time(NULL)), mStart(std::clock()), mStart2(std::chrono::steady_clock::now()) {
}
void TimerClass::Output() {
	//std::clock_t end=std::clock();
	//std::clock_t elapsed=end-start;
	std::time_t end_sec = time(NULL);
	double elapsed_sec = end_sec - mStartSec;
	double elapsed_min = elapsed_sec / 60;
	double elapsed_hour = elapsed_sec / 3600;
	//CLOCKS_PER_SEC
	cout << endl << setprecision(2) << "Elapsed time: h: " << floor(elapsed_hour) << " m: " << floor(elapsed_min) << " s: " << floor(elapsed_sec) << endl;
}

double TimerClass::getElapsed(){
	auto end = std::chrono::steady_clock::now();
	//chrono::duration <double, milli> (diff).count()
   std::chrono::milliseconds elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - mStart2);
	return elapsed.count();
}

//------------------------------------------------------------------------------------------------------------------------
ProgressBar::ProgressBar(unsigned aStep) :
		mStep(aStep), mCounter(0) {
}
ProgressBar::~ProgressBar() {
	if (mCounter > 1) {
		cout << endl << "Counted " << mCounter << " times." << endl;
		mTimer.Output();
	}
}
void ProgressBar::Begin() {
	mCounter = 0;
}

void ProgressBar::Count() {
	std::lock_guard<std::mutex> lk(mut);
	mCounter++;
	double t = 1000/mStep;
	if (mCounter % mStep == 0)
		cout << "." << flush;
	if (mCounter % (mStep*10) == 0)
		cout << mCounter / (mStep*t) << "K" << flush;
}

void ProgressBar::Count(int max) {
	std::lock_guard<std::mutex> lk(mut);
	double t = 1000/mStep;
	for (int i=mCounter; i<max; i++){
		mCounter++;
		if (mCounter % mStep == 0)
			cout << "." << flush;
		if (mCounter % (mStep*10) == 0)
			cout << mCounter / (mStep*t) << "K" << flush;
	}
}


double ProgressBar::getElapsed() {
	return mTimer.getElapsed();
}

unsigned ProgressBar::End() {
	return mCounter;
}

void ProgressBar::PrintElapsed(){
	mTimer.Output();
}

//------------------------------------------------------------------------------------------------------------------------
//ostream& operator<<(ostream& out, const VectorClass& aV) {
//	aV.Output(out);
//	return out;
//}
//VectorClass::VectorClass() {
//}
//VectorClass::VectorClass(unsigned aSize) {
//	Init(aSize);
//}
//void VectorClass::operator=(const VectorClass& aVector) {
//	mV = aVector.mV;
//}
//VectorClass::VectorClass(const VectorClass& aVector) {
//	(*this) = aVector;
//}
//VectorClass::VectorClass(const vector<double>& aVector) {
//	mV = aVector;
//}
//VectorClass::VectorClass(const vector<int>& aVector) {
//	mV.clear();
//	for (unsigned i=0;i<aVector.size();i++)
//	mV.push_back((double)aVector[i]);
//}
//void VectorClass::Init(unsigned aSize) {
//	mV.clear();
//	for (unsigned i = 0; i < aSize; ++i)
//		mV.push_back(-1);
//}
//void VectorClass::Import(const string& aFileName) {
//	mV.clear();
//	ifstream fin;
//	fin.open(aFileName.c_str());
//	if (!fin) {
//		cerr << "Cannot open file: " << aFileName << endl;
//		throw exception();
//	}
//	string line;
//
//	//read size
//	while (getline(fin, line)) {
//		if (line != "") {
//			stringstream ss;
//			ss << line;
//			while (ss.good()) {
//				string value_str;
//				ss >> value_str;
//				if (value_str != "")
//					mV.push_back(stream_cast<double>(value_str));
//			}
//		}
//	}
//	fin.close();
//}
//void VectorClass::Clear() {
//	mV.clear();
//}
//unsigned VectorClass::Size() const {
//	return mV.size();
//}
//ostream& VectorClass::Output(ostream& out) const {
//	for (unsigned i = 0; i < Size(); i++)
//		out << mV[i] << " ";
//	return out;
//}
//void VectorClass::PushBack(double aValue) {
//	mV.push_back(aValue);
//}
//double& VectorClass::operator[](unsigned i) {
//	assert(i<Size());
//	return mV[i];
//}
//double VectorClass::operator[](unsigned i) const {
//	assert(i<Size());
//	return mV[i];
//}
//double VectorClass::Prod() const {
//	double prod = 1;
//	for (unsigned i = 0; i < Size(); i++)
//		prod *= mV[i];
//	return prod;
//}
//double VectorClass::Sum() const {
//	double sum = 0;
//	for (unsigned i = 0; i < Size(); i++)
//		sum += mV[i];
//	return sum;
//}
//double VectorClass::Mean() const {
//	return Sum() / Size();
//}
//double VectorClass::StandardDeviation() const {
//	double avg = Mean();
//	double sd = 0;
//	for (unsigned i = 0; i < Size(); i++)
//		sd += (mV[i] - avg) * (mV[i] - avg);
//	sd = sd / (Size() - 1);
//	sd = sqrt(sd);
//	return sd;
//}
//double VectorClass::Order(double aOrder) const {
//	vector<double> v(mV);
//	sort(v.begin(), v.end());
//	if (aOrder == 0)
//		return v[0];
//	else if (aOrder == 1)
//		return v[v.size() - 1];
//	else {
//		unsigned id=(unsigned) ceil((double) Size() * aOrder);
//		return v[id];
//	}
//}
//double VectorClass::Median() const {
//	return Order(.5);
//}
//double VectorClass::MedianAbsoluteDifference() const {
//	double median = Median();
//	VectorClass v;
//	for (unsigned i = 0; i < Size(); i++)
//		v.PushBack(fabs(mV[i] - median));
//	return v.Median();
//}
//double VectorClass::Min() const {
//	return Order(0);
//}
//double VectorClass::Max() const {
//	return Order(1);
//}
//VectorClass VectorClass::RemoveNulls() {
//	VectorClass v;
//	for (unsigned i = 0; i < Size(); i++)
//		if (mV[i] != -1)
//			v.PushBack(mV[i]);
//	return v;
//}
//ostream& VectorClass::OutputStatistics(ostream& out) {
//	VectorClass v = RemoveNulls();
//	if (v.Size() > 0) {
//		out << "num: " << v.Size() << " sum: " << v.Sum() << " avg: " << v.Mean() << " sd: " << v.StandardDeviation();
//		out << " min: " << v.Min() << " Q.05: " << v.Order(.05) << " Q.25: " << v.Order(.25) << " Q.5: " << v.Median() << " Q.75: " << v.Order(.75) << " Q.95: " << v.Order(.95) << " max: " << v.Max();
//	} else {
//	}
//	return out;
//}

//--------------------------------------------------------------------------------------
OutputManager::OutputManager(string aFileName, string aDirectoryPath) :
		mFileName(aFileName), mDirectoryPath(aDirectoryPath) {
	if (mDirectoryPath != "") {
		//check that directory exists otherwise make one
		int status;
		status = mkdir(mDirectoryPath.c_str(), 0777);
		if (status == -1) {
			if (errno == EEXIST) {
			} //ignore the issue if the directory exists already
			else {
				//otherwise bail out
				string error_code=stream_cast<string>(strerror(errno));
				string msg="ERROR " + error_code + ": Cannot access directory:" + mDirectoryPath;
				throw range_error(msg);
			}
		}
	}

	//check that file can be opened in w mode
	string output_filename = GetFullPathFileName();
	mOut.open(output_filename.c_str(),std::ofstream::trunc);
	if (!mOut)
		throw range_error("ERROR: Cannot open file:" + output_filename);
}

string OutputManager::GetFullPathFileName() {
	string output_filename;
	if (mDirectoryPath != "")
		output_filename = mDirectoryPath + "/" + mFileName;
	else output_filename = mFileName;
	return output_filename;
}



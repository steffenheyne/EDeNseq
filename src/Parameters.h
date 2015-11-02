/* -*- mode:c++ -*- */
#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "Utility.h"

using namespace std;

const string PROG_VERSION="1.0";
const string PROG_DATE="Autumn 2015";
const string PROG_NAME = "EDeNseq (Explicit Decomposition with Neighborhoods for Biological Sequences)";
const string PROG_CREDIT = "Author: Steffen Heyne & Fabrizio Costa     Email:heyne@ie-freiburg.mpg.de";
const string CITATIONS = "For the NSPDK graph kernels see Fabrizio Costa, Kurt De Grave, ''Fast Neighborhood Subgraph Pairwise Distance Kernel'', Proceedings of the 27th International Conference on Machine Learning (ICML-2010), Haifa, Israel, 2010..\n";


enum ActionType {
	NULL_ACTION, CLUSTER, CLASSIFY, TEST
};

enum InputFileType {
	STRINGSEQ, FASTA
};

//------------------------------------------------------------------------------------------------------------------------
enum OptionsType {
	FLAG, LIST, REAL, INTEGER, POSITIVE_INTEGER, STRING
};

//------------------------------------------------------------------------------------------------------------------------
class ParameterType {
public:
	string mShortSwitch;
	string mLongSwitch;
	string mShortDescription;
	string mLongDescription;
	OptionsType mTypeCode;
	string mValue;
	string mMinValue;
	string mMaxValue;
	string mNumStepsValue;
	bool mExponentialStepIncrease;
	bool mIsSet;
	vector<string> mCloseValuesList;

public:
	ParameterType();
public:
	void Parse(vector<string>& aParameterList);
	void OutputCompact(ostream& out) const;
	void OutputExtended(ostream& out) const;
};

//------------------------------------------------------------------------------------------------------------------------
class Parameters {
public:
	map<string, ParameterType> mOptionList;
	map<ActionType, vector<ParameterType*> > mActionOptionList;
	map<ActionType, string> mActionSummary;
	map<ActionType, string> mActionReferences;

	string mProgramVersion;
	string mAction;
	ActionType mActionCode;
	string mInputDataFileName;
	string mFileType;
	InputFileType mFileTypeCode;
	unsigned mRadius;
	unsigned mDistance;
	unsigned mHashBitSize;
	unsigned mRandomSeed;
	string mKernelType;
	string mSuffix;

	bool mVerbose;
	unsigned mNumThreads;
	unsigned mNumIndexThreads;

	unsigned mNumHashFunctions;
	unsigned mNumRepeatsHashFunction;
	double mMaxSizeBin;
	string mClusterType;
	unsigned mNumHashShingles;
	double mPureApproximateSim;

	string mDirectoryPath;

	// MetaGenome
	string mIndexBedFile;
	string mIndexSeqFile;
	bool mNoIndexCacheFile;
	unsigned mSeqWindow;
	unsigned mIndexSeqShift;
	unsigned mSeqShift;

	unsigned mSeqClip;
	unsigned mMinRadius;
	unsigned mMinDistance;
	string mDenseCenterNamesFile;
	bool mWriteApproxNeighbors;

public:
	Parameters();
	void SetupOptions();
	void Usage(string aCommandName, string aCompactOrExtended);
	void Init(int argc, const char** argv);
};

#endif /* PARAMETERS_H */

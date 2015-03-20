/* -*- mode:c++ -*- */
#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "Utility.h"

using namespace std;

const string PROG_VERSION="1.3.7";
const string PROG_DATE="30 May 2014";
const string PROG_NAME = "EDeN (Explicit Decomposition with Neighborhoods)";
const string PROG_CREDIT = "Author: Fabrizio Costa     Email:costa@informatik.uni-freiburg.de";
const string CITATIONS = "For the Decompositional DAG graph kernel see G. Da San Martino, N. Navarin and A. Sperduti , ''A tree-based kernel for graphs'', Proceedings of the Twelfth SIAM International Conference on Data Mining, Anaheim, California, April 26 - 28, 2012, p. 975-986.\n"
		"For the NSPDK graph kernels see Fabrizio Costa, Kurt De Grave, ''Fast Neighborhood Subgraph Pairwise Distance Kernel'', Proceedings of the 27th International Conference on Machine Learning (ICML-2010), Haifa, Israel, 2010..\n"
		"The code for Stochastic Gradient Descent SVM is adapted from http://leon.bottou.org/projects/sgd. Léon Bottou and Yann LeCun, ''Large Scale Online Learning'', Advances in Neural Information Processing Systems 16, Edited by Sebastian Thrun, Lawrence Saul and Bernhard Schölkopf, MIT Press, Cambridge, MA, 2004.\n"
		"The embedding method is adapted from L.Chen, A.Buja ''Local Multidimensional Scaling for Nonlinear Dimension Reduction, Graph Drawing, and Proximity Analysis'', Journal of the American Statistical Association, 2009.";

enum ActionType {
	NULL_ACTION, TRAIN, TEST, PARAMETERS_OPTIMIZATION, CROSS_VALIDATION, BIAS_VARIANCE_DECOMPOSITION,
	LEARNING_CURVE, TEST_PART, FEATURE, FEATURE_PART, FEATURE_SCALED, MATRIX, EMBED, TARGET_ALIGNMENT,
	CLUSTER, NEAREST_NEIGHBOR, MIN_HASH, SEMI_SUPERVISED, PRE_CLUSTER_FILTER, DENSITY, TRAIN_MULTILABEL, TEST_MULTILABEL,
	GENOME
};

enum InputFileType {
	SPARSE_VECTOR, GRAPH, MOLECULAR_GRAPH, SEQUENCE, STRINGSEQ, FASTA
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
	string mTargetFileName;
	string mRowIndexFileName;
	string mColIndexFileName;
	string mModelFileName;
	string mFileType;
	string mEmbedFileName;
	bool mBinaryFormat;
	string mOpenBabelFormat;
	InputFileType mFileTypeCode;
	bool mKernelNoNormalization;
	bool mMinKernel;
	bool mUseRealVectorInformation;
	unsigned mRadius;
	unsigned mDistance;
	unsigned mVertexDegreeThreshold;
	double mLambda;
	unsigned mEpochs;
	unsigned mHashBitSize;
	unsigned mCrossValidationNumFolds;
	unsigned mNumPoints;
	unsigned mRandomSeed;
	string mKernelType;
	string mGraphType;
	unsigned mSemiSupervisedNumIterations;
	double mSemiSupervisedThreshold;
	bool mSemiSupervisedInduceOnlyPositive;
	bool mSemiSupervisedInduceOnlyNegative;
	string mSuffix;
	unsigned mSequenceDegree;
	unsigned mLMDSNumRandomRestarts;
	unsigned mLMDSNumIterations;
	unsigned mLMDSDimensionality;
	double mLMDSIterationEpsilon;
	unsigned mLMDSNeighborhoodSize;
	unsigned mLMDSNonNeighborhoodSize;
	unsigned mLMDSNeighborhoodSizeRange;
	double mLMDSTau;
	unsigned mLMDSTauExponentRange;
	unsigned mLMDSRefineNumNearestNeighbor;

	bool mVerbose;
	bool mMinimalOutput;
	bool mSequenceToken;
	bool mSequenceMultiLine;
	bool mSequencePairwiseInteraction;
	unsigned mSparsificationNumIterations;
	unsigned mTopologicalRegularizationNumNeighbors;
	double mTopologicalRegularizationRate;
	unsigned mNumLineSearchIterations;
	unsigned mNumThreads;

	unsigned mNumHashFunctions;
	unsigned mNumRepeatsHashFunction;
	double mMaxSizeBin;
	double mEccessNeighborSizeFactor;
	unsigned mSampleSize;
	unsigned mNumNearestNeighbors;
	bool mSharedNeighborhood;
	double mFractionCenterScan;
	unsigned mNeighborhoodIntersectionSize;
	double mClusterThreshold;
	string mClusterType;
	bool mNoNeighborhoodCache;
	bool mNoMinHashCache;
	bool mUseApproximate;
	double mPureApproximateSim;
	unsigned mNumHashShinglets;

	double mSemiSupervisedAlpha;

	//DD kernel family additional parameters
	double mTreeLambda;
	double mRadiusTwo;

	bool mSmooth;
	double mSmootherParam;

	string mDirectoryPath;
	bool mExtendedMatrixInformation;

	unsigned mNumOutputInstances;

	unsigned mSizeRandomDataset;
	unsigned mNumReshuffles;
	double mMaxFractionOfDataset;
	double mStdDevRandomSampling;

	unsigned mNumRandomProjections;

	// MetaGenome
	string mIndexDataList;
	unsigned mSeqWindow;
	double mSeqShift;
	unsigned mMinRadius;
	unsigned mMinDistance;

public:
	Parameters();
	void SetupOptions();
	void Usage(string aCommandName, string aCompactOrExtended);
	void Init(int argc, const char** argv);
};

#endif /* PARAMETERS_H */

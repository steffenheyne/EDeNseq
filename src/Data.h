/* -*- mode:c++ -*- */
#ifndef DATA_H
#define DATA_H

#include "Utility.h"
#include "Parameters.h"
//#include "Kernel.h"
#include "gzstream.h"
#include "GraphClass.h"

using namespace std;

const unsigned BUFFER_SIZE = 100;

struct SeqDataSet {
	string filename;
	unsigned idx;
	string desc;
	InputFileType filetype;
};

class Data {
	public:
		Parameters* mpParameters;
		//Kernel mKernel;
		vector<double> mTargetList;
		vector<vector<double> > mMultilabelTargetList;
		vector<SVector> mVectorList;
		map<unsigned, SVector> mFeatureCorrelationMatrix;
		vector<unsigned> mRowIndexList;
		vector<unsigned> mColIndexList;
		bool mDataIsLoaded;
	protected:
		bool mMultilabelTargetIsLoaded;
		bool mTargetIsLoaded;
		bool mIndexIsLoaded;
		unsigned mDataSize;
	public:
		Data();
		void Init(Parameters* apParameters);
		double ComputeKernel(unsigned i, unsigned j);
		double ComputeKernel(SVector& x, SVector& z);
		void ComputeFeatureCorrelationMatrix();
		void LoadIndex();
		bool IsIndexLoaded();
		void LoadMultilabelTarget();
		bool IsMultilabelTargetLoaded();
		void LoadMultilabelRealList(string aFileName, vector<vector<double> >& oList);
		unsigned MultilabelTargetDimension();
		void LoadTarget();
		bool IsTargetLoaded();
		void LoadRealList(string aFileName, vector<double>& oList);
		void LoadUnsignedList(string aFileName, vector<unsigned>& oList);
		void LoadGspanList(vector<GraphClass>& aGraphList);
		void LoadData(bool aLoadIndex, bool aLoadTarget, bool aLoadMultilabelTarget);
		void LoadData(istream& fin);
		bool IsDataLoaded();
		void SetGraphFromFile(istream& in, GraphClass& oG);
#ifdef USEOBABEL
		void SetGraphFromGraphOpenBabelFile(istream& in, GraphClass& oG, string aFormat);
#endif
		bool SetVectorFromSparseVectorAsciiFile(istream& in, SVector& aX);
		bool SetVectorFromSparseVectorBinaryFile(istream& in, SVector& aX);
		bool SetGraphFromStringFile(istream& in, GraphClass& oG);
		bool SetGraphFromFASTAFile(istream& in, GraphClass& oG);
		void SetGraphFromSequenceFile(istream& in, GraphClass& oG);
		void SetGraphFromSequenceTokenFile(istream& in, GraphClass& oG);
		void SetGraphFromSequenceMultiLineFile(istream& in, GraphClass& oG);
		void SetGraphFromSequenceMultiLineTokenFile(istream& in, GraphClass& oG);
		void ManageDirectedAndMultiComponent(GraphClass& oG, vector<vector<unsigned> >& aVertexComponentList);
		void AddVertexAndEdgesForSequence(GraphClass& oG, string aLabel, unsigned aVertexCounter, bool aGraphDisconnect);
		void AddAbstractConnections(GraphClass& oG, vector<vector<unsigned> >& aVertexComponentList);
		void SetGraphFromGraphGspanFile(istream& in, GraphClass& oG);
		void AddReverseGraph(GraphClass& oG);
		unsigned Size();
		void SetDataSize(unsigned aSize);
		void Bootstrap(vector<unsigned>& oBootstrappedSample, vector<unsigned>& oOutOfBagSample);
		void RandomPartition(double aPartitionRatio, vector<unsigned>& oFirstSample, vector<unsigned>& oSecondSample);
		void RandomPartition(double aPartitionRatio, vector<bool>& oBitVector);
		void GenerateRandomPermutationInstance(SVector& oX);
		void GenerateRandomPermutationDataset(unsigned aSize, vector<SVector>& oDataset);
		void GenerateRandomSparseInstance(SVector& oX, const SVector& aRef, double aStdDev);
		void GenerateRandomInstance(SVector& oX, const SVector& aRef, double aStdDev);
		void GenerateRandomInstance(SVector& oX);
		void GenerateRandomDataset(unsigned aSize, vector<SVector>& oDataset);
};

#endif /* DATA_H */

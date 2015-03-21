/* -*- mode:c++ -*- */
#ifndef DATA_H
#define DATA_H

#include "Utility.h"
#include "Parameters.h"
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
		vector<double> mTargetList;
		vector<vector<double> > mMultilabelTargetList;
		vector<SVector> mVectorList;
		map<unsigned, SVector> mFeatureCorrelationMatrix;
		vector<unsigned> mRowIndexList;
		vector<unsigned> mColIndexList;
		bool mDataIsLoaded;
	protected:
		unsigned mDataSize;
	public:
		Data();
		void Init(Parameters* apParameters);
		void LoadMultilabelTarget();
		void LoadMultilabelRealList(string aFileName, vector<vector<double> >& oList);
		void LoadTarget();
		void LoadRealList(string aFileName, vector<double>& oList);
		void LoadUnsignedList(string aFileName, vector<unsigned>& oList);
		void LoadData(bool aLoadIndex, bool aLoadTarget, bool aLoadMultilabelTarget);
		void LoadData(istream& fin);
		bool IsDataLoaded();
		void SetGraphFromFile(istream& in, GraphClass& oG);
		bool SetGraphFromStringFile(istream& in, GraphClass& oG);

		unsigned Size();
		void SetDataSize(unsigned aSize);
};

#endif /* DATA_H */

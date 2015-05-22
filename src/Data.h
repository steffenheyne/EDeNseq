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
	unsigned idx;
	string filename;
	string desc;
	InputFileType filetype;
	unsigned numSequences;
	bool updateIndex;
	bool updateSigCache;
	//vector<vector<unsigned> > sigCache;
};

class Data {
public:
	Parameters* mpParameters;

protected:
	unsigned mDataSize;
public:
	Data();
	void Init(Parameters* apParameters);
	vector<SeqDataSet> LoadIndexDataList(string filename);
	void SetGraphFromFile(istream& in, GraphClass& oG);
	bool SetGraphFromFASTAFile(istream& in, GraphClass& oG, string& currSeq, uint& pos, string& name);
	bool SetGraphFromSeq(string& seq, GraphClass& oG);
	bool SetGraphFromStringFile(istream& in, GraphClass& oG);
	string GetNextFastaSeq(istream& in,string& header);
	bool SetGraphFromSeq(GraphClass& oG, string& currSeq);
	unsigned Size();
	void SetDataSize(unsigned aSize);
};

#endif /* DATA_H */

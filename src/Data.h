/* -*- mode:c++ -*- */
#ifndef DATA_H
#define DATA_H

#include "Utility.h"
#include "Parameters.h"
#include "gzstream.h"
#include "GraphClass.h"

using namespace std;

class Data {
	
public:
	typedef vector<vector<unsigned> > SigCacheT;
	typedef std::shared_ptr<SigCacheT> SigCacheP;

	struct BEDentryS {
		// cols 1-6 from file
		string SEQ;
		uint START;
		uint END;
		string NAME;
		double SCORE;
		char	STRAND;
		// col7 and beyond
		string INFO;
	};

	typedef BEDentryS BEDentryT;
	typedef multimap<string, BEDentryT> BEDdataT;
	typedef BEDdataT::iterator BEDdataIt;
	typedef std::shared_ptr<BEDdataT> BEDdataP;


	Parameters* mpParameters;

public:
	Data() {};
	Data(Parameters* apParamters);
	void Init(Parameters* apParameters);
	//	vector<SeqDataSet> LoadIndexDataList(string filename);
	BEDdataP 			LoadBEDfile(string filename);

	//	void SetGraphFromFile(istream& in, GraphClass& oG);
	//	bool SetGraphFromFASTAFile(istream& in, GraphClass& oG, string& currSeq, unsigned& pos, string& name);
	bool SetGraphFromSeq(string& seq, GraphClass& oG);
	bool SetGraphFromSeq2(GraphClass& oG, string& currSeq, unsigned& pos);
	//	bool SetGraphFromStringFile(istream& in, GraphClass& oG);
	void GetNextFastaSeq(istream& in,string& currSeq, string& header);
	void GetNextStringSeq(istream& in,string& currSeq, string& header);
	//	bool SetGraphFromSeq(GraphClass& oG, string& currSeq);

	void LoadStringList(string aFileName, vector<string>& oList, uint numTokens);
};

#endif /* DATA_H */

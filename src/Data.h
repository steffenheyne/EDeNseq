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
		vector<string> COLS;
	};
	typedef BEDentryS BEDentryT;
	typedef std::shared_ptr<BEDentryS> BEDentryP;
	typedef multimap<string, BEDentryP> BEDdataT;
	typedef BEDdataT::iterator BEDdataIt;
	typedef std::shared_ptr<BEDdataT> BEDdataP;


	Parameters* mpParameters;

public:
	Data() {};
	Data(Parameters* apParamters);

	void Init(Parameters* apParameters);

	BEDdataP	LoadBEDfile(string filename);
	bool SetGraphFromSeq2(GraphClass& oG, string& currSeq, unsigned& pos, bool& lastGr, string& seq);
	bool SetGraphFromSeq(string& seq, GraphClass& oG);
	void GetNextFastaSeq(istream& in,string& currSeq, string& header);
	void GetNextStringSeq(istream& in,string& currSeq);
	void LoadStringList(string aFileName, vector<string>& oList, uint numTokens);

	//	vector<SeqDataSet> LoadIndexDataList(string filename);
	//	void SetGraphFromFile(istream& in, GraphClass& oG);
	//	bool SetGraphFromFASTAFile(istream& in, GraphClass& oG, string& currSeq, unsigned& pos, string& name);
	//bool SetGraphFromSeq(GraphClass& oG, string& currSeq);
	//	bool SetGraphFromStringFile(istream& in, GraphClass& oG);
};

#endif /* DATA_H */

/* -*- mode:c++ -*- */
#ifndef DATA_H
#define DATA_H

#include "Utility.h"
#include "Parameters.h"
#include "gzstream.h"

using namespace std;

class Data {

public:

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
	bool GetNextWinFromSeq(string& currSeq, unsigned& pos, bool& lastGr, string& seq);
	//bool SetGraphFromSeq2(GraphClass& oG, string& currSeq, unsigned& pos, bool& lastGr, string& seq);
	//bool SetGraphFromSeq(string& seq, GraphClass& oG);
	void GetRevComplSeq(string& in_seq,string& out_seq);
	void GetNextFastaSeq(istream& in,string& currSeq, string& header);
	void GetNextStringSeq(istream& in,string& currSeq);
	void LoadStringList(string aFileName, vector<string>& oList, uint numTokens);

	void GetNextLargeWinFromSeq(string& currSeq, unsigned& pos, bool& lastGr, string& seq, unsigned win_size_large, unsigned win_size_small, unsigned shift);
	//	vector<SeqDataSet> LoadIndexDataList(string filename);
	//	void SetGraphFromFile(istream& in, GraphClass& oG);
	//	bool SetGraphFromFASTAFile(istream& in, GraphClass& oG, string& currSeq, unsigned& pos, string& name);
	//bool SetGraphFromSeq(GraphClass& oG, string& currSeq);
	//	bool SetGraphFromStringFile(istream& in, GraphClass& oG);
};

#endif /* DATA_H */

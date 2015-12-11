#include "Data.h"

Data::Data(Parameters* apParameters) :
mpParameters(apParameters) {
}

void Data::Init(Parameters* apParameters){
	mpParameters = apParameters;
}

Data::BEDdataP Data::LoadBEDfile(string filename){

	Data::BEDdataP myBED = std::make_shared<BEDdataT>();
	igzstream fin;
	string line;

	fin.open(filename.c_str());
	if (!fin)
		throw range_error("ERROR LoadData: Cannot open index data file: " + filename);

	while (!fin.eof()) {

		char line2[2048];
		fin.getline(line2,2048);
		string tmp(line2);
		istringstream iss2(line2,istringstream::in);

		BEDentryP myEnt = std::make_shared<BEDentryT>();

		if (iss2 >> myEnt->SEQ >> myEnt->START >> myEnt->END >> myEnt->NAME){

			string col;
			double sc;
			string str;

			if (iss2 >> sc ){
				myEnt->SCORE = sc;
			} else if (!iss2.eof()) {
				string tmp(line2);
				cout << "BED file line:\n" << tmp << endl;
				throw range_error("ERROR BEDfile format col 5 error. Expected numerical value");
			} else
				myEnt->SCORE = 0.0;

			if (iss2 >> str && (str=="." || str=="-" || str =="+") ){
				myEnt->STRAND = str[0];
			} else if (!iss2.eof()) {
				string tmp(line2);
				cout << "BED file line:\n" << tmp << endl;
				throw range_error("ERROR BEDfile format col 6 error. Expected strand symbol: + - . ");
			} else {
				myEnt->STRAND = '.';
			}

			while (iss2 >> col) {
				myEnt->COLS.push_back(col);
			}

			myBED->insert(make_pair(myEnt->SEQ,myEnt));

		} else if (tmp.size()>0 && tmp[0]!='#') {
			cout << endl << "Wrong BED line found! Expect at least 4 columns: <STRING> <INT> <INT> <STRING>! Comment lines start with '#'! Line found:" << endl;
			cout << tmp << endl;
			throw range_error("ERROR: BEDfile format error.\n");
		}
	}
	fin.close();
	cout << "BED file loaded with " <<  myBED->size() << " entries" << endl;
	if (!myBED->size())
		throw range_error("ERROR LoadIndexData: No data found in " + filename + "!");

	return myBED;
}


void Data::GetNextFastaSeq(istream& in,string& currSeq, string& header) {

	in >> std::ws;

	char c = in.peek();
	currSeq.clear();
	header.clear();

	if (!in.eof() && c != EOF && c=='>' ){
		getline(in, header,'>');
		getline(in, header);
		getline(in, currSeq,'>');
		currSeq.erase(std::remove(currSeq.begin(), currSeq.end(), '\n'),currSeq.end());
		currSeq.erase(std::remove(currSeq.begin(), currSeq.end(), ' '),currSeq.end());
		std::transform(currSeq.begin(), currSeq.end(), currSeq.begin(), ::toupper);

		//string seq = currSeq.substr(mpParameters->mSeqClip,currSeq.size()-(2*mpParameters->mSeqClip));
		//currSeq=seq;

		const unsigned pos = header.find_first_of(" ");
		if (std::string::npos != pos)
			header = header.substr(0,pos);
		in.unget();
		if (currSeq.size()==0 || header.size()==0)
			throw range_error("ERROR FASTA reader - empty Sequence or header found! Header:"+header);
		//cout << " found seq " << header << " " << currSeq.size() << " length" << " EOF? "<< in.eof() << endl;
	} else if (c != '>' && c != EOF && c!= '\n') {
		throw range_error("ERROR FASTA format error  -2-!");
	}
}

void Data::GetNextLargeWinFromSeq(string& currSeq, unsigned& pos, bool& last, string& seq, unsigned win_size_large, unsigned win_size_small, unsigned shift){

	// pos is 0 based

	// estimate number of shifts that fit into win_size_large
	unsigned num_shifts = std::ceil( (double)( win_size_large - win_size_small) / (double)shift);

	// minimal shifts that should fit is the fragment fter the current one
	unsigned min_shifts_left = 50;

	// we either return a large window that hold exactly num_shift small windows
	// or the rest of the sequence
	if (currSeq.size()-pos <= (shift * (num_shifts+min_shifts_left) + win_size_small) ){
		seq = currSeq.substr(pos,currSeq.size()-pos);
		pos = currSeq.size();
		last = true;
	} else {
		seq = currSeq.substr(pos, (shift * num_shifts + win_size_small) );
		pos += shift * (num_shifts+1);
		last = false;
	}
	//cout << "GetNext " << seq.size()<< " pos " << pos << " last " << last << " numshifts " << num_shifts << endl;

}

// 1..100       0*20+1
//  21..120     1*20+1
//   41..140    2*20+1
//    61..160
//            lws =   shift*x + ws

bool Data::GetNextWinFromSeq(string& currSeq, unsigned& pos, bool& lastGr, string& seq){

	bool success_status = false;
	lastGr = false;
	unsigned win=mpParameters->mSeqWindow;
	unsigned shift = mpParameters->mSeqShift;

	if (currSeq.size() > pos ) {

		// default case for window/shift
		unsigned currSize = win;
		// case no window/shift
		if (win==0){
			unsigned clipSize = mpParameters->mSeqClip;
			currSize = currSeq.size()-(2*clipSize);
			pos = clipSize;
		} else if (win>currSeq.size()-pos) {
			// case seq left is smaller than win
			// then we take a full window from the end
			pos = std::max((int)0,((int)currSeq.size()-(int)win));
			currSize=currSeq.size()-pos;
		}

		if (currSize>=1){
			seq = currSeq.substr(pos,currSize);
			//cout << currSeq.size() << " " << currSize  << " pos " << pos << " win " << win << " " << seq << endl;
		} else
			throw range_error("ERROR FASTA reader! Too short sequence found. " + std::to_string(pos));

		if ((win>0) && (currSeq.size()-pos-shift>=win)){
			pos += shift;
		} else if ((currSeq.size()-shift-pos<win ) || (win == 0))	{
			currSeq="";
			pos = 0;
			lastGr=true;
		}
		success_status = true;
	}
	return success_status;
}

//
//bool Data::SetGraphFromSeq(string& seq, GraphClass& oG) {
//	vector<bool> vertex_status(5, false);
//	vertex_status[0] = true; //kernel point
//	vertex_status[1] = true; //kind
//	vertex_status[2] = true; //viewpoint
//	vertex_status[3] = false; //dead
//	vertex_status[4] = false; //abstraction
//
//	bool success_status = false;
//	//cout << graphSeq << " " << currSeq.size() << " " << currSize << " " << in.eof() << endl;
//	unsigned real_vertex_index = oG.InsertVertex();
//	vector<string> vertex_symbolic_attribute_list(1);
//	vertex_symbolic_attribute_list[0] = seq;
//	oG.SetVertexSymbolicAttributeList(real_vertex_index, vertex_symbolic_attribute_list);
//	oG.SetVertexStatusAttributeList(real_vertex_index, vertex_status);
//	success_status = true;
//
//	return success_status;
//}


void Data::GetRevComplSeq(string& in_seq,string& out_seq){

	out_seq = "";
	for (std::string::reverse_iterator rit=in_seq.rbegin(); rit!=in_seq.rend(); ++rit) {
		switch (*rit) {
			case 'A': out_seq.push_back('T');
		 		break;
		 	case 'T': out_seq.push_back('A');
		 		break;
		 	case 'G': out_seq.push_back('C');
		 		break;
		 	case 'C': out_seq.push_back('G');
		 		break;
		 	default:
		 		out_seq.push_back('N');
		 		break;
		}
	}
}

void Data::GetNextStringSeq(istream& in,string& currSeq) {

	currSeq.clear();
	getline(in, currSeq);
}


void Data::LoadStringList(string aFileName, vector<string>& oList, uint numTokens) {
	oList.clear();
	cout << endl << "Reading file: " << aFileName << " ..";
	ifstream fin;
	fin.open(aFileName.c_str());
	if (!fin)
		throw range_error("ERROR Data::LoadStringList: Cannot open file:" + aFileName);
	while (!fin.eof()) {
		string line;
		getline(fin, line);
		stringstream ss;
		ss << line << endl;
		uint tok = 0;
		while (!ss.eof() && tok<numTokens) {
			string value;
			ss >> value;
			if (ss.good()) {
				oList.push_back(value);
			}
			tok++;
		}
	}
	fin.close();
	cout << ".. read: " << oList.size() << " values." << endl;
}

//bool Data::SetGraphFromSeq(GraphClass& oG, string& currSeq) {
//	vector<bool> vertex_status(5, false);
//	vertex_status[0] = true; //kernel point
//	vertex_status[1] = true; //kind
//	vertex_status[2] = true; //viewpoint
//	vertex_status[3] = false; //dead
//	vertex_status[4] = false; //abstraction
//
//	bool success_status = false;
//
//	unsigned win=mpParameters->mSeqWindow;
//	unsigned shift = (unsigned)((double)win*mpParameters->mSeqShift);
//
//	if (currSeq.size() > 0 ) {
//
//		// default case for window/shift
//		unsigned currSize = win;
//		// case now window/shift
//		if (win==0){
//			currSize = currSeq.size();
//		} else if (win>currSeq.size()) {
//			// case seq left is smaller than win
//			currSize=currSeq.size();
//		}
//
//		if (currSize>=win){
//
//			string graphSeq = currSeq.substr(0,currSize);
//			//cout << graphSeq << " " << currSeq.size() << " " << currSize << " " << in.eof() << endl;
//			unsigned real_vertex_index = oG.InsertVertex();
//			vector<string> vertex_symbolic_attribute_list(1);
//			vertex_symbolic_attribute_list[0] = graphSeq;
//			oG.SetVertexSymbolicAttributeList(real_vertex_index, vertex_symbolic_attribute_list);
//			oG.SetVertexStatusAttributeList(real_vertex_index, vertex_status);
//		}
//		if (win>0 && currSeq.size()-shift>=win){
//			currSeq.erase(0,shift);
//		} else if (currSeq.size()-shift<win || win == 0)
//			currSeq="";
//		success_status = true;
//	}
//	return success_status;
//}




//bool Data::SetGraphFromFASTAFile(istream& in, GraphClass& oG, string& currSeq, unsigned& pos, string& name) {
//
//	bool success_status = false;
//
//	unsigned win=mpParameters->mSeqWindow;
//	unsigned shift = (unsigned)((double)win*mpParameters->mSeqShift);
//
//	if (currSeq.size() == 0){
//		in >> std::ws;
//	}
//
//	char c = in.peek();
//	string header;
//	bool newSeq=false;
//	//cout << "here "<< " " << currSeq.size() << " " << in.eof() << " c "  << c << endl;
//	if (!in.eof() && c != EOF && c=='>' && currSeq.size() == 0 ){
//		getline(in, header,'>');
//		getline(in, header);
//		getline(in, currSeq,'>');
//
//		name = header;
//		currSeq.erase(std::remove(currSeq.begin(), currSeq.end(), '\n'),currSeq.end());
//		currSeq.erase(std::remove(currSeq.begin(), currSeq.end(), ' '),currSeq.end());
//		std::transform(currSeq.begin(), currSeq.end(), currSeq.begin(), ::toupper);
//		in.unget();
//		if (currSeq.size()==0 || header.size()==0)
//			throw range_error("ERROR FASTA reader - empty Sequence or header found! Header:"+header);
//		//cout << " found seq " << header << " " << currSeq.size() << " length" << " EOF? "<< in.eof() << endl;
//		newSeq=true;
//		pos=0;
//	} else if (c != '>' && c != EOF && c!= '\n') {
//		throw range_error("ERROR FASTA format error!");
//	}
//
//	if (currSeq.size() > pos ) {
//
//		// default case for window/shift
//		unsigned currSize = win;
//		// case no window/shift
//		if (win==0){
//			uint clipSize = 0;
//			currSize = currSeq.size()-(2*clipSize);
//			pos = clipSize;
//		} else if (win>currSeq.size()-pos) {
//			// case seq left is smaller than win
//			currSize=currSeq.size()-pos;
//		}
//
//		if (currSize>=win){
//			string seq = currSeq.substr(pos,currSize);
//			SetGraphFromSeq( seq ,oG);
//			//cout << currSeq.size() << " " << currSize << " eof? " << in.eof() << " pos " << pos << " win " << win << " " << seq << endl;
//		} else if (newSeq){
//			throw range_error("ERROR FASTA reader! Too short sequence found. Either use win=0 or increase window size!");
//		}
//
//
//		if ((win>0) && (currSeq.size()-pos-shift>=win)){
//			pos += shift;
//		} else if ((currSeq.size()-shift-pos<win) || (win == 0))
//			{
//				currSeq="";
//				pos = 0;
//			}
//		success_status = true;
//	}
//	return success_status;
//}

//bool Data::SetGraphFromStringFile(istream& in, GraphClass& oG) {
//
//	bool success_status = false;
//	string line;
//	getline(in, line);
//	if (line == "")
//		return false;
//
//	SetGraphFromSeq(line,oG);
//
//	success_status = true;
//	return success_status;
//}

//void Data::SetGraphFromFile(istream& in, GraphClass& oG) {
//	switch (mpParameters->mFileTypeCode) {
//	case STRINGSEQ: {
//		SetGraphFromStringFile(in, oG);
//		break;
//	}
//	default:
//		throw range_error("ERROR Data::SetGraphFromFile: file type not recognized: " + mpParameters->mFileType);
//	}
//	//******************************************************************************************************
//}

//vector<SeqDataSet> Data::LoadIndexDataList(string filename){
//
//	vector<SeqDataSet> myList;
//	bool valid_input = true;
//	igzstream fin;
//	string line;
//
//	fin.open(filename.c_str());
//	if (!fin)
//		throw range_error("ERROR LoadData: Cannot open index data file: " + filename);
//	while (!fin.eof() && valid_input) {
//		SeqDataSet mySet;
//		mySet.filetype=FASTA;
//		if (fin >> mySet.uIdx >> mySet.filename >> mySet.desc){
//			cout << "found file idx " << mySet.uIdx << "\t" << mySet.filename << "\t" << mySet.desc << endl;
//			mySet.idx=0;
//			mySet.updateIndex=true;
//			mySet.updateSigCache=false;
//			myList.push_back(mySet);
//		}
//		getline(fin, line);
//	}
//	fin.close();
//	if (!myList.size())
//		throw range_error("ERROR LoadIndexData: No data found in " + filename + "!");
//
//	return myList;
//}


//bool Data::SetGraphFromSeq2(GraphClass& oG, string& currSeq, unsigned& pos, bool& lastGr, string& seq) {
//
//	bool success_status = false;
//	lastGr = false;
//	unsigned win=mpParameters->mSeqWindow;
//	unsigned shift = std::max((double)1,(double)win*mpParameters->mSeqShift);
//
//	if (currSeq.size() > pos ) {
//
//		// default case for window/shift
//		unsigned currSize = win;
//		// case no window/shift
//		if (win==0){
//			unsigned clipSize = mpParameters->mSeqClip;
//			currSize = currSeq.size()-(2*clipSize);
//			pos = clipSize;
//		} else if (win>currSeq.size()-pos) {
//			// case seq left is smaller than win
//			// then we take a full window from the end
//			pos = std::max((int)0,((int)currSeq.size()-(int)win));
//			currSize=currSeq.size()-pos;
//		}
//
//		if (currSize>=1){
//			seq = currSeq.substr(pos,currSize);
//			SetGraphFromSeq( seq ,oG);
//			//cout << currSeq.size() << " " << currSize << " eof? " << in.eof() << " pos " << pos << " win " << win << " " << seq << endl;
//		} else
//			throw range_error("ERROR FASTA reader! Too short sequence found. " + std::to_string(pos));
//
//		if ((win>0) && (currSeq.size()-pos-shift>=win)){
//			pos += shift;
//		} else if ((currSeq.size()-shift-pos<win ) || (win == 0))
//		{
//			currSeq="";
//			pos = 0;
//			lastGr=true;
//		}
//		success_status = true;
//	}
//	return success_status;
//}

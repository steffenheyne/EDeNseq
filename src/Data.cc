#include "Data.h"

Data::Data(Parameters* apParameters) :
mpParameters(apParameters) {
}

void Data::Init(Parameters* apParameters){
	mpParameters = apParameters;
}

Data::BEDdataP Data::LoadBEDfile(string filename){

	Data::BEDdataP myBED = std::make_shared<BEDdataT>();
	bool valid_input = true;
	igzstream fin;
	string line;

	fin.open(filename.c_str());
	if (!fin)
		throw range_error("ERROR LoadData: Cannot open index data file: " + filename);

	while (!fin.eof() && valid_input) {
		BEDentryP myEnt = std::make_shared<BEDentryT>();
		if (fin >> myEnt->SEQ >> myEnt->START >> myEnt->END >> myEnt->NAME >> myEnt->SCORE >> myEnt->STRAND){
			cout << "BED " << myEnt->SEQ << "\t" <<  myEnt->START << "\t" <<  myEnt->END << "\t" <<  myEnt->NAME << "\t" <<  myEnt->SCORE << "\t" <<  myEnt->STRAND << endl;

			getline(fin, line);
			istringstream iss(line,istringstream::in);
			string col;
			while (iss >> col) {
				myEnt->COLS.push_back(col);
			}
			myBED->insert(make_pair(myEnt->SEQ,myEnt));
		}
	}
	fin.close();

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

bool Data::SetGraphFromSeq2(GraphClass& oG, string& currSeq, unsigned& pos, bool& lastGr) {

	bool success_status = false;
	lastGr = false;
	unsigned win=mpParameters->mSeqWindow;
	unsigned shift = std::max((double)1,(double)win*mpParameters->mSeqShift);

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
			string seq = currSeq.substr(pos,currSize);
			SetGraphFromSeq( seq ,oG);
			//cout << currSeq.size() << " " << currSize << " eof? " << in.eof() << " pos " << pos << " win " << win << " " << seq << endl;
		} else
			throw range_error("ERROR FASTA reader! Too short sequence found. " + std::to_string(pos));

		if ((win>0) && (currSeq.size()-pos-shift>=win)){
			pos += shift;
		} else if ((currSeq.size()-shift-pos<win ) || (win == 0))
		{
			currSeq="";
			pos = 0;
			lastGr=true;
		}
		success_status = true;
	}
	return success_status;
}


bool Data::SetGraphFromSeq(string& seq, GraphClass& oG) {
	vector<bool> vertex_status(5, false);
	vertex_status[0] = true; //kernel point
	vertex_status[1] = true; //kind
	vertex_status[2] = true; //viewpoint
	vertex_status[3] = false; //dead
	vertex_status[4] = false; //abstraction

	bool success_status = false;
	//cout << graphSeq << " " << currSeq.size() << " " << currSize << " " << in.eof() << endl;
	unsigned real_vertex_index = oG.InsertVertex();
	vector<string> vertex_symbolic_attribute_list(1);
	vertex_symbolic_attribute_list[0] = seq;
	oG.SetVertexSymbolicAttributeList(real_vertex_index, vertex_symbolic_attribute_list);
	oG.SetVertexStatusAttributeList(real_vertex_index, vertex_status);
	success_status = true;

	return success_status;
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


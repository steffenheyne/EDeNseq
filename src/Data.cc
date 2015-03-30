#include "Data.h"

Data::Data() :
mpParameters(0) {
}

void Data::Init(Parameters* apParameters) {
	mpParameters = apParameters;
	//mKernel.Init(mpParameters);
	if (mpParameters->mNumThreads>0 ){
		omp_set_num_threads(mpParameters->mNumThreads);
	} else {
		omp_set_num_threads(std::thread::hardware_concurrency());
	}
}

vector<SeqDataSet> Data::LoadIndexDataList(string filename){

	vector<SeqDataSet> myList;
	bool valid_input = true;
	igzstream fin;
	string line;

	fin.open(filename.c_str());
	if (!fin)
		throw range_error("ERROR LoadData: Cannot open index data file: " + filename);
	while (!fin.eof() && valid_input) {
		SeqDataSet mySet;
		mySet.filetype=FASTA;
		if (fin >> mySet.idx >> mySet.filename >> mySet.desc){
			cout << "found file idx " << mySet.idx << "\t" << mySet.filename << "\t" << mySet.desc << endl;
			mySet.updateIndex=true;
			mySet.updateSigCache=false;
			myList.push_back(mySet);
		}
		getline(fin, line);
	}
	fin.close();
	if (!myList.size())
		throw range_error("ERROR LoadIndexData: No data found in " + filename + "!");
	return myList;
}


void Data::SetGraphFromFile(istream& in, GraphClass& oG) {
	switch (mpParameters->mFileTypeCode) {
	case STRINGSEQ: {
		SetGraphFromStringFile(in, oG);
		break;
	}
	default:
		throw range_error("ERROR Data::SetGraphFromFile: file type not recognized: " + mpParameters->mFileType);
	}
	//******************************************************************************************************
}

bool Data::SetGraphFromFASTAFile(istream& in, GraphClass& oG, string& currSeq) {
	vector<bool> vertex_status(5, false);
	vertex_status[0] = true; //kernel point
	vertex_status[1] = true; //kind
	vertex_status[2] = true; //viewpoint
	vertex_status[3] = false; //dead
	vertex_status[4] = false; //abstraction

	bool success_status = false;

	unsigned win=mpParameters->mSeqWindow;
	unsigned shift = (unsigned)((double)win*mpParameters->mSeqShift);

	if (currSeq.size() == 0){
		in >> std::ws;
	}

	char c = in.peek();
	string header;
	bool newSeq=false;
	//cout << "here "<< " " << currSeq.size() << " " << in.eof() << " c "  << c << endl;
	if (!in.eof() && c != EOF && c=='>' && currSeq.size() == 0 ){
		getline(in, header,'>');
		getline(in, header);
		getline(in, currSeq,'>');

		currSeq.erase(std::remove(currSeq.begin(), currSeq.end(), '\n'),currSeq.end());
		currSeq.erase(std::remove(currSeq.begin(), currSeq.end(), ' '),currSeq.end());
		std::transform(currSeq.begin(), currSeq.end(), currSeq.begin(), ::toupper);
		in.unget();
		if (currSeq.size()==0 || header.size()==0)
			throw range_error("ERROR FASTA reader - empty Sequence or header found! Header:"+header);
		//cout << " found seq " << header << " " << currSeq.size() << " length" << " EOF? "<< in.eof() << endl;
		newSeq=true;
	} else if (c != '>' && c != EOF && c!= '\n') {
		throw range_error("ERROR FASTA format error!");
	}

	if (currSeq.size() > 0 ) {

		// default case for window/shift
		unsigned currSize = win;
		// case now window/shift
		if (win==0){
			currSize = currSeq.size();
		} else if (win>currSeq.size()) {
			// case seq left is smaller than win
			currSize=currSeq.size();
		}

		if (currSize>=win){

			string graphSeq = currSeq.substr(0,currSize);
			//cout << graphSeq << " " << currSeq.size() << " " << currSize << " " << in.eof() << endl;
			unsigned real_vertex_index = oG.InsertVertex();
			vector<string> vertex_symbolic_attribute_list(1);
			vertex_symbolic_attribute_list[0] = graphSeq;
			oG.SetVertexSymbolicAttributeList(real_vertex_index, vertex_symbolic_attribute_list);
			oG.SetVertexStatusAttributeList(real_vertex_index, vertex_status);
		} else if (newSeq){
			throw range_error("ERROR FASTA reader! Too short sequence found. Either use win=0 or increase window size!");
		}

		if (win>0 && currSeq.size()-shift>=win){
			currSeq.erase(0,shift);
		} else if (currSeq.size()-shift<win || win == 0)
			currSeq="";
		success_status = true;
	}
	return success_status;
}

bool Data::SetGraphFromStringFile(istream& in, GraphClass& oG) {
	vector<bool> vertex_status(5, false);
	vertex_status[0] = true; //kernel point
	vertex_status[1] = true; //kind
	vertex_status[2] = true; //viewpoint
	vertex_status[3] = false; //dead
	vertex_status[4] = false; //abstraction

	bool success_status = false;
	string line;
	getline(in, line);
	if (line == "")
		return false;
	unsigned real_vertex_index = oG.InsertVertex();
	vector<string> vertex_symbolic_attribute_list(1);
	vertex_symbolic_attribute_list[0] = line;
	oG.SetVertexSymbolicAttributeList(real_vertex_index, vertex_symbolic_attribute_list);
	oG.SetVertexStatusAttributeList(real_vertex_index, vertex_status);
	success_status = true;

	return success_status;
}


unsigned Data::Size() {
	//return mVectorList.size();
	return mDataSize;
}

void Data::SetDataSize(unsigned aSize){
	mDataSize=aSize;
}

void Data::writeBinaryIndex(ostream &out, const vector<umap_uint_vec_uint> &index) {
	// create binary reverse index representation
	// format:
	unsigned numHashFunc = index.size();
	out.write((const char*) &numHashFunc, sizeof(unsigned));
	for (vector<umap_uint_vec_uint>::const_iterator it = index.begin(); it!= index.end(); it++){
		unsigned numBins = it->size();
		out.write((const char*) &numBins, sizeof(unsigned));
		for (umap_uint_vec_uint::const_iterator itBin = it->begin(); itBin!=it->end(); itBin++){
			unsigned binId = itBin->first;
			unsigned numBinEntries = itBin->second.size();
			out.write((const char*) &binId, sizeof(unsigned));
			out.write((const char*) &numBinEntries, sizeof(unsigned));
			for (vector<unsigned>::const_iterator binEntry = itBin->second.begin(); binEntry != itBin->second.end(); binEntry++){
				unsigned short int t= *binEntry;
				out.write((const char*) &(t), sizeof(unsigned short int));
			}
		}
	}
}

bool Data::readBinaryIndex(string filename, vector<umap_uint_vec_uint> &index){
	igzstream fin;
	fin.open(filename.c_str());
	unsigned numHashFunc = 0;
	fin.read((char*) &numHashFunc, sizeof(unsigned));
	if (numHashFunc <= 0)
		fin.setstate(std::ios::badbit);
	if (!fin.good())
		return false;
	//cout << "here1"<< endl;
	index.resize(numHashFunc);
	for (unsigned  hashFunc = 0; hashFunc < numHashFunc; hashFunc++){

		unsigned numBins = 0;
		fin.read((char*) &numBins, sizeof(unsigned));
		if (numBins < 0)
			fin.setstate(std::ios::badbit);
		if (!fin.good())
			return false;
		//cout << "here2"<< endl;
		for (unsigned  bin = 0; bin < numBins; bin++){

			unsigned binId = 0;
			unsigned numBinEntries = 0;
			fin.read((char*) &binId, sizeof(unsigned));
			fin.read((char*) &numBinEntries, sizeof(unsigned));
			if (numBinEntries < 0 || binId <= 0)
				fin.setstate(std::ios::badbit);
			if (!fin.good())
				return false;
			//cout << "here3"<< endl;
			vector<unsigned> tmp(numBinEntries);
			index[hashFunc].insert(make_pair(binId,tmp));

			for (unsigned entry = 0; entry < numBinEntries; entry++ ){

				unsigned short int t = 0;
				fin.read((char*) &t, sizeof(unsigned short int));
				if (t < 0)
					fin.setstate(std::ios::badbit);
				if (!fin.good())
					return false;
				//cout << "here4"<< endl;
				index[hashFunc][binId].push_back(t);
			}
		}
	}
	fin.close();
	return true;
}

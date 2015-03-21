#include "Data.h"

Data::Data() :
mpParameters(0), mDataIsLoaded(false) {
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

void Data::LoadMultilabelRealList(string aFileName, vector<vector<double> >& oList) {
	oList.clear();
	if (mpParameters->mVerbose)
		cout << endl << "Reading multilabel file: " << aFileName << " ...";
	ifstream fin;
	fin.open(aFileName.c_str());
	if (!fin)
		throw range_error("ERROR Data::LoadMultilabelRealList: Cannot open file:" + aFileName);
	while (!fin.eof()) {
		string line;
		getline(fin, line);
		stringstream ss;
		ss << line << endl;
		if (line.size() > 1) {
			vector<double> multilabel_list;
			while (!ss.eof()) {
				double value(0);
				ss >> value;
				if (ss.good()) {
					multilabel_list.push_back(value);
				}
			}
			oList.push_back(multilabel_list);
		}
	}
	fin.close();
}


void Data::LoadRealList(string aFileName, vector<double>& oList) {
	oList.clear();
	if (mpParameters->mVerbose)
		cout << endl << "Reading file: " << aFileName << " ..";
	ifstream fin;
	fin.open(aFileName.c_str());
	if (!fin)
		throw range_error("ERROR Data::LoadRealList: Cannot open file:" + aFileName);
	while (!fin.eof()) {
		string line;
		getline(fin, line);
		stringstream ss;
		ss << line << endl;
		while (!ss.eof()) {
			double value(0);
			ss >> value;
			if (ss.good()) {
				oList.push_back(value);
			}
		}
	}
	fin.close();
	if (mpParameters->mVerbose)
		cout << ".. read: " << oList.size() << " values." << endl;
}

void Data::LoadUnsignedList(string aFileName, vector<unsigned>& oList) {
	oList.clear();
	if (!mpParameters->mVerbose)
		cout << endl << "Reading file: " << aFileName << " ..";
	ifstream fin;
	fin.open(aFileName.c_str());
	if (!fin)
		throw range_error("ERROR Data::LoadUnsignedList: Cannot open file:" + aFileName);
	while (!fin.eof()) {
		string line;
		getline(fin, line);
		stringstream ss;
		ss << line << endl;
		while (!ss.eof()) {
			unsigned value(0);
			ss >> value;
			if (ss.good()) {
				oList.push_back(value);
			}
		}
	}
	fin.close();
	if (mpParameters->mVerbose)
		cout << ".. read: " << oList.size() << " values." << endl;
}

void Data::LoadData(istream& fin) {
	ProgressBar pb;
	bool valid_input = true;
	unsigned instance_counter = 0;
	while (!fin.eof() && valid_input) {
		switch (mpParameters->mFileTypeCode) {
		case STRINGSEQ:{
			vector<GraphClass> g_list(BUFFER_SIZE);
			unsigned i = 0;
			while (i < BUFFER_SIZE && !fin.eof() && valid_input) {
				SetGraphFromFile(fin, g_list[i]);
				if (g_list[i].IsEmpty()) {
					valid_input = false;
				} else {
					i++;
					instance_counter++;
				}
			}
		}
		break;
		default:
			throw range_error("ERROR Data::LoadData: file type not recognized: " + mpParameters->mFileType);
		}
	}

	if (instance_counter > 0)
		mDataIsLoaded = true;
	else
		throw range_error("ERROR Data::LoadData: something went wrong: data file was expected but no data is available");

	SetDataSize(instance_counter);
}

void Data::LoadData(bool aLoadIndex = false, bool aLoadTarget = false, bool aLoadMultilabelTarget = false) {
	mVectorList.clear();
	//	mKernel.ParametersSetup();
	igzstream fin;
	fin.open(mpParameters->mInputDataFileName.c_str());
	if (!fin)
		throw range_error("ERROR Data::LoadData: Cannot open file: " + mpParameters->mInputDataFileName);
	if (mpParameters->mVerbose)
		cout << "Processing file: " << mpParameters->mInputDataFileName << endl;
	LoadData(fin);
}

bool Data::IsDataLoaded() {
	return mDataIsLoaded;
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


#include "Data.h"

Data::Data() :
		mpParameters(0), mDataIsLoaded(false), mTargetIsLoaded(false), mIndexIsLoaded(false) {
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


void Data::LoadIndex() {
	if (mpParameters->mRowIndexFileName != "")
		if (mRowIndexList.size() == 0) {
			LoadUnsignedList(mpParameters->mRowIndexFileName, mRowIndexList);
			mIndexIsLoaded = true;
		}
	if (mpParameters->mColIndexFileName != "")
		if (mColIndexList.size() == 0) {
			LoadUnsignedList(mpParameters->mColIndexFileName, mColIndexList);
			mIndexIsLoaded = true;
		}

	if (IsDataLoaded() == false)
		throw range_error("ERROR Data::LoadIndex: Cannot assign indices before reading data file.");

	if (mRowIndexList.size() == 0) {
		if (!mpParameters->mMinimalOutput)
			cout << endl << "No row index list specified. Assuming all " << Size() << " row indices as valid." << endl;
		for (unsigned i = 0; i < Size(); ++i)
			mRowIndexList.push_back(i);
	}
	if (mColIndexList.size() == 0) {
		if (!mpParameters->mMinimalOutput)
			cout << endl << "No col index list specified. Assuming all " << Size() << " col indices as valid." << endl;
		for (unsigned i = 0; i < Size(); ++i)
			mColIndexList.push_back(i);
	}
}

bool Data::IsIndexLoaded() {
	return mIndexIsLoaded;
}

void Data::LoadMultilabelTarget() {
	if (mpParameters->mTargetFileName != "") {
		mTargetList.clear();
		LoadMultilabelRealList(mpParameters->mTargetFileName, mMultilabelTargetList);
		mMultilabelTargetIsLoaded = true;
	}
}

bool Data::IsMultilabelTargetLoaded() {
	return mMultilabelTargetIsLoaded;
}

void Data::LoadMultilabelRealList(string aFileName, vector<vector<double> >& oList) {
	oList.clear();
	if (!mpParameters->mMinimalOutput)
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

	//check all targets with same dimension
	unsigned dimension = MultilabelTargetDimension();
	for (unsigned i = 0; i < mMultilabelTargetList.size(); ++i) {
		if (mMultilabelTargetList[i].size() != dimension)
			throw range_error("ERROR Data::LoadMultilabelRealList: target row:" + stream_cast<string>(i) + " has dimension: " + stream_cast<string>(mMultilabelTargetList[i].size()) + " instead of: " + stream_cast<string>(dimension));
	}

	if (!mpParameters->mMinimalOutput)
		cout << "... read: " << oList.size() << " vectors each one of size " << MultilabelTargetDimension() << "." << endl;
}

unsigned Data::MultilabelTargetDimension() {
	return mMultilabelTargetList[0].size();
}

void Data::LoadTarget() {
	if (mpParameters->mTargetFileName != "") {
		mTargetList.clear();
		LoadRealList(mpParameters->mTargetFileName, mTargetList);
		mTargetIsLoaded = true;
	}
}

bool Data::IsTargetLoaded() {
	return mTargetIsLoaded;
}

void Data::LoadRealList(string aFileName, vector<double>& oList) {
	oList.clear();
	if (!mpParameters->mMinimalOutput)
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
	if (!mpParameters->mMinimalOutput)
		cout << ".. read: " << oList.size() << " values." << endl;
}

void Data::LoadUnsignedList(string aFileName, vector<unsigned>& oList) {
	oList.clear();
	if (!mpParameters->mMinimalOutput)
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
	if (!mpParameters->mMinimalOutput)
		cout << ".. read: " << oList.size() << " values." << endl;
}

void Data::LoadData(istream& fin) {
	ProgressBar pb;
	bool valid_input = true;
	unsigned instance_counter = 0;
	while (!fin.eof() && valid_input) {
		switch (mpParameters->mFileTypeCode) {
			case GRAPH:
			case STRINGSEQ:
			case SEQUENCE: {
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
			case SPARSE_VECTOR: {
				bool success_status = false;
				SVector x(pow(2, mpParameters->mHashBitSize));
				if (mpParameters->mBinaryFormat)
					success_status = SetVectorFromSparseVectorBinaryFile(fin, x);
				else
					success_status = SetVectorFromSparseVectorAsciiFile(fin, x);
				if (success_status) {
					mVectorList.push_back(x);
					instance_counter++;
					if (!mpParameters->mMinimalOutput)
						pb.Count();
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
	if (!mpParameters->mMinimalOutput)
		cout << "Processing file: " << mpParameters->mInputDataFileName << endl;
	LoadData(fin);

	//additional files to load if required
	if (IsIndexLoaded() == false && aLoadIndex == true)
		LoadIndex();

	if (IsTargetLoaded() == false && aLoadTarget == true)
		LoadTarget();

	if (IsMultilabelTargetLoaded() == false && aLoadMultilabelTarget == true)
		LoadMultilabelTarget();

}

bool Data::IsDataLoaded() {
	return mDataIsLoaded;
}

void Data::SetGraphFromFile(istream& in, GraphClass& oG) {
	switch (mpParameters->mFileTypeCode) {
		case GRAPH:
			SetGraphFromGraphGspanFile(in, oG);
			break;
		case SEQUENCE: {
			if (mpParameters->mSequenceToken && mpParameters->mSequenceMultiLine)
				SetGraphFromSequenceMultiLineTokenFile(in, oG);
			if (mpParameters->mSequenceToken && !mpParameters->mSequenceMultiLine)
				SetGraphFromSequenceTokenFile(in, oG);
			if (!mpParameters->mSequenceToken && mpParameters->mSequenceMultiLine)
				SetGraphFromSequenceMultiLineFile(in, oG);
			if (!mpParameters->mSequenceToken && !mpParameters->mSequenceMultiLine)
				SetGraphFromSequenceFile(in, oG);
			break;
		}
		case STRINGSEQ: {
			SetGraphFromStringFile(in, oG);
			break;
		}
		case FASTA: {
			SetGraphFromFASTAFile(in, oG);
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

bool Data::SetGraphFromFASTAFile(istream& in, GraphClass& oG) {
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

bool Data::SetVectorFromSparseVectorAsciiFile(istream& in, SVector& aX) {
	bool success_status = false;
	string line;
	getline(in, line);
	if (line == "")
		return false;
	stringstream ss;
	ss << line << endl;
	while (!ss.eof() && ss.good()) {
		string key_value;
		ss >> key_value;
		size_t limit = key_value.find_first_of(":", 0);
		if (limit != string::npos) { //if the delimiter ':' is found then proceed
			string key = key_value.substr(0, limit);
			string value = key_value.substr(limit + 1, key_value.size());
			unsigned key_int = stream_cast<unsigned>(key);
			double val_real = stream_cast<double>(value);
			aX.insert(key_int) = val_real;
			success_status = true;
		}
	}
	return success_status;
}

bool Data::SetVectorFromSparseVectorBinaryFile(istream& in, SVector& aX) {
	return loadSparseBinary(in, aX);
}

void Data::SetGraphFromSequenceFile(istream& in, GraphClass& oG) {
	vector<vector<unsigned> > vertex_component_list;
	vector<unsigned> vertex_component;
	bool graph_disconnect = true;
	string line;
	getline(in, line);
	if (line == "")
		return;
	unsigned vertex_counter = 0;
	for (unsigned position_counter = 0; position_counter < line.length(); position_counter++) {
		char char_label = line.at(position_counter);
		string label = stream_cast<string>(char_label);
		if (label == "|") {
			graph_disconnect = true;
			vertex_component_list.push_back(vertex_component);
			vertex_component.clear();
		} else {
			vertex_component.push_back(vertex_counter);
			AddVertexAndEdgesForSequence(oG, label, vertex_counter, graph_disconnect);
			vertex_counter++;
			graph_disconnect = false;
		}
	}
	vertex_component_list.push_back(vertex_component);
	ManageDirectedAndMultiComponent(oG, vertex_component_list);
}

void Data::SetGraphFromSequenceTokenFile(istream& in, GraphClass& oG) {
	vector<vector<unsigned> > vertex_component_list;
	vector<unsigned> vertex_component;
	bool graph_disconnect = true;
	string line;
	getline(in, line);
	if (line == "")
		return;
	stringstream ss;
	ss << line << endl;
	unsigned vertex_counter = 0;
	while (!ss.eof() && ss.good()) {
		//add vertex
		string label;
		ss >> label;
		if (label != "") {
			if (label == "|") {
				graph_disconnect = true;
				vertex_component_list.push_back(vertex_component);
				vertex_component.clear();
			} else {
				vertex_component.push_back(vertex_counter);
				AddVertexAndEdgesForSequence(oG, label, vertex_counter, graph_disconnect);
				vertex_counter++;
				graph_disconnect = false;
			}
		}
	}
	vertex_component_list.push_back(vertex_component);
	ManageDirectedAndMultiComponent(oG, vertex_component_list);
}

void Data::SetGraphFromSequenceMultiLineFile(istream& in, GraphClass& oG) {
	vector<vector<unsigned> > vertex_component_list;
	vector<unsigned> vertex_component;
	bool graph_disconnect = true;
	vector<string> line(mpParameters->mSequenceDegree);
	for (unsigned i = 0; i < mpParameters->mSequenceDegree; ++i) {
		getline(in, line[i]);
	}
	unsigned sequence_length = line[0].length();
	unsigned vertex_counter = 0;
	for (unsigned j = 0; j < sequence_length; j++) {
		for (unsigned i = 0; i < mpParameters->mSequenceDegree; ++i) {
			if (line[i] == "")
				return;
			char char_label = line[i].at(j);
			string label = stream_cast<string>(char_label);
			if (label == "|") {
				if (i == 0) { //process the first | symbol
					graph_disconnect = true;
					vertex_component_list.push_back(vertex_component);
					vertex_component.clear();
				} else { //skip the remaining | symbols
					//do nothing
				}
			} else {
				vertex_component.push_back(vertex_counter);
				AddVertexAndEdgesForSequence(oG, label, vertex_counter, graph_disconnect);
				vertex_counter++;
				graph_disconnect = false;
			}
		}
	}
	//at last get the last vertex_component in
	vertex_component_list.push_back(vertex_component);
	ManageDirectedAndMultiComponent(oG, vertex_component_list);
}

void Data::SetGraphFromSequenceMultiLineTokenFile(istream& in, GraphClass& oG) {
	vector<vector<unsigned> > vertex_component_list;
	vector<unsigned> vertex_component;
	bool graph_disconnect = true;
	vector<vector<string> > multi_line;
	for (unsigned i = 0; i < mpParameters->mSequenceDegree; ++i) {
		vector<string> single_line;
		string line;
		getline(in, line);
		stringstream ss;
		ss << line << endl;
		while (!ss.eof() && ss.good()) {
			string label;
			ss >> label;
			if (label != "") {
				single_line.push_back(label);
			}
		}
		if (single_line.size() > 0) {
			multi_line.push_back(single_line);
		}
	}

	if (multi_line.size() > 0) {
		unsigned sequence_length = multi_line[0].size();
		for (unsigned i = 0; i < mpParameters->mSequenceDegree; ++i)
			assert(multi_line[i].size() == sequence_length);

		unsigned vertex_counter = 0;
		for (unsigned j = 0; j < sequence_length; j++) {
			for (unsigned i = 0; i < mpParameters->mSequenceDegree; ++i) {
				string label = multi_line[i][j];
				if (label == "|") {
					if (i == 0) { //process the first | symbol
						graph_disconnect = true;
						vertex_component_list.push_back(vertex_component);
						vertex_component.clear();
					} else { //skip the remaining | symbols
						//do nothing
					}
				} else {
					vertex_component.push_back(vertex_counter);
					AddVertexAndEdgesForSequence(oG, label, vertex_counter, graph_disconnect);
					vertex_counter++;
					graph_disconnect = false;
				}
			}
		}
		//at last get the last vertex_component in
		vertex_component_list.push_back(vertex_component);
		ManageDirectedAndMultiComponent(oG, vertex_component_list);
	}
}

void Data::ManageDirectedAndMultiComponent(GraphClass& oG, vector<vector<unsigned> >& aVertexComponentList) {
	if (mpParameters->mGraphType == "DIRECTED" && aVertexComponentList.size() > 1) {
		//create the components for the reversed direction graph
		vector<vector<unsigned> > reverse_vertex_component_list;
		for (unsigned u = 0; u < aVertexComponentList.size(); ++u) {
			vector<unsigned> reverse_vertex_component;
			for (unsigned t = 0; t < aVertexComponentList[u].size(); ++t) {
				unsigned id = aVertexComponentList[u][t];
				reverse_vertex_component.push_back(id + oG.VertexSize());
			}
			reverse_vertex_component_list.push_back(reverse_vertex_component);
		}
		//add the reverse direction components to the component list
		aVertexComponentList.insert(aVertexComponentList.end(), reverse_vertex_component_list.begin(), reverse_vertex_component_list.end());
	}

	if (mpParameters->mGraphType == "DIRECTED")
		AddReverseGraph(oG);

	if (mpParameters->mSequencePairwiseInteraction)
		if (aVertexComponentList.size() > 1)
			AddAbstractConnections(oG, aVertexComponentList);
}

void Data::AddVertexAndEdgesForSequence(GraphClass& oG, string aLabel, unsigned aVertexCounter, bool aGraphDisconnect) {
	//set (once) boolean status vectors
	static vector<bool> vertex_status(5, false);
	vertex_status[0] = true; //kernel point
	vertex_status[1] = true; //kind
	vertex_status[2] = true; //viewpoint
	//vertex_status[3] = false; //dead
	//vertex_status[4] = false; //abstraction

	static vector<bool> edge_status(3, false);
	//edge_status[0] = false; //edge dead
	//edge_status[1] = false; //edge abstraction_of
	//edge_status[2] = false; //edge part_of

	unsigned real_vertex_index = oG.InsertVertex();
	vector<string> vertex_symbolic_attribute_list;
	vertex_symbolic_attribute_list.push_back(aLabel);
	oG.SetVertexSymbolicAttributeList(real_vertex_index, vertex_symbolic_attribute_list);
	oG.SetVertexStatusAttributeList(real_vertex_index, vertex_status);
	if (aVertexCounter % mpParameters->mSequenceDegree != 0) {
		oG.SetVertexViewPoint(real_vertex_index, false);
	}
	if (aGraphDisconnect == false) {
		//add edge
		unsigned real_src_index;
		unsigned real_dest_index;
		if (aVertexCounter % mpParameters->mSequenceDegree != 0) { //add edge to preceding vertex
			real_src_index = aVertexCounter;
			real_dest_index = aVertexCounter - 1;
		} else { //vertices with position multiple than mSequenceDegree are connected in sequence
			real_src_index = aVertexCounter;
			real_dest_index = aVertexCounter - mpParameters->mSequenceDegree;
		}
		vector<string> edge_symbolic_attribute_list;
		edge_symbolic_attribute_list.push_back("-"); //NOTE: default edge label is '-'
		unsigned edge_index = oG.InsertEdge(real_src_index, real_dest_index);
		oG.SetEdgeSymbolicAttributeList(edge_index, edge_symbolic_attribute_list);
		oG.SetEdgeStatusAttributeList(edge_index, edge_status);
		if (mpParameters->mGraphType == "UNDIRECTED" || aVertexCounter % mpParameters->mSequenceDegree != 0) { //add reverse edge in case the graph is undirected or for the attribute vertices, so that in both directions these are accessible
			unsigned reverse_edge_index = oG.InsertEdge(real_dest_index, real_src_index);
			oG.SetEdgeSymbolicAttributeList(reverse_edge_index, oG.GetEdgeSymbolicAttributeList(edge_index));
			oG.SetEdgeStatusAttributeList(reverse_edge_index, oG.GetEdgeStatusAttributeList(edge_index));
		}
	}
}

void Data::AddAbstractConnections(GraphClass& oG, vector<vector<unsigned> >& aVertexComponentList) {
	//set (once) boolean status vectors
	static vector<bool> vertex_status(5, false);
	//vertex_status[0] = false; //kernel point
	//vertex_status[1] = false; //kind
	//vertex_status[2] = false; //viewpoint
	//vertex_status[3] = false; //dead
	vertex_status[4] = true; //abstraction

	static vector<bool> edge_status(3, false);
	//edge_status[0] = false; //edge dead
	//edge_status[1] = false; //edge abstraction_of
	//edge_status[2] = false; //edge part_of

	//for all pairs of components
	for (unsigned i = 0; i < aVertexComponentList.size(); ++i) {
		for (unsigned j = i + 1; j < aVertexComponentList.size(); ++j) {
			//join all vertices in one component to all other vertices in the other component
			//add 1 abstract vertex
			unsigned real_vertex_index = oG.InsertVertex();
			vector<string> vertex_symbolic_attribute_list;
			vertex_symbolic_attribute_list.push_back("^L");
			oG.SetVertexSymbolicAttributeList(real_vertex_index, vertex_symbolic_attribute_list);
			oG.SetVertexStatusAttributeList(real_vertex_index, vertex_status);

			for (unsigned ii = 0; ii < aVertexComponentList[i].size(); ii++) { //add part_of edges
				unsigned real_src_index = real_vertex_index;
				unsigned real_dest_index = aVertexComponentList[i][ii];
				//only add edges to point vertices
				if (oG.GetVertexViewPoint(real_dest_index)) {
					vector<string> edge_symbolic_attribute_list;
					edge_symbolic_attribute_list.push_back("@-");
					unsigned edge_index = oG.InsertEdge(real_src_index, real_dest_index);
					oG.SetEdgeSymbolicAttributeList(edge_index, edge_symbolic_attribute_list);
					oG.SetEdgeStatusAttributeList(edge_index, edge_status);
					oG.SetEdgePartOf(edge_index, true);
				}
			}

			for (unsigned jj = 0; jj < aVertexComponentList[j].size(); jj++) { //add abstraction_of edges
				unsigned real_src_index = real_vertex_index;
				unsigned real_dest_index = aVertexComponentList[j][jj];
				if (oG.GetVertexViewPoint(real_dest_index)) {
					vector<string> edge_symbolic_attribute_list;
					edge_symbolic_attribute_list.push_back("^-");
					unsigned edge_index = oG.InsertEdge(real_src_index, real_dest_index);
					oG.SetEdgeSymbolicAttributeList(edge_index, edge_symbolic_attribute_list);
					oG.SetEdgeStatusAttributeList(edge_index, edge_status);
					oG.SetEdgeAbstractionOf(edge_index, true);
				}
			}
		}
	}
}

void Data::SetGraphFromGraphGspanFile(istream& in, GraphClass& oG) {
	//status
	vector<bool> vertex_status(5, false);
	vertex_status[0] = true; //kernel point
	vertex_status[1] = true; //kind
	vertex_status[2] = true; //viewpoint
	vertex_status[3] = false; //dead
	vertex_status[4] = false; //abstraction

	vector<bool> edge_status(3, false);
	edge_status[0] = false; //edge dead
	edge_status[1] = false; //edge abstraction_of
	edge_status[2] = false; //edge part_of

	map<string, int> index_map_nominal_to_real;
	//string gid;
	unsigned line_counter = 0;
	bool instance_started = false;
	do {
		char c = in.peek();
		if (c == 't' && instance_started == true)
			break;

		line_counter++;
		string line;
		getline(in, line);
		if (line == "")
			break;
		stringstream ss;
		ss << line << endl;
		char code;
		ss >> code;
		if (code == 't') {
			instance_started = true;
		} else if (code == 'v' || code == 'V' || code == 'W') {
			//extract vertex id and make map nominal_id -> real_id
			string nominal_vertex_index;
			ss >> nominal_vertex_index;
			unsigned real_vertex_index = oG.InsertVertex();
			index_map_nominal_to_real[nominal_vertex_index] = real_vertex_index;
			//label
			vector<string> vertex_symbolic_attribute_list;
			string label;
			ss >> label;
			if (label == "")
				throw range_error("ERROR Data::SetGraphFromGraphGspanFile: vertex label doies not exist for vertex with index: " + nominal_vertex_index + " at line: " + stream_cast<string>(line_counter));
			vertex_symbolic_attribute_list.push_back(label);
			oG.SetVertexSymbolicAttributeList(real_vertex_index, vertex_symbolic_attribute_list);
			oG.SetVertexStatusAttributeList(real_vertex_index, vertex_status);
			if (code == 'V') //non-viewpoint but kernel point
				oG.SetVertexViewPoint(real_vertex_index, false);
			if (code == 'W') { //non kernel point
				oG.SetVertexViewPoint(real_vertex_index, false);
				oG.SetVertexKernelPoint(real_vertex_index, false);
			}
			char vertex_abstraction_code = label.at(0);
			if (vertex_abstraction_code == '^') {
				oG.SetVertexAbstraction(real_vertex_index, true);
				oG.SetVertexViewPoint(real_vertex_index, false);
				oG.SetVertexKernelPoint(real_vertex_index, false);
			}

			std::size_t found = line.find("vec:");
			if (mpParameters->mUseRealVectorInformation && found != std::string::npos) { //check that the key string "vec" is found in the line
			//advance until the string "vec:" is found, then read in all values as real components
				string var = "";
				string vec_code = "vec:";
				do {
					ss >> var;
				} while (var != vec_code);
				vector<double> vec;
				while (!ss.fail()) {
					double val;
					ss >> val;
					if (!ss.fail())
						vec.push_back(val);
				}
				oG.SetVertexNumericAttributeList(real_vertex_index, vec);
			}
		} else if (code == 'e') {
			//extract src and dest vertex id
			string nominal_src_index, nominal_dest_index;
			string label;
			ss >> nominal_src_index >> nominal_dest_index >> label;
			if (index_map_nominal_to_real.count(nominal_src_index) == 0)
				throw range_error("ERROR Data::SetGraphFromGraphGspanFile: Did not find nominal src index: " + nominal_src_index + " at line: " + stream_cast<string>(line_counter));
			if (index_map_nominal_to_real.count(nominal_dest_index) == 0)
				throw range_error("ERROR Data::SetGraphFromGraphGspanFile: Did not find nominal dest index: " + nominal_dest_index + " at line: " + stream_cast<string>(line_counter));
			if (label == "")
				throw range_error("ERROR Data::SetGraphFromGraphGspanFile: edge label does not exist for edge with src: " + nominal_src_index + " and dest: " + nominal_dest_index + " at line: " + stream_cast<string>(line_counter));
			vector<string> edge_symbolic_attribute_list;
			edge_symbolic_attribute_list.push_back(label);
			unsigned real_src_index = index_map_nominal_to_real[nominal_src_index];
			unsigned real_dest_index = index_map_nominal_to_real[nominal_dest_index];
			unsigned edge_index = oG.InsertEdge(real_src_index, real_dest_index);
			oG.SetEdgeSymbolicAttributeList(edge_index, edge_symbolic_attribute_list);
			oG.SetEdgeStatusAttributeList(edge_index, edge_status);

			char edge_abstraction_code = label.at(0);
			if (edge_abstraction_code == '^')
				oG.SetEdgeAbstractionOf(edge_index, true);
			if (edge_abstraction_code == '@')
				oG.SetEdgePartOf(edge_index, true);

			if (mpParameters->mGraphType == "UNDIRECTED" || edge_abstraction_code == '^' || edge_abstraction_code == '@') { //NOTE: edges that are part of the abstraction mechanism should be treated as undirected
				unsigned reverse_edge_index = oG.InsertEdge(real_dest_index, real_src_index);
				oG.SetEdgeSymbolicAttributeList(reverse_edge_index, oG.GetEdgeSymbolicAttributeList(edge_index));
				oG.SetEdgeStatusAttributeList(reverse_edge_index, oG.GetEdgeStatusAttributeList(edge_index));
			}
		} else {
		} //NOTE: ignore other markers
	} while (!in.eof() && in.good());
	if (mpParameters->mGraphType == "DIRECTED")
		AddReverseGraph(oG);
}

void Data::AddReverseGraph(GraphClass& oG) {
	unsigned vsize = oG.VertexSize();
	//add a copy of all vertices
	for (unsigned i = 0; i < vsize; i++) {
		unsigned real_vertex_index = oG.InsertVertex();
		assert(real_vertex_index == i + vsize);
		vector<string> vertex_symbolic_attribute_list = oG.GetVertexSymbolicAttributeList(i);
		for (unsigned t = 0; t < vertex_symbolic_attribute_list.size(); t++) //prepend a prefix to mark the reverse direction
			vertex_symbolic_attribute_list[t] = "r." + vertex_symbolic_attribute_list[t];
		oG.SetVertexSymbolicAttributeList(real_vertex_index, vertex_symbolic_attribute_list);
		oG.SetVertexNumericAttributeList(real_vertex_index, oG.GetVertexNumericAttributeList(i)); //assign original numerical vector
		oG.SetVertexStatusAttributeList(real_vertex_index, oG.GetVertexStatusAttributeList(i)); //assign original status vector
	}
	//copy all edges swapping src with dest
	for (unsigned i = 0; i < vsize; i++) {
		//get all edges
		vector<unsigned> adj = oG.GetVertexAdjacentList(i);
		for (unsigned j = 0; j < adj.size(); j++) {
			unsigned orig_src = i;
			unsigned orig_dest = adj[j];
			unsigned reverse_src = orig_dest + vsize;
			unsigned reverse_dest = orig_src + vsize;
			unsigned edge_index = oG.InsertEdge(reverse_src, reverse_dest);
			oG.SetEdgeSymbolicAttributeList(edge_index, oG.GetEdgeSymbolicAttributeList(orig_src, orig_dest));
			oG.SetEdgeStatusAttributeList(edge_index, oG.GetEdgeStatusAttributeList(orig_src, orig_dest));
		}
	}
}

unsigned Data::Size() {
	//return mVectorList.size();
	return mDataSize;
}

void Data::SetDataSize(unsigned aSize){
	mDataSize=aSize;
}

void Data::RandomPartition(double aPartitionRatio, vector<unsigned>& oFirstPart, vector<unsigned>& oSecondPart) {
	oFirstPart.clear();
	oSecondPart.clear();
	vector<unsigned> shuffled_index_list;
	MakeShuffledDataIndicesList(shuffled_index_list, Size());
	for (unsigned i = 0; i < Size(); ++i) {
		if (i <= (unsigned) (Size() * aPartitionRatio))
			oFirstPart.push_back(shuffled_index_list[i]);
		else
			oSecondPart.push_back(shuffled_index_list[i]);
	}
}

void Data::RandomPartition(double aPartitionRatio, vector<bool>& oBitVector) {
	if (aPartitionRatio < 1e-8)
		throw range_error("ERROR Data::RandomPartition: aPartitionRatio value is suspiciously low: " + stream_cast<string>(aPartitionRatio));
	oBitVector = vector<bool>(Size(), false);
	//shuffle the index list
	vector<unsigned> shuffled_index_list;
	MakeShuffledDataIndicesList(shuffled_index_list, Size());
	//mark as true the bit corresponding to the index selected in the shuffling
	unsigned part = (unsigned) ((double) shuffled_index_list.size() * aPartitionRatio);
	for (unsigned i = 0; i < part && i < shuffled_index_list.size(); ++i)
		oBitVector[shuffled_index_list[i]] = true;
}

void Data::Bootstrap(vector<unsigned>& oBootstrappedSample, vector<unsigned>& oOutOfBagSample) {
	if (!IsDataLoaded())
		throw range_error("ERROR Data::GetBootstrappedSampleIDList: IsDataLoaded() failed.");

	//sample with replacement Size() times from data index id range
	set<unsigned> sample_set;
	for (unsigned i = 0; i < Size(); i++) {
		unsigned id = randomUnsigned(Size());
		oBootstrappedSample.push_back(id);
		sample_set.insert(id);
	}

	//compute the out of bag sample
	for (unsigned i = 0; i < Size(); i++) {
		if (sample_set.count(i) == 0)
			oOutOfBagSample.push_back(i);
	}
}

void Data::GenerateRandomSparseInstance(SVector& oX, const SVector& aRef, double aStdDev) {
	const double GUARD_LIMIT = 10000;
	unsigned size = pow(2, mpParameters->mHashBitSize);
	unsigned nnz_size = aRef.nonZeros();
	if (nnz_size == 0)
		throw range_error("ERROR Data::GenerateRandomSparseInstance: size of  reference vector is null.");
	unsigned guard_reject = 0;
	double rel;
	do {
		unsigned guard_dim = 0;
		do {
			double value = randomn();
			int index = randomUnsigned(size);
			oX.coeffRef(index) = value;
			guard_dim++;
		} while ((unsigned) oX.nonZeros() < nnz_size && guard_dim < GUARD_LIMIT); //continue to add random features until the same sparsity coefficient is reached
		oX = oX / oX.norm();
		rel = oX.dot(aRef);
		guard_reject++;
	} while (random01() > rel && guard_reject < GUARD_LIMIT); //reject instances if they are not parallel to the reference vector

	oX = oX * aStdDev;
	oX += aRef;
	oX = oX / oX.norm();
}

void Data::GenerateRandomInstance(SVector& oX, const SVector& aRef, double aStdDev) {
	oX = aRef;
	unsigned size = pow(2, mpParameters->mHashBitSize);
	for (unsigned i = 0; i < size; ++i)
		oX.coeffRef(i) += randomn() / aStdDev;
	oX = oX / oX.norm();
}

void Data::GenerateRandomInstance(SVector& oX) {
	const bool SPARSE_FLAG = true;

	unsigned size = pow(2, mpParameters->mHashBitSize);
	if (SPARSE_FLAG) {
		unsigned id = randomUnsigned(Size());
		//sample a sparsity ratio, i.e. the size of the sparse vector
		unsigned nnz_size = mVectorList[id].nonZeros();
		if (nnz_size == 0)
			throw range_error("ERROR Data::GenerateRandomPermutationInstance: size of  mVectorList[id] is null.");
		unsigned guard = 0;
		do {
			double value = randomn();
			int index = randomUnsigned(size);
			oX.coeffRef(index) = value;
			guard++;
		} while ((unsigned) oX.nonZeros() < nnz_size && guard < 10000);
	} else {
		for (unsigned i = 0; i < size; ++i)
			oX.insert(i) = randomn();
	}
	oX = oX / oX.norm();
}

void Data::GenerateRandomDataset(unsigned aSize, vector<SVector>& oDataset) {
	for (unsigned i = 0; i < aSize; i++) {
		SVector x(pow(2, mpParameters->mHashBitSize));
		GenerateRandomInstance(x);
		oDataset.push_back(x);
	}
}

void Data::GenerateRandomPermutationInstance(SVector& oX) {
	unsigned id = randomUnsigned(Size());
	//sample a sparsity ratio, i.e. the size of the sparse vector
	unsigned size = mVectorList[id].nonZeros();
	if (size == 0)
		throw range_error("ERROR Data::GenerateRandomPermutationInstance: size of  mVectorList[id] is null.");
	unsigned guard = 0;
	do {
		//sample an instance at random
		unsigned idx = randomUnsigned(Size());
		SVector& z = mVectorList[idx];
		//sample a feature at random
		unsigned z_size = z.nonZeros();
		unsigned id_z = randomUnsigned(z_size);
		SVector::InnerIterator it(z);
		for (unsigned counter = 0; it && counter < id_z; ++counter) {
			++it;
		}
		//add it to the current instance
		oX.coeffRef(it.index()) = it.value();
		guard++;
	} while ((unsigned) oX.nonZeros() < size && guard < 10000);
	oX = oX / oX.norm();
}

void Data::GenerateRandomPermutationDataset(unsigned aSize, vector<SVector>& oDataset) {
	for (unsigned i = 0; i < aSize; i++) {
		SVector x(pow(2, mpParameters->mHashBitSize));
		GenerateRandomPermutationInstance(x);
		oDataset.push_back(x);
	}
}


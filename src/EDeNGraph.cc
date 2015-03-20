#include "EDeNGraph.h"

void EDeN::Graph::GetNeighbors(unsigned aSrcVertexID, unsigned aDistance, list<unsigned>& oNeighborsList) const {
	if (mNeighborsCache[aSrcVertexID][aDistance].size() == 0)
		throw std::logic_error("ERROR:EDeN::Graph::GetNeighbors: Something went wrong: no values are present in cache for vertex id=" + stream_cast<string>(aSrcVertexID) + " and distance=" + stream_cast<string>(aDistance));
	oNeighborsList = mNeighborsCache[aSrcVertexID][aDistance];
}

void EDeN::Graph::GetNeighbors(unsigned aSrcVertexID, unsigned aDistance, vector<list<unsigned> >& oNeighborDistanceList) const {
	//if not cached then compute answer
	if (mNeighborsCache[aSrcVertexID][aDistance].size() == 0) {
		//compute distance from source to all vertices up to aDistance and fill cache
		SingleVertexBoundedBreadthFirstVisit(aSrcVertexID, aDistance);
	}
	oNeighborDistanceList[aDistance] = mNeighborsCache[aSrcVertexID][aDistance];
}

void EDeN::Graph::SingleVertexBoundedBreadthFirstVisit(unsigned aRootVertexIndex, unsigned aDistance) const {
	//breadth first visit over allowed vertices
	mDestMaptoDistance[aRootVertexIndex] = 0;
	map<unsigned, bool> already_explored; //NOTE: we use a map instead of a vector since the limited depth breadth first visit can visit fewer vertices than there are vertices in total
	already_explored[aRootVertexIndex] = true;
	queue<unsigned> q;
	q.push(aRootVertexIndex); //initialize queue with the root vertex
	while (q.empty() == false) {
		int u = q.front();
		for (list<unsigned>::const_iterator it = mAdjacencyList[u].begin(); it != mAdjacencyList[u].end(); it++) {
			int v = *it;
			if (already_explored[v] == true || IsTraversable(v) == false) {
				//do nothing, ignore the vertex
			} else {
				if (mDestMaptoDistance[u] + 1 <= aDistance) {
					unsigned dist = mDestMaptoDistance[u] + 1;
					mDestMaptoDistance[v] = dist;
					already_explored[v] = true;
					q.push(v);
					//update cache
					mNeighborsCache[aRootVertexIndex][dist].push_back(v);
				}
			}
		}
		q.pop();
	}
}

bool EDeN::Graph::IsTraversable(unsigned aVertexID) const {
	//TODO: find more general solution
	bool value = true;
	string att = "is_nesting_rel";
	if (mNodes[aVertexID].mBool.count(att) > 0) {
		value = mNodes[aVertexID].mBool.find(att)->second;
	} else
		value = true; //default
	return value;
}

double EDeN::Graph::GetWeight(unsigned aVertexID, const string aKey) const {
	double value = 1; //default
	if (mNodes[aVertexID].mNumber.count(aKey) > 0)
		value = mNodes[aVertexID].mNumber.find(aKey)->second;
	return value;
}

string EDeN::Graph::GetLabel(unsigned aVertexID, const string aKey) const {
	string value = "."; //default
	if (mNodes[aVertexID].mDiscrete.count(aKey) > 0)
		value = mNodes[aVertexID].mDiscrete.find(aKey)->second;
	return value;
}

bool EDeN::Graph::HasSVector(unsigned aVertexID, const string aKey) const {
	bool value = false; //default
	if (mNodes[aVertexID].mSVector.count(aKey) > 0)
		value = true;
	return value;
}

bool EDeN::Graph::HasDVector(unsigned aVertexID, const string aKey) const {
	bool value = false; //default
	if (mNodes[aVertexID].mDVector.count(aKey) > 0)
		value = true;
	return value;
}

unsigned EDeN::Graph::SizeSVector(unsigned aVertexID, const string aKey) const {
	unsigned value = 0; //default
	if (mNodes[aVertexID].mSVector.count(aKey) > 0) {
		const SVector& svec = mNodes[aVertexID].mSVector.find(aKey)->second;
		value = svec.size();
	}
	return value;
}

unsigned EDeN::Graph::SizeDVector(unsigned aVertexID, const string aKey) const {
	unsigned value = 0; //default
	if (mNodes[aVertexID].mDVector.count(aKey) > 0) {
		const FVector& dvec = mNodes[aVertexID].mDVector.find(aKey)->second;
		value = dvec.size();
	}
	return value;
}

SVector EDeN::Graph::GetSVector(unsigned aVertexID, const string aKey) const {
	if (mNodes[aVertexID].mSVector.count(aKey) > 0) {
		unsigned size = SizeSVector(aVertexID, aKey);
		SVector svec(size);
		svec = mNodes[aVertexID].mSVector.find(aKey)->second;
		return svec;
	}
	return SVector();
}

FVector EDeN::Graph::GetDVector(unsigned aVertexID, const string aKey) const {
	if (mNodes[aVertexID].mDVector.count(aKey) > 0) {
		unsigned size = SizeDVector(aVertexID, aKey);
		FVector dvec(size);
		dvec = mNodes[aVertexID].mDVector.find(aKey)->second;
		return dvec;
	}
	return FVector();
}

void EDeN::Graph::Parse(istream& aIn, const string& aFormat) {
	//read one graph at a time
	//read one line at a time
	//check if it is a vertex or an edge
	//get the reference id
	//manage the comments
	if (mParser.ParseRoot(aIn, mRoot)) {
		while (mParser.IsValid(aIn)) {
			if (mParser.TestForChar(aIn, 'v')) {
				Node node;
				if (mParser.ParseVertex(aIn, node))
					mNodes.push_back(node);
			}
			if (mParser.TestForChar(aIn, 'e')) {
				Node node;
				if (mParser.ParseEdge(aIn, node)) {
					//UpdateAdjacencyList(node);
					//mNodes.push_back(node);
				}
			}
		}
	} else
		throw logic_error("ERROR: the format is not recognized. \n in code file: " + stream_cast<string>(__FILE__) + " line: " + stream_cast<string>(__LINE__) + " function: " + stream_cast<string>(__FUNCTION__));
}

//---------------------------------------------------------------------------------------------------------------------------------------------
istream& EDeN::Parser::TestForChar(istream& in, char aChar) {
	if ((in >> std::ws).peek() != std::char_traits<char>::to_int_type(aChar)) {
		in.setstate(std::ios_base::failbit);
	}
	return in;
}

bool EDeN::Parser::IsValid(istream& aIn) {
	//skip any number of comment lines, i.e. lines that start with a '/' symbol
	while (TestForChar(aIn, '/')) {
		string line;
		getline(aIn, line);
	}
	//check for vertices or edges
	if (TestForChar(aIn, 'v') || TestForChar(aIn, 'e'))
		return true;
	else
		return false;
}

bool EDeN::Parser::ParseRoot(istream& aIn, Node& oRoot) {
	//a valid graph is introduced by a '#' token
	if (TestForChar(aIn, '#')) {
		string symbol;
		aIn >> symbol;
		if (symbol != "#")
			throw logic_error("ERROR: expecting the token '#' instead of '" + symbol + "' \n in file: " + stream_cast<string>(__FILE__) + " line: " + stream_cast<string>(__LINE__) + " function: " + stream_cast<string>(__FUNCTION__));
		//the first token after the # symbol is the graph identifier
		string graph_identifier;
		aIn >> graph_identifier;
		oRoot.mDiscrete.insert(std::pair<string, string>(GRAPH_IDENTIFIER_KEY, graph_identifier));
		//the rest of tokens are attributes
		while (ParseAttributes(aIn, oRoot)) {
		}
	} else
		throw logic_error("ERROR: the '#' symbol is missing \n in code file: " + stream_cast<string>(__FILE__) + " line: " + stream_cast<string>(__LINE__) + " function: " + stream_cast<string>(__FUNCTION__));
	return true;
}

bool EDeN::Parser::ParseVertex(istream& aIn, Node& oNode) {
	//a valid vertex is introduced by a 'v' token
	if (TestForChar(aIn, 'v')) {
		string symbol;
		aIn >> symbol;
		if (symbol != "v")
			throw logic_error("ERROR: expecting the token 'v' instead of '" + symbol + "' \n in file: " + stream_cast<string>(__FILE__) + " line: " + stream_cast<string>(__LINE__) + " function: " + stream_cast<string>(__FUNCTION__));
		//mark the node as a vertex
		oNode.mDiscrete.insert(std::pair<string, string>(NODE_TYPE_KEY, VERTEX_TYPE_VALUE));
		//the first token after the v symbol is the vertex identifier
		string vertex_identifier;
		aIn >> vertex_identifier;
		oNode.mDiscrete.insert(std::pair<string, string>(VERTEX_IDENTIFIER_KEY, vertex_identifier));
		//the rest of tokens are attributes
		while (ParseAttributes(aIn, oNode)) {
		}
	} else
		throw logic_error("ERROR: the 'v' symbol is missing \n in code file: " + stream_cast<string>(__FILE__) + " line: " + stream_cast<string>(__LINE__) + " function: " + stream_cast<string>(__FUNCTION__));
	return true;
}

bool EDeN::Parser::ParseEdge(istream& aIn, Node& oNode) {
	//a valid edge is introduced by a 'e' token
	if (TestForChar(aIn, 'e')) {
		string symbol;
		aIn >> symbol;
		if (symbol != "e")
			throw logic_error("ERROR: expecting the token 'e' instead of '" + symbol + "' \n in file: " + stream_cast<string>(__FILE__) + " line: " + stream_cast<string>(__LINE__) + " function: " + stream_cast<string>(__FUNCTION__));
		//mark the node as an edge
		oNode.mDiscrete.insert(std::pair<string, string>(NODE_TYPE_KEY, EDGE_TYPE_VALUE));
		//the first token after the e symbol is the source vertex identifier, the second is the destination vertex identifier
		string edge_src_identifier;
		string edge_dest_identifier;
		aIn >> edge_src_identifier >> edge_dest_identifier;
		oNode.mDiscrete.insert(std::pair<string, string>(EDGE_SRC_IDENTIFIER_KEY, edge_src_identifier));
		oNode.mDiscrete.insert(std::pair<string, string>(EDGE_DEST_IDENTIFIER_KEY, edge_dest_identifier));
		//the rest of tokens are attributes
		while (ParseAttributes(aIn, oNode)) {
		}
	} else
		throw logic_error("ERROR: the 'e' symbol is missing \n in code file: " + stream_cast<string>(__FILE__) + " line: " + stream_cast<string>(__LINE__) + " function: " + stream_cast<string>(__FUNCTION__));
	return true;
}

bool EDeN::Parser::ParseAttributes(istream& aIn, Node& oNode) {
	while(ParseAttributeValuePair(aIn,oNode)){}
	return true;
}

bool EDeN::Parser::ParseAttributeValuePair(istream& aIn, Node& oNode) {
	//extract the key
	string key;
	aIn >> key;
	//extract the '=' symbol
	string key_value_separator;
	aIn >> key_value_separator;
	if (key_value_separator != "=")
		throw logic_error("ERROR: expecting the token '=' instead of '" + key_value_separator + "' \n in file: " + stream_cast<string>(__FILE__) + " line: " + stream_cast<string>(__LINE__) + " function: " + stream_cast<string>(__FUNCTION__));
	//determine the data type for the value that follows
	string next_token;
	aIn >> next_token;
	//if the token is "true" or "false" then it is boolean
	if (next_token == "true")
		oNode.mBool.insert(std::pair<string, bool>(key, true));
	if (next_token == "false")
		oNode.mBool.insert(std::pair<string, bool>(key, false));
	if (isdigit(next_token.at(0)) == false) //if the first char is not a digit then it is discrete (i.e. a string)
		oNode.mDiscrete.insert(std::pair<string, string>(key, next_token));
	//else if the token contains the symbol ':' it is a sparse vector
	std::size_t found = next_token.find(":");
	if (found != std::string::npos) {
		//determine the number of tokens
		//initialize a sparse vector
	}
	if (isdigit(next_token.at(0)) == true) {
		//else if the next token starts with a digit then it is a dense vector
		//else it is a number
	}

	return false;//silently bail out if the format is not recognised
}

//---------------------------------------------------------------------------------------------------------------------------------------------
inline
unsigned EDeN::Kernel::HashFunc(const vector<unsigned>& aList, unsigned aBitMask) const {
	unsigned int hash = 0xAAAAAAAA;
	for (std::size_t i = 0; i < aList.size(); i++) {
		hash ^= ((i & 1) == 0) ? ((hash << 7) ^ aList[i] * (hash >> 3)) : (~(((hash << 11) + aList[i]) ^ (hash >> 5)));
	}
	return hash & aBitMask;
}

inline vector<unsigned> EDeN::Kernel::HashFuncVec(const vector<unsigned>& aList, unsigned aBitMask) const {
	vector<unsigned> hash_vec(aList.size());
	unsigned int hash = 0xAAAAAAAA;
	for (std::size_t i = 0; i < aList.size(); i++) {
		hash ^= ((i & 1) == 0) ? ((hash << 7) ^ aList[i] * (hash >> 3)) : (~(((hash << 11) + aList[i]) ^ (hash >> 5)));
		hash_vec[i] = hash & aBitMask;
	}
	return hash_vec;
}

inline
unsigned EDeN::Kernel::HashFuncString(const string& aString, unsigned aBitMask) const { //NOTE: extract the least significant bits from the hash
	unsigned int hash = 0xAAAAAAAA;
	for (std::size_t i = 0; i < aString.length(); i++) {
		hash ^= ((i & 1) == 0) ? ((hash << 7) ^ aString[i] * (hash >> 3)) : (~(((hash << 11) + aString[i]) ^ (hash >> 5)));
	}
	return hash & aBitMask;
}

double EDeN::Kernel::ComputeKernel(const SVector& aX, const SVector& aZ) {
	double k = 0;
	k = aX.dot(aZ);
	return k;
}

double EDeN::Kernel::ComputeKernel(const Graph& aG, const Graph& aM) {
	SVector xg(pow(2, mHashBitSize));
	SVector xm(pow(2, mHashBitSize));
	GenerateFeatureVector(aG, xg);
	GenerateFeatureVector(aM, xm);
	return ComputeKernel(xg, xm);
}

void EDeN::Kernel::GenerateFeatureVector(const Graph& aG, SVector& oX) {
//TODO insert a type: graph or string; for string use a specialized code

//extract vertices of interest
	list<unsigned> vertex_list;
	SelectVertexListPolicy(aG, vertex_list);

//for each vertex
	for (list<unsigned>::iterator it = vertex_list.begin(); it != vertex_list.end(); ++it) {
		unsigned id = *it;
		//for each attribute
		ComputeFeatures(aG, id);
	}

//normalization policy
	NormalizationPolicy(oX);
}

void EDeN::Kernel::SelectVertexListPolicy(const Graph& aG, list<unsigned>& oVertexList) const {
//check attribute weight=0 to exclude vertices
	for (unsigned i = 0; i < aG.Size(); i++) {
		if (aG.GetWeight(i) > 0)
			oVertexList.push_back(i);
	}
}

void EDeN::Kernel::NormalizationPolicy(SVector& oX) const {
//extract all different keys
	set<string> keys;
	for (multimap<string, SVector>::const_iterator it = mAttributeFeatureVectorDictionary.begin(); it != mAttributeFeatureVectorDictionary.end(); ++it)
		keys.insert(it->first);

//sum all feature vectors from same attribute (e.g. same radius and distance) and normalize
	for (set<string>::const_iterator it = keys.begin(); it != keys.end(); ++it) {
		string key = (*it);
		SVector x(pow(2, mHashBitSize));
		AggregateFeatures(key, x);
		x /= x.norm();
		oX += x;
	}
//normalize the final vector
	oX /= oX.norm();
}

void EDeN::Kernel::AggregateFeatures(string aKey, SVector& oX) const {
//compute aggregate feature and cache it
	string aggregate_key = EDeN::AGGREGATE_KEY + aKey;
	SVector x(pow(2, mHashBitSize));
	if (mAttributeFeatureVectorDictionary.count(aggregate_key) == 0) {
		pair<multimap<string, SVector>::const_iterator, multimap<string, SVector>::const_iterator> ret;
		ret = mAttributeFeatureVectorDictionary.equal_range(aKey);
		for (multimap<string, SVector>::const_iterator it = ret.first; it != ret.second; ++it) {
			x += it->second;
		}
		mAttributeFeatureVectorDictionary.insert(std::pair<string, SVector>(aggregate_key, x));
	}
	oX = mAttributeFeatureVectorDictionary.find(aggregate_key)->second;
}

void EDeN::Kernel::ComputeFeatures(const Graph& aG, unsigned aID) const {
//TODO: manage string kernel

//TODO:for each type of attribute present in the vertex attribute list
//for now: just all possible radius and distance values
	ComputeFeaturesRadiusDistance(aG, aID);
	ComputeFeaturesNonRootRadiusDistance(aG, aID);
	ComputeFeaturesDenseRealVector(aG, aID);

//TODO: add the nested features with abstract relations
}

string EDeN::Kernel::GenerateKeyRadiusDistance(unsigned aRadius, unsigned aDistance) const {
	string key = "r=" + stream_cast<string>(aRadius) + "," + "d=" + stream_cast<string>(aDistance);
	return key;
}

void EDeN::Kernel::ComputeFeaturesRadiusDistance(const Graph& aG, unsigned aID) const {
	for (unsigned r = 0; r <= mRadius; r++) {
		for (unsigned d = 0; d <= mDistance; d++) {
			string key = GenerateKeyRadiusDistance(r, d);
			SVector x(pow(2, mHashBitSize));
			ComputeFeaturesRadiusDistance(aG, aID, r, d, x);
			mAttributeFeatureVectorDictionary.insert(std::pair<string, SVector>(key, x));
		}
	}
}

void EDeN::Kernel::ComputeFeaturesRadiusDistance(const Graph& aG, unsigned aID, unsigned aRadius, unsigned aDistance, SVector& oX) const {
//extract all vertices at distance, compose hash code and build sparse vector
	const unsigned feature_namespace = 1;

//get code for aID
	unsigned id_code = ComputeNeighborhoodGraphHashCode(aG, aID, aRadius);

//create a vector to store dist, radius, code src, code dest to generate final feature hash code
	vector<unsigned> code_vec(5);
	code_vec[4] = feature_namespace;
	code_vec[3] = aDistance;
	code_vec[2] = aRadius;
	code_vec[1] = id_code;

//extract all distant vertices
	list<unsigned> neighbors_list;
	aG.GetNeighbors(aID, aDistance, neighbors_list);
	for (list<unsigned>::const_iterator it = neighbors_list.begin(); it != neighbors_list.end(); ++it) {
		unsigned v = *it;
		//get code for distant vertices
		unsigned dist_code = ComputeNeighborhoodGraphHashCode(aG, v, aRadius);
		code_vec[0] = dist_code;
		//compute feature code: dist, radius, code src, code dest
		unsigned hash_subgraph_code = HashFunc(code_vec, pow(2, mHashBitSize) - 1);
		oX.coeffRef(hash_subgraph_code) += 1;
	}
//multiply counts by the weight if given, otherwise resort to default=1
	double w = aG.GetWeight(aID);
	if (w != 1)
		oX *= w;
}

void EDeN::Kernel::ComputeFeaturesNonRootRadiusDistance(const Graph& aG, unsigned aID) const {
	for (unsigned r = 0; r <= mRadius; r++) {
		for (unsigned d = 0; d <= mDistance; d++) {
			string key = "non_root,r=" + stream_cast<string>(r) + "," + "d=" + stream_cast<string>(d);
			SVector x(pow(2, mHashBitSize));
			ComputeFeaturesNonRootRadiusDistance(aG, aID, r, d, x);
			mAttributeFeatureVectorDictionary.insert(std::pair<string, SVector>(key, x));
		}
	}
}

void EDeN::Kernel::ComputeFeaturesNonRootRadiusDistance(const Graph& aG, unsigned aID, unsigned aRadius, unsigned aDistance, SVector& oX) const {
//extract all vertices at distance, compose hash code and build sparse vector
	const unsigned feature_namespace = 2;

//assign fix code for root node
	const unsigned id_code = 1;

//create a vector to store dist, radius, code src, code dest to generate final feature hash code
	vector<unsigned> code_vec(5);
	code_vec[4] = feature_namespace;
	code_vec[3] = aDistance;
	code_vec[2] = aRadius;
	code_vec[1] = id_code;

//extract all distant vertices
	list<unsigned> neighbors_list;
	aG.GetNeighbors(aID, aDistance, neighbors_list);
	for (list<unsigned>::const_iterator it = neighbors_list.begin(); it != neighbors_list.end(); ++it) {
		unsigned v = *it;
		//get code for distant vertices
		unsigned dist_code = ComputeNeighborhoodGraphHashCode(aG, v, aRadius);
		code_vec[0] = dist_code;
		//compute feature code: dist, radius, code src, code dest
		unsigned hash_subgraph_code = HashFunc(code_vec, pow(2, mHashBitSize) - 1);
		oX.coeffRef(hash_subgraph_code) += 1;
	}
//multiply counts by the weight if given, otherwise resort to default=1
	double w = aG.GetWeight(aID);
	if (w != 1)
		oX *= w;
}

void EDeN::Kernel::ComputeFeaturesDenseRealVector(const Graph& aG, unsigned aID) const {
	for (unsigned r = 0; r <= mRadius; r++) {
		for (unsigned d = 0; d <= mDistance; d++) {
			string key = "dvec,r=" + stream_cast<string>(r) + "," + "d=" + stream_cast<string>(d);
			SVector x(pow(2, mHashBitSize));
			ComputeFeaturesDenseRealVector(aG, aID, r, d, x);
			if (x.nonZeros() > 0)
				mAttributeFeatureVectorDictionary.insert(std::pair<string, SVector>(key, x));
		}
	}
}

void EDeN::Kernel::ComputeFeaturesDenseRealVector(const Graph& aG, unsigned aID, unsigned aRadius, unsigned aDistance, SVector& oX) const {
//bail out if no vector information is associated to the node
	if (!aG.HasDVector(aID))
		return;
//retrieve features for plain case
	string key = GenerateKeyRadiusDistance(aRadius, aDistance);
	SVector disc_x(pow(2, mHashBitSize));
	AggregateFeatures(key, disc_x);

//extract dense vector for vertex aID
	unsigned dim = aG.SizeDVector(aID);
	FVector dvec(dim);
	dvec = aG.GetDVector(aID);

//for a number of times repeat
//if there exist some information in the real vector
//for each feature in x induce a random projection of the real_vector
//the resulting value is used as the value for the corresponding feature in a new rehashed feature
//this is repeated k=mNumRandProjections times
//the final num of features is then k+1 times the original discrete case

	const unsigned feature_namespace = 3;
//create a vector to store dist, radius, code src, code dest to generate final feature hash code
	vector<unsigned> code_vec(5);
	code_vec[4] = feature_namespace;
	code_vec[3] = aDistance;
	code_vec[2] = aRadius;

	for (unsigned k = 0; k < mNumProjections; ++k) {
		code_vec[1] = k;
		for (SVector::InnerIterator it(disc_x); it; ++it) {
			unsigned original_feature_id = it.index();
			//double original_feature_val = it.value();
			code_vec[0] = original_feature_id;
			unsigned code = HashFunc(code_vec, pow(2, mHashBitSize) - 1);
			double random_projection_value = MinHashProjection(dvec, code);
			if (random_projection_value)
				oX.coeffRef(code) += random_projection_value;
		}
	}
//multiply counts by the weight if given, otherwise resort to default=1
	double w = aG.GetWeight(aID);
	if (w != 1)
		oX *= w;
}

void EDeN::Kernel::ComputeFeaturesSparseRealVector(const Graph& aG, unsigned aID) const {
	for (unsigned r = 0; r <= mRadius; r++) {
		for (unsigned d = 0; d <= mDistance; d++) {
			string key = "dvec,r=" + stream_cast<string>(r) + "," + "d=" + stream_cast<string>(d);
			SVector x(pow(2, mHashBitSize));
			ComputeFeaturesSparseRealVector(aG, aID, r, d, x);
			if (x.nonZeros() > 0)
				mAttributeFeatureVectorDictionary.insert(std::pair<string, SVector>(key, x));
		}
	}
}

void EDeN::Kernel::ComputeFeaturesSparseRealVector(const Graph& aG, unsigned aID, unsigned aRadius, unsigned aDistance, SVector& oX) const {
//bail out if no vector information is associated to the node
	if (!aG.HasSVector(aID))
		return;
//retrieve features for plain case
	string key = GenerateKeyRadiusDistance(aRadius, aDistance);
	SVector disc_x(pow(2, mHashBitSize));
	AggregateFeatures(key, disc_x);

//extract dense vector for vertex aID
	unsigned dim = aG.SizeSVector(aID);
	SVector svec(dim);
	svec = aG.GetSVector(aID);

//for a number of times repeat
//if there exist some information in the real vector
//for each feature in x induce a random projection of the real_vector
//the resulting value is used as the value for the corresponding feature in a new rehashed feature
//this is repeated k=mNumRandProjections times
//the final num of features is then k+1 times the original discrete case

	const unsigned feature_namespace = 4;
//create a vector to store dist, radius, code src, code dest to generate final feature hash code
	vector<unsigned> code_vec(5);
	code_vec[4] = feature_namespace;
	code_vec[3] = aDistance;
	code_vec[2] = aRadius;

	for (unsigned k = 0; k < mNumProjections; ++k) {
		code_vec[1] = k;
		for (SVector::InnerIterator it(disc_x); it; ++it) {
			unsigned original_feature_id = it.index();
			//double original_feature_val = it.value();
			code_vec[0] = original_feature_id;
			unsigned code = HashFunc(code_vec, pow(2, mHashBitSize) - 1);
			double random_projection_value = MinHashProjection(svec, code);
			if (random_projection_value)
				oX.coeffRef(code) += random_projection_value;
		}
	}
//multiply counts by the weight if given, otherwise resort to default=1
	double w = aG.GetWeight(aID);
	if (w != 1)
		oX *= w;
}

unsigned EDeN::Kernel::ComputeNeighborhoodGraphHashCode(const Graph& aG, unsigned aID, unsigned aRadius) const {
//if not cached then compute answer and cache it
	if (mIdRadiusHashCodeCache[aID][aRadius] == 0) {
		CacheNeighborhoodGraphHashCode(aG, aID, aRadius);
	}
//always return the cached answer
	return mIdRadiusHashCodeCache[aID][aRadius];
}

void EDeN::Kernel::CacheNeighborhoodGraphHashCode(const Graph& aG, unsigned aID, unsigned aMaxRadius) const {
//compute the incremental encoding of all neigh. vertices for increasing radii and fill cache

//get all ids at increasing distances
	vector<list<unsigned> > neighbor_distance_list;
	aG.GetNeighbors(aID, aMaxRadius, neighbor_distance_list);
	unsigned max_radius = neighbor_distance_list.size();
	vector<unsigned> code_vec(max_radius);
	for (unsigned r = 0; r <= max_radius; r++) {
		//get all vertices at distance r
		list<unsigned>& neighbors = neighbor_distance_list[r];
		vector<unsigned> code_list(neighbors.size());
		unsigned i = 0;
		//compute hash code for each label
		for (list<unsigned>::const_iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
			unsigned dist_id = *it;
			string label = aG.GetLabel(dist_id);
			code_list[i] = HashFuncString(label);
		}
		//sort codes to obtain canonical representation within fixed distance
		sort(code_list.begin(), code_list.end());
		//hash the list of codes into a single code representing all vertices at distance r
		unsigned code = HashFunc(code_list);
		code_vec[r] = code;
	}
//extract incrementally a code for distances 0, 0,1, then 0,1,2 etc
	vector<unsigned> hash_code_vec = HashFuncVec(code_vec);

//cache the resulting codes as these are exacly the codes for neighborhoods with the desired radius
	for (unsigned i = 0; i < hash_code_vec.size(); ++i)
		mIdRadiusHashCodeCache[aID][i] = hash_code_vec[i];
}

inline double EDeN::Kernel::PositiveOrthantRandomProjection(const FVector& aRealVector, unsigned aCode) const {
	double projection = 0;
	srand(aCode); //Initialize the random number generator with aCode so to be deterministic
	for (unsigned i = 0; i < aRealVector.size(); i++) {
		if (random01() < 0.5)
			projection += aRealVector(i);
	}
	return projection;
}

inline double EDeN::Kernel::RandomProjection(const FVector& aRealVector, unsigned aCode) const {
//refrerence: Achlioptas, Dimitris. "Database-friendly random projections." Proceedings of the twentieth ACM SIGMOD-SIGACT-SIGART symposium on Principles of database systems. ACM, 2001.
	const double THRESHOLD = .3333333;
	const double SQRT3 = 1.73205080756887729352;
	double projection = 0;
	srand(aCode); //Initialize the random number generator with aCode so to be deterministic
	for (unsigned i = 0; i < aRealVector.size(); i++) {
		if (random01() <= THRESHOLD) {
			if (random01() < 0.5)
				projection += SQRT3 * aRealVector(i);
			else
				projection -= SQRT3 * aRealVector(i);
		}
	}
	return projection;
}

inline double EDeN::Kernel::PositiveOrthantRandomProjection(const SVector& aRealVector, unsigned aCode) const {
	double projection = 0;
	vector<unsigned> code_vec(2);
	code_vec[1] = aCode;
	for (SVector::InnerIterator it(aRealVector); it; ++it) {
		unsigned original_feature_id = it.index();
		double original_feature_value = it.value();
		code_vec[0] = original_feature_id;
		unsigned code = HashFunc(code_vec, pow(2, mHashBitSize) - 1);
		srand(code); //Initialize the random number generator with aCode so to be deterministic
		if (random01() < 0.5)
			projection += original_feature_value;
	}
	return projection;
}

inline double EDeN::Kernel::RandomProjection(const SVector& aRealVector, unsigned aCode) const {
//refrerence: Achlioptas, Dimitris. "Database-friendly random projections." Proceedings of the twentieth ACM SIGMOD-SIGACT-SIGART symposium on Principles of database systems. ACM, 2001.
	const double THRESHOLD = .3333333;
	const double SQRT3 = 1.73205080756887729352;
	double projection = 0;
	vector<unsigned> code_vec(2);
	code_vec[1] = aCode;
	for (SVector::InnerIterator it(aRealVector); it; ++it) {
		unsigned original_feature_id = it.index();
		double original_feature_value = it.value();
		code_vec[0] = original_feature_id;
		unsigned code = HashFunc(code_vec, pow(2, mHashBitSize) - 1);
		srand(code); //Initialize the random number generator with aCode so to be deterministic
		if (random01() <= THRESHOLD) {
			if (random01() < 0.5)
				projection += SQRT3 * original_feature_value;
			else
				projection -= SQRT3 * original_feature_value;
		}
	}
	return projection;
}

//TODO: unify the code, this is a generic function, it should just have an iterator over the container, be it dense or sparse vector

inline double EDeN::Kernel::MinHashProjection(const SVector& aRealVector, unsigned aCode) const {
	vector<unsigned> code_vec(2);
//initialize with first value
	code_vec[1] = aCode;
	SVector::InnerIterator it(aRealVector);
	unsigned original_feature_id = it.index();
	double original_feature_value = it.value();
	code_vec[0] = original_feature_id;
	unsigned code = HashFunc(code_vec, pow(2, mHashBitSize) - 1);
	unsigned min_code = code;
	double min_value = original_feature_value;
//for all feature values
	for (SVector::InnerIterator it(aRealVector); it; ++it) {
		unsigned original_feature_id = it.index();
		double original_feature_value = it.value();
		code_vec[0] = original_feature_id;
		unsigned code = HashFunc(code_vec, pow(2, mHashBitSize) - 1);
		if (code < min_code) {
			min_code = code;
			min_value = original_feature_value;
		}
	}
	return min_value;
}

inline double EDeN::Kernel::MinHashProjection(const FVector& aRealVector, unsigned aCode) const {
	const unsigned FEATURE_ID_OFFSET = 1;
	vector<unsigned> code_vec(2);
//initialize with first value
	code_vec[1] = aCode;
	code_vec[0] = FEATURE_ID_OFFSET;
	unsigned code = HashFunc(code_vec, pow(2, mHashBitSize) - 1);
	double original_feature_value = aRealVector(0);
	unsigned min_code = code;
	double min_value = original_feature_value;
//for all feature values
	for (unsigned i = 1; i < aRealVector.size(); i++) {
		unsigned original_feature_id = i + FEATURE_ID_OFFSET;
		double original_feature_value = aRealVector(i);
		code_vec[0] = original_feature_id;
		unsigned code = HashFunc(code_vec, pow(2, mHashBitSize) - 1);
		if (code < min_code) {
			min_code = code;
			min_value = original_feature_value;
		}
	}
	return min_value;
}

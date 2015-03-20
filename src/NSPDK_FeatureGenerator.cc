#include "NSPDK_FeatureGenerator.h"

//----------------------------------------------------------------------------------------------------------------------------
void DebugClass::Clear() {
	mHashToFeatureMap.clear();
	mHashToGraphMap.clear();
}

void DebugClass::Output(ostream& out) const {
	OutputFeatureEncoding(out);
}

void DebugClass::OutputFeatureEncoding(ostream& out) const {
	out << "#Feature encodings: [" << mHashToFeatureMap.size() << "]" << endl;
	for (map<unsigned, string>::const_iterator it = mHashToFeatureMap.begin(); it != mHashToFeatureMap.end(); ++it)
		out << it->first << " -> " << it->second << endl;
}

void DebugClass::StoreFeatureCodeToFeatureInfo(unsigned aFeatureCode, vector<unsigned>& aDetailsList) {
	if (mHashToFeatureMap.count(aFeatureCode) == 0) {
		string feature_information = "r:" + stream_cast<string>(aDetailsList[0]) + " d:" + stream_cast<string>(aDetailsList[1]);
		feature_information += " [g: " + mHashToGraphMap[aDetailsList[2]] + "]";
		feature_information += " [g: " + mHashToGraphMap[aDetailsList[3]] + "]";
		mHashToFeatureMap.insert(make_pair(aFeatureCode, feature_information));
	}
}

void DebugClass::SerializedRootedGraphCanonicalFormEncoding(unsigned aDiscreteEncoding, int aRootVertexIndex, const GraphClass& aG, int aRadius) {
	string value = "";
	if (aRadius > 0) {
		//extract set of vertices in the ball of radius aMaxDepth
		set<unsigned> ball;
		for (int r = 0; r <= aRadius; r++) {
			vector<unsigned> dest_id_list = aG.GetFixedDistanceVertexIDList(aRootVertexIndex, r);
			ball.insert(dest_id_list.begin(), dest_id_list.end());
		}
		if (ball.size() == 0)
			throw std::logic_error("ERROR:DebugClass::SerializedRootedGraphCanonicalFormEncoding: Something went wrong in SerializedRootedGraphCanonicalFormEncoding: cannot generate features over an empty neighbourhood graph!");
		//induce the subgraph from the ball and return the new index for the root vertex
		GraphClass gal;
		unsigned root = aG.GetVertexInducedRootedSubGraph(ball, aRootVertexIndex, gal);
		string root_label = gal.GetVertexLabel(root);
		root_label += "*";
		gal.SetVertexLabel(root, root_label);
		value = gal.Serialize();
	} else
		value = "*" + aG.GetVertexLabelConcatenated(aRootVertexIndex);
	mHashToGraphMap.insert(make_pair(aDiscreteEncoding, value));
}

//----------------------------------------------------------------------------------------------------------------------------------------------
NSPDK_FeatureGenerator::NSPDK_FeatureGenerator(const std::string& id) :
		FeatureGenerator(id)  {
	mRadius = 0;
	mDistance = 0;
	mMatchType = "hard";
	mHashBitSize = (unsigned) (numeric_limits<unsigned>::digits - 1);
	mHashBitMask = numeric_limits<unsigned>::max() >> 1;
	mMinKernel = false;
	mNormalization = true;
	mUseRealVectorInformation = false;
	mNumRandProjections = 1;
	mDebugVerbosity = 0;
	mVertexDegreeThreshold = 10;
//	new_flag(&mVertexDegreeThreshold, "vertex_degree_threshold", "(unsigned)\nThreshold vertex degree above which features for a specific vertex are not generated");
//	new_flag(&mRadius, "radius", "(unsigned)\nMax radius of kernel neighborhoods");
//	new_flag(&mDistance, "distance", "(unsigned)\nMax distance between pairs of neighborhoods");
//	new_flag(&mMatchType, "match_type", "(string)\nHow to match neighborhoods: soft, hard");
//	new_flag(&mNormalization, "normalization", "(bool)\nNormalize feature vectors");
//	new_flag(&mMinKernel, "min_kernel", "(bool)\nApply the min-kernel on top of the generated features");
//	new_flag(&mNumRandProjections, "num_rand_projections", "(unsigned)\nNumber of hash funcions in the minhash signature");
//
//	new_flag(&mUseRealVectorInformation, "use_real_vector_information", "(bool)\nUse information encoded as a real vector per vertex");
//	new_flag(&mHashBitSize, "hash_bit_size", "(unsigned)\nNumber of bits for hash values"); //FIXME: since parameter variables are accessed directly the bit size is useless as setting it cannot trigger automatically the computation of the bit mask; the only solution is to set directly the bit mask itself
//	mHashBitMask = (2 << (mHashBitSize - 1)) - 1;
//	new_flag(&mHashBitMask, "hash_bit_mask", "(unsigned)\nMask for hash values");
//	new_flag(&mDebugVerbosity, "verbosity", "(unsigned)\nNumber for debug verbosity level");
}

void NSPDK_FeatureGenerator::OutputParameters(ostream& out) const {
	out << "Radius: " << mRadius << endl;
	out << "Distance: " << mDistance << endl;
	out << "Match_Type: " << mMatchType << endl;
	out << "Hash_Bit_Size: " << mHashBitSize << endl;
	out << "Hash_Bit_mask: " << mHashBitMask << endl;
	out << "Min_Kernel: " << mMinKernel << endl;
	out << "Normalization: " << mNormalization << endl;
	out << "Use_Real_Vector_Information: " << mUseRealVectorInformation << endl;
	out << "Num_rand_projections: " << mNumRandProjections << endl;
	out << "Vertex Degree Threshold: " << mVertexDegreeThreshold << endl;
	out << "Debug_Verbosity: " << mDebugVerbosity << endl;
}

void NSPDK_FeatureGenerator::OutputFeatureMap(ostream& out) const {
	mDebugInfo.OutputFeatureEncoding(out);
}

unsigned NSPDK_FeatureGenerator::GenerateGraphHashCode(const GraphClass& aG, const vector<unsigned>& aFirstEndpointList) {
	SVector x(pow(2, mHashBitSize));
	generate_feature_vector(aG, x, aFirstEndpointList);
	unsigned code = GenerateVectorHashCode(x);
	return code;
}

void NSPDK_FeatureGenerator::InitFeatureCache(const GraphClass& aG, unsigned aRadius) {
	mFeatureCache.clear();
	unsigned vertex_size = aG.VertexSize();
	for (unsigned i = 0; i <= aRadius; ++i) {
		mFeatureCache.push_back(vector<unsigned>(vertex_size, 0));
	}
}

void NSPDK_FeatureGenerator::generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList) {
	vector<unsigned> first_endpoint_list = aFirstEndpointList;
	GetFirstEndpoints(aG, first_endpoint_list);
	if (first_endpoint_list.size() == 0)
		throw std::logic_error("ERROR: NSPDK_FeatureGenerator::generate_feature_vector: Something went wrong: cannot generate features over an empty set of first endpoints!");
	int horizon = max(mDistance, mRadius);
	aG.ComputePairwiseDistanceInformation(horizon, first_endpoint_list);
	if (aG.Check() == false)
		throw logic_error("ERROR: NSPDK_FeatureGenerator::generate_feature_vector: the graph data structure is not sound and it has not passed the checkup procedure"); //check graph data structure soundness

	InitFeatureCache(aG, mRadius);
	for (unsigned r = 0; r <= mRadius; r++) {
		for (unsigned d = 0; d <= mDistance; d++) {
			SVector z(pow(2, mHashBitSize));
			for (unsigned i = 0; i < first_endpoint_list.size(); i++) {
				unsigned src_id = first_endpoint_list[i];
				if (aG.GetVertexViewPoint(src_id) && aG.GetVertexKernelPoint(src_id) && aG.GetVertexAlive(src_id)) { //proceed to extract features only if the *src* vertex is a kernel point and is alive
					SVector zv(pow(2, mHashBitSize));
					GenerateVertexFeatures(src_id, aG, r, d, zv);
					z += zv;
				}
			}
			if (mNormalization)
				z /= z.norm();
			x += z;
		}
	}

	if (mNormalization)
		x /= x.norm();
	if (mMinKernel)
		ConvertSparseVectorToMinFeatureVector(x);
	if (mDebugVerbosity > 0) {
		cout << x;
		cout << endl;
		OutputFeatureMap(cout);
		aG.Output(cout);
	}
}

void NSPDK_FeatureGenerator::generate_vertex_feature_vector(const GraphClass& aG, vector<SVector>& x_list, const vector<unsigned>& aFirstEndpointList) {
	vector<unsigned> first_endpoint_list = aFirstEndpointList;
	GetFirstEndpoints(aG, first_endpoint_list);
	if (first_endpoint_list.size() == 0)
		throw std::logic_error("ERROR:NSPDK_FeatureGenerator::generate_vertex_feature_vector: Something went wrong: cannot generate features over an empty set of first endpoints!");
	int horizon = max(mDistance, mRadius);
	aG.ComputePairwiseDistanceInformation(horizon, first_endpoint_list);
	if (aG.Check() == false)
		throw logic_error("ERROR:NSPDK_FeatureGenerator::generate_vertex_feature_vector: the graph data structure is not sound and it has not passed the checkup procedure"); //check graph data structure soundness

	InitFeatureCache(aG, mRadius);
	for (unsigned i = 0; i < first_endpoint_list.size(); i++) {
		SVector z(pow(2, mHashBitSize));
		unsigned src_id = first_endpoint_list[i];
		if (aG.GetVertexViewPoint(src_id) && aG.GetVertexKernelPoint(src_id) && aG.GetVertexAlive(src_id)) { //proceed to extract features only if the *src* vertex is a kernel point and is alive
			for (unsigned r = 0; r <= mRadius; r++) {
				for (unsigned d = 0; d <= mDistance; d++) {
					SVector zv(pow(2, mHashBitSize));
					GenerateVertexFeatures(src_id, aG, r, d, zv);
					z += zv;
				}
			}
		}
		if (mMinKernel)
			ConvertSparseVectorToMinFeatureVector(z);
		x_list.push_back(z);
	}
	if (mDebugVerbosity > 0) {
		for (unsigned i = 0; i < x_list.size(); ++i) {
			cout << i << " " << x_list[i];
			cout << endl;
		}
		OutputFeatureMap(cout);
		aG.Output(cout);
	}
}

unsigned NSPDK_FeatureGenerator::GenerateVectorHashCode(SVector& x) {
	vector<unsigned> hash_vec;
	for (SVector::InnerIterator it(x); it; ++it) {
		int key = it.index();
		double val = it.value();
		hash_vec.push_back((unsigned) key);
		hash_vec.push_back((unsigned) (val * 10000)); //Note: in general the value associated to a feature is a real number and it is therefore converted into an unsigned via scaling
	}
	unsigned code = HashFunc(hash_vec, mHashBitMask);
	return code;
}

unsigned NSPDK_FeatureGenerator::GenerateVertexHashCode(unsigned aSrcID, const GraphClass& aG) {
	SVector x(pow(2, mHashBitSize));
	GenerateVertexFeatures(aSrcID, aG, x);
	unsigned code = GenerateVectorHashCode(x);
	return code;
}

void NSPDK_FeatureGenerator::GenerateVertexFeatures(unsigned aSrcID, const GraphClass& aG, SVector& x) {
	for (unsigned r = 0; r <= mRadius; r++) {
		for (unsigned d = 0; d <= mDistance; d++) {
			SVector zv(pow(2, mHashBitSize));
			GenerateVertexFeatures(aSrcID, aG, r, d, zv);
			x += zv;
		}
	}
}

void NSPDK_FeatureGenerator::GenerateVertexFeatures(unsigned aSrcID, const GraphClass& aG, unsigned aRadius, unsigned aDistance, SVector& x) {
	vector<unsigned> endpoint_list(4);
	endpoint_list[0] = aRadius;
	endpoint_list[1] = aDistance;

	SVector xl(pow(2, mHashBitSize));
	unsigned src_code = GenerateVertexNeighbourhoodHashCode(aSrcID, aG, aRadius);
	vector<unsigned> dest_id_list = aG.GetFixedDistanceVertexIDList(aSrcID, aDistance);
	for (unsigned dest_j = 0; dest_j < dest_id_list.size(); dest_j++) {
		unsigned dest_id = dest_id_list[dest_j];
		unsigned dest_code = 0;
		if (aG.GetVertexKernelPoint(dest_id) && aG.GetVertexAlive(dest_id)) { //proceed to extract features only if the *dest* vertex is a kernel point and is alive
			dest_code = GenerateVertexNeighbourhoodHashCode(dest_id, aG, aRadius);
			//impose canonical order for pair: i.e. A-B and B-A must generate the same feature
			if (src_code < dest_code) {
				endpoint_list[2] = src_code;
				endpoint_list[3] = dest_code;
			} else {
				endpoint_list[2] = dest_code;
				endpoint_list[3] = src_code;
			}
			unsigned code = HashFunc(endpoint_list, mHashBitMask);
			if (mDebugVerbosity > 0)
				mDebugInfo.StoreFeatureCodeToFeatureInfo(code, endpoint_list);
			x.coeffRef(code) += 1;

			//add feature that ignore src endpoint
			//NOTE: this is important when the label of src vertex is considered noisy, in this way, it is only the context that will define the features
			endpoint_list[2] = 0;
			endpoint_list[3] = dest_code;
			unsigned nosrc_code = HashFunc(endpoint_list, mHashBitMask);
			if (mDebugVerbosity > 0)
				mDebugInfo.StoreFeatureCodeToFeatureInfo(code, endpoint_list);
			xl.coeffRef(nosrc_code) += 1;
		}
	}
	x += xl;
	if (mUseRealVectorInformation) {
		SVector z(pow(2, mHashBitSize));
		GenerateVertexRealVectorFeatures(aSrcID, aG, xl, z);
		x += z;
	}
}

void NSPDK_FeatureGenerator::GenerateVertexRealVectorFeatures(unsigned aSrcID, const GraphClass& aG, const SVector& x, SVector& z) {
	const unsigned FEATURE_SPACE_SCRAMBLER = 7919; //use a const to rehash the feature ids
	//get real vector per vertex
	vector<double> real_vector = aG.GetVertexNumericAttributeList(aSrcID);
	if (real_vector.size() != 0) { //if there exist some information in the real vector
		//for each feature in x induce a random projection of the real_vector
		//the resulting value is used as the value for the corresponding feature in a new rehashed feature
		//this is repeated k=mNumRandProjections times
		//the final num of features is then k+1 times the original discrete case
		vector<unsigned> convolution_feature_id(3);
		convolution_feature_id[0] = FEATURE_SPACE_SCRAMBLER;
		for (unsigned k = 0; k < mNumRandProjections; ++k) {
			convolution_feature_id[1] = k;
			for (SVector::InnerIterator it(x); it; ++it) {
				unsigned original_feature_id = it.index();
				//double original_feature_val = it.value();
				convolution_feature_id[2] = original_feature_id;
				unsigned code = HashFunc(convolution_feature_id, mHashBitMask);
				double random_projection_value = RandomProjection(real_vector, code);
				if (random_projection_value)
					z.coeffRef(code) += random_projection_value;
			}
		}
	}
}

//inline double NSPDK_FeatureGenerator::RandomProjection(const vector<double>& aRealVector, unsigned aCode) {
//	double projection = 0;
//	srand(aCode); //Initialize the random number generator with aCode so to be deterministic
//	for (unsigned i = 0; i < aRealVector.size(); ++i) {
//		if (random01() < 0.5)
//			projection += aRealVector[i];
//	}
//	return projection;
//}

inline double NSPDK_FeatureGenerator::RandomProjection(const vector<double>& aRealVector, unsigned aCode) {
	//refrerence: Achlioptas, Dimitris. "Database-friendly random projections." Proceedings of the twentieth ACM SIGMOD-SIGACT-SIGART symposium on Principles of database systems. ACM, 2001.
	const double THRESHOLD = .3333333;
	const double SQRT3 = 1.73205080756887729352;
	double projection = 0;
	srand(aCode); //Initialize the random number generator with aCode so to be deterministic
	for (unsigned i = 0; i < aRealVector.size(); ++i) {
		if (random01() <= THRESHOLD) {
			if (random01() < 0.5)
				projection += SQRT3 * aRealVector[i];
			else
				projection -= SQRT3 * aRealVector[i];
		}
	}
	return projection;
}

unsigned NSPDK_FeatureGenerator::GenerateVertexNeighbourhoodHashCode(unsigned aSrcID, const GraphClass& aG, unsigned aRadius) {
	unsigned src_code = 1;
	if (mFeatureCache[aRadius][aSrcID] == 0) {
		src_code = RootedGraphCanonicalFormEncoding(aSrcID, aG, aRadius);
		mFeatureCache[aRadius][aSrcID] = src_code;
	} else {
		src_code = mFeatureCache[aRadius][aSrcID];
	}
	return src_code;
}

void NSPDK_FeatureGenerator::ConvertSparseVectorToMinFeatureVector(SVector& x) {
	vector<unsigned> hash_vec(2, 0);
	SVector z(pow(2, mHashBitSize));
	for (SVector::InnerIterator it(x); it; ++it) {
		int key = it.index();
		double val = it.value();
		hash_vec[0] = (unsigned) (key);
		for (unsigned j = 0; j < val; ++j) {
			hash_vec[1] = j;
			unsigned code = HashFunc(hash_vec, mHashBitMask);
			z.coeffRef(code) = 1; //NOTE: unresolved issues in case of collisions
		}
	}
	x = z;
}

unsigned NSPDK_FeatureGenerator::RootedGraphCanonicalFormEncoding(int aRootVertexIndex, const GraphClass& aG, int aRadius) {
	//NOTE:for efficiency reasons case radius=0 and radius=1 are treated as special cases
	unsigned discrete_encoding = 1;
	unsigned root_degree = aG.VertexAdjacentListSize(aRootVertexIndex);
	if (root_degree > mVertexDegreeThreshold)
		return discrete_encoding;

	if (aRadius == 0) { //return root label
		discrete_encoding = Radius0RootedGraphCanonicalFormEncoding(aRootVertexIndex, aG);
	} else if (aRadius == 1) { //return the sorted sequence of root's children
		discrete_encoding = Radius1RootedGraphCanonicalFormEncoding(aRootVertexIndex, aG);
	} else { //general case
		discrete_encoding = RadiusKRootedGraphCanonicalFormEncoding(aRootVertexIndex, aG, aRadius);
	}
	if (mDebugVerbosity > 0)
		mDebugInfo.SerializedRootedGraphCanonicalFormEncoding(discrete_encoding, aRootVertexIndex, aG, aRadius);
	return discrete_encoding;
}

unsigned NSPDK_FeatureGenerator::Radius0RootedGraphCanonicalFormEncoding(int aRootVertexIndex, const GraphClass& aG) {
	string encoding = aG.GetVertexLabelConcatenated(aRootVertexIndex);
	unsigned hash_subgraph_code = HashFunc(encoding, mHashBitMask);
	return hash_subgraph_code;
}

unsigned NSPDK_FeatureGenerator::Radius1RootedGraphCanonicalFormEncoding(int aRootVertexIndex, const GraphClass& aG) {
	unsigned hash_subgraph_code = 1;
	string encoding;
	encoding = aG.GetVertexLabelConcatenated(aRootVertexIndex) + ":";
	vector<pair<string, unsigned> > vertex_label_id_list;
	vector<unsigned> vertex_adjacency_list = aG.GetVertexAdjacentList(aRootVertexIndex);
	vector<unsigned> edge_adjacency_list = aG.GetEdgeAdjacentList(aRootVertexIndex);
	for (unsigned i = 0; i < vertex_adjacency_list.size(); ++i) {
		unsigned child_vertex_id = vertex_adjacency_list[i];
		unsigned child_edge_id = edge_adjacency_list[i];
		string child_label = aG.GetVertexLabelConcatenated(child_vertex_id) + "-" + aG.GetEdgeLabelConcatenated(child_edge_id);
		vertex_label_id_list.push_back(make_pair(child_label, child_vertex_id));
	}
	sort(vertex_label_id_list.begin(), vertex_label_id_list.end());
	if (vertex_label_id_list.size() > 0)
		encoding += vertex_label_id_list[0].first;
	for (unsigned i = 1; i < vertex_label_id_list.size(); i++)
		encoding += "." + vertex_label_id_list[i].first;
	hash_subgraph_code = HashFunc(encoding, mHashBitMask);
	return hash_subgraph_code;
}

unsigned NSPDK_FeatureGenerator::RadiusKRootedGraphCanonicalFormEncoding(int aRootVertexIndex, const GraphClass& aG, int aRadius) {
	unsigned hash_subgraph_code = 1;
	//extract set of vertices in the ball of radius aRadius
	set<unsigned> ball;
	for (int r = 0; r <= aRadius; r++) {
		vector<unsigned> dest_id_list = aG.GetFixedDistanceVertexIDList(aRootVertexIndex, r);
		ball.insert(dest_id_list.begin(), dest_id_list.end());
	}

	//induce the subgraph from the ball and return the new index for the root vertex
	GraphClass gal;
	unsigned root = aG.GetVertexInducedRootedSubGraph(ball, aRootVertexIndex, gal);
	gal.ComputePairwiseDistanceInformation(gal.VertexSize());

	hash_subgraph_code = RootedGraphCanonicalFormEncoding(gal, root);
	return hash_subgraph_code;
}

void NSPDK_FeatureGenerator::CanonicalGraphVertexList(const GraphClass& aG, vector<unsigned>& oVertexList) {
	oVertexList.clear();
	GraphClass g(aG);
	g.ComputePairwiseDistanceInformation(g.VertexSize());
	for (unsigned i = 0; i < g.VertexSize(); ++i) {
		oVertexList.push_back(RootedGraphCanonicalFormEncoding(g, i));
	}
}

unsigned NSPDK_FeatureGenerator::MultiRootedGraphCanonicalFormEncoding(vector<int>& aRootVertexIndexList, const GraphClass& aG) {
	GraphClass g(aG); //make a copy of the graph
	g.ComputePairwiseDistanceInformation(g.VertexSize()); //compute complete breadth first visit

	//for all vertices extract the vertex's signature: sorted distance from all the other vertices + their vertex label + distance from labeled root vertices
	vector<unsigned> vertex_encoding_list;
	for (unsigned i = 0; i < g.VertexSize(); ++i) {
		vector<unsigned> vertex_encoding;
		vector<unsigned> vertex_distance_list;

		for (unsigned t = 0; t < aRootVertexIndexList.size(); ++t) {
			int dist = aG.PairwiseDistance(aRootVertexIndexList[t], i);
			string distance_label = stream_cast<string>(dist) + aG.GetVertexLabelConcatenated(aRootVertexIndexList[t]);
			unsigned hash_distance_label = HashFunc(distance_label, mHashBitMask);
			vertex_distance_list.push_back(hash_distance_label);
		}

		for (unsigned j = 0; j < g.VertexSize(); ++j) {
			int dist = g.PairwiseDistance(i, j);
			string distance_label = stream_cast<string>(dist) + g.GetVertexLabelConcatenated(j);
			unsigned hash_distance_label = HashFunc(distance_label, mHashBitMask);
			vertex_distance_list.push_back(hash_distance_label);
		}
		sort(vertex_distance_list.begin(), vertex_distance_list.end());
		unsigned hash_encoding = HashFunc(vertex_distance_list, mHashBitMask);
		vertex_encoding_list.push_back(hash_encoding);
	}
	//extract list of all edge's signatures in the induced graph: v-signature,u-signature,label(uv)
	vector<unsigned> edge_list;
	vector<unsigned> edge_encoding(3);
	for (unsigned u = 0; u < g.VertexSize(); ++u) {
		//get all edges of vertex u
		vector<unsigned> vertex_adjacency_list = g.GetVertexAdjacentList(u);
		vector<unsigned> edge_adjacency_list = g.GetEdgeAdjacentList(u);
		if (vertex_adjacency_list.size() == 0) { //NOTE: if a vertex is isolated use its encoding as edge encoding
			edge_list.push_back(vertex_encoding_list[u]);
		}
		for (unsigned j = 0; j < vertex_adjacency_list.size(); ++j) {
			unsigned v = vertex_adjacency_list[j];

			if (vertex_encoding_list[u] < vertex_encoding_list[v]) {
				edge_encoding[0] = vertex_encoding_list[u];
				edge_encoding[1] = vertex_encoding_list[v];
			} else {
				edge_encoding[0] = vertex_encoding_list[v];
				edge_encoding[1] = vertex_encoding_list[u];
			}
			unsigned e = edge_adjacency_list[j];
			string edge_label = g.GetEdgeLabelConcatenated(e);
			unsigned hash_edge_label = HashFunc(edge_label, mHashBitMask);
			edge_encoding[2] = hash_edge_label;
			unsigned hash_edge_encoding = HashFunc(edge_encoding, mHashBitMask);
			edge_list.push_back(hash_edge_encoding);
		}
	}
	//the graph encoding is the sorted list of edge encodings
	sort(edge_list.begin(), edge_list.end());
	unsigned hash_subgraph_code = HashFunc(edge_list, mHashBitMask);
	return hash_subgraph_code;
}

unsigned NSPDK_FeatureGenerator::GraphCanonicalFormEncoding(const GraphClass& aG) {
	GraphClass g(aG); //make a copy of the graph
	g.ComputePairwiseDistanceInformation(g.VertexSize()); //compute complete breadth first visit

	//for all vertices extract the vertex's signature: sorted distance from all the other vertices + their vertex label
	vector<unsigned> vertex_encoding_list;
	for (unsigned i = 0; i < g.VertexSize(); ++i) {
		vector<unsigned> vertex_encoding;
		vector<unsigned> vertex_distance_list;
		for (unsigned j = 0; j < g.VertexSize(); ++j) {
			int dist = g.PairwiseDistance(i, j);
			string distance_label = stream_cast<string>(dist) + g.GetVertexLabelConcatenated(j);
			unsigned hash_distance_label = HashFunc(distance_label, mHashBitMask);
			vertex_distance_list.push_back(hash_distance_label);
		}
		sort(vertex_distance_list.begin(), vertex_distance_list.end());
		unsigned hash_encoding = HashFunc(vertex_distance_list, mHashBitMask);
		vertex_encoding_list.push_back(hash_encoding);
	}
	//extract list of all edge's signatures in the induced graph: v-signature,u-signature,label(uv)
	vector<unsigned> edge_list;
	vector<unsigned> edge_encoding(3);
	for (unsigned u = 0; u < g.VertexSize(); ++u) {
		//get all edges of vertex u
		vector<unsigned> vertex_adjacency_list = g.GetVertexAdjacentList(u);
		vector<unsigned> edge_adjacency_list = g.GetEdgeAdjacentList(u);
		if (vertex_adjacency_list.size() == 0) { //NOTE: if a vertex is isolated use its encoding as edge encoding
			edge_list.push_back(vertex_encoding_list[u]);
		}
		for (unsigned j = 0; j < vertex_adjacency_list.size(); ++j) {
			unsigned v = vertex_adjacency_list[j];

			if (vertex_encoding_list[u] < vertex_encoding_list[v]) {
				edge_encoding[0] = vertex_encoding_list[u];
				edge_encoding[1] = vertex_encoding_list[v];
			} else {
				edge_encoding[0] = vertex_encoding_list[v];
				edge_encoding[1] = vertex_encoding_list[u];
			}
			unsigned e = edge_adjacency_list[j];
			string edge_label = g.GetEdgeLabelConcatenated(e);
			unsigned hash_edge_label = HashFunc(edge_label, mHashBitMask);
			edge_encoding[2] = hash_edge_label;
			unsigned hash_edge_encoding = HashFunc(edge_encoding, mHashBitMask);
			edge_list.push_back(hash_edge_encoding);
		}
	}
	//the graph encoding is the sorted list of edge encodings
	sort(edge_list.begin(), edge_list.end());
	unsigned hash_subgraph_code = HashFunc(edge_list, mHashBitMask);
	return hash_subgraph_code;
}

unsigned NSPDK_FeatureGenerator::RootedGraphCanonicalFormEncoding(const GraphClass& aG, unsigned aRootID) {
	//for all vertices extract the vertex's signature: distance from root-sorted distance from all the other vertices + their vertex label
	vector<unsigned> vertex_encoding_list;
	for (unsigned i = 0; i < aG.VertexSize(); ++i) {
		vector<unsigned> vertex_distance_list;
		int dist = aG.PairwiseDistance(aRootID, i);
		string distance_label = stream_cast<string>(dist) + aG.GetVertexLabelConcatenated(aRootID);
		unsigned hash_distance_label = HashFunc(distance_label, mHashBitMask);
		vertex_distance_list.push_back(hash_distance_label);

		for (unsigned j = 0; j < aG.VertexSize(); ++j) {
			int dist = aG.PairwiseDistance(i, j);
			string distance_label = stream_cast<string>(dist) + aG.GetVertexLabelConcatenated(j);
			unsigned hash_distance_label = HashFunc(distance_label, mHashBitMask);
			vertex_distance_list.push_back(hash_distance_label);
		}
		sort(vertex_distance_list.begin(), vertex_distance_list.end());
		unsigned hash_encoding = HashFunc(vertex_distance_list, mHashBitMask);
		vertex_encoding_list.push_back(hash_encoding);
	}
	//extract list of all edge's signatures in the induced graph: v-signature,u-signature,label(uv)
	vector<unsigned> edge_list;
	vector<unsigned> edge_encoding(3);
	for (unsigned u = 0; u < aG.VertexSize(); ++u) {
		//get all edges of vertex u
		vector<unsigned> vertex_adjacency_list = aG.GetVertexAdjacentList(u);
		vector<unsigned> edge_adjacency_list = aG.GetEdgeAdjacentList(u);
		if (vertex_adjacency_list.size() == 0) { //NOTE: if a vertex is isolated use its encoding as edge encoding
			edge_list.push_back(vertex_encoding_list[u]);
		}
		for (unsigned j = 0; j < vertex_adjacency_list.size(); ++j) {
			unsigned v = vertex_adjacency_list[j];

			if (vertex_encoding_list[u] < vertex_encoding_list[v]) {
				edge_encoding[0] = vertex_encoding_list[u];
				edge_encoding[1] = vertex_encoding_list[v];
			} else {
				edge_encoding[0] = vertex_encoding_list[v];
				edge_encoding[1] = vertex_encoding_list[u];
			}
			unsigned e = edge_adjacency_list[j];
			string edge_label = aG.GetEdgeLabelConcatenated(e);
			unsigned hash_edge_label = HashFunc(edge_label, mHashBitMask);
			edge_encoding[2] = hash_edge_label;
			unsigned hash_edge_encoding = HashFunc(edge_encoding, mHashBitMask);
			edge_list.push_back(hash_edge_encoding);
		}
	}
	//the graph encoding is the sorted list of edge encodings
	sort(edge_list.begin(), edge_list.end());
	unsigned hash_subgraph_code = HashFunc(edge_list, mHashBitMask);
	return hash_subgraph_code;
}

void NSPDK_FeatureGenerator::GetFirstEndpoints(const GraphClass& aG, vector<unsigned>& oFirstEndpointList) const {
	//insert additional vertices
	if (oFirstEndpointList.size() == 0) //if oFirstEndpointList is empty then fill it with all vertices that are viewpoints, otherwise do nothing i.e. use the given list
		for (unsigned i = 0; i < aG.VertexSize(); ++i)
			if (aG.GetVertexViewPoint(i) || aG.GetVertexAbstraction(i))
				oFirstEndpointList.push_back(i);
	if (mDebugVerbosity > 0) {
		cout << "First endpoint id list [" << oFirstEndpointList.size() << "]:" << endl;
		for (unsigned i = 0; i < oFirstEndpointList.size(); i++)
			cout << oFirstEndpointList[i] << " ";
		cout << endl;
	}
}

inline
unsigned NSPDK_FeatureGenerator::HashFunc(const string& aString, unsigned aBitMask) { //NOTE: extract the least significant bits from the hash
	unsigned int hash = 0xAAAAAAAA;
	for (std::size_t i = 0; i < aString.length(); i++) {
		hash ^= ((i & 1) == 0) ? ((hash << 7) ^ aString[i] * (hash >> 3)) : (~(((hash << 11) + aString[i]) ^ (hash >> 5)));
	}
	return hash & aBitMask;
}

inline
unsigned NSPDK_FeatureGenerator::HashFunc(const vector<unsigned>& aList, unsigned aBitMask) {
	unsigned int hash = 0xAAAAAAAA;
	for (std::size_t i = 0; i < aList.size(); i++) {
		hash ^= ((i & 1) == 0) ? ((hash << 7) ^ aList[i] * (hash >> 3)) : (~(((hash << 11) + aList[i]) ^ (hash >> 5)));
	}
	return hash & aBitMask;
}

//----------------------------------------------------------------------------------------------------------------------------------------------
String_FeatureGenerator::String_FeatureGenerator(const std::string& id) :
		NSPDK_FeatureGenerator(id) {
}

void String_FeatureGenerator::generate_feature_vector(const GraphClass& aG, SVector& x, const vector<unsigned>& aFirstEndpointList) {
	//assume 1 vertex with all info on the label
	string seq = aG.GetVertexLabel(0);
	unsigned size = seq.size();

	if (aFirstEndpointList.size() == 0) {
		InitFeatureCache(size, mRadius);

		//create neighborhood features
		for (unsigned start = 0; start < size; ++start)
			mFeatureCache[start] = HashFunc(seq, start, mRadius, mHashBitMask);

		vector<unsigned> endpoint_list(4);
		for (unsigned r = 0; r <= mRadius; r++) {
			endpoint_list[0] = r;
			for (unsigned d = 0; d <= mDistance; d++) {
				endpoint_list[1] = d;
				SVector z(pow(2, mHashBitSize));
				for (unsigned start = 0; start < size; ++start) {
					unsigned src_code = mFeatureCache[start][r];
					unsigned effective_dest = min(start + d, size - 1);
					unsigned dest_code = mFeatureCache[effective_dest][r];
					//impose canonical order for pair: i.e. A-B and B-A must generate the same feature
					if (src_code < dest_code) {
						endpoint_list[2] = src_code;
						endpoint_list[3] = dest_code;
					} else {
						endpoint_list[2] = dest_code;
						endpoint_list[3] = src_code;
					}
					unsigned code = HashFunc(endpoint_list, mHashBitMask);
					z.coeffRef(code) += 1;

					//add feature that ignore src endpoint
					//NOTE: this is important when the label of src vertex is considered noisy, in this way, it is only the context that will define the features
					endpoint_list[2] = 0;
					endpoint_list[3] = dest_code;
					unsigned nosrc_code = HashFunc(endpoint_list, mHashBitMask);
					z.coeffRef(nosrc_code) += 1;
				}
				if (mNormalization)
					z /= z.norm();
				x += z;
			}
		}
		if (mNormalization)
			x /= x.norm();
	} else
		throw range_error("ERROR String_FeatureGenerator::generate_feature_vector: expecting an empty end point list");
}

void String_FeatureGenerator::InitFeatureCache(unsigned aSize, unsigned aRadius) {
	mFeatureCache.clear();
	for (unsigned i = 0; i < aSize; ++i) {
		mFeatureCache.push_back(vector<unsigned>(aRadius, 0));
	}
}

void String_FeatureGenerator::generate_vertex_feature_vector(const GraphClass& aG, vector<SVector>& x_list, const vector<unsigned>& aFirstEndpointList) {
	//assume 1 vertex with all info on the label
	string seq = aG.GetVertexLabel(0);
	unsigned size = seq.size();

	if (aFirstEndpointList.size() != 0)
		throw range_error("ERROR String_FeatureGenerator::generate_feature_vector: expecting an empty end point list");

	InitFeatureCache(size, mRadius);

	//create neighborhood features
	for (unsigned start = 0; start < size; ++start) {
		mFeatureCache[start] = HashFunc(seq, start, mRadius, mHashBitMask);
	}

	vector<unsigned> endpoint_list(4);
	for (unsigned start = 0; start < size; ++start) {
		SVector z(pow(2, mHashBitSize));
		for (unsigned r = 0; r <= mRadius; r++) {
			endpoint_list[0] = r;
			for (unsigned d = 0; d <= mDistance; d++) {
				endpoint_list[1] = d;
				unsigned src_code = mFeatureCache[start][r];
				unsigned effective_dest = min(start + d, size - 1);
				unsigned dest_code = mFeatureCache[effective_dest][r];
				//impose canonical order for pair: i.e. A-B and B-A must generate the same feature
				if (src_code < dest_code) {
					endpoint_list[2] = src_code;
					endpoint_list[3] = dest_code;
				} else {
					endpoint_list[2] = dest_code;
					endpoint_list[3] = src_code;
				}
				unsigned code = HashFunc(endpoint_list, mHashBitMask);
				z.coeffRef(code) += 1;

				//add feature that ignore src endpoint
				//NOTE: this is important when the label of src vertex is considered noisy, in this way, it is only the context that will define the features
				endpoint_list[2] = 0;
				endpoint_list[3] = dest_code;
				unsigned nosrc_code = HashFunc(endpoint_list, mHashBitMask);
				z.coeffRef(nosrc_code) += 1;
			}
		}
		x_list.push_back(z);
	}
}

inline vector<unsigned> String_FeatureGenerator::HashFunc(const string& aString, unsigned aStart, unsigned aMaxRadius, unsigned aBitMask) {
	unsigned int hash = 0xAAAAAAAA;
	unsigned effective_end = min((unsigned) aString.size() - 1, aStart + aMaxRadius);
	unsigned radius = 0;
	vector<unsigned> code_list(aMaxRadius + 1, 0);
	for (std::size_t i = aStart; i <= effective_end; i++) {
		hash ^= ((i & 1) == 0) ? ((hash << 7) ^ aString[i] * (hash >> 3)) : (~(((hash << 11) + aString[i]) ^ (hash >> 5)));
		code_list[radius] = hash & aBitMask;
		radius++;
	}
	return code_list;
}

unsigned String_FeatureGenerator::HashFunc(const vector<unsigned>& aList, unsigned aBitMask) {
	unsigned int hash = 0xAAAAAAAA;
	for (std::size_t i = 0; i < aList.size(); i++) {
		hash ^= ((i & 1) == 0) ? ((hash << 7) ^ aList[i] * (hash >> 3)) : (~(((hash << 11) + aList[i]) ^ (hash >> 5)));
	}
	return hash & aBitMask;
}

//----------------------------------------------------------------------------------------------------------------------------------------------

#ifndef EDeNGRAPH_H
#define EDeNGRAPH_H

#include "Utility.h"

namespace EDeN {

	const string GRAPH_IDENTIFIER_KEY = "__gid__";
	const string VERTEX_IDENTIFIER_KEY = "__vid__";
	const string EDGE_SRC_IDENTIFIER_KEY = "__srcid__";
	const string EDGE_DEST_IDENTIFIER_KEY = "__destid__";
	const string NODE_TYPE_KEY = "__type__";
	const string VERTEX_TYPE_VALUE = "__vertex__";
	const string EDGE_TYPE_VALUE = "__edge__";

	const string LABEL_KEY = "L";
	const string SVECTOR_KEY = "svec";
	const string DVECTOR_KEY = "dvec";
	const string WEIGHT_KEY = "w";
	const string AGGREGATE_KEY = "aggregate";

	const string DEFAULT_VERTEX_LABEL = "__v__";
	const string DEFAULT_EDGE_LABEL = "__e__";
	const string DEFAULT_VERTEX_ID_NAMESPACE = "__g__";

	///
	/// The vertex data structure is a dictionary of different primitive data types.
	/// These include: numbers, discrete (as strings), booleans and sparse or dense vectors
	///
	struct Node {
			map<string, double> mNumber;
			map<string, string> mDiscrete;
			map<string, bool> mBool;
			map<string, SVector> mSVector;
			map<string, FVector> mDVector;
	};

	///
	/// Design choices:
	/// 1) the parser class knows only the Node structure
	///
	class Parser {
		public:
			Parser() {
			}
			istream& TestForChar(istream& in, char aChar);
			bool IsValid(istream& in);
			bool ParseRoot(istream& in, Node& oRoot);
			bool ParseVertex(istream& in, Node& oNode);
			bool ParseEdge(istream& in, Node& oNode);
			bool ParseAttributes(istream& in, Node& oNode);
			bool ParseAttributeValuePair(istream& aIn, Node& oNode);
	};

	///
	/// Design choices:
	/// 1) Edges and hyperedges are transformed into nodes (internally) so that all labels are on vertices
	///
	class Graph {
		protected:
			Node mRoot;
			vector<Node> mNodes;
			vector<list<unsigned> > mAdjacencyList;
			mutable vector<vector<list<unsigned> > > mNeighborsCache;
			mutable vector<unsigned> mDestMaptoDistance;
			Parser mParser;
		public:
			Graph() {
			}
			//convert the expected format in a JSON format and then use the JSON parser
			//group all vertices in a vector of vertices for faster access
			//same with edges
			//eg. string format
			Graph(istream& aI, const string& aFormat) {
				Parse(aI, aFormat);
			}
			bool IsGood() const {
				return mNodes.size() > 0;
			}
			unsigned Size() const {
				return mNodes.size();
			}
			//have specialized getters (i.e. they check and extract specific attribute values)
			string GetLabel(unsigned aVertexID, const string aKey = LABEL_KEY) const;
			double GetWeight(unsigned aVertexID, const string aKey = WEIGHT_KEY) const;
			bool HasSVector(unsigned aVertexID, const string aKey = SVECTOR_KEY) const;
			bool HasDVector(unsigned aVertexID, const string aKey = DVECTOR_KEY) const;
			unsigned SizeSVector(unsigned aVertexID, const string aKey = SVECTOR_KEY) const;
			unsigned SizeDVector(unsigned aVertexID, const string aKey = DVECTOR_KEY) const;
			SVector GetSVector(unsigned aVertexID, const string aKey = SVECTOR_KEY) const;
			FVector GetDVector(unsigned aVertexID, const string aKey = DVECTOR_KEY) const;

			//return the list of all ids of vertices that are at distance 0,1,2,.. in a vector
			//cache the returned lists under lazy style framework
			void GetNeighbors(unsigned aSrcVertexID, unsigned aDistance, vector<list<unsigned> >& oNeighborDistanceList) const;
			void GetNeighbors(unsigned aSrcVertexID, unsigned aDistance, list<unsigned>& oNeighborsList) const;
		protected:
			void Parse(istream& aI, const string& aFormat);
			void UpdateAdjacencyList(const Node& aEdgeNode);//TODO
			bool IsTraversable(unsigned aVertexID) const;//TODO
			void SingleVertexBoundedBreadthFirstVisit(unsigned aSrcVertexID, unsigned aDistance) const;
	};

	///
	/// Design choices: neighborhood extractions are specialized for specific values of radius and distances to increase efficiency.
	/// The kernel for string data type is specialized.
	/// The design has template policies to manage different attribute types and provide compile time type checking
	/// Note that having both the graph class and the kernel class are in the same translation unit
	/// ensures ease of maintenance.
	///
	class Kernel {
		protected:
			unsigned mRadius;
			unsigned mDistance;
			unsigned mHashBitSize;
			unsigned mNumProjections;

			mutable multimap<string, SVector> mAttributeFeatureVectorDictionary;
			mutable vector<vector<unsigned> > mIdRadiusHashCodeCache;
		public:
			Kernel() :
					mRadius(2), mDistance(5), mHashBitSize(15), mNumProjections(1) {
			}
			void GenerateFeatureVector(const Graph& aG, SVector& oX);
			void GenerateVertexFeatureVector(const Graph& aG, vector<SVector>& oXList); //TODO
			double ComputeKernel(const Graph& aG, const Graph& aM);
			double ComputeKernel(const SVector& aX, const SVector& aZ);
		protected:
			void AggregateFeatures(string aKey, SVector& oX) const;
			unsigned HashFunc(const vector<unsigned>& aList, unsigned aBitMask = 2147483647) const;
			vector<unsigned> HashFuncVec(const vector<unsigned>& aList, unsigned aBitMask = 2147483647) const;
			unsigned HashFuncString(const string& aString, unsigned aBitMask = 2147483647) const;
			void SelectVertexListPolicy(const Graph& aG, list<unsigned>& oVetexList) const;
			void NormalizationPolicy(SVector& oX) const;
			void ComputeFeatures(const Graph& aG, unsigned aID) const;
			void CacheNeighborhoodGraphHashCode(const Graph& aG, unsigned aID, unsigned aMaxRadius) const;
			unsigned ComputeNeighborhoodGraphHashCode(const Graph& aG, unsigned aID, unsigned aRadius) const;
			string GenerateKeyRadiusDistance(unsigned aRadius, unsigned aDistance) const;
			void ComputeFeaturesRadiusDistance(const Graph& aG, unsigned aID) const;
			void ComputeFeaturesRadiusDistance(const Graph& aG, unsigned aID, unsigned aRadius, unsigned aDistance, SVector& oX) const;
			void ComputeFeaturesNonRootRadiusDistance(const Graph& aG, unsigned aID) const;
			void ComputeFeaturesNonRootRadiusDistance(const Graph& aG, unsigned aID, unsigned aRadius, unsigned aDistance, SVector& oX) const;
			void ComputeFeaturesDenseRealVector(const Graph& aG, unsigned aID) const;
			void ComputeFeaturesDenseRealVector(const Graph& aG, unsigned aID, unsigned aRadius, unsigned aDistance, SVector& oX) const;
			void ComputeFeaturesSparseRealVector(const Graph& aG, unsigned aID) const;
			void ComputeFeaturesSparseRealVector(const Graph& aG, unsigned aID, unsigned aRadius, unsigned aDistance, SVector& oX) const;
			double PositiveOrthantRandomProjection(const FVector& aRealVector, unsigned aCode) const;
			double RandomProjection(const FVector& aRealVector, unsigned aCode) const;
			double PositiveOrthantRandomProjection(const SVector& aRealVector, unsigned aCode) const;
			double RandomProjection(const SVector& aRealVector, unsigned aCode) const;
			double MinHashProjection(const SVector& aRealVector, unsigned aCode) const;
			double MinHashProjection(const FVector& aRealVector, unsigned aCode) const;
	};
}
#endif

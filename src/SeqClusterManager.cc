#include "SeqClusterManager.h"
#include <ostream>


//---------------------------------------------------------------------------------------------------
// Sequence clustering via approximate neighborhood density
SeqClusterManager::SeqClusterManager(Parameters* apParameters, Data* apData)
: NeighborhoodIndex(apParameters, apData)
{
	mMinHashCache = std::make_shared<SigCacheT>();
	mIndexDataSet = std::make_shared<SeqFileT>();

	mIndexDataSet->filename = mpParameters->mInputDataFileName.c_str();
	mIndexDataSet->filename_BED = "";
	mIndexDataSet->filetype = mpParameters->mFileTypeCode;
	mIndexDataSet->groupGraphsBy=SEQ_WINDOW;
	mIndexDataSet->sigCache = mMinHashCache;
	mIndexDataSet->lastMetaIdx = 0;
	mIndexDataSet->checkUniqueSeqNames=true;
	mIndexDataSet->signatureAction	= INDEX_SIGCACHE;
	mIndexDataSet->strandType			= FWD;

	SeqFilesT myList;
	myList.push_back(mIndexDataSet);
	LoadData_Threaded(myList);

	NeighborhoodCacheReset();
}

void SeqClusterManager::finishUpdate(ChunkP& myData) {

	switch ((*myData)[0].seqFile->signatureAction){
	case INDEX_SIGCACHE: {
		unsigned chunkSize 	= myData->size();
		unsigned offset 		= (*myData)[0].idx-1;	//first chunk has idx=1
		//cout << "offset " << offset << " chunk " << chunkSize << endl;
		if ((*myData)[0].seqFile->sigCache->size() < offset + chunkSize) {
			(*myData)[0].seqFile->sigCache->resize( offset + chunkSize);
			idx2nameMap.resize(offset + chunkSize);
		}

		for (unsigned j = 0; j < chunkSize; j++) {
			(*myData)[0].seqFile->sigCache->at(offset+j) = (*myData)[j].sig;
			name2idxMap.insert(make_pair((*myData)[j].name, offset+j));
			idx2nameMap.at(offset+j) = (*myData)[j].name;
		}
	}
	case INDEX: {
			for (unsigned j = 0; j < myData->size(); j++) {
				UpdateInverseIndex((*myData)[j].sig, (*myData)[j].idx);
			}
			break;
	}
	default:
		break;
	}
}


void SeqClusterManager::Exec() {
	ProgressBar pb(1000);

	DenseCluster();

	pb.Count(1);
	cout << endl << "Total clustering time:" << endl;
}

void SeqClusterManager::DenseCluster() {
	vector<unsigned> mDenseCenterIdxList;
	if (mpParameters->mDenseCenterNamesFile != ""){
		vector<string> centerNames;
		mpData->LoadStringList(mpParameters->mDenseCenterNamesFile,centerNames,1);
		for (uint i=0;i<centerNames.size();i++){
			//cout << "name " << centerNames[i] << " map: "<< mNearestNeighbor.mMinHashEncoder.name2idxMap[centerNames[i]] << endl;
			mDenseCenterIdxList.push_back(name2idxMap[centerNames[i]]);
		}

	} else {

		for (uint i=0; i<GetLoadedInstances(); i++){
			mDenseCenterIdxList.push_back(i);
		}
	}

	cout << endl << "Compute neighborhood and density for selected " << mDenseCenterIdxList.size() << " instances." << endl;
	vector<pair<double, unsigned> > DensityList(mDenseCenterIdxList.size());
	unsigned fullCollisions = 0;

	{
		ProgressBar ppb(1000);
#ifdef USEMULTITHREAD
#pragma omp parallel for schedule(dynamic,100)
#endif
		for (unsigned i = 0; i < mDenseCenterIdxList.size(); ++i) {

			uint ii = mDenseCenterIdxList[i];
			//compute neighbors
			vector<unsigned> neighborhood_list;
			unsigned collisions = 0;
			double density = 0;
			neighborhood_list = ComputeNeighborhood(ii,collisions,density);
			double score = density * neighborhood_list.size();
			DensityList[i]= make_pair(-score,ii);
			if (collisions >= mpParameters->mNumHashFunctions) fullCollisions++;
			//cout << i << " dens " << density << " size " << neighborhood_list.size() << " coll " << collisions <<  endl;
			ppb.Count();
		}
	}
	cout << endl << " Instances with complete signature collision (MaxSizeBin): " << fullCollisions << endl;

	cout << endl << " *** Approx-Density Clustering *** " << endl;
	cout << "sort instances according to neighborhood density score..." << endl;

	sort(DensityList.begin(), DensityList.end());

	vector<unsigned> neighborhood_list;
	unsigned collisions = 0;
	double density = 0;
	vector<long int> adjacencyList(GetLoadedInstances());

	for (unsigned k = 0; k < adjacencyList.size(); ++k) {
		adjacencyList[k] = -1;
	}

	cout << "make adjacency list..." << endl;
	for (	unsigned i = 0; i < DensityList.size(); ++i) {
		unsigned ii = DensityList[i].second;
		// cluster seed seq is not yet (-1) part of any cluster
		if ( adjacencyList[ii] == -1 ) {
			neighborhood_list = ComputeNeighborhood(ii,collisions,density);

			// count not clustered instances of neighborhood(ii)
			// and also the cluster-frequency of already clustered instances
			umap_uint_int clustered;
			unsigned notClustered = 0;
			for (vector<unsigned>::const_iterator it = neighborhood_list.begin(); it != neighborhood_list.end(); ++it) {
				if ( adjacencyList[*it] != -1 ) {
					if (clustered.count(adjacencyList[*it]) == 0 ){
						clustered.insert(make_pair(adjacencyList[*it],1));
					} else {
						clustered[adjacencyList[*it]]++;
					}
				} else {
					notClustered++;
				}
			}

			unsigned newClusId = ii;
			if (notClustered<neighborhood_list.size()){
				vector<pair<signed,unsigned> > sortedAb;
				for (umap_uint_int::iterator it = clustered.begin(); it != clustered.end(); it++){
					sortedAb.push_back(make_pair(-it->second,it->first));
				}

				sort(sortedAb.begin(),sortedAb.end());

				// merge not clustered instances of nh(ii) into most abundant cluster (sortedAb[0].second) of nh(ii) if:
				// a) all already clustered instances belong to the same cluster
				// b) neighborhood size is less than half of the instances that are not yet clustered
				if (-sortedAb.front().first >= (signed)(neighborhood_list.size() - notClustered) && (notClustered+2)*2 <= neighborhood_list.size()){
					newClusId = sortedAb[0].second;
				}
			}

			// assign cluster id to all not clustered instances
			for (vector<unsigned>::const_iterator it = neighborhood_list.begin(); it != neighborhood_list.end(); ++it) {
				if ( adjacencyList[*it] == -1 ) {
					adjacencyList[*it] = newClusId;
				}
			}
			adjacencyList[ii] = newClusId;
		}
	}

	// build clusters
	umap_uint_vec_uint clusters;
	for (	unsigned i = 0; i < adjacencyList.size(); ++i) {
		if ( adjacencyList[i] != -1 ) {
			clusters[adjacencyList[i]].push_back(i);
		} else if (mpParameters->mDenseCenterNamesFile == ""){
			clusters[i].push_back(i);
		}
	}

	cout << "Found " << clusters.size() << " clusters - output them now ..." << endl;

	// OUTPUT ---------


	OutputManager out_c("approx_dense_cluster", mpParameters->mDirectoryPath);

	for (umap_uint_vec_uint::const_iterator it = clusters.begin(); it != clusters.end(); ++it) {
		out_c.mOut << it->first << "    ";
		for (vector<unsigned>::const_iterator mIt = it->second.begin(); mIt != it->second.end(); ++mIt) {
			out_c.mOut << *mIt << " ";
		}
		out_c.mOut << endl;
	}
	cout << endl << "Approx dense cluster results written in file " << out_c.GetFullPathFileName() << endl;

	if (mpParameters->mWriteApproxNeighbors){
		//output all info: id of neighbour, similarity of neighbor, target of neighbor
		OutputManager out_n("approx_dense_centers", mpParameters->mDirectoryPath);
		cout << "Output all ranked Neighborhoods..." << endl;
		for (	unsigned i = 0; i < DensityList.size(); ++i) {
			unsigned ii = DensityList[i].second;
			out_n.mOut << ii << ":" << -DensityList[i].first <<  "    ";
			neighborhood_list = ComputeNeighborhood(ii,collisions,density);
			for (vector<unsigned>::const_iterator it = neighborhood_list.begin(); it != neighborhood_list.end(); ++it) {
				out_n.mOut << *it << " ";
			}
			out_n.mOut << endl;
		}
		cout << "Approx dense center results written in file " << out_n.GetFullPathFileName() << endl;
	}
}

//---------------------------------------------------------------------------------------------------

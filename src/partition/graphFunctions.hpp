/*
  Copyright 2018 Andreas Thune.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_GRAPHFUNCTIONS_HEADER_INCLUDED
#define OPM_GRAPHFUNCTIONS_HEADER_INCLUDED

#endif // OPM_GRAPHFUNCTIONS_HEADER_INCLUDED

#include <algorithm>
#include <map>
#include <queue>

#include "transWellgraph.hpp"
#include "hyperGraphFunctions.hpp"
#include "parmetisPart.hpp"

// Num object NULL
int getNullNumCells(void* graphPointer, int* err)
{
    (void) graphPointer;
    *err = ZOLTAN_OK;
    return 0;
}

// Object List NULL
void getNullVertexList(void* graphPointer, int numGlobalIdEntries,
                       int numLocalIdEntries, ZOLTAN_ID_PTR gids,
                       ZOLTAN_ID_PTR lids, int wgtDim,
                       float *objWgts, int *err)
{
    (void) graphPointer; (void) numGlobalIdEntries;
    (void) numLocalIdEntries; (void) gids; (void) lids; (void) objWgts;
    (void) wgtDim;
    // We do nothing as we pretend to not have any grid cells.
    *err = ZOLTAN_OK;
}

// Num Edges NULL
void getNullNumEdgesList(void *cpGridPointer, int sizeGID, int sizeLID,
			 int numCells,
			 ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
			 int *numEdges, int *err)
{
    (void) sizeGID; (void) sizeLID; (void) numCells; (void) globalID;
    (void) localID; (void) numEdges; (void) cpGridPointer;
    // Pretend that there are no edges
    numEdges = 0;
    *err = ZOLTAN_OK;
}

// Edges List NULL
void getNullEdgeList(void *cpGridPointer, int sizeGID, int sizeLID,
		     int numCells, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
		     int *numEdges,
		     ZOLTAN_ID_PTR nborGID, int *nborProc,
		     int wgtDim, float *ewgts, int *err)
{
    (void) cpGridPointer; (void) sizeGID; (void) sizeLID; (void) numCells;
    (void) globalID; (void) localID; (void) numEdges; (void) nborGID;
    (void) nborProc; (void) wgtDim; (void) ewgts;
    *err = ZOLTAN_OK;
}

// Num objects 
int getMatNumCells(void* graphPointer, int* err)
{
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> Mat;
    typedef TransWellGraph<Mat> Graph;

    const Graph& graph = *static_cast<const Graph*>(graphPointer);
    const Mat& transMat = graph.edgeWgt_;

    *err = ZOLTAN_OK;
    return transMat.N();
}


// Num objects 2 Coarse graph
int getMatNumCellsCoarse(void* graphPointer, int* err)
{
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> Mat;
    typedef TransWellGraph<Mat> Graph;

    const Graph& graph = *static_cast<const Graph*>(graphPointer);
    const std::vector<std::vector<int>>& nodes = graph.courseNodes_;

    *err = ZOLTAN_OK;
    return nodes.size();
}

// Object List 
void getMatVertexList(void* graphPointer, int numGlobalIdEntries,
		      int numLocalIdEntries, ZOLTAN_ID_PTR gids,
		      ZOLTAN_ID_PTR lids, int wgtDim,
		      float *objWgts, int *err)
{
    //(void) wgtDim; (void) objWgts;
    
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> Mat;
    typedef TransWellGraph<Mat> Graph;

    const Graph& graph = *static_cast<const Graph*>(graphPointer);
    const Mat& transMat = graph.edgeWgt_;
    const std::vector<int>& rs = graph.row_size;

    for (int idx = 0; idx < transMat.N(); ++idx)
    {
        gids[idx] = idx;
        lids[idx] = idx;
	if (wgtDim == 1)
	    objWgts[idx] = rs[idx];
	if (wgtDim == 2) {
	    objWgts[2 * idx ]     = 1;
	    objWgts[2 * idx + 1 ] = rs[idx];
	}
    }

    *err = ZOLTAN_OK;
}

// Object List 2 coarse graph
void getMatVertexListCoarse(void* graphPointer, int numGlobalIdEntries,
		      int numLocalIdEntries, ZOLTAN_ID_PTR gids,
		      ZOLTAN_ID_PTR lids, int wgtDim,
		      float *objWgts, int *err)
{
    //(void) wgtDim; (void) objWgts;

    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> Mat;    
    typedef TransWellGraph<Mat> Graph;

    const Graph& graph = *static_cast<const Graph*>(graphPointer);
    const std::vector<std::vector<int>>& nodes = graph.courseNodes_;
    

    for (int idx = 0; idx < nodes.size(); ++idx)
    {
	//std::cout << "in getMatVertexListCoarse idx: " << idx << std::endl;
        gids[idx] = idx;
        lids[idx] = idx;
	if (wgtDim == 1)
	    objWgts[idx] = nodes[idx].size();
    }

    *err = ZOLTAN_OK;
}

// Num edges
void getMatNumEdgesList(void *graphPointer, int sizeGID, int sizeLID,
			int numCells,
			ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
			int *numEdges, int *err)
{
    (void) globalID;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> Mat;
    typedef TransWellGraph<Mat> Graph;

    const Graph& graph = *static_cast<const Graph*>(graphPointer);
    const Mat& transMat = graph.edgeWgt_;

    for (auto row = transMat.begin(); row != transMat.end(); ++row)
    {
	int edges = 0;
	
	auto col = row->begin();
	for (; col!=row->end(); ++col)
	{
	    edges++;
	}
	numEdges[row.index()] = edges;
    }
    
    *err = ZOLTAN_OK;
}

// Num edges
void getMatWellNumEdgesList(void *graphPointer, int sizeGID, int sizeLID,
			    int numCells,
			    ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
			    int *numEdges, int *err)
{
    (void) globalID;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> Mat;
    typedef TransWellGraph<Mat> Graph;

    const Graph& graph = *static_cast<const Graph*>(graphPointer);
    const Mat& transMat = graph.edgeWgt_;
    const Mat& wellsMat = graph.wells;
    
    auto trow = transMat.begin();
    auto wrow = wellsMat.begin();
    for (; trow != transMat.end(); trow++, wrow++)
    {
	int edges = 0;
	
	auto wcol = wrow->begin();
	for (; wcol != wrow->end(); ++wcol)
	{
	    edges++;
	}
	auto tcol = trow->begin();
	for (; tcol != trow->end(); ++tcol)
	{
	    edges++;
	}

	numEdges[trow.index()] = edges;
    }
    
    *err = ZOLTAN_OK;
}

// Num edges 3 coarse graph
void getMatNumEdgesListCoarse(void *graphPointer, int sizeGID, int sizeLID,
			      int numCells, ZOLTAN_ID_PTR globalID,
			      ZOLTAN_ID_PTR localID, int *numEdges, int *err)
{
    (void) globalID;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> Mat;
    typedef TransWellGraph<Mat> Graph;

    const Graph& graph = *static_cast<const Graph*>(graphPointer);
    const std::vector<std::map<int, double> >& edges = graph.cedges;
    

    for (int idx = 0; idx < edges.size(); ++idx)
    {

	//std::cout << "in getMatNumEdgesListCoarse idx: " << idx << " "<< edges.size() << std::endl;
	numEdges[idx] = edges[idx].size();
    }
    
    *err = ZOLTAN_OK;
}

// Edge List 1, No Weights
void getMatEdgeList(void *graphPointer, int sizeGID, int sizeLID,
		    int numCells, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
		    int *numEdges,
		    ZOLTAN_ID_PTR nborGID, int *nborProc,
		    int wgtDim, float *ewgts, int *err)
{
    (void) wgtDim; (void) globalID; (void) numEdges; (void) ewgts;

    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> Mat;
    typedef TransWellGraph<Mat> Graph;

    const Graph& graph = *static_cast<const Graph*>(graphPointer);
    const Mat& transMat = graph.edgeWgt_;
    const Mat& wellsMat = graph.wells;

    int idx = 0;    
    auto trow = transMat.begin();
    for (; trow != transMat.end(); ++trow)
    {
	auto tcol = trow->begin();
	for (; tcol != trow->end(); ++tcol)
	{    
	    if (tcol.index() != trow.index()) {
		nborGID[idx++] = tcol.index();
	    }
	}	
    }
    for ( int i = 0; i < idx; ++i )
    {
        nborProc[i] = 0;
    }
}

// Edge List 2, wells + uniform
void getWellMatEdgeList(void *graphPointer, int sizeGID, int sizeLID,
		    int numCells, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
		    int *numEdges,
		    ZOLTAN_ID_PTR nborGID, int *nborProc,
		    int wgtDim, float *ewgts, int *err)
{
    (void) wgtDim; (void) globalID; (void) numEdges; (void) ewgts;

    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> Mat;
    typedef TransWellGraph<Mat> Graph;
    
    const Graph& graph = *static_cast<const Graph*>(graphPointer);
    const Mat& transMat = graph.edgeWgt_;
    const Mat& wellsMat = graph.wells;

    int idx = 0;
    auto trow = transMat.begin();
    auto wrow = wellsMat.begin();
    for (; trow != transMat.end(); ++trow, ++wrow)
    {	
	auto wcol = wrow->begin();
	for (; wcol != wrow->end(); ++wcol)
	{
	    nborGID[idx] = wcol.index();
	    ewgts[idx++] = std::numeric_limits<float>::max();//1.0e10;//std::numeric_limits<float>::max();
	}

	auto tcol = trow->begin();
	for (; tcol != trow->end(); ++tcol)
	{    
	    if (tcol.index() != trow.index()) {
		
		nborGID[idx] = tcol.index();
		ewgts[idx++] = 1.0;
	    }
	}	
    }
    for ( int i = 0; i < idx; ++i )
    {
        nborProc[i] = 0;
    }
}

// Edge List 3, trans weights no wells 
void getWeightMatEdgeList(void *graphPointer, int sizeGID, int sizeLID,
			  int numCells, ZOLTAN_ID_PTR globalID, 
			  ZOLTAN_ID_PTR localID, int *numEdges,
			  ZOLTAN_ID_PTR nborGID, int *nborProc,
			  int wgtDim, float *ewgts, int *err)
{
    (void) wgtDim; (void) globalID; (void) numEdges; (void) ewgts;
    
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> Mat;
    typedef TransWellGraph<Mat> Graph;

    const Graph& graph = *static_cast<const Graph*>(graphPointer);
    const Mat& transMat = graph.edgeWgt_;
    double scaler = graph.scaler_;

    int idx = 0;
    
    for (auto row = transMat.begin(); row != transMat.end(); ++row)
    {
	auto col = row->begin();
	for (; col != row->end(); ++col)
	{    
	    if (col.index() != row.index()) {
		
		nborGID[idx] = col.index();
		ewgts[idx++] = scaler*(*col);
	    }
	}	
    }
    for ( int i = 0; i < idx; ++i )
    {
        nborProc[i] = 0;
    }
}

// Edge list 4, trans and well weights
void getWellWeightMatEdgeList(void *graphPointer, int sizeGID, int sizeLID,
			      int numCells, ZOLTAN_ID_PTR globalID, 
			      ZOLTAN_ID_PTR localID, int *numEdges,
			      ZOLTAN_ID_PTR nborGID, int *nborProc,
			      int wgtDim, float *ewgts, int *err)
{
    (void) wgtDim; (void) globalID; (void) numEdges; (void) ewgts;
    
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> Mat;
    typedef TransWellGraph<Mat> Graph;
    
    const Graph& graph = *static_cast<const Graph*>(graphPointer);
    const Mat& transMat = graph.edgeWgt_;
    const Mat& wellsMat = graph.wells;
    double scaler = graph.scaler_;

    int idx = 0;
    auto trow = transMat.begin();
    auto wrow = wellsMat.begin();
    for (; trow != transMat.end(); ++trow, ++wrow)
    {	
	auto wcol = wrow->begin();
	for (; wcol != wrow->end(); ++wcol)
	{
	    nborGID[idx] = wcol.index();
	    ewgts[idx++] = 1.0e10;//std::numeric_limits<float>::max();
	}

	auto tcol = trow->begin();
	for (; tcol != trow->end(); ++tcol)
	{    
	    if (tcol.index() != trow.index()) {
		
		nborGID[idx] = tcol.index();
		ewgts[idx++] = scaler*(*tcol);
	    }
	}
    }
    for ( int i = 0; i < idx; ++i )
    {
        nborProc[i] = 0;
    }
}

// Edge List 5, coarse graph
void getWeightMatEdgeListCoarse(void *graphPointer, int sizeGID, int sizeLID,
				int numCells, ZOLTAN_ID_PTR globalID, 
				ZOLTAN_ID_PTR localID, int *numEdges,
				ZOLTAN_ID_PTR nborGID, int *nborProc,
				int wgtDim, float *ewgts, int *err)
{
    (void) wgtDim; (void) globalID; (void) numEdges;
    
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> Mat;
    typedef TransWellGraph<Mat> Graph;
    

    const Graph& graph = *static_cast<const Graph*>(graphPointer);
    const std::vector<std::map<int, double> >& edges = graph.cedges;

    int idx = 0;
    
    for (const auto& node : edges)
    {
	for (const auto& edge : node) {
	    //std::cout << "in getWeightMatEdgeListCoarse idx: " << idx << " "<< edge.first << std::endl;
	    nborGID[idx] = edge.first;
	    ewgts[idx++] = edge.second;
	}
    }
    for ( int i = 0; i < idx; ++i )
    {
        nborProc[i] = 0;
    }
}

template<class G>
void setMatZoltanGraphFunctions(Zoltan_Struct *zz, const G& graph,
				bool pretendNull, bool weights,
				bool wells, bool coarse)
{
    G *graphPointer = const_cast<G*>(&graph);
    if ( pretendNull )
    {
        Zoltan_Set_Num_Obj_Fn(zz, getNullNumCells, graphPointer);
        Zoltan_Set_Obj_List_Fn(zz, getNullVertexList, graphPointer);
        Zoltan_Set_Num_Edges_Multi_Fn(zz, getNullNumEdgesList, graphPointer);
        Zoltan_Set_Edge_List_Multi_Fn(zz, getNullEdgeList, graphPointer);
    }
    else
    {
	if (coarse) {
	    Zoltan_Set_Num_Obj_Fn(zz, getMatNumCellsCoarse, graphPointer);
	    Zoltan_Set_Obj_List_Fn(zz, getMatVertexListCoarse, graphPointer);
	    Zoltan_Set_Num_Edges_Multi_Fn(zz, getMatNumEdgesListCoarse, graphPointer);
	    Zoltan_Set_Edge_List_Multi_Fn(zz, getWeightMatEdgeListCoarse, graphPointer);

	}
	else {
	    Zoltan_Set_Num_Obj_Fn(zz, getMatNumCells, graphPointer);
	    Zoltan_Set_Obj_List_Fn(zz, getMatVertexList, graphPointer);
	    if ( wells )
	    {
		Zoltan_Set_Num_Edges_Multi_Fn(zz, getMatWellNumEdgesList, graphPointer);
	    
		if ( weights )
		    Zoltan_Set_Edge_List_Multi_Fn(zz, getWellWeightMatEdgeList, graphPointer);
		else
		Zoltan_Set_Edge_List_Multi_Fn(zz, getWellMatEdgeList, graphPointer);
	    }
	    else 
	    {
		Zoltan_Set_Num_Edges_Multi_Fn(zz, getMatNumEdgesList, graphPointer);
		if ( weights )
		    Zoltan_Set_Edge_List_Multi_Fn(zz, getWeightMatEdgeList, graphPointer);
		else
		    Zoltan_Set_Edge_List_Multi_Fn(zz, getMatEdgeList, graphPointer);
	    }
	}
    }
}

template<class G>
void setMatZoltanHyperGraphFunctions(Zoltan_Struct *zz, const G& graph, bool pretendNull)
{
    G *graphPointer = const_cast<G*>(&graph);
    if ( pretendNull ) {
	Zoltan_Set_Num_Obj_Fn(zz, getNullNumCells, graphPointer);
	Zoltan_Set_Obj_List_Fn(zz, getNullVertexList, graphPointer);
	Zoltan_Set_HG_Size_CS_Fn(zz, getNullHyperGraphSize, graphPointer);
	Zoltan_Set_HG_CS_Fn(zz, getNullHyperGraphList, graphPointer);
	Zoltan_Set_HG_Size_Edge_Wts_Fn(zz, getNullHyperGraphWgtSize, graphPointer);
	Zoltan_Set_HG_Edge_Wts_Fn(zz, getNullHyperGraphWgtVal, graphPointer);
    }

    else {
	Zoltan_Set_Num_Obj_Fn(zz, getMatNumCellsCoarseHyper, graphPointer);
	Zoltan_Set_Obj_List_Fn(zz, getMatVertexListCoarseHyper, graphPointer);
	Zoltan_Set_HG_Size_CS_Fn(zz, getCpGridHyperGraphSize, graphPointer);
	Zoltan_Set_HG_CS_Fn(zz, getCpGridHyperGraphList, graphPointer);
	Zoltan_Set_HG_Size_Edge_Wts_Fn(zz, getCpGridHyperGraphWgtSize, graphPointer);
	Zoltan_Set_HG_Edge_Wts_Fn(zz, getCpGridHyperGraphWgtVal, graphPointer);
    }
}


template<class Comm, class M, class D>
void zoltanPartitionFunction(std::vector<int>& mpirank, M& g , M& wells, Comm comm, D dr, std::vector<int>& row_size, int numGlobParts=-1, double coarsenGraph=-1)
{
    int rank = comm.rank();

    bool pretendNull = rank!=0;
    
    int wgtType = std::stoi(dr.dict[0]);
    bool useWeights = wgtType != 0;
    bool useWells = std::stoi(dr.dict[4]) == 1;
    int objWgtMet = std::stoi(dr.dict[7]);
    bool useObjWeights = objWgtMet  > 0;
    double logBase = std::exp(std::stod(dr.dict[11]));    
    bool partCoarseGraph = std::stoi(dr.dict[13]) == 1;
    int maxNodeSize = std::stoi(dr.dict[15]);
    double dictCoarsenGraph = std::stod(dr.dict[16]);
    bool useParMetis = std::stoi(dr.dict[17]) == 1;
    bool useHyper = std::stoi(dr.dict[17]) == 2;
    bool useMetis = std::stoi(dr.dict[17]) == 3;
    bool useAMG = std::stoi(dr.dict[17]) == 4;
    int amgLevel = std::stoi(dr.dict[18]);
    
    TransWellGraph<M> twg(g, wells, row_size, wgtType, logBase, rank==0);

    if (!useMetis) {
	int rc = ZOLTAN_OK - 1;
	float ver= 0;
	int argcc = 0;
	char ** argvv = 0;
	struct Zoltan_Struct *zz;
    
	int changes, numGidEntries,numLidEntries,numImport,numExport;
	ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids,exportLocalGids;
	int *importProcs, *importToPart, *exportProcs,*exportToPart;

	//MPI_Init(&argc,&argv);
	rc = Zoltan_Initialize(argcc, argvv, &ver);
	zz = Zoltan_Create(comm);
    
	Zoltan_Set_Param(zz,"DEBUG_LEVEL",dr.dict[5].data());
    
	if (!useHyper)
	    Zoltan_Set_Param(zz,"LB_METHOD","GRAPH");
	else
	    Zoltan_Set_Param(zz,"LB_METHOD","HYPERGRAPH");
    
	if (useParMetis) {
	    Zoltan_Set_Param(zz,"GRAPH_PACKAGE","Parmetis");
	    Zoltan_Set_Param(zz,"PARMETIS_METHOD","PartKway");
	    Zoltan_Set_Param(zz,"PARMETIS_OUTPUT_LEVEL","2");
	    Zoltan_Set_Param(zz,"GRAPH_SYMMETRIZE","TRANSPOSE");
	    Zoltan_Set_Param(zz,"GRAPH_SYM_WEIGHT","MAX");
	    Zoltan_Set_Param(zz,"PARMETIS_COARSE_ALG","1");
	}
    
	Zoltan_Set_Param(zz,"LB_APPROACH","PARTITION");
	Zoltan_Set_Param(zz,"NUM_GID_ENTRIES","1");
	Zoltan_Set_Param(zz,"NUM_LID_ENTRIES","1");
	Zoltan_Set_Param(zz,"RETURN_LISTS","ALL");
	Zoltan_Set_Param(zz,"CHECK_GRAPH","2");
	Zoltan_Set_Param(zz,"EDGE_WEIGHT_DIM","0");    
	Zoltan_Set_Param(zz,"OBJ_WEIGHT_DIM","0");
	Zoltan_Set_Param(zz,"PHG_EDGE_SIZE_THRESHOLD",".35");
	Zoltan_Set_Param(zz,"IMBALANCE_TOL",dr.dict[3].data());
	Zoltan_Set_Param(zz,"PHG_USE_TIMERS","0");
	Zoltan_Set_Param(zz,"PHG_COARSEPARTITION_METHOD", "GREEDY");
	//Zoltan_Set_Param(zz,"PHG_REFINEMENT_QUALITY",dr.dict[0].data());
	//Zoltan_Set_Param(zz,"PHG_COARSENING_LIMIT",dr.dict[5].data());
	//Zoltan_Set_Param(zz,"PHG_COARSENING_NCANDIDATE",dr.dict[6].data());
	//Zoltan_Set_Param(zz,"PHG_COARSENING_METHOD",dr.dict[2].data()); 
	//Zoltan_Set_Param(zz,"PHG_REFINEMENT_LOOP_LIMIT",dr.dict[8].data()); 

	if (numGlobParts != -1) {
	    Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", std::to_string(numGlobParts).c_str());
	    Zoltan_Set_Param(zz, "RETURN_LISTS", "PARTS");
	    Zoltan_Set_Param(zz,"DEBUG_LEVEL","0");
	}
    
    
	if (useWeights || useWells)
	    Zoltan_Set_Param(zz,"EDGE_WEIGHT_DIM","1");
    
	if (objWgtMet == 1)
	    Zoltan_Set_Param(zz,"OBJ_WEIGHT_DIM","1");
	if (objWgtMet == 2)
	    Zoltan_Set_Param(zz,"OBJ_WEIGHT_DIM","2");


	if (coarsenGraph != -1) {
	    if (maxNodeSize == -1)
		twg.coarsenGraph(coarsenGraph);
	    else
		twg.coarsenGraphMaxNodeSize(coarsenGraph, maxNodeSize, rank);
	    if (partCoarseGraph)
		Zoltan_Set_Param(zz,"OBJ_WEIGHT_DIM","1");
	}
	else {
	    if (dictCoarsenGraph != -1) {

		Zoltan_Set_Param(zz,"OBJ_WEIGHT_DIM","1");
		double ctv = 0;
		if (rank == 0)
		    ctv = twg.sortTransFindThreshold(dictCoarsenGraph);
		if (maxNodeSize == -1)
		    twg.coarsenGraph(ctv);
		else
		    twg.coarsenGraphMaxNodeSize(ctv, maxNodeSize, rank);
	    } else {
		partCoarseGraph = false;
	    }
	}
	if (useAMG) {
	    Zoltan_Set_Param(zz,"OBJ_WEIGHT_DIM","1");
	    partCoarseGraph = true;
	    typedef Dune::OwnerOverlapCopyCommunication<int,int> DummyComm;
	    std::shared_ptr<DummyComm> dummyComm(new DummyComm(comm));
	    dummyComm->remoteIndices().template rebuild<false>();
	    twg.createAmgGraph(*dummyComm, amgLevel);
	    if (rank == 0) {
		std::cout << "Created the AMG graph" << std::endl;
	    }
	}
	if (useHyper)
	    setMatZoltanHyperGraphFunctions(zz, twg, pretendNull);
	else
	    setMatZoltanGraphFunctions(zz, twg, pretendNull, useWeights, useWells, partCoarseGraph);

	rc = Zoltan_LB_Partition(zz,
				 &changes, /* 1 if partitioning was changed, 0 otherwise */
				 &numGidEntries,
				 &numLidEntries,
				 &numImport,
				 &importGlobalGids,
				 &importLocalGids,
				 &importProcs,
				 &importToPart,
				 &numExport, /* Number of vertices I must send to other processes*/
				 &exportGlobalGids, /* Global IDs of the vertices I must send */
				 &exportLocalGids, /* Local IDs of the vertices I must send */
				 &exportProcs, /* Process to which I send each of the vertices */
				 &exportToPart); /* Partition to which each vertex will belong */
  
	if (rc!=ZOLTAN_OK)
	    std::cout << "Error occured" << std::endl;

	if (partCoarseGraph) {

	    if (rank == 0) {
		std::vector<int> coarsePartRes(twg.courseNodes_.size());
		for (int i = 0; i < numExport; ++i) {
		    coarsePartRes[exportLocalGids[i]] = exportProcs[i];
		}

		for (int i = 0; i < mpirank.size(); ++i) {
	
		    mpirank[i] = coarsePartRes[twg.f2c[i]];
		}
	    }
	}  else {
	    for (int i = 0; i < numExport; ++i){
		mpirank[exportLocalGids[i]] = exportProcs[i];
	    }
	}

	if (numGlobParts != -1) {

	    if (partCoarseGraph) {

		for (int i = 0; i < mpirank.size(); ++i) {
		    mpirank[i] = exportToPart[twg.f2c[i]];
		}

	    } else {
		for (int i = 0; i < mpirank.size(); ++i) {
		    mpirank[i] = exportToPart[i];
		}
	    }

	} else {

	    std::vector<int> rankIsZero(comm.size(), 0);

	    for (int i = 0; i < mpirank.size(); ++i) {
	
		rankIsZero[mpirank[i]] +=1;
	    }

	    bool zeroRankPresent = false;
	    for (int r = 0; r < rankIsZero.size(); ++r) {
		if (rankIsZero[r] == 0) {
		    zeroRankPresent = true;
		}
	    }
    

	    if (rank ==0) {

		for (int r = 0; r < rankIsZero.size(); ++r) {
		    std::cout << r << ":" << rankIsZero[r]<<std::endl;
		}
		std::cout <<std::endl;

		if (zeroRankPresent) {
		    std::cout << "Zero partitions present: ";
		    for (int r = 0; r < rankIsZero.size(); ++r) {
			if (rankIsZero[r] == 0)
			    std::cout << r << " ";
		    }
		    std::cout <<std::endl;
		}
	    }
	}
	Zoltan_Destroy(&zz);  
    }
    else {
	if (coarsenGraph != -1) {
	    if (maxNodeSize == -1)
		twg.coarsenGraph(coarsenGraph);
	    else
		twg.coarsenGraphMaxNodeSize(coarsenGraph, maxNodeSize, rank);
	}
	else {
	    if (dictCoarsenGraph != -1) {
		double ctv = 0;
		if (rank == 0)
		    ctv = twg.sortTransFindThreshold(dictCoarsenGraph);
		if (maxNodeSize == -1)
		    twg.coarsenGraph(ctv);
		else
		    twg.coarsenGraphMaxNodeSize(ctv, maxNodeSize, rank);
	    } else {
		partCoarseGraph = false;
	    }
	}
	std::vector<int> partC;
	if (rank == 0) {
	    partWithMetis(twg, dr, partC, comm.size());

	    for (int i = 0; i < mpirank.size(); ++i) {
	
		mpirank[i] = partC[twg.f2c[i]];
	    }
	}
    }	  
    comm.broadcast(&mpirank[0], mpirank.size(), 0);
}

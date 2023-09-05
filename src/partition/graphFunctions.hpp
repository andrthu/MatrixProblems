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

template <class Mat>
class TransWellGraph
{
public:

    TransWellGraph(Mat T, Mat W, std::vector<int> rs, int useNormal, double logBase, bool isRoot)
    {
	if (isRoot) {
	    trans = T;
	    wells = W;
	    row_size = rs;
	    base_ = std::exp(logBase);
	    if (useNormal < 2)
	    {
		edgeWgt_ = T;
		scaler_ = 1.0e18;
		sortTrans();
	    }
	    else if (useNormal == 2)
	    {
		findMin();
		createLogWeights();
		scaler_ = 1.0;
	    }
	    else if (useNormal == 4)
	    {
		sortTrans();
		createCatWeights();
		scaler_ = 1.0;
	    }
	}
    }

    void findMin()
    {
	double minVal = 1.0e18;
	for (auto trow = trans.begin(); trow != trans.end(); trow++)
	{	    		    	
	    auto tcol = trow->begin();
	    for (; tcol != trow->end(); ++tcol)
	    {
		if (*tcol != 0.0)
		{
		    if (minVal > *tcol)
		    {
			minVal = *tcol;
		    }
		}
	    }
	}
	minLogWgt = std::log(minVal);
	//std::cout << minLogWgt << " "<< minVal << std::endl; 
    }
    
    void createLogWeights()
    {
	Mat logWgt(trans);
	for (auto trow = trans.begin(); trow != trans.end(); trow++)
	{	    		    	
	    auto tcol = trow->begin();
	    int rid = trow.index();
	    for (; tcol != trow->end(); ++tcol)
	    {
		int cid = tcol.index();
		if (*tcol != 0.0)
		{
		    logWgt[rid][cid] = 1.0 + (std::log(*tcol) - minLogWgt)/std::log(base_);
		}
		else
		{
		    logWgt[rid][cid] = 0.0;
		}
		//if (logWgt[rid][cid] < 0)
		//std::cout << rid << " " << cid<< " "<< logWgt[rid][cid] << " " << *tcol<<" "<< trans[rid][cid]<<std::endl; 
	    }
	}	
	edgeWgt_ = logWgt;
    }

    void sortTrans()
    {
	unsigned nnz = trans.nonzeroes();
	std::vector<double>trans_list(nnz, 0.0);

	int idx = 0;
	for (auto row = trans.begin(); row != trans.end(); ++row)
	{
	    auto col = row->begin();
	    for (; col!=row->end(); ++col)
	    {
		trans_list[idx] = *col;	
	    }
	}
	
	std::sort(trans_list.begin(), trans_list.end());

	trans_bound_.push_back(trans_list[nnz/4]);
	trans_bound_.push_back(trans_list[3*(nnz/4)]);

    }

    void createCatWeights()
    {
	Mat catWgt(trans);
	for (auto trow = trans.begin(); trow != trans.end(); trow++)
	{	    		    	
	    auto tcol = trow->begin();
	    int rid = trow.index();
	    for (; tcol != trow->end(); ++tcol)
	    {
		int cid = tcol.index();
		double t = *tcol;
		if (t == 0.0)
		{
		    catWgt[rid][cid] = 0.1;
		}
		else if (t < trans_bound_[0])
		{
		    catWgt[rid][cid] = 1.0;
		}
		else if ( t < trans_bound_[1] && t > trans_bound_[0])
		{
		    catWgt[rid][cid] = 10.0;
		}
		else
		{
		    catWgt[rid][cid] = 100.0;
		}
	    }
	}	
	edgeWgt_ = catWgt;
    }

    Mat trans;
    Mat wells;
    Mat edgeWgt_;
    double minLogWgt;
    double scaler_;
    double base_;

    std::vector<double> trans_bound_;
    std::vector<int> row_size;
};

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

template<class G>
void setMatZoltanGraphFunctions(Zoltan_Struct *zz, const G& graph,
				bool pretendNull, bool weights, bool wells)
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


template<class Comm, class M, class D>
void zoltanPartitionFunction(std::vector<int>& mpirank, M& g , M& wells, Comm comm, D dr, std::vector<int>& row_size)
{
    int rank = comm.rank();
    
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
    Zoltan_Set_Param(zz,"LB_METHOD","GRAPH");

    //Zoltan_Set_Param(zz,"GRAPH_PACKAGE","Parmetis");
    //Zoltan_Set_Param(zz,"PARMETIS_METHOD","PartKway");
    //Zoltan_Set_Param(zz,"PARMETIS_OUTPUT_LEVEL","2");

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
    //Zoltan_Set_Param(zz,"PHG_REFINEMENT_QUALITY",dr.dict[0].data());
    Zoltan_Set_Param(zz,"PHG_COARSEPARTITION_METHOD", "GREEDY");
    //Zoltan_Set_Param(zz,"PHG_COARSENING_LIMIT",dr.dict[5].data());
    //Zoltan_Set_Param(zz,"PHG_COARSENING_NCANDIDATE",dr.dict[6].data());
    //Zoltan_Set_Param(zz,"PHG_COARSENING_METHOD",dr.dict[2].data()); 
    //Zoltan_Set_Param(zz,"PHG_REFINEMENT_LOOP_LIMIT",dr.dict[8].data()); 
  
    bool pretendNull = rank!=0;
    
    int wgtType = std::stoi(dr.dict[0]);
    bool useWeights = wgtType != 0;
    bool useWells = std::stoi(dr.dict[4]) == 1;
    bool useObjWeights = std::stoi(dr.dict[7]) == 1;
    double logBase = std::exp(std::stod(dr.dict[11]));    
    
    if (useWeights || useWells)
	Zoltan_Set_Param(zz,"EDGE_WEIGHT_DIM","1");
    
    if (useObjWeights)
	Zoltan_Set_Param(zz,"OBJ_WEIGHT_DIM","1");
    
    TransWellGraph<M> twg(g, wells, row_size, wgtType, logBase, rank==0);
    //twg.findMin();
    //twg.createLogWeights();

    setMatZoltanGraphFunctions(zz, twg, pretendNull, useWeights, useWells);

    rc = Zoltan_LB_Partition(zz,
			     &changes,
			     &numGidEntries,
			     &numLidEntries,
			     &numImport,
			     &importGlobalGids,
			     &importLocalGids,
			     &importProcs,
			     &importToPart,
			     &numExport,
			     &exportGlobalGids,
			     &exportLocalGids,
			     &exportProcs,
			     &exportToPart);
  
    if (rc!=ZOLTAN_OK)
	std::cout << "Error occured" << std::endl;

    for (int i = 0; i < numExport; ++i)
    {
	mpirank[exportLocalGids[i]] = exportProcs[i];
    }

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
    
    comm.broadcast(&mpirank[0], mpirank.size(), 0);
    
    Zoltan_Destroy(&zz);  
}

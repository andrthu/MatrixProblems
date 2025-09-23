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

template <class Mat>
class TransWellGraph
{
public:

    TransWellGraph(Mat T, Mat W, std::vector<int> rs, int useNormal, double logBase, bool isRoot)
    {
	if (isRoot) {
	    wgtType_ = useNormal;
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

    template<class R>
    void dps(R row, int v, int master, double w, std::vector<bool>& visited,
	     std::vector<int>& f2c, std::vector<int>& cnode, std::vector<std::tuple<int,int,double> >& edges) {

	visited[v] = true;
	f2c[v] = master;
	cnode.push_back(v);
	
	auto col = row.begin();
	for (; col != row.end(); ++col) {
	    int nab = col.index();

	    if (trans[v][nab] > w) {
		if (!visited[nab]) {
		    dps(trans[nab],nab,master,w,visited,f2c,cnode,edges);
		} else {
		    if (f2c[v]!=f2c[nab]) {
			std::cout << "Problem " << nab << " " << v <<
			    " " << f2c[v] << " " << f2c[nab] <<std::endl; 
		    }
		}
	    }
	}
	col = row.begin();
	for (; col != row.end(); ++col) {
	    int nab = col.index();
	    if (f2c[v]!=f2c[nab]) {
		edges.push_back({v,nab,trans[v][nab]});
	    }
	}
    }

    void createCoarseEdges(std::vector<std::vector<std::tuple<int,int,double> >> gEdges,
			   std::vector<int> f2c) {

	//std::vector<std::map<int, double> > cedges;

	for (std::vector<std::tuple<int,int,double> > es : gEdges ) {

	    std::map<int, double> ce;

	    for (std::tuple<int,int,double> fe : es) {

		int coarseNab = f2c[std::get<1>(fe)];
		double weight = wgtType_ == 0 ? 1.0 : std::get<2>(fe);
		if ( ce.count(coarseNab) == 1 ) {
		    ce[coarseNab] += weight;
		} else {
		    ce.insert({coarseNab,weight});
		}
	    }
	    cedges.push_back(ce);
	}
    }
    
    void coarsenGraph(double w) {

	int N = trans.N();
	std::vector<bool> visited(N, false);

	f2c.resize(N, 0);
	std::vector<int> c2f;

	int biggest = 0;
	int single = 0;

	std::vector<std::vector<std::tuple<int,int,double> >> gEdges;
	int newV = 0;
	for (int v = 0; v < N; ++v) {

	    if (!visited[v]) {

		//f2c[v] = v;
		c2f.push_back(v);
		std::vector<int> cnode;
		std::vector<std::tuple<int,int,double> > edges;		
		dps(trans[v],v,newV,w,visited,f2c,cnode,edges);
		newV++;
		if (cnode.size() > biggest)
		    biggest = cnode.size();
		if (cnode.size() == 1)
		    single++;
		gEdges.push_back(edges);
		courseNodes_.push_back(cnode);
	    }
	}
	createCoarseEdges(gEdges,f2c);
	std::cout << "Coarse graph size: " << c2f.size()<< " "<< biggest<< " "<< single<< " " << N << std::endl;
    }

    Mat trans;
    Mat wells;
    Mat edgeWgt_;
    double minLogWgt;
    double scaler_;
    double base_;
    int wgtType_;

    
    std::vector<std::vector<int>> courseNodes_;
    std::vector<std::map<int, double> > cedges;
    std::vector<int> f2c;

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


template<class Comm, class M, class D>
void zoltanPartitionFunction(std::vector<int>& mpirank, M& g , M& wells, Comm comm, D dr, std::vector<int>& row_size, int numGlobParts=-1, double coarsenGraph=-1)
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

    if (numGlobParts != -1) {
	Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", std::to_string(numGlobParts).c_str());
	Zoltan_Set_Param(zz, "RETURN_LISTS", "PARTS");
	Zoltan_Set_Param(zz,"DEBUG_LEVEL","0");
    }
    
    bool pretendNull = rank!=0;
    
    int wgtType = std::stoi(dr.dict[0]);
    bool useWeights = wgtType != 0;
    bool useWells = std::stoi(dr.dict[4]) == 1;
    int objWgtMet = std::stoi(dr.dict[7]);
    bool useObjWeights = objWgtMet  > 0;
    double logBase = std::exp(std::stod(dr.dict[11]));    
    bool partCoarseGraph = std::stoi(dr.dict[13]) == 1;
    
    if (useWeights || useWells)
	Zoltan_Set_Param(zz,"EDGE_WEIGHT_DIM","1");
    
    if (objWgtMet == 1)
	Zoltan_Set_Param(zz,"OBJ_WEIGHT_DIM","1");
    if (objWgtMet == 2)
	Zoltan_Set_Param(zz,"OBJ_WEIGHT_DIM","2");
    
    TransWellGraph<M> twg(g, wells, row_size, wgtType, logBase, rank==0);
    //twg.findMin();
    //twg.createLogWeights();

    if (coarsenGraph != -1) {
	twg.coarsenGraph(coarsenGraph);
	if (partCoarseGraph)
	    Zoltan_Set_Param(zz,"OBJ_WEIGHT_DIM","1");
    }
    else
	partCoarseGraph = false;
    
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

    
    for (int i = 0; i < numExport; ++i)
    {
	mpirank[exportLocalGids[i]] = exportProcs[i];
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
    
    comm.broadcast(&mpirank[0], mpirank.size(), 0);
    
    Zoltan_Destroy(&zz);  
}

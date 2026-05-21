/*
  Copyright 2026 Andreas Thune.

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

 
int getMatNumCellsCoarseHyper(void* graphPointer, int* err)
{
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> Mat;
    typedef TransWellGraph<Mat> Graph;

    const Graph& graph = *static_cast<const Graph*>(graphPointer);
    const std::vector<std::vector<int>>& nodes = graph.courseNodes_;

    *err = ZOLTAN_OK;
    return nodes.size();
}

void getMatVertexListCoarseHyper(void* graphPointer, int numGlobalIdEntries,
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
        gids[idx] = idx;
        lids[idx] = idx;
	if (wgtDim == 1)
	    objWgts[idx] = nodes[idx].size();
    }

    *err = ZOLTAN_OK;
}

void getCpGridHyperGraphSize(void *graphPointer, int* num_lists,
			     int *num_pins, int *format, int *err)
{
    *format = ZOLTAN_COMPRESSED_EDGE;
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> Mat;
    typedef TransWellGraph<Mat> Graph;
    const Graph& graph = *static_cast<const Graph*>(graphPointer);

    const std::vector<std::map<int, double> >& edges = graph.cedges;

    *num_lists = edges.size();

    int numEdge = 0;
    for (const auto& neighbors : edges) {
        numEdge += neighbors.size();
    }
    *num_pins = numEdge;
    *err = ZOLTAN_OK;
}

void getCpGridHyperGraphList(void *graphPointer, int num_gid_entries, 
			     int num_vtx_edge, int num_pins, int format,
			     ZOLTAN_ID_PTR vtxedge_GID, int *vtxedge_ptr, 
			     ZOLTAN_ID_PTR pin_GID, int *err)
{
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> Mat;
    typedef TransWellGraph<Mat> Graph;
    const Graph& graph = *static_cast<const Graph*>(graphPointer);
    const std::vector<std::map<int, double> >& edges = graph.cedges;

    int hyperEdgeIdx = 0;
    int edgeIdx = 0;


    for( int i = 0; i < edges.size();  i++ ) {

        vtxedge_GID[hyperEdgeIdx] = hyperEdgeIdx;
	vtxedge_ptr[hyperEdgeIdx] = edgeIdx;
	hyperEdgeIdx++;

	for (const auto& edge : edges[i]) {
	    pin_GID[edgeIdx] = edge.first;
	    edgeIdx++;
	}
    }
    *err = ZOLTAN_OK;
}

void getCpGridHyperGraphWgtSize(void *graphPointer, int *num_edges, int *err)
{
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> Mat;
    typedef TransWellGraph<Mat> Graph;
    const Graph& graph = *static_cast<const Graph*>(graphPointer);
    const std::vector<std::vector<int>>& nodes = graph.courseNodes_;

    *num_edges = nodes.size();
    *err = ZOLTAN_OK;
}

void getCpGridHyperGraphWgtVal(void *graphPointer, int num_gid_entries, 
			       int num_lid_entries, int num_edges, 
			       int edge_weight_dim, ZOLTAN_ID_PTR edge_GID, 
			       ZOLTAN_ID_PTR edge_LID, float  *edge_weight, int *err)
{
    typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1>> Mat;
    typedef TransWellGraph<Mat> Graph;
    const Graph& graph = *static_cast<const Graph*>(graphPointer);
    const std::vector<std::vector<int>>& nodes = graph.courseNodes_;


    int edgeIdx = 0;
    for( int i = 0; i < nodes.size();  i++ ) {

        edge_GID[edgeIdx] = edgeIdx;
	edge_LID[edgeIdx] = edgeIdx;
	edge_weight[edgeIdx] = nodes[i].size();
	edgeIdx++;
    }
    *err = ZOLTAN_OK;
    
}

void getNullHyperGraphSize(void *graphPointer, int* num_lists,
			   int *num_pins, int *format, int *err)
{
    (void) graphPointer;
    *num_lists = 0;
    *num_pins = 0;
    *format = ZOLTAN_COMPRESSED_EDGE;

    *err = ZOLTAN_OK;
}

void getNullHyperGraphList(void *graphPointer, int num_gid_entries, 
			   int num_vtx_edge, int num_pins, int format,
			   ZOLTAN_ID_PTR vtxedge_GID, int *vtxedge_ptr, 
			   ZOLTAN_ID_PTR pin_GID, int *err)
{
    (void) graphPointer; (void) num_gid_entries; (void) num_vtx_edge; 
    (void) num_pins; (void) format; (void) vtxedge_GID;
    (void) vtxedge_ptr; (void) pin_GID;
    
    *err = ZOLTAN_OK;
}

void getNullHyperGraphWgtSize(void *graphPointer, int *num_edges, int *err)
{
    (void) graphPointer;
    *num_edges = 0;

    *err = ZOLTAN_OK;
}

void getNullHyperGraphWgtVal(void *graphPointer, int num_gid_entries, 
			     int num_lid_entries, int num_edges, 
			     int edge_weight_dim, ZOLTAN_ID_PTR edge_GID, 
			     ZOLTAN_ID_PTR edge_LID, float  *edge_weight, int *err)
{
    (void) graphPointer; (void) num_gid_entries; (void) num_lid_entries; 
    (void) num_edges; (void) edge_weight_dim; (void) edge_GID;
    (void) edge_LID; (void) edge_weight;

    *err= ZOLTAN_OK;
}


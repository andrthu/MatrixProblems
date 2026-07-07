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

#include <metis.h>

template<class Graph, class D>
void partWithMetis(const Graph& graph, D dr, std::vector<int>& mpirank, int numPart)
{

    // 1. Prepare CSR Arrays
    std::vector<idx_t> xadj;
    std::vector<idx_t> adjncy;
    std::vector<idx_t> vwgt;
    std::vector<idx_t> adjwgt;

    const std::vector<std::map<int, double> >& edges = graph.cedges;
    const std::vector<std::vector<int>>& nodes = graph.courseNodes_;

    int n_nodes = nodes.size();
    mpirank.resize(n_nodes,0);
    std::vector<idx_t> dummy(n_nodes,0);
    xadj.push_back(0);
    for (int i = 0; i < n_nodes; ++i) {
        // Add vertex weight
        vwgt.push_back(nodes[i].size());

	for (const auto& edge : edges[i]) {

	    adjncy.push_back(edge.first);
	    adjwgt.push_back(edge.second);
	}
	// Mark the end of this node's list in adjncy
	xadj.push_back(adjncy.size());
    }

    // 2. METIS Parameters
    idx_t nvtxs = n_nodes;     // Number of vertices
    idx_t ncon = 1;            // Number of balancing constraints
    idx_t nparts = numPart;    // Number of partitions you want
    idx_t objval;              // Output: edge-cut value

    std::vector<real_t> tpwgts(nparts * ncon);
    for (int i = 0; i < nparts; ++i) {
        tpwgts[i] = 1.0 / nparts;
    }
    std::vector<real_t> ubvec(ncon);
    ubvec[0] = std::stod(dr.dict[3].data());

    // METIS options (0 sets defaults)
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    options[METIS_OPTION_NUMBERING] = 0; // C-style 0-based indexing

    // 3. Call METIS
    int result = METIS_PartGraphKway(&nvtxs, &ncon, xadj.data(), adjncy.data(), 
                                     vwgt.data(), NULL, adjwgt.data(), &nparts, tpwgts.data(),
				     ubvec.data(), options, &objval, dummy.data());

    for (int i=0; i<n_nodes;++i)
        mpirank[i] = dummy[i]; 
    if (result == METIS_OK) {
        std::cout << "Partitioning successful. Edge cut: " << objval << std::endl;
	//for (int i = 0; i < nvtxs; ++i) {
	//    std::cout << "Node " << i << " is in partition " << mpirank[i] << std::endl;
        //}
    }
}

/*
  Copyright 2022 Andreas Thune.

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

#include "config.h"
#include <array>
#include <vector>
#include <cstdlib>
#include <memory>
#include <cmath>
#include <tuple>
#include <limits>
#include <chrono>
#include <unordered_set>

#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/parallel/indexset.hh>
#include <dune/common/timer.hh>
#include <dune/common/shared_ptr.hh>

#include <dune/istl/matrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/matrixindexset.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/repartition.hh>
#include <dune/istl/matrixredistribute.hh>
#include <dune/istl/schwarz.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <dune/common/function.hh>
#include <boost/property_tree/json_parser.hpp>

#include <opm/simulators/linalg/ParallelOverlappingILU0.hpp>
#include <opm/simulators/linalg/PreconditionerFactory.hpp>
#include <opm/simulators/linalg/WellOperators.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/amgcpr.hh>
#include <opm/simulators/linalg/PressureSolverPolicy.hpp>
#include <opm/simulators/linalg/PressureTransferPolicy.hpp>
#include <opm/simulators/linalg/OwningTwoLevelPreconditioner.hpp>
#include <opm/simulators/linalg/getQuasiImpesWeights.hpp>
//#include <opm/simulators/linalg/matrixblock.hh>
//#include <opm/simulators/linalg/MatrixBlock.hpp>


#include <opm/io/eclipse/EGrid.hpp>
#include <opm/io/eclipse/ERft.hpp>
#include <opm/io/eclipse/ERst.hpp>
#include <opm/io/eclipse/ESmry.hpp>
#include <opm/io/eclipse/ERsm.hpp>
#include <opm/io/eclipse/RestartFileView.hpp>
#include <opm/io/eclipse/rst/state.hpp>

#include <opm/output/eclipse/RestartIO.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/utility/numeric/cmp.hpp>

#include <mpi.h>
#include <zoltan.h>

template<class Vec>
void simple_vec_diffs(const Vec& a1, const Vec& a2 )
{

    typedef typename Vec::value_type X;
    
    int n = a1.size();

    X sum = 0;
    X maxE = 0;
    X relMax = 0;
    X relSum = 0;
    
    for (int i = 0; i < n; ++i) {

	auto e = std::abs(a1[i]-a2[i]);
	auto rel_e = e;
	if (std::max(a1[i], a2[i]) > 0)
	    rel_e = e/std::max(a1[i], a2[i]);
	sum += e;
	relSum += rel_e;
	
	if (e > maxE)
	    maxE = e;
	if (rel_e > relMax)
	    relMax = rel_e;
    }

    std::cout << "Max and average error " << maxE << " " << sum/n << std::endl;
    std::cout << "Relative max and average error " << relMax << " " << relSum/n << std::endl;
    std::cout << std::endl;
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    std::string rootName1 = std::string(argv[1]);
    std::string rootName2 = std::string(argv[2]);
    
    Opm::EclIO::EGrid* grid1 = new Opm::EclIO::EGrid(rootName1+".EGRID");
    auto rst1 = std::make_shared<Opm::EclIO::ERst>(rootName1+".UNRST");

    Opm::EclIO::EGrid* grid2 = new Opm::EclIO::EGrid(rootName2+".EGRID");
    auto rst2 = std::make_shared<Opm::EclIO::ERst>(rootName2+".UNRST");


    std::cout << "Active cells: " << grid1->activeCells() << std::endl;
    std::vector<int> seqnums1 = rst1->listOfReportStepNumbers();

    for (int i=0; i < seqnums1.size(); ++i) { 
	int seqn = seqnums1[i];
	std::cout << seqnums1[i] << " "<< seqnums1.size()<<std::endl;
	rst1->loadReportStepNumber(seqn);

	auto arrays = rst1->listOfRstArrays(seqn);

	auto P1 = rst1->getRestartData<float>("PRESSURE", seqn, 0);
	//std::cout << "Pressure vec size "<<P1.size()<< std::endl;
	auto P2 = rst2->getRestartData<float>("PRESSURE", seqn, 0);

	std::cout << "PRESSURE diff" <<std::endl;
	simple_vec_diffs(P1,P2);

	std::cout << "SGAS diff" <<std::endl;
	auto G1 = rst1->getRestartData<float>("SGAS", seqn, 0);
	auto G2 = rst2->getRestartData<float>("SGAS", seqn, 0);
	simple_vec_diffs(G1,G2);

	std::cout << "SWAT diff" <<std::endl;
	auto W1 = rst1->getRestartData<float>("SWAT", seqn, 0);
	auto W2 = rst2->getRestartData<float>("SWAT", seqn, 0);
	simple_vec_diffs(W1,W2);
    }
    return 0;
    //./src/eclipseIO ../../../fork_opm/test/caseNew/epycMilan/newCase/amg/cpr-gl/trans/tasks8nodes1/run1/FLOW ../../../fork_opm/test/caseNew/epycMilan/newCase/amg/cpr-gl/trans/tasks16nodes1/run1/FLOW
}

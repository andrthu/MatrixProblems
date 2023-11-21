/*
  Copyright 2023 Andreas Thune.

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
#include <dune/istl/umfpack.hh>
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

#include <mpi.h>
#include <zoltan.h>

#include "io/inputRead.hpp"
#include "solveParallel/ghostLastOperations.hpp"
#include "partition/graphFunctions.hpp"
#include "solveParallel/buildLocalMatrix.hpp"
#include "partition/overlapCreation.hpp"
#include "solveParallel/reorder.hpp"
#include "solveParallel/nonQuadraticLocalMatrix.hpp"
#include "partition/setIndexSet.hpp"
#include "partition/evaluatePartition.hpp"
#include "io/dictRead.hpp"
#include "partition/meassureCommunication.hpp"
#include "partition/meassureOperationsMinLoop.hpp"
#include "partition/timeMultipleTimes.hpp"
#include "partition/flexibleSolverTimer.hpp"
#include "solveParallel/amgSetup.hpp"
#include "generalDim/amgCpr.hpp"
#include "generalDim/multJsonSolveGeneral.hpp"

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    if (help_detect_mult(argc, argv)) {
	return 0;
    }
    auto systemDirs = parse_multiple_systems(argc, argv);
    
    typedef Dune::FieldMatrix<double,3,3> BlockMat3;
    typedef Dune::BCRSMatrix<BlockMat3> Mat;
    typedef Dune::FieldVector<double,3> BlockVec;
    typedef Dune::BlockVector<BlockVec> Vec;
    gen_dim_jsonSolve_mult_sys<Mat,Vec>(systemDirs);

    return 0;
}

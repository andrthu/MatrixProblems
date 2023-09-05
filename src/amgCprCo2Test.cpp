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
//#include <opm/simulators/linalg/PreconditionerFactory.hpp>
#include <opm/simulators/linalg/WellOperators.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/amgcpr.hh>
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

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    typedef Dune::FieldMatrix<double,2,2> BlockMat2;
    typedef Dune::BCRSMatrix<BlockMat2> Mat;
    typedef Dune::FieldVector<double,2> BlockVec;
    typedef Dune::BlockVector<BlockVec> Vec;
    typedef Dune::MPIHelper::MPICommunicator MPICommunicator;
    typedef Dune::CollectiveCommunication<MPICommunicator> CollectiveCommunication;
    typedef Dune::BiCGSTABSolver<Vec> Solver;
    typedef Dune::InverseOperatorResult Stat;

    typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
    typedef Dune::OverlappingSchwarzScalarProduct<Vec,Comm> ScalarProduct; 
    typedef Dune::OverlappingSchwarzOperator<Mat,Vec,Vec,Comm> Operator;
    typedef GhostLastMatrixAdapter<Mat,Vec,Vec,Comm> GLO;
    typedef Opm::ParallelOverlappingILU0<Mat,Vec,Vec,Comm> ILU;
    typedef Dune::Amg::AMGCPR<GLO, Vec, ILU, Comm> AMGCPR;
    typedef Dune::FlexibleSolver<Mat, Vec> FlexibleSolverType;

    using Smoother = ILU;
    using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
    using CriterionBase
	= Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricDependency<Mat, Dune::Amg::FirstDiagonal>>;
    using Criterion = Dune::Amg::CoarsenCriterion<CriterionBase>;
    
    typedef Dune::SeqILU<Mat,Vec,Vec> DuneSmoother;
    typedef Dune::BlockPreconditioner<Vec, Vec, Comm, DuneSmoother> ParSmoother;
    typedef typename Dune::Amg::SmootherTraits<ParSmoother>::Arguments Smoother1Args;
    typedef Dune::Amg::AMGCPR<Operator, Vec, ParSmoother, Comm> AMGCPR_DUNE;

    // ------------------------------------------------ 
    // --- start reading input 
    CollectiveCommunication cc(MPI_COMM_WORLD);
    int rank = cc.rank();

    Mat A_loc;
    Vec rhs_loc;

    DictRead DR;
    Comm comm(cc);
    std::shared_ptr<Comm> parComm(new(Comm));
    readMatOnRootAndDist(argc, argv, A_loc, rhs_loc, DR, comm, parComm, cc);
    // --- Complete reading input

    // ------------------------------------------------
    // --- Set up operators, preconditioners and solvers 
    ScalarProduct sp(*parComm);
    Operator linOp(A_loc, *parComm);
    GLO glLinOp(A_loc, *parComm);
    
    Opm::FlowLinearSolverParameters flsp_amg;
    Opm::PropertyTree prm_amg = setupAMG(std::string("amg"), flsp_amg);
    prm_amg.put("preconditioner.verbosity", 10);

    Criterion criterion(15, prm_amg.get<int>("coarsenTarget", 1200));
    setCrit(criterion, prm_amg);

    SmootherArgs smootherArgs;
    setOpmILU0args(smootherArgs, prm_amg);

    Smoother1Args smoother1Args;
    smoother1Args.iterations = prm_amg.get<int>("iterations", 1);
    
    auto amg = std::make_shared<AMGCPR>(glLinOp, criterion, smootherArgs, *parComm);
    // --- Complete set up
    
    int verb=0;
    if (rank==0) {verb=1;}
    Solver solver(linOp, sp, *amg, 0.01, 30, verb);

    Vec x(rhs_loc.size());
    Vec crhs1(rhs_loc);
    Vec crhs2(rhs_loc);
    x=0;
    Dune::InverseOperatorResult stat;
    solver.apply(x, crhs1,stat);

    auto amg_dune = std::make_shared<AMGCPR_DUNE>(linOp, criterion, smoother1Args, *parComm);

    Solver solver2(linOp, sp, *amg_dune, 0.01, 30, verb);
    solverAndPreconTimer(solver, *amg, rhs_loc, cc, std::string("AMG"), 0.01);
    solverAndPreconTimer(solver2, *amg_dune, rhs_loc, cc, std::string("AMG_DUNE"), 0.01);

    amgHiaInfo(*amg, rhs_loc, cc); // in partition/flexibleSolverTimer.hpp
    return 0;
}

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

#include <mpi.h>
#include <zoltan.h>
#include <boost/property_tree/json_parser.hpp>

//#include <opm/autodiff/ParallelOverlappingILU0.hpp>
#include <opm/simulators/linalg/ParallelOverlappingILU0.hpp>
#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/getQuasiImpesWeights.hpp>
//#include <opm/simulators/linalg/matrixblock.hh>
//#include <opm/simulators/linalg/MatrixBlock.hpp>

#include "io/inputRead.hpp"
#include "solveParallel/ghostLastOperations.hpp"
#include "partition/graphFunctions.hpp"
#include "partition/overlapCreation.hpp"
#include "solveParallel/reorder.hpp"
#include "solveParallel/buildLocalMatrix.hpp"
#include "solveParallel/nonQuadraticLocalMatrix.hpp"
#include "partition/setIndexSet.hpp"
#include "partition/evaluatePartition.hpp"
#include "io/dictRead.hpp"
#include "partition/meassureCommunication.hpp"
#include "partition/meassureOperationsMinLoop.hpp"
#include "partition/timeMultipleTimes.hpp"
#include "partition/flexibleSolverTimer.hpp"

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    typedef Dune::FieldMatrix<double,3,3> BlockMat3;
    typedef Dune::FieldMatrix<double,1,1> BlockMat1;
    typedef Dune::BCRSMatrix<BlockMat3> Mat;
    typedef Dune::BCRSMatrix<BlockMat1> Mat1;
    typedef Dune::FieldVector<double,3> BlockVec;
    typedef Dune::BlockVector<BlockVec> Vec;
    typedef Dune::MPIHelper::MPICommunicator MPICommunicator;
    typedef Dune::CollectiveCommunication<MPICommunicator> CollectiveCommunication;
    typedef Dune::OwnerOverlapCopyAttributeSet::AttributeSet AttributeSet;
    typedef Dune::ParallelIndexSet<int,Dune::ParallelLocalIndex<AttributeSet>> PIS;
    typedef Dune::RemoteIndices<PIS> RIS;
    typedef Dune::BiCGSTABSolver<Vec> Solver;
    typedef Dune::InverseOperatorResult Stat;

    typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
    typedef Dune::OverlappingSchwarzScalarProduct<Vec,Comm> ScalarProduct; 
    typedef Dune::OverlappingSchwarzOperator<Mat,Vec,Vec,Comm> Operator;
    typedef Opm::ParallelOverlappingILU0<Mat,Vec,Vec,Comm> ILU;
    typedef Dune::FlexibleSolver<Mat, Vec> FlexibleSolverType;

    Dune::FMatrixPrecision<double>::set_absolute_limit(1.e-100);
    
    CollectiveCommunication cc(MPI_COMM_WORLD);
    int rank = cc.rank();

    Mat A;
    Mat A_loc;
    Mat1 trans, wells;
    Vec rhs, rhs_loc;

    DictRead DR;

    // Read matrices
    handleMatrixSystemInputSomeRanks(argc, argv, A, trans, wells, rhs, DR, cc, rank==0);
    DR.dict[5] = "2";

    int N;
    if (rank == 0)
	N = A.N();
    cc.broadcast(&N, 1, 0);
    
    std::vector<int> row_size(N);
    storeRowSizeFromRoot(A, row_size, cc);
    
    //partition matrix
    std::vector<int> mpivec(N, rank);
    zoltanPartitionFunction(mpivec, trans, wells, cc, DR, row_size);
    if (rank == 0) {std::cout << "After zoltan "<<std::endl;}
    
    Comm comm(cc);
    std::shared_ptr<Comm> parComm(new(Comm));
    Mat A_loc_;

    //Distribute Matrix and rhs vector from root. A_loc_ incudes off-diagonal NNZ on ghost rows.
    constructLocalFromRoot(A, A_loc_, N, rhs, rhs_loc, mpivec, comm, parComm, cc);
    //constructLocalFromRoot(A, A_loc, N, rhs, rhs_loc, mpivec, comm, parComm, cc);
    parComm->remoteIndices().rebuild<false>();

    std::vector<int> rowType(A_loc_.N(), 0);
    //std::vector<int> rowType(A_loc_.N(), 0);
    std::vector<int> comTab;
    getIndexSetInfo(parComm, cc, comTab, rowType);
    printComTabOnRoot(cc, comTab);

    //remove off-diagonal NNZ on ghost rows.
    buildLocalMatrixFromLoc(A_loc_, A_loc, rowType);
    printNumCells(cc, A_loc.nonzeroes(), 2);

    ScalarProduct sp(*parComm);
    Operator linOp(A_loc, *parComm);
    
    Opm::FlowLinearSolverParameters flsp;
    Opm::FlowLinearSolverParameters flsp_amg;
    Opm::FlowLinearSolverParameters flsp_cpr;
    Opm::FlowLinearSolverParameters flsp_json;
    Opm::FlowLinearSolverParameters flsp_json2;
    flsp_json.linsolver_ = std::string("/home/andreast/fork_opm/test/linear_solver_config_files/2021/cpr/sleipner_cpr_verbose_quasiimpes_tol3.json");
    flsp_json2.linsolver_ = std::string("/home/andreast/fork_opm/test/linear_solver_config_files/2021/cpr/fieldCase1.json");
    
    Opm::PropertyTree prm_ilu = setupILU(std::string("ILU"), flsp);
    
    Opm::PropertyTree prm_amg = setupAMG(std::string("amg"), flsp_amg);
    
    Opm::PropertyTree prm_json(flsp_json.linsolver_);
    Opm::PropertyTree prm_json2(flsp_json2.linsolver_);

    int verbose = 10;
    prm_ilu.put("verbosity",verbose);
    prm_json2.put("verbosity",verbose);
    prm_json.put("verbosity",verbose);
    prm_amg.put("verbosity",verbose);
    prm_amg.put("preconditioner.verbosity", verbose);

    int maxI = 1;
    prm_ilu.put("maxiter", maxI);
    prm_amg.put("maxiter", maxI);
    prm_json.put("maxiter", maxI);
    prm_json2.put("maxiter", maxI);

    std::function<Vec()> weightsCalculator;
    std::function<Vec()> quasi;

    if (rank == 0) {std::cout << "Create ILU "<< std::endl;}
    auto fs = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_ilu, weightsCalculator, 0);
    if (rank == 0) {std::cout << std::endl;}
    if (rank == 0) {std::cout << std::endl;}
    if (rank == 0) {std::cout << "Create AMG "<< std::endl;}
    auto fs_amg = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_amg, weightsCalculator, 0);

    if (rank == 0) {std::cout << std::endl;}
    if (rank == 0) {std::cout << std::endl;}
    if (rank == 0) {std::cout << "Create second AMG "<< std::endl;}
    auto fs_amg2 = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_amg, weightsCalculator, 0);

    int pIdx = 1;
    quasi = [A_loc, pIdx]() {
	return Opm::Amg::getQuasiImpesWeights<Mat, Vec>(A_loc, pIdx, false);
    };
    if (rank == 0) {std::cout << std::endl;}
    if (rank == 0) {std::cout << std::endl;}
    if (rank == 0) {std::cout << "Create CPR sleipner"<< std::endl;}
    auto fs_sleip = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_json, quasi, pIdx);
    if (rank == 0) {std::cout << std::endl;}
    if (rank == 0) {std::cout << std::endl;}
    if (rank == 0) {std::cout << "Create second CPR sleipner"<< std::endl;}
    auto fs_sleip2 = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_json, quasi, pIdx);
    if (rank == 0) {std::cout << std::endl;}
    if (rank == 0) {std::cout << std::endl;}
    if (rank == 0) {std::cout << "Create CPR fieldCase "<< std::endl;}
    auto fs_field = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_json2, quasi, pIdx);

    typedef std::numeric_limits< double > dbl;
    std::cout.precision(dbl::max_digits10);
    double tol = 0.00001;

    //solveTwice(*fs, rhs_loc, cc, sp, std::string("ILU"), tol);
    double amg_norm =solveTwice(*fs_amg, *fs_amg2, rhs_loc, cc, sp, std::string("AMG"), tol);
    if (amg_norm > 0)
	applyTwice(*fs_amg, *fs_amg2, rhs_loc, cc, sp, std::string("AMG"));

    //solveTwice(*fs_field, rhs_loc, cc, sp, std::string("CPR_FIELDCASE"), tol);
    double sleip_norm = solveTwice(*fs_sleip, *fs_sleip2, rhs_loc, cc, sp, std::string("CPR_SLEIPNER"), tol);
    if (sleip_norm > 0)
	applyTwice(*fs_sleip, *fs_sleip2, rhs_loc, cc, sp, std::string("CPR_SLEIPNER"));
    
    //solveTwice(*fs_sleip2, rhs_loc, cc, sp, std::string("CPR_SLEIPNER"), tol);
    
    return 0;
}


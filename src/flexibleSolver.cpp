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
#include "solveParallel/reorder.hpp"
#include "solveParallel/buildLocalMatrix.hpp"
#include "partition/overlapCreation.hpp"
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
    typedef Dune::BCRSMatrix<BlockMat3> Mat;
    typedef Dune::FieldVector<double,3> BlockVec;
    typedef Dune::BlockVector<BlockVec> Vec;
    typedef Dune::MPIHelper::MPICommunicator MPICommunicator;
    typedef Dune::CollectiveCommunication<MPICommunicator> CollectiveCommunication;
    typedef Dune::BiCGSTABSolver<Vec> Solver;
    typedef Dune::InverseOperatorResult Stat;
    
    typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
    typedef Dune::OverlappingSchwarzScalarProduct<Vec,Comm> ScalarProduct;
    typedef GhostLastMatrixAdapter<Mat,Vec,Vec,Comm> GLO;                   // solveParallel/ghostLastOperations.hpp
    typedef Dune::OverlappingSchwarzOperator<Mat,Vec,Vec,Comm> Operator;
    typedef Opm::ParallelOverlappingILU0<Mat,Vec,Vec,Comm> ILU;
    typedef Dune::FlexibleSolver<Mat, Vec> FlexibleSolverType;
    //typedef Dune::FlexibleSolver<GLO> FlexibleSolverType;
    
    CollectiveCommunication cc(MPI_COMM_WORLD);
    int rank = cc.rank();

    Mat A_loc, A_loc_;
    Vec rhs, rhs_loc;

    DictRead DR;
    Comm comm(cc);
    std::shared_ptr<Comm> parComm(new(Comm));
    readMatOnRootAndDist(argc, argv, A_loc, rhs_loc, DR, comm, parComm, cc); // in partition/overlapCreation.hpp
    
    ScalarProduct sp(*parComm);
    GLO linOp(A_loc, *parComm);
    
    Opm::FlowLinearSolverParameters flsp;
    Opm::FlowLinearSolverParameters flsp_amg;
    Opm::FlowLinearSolverParameters flsp_cpr;
    Opm::FlowLinearSolverParameters flsp_json;
    Opm::FlowLinearSolverParameters flsp_json2;
    flsp_json.linsolver_ = std::string("/home/andreast/fork_opm/test/linear_solver_config_files/2021/cpr/sleipner_cpr_verbose_quasiimpes_tol3.json");
    flsp_json2.linsolver_ = DR.dict[12];
    
    Opm::PropertyTree prm_ilu = setupILU(std::string("ILU"), flsp);
    prm_ilu.put("verbosity", 2);
    prm_ilu.put("maxiter", 200);
    
    Opm::PropertyTree prm_amg = setupAMG(std::string("amg"), flsp_amg);
    prm_amg.put("verbosity", 2);
    prm_amg.put("maxiter", 40);
    prm_amg.put("preconditioner.verbosity", 10);
    Opm::PropertyTree prm_cpr = setupCPR(std::string("cpr_quasiimpes"), flsp_amg);
    prm_cpr.put("verbosity", 2);
    prm_cpr.put("preconditioner.verbosity", 10);

    Opm::PropertyTree prm_json(flsp_json.linsolver_);
    Opm::PropertyTree prm_json2(flsp_json2.linsolver_);
    prm_json2.put("verbosity", 2);
    prm_json.put("maxiter", 40);
    prm_json2.put("maxiter", 40);
    
    std::function<Vec()> weightsCalculator;
    std::function<Vec()> quasi;

    Dune::Timer ILUtimer;
    if (rank == 0) {std::cout << "Create ILU "<< std::endl;}
    cc.barrier();
    ILUtimer.start();
    auto fs = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_ilu, weightsCalculator, 0);
    cc.barrier();
    double iluTime = ILUtimer.stop();
    if (rank == 0) {std::cout << "Setup ILU time: "<< iluTime<< std::endl;}

    if (rank == 0) {std::cout << "Create AMG "<< std::endl;}
    cc.barrier();
    ILUtimer.reset();
    ILUtimer.start();
    auto fs_amg = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_amg, weightsCalculator, 0);
    cc.barrier();
    double amgTime = ILUtimer.stop();
    if (rank == 0) {std::cout << "Setup AMG time: "<< amgTime<< std::endl;}
    
    int pidx = 1;

    cc.barrier();
    ILUtimer.reset();
    ILUtimer.start();
    quasi = [A_loc, pidx]() {
	return Opm::Amg::getQuasiImpesWeights<Mat, Vec>(A_loc, pidx, false);
    };
    cc.barrier();
    double wgtTime = ILUtimer.stop();
    if (rank == 0) {std::cout << "Setup quasiWgt time: "<< wgtTime<< std::endl;}

    if (rank == 0) {std::cout << std::endl;}
    if (rank == 0) {std::cout << std::endl;}
    if (rank == 0) {std::cout << "Create CPR default "<< std::endl;}
    auto fs_cpr = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_cpr, quasi, pidx);
    if (rank == 0) {std::cout << std::endl;}
    if (rank == 0) {std::cout << std::endl;}
    if (rank == 0) {std::cout << "Create CPR sleipner"<< std::endl;}
    cc.barrier();
    ILUtimer.reset();
    ILUtimer.start();
    auto fs_json = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_json, quasi, pidx);
    cc.barrier();
    double fieldTime = ILUtimer.stop();
    if (rank == 0) {std::cout << "Setup CPR_SLEIPNER time: "<< fieldTime<< std::endl;}
    
    if (rank == 0) {std::cout << std::endl;}
    if (rank == 0) {std::cout << std::endl;}
    if (rank == 0) {std::cout << "Create CPR fieldCase "<< std::endl;}
    cc.barrier();
    ILUtimer.reset();
    ILUtimer.start();
    auto fs_json2 = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_json2, quasi, pidx);
    cc.barrier();
    double sleipTime = ILUtimer.stop();
    if (rank == 0) {std::cout << "Setup CPR_FIELDCASE time: "<< sleipTime<< std::endl;}

    flexibleSolverAndPreconTimer(*fs_cpr, rhs_loc, cc, std::string("CPR_DEFAULT"), 0.01);
    flexibleSolverAndPreconTimer(*fs, rhs_loc, cc, std::string("ILU"), 0.01, true);
    flexibleSolverAndPreconTimer(*fs_amg, rhs_loc, cc, std::string("AMG"), 0.01, true);
    flexibleSolverAndPreconTimer(*fs_json2, rhs_loc, cc, std::string("CPR_FIELDCASE"), 0.01, true);
    flexibleSolverAndPreconTimer(*fs_json, rhs_loc, cc, std::string("CPR_SLEIPNER"), 0.01, true);

    return 0;
}

/*
  Copyright 2019 Andreas Thune.

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
#include <dune/istl/schwarz.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>


#include <dune/common/function.hh>

#include <mpi.h>
#include <zoltan.h>

//#include <opm/autodiff/ParallelOverlappingILU0.hpp>
#include <opm/simulators/linalg/ParallelOverlappingILU0.hpp>

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
    //typedef Opm::GhostLastParallelOverlappingILU0<Mat,Vec,Vec,Comm> Pre;
    typedef Opm::ParallelOverlappingILU0<Mat,Vec,Vec,Comm> Pre;

    typedef GhostLastOperator<Mat,Vec,Vec> GLOp;
    typedef GhostLastScalarProduct<Vec,CollectiveCommunication> GLSP;

    typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<Mat,Dune::Amg::FirstDiagonal> > Criterion;
    typedef Dune::SeqILU<Mat,Vec,Vec> Smoother;
    typedef Dune::SeqJac<Mat,Vec,Vec> JacSmoother;
    typedef Dune::BlockPreconditioner<Vec, Vec, Comm, Smoother> ParSmoother;
    typedef Dune::BlockPreconditioner<Vec, Vec, Comm, JacSmoother> ParJacSmoother;
    typedef typename Dune::Amg::SmootherTraits<ParSmoother>::Arguments Smoother1Args;
    typedef typename Dune::Amg::SmootherTraits<ParJacSmoother>::Arguments Smoother2Args;
    typedef typename Dune::Amg::SmootherTraits<Pre>::Arguments Smoother3Args;

    typedef Dune::Amg::AMG<Operator,Vec,ParSmoother,Comm> AMGILU0;
    typedef Dune::Amg::AMG<Operator,Vec,ParJacSmoother,Comm> AMGJac;
    typedef Dune::Amg::AMG<Operator,Vec,Pre,Comm> AMG_GLILU;
    
    CollectiveCommunication cc(MPI_COMM_WORLD);
    int rank = cc.rank();

    Mat A, A_loc, A_rec, A_all;
    Mat1 trans, wells;
    Vec rhs, rhs_loc;

    DictRead DR;

    // Read matrices
    handleMatrixSystemInput(argc, argv, A, trans, wells, rhs, DR, rank);
    DR.dict[5] = "2";
    std::cout << std::numeric_limits<float>::max() << std::endl;

    std::vector<int> row_size(A.N());
    storeRowSize(A, row_size);

    //partition matrix
    std::vector<int> mpivec(A.N(), rank);
    zoltanPartitionFunction(mpivec, trans, wells, cc, DR, row_size);
    evaluatePartition(trans, wells, mpivec, cc);
    
    std::vector<std::set<int>> overlap(A.N());
    std::vector<int> overlapMap, local2global, global2local, reorder, l2r, r2l;

    //Add overlap
    int partSize = addOverlap(overlap, overlapMap, mpivec, A, cc);
    myIds(local2global, global2local, overlapMap, cc);
    myRecIds(local2global, global2local, overlapMap, l2r, r2l);
    do_reorder(A, reorder, local2global, global2local, overlapMap, DR);

    //Set up parallel index set stuff
    PIS indexSet;
    setParallelLocalIndex(indexSet, local2global, overlapMap);
    RIS remoteIndexSet(indexSet, indexSet, cc);
    remoteIndexSet.rebuild<true>();
    auto info = getParallelInfoReorder(indexSet, remoteIndexSet, reorder);

    bool assembleGhost = std::stoi(DR.dict[9]) == 1;
    //build local matrix and vector
    buildLocalMatrixReorder(A, A_loc, overlapMap, local2global, global2local, rank, reorder, false);
    buildLocalMatrixReorder(A, A_all, overlapMap, local2global, global2local, rank, reorder, true);
    buildLocalVectorReorder(rhs, rhs_loc, overlapMap, local2global, rank, reorder);
    buildLocalMatrixRectangular(A, A_rec, overlapMap, local2global, global2local, rank, reorder, partSize, l2r, r2l);

    printNumCells(cc, A_rec.nonzeroes(), 3);
    printNumCells(cc, A_loc.nonzeroes(), 2);

    printNumCells(cc, local2global.size(), 0);
    printNumCells(cc, r2l.size(), 1);

    Comm comm(info, cc);
    
    std::string use_ilu("ILU");
    auto ilu_help = Opm::convertString2Milu(use_ilu);
    std::unique_ptr<Pre> preCon(new Pre(A_loc, comm, 0.99, ilu_help, local2global.size(), false, false) );
    std::unique_ptr<Pre> GLpre(new Pre(A_loc, comm, 0.99, ilu_help, r2l.size(), false, false) );

    ScalarProduct sp(comm);
    Operator linOp(A_loc,comm);
    Operator linOpAll(A_all,comm);
    
    Solver oneIterBicg(linOp, sp, *preCon, 1e-5, 1, 0);
    Solver hundredIterBiCG(linOp,sp,*preCon,1e-40, 100, 0);
    
    Stat statistics;

    //GLOp ghostLastOp(A_loc, r2l.size(), comm);
    GLSP ghostLastSp(cc, r2l.size());
    //Solver ghostLastSolver(ghostLastOp, ghostLastSp, *GLpre, 1e-5, 1, 0);
    //Solver hundredIterGLSolver(ghostLastOp, ghostLastSp, *GLpre, 1e-40, 100, 0);

    // ---------- Set up AMG ----------
    
    Smoother1Args smoother1Args;
    Smoother2Args smoother2Args;
    Smoother3Args smoother3Args;

    
    smoother1Args.iterations = 1;
    Criterion criterion(15, 100);
    criterion.setDefaultValuesIsotropic(2);

    if (rank == 0) {std::cout << std::endl;}
    std::unique_ptr<AMGILU0> amg1( new AMGILU0(linOp, criterion, smoother1Args, comm) );

    if (rank == 0) {std::cout << std::endl;}
    std::unique_ptr<AMGILU0> amg1_all( new AMGILU0(linOpAll, criterion, smoother1Args, comm) );

    if (rank == 0) {std::cout << std::endl;}
    std::unique_ptr<AMGJac> amg2( new AMGJac(linOp, criterion, smoother2Args, comm) );

    if (rank == 0) {std::cout << std::endl;}
    std::unique_ptr<AMGJac> amg2_all( new AMGJac(linOpAll, criterion, smoother2Args, comm) );
    
    if (rank == 0) {std::cout << std::endl;}
    //AMG_GLILU amg3(linOp, criterion, smoother3Args, comm);
    
    if (rank == 0) {std::cout << std::endl;}
    
    // --------------------------------

    
    int numTimeOp = std::stoi(DR.dict[10]);
    
    Vec x(rhs_loc.size());
    Vec crhs(rhs_loc);

    //Solver amgJacSolver(linOp, ghostLastSp, *amg1, 1e-3, 200, 2);
    int verbose=0;
    if (rank==0) {verbose=2;}
    Solver amgJacSolver(linOp, sp, *amg1, 1e-3, 200, verbose);
    x=0;
    Dune::InverseOperatorResult stat;
    amgJacSolver.apply(x,crhs,stat);
    x=0;
    amg1->pre(x,rhs_loc);
    x=0;
    amg2->pre(x,rhs_loc);
    x=0;
    amg1_all->pre(x,rhs_loc);
    x=0;
    amg2_all->pre(x,rhs_loc);

    //time operator.apply
    cc.barrier();
    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimeSpMV(cc, linOp, rhs_loc, numTimeOp, false);


    std::cout << A_loc.nonzeroes() << " " << A_all.nonzeroes() << std::endl;
    amg1->apply(rhs_loc, x);
    cc.barrier();
    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimePre(cc, *amg1, rhs_loc, numTimeOp);

    cc.barrier();
    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimePre(cc, *amg1_all, rhs_loc, numTimeOp);

    cc.barrier();
    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimePre(cc, *amg2, rhs_loc, numTimeOp, true, false);

    cc.barrier();
    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimePre(cc, *amg2_all, rhs_loc, numTimeOp);

    return 0;
    
}

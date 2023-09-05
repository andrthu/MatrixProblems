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
#include <string>

#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/parallel/indexset.hh>
#include <dune/common/timer.hh>
#include <dune/common/precision.hh>

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
#include "io/wells.hpp"
#include "partition/graphFunctions.hpp"
#include "partition/overlapCreation.hpp"
#include "solveParallel/reorder.hpp"
#include "solveParallel/buildLocalMatrix.hpp"
#include "partition/setIndexSet.hpp"
#include "partition/evaluatePartition.hpp"
#include "solveParallel/solveSystem.hpp"
#include "solveParallel/wellSolveSystem.hpp"
#include "io/dictRead.hpp"
#include "solveParallel/compareSolutions.hpp"
#include "partition/meassureCommunication.hpp"


int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    Dune::FMatrixPrecision<double> pres;
    //pres.set_singular_limit(1E-20);

    typedef Dune::FieldMatrix<double,3,3> BlockMat3;
    typedef Dune::FieldMatrix<double,4,4> WellBlockD;
    typedef Dune::FieldMatrix<double,4,3> WellBlockOff;
    typedef Dune::BCRSMatrix<WellBlockOff> WellMat;
    typedef Dune::FieldMatrix<double,1,1> BlockMat1;
    typedef Dune::BCRSMatrix<BlockMat3> Mat;
    typedef Dune::BCRSMatrix<BlockMat1> Mat1;
    typedef Dune::FieldVector<double,3> BlockVec;
    typedef Dune::FieldVector<double,4> WellVec;
    typedef Dune::BlockVector<BlockVec> Vec;
    typedef Dune::MPIHelper::MPICommunicator MPICommunicator;
    typedef Dune::CollectiveCommunication<MPICommunicator> CollectiveCommunication;
    typedef Dune::OwnerOverlapCopyAttributeSet::AttributeSet AttributeSet;
    typedef Dune::ParallelIndexSet<int,Dune::ParallelLocalIndex<AttributeSet>> PIS;
    typedef Dune::RemoteIndices<PIS> RIS;
    typedef Wells<WellBlockD,WellMat,WellBlockOff,WellVec> MyWells;
    typedef Dune::BiCGSTABSolver<Vec> Solver;
    typedef Dune::InverseOperatorResult Stat;
    
    typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
    typedef Dune::OverlappingSchwarzScalarProduct<Vec,Comm> ScalarProduct; 
    typedef WellModelMatrixAdapter<Mat,Vec,Vec,MyWells,Comm> Operator;
    typedef Opm::ParallelOverlappingILU0<Mat,Vec,Vec,Comm> Pre;

    typedef GhostLastOperator<Mat,Vec,Vec> GLOp;
    typedef GhostLastScalarProduct<Vec,CollectiveCommunication> GLSP;

    CollectiveCommunication cc(MPI_COMM_WORLD);
    int rank = cc.rank();
    
    Mat A, A_loc;    
    Mat1 trans, wellAdj;
    Vec rhs, rhs_loc, x, x_seq;
    MyWells wells, wells_loc;
    

    DictRead DR;
    
    handleMatrixSystemInput(argc, argv, A, trans, wellAdj, rhs, DR, rank);
    
    std::vector<int> row_size(A.N());
    storeRowSize(A, row_size);

    wells.readWellsFromFiles(argv[1],A.N());
    if (rank==0)
	std::cout << "Active wells: " << wells.numWells() << std::endl;
    
    std::vector<int> mpivec(A.N(), rank);
    zoltanPartitionFunction(mpivec, trans, wellAdj, cc, DR, row_size);

    std::vector<std::set<int>> overlap(A.N());
    std::vector<int> overlapMap;
    std::vector<int> local2global, global2local, l2r, r2l;
    
    evaluatePartition(trans, wellAdj, mpivec, cc);
    addOverlap(overlap, overlapMap, mpivec, A, cc);
    myIds(local2global, global2local, overlapMap, cc);
    myRecIds(local2global, global2local, overlapMap, l2r, r2l);
    
    std::cout << "After partitioning rank " << rank << " has " << local2global.size() << " Blocks." << std::endl;
    std::cout << "Rank "<< rank << " InteriorSize: " << r2l.size()<<std::endl;

    bool assembleGhost = std::stoi(DR.dict[9]) == 1;
    buildLocalMatrix(A, A_loc, overlapMap, local2global, global2local, rank, assembleGhost);
    buildLocalVector(rhs, rhs_loc, overlapMap, local2global, rank);

    buildLocalWells(wells, wells_loc, mpivec, local2global, global2local, rank);
    std::cout << "WellsOnRank "<< rank <<": "<< wells_loc.numWells()<< std::endl;
    std::cout << "Rank "<< rank << " nnz: " << A_loc.nonzeroes() << std::endl;


    PIS indexSet;
    setParallelLocalIndex(indexSet, local2global, overlapMap);

    RIS remoteIndexSet(indexSet, indexSet,cc);
    remoteIndexSet.rebuild<true>();
    auto info = getParallelInfo(indexSet, remoteIndexSet);

    Comm comm(info, cc);
    ScalarProduct sp(comm);
    Operator linOp(A_loc, wells_loc, comm);
    
    std::string use_ilu("ILU");
    std::unique_ptr<Pre> preCon(new Pre(A_loc, comm, 0.99, Opm::convertString2Milu(use_ilu)));

    Solver oneIterBicg(linOp, sp, *preCon, 1e-5, 1, 0);

    int numTimeOp = std::stoi(DR.dict[10]);

    //time copyOwnerToAll call
    cc.barrier();
    multipleMeassureComm(cc, info, local2global.size(), numTimeOp);    
    
    //time A.apply + well.apply
    cc.barrier();
    if (rank == 0) {std::cout << std::endl;}
    multipleTimeWellSpMV(cc, info, linOp, rhs_loc, numTimeOp);

    //time pre.apply
    cc.barrier();
    if (rank == 0) {std::cout << std::endl;}
    multipleTimePre(cc, *preCon, rhs_loc, numTimeOp);
    
    //time well.apply
    cc.barrier();
    if (rank == 0) {std::cout << std::endl;}
    multipleTimeWell(cc, wells_loc, rhs_loc, numTimeOp);
    
    //time scalar product
    cc.barrier();
    if (rank == 0) {std::cout << std::endl;}    
    multipleTimeScalarProduct(cc, info, rhs_loc, numTimeOp);

    cc.barrier();
    if (rank == 0) {std::cout << std::endl;}    
    multipleTimeSolver(cc, oneIterBicg, rhs_loc, 10);

    //wellSolveSystem(A_loc, wells_loc, rhs_loc, x, cc, info, DR);

    /*
    Vec y(A.N());
    y=0;
    wells.apply(rhs,y);
    */

    return 0;
}

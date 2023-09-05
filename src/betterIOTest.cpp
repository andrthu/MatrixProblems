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

    typedef Dune::FieldMatrix<double,2,2> BlockMat3;
    typedef Dune::FieldMatrix<double,1,1> BlockMat1;
    typedef Dune::BCRSMatrix<BlockMat3> Mat;
    typedef Dune::BCRSMatrix<BlockMat1> Mat1;
    typedef Dune::FieldVector<double,2> BlockVec;
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
    //typedef Dune::FlexibleSolver<Operator> FlexibleSolverType;
    
    CollectiveCommunication cc(MPI_COMM_WORLD);
    int rank = cc.rank();

    //std::unique_ptr<Mat> A(new Mat);
    Mat* A = new Mat;
    Mat A_loc;
    Mat1 trans, wells;
    Vec* rhs = new Vec;
    Vec rhs_loc;

    DictRead DR;

    // Read matrices
    handleMatrixSystemInputSomeRanks(argc, argv, *A, trans, wells, *rhs, DR, cc, rank%2==0);
    DR.dict[5] = "2";

    int N;
    if (rank == 0)
	N = A->N();
    cc.broadcast(&N, 1, 0);
    std::cout << rank << " num rows: "<< N << std::endl;
    std::vector<int> row_size(N);
    storeRowSizeFromRoot(*A, row_size, cc);

    std::cout << rank << " "<< row_size[0] << std::endl;
    //partition matrix
    std::vector<int> mpivec(N, rank);
    zoltanPartitionFunction(mpivec, trans, wells, cc, DR, row_size);
    
    evaluatePartition(trans, wells, mpivec, cc);
    
    std::vector<std::set<int>> overlap(N);
    std::vector<int> overlapMap, local2global, global2local, reorder, l2r, r2l;
    //bool assembleGhost = std::stoi(DR.dict[9]) == 1;
    if (rank%8==0) {
    
	//Add overlap
	int partSize = addOverlap(overlap, overlapMap, mpivec, *A, cc);
	myIds(local2global, global2local, overlapMap, cc);
	myRecIds(local2global, global2local, overlapMap, l2r, r2l);
	do_reorder(*A, reorder, local2global, global2local, overlapMap, DR);

	buildLocalMatrixReorder(*A, A_loc, overlapMap, local2global, global2local, rank, reorder, false);
	buildLocalVectorReorder(*rhs, rhs_loc, overlapMap, local2global, rank, reorder);
	//Mat* A_ = A.release();
	delete A;
	delete rhs;
	A = NULL;
	rhs = NULL;
    }

    cc.barrier();
    for (int modRank = 1; modRank < 8; ++modRank) {
	if (rank == 0) {std::cout<< "modrank "<< modRank <<std::endl;} 
	if (rank%8 == modRank) {
	    handleMatrixSystemInputSomeRanks(argc, argv, *A, trans, wells, *rhs, DR, cc, rank%2!=0);

	    //Add overlap
	    int partSize = addOverlap(overlap, overlapMap, mpivec, *A, cc);
	    myIds(local2global, global2local, overlapMap, cc);
	    myRecIds(local2global, global2local, overlapMap, l2r, r2l);
	    do_reorder(*A, reorder, local2global, global2local, overlapMap, DR);
	

	    buildLocalMatrixReorder(*A, A_loc, overlapMap, local2global, global2local, rank, reorder, false);
	    buildLocalVectorReorder(*rhs, rhs_loc, overlapMap, local2global, rank, reorder);
	    //Mat* A_ = A.release();
	    delete A;
	    delete rhs;
	    A = NULL;
	    rhs = NULL;
	}
	cc.barrier();
    }
    
    PIS indexSet;
    setParallelLocalIndex(indexSet, local2global, overlapMap);
    RIS remoteIndexSet(indexSet, indexSet, cc);
    remoteIndexSet.rebuild<true>();
    auto info = getParallelInfoReorder(indexSet, remoteIndexSet, reorder);
    printNumCells(cc, A_loc.nonzeroes(), 2);
    printNumCells(cc, local2global.size(), 0);
    printNumCells(cc, r2l.size(), 1);
    
    return 0;
}

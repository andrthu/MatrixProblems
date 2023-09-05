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
#include <dune/istl/repartition.hh>
#include <dune/istl/matrixredistribute.hh>
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


template<class M, class M1>
void createMat1(const M& orig, M1& newM)
{
    
    auto createIter =  newM.createbegin();
    for (const auto& row : orig) {
        for (auto col = row.begin(), cend = row.end(); col != cend; ++col) {
            createIter.insert(col.index());
        }
        ++createIter;
    }
    auto rowCoarse = newM.begin();
    for (auto row = orig.begin(), rowEnd = orig.end(); row != rowEnd; ++row, ++rowCoarse) {
	auto entryCoarse = newM.begin();
	for (auto entry = row->begin(), entryEnd = row->end(); entry != entryEnd; ++entry, ++entryCoarse) {
	    double matrix_el = (*entry)[0][0];
	    (*entryCoarse) = matrix_el;
	}
    }
}


int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    typedef Dune::FieldMatrix<double,3,3> BlockMat3;
    typedef Dune::FieldMatrix<double,1,1> BlockMat1;
    typedef Dune::BCRSMatrix<BlockMat3> Mat;
    typedef Dune::BCRSMatrix<BlockMat1> Mat1;
    typedef Dune::FieldVector<double,3> BlockVec;
    typedef Dune::FieldVector<double,3> BlockVec1;
    typedef Dune::BlockVector<BlockVec> Vec;
    typedef Dune::BlockVector<BlockVec1> Vec1;
    typedef Dune::MPIHelper::MPICommunicator MPICommunicator;
    typedef Dune::CollectiveCommunication<MPICommunicator> CollectiveCommunication;
    typedef Dune::OwnerOverlapCopyAttributeSet::AttributeSet AttributeSet;
    typedef Dune::ParallelIndexSet<int,Dune::ParallelLocalIndex<AttributeSet>> PIS;
    typedef Dune::RemoteIndices<PIS> RIS;
    typedef Dune::BiCGSTABSolver<Vec> Solver;
    typedef Dune::InverseOperatorResult Stat;
    typedef Dune::Amg::MatrixGraph<Mat> MatrixGraph;
    
    typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
    typedef Dune::OverlappingSchwarzScalarProduct<Vec,Comm> ScalarProduct; 
    typedef Dune::OverlappingSchwarzOperator<Mat,Vec,Vec,Comm> Operator;
    typedef Dune::OverlappingSchwarzOperator<Mat1,Vec1,Vec1,Comm> Operator1;
    typedef Opm::ParallelOverlappingILU0<Mat,Vec,Vec,Comm> ILU;
    typedef Opm::ParallelOverlappingILU0<Mat1,Vec1,Vec1,Comm> ILU1;
    typedef Dune::FlexibleSolver<Mat, Vec> FlexibleSolverType;
    
    CollectiveCommunication cc(MPI_COMM_WORLD);
    int rank = cc.rank();

    //std::unique_ptr<Mat> A(new Mat);
    Mat A;
    Mat A_loc;
    Mat1 trans, wells;
    Vec rhs;
    
    DictRead DR;

    // Read matrices
    handleMatrixSystemInputSomeRanks(argc, argv, A, trans, wells, rhs, DR, cc, rank==0);
    DR.dict[5] = "2";
    
    int N;
    if (rank == 0)
	N = A.N();
    cc.broadcast(&N, 1, 0);
    std::cout << rank << " num rows: "<< N << std::endl;
    std::vector<int> row_size(N);
    storeRowSizeFromRoot(A, row_size, cc);

    std::cout << rank << " "<< row_size[0] << std::endl;
    //partition matrix
    std::vector<int> mpivec(N, rank);
    zoltanPartitionFunction(mpivec, trans, wells, cc, DR, row_size);

    Comm comm(cc);

    std::shared_ptr<Comm> parComm(new(Comm));
    Vec rhs_loc;
    constructLocalFromRoot(A, A_loc, N, rhs, rhs_loc, mpivec, comm, parComm, cc);
    parComm->remoteIndices().rebuild<false>();
    
    
    Mat1 A1(A_loc.N(), A_loc.N(), Mat1::row_wise);
    createMat1(A_loc, A1);

    Vec1 rhs1(rhs_loc.size());
    for (int i = 0; i < rhs_loc.size(); ++i)
	rhs1[i] = rhs_loc[i][1];

    
    Operator linOp(A_loc, *parComm);
    Operator1 linOp1(A1, *parComm);
    //Operator1 linOpT(trans, comm);
    
    multipleMinLoopTimeSpMV(cc, linOp, rhs_loc, 3);
    multipleMinLoopTimeSpMV(cc, linOp1, rhs1, 3);
    //multipleMinLoopTimeSpMV(cc, linOpT, rhs1, 3);

    std::string use_ilu("ILU");
    auto ilu_help = Opm::convertString2Milu(use_ilu);
    std::unique_ptr<ILU> preCon(new ILU(A, *parComm, 0.99, ilu_help, rhs_loc.size(), false, false) );
    std::unique_ptr<ILU1> preCon1(new ILU1(A1, *parComm, 0.99, ilu_help, rhs_loc.size(), false, false) );

    multipleMinLoopTimePre(cc, *preCon, rhs_loc, 3);
    multipleMinLoopTimePre(cc, *preCon1, rhs1, 3);
    
    return 0;
}

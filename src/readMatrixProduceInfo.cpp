/*
  Copyright 2018 Andreas Thune.

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

#include "io/inputRead.hpp"
#include "partition/graphFunctions.hpp"
#include "partition/overlapCreation.hpp"
#include "solveParallel/reorder.hpp"
#include "solveParallel/buildLocalMatrix.hpp"
#include "partition/setIndexSet.hpp"
#include "partition/evaluatePartition.hpp"
#include "solveParallel/solveSystem.hpp"
#include "io/dictRead.hpp"
#include "solveParallel/compareSolutions.hpp"
#include "partition/meassureCommunication.hpp"


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
    typedef Dune::MatrixIndexSet AdjecencyPattern;
    
    CollectiveCommunication cc(MPI_COMM_WORLD);
    int rank = cc.rank();

    Mat A, A_loc, T;
    Mat1 trans, wells;
    Vec rhs, rhs_loc, x, x_seq;
    
    DictRead DR;
    
    handleMatrixSystemInput(argc, argv, A, trans, wells, rhs, DR, rank);
    
    //printMatrixAndInfo(A,trans);
    //writeMatrixRowSizeHist(A);
    std::vector<int> row_size(A.N());
    storeRowSize(A, row_size);

    std::vector<int> mpivec(A.N(), rank);
    zoltanPartitionFunction(mpivec, trans, wells, cc, DR, row_size);
    
    std::vector<std::set<int>> overlap(A.N());
    std::vector<int> overlapMap;
    std::vector<int> local2global, global2local, reorder;
    
    evaluatePartition(trans, wells, mpivec, cc);
    addOverlap(overlap, overlapMap, mpivec, A, cc);
    myIds(local2global, global2local, overlapMap, cc);
    do_reorder(A, reorder, local2global, global2local, overlapMap, DR);

    std::cout << "After partitioning rank " << rank << " has " << local2global.size() << " Blocks." << std::endl;

    bool assembleGhost = std::stoi(DR.dict[9]) == 1;
    //build local matrices and vectors
    buildLocalMatrixReorder(A, A_loc, overlapMap, local2global, global2local, rank, reorder, assembleGhost);
    buildLocalVectorReorder(rhs, rhs_loc, overlapMap, local2global, rank, reorder);
    
    PIS indexSet;
    setParallelLocalIndex(indexSet, local2global, overlapMap);

    RIS remoteIndexSet(indexSet, indexSet,cc);
    remoteIndexSet.rebuild<true>();
    auto info = getParallelInfoReorder(indexSet, remoteIndexSet, reorder);
    
    //solveSystem(A_loc, rhs_loc, x, cc, info, DR);

    if (cc.size()>1 && std::stoi(DR.dict[6])==1)
    {
	solveSystem(A, rhs, x_seq, cc, DR);
	compareSolutions(x,x_seq,cc,local2global,overlapMap);
    }



    std::cout << A_loc.nonzeroes() << std::endl;

    //A_loc.j_ = NULL;


    if (rank == 0) {
	std::string test("start");

	test.push_back('_');
	test+=std::to_string(A.N());
	test.push_back('_');
	test += std::to_string(cc.size());
	test += std::string(".txt");
	
	std::cout <<test.data()<<std::endl;

	std::ofstream ofile(test.data());

	ofile << "hihi \n";

	ofile.close();

	Mat1 t1,t2;

	size_t size = 5;
	AdjecencyPattern adjDiag, adjt2;
	adjDiag.resize(size,size);
	adjt2.resize(size,size);
	for (int idx = 0; idx < size; ++idx) {
	    adjDiag.add(idx, idx);
	    adjt2.add(idx, idx);
	    if (idx < size -1)
		adjt2.add(idx, idx + 1);
	}

	adjDiag.exportIdx(t2);
	adjt2.exportIdx(t1);

	t1=0;
	t2=1;

	size_t teller = 0;
	auto rowt2 = t2.begin();
	for (auto rowt = t1.begin(); rowt!= t1.end(); ++rowt, ++rowt2) {

	    if (teller % 2 == 0) {
		*rowt = *rowt2; 
	    }
	    teller += 1;
	}

	t1[0]=2;
	t1[0][0]=3;

	for (auto rowt = t1.begin(); rowt!= t1.end(); ++rowt, ++rowt2) {

	    for (auto col = rowt->begin(); col!=rowt->end(); ++col)
		std::cout << *col << " ";
	    std::cout <<std::endl;
	}
	
    }
    return 0;
}


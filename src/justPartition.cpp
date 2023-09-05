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

#include "io/inputRead.hpp"
#include "partition/graphFunctions.hpp"
#include "partition/overlapCreation.hpp"
#include "partition/setIndexSet.hpp"
#include "partition/evaluatePartition.hpp"
#include "io/dictRead.hpp"

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    typedef Dune::MPIHelper::MPICommunicator MPICommunicator;
    typedef Dune::CollectiveCommunication<MPICommunicator> CollectiveCommunication;
    
    typedef Dune::FieldMatrix<double,1,1> BlockMat3;
    typedef Dune::BCRSMatrix<BlockMat3> Mat;

    CollectiveCommunication cc(MPI_COMM_WORLD);
    
    Mat A;
    Mat dummy;
    
    if (argc == 2)
	readMatMarketObject(A, argv[1]);
    else
	return 1;

    std::cout << A.N() <<std::endl;

    DictRead DR;
    DR.dict[0] = "0";
    DR.dict[4] = "0";
    DR.dict[7] = "0";
    
    int rank = cc.rank();
    std::vector<int> row_size(A.N(), 1);
    std::vector<int> mpivec(A.N(), rank);
    zoltanPartitionFunction(mpivec, A, dummy, cc, DR, row_size);
    evaluatePartition(A,dummy,mpivec,cc);
    
    return 0;
}

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

#ifndef OPM_PROBLEM_HEADER_INCLUDED
#define OPM_PROBLEM_HEADER_INCLUDED

#endif // OPM_PROBLEM_HEADER_INCLUDED

class Problem
{
public:

    //Generic types
    typedef Dune::FieldMatrix<double,1,1> Bloc;
    typedef Dune::BCRSMatrix<BlockMat1> Graph;
    typedef Dune::MPIHelper::MPICommunicator MPICommunicator;
    typedef Dune::CollectiveCommunication<MPICommunicator> CollectiveCommunication;
    typedef Dune::OwnerOverlapCopyAttributeSet::AttributeSet AttributeSet;
    typedef Dune::ParallelIndexSet<int,Dune::ParallelLocalIndex<AttributeSet>> PIS;
    typedef Dune::RemoteIndices<PIS> RIS;
    typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;

    //Black-oil types
    typedef Dune::FieldMatrix<double,3,3> BlockMat3;
    typedef Dune::BCRSMatrix<BlockMat3> Mat;
    typedef Dune::FieldVector<double,3> BlockVec;
    typedef Dune::BlockVector<BlockVec> Vec;
    typedef Dune::BiCGSTABSolver<Vec> Solver;
    typedef Dune::OverlappingSchwarzScalarProduct<Vec,Comm> ScalarProduct;
    typedef Dune::OverlappingSchwarzOperator<Mat,Vec,Vec,Comm> Operator;
    typedef GhostLastMatrixAdapter<Mat,Vec,Vec,Comm> GLO;
    typedef Opm::ParallelOverlappingILU0<Mat,Vec,Vec,Comm> ILU;
    typedef Dune::Amg::AMGCPR<GLO, Vec, ILU, Comm> AMGCPR;
    typedef Dune::FlexibleSolver<Mat, Vec> FlexibleSolverType;  
    
    //Co2 types
    typedef Dune::FieldMatrix<double,2,2> BlockMat2;
    typedef Dune::BCRSMatrix<BlockMat2> MatCO2;
    typedef Dune::FieldVector<double,2> BlockVecCO2;
    typedef Dune::BlockVector<BlockVecCO2> VecCO2;
    typedef Dune::BiCGSTABSolver<VecCO2> SolverCO2;
    typedef Dune::OverlappingSchwarzScalarProduct<VecCO2,Comm> ScalarProductCO2;
    typedef Dune::OverlappingSchwarzOperator<MatCO2,VecCO2,VecCO2,Comm> OperatorCO2;
    typedef GhostLastMatrixAdapter<MatCO2,VecCO2,VecCO2,Comm> GLOCO2;
    typedef Opm::ParallelOverlappingILU0<MatCO2,VecCO2,VecCO2,Comm> ILUCO2;
    typedef Dune::Amg::AMGCPR<GLOCO2, VecCO2, ILUCO2, Comm> AMGCPRCO2;
    typedef Dune::FlexibleSolver<MatCO2, VecCO2> FlexibleSolverTypeCO2;  
    
}

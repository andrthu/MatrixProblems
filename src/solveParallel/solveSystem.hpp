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

#ifndef OPM_SOLVESYSTEM_HEADER_INCLUDED
#define OPM_SOLVESYSTEM_HEADER_INCLUDED

#endif // OPM_SOLVESYSTEM_HEADER_INCLUDED

template<class Mat, class Vec, class C, class PI, class D>
void solveSystem(Mat& A, Vec& rhs, Vec& x, C cc, PI info, D dr)
{
    typedef Dune::OwnerOverlapCopyCommunication<int,int>        Comm;
    typedef Dune::OverlappingSchwarzScalarProduct<Vec,Comm>     ScalarProduct; 
    typedef Dune::OverlappingSchwarzOperator<Mat,Vec,Vec,Comm>  Operator;
    typedef Dune::SeqILU<Mat,Vec,Vec>                          ILU; 
    typedef Dune::BlockPreconditioner<Vec,Vec,Comm,ILU>         Pre;
    typedef Dune::BiCGSTABSolver<Vec>                           Solver;
    typedef Dune::InverseOperatorResult                         Stat;

    int linItr = std::stoi(dr.dict[2]);
    double tol = std::stod(dr.dict[1]);

    int rank = cc.rank();
    int verbose = 0;
    if (rank==0)
	verbose = 2;

    Comm           comm(info, cc);
    ScalarProduct  sp(comm);
    Operator       linOp(A, comm);
    ILU            SILU(A, 0.99);
    Pre            BJILU(SILU, comm);
    Solver         bicg(linOp, sp, BJILU, tol, linItr, verbose);
    Stat           statistics;
    
    x.resize(A.N());
    x = 0;
    
    bicg.apply(x, rhs, statistics);
    
    if (rank==0)
    {
	double par_time = statistics.elapsed;
	int par_iter = statistics.iterations;
	double par_time_per_iter = par_time/par_iter;

	std::cout << "Iterations: " << par_iter << std::endl;
	std::cout << "SolveTime: " << par_time << std::endl;
	std::cout << "TimePerIter: " << par_time_per_iter << std::endl;
    }
}

template<class Mat, class Vec, class C, class D>
void solveSystem(Mat& A, Vec& rhs, Vec& x, C cc, D dr)
{
    typedef Dune::MatrixAdapter<Mat,Vec,Vec>                    Operator; 
    typedef Dune::SeqILU<Mat,Vec,Vec>                          ILU; 
    typedef Dune::BiCGSTABSolver<Vec>                           Solver;
    typedef Dune::InverseOperatorResult                         Stat;

    int linItr = std::stoi(dr.dict[2]);
    double tol = std::stod(dr.dict[1]);

    int rank = cc.rank();
    int verbose = 0;

    Operator       linOp(A);
    ILU            ilu(A, 0.99);
    Solver         bicg(linOp, ilu, tol, linItr, verbose);
    Stat           statistics;
       
    x.resize(A.N());
    x=0;
    
    bicg.apply(x, rhs, statistics);
    
    if (rank==0)
    {
	double mes_time = statistics.elapsed;
	int mes_iter = statistics.iterations;
	double time_per_iter = mes_time/mes_iter;
	std::cout << std::endl;
	std::cout << "Seq_iterations: " << mes_iter << std::endl;
	std::cout << "Seq_time: " << mes_time << std::endl;
	std::cout << "Seq_timePerIter: " << time_per_iter << std::endl;
    }
}

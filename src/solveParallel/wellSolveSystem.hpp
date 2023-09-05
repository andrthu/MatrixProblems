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

#ifndef OPM_WELLSOLVESYSTEM_HEADER_INCLUDED
#define OPM_WELLSOLVESYSTEM_HEADER_INCLUDED

#endif // OPM_WELLSOLVESYSTEM_HEADER_INCLUDED

template<class M, class X, class Y, class Well, class Comm>
class WellModelMatrixAdapter : public Dune::AssembledLinearOperator<M,X,Y>
{
public:
    
    enum {
	//! \brief The solver category.
	category =Dune::SolverCategory::overlapping};

    typedef M matrix_type;
    typedef X domain_type;
    typedef Y range_type;
    typedef typename X::field_type field_type;
    typedef Comm comunication_type;
    
    WellModelMatrixAdapter (const M& A, const Well& wellMod, const Comm& comm)
	: A_( A ), well_( wellMod ), comm_( comm )
    {}
    
    virtual void apply( const X& x, Y& y ) const
    {
	A_.mv( x, y );
	
	// add well model modification to y
	well_.apply(x, y );
	comm_.project( y );
    }

    virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
    {
	A_.usmv(alpha,x,y);
	
	// add scaled well model modification to y
	//well_.applyScaleAdd( alpha, x, y );
	comm_.project( y );
    }
    
    virtual const matrix_type& getmat() const { return A_; }

protected:
    const matrix_type& A_ ;
    //const matrix_type& A_for_precond_ ;
    const Well& well_;
    const Comm& comm_;
};

template<class Mat, class W, class Vec, class C, class PI, class D>
void wellSolveSystem(Mat& A, W& w, Vec& rhs, Vec& x, C cc, PI info, D dr)
{
    typedef Dune::OwnerOverlapCopyCommunication<int,int>        Comm;
    typedef Dune::OverlappingSchwarzScalarProduct<Vec,Comm>     ScalarProduct; 
    //typedef Dune::OverlappingSchwarzOperator<Mat,Vec,Vec,Comm>  Operator;
    typedef WellModelMatrixAdapter<Mat,Vec,Vec,W,Comm>          Operator;
    //typedef Dune::SeqILU0<Mat,Vec,Vec>                          ILU; 
    typedef Dune::SeqJac<Mat,Vec,Vec> ILU;
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
    Operator       linOp(A, w, comm);
    //ILU            SILU(A, 0.99);
    ILU            SILU(A, 1, 0.99);
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

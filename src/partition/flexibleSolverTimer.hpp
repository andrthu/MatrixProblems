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

#ifndef OPM_FLEXIBLESOLVERTIMER_HEADER_INCLUDED
#define OPM_FLEXIBLESOLVERTIMER_HEADER_INCLUDED

#endif // OPM_FLEXIBLESOLVERTIMER_HEADER_INCLUDED

#include <dune/istl/matrixmarket.hh>

template<class FS, class Vec, class C, class O, class SP>
void mockSolver(FS& fs, Vec rhs, const C& cc, O o, SP sp)
{
    Vec x(rhs.size());
    x=0;

    Vec& r = rhs;
    Vec p(x);
    Vec v(x);
    Vec t(x);
    Vec y(x);
    Vec y2(x);
    Vec y3(x);
    Vec rt(x);

    o.applyscaleadd(-1,x,r);

    rt = r;
    double norm = sp.norm(r);

    std::cout << "first norm "<< norm <<std::endl;
    //first iter

    double rho_new = sp.dot(rt,r);
    p = r;
    
    
    y = 0;
    fs.preconditioner().apply(y,p);

    y2=0;
    y3=0;
    fs.preconditioner().apply(y2,y3);
    std::cout << "prec y2 norm "<< sp.norm(y2) <<std::endl;
    bool first = true;
    for (int i = 0; i < y.size(); ++i) {
	for (int k = 0; k<3; ++k) {

	    
	    if (std::isnan(y[i][k])) {
		if (first) {
		    std::cout << "y idx: "<< i<<" y val "<< y[i][k] << " " << p[i][k]<<std::endl;
		    first = false;
		}
	    }
	}
    }
    std::cout << "prec y norm "<< sp.norm(y) <<std::endl;
    
    o.apply(y,v);
    std::cout << "mat v norm "<< sp.norm(v) <<std::endl;
    double h = sp.dot(rt,v);
    double alpha = rho_new / h;

    x.axpy(alpha,y);
    r.axpy(-alpha,v);
    norm = sp.norm(r);
    std::cout << "second norm "<< norm <<std::endl;
}

template<class FS, class Vec, class C>
void flexibleSolverAndPreconTimer(FS& fs, Vec rhs, const C& cc, std::string pc, double tol, bool update=false,
				  bool applyMeas=true)
{
    int rank = cc.rank();
    Vec x(rhs.size());
    Vec crhs(rhs);

    if (rank == 0) {std::cout << pc.c_str() << std::endl;}
    Dune::InverseOperatorResult stat;
    fs.apply(x, crhs, tol, stat);
    if (rank == 0) {std::cout << pc.c_str() << " result: iter: " << stat.iterations <<" elapsed " << stat.elapsed << std::endl;}

    if (applyMeas) {
	multipleMinLoopTimeFS(cc, fs, x, pc, 5);

    }
    if (update) {
	Dune::Timer timer;

	cc.barrier();
	timer.start();
	fs.preconditioner().update();
	cc.barrier();
	double utime = timer.stop();
	if (rank == 0) {std::cout << "Update "<< pc.c_str() << " time "<< utime << std::endl;}
    }
    if (rank == 0) {std::cout << std::endl;}
}

template<class S, class P, class Vec, class C>
void solverAndPreconTimer(S& s, P& p, Vec rhs, const C& cc, std::string pc, double tol, bool update=false)
{
    int rank = cc.rank();
    Vec x(rhs.size());
    Vec crhs(rhs);

    if (rank == 0) {std::cout << pc.c_str() << std::endl;}
    Dune::InverseOperatorResult stat;
    s.apply(x, crhs, tol, stat);
    if (rank == 0) {std::cout << pc.c_str() << " result: iter: " << stat.iterations <<" elapsed " << stat.elapsed << std::endl;}

    multipleMinLoopTimePre(cc, p, x, pc, 5);
    if (rank == 0) {std::cout << std::endl;}

    if (update) {
	Dune::Timer timer;

	cc.barrier();
	timer.start();
	p.update();
	cc.barrier();
	double utime = timer.stop();
	if (rank == 0) {std::cout << "Update "<< pc.c_str() << " time "<< utime << std::endl;}
    }
    
}

template<class FS, class FS2, class Vec, class C, class SP>
double solveTwice(FS& fs, FS2& fs2,Vec rhs, const C& cc, const SP& sp, std::string pc, double tol)
{
    int rank = cc.rank();
    Vec x1(rhs.size());
    Vec x2(rhs.size());
    Vec x3(rhs.size());
    Vec crhs1(rhs);
    Vec crhs2(rhs);
    Vec crhs3(rhs);
    Dune::InverseOperatorResult stat1;
    Dune::InverseOperatorResult stat2;
    fs.apply(x1, crhs1, tol, stat1);
    fs2.apply(x2, crhs2, tol, stat2);

    auto diff = x1;
    diff.axpy(-1,x2);

    double norm = sp.norm(diff);
    
    if (rank == 0) {std::cout << pc.c_str() << " norm diff "<< norm << std::endl;}

    fs.preconditioner().update();
    Dune::InverseOperatorResult stat3;
    fs.apply(x3, crhs3, tol, stat3);

    auto diff2 = x1;
    diff2.axpy(-1,x3);
    double norm2 = sp.norm(diff2);
    if (rank == 0) {std::cout << pc.c_str() << " update norm diff "<< norm2 << std::endl;}
    
    return norm;
}

template<class FS, class FS2, class Vec, class C, class SP>
void applyTwice(FS& fs, FS2& fs2,Vec rhs, const C& cc, const SP& sp, std::string pc)
{
    int rank = cc.rank();
    Vec y1(rhs.size());
    Vec y2(rhs.size());
    Vec x1(rhs.size());
    Vec x2(rhs.size());

    y1 = 0;
    y2 = 0;
    x1 = 0.0001;
    x2 = 0.0001;

    fs.preconditioner().pre(y1,x1);
    fs2.preconditioner().pre(y2,x2);
    
    fs.preconditioner().apply(y1,x1);
    fs2.preconditioner().apply(y2,x2);

    int numDiff = 0;
    for (int i = 0; i < y1.size(); ++i) {
	for (int j = 0; j<3; ++j) {
	    if (y1[i][j] != y2[i][j]) {

		double rel = (y1[i][j]-y2[i][j])/y1[i][j];
		numDiff++;
		//std::cout << rank <<" " << i<< " " << j <<" "<< y1[i][j] << " " <<y2[i][j]<< " "<<rel<<std::endl;
	    }
	}
    }

    std::cout << rank <<" " << numDiff << " "<< 3*y1.size()<< " "<< (double)numDiff/(3*y1.size()) <<std::endl;
    
    double normy = sp.norm(y1);
    auto diff = y1;
    diff.axpy(-1, y2);

    
    double norm = sp.norm(diff);

    if (rank == 0) {std::cout << pc.c_str() << " ynorm "<< normy <<" norm diff "<< norm << std::endl;}
    
}


template<class B>
void printBlock(B b)
{
    const int bs = B::rows;
    for (size_t i = 0; i < bs; ++i) {
	for (size_t k = 0; k < bs; ++k) {
	    std::cout << b[i][k]<<" ";
	}
	std::cout <<std::endl;
    }
}

template<class M>
void mock_ghost_last_bilu0_decomposition (M& A, size_t interiorSize)
{
    // iterator types
    typedef typename M::RowIterator rowiterator;
    typedef typename M::ColIterator coliterator;
    typedef typename M::block_type block;
    
    // implement left looking variant with stored inverse
    for (rowiterator i = A.begin(); i.index() < interiorSize; ++i)
    {
	// coliterator is diagonal after the following loop
	coliterator endij=(*i).end();           // end of row i
	coliterator ij;

	// eliminate entries left of diagonal; store L factor
	for (ij=(*i).begin(); ij.index()<i.index(); ++ij)
        {
	    // find A_jj which eliminates A_ij
	    coliterator jj = A[ij.index()].find(ij.index());

	    if (i.index() == 604143) {
		std::cout << std::endl;
		std::cout << i.index()<< " "<< ij.index()<< std::endl;
		printBlock(*jj);
		printBlock(*ij);
		std::cout << std::endl;
	    }
	    // compute L_ij = A_jj^-1 * A_ij
	    (*ij).rightmultiply(*jj);

	    if (i.index() == 604143) {
		printBlock(*ij);
		std::cout << std::endl;
	    }
	    
	    // modify row
	    coliterator endjk=A[ij.index()].end();    // end of row j
	    coliterator jk=jj; ++jk;
	    coliterator ik=ij; ++ik;
	    while (ik!=endij && jk!=endjk)
		if (ik.index()==jk.index())
                {
		    
		    block B(*jk);
		    if (i.index() == 604143) {
			std::cout << "B " << i.index()<< " "<< jk.index()<< " "<< ij.index()<<std::endl;
			printBlock(B);
			//std::cout << std::endl;
		    }
		    B.leftmultiply(*ij);
		    if (i.index() == 604143) {
			std::cout << "B2 ik " << i.index()<< " "<< jk.index()<< std::endl;
			printBlock(B);
			printBlock(*ik);
			std::cout << std::endl;
		    }
		    *ik -= B;
		    if (i.index() == 604143) {
			printBlock(*ik);
			std::cout << std::endl;
		    }
		    ++ik; ++jk;
		}
		else
                {
		    if (ik.index()<jk.index())
			++ik;
		    else
			++jk;
		}
	}

	// invert pivot and store it in A
	if (ij.index()!=i.index())
	    DUNE_THROW(Dune::ISTLError,"diagonal entry missing");
	try {
	    if (ij.index() == 604143) {
		printBlock(*ij);
		(*ij)[1][1] = 1.91483e-05;
	    }
	    (*ij).invert();   // compute inverse of diagonal block
	}
	catch (Dune::FMatrixError & e) {
	    DUNE_THROW(Dune::ISTLError,"ILU failed to invert matrix block");
	}
    }
}

template<class Mat, class C>
void compMat(const Mat& A, const Mat& B, const C& cc, int level)
{
    int rank = cc.rank();

    typedef typename Mat::RowIterator rowiterator;
    typedef typename Mat::ColIterator coliterator;
    typedef typename Mat::block_type block;
    const int bs = block::rows;
    
    auto brow = B.begin();

    int numDiff = 0;
    for (auto arow = A.begin(); arow != A.end(); ++arow) {

	auto bcol = brow->begin();

	std::vector<int> aids;
	std::vector<int> bids;
	for (auto acol = (*arow).begin(); acol!=(*arow).end(); ++acol) {

	    block a = *acol;
	    block b = *bcol;

	    auto aidx = acol.index();
	    auto bidx = bcol.index();

	    aids.push_back(aidx);
	    bids.push_back(bidx);
	    //if (aidx!= bidx)
	    //std::cout << rank << " " << level << " "<< arow.index() << " "<< aidx << " "<< bidx <<std::endl;
	    
	    for (int r = 0; r < bs; ++r){
		for (int c = 0; c < bs; ++c) {

		    double diff = std::abs(a[r][c]-b[r][c]);

		    if (diff > 0) {
			//std::cout << rank << " " << level << " "<< arow.index()<<" "<< acol.index()<< " "<< " "<< r <<" "<< c <<" " << diff << std::endl;
			numDiff++;
 		    }
		    
		}
	    }

	
	    bcol++;
	}
	if (false) {//rank == 0) {
	    std::cout << "A Cols: "<< rank<< " "<< level <<" "<< arow.index()<<" ";
	    for (int k = 0; k < aids.size(); ++k)
		std::cout << aids[k] << " ";
	    std::cout <<std::endl;
	    std::cout << "B Cols: "<< rank<< " "<< level <<" "<< brow.index()<<" ";
	    for (int k = 0; k < bids.size(); ++k)
		std::cout << bids[k] << " ";
	    std::cout <<std::endl;
	}
	brow++;
    }

    std::cout << "Num diff vals " << rank <<" "<<level<<" "<< numDiff<<std::endl; 
}

template<class Mat>
void countMatVals(const Mat& mat)
{
    const auto block_size = Mat::block_type::rows;

    int numOff = 0;
    int numDiag = 0;
    
    int numOffNeg = 0;
    int numOffZero = 0;
    int numDiagNeg = 0;
    int numDiagZero = 0;

    for ( auto row = mat.begin(); row!=mat.end(); ++row) {

	auto rid = row.index();
	for (auto col = row->begin(); col != row->end(); ++col) {

	    auto cid = col.index();

	    auto a = *col;
	    for (int i = 0; i < block_size; ++i) {
		for (int j = 0; j < block_size; ++j) {
		    if (rid != cid) {
			numOff++;

			if (a[i][j] < 0)
			    numOffNeg++;
			else if (a[i][j] == 0)
			    numOffZero++;
		    }
		    else {

			if (i!=j) {
			    numOff++;
			    if (a[i][j] < 0)
				numOffNeg++;
			    else if (a[i][j] == 0)
				numOffZero++;
			}
			else {
			    numDiag++;
			    if (a[i][j] < 0)
				numDiagNeg++;
			    else if (a[i][j] == 0)
				numDiagZero++;
			}
		    }
		}
	    }
	}
    }

    int numOffPos = numOff-numOffNeg-numOffZero;
    //std::cout << "Negative off count " << numOff << " " << numOffNeg<< " " << numOffZero << " " << numOffPos
    //<< " " << numDiag << " " << numDiagNeg << " "<< numDiagZero << std::endl;
    
}

template<class AMG, class C, class Vec>
void amgHiaInfo(AMG& amg, Vec rhs_loc, const C& cc )
{

    typedef typename AMG::ParallelInformation Comm;
    typedef typename AMG::Operator::matrix_type Mat;
    typedef typename AMG::OperatorHierarchy::AggregatesMap::AggregateDescriptor Agg;
    typedef Opm::ParallelOverlappingILU0<Mat,Vec,Vec,Comm> ILU;

    typedef Dune::Amg::Transfer<Agg,Vec,Dune::Amg::SequentialInformation> Transfer;
    
    int rank = cc.rank();
    
    auto opHir1 = amg.operatorHirarchyList();
    auto mats1 = opHir1->matrices();
    auto parHir1 = opHir1->parallelInformation().finest();
    auto fineMat1 = mats1.finest();
    auto coarseMat = mats1.finest();
    auto aggs = opHir1->aggregatesMaps().begin();
    coarseMat++;
    auto cmat1 = mats1.coarsest();
    //auto redistInfoL = opHir1->redistributeInformation();

    int redistVecSizeTo = 0;
    int redistVecSizeFrom = 0;
    
    int level = 0;
    bool movedToProc0 = false;
    bool writeMatrix = true;

    for (;fineMat1!=cmat1;) {
	if (rank==0) {
	    std::cout <<std::endl;
	    std::cout << "AMGLevelHiarchy " << level<<std::endl;
	}
	//double fnorm1 = fineMat1->getmat().frobenius_norm();
	int numRow = fineMat1->getmat().N();

	if (!fineMat1.isRedistributed() && !movedToProc0) {
	    std::vector<int> ct1;
	    std::vector<int> rt(numRow, 0);
	    getIndexSetInfo(parHir1, cc, ct1, rt, true);
	    printComTabOnRoot(cc, ct1);

	    countMatVals(fineMat1->getmat());
	    ILU ilu(fineMat1->getmat(), *parHir1, 0.99, Opm::convertString2Milu(std::string("ILU")), false, false);
	    timeLevel(*fineMat1, rhs_loc, ilu, *parHir1, cc); //solveParallel/amgSetup.hpp

	    Vec fineVec(fineMat1->getmat().N());
	    Vec coarseVec(coarseMat->getmat().N());

	    fineVec=1.0;
	    coarseVec=0;
	    Dune::Timer transferTime;

	    cc.barrier();
	    transferTime.start();
	    Transfer::restrictVector(*(*aggs),coarseVec,fineVec,Dune::Amg::SequentialInformation());
	    cc.barrier();
	    double tTime= transferTime.stop();
	    if (rank==0) {
		std::cout <<std::endl;
		std::cout << "MoveToCoarse " << tTime <<std::endl;
	    }

	    //aggs++;
	    transferTime.reset();
	    cc.barrier();
	    transferTime.start();
	    Transfer::prolongateVector(*(*aggs),coarseVec,fineVec,1.6,Dune::Amg::SequentialInformation());
	    double tTime2= transferTime.stop();
	    if (rank==0) {
		std::cout <<std::endl;
		std::cout << "MoveToFine " << tTime2 <<std::endl;
	    }
	}
	else {
	    if (!movedToProc0) {
		redistVecSizeTo = fineMat1->getmat().N();
		if (fineMat1.isRedistributed()) {
		    auto rmat = fineMat1.getRedistributed().getmat();
		    redistVecSizeFrom = rmat.N();
		}
	    }
	    if (rank==0) {
		std::cout << "Redist on non-coarsest level "<< std::endl;
		movedToProc0 = true;
	    }
	}

	if (writeMatrix) {
	    if (rank == 0) {
		std::string fname = "coarseLevelMat_" + std::to_string(level);
		std::ofstream filename(fname.c_str());
		Dune::writeMatrixMarket(fineMat1->getmat(), filename);
	    }
	}
	
	parHir1++;
	fineMat1++;
	coarseMat++;
	aggs++;
	level++;
    }
    if (fineMat1 == cmat1) {

	auto cmat = fineMat1->getmat();
	if (redistVecSizeTo == 0)
	    redistVecSizeTo = fineMat1->getmat().N();
	auto rmat = cmat;
	if (cc.size()>1)
	    if (fineMat1.isRedistributed()) {
		rmat = fineMat1.getRedistributed().getmat();
		redistVecSizeFrom = rmat.N();
	    }
	auto rN = rmat.N();
	if (rN > 0) {

	    Vec re_rhs(rN);
	    for (size_t i = 0; i < rN; ++i) {re_rhs[i]=rhs_loc[i];}
	    //std::cout << rank<< " N: "<< cmat.N() <<" "<< rmat.N()<<std::endl;
	    //timeCoarseSolverOnRoot();
	    
	    Dune::Timer coarsetimer;

	    coarsetimer.start();
	    Dune::UMFPack<Mat> coarseSolver(rmat, 1);
	    double decompTime = coarsetimer.stop();
	    std::cout << "DecompTime "<< decompTime <<std::endl;
	    Dune::InverseOperatorResult res;
	    Vec x(rN);
	    x=0;

	    coarsetimer.reset();
	    coarsetimer.start();
	    coarseSolver.apply(x, re_rhs, res);
	    double coarseApplyTime = coarsetimer.stop();
	    std::cout << "CoarseSolverTime " << coarseApplyTime << std::endl;

	    if (writeMatrix) {
		if (rank == 0) {
		    std::string fname = "coarseMat";
		    std::ofstream filename(fname.c_str());
		    Dune::writeMatrixMarket(cmat, filename);
		}
	    }
	}
    }
    auto redistInfo = opHir1->redistributeInformation().begin();
    auto redistInfoEnd = opHir1->redistributeInformation().end();
    for (auto ri = redistInfo; ri != redistInfoEnd; ++ri) {
	if (ri->isSetup()) {
	    if (false)//rank==0)
		std::cout << "redist is setup"<< std::endl;

	    Vec to(redistVecSizeTo);
	    Vec from(redistVecSizeFrom);
	    //std::cout << "resdistBack "<< rank<<" "<<redistVecSizeTo<<" "<<redistVecSizeFrom<<std::endl;
	    Dune::Timer redistTimer;

	    cc.barrier();
	    redistTimer.start();
	    ri->redistributeBackward(to, from);
	    cc.barrier();
	    double rdt = redistTimer.stop();

	    if (rank==0) {
		std::cout << std::endl;
		std::cout << "RedistTime " << rdt << std::endl;
	    }
	}
	else {
	    if (false)//rank==0)
		std::cout << "redist is not setup"<< std::endl;
	}
    }
	
}

template<class O, class C, class Vec, class CT, class LTP>
void timeFineLevel(O& op, Vec rhs_loc, const C& cc, CT comm, LTP ltp)
{

    typedef typename O::matrix_type Mat;
    typedef typename O::communication_type Comm;
    typedef Opm::ParallelOverlappingILU0<Mat,Vec,Vec,Comm> ILU;
    
    //const Comm comm = op.comm();
    auto mat = op.getmat();
    
    int rank = cc.rank();

    if (rank==0) {
	std::cout <<std::endl;
	std::cout << "AMGLevelHiarchy -1" <<std::endl;
    }

    std::vector<int> ct1;
    int numRow = mat.N();
    std::vector<int> rt(numRow, 0);
    
    getIndexSetInfo(comm, cc, ct1, rt, true);
    printComTabOnRoot(cc, ct1);

    ILU ilu(mat, *comm, 0.99, Opm::convertString2Milu(std::string("ILU")), false, false);
    timeLevel(op, rhs_loc, ilu, *comm, cc);


    ltp.createCoarseLevelSystem(op);
    if (rank == 0) {std::cout << std::endl;}
    cc.barrier();
    multipleMinLoopTimeTLP(cc, ltp, rhs_loc, 3);

    if (rank == 0) {std::cout << std::endl;}
    cc.barrier();
    multipleMinLoopTimeTLPC(cc, ltp, rhs_loc, 3);
}

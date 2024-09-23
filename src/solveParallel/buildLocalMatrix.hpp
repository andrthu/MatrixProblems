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

#ifndef OPM_BUILDLOCALMATRIX_HEADER_INCLUDED
#define OPM_BUILDLOCALMATRIX_HEADER_INCLUDED

#endif // OPM_BUILDLOCALMATRIX_HEADER_INCLUDED


///Construct adjecency pattern of local matrix based on global adjecency pattern.
template<class Mat>
Dune::MatrixIndexSet getAdjecency(Mat& A, std::vector<int>& overlap, 
				  std::vector<int>& local2global, 
				  std::vector<int>& global2local, int rank, bool assembleGhost=false)
{
    Dune::MatrixIndexSet op;
    op.resize(local2global.size(), local2global.size());
    
    for (auto row = A.begin(); row != A.end(); ++row)
    {
	int d = row.index();
	if (overlap[d] > -1)
	{
	    op.add(global2local[d], global2local[d]);
	    auto col = row->begin();
	    if (overlap[d] == 0 || assembleGhost)
	    {
		for (; col != row->end(); ++col)
		{
		    int nab = col.index();		    
		    if (overlap[nab] > -1)
		    {
			op.add(global2local[d], global2local[nab]);
		    }
		}	
	    }
	}
    }
    return op;
}

///Build local matrix fromm global matrix.
template<class Mat> 
void buildLocalMatrix(Mat& A, Mat& loc, std::vector<int>& overlap, 
		      std::vector<int>& local2global, 
		      std::vector<int>& global2local, int rank, bool assembleGhost=false)
{
    auto op = getAdjecency(A, overlap, local2global, global2local, rank, assembleGhost);
    op.exportIdx(loc);
    loc = 1;
    
    typename Mat::block_type diag(0.0);
    for (int el = 0; el < diag.size(); el++)
	diag[el][el] = 1.0e100;

    for (auto row = loc.begin(); row != loc.end(); ++row)
    {
	int d = row.index();
	
	int gd = local2global[d];
	auto block = A[gd][gd];
	if (overlap[gd] == 1)
	{
	    block = diag;
	}
	
	loc[d][d] = block;
	
	auto col = row->begin();
	for (; col != row->end(); ++col)
	{
	    int nab = col.index();
	    int gnab = local2global[nab];
	    if (overlap[gd] == 0)
		loc[d][nab] = A[gd][gnab];
	}	
    }
}

/// Build local vector based on global vector.
template<class Vec>
void buildLocalVector(Vec& v, Vec& loc, std::vector<int>& overlap, 
		      std::vector<int>& local2global, int rank)
{
    loc.resize(local2global.size());
    loc = 0;
    for (auto el = loc.begin(); el != loc.end(); ++el)
    {
	int idx = el.index();
	if (overlap[local2global[idx]] == 0)
	    loc[idx] = v[local2global[idx]];
    }
}

///Same functions with reordering 
///Construct adjecency pattern of local matrix based on global adjecency pattern.
template<class Mat>
Dune::MatrixIndexSet getAdjecencyReorder(Mat& A, std::vector<int>& overlap, 
					 std::vector<int>& local2global, 
					 std::vector<int>& global2local, int rank,
					 std::vector<int>& reorder, bool assembleGhost=false)
{
    Dune::MatrixIndexSet op;
    op.resize(local2global.size(), local2global.size());
    
    for (auto row = A.begin(); row != A.end(); ++row)
    {
	int d = row.index();
	if (overlap[d] > -1)
	{
	    int rd = reorder[global2local[d]];
	    op.add(rd, rd);
	    auto col = row->begin();
	    //std::cout << rd << " " << d << " " << overlap[d] << std::endl;
	    if (overlap[d] == 0 || assembleGhost)
	    {
		for (; col != row->end(); ++col)
		{
		    int nab = col.index();		
		    if (overlap[nab] > -1)		    
		    {
			op.add(rd, reorder[global2local[nab]]);
		    }
		}	
	    }
	}
    }
    return op;
}

///Build local matrix fromm global matrix.
template<class Mat> 
void buildLocalMatrixReorder(Mat& A, Mat& loc, std::vector<int>& overlap, 
			     std::vector<int>& local2global, 
			     std::vector<int>& global2local, int rank,
			     std::vector<int>& reorder, bool assembleGhost=false)
{        
    auto op = getAdjecencyReorder(A, overlap, local2global, global2local, rank, reorder, assembleGhost);    
    auto op_d = getAdjecency(A, overlap, local2global, global2local, rank, assembleGhost);
    
    op.exportIdx(loc);
    loc = 0;   
    
    Mat dummy;
    op_d.exportIdx(dummy);

    typename Mat::block_type diag(0.0);
    for (int el = 0; el < diag.size(); el++)
	diag[el][el] = 1.0e100;
        
    for (auto row = dummy.begin(); row != dummy.end(); ++row)
    {
	int d = row.index();
	int rd = reorder[d];
	int gd = local2global[d];
	auto block = A[gd][gd]; 
	if (overlap[gd] == 1)
	{
	    block = diag;
	}
	
	loc[rd][rd] = block;
	
	auto col = row->begin();
	for (; col != row->end(); ++col)
	{
	    int nab = col.index();
	    int rnab = reorder[nab];
	    int gnab = local2global[nab];
	    if (overlap[gd] == 0)
		loc[rd][rnab] = A[gd][gnab];
	}
    }
}

/// Build local vector based on global vector.
template<class Vec>
void buildLocalVectorReorder(Vec& v, Vec& loc, std::vector<int>& overlap, 
			     std::vector<int>& local2global, int rank,
			     std::vector<int>& reorder)
{
    loc.resize(local2global.size());
    loc = 0;
    for (auto el = loc.begin(); el != loc.end(); ++el)
    {
	int idx = el.index();
	if (overlap[local2global[idx]] == 0)
	    loc[reorder[idx]] = v[local2global[idx]];
    }
}

template<class Mat, class C>
void setGhostToUniform(Mat& loc, const std::vector<int>& rowType, const C& cc)
{
    for (auto row = loc.begin(); row != loc.end(); ++row) {

	int d = row.index();

	if (rowType[d] == 1) {
	    auto col = row->begin();

	    int numCol = 0;
	    for (; col != row->end(); ++col) {
		numCol+=1;
		//std::cout << loc[d][col.index()].determinant()<< " ";
	    }
	    //std::cout << cc.rank()<<" "<<d<< " "<< numCol<<std::endl;
	}	    
    }
}

template<class Mat>
Dune::MatrixIndexSet getAdjecencyFromLoc(Mat& A, std::vector<int>& overlap)
{

    Dune::MatrixIndexSet op;
    op.resize(A.N(), A.N());

    for (auto row = A.begin(); row != A.end(); ++row) {
	
	int d = row.index();

	if (overlap[d] == 0) {

	    auto col = row->begin();
	    for (; col != row->end(); ++col) {
		int nab = col.index();
		op.add(d, nab);
	    }
	}
	else {
	    op.add(d,d);
	}
    }
    return op;
}
template<class Mat, class Mat2>
void buildLocalMatrixFromLoc(Mat& A, Mat2& newLoc, std::vector<int>& overlap)
{
    auto op = getAdjecencyFromLoc(A, overlap);
    op.exportIdx(newLoc);

    for (auto row = A.begin(); row != A.end(); ++row) {

	int d = row.index();
	if (overlap[d] == 0) {
	    
	    auto col = row->begin();
	    for (; col != row->end(); ++col) {
		int nab = col.index();

		newLoc[d][nab] = A[d][nab];
	    }
	}
	else {
	    newLoc[d][d] = A[d][d];
	}
    }

    //std::cout << "Nonzeros for before and after adjecency fix "<< A.nonzeroes() << " " << newLoc.nonzeroes() << " " << A.nonzeroes()-newLoc.nonzeroes() << std::endl;
}

template<class Mat, class Mat2, class Comm>
void buildLocalMatrixFromLocAndComm(Mat& A, Mat2& newLoc, const Comm& comm)
{
    std::vector<int> overlap(A.N(), 0);

    auto indexSet = comm.indexSet();
    int numg = 0;
    for (auto idx = indexSet.begin(); idx!=indexSet.end(); ++idx) {

	if (idx->local().attribute()!=1) {
	    //std::cout << idx->local().local() << " ";
	    overlap[idx->local().local()] = 1;
	    numg += 1;
	}
	
    }
    //std::cout << std::endl;
    auto op = getAdjecencyFromLoc(A, overlap);
    op.exportIdx(newLoc);

    for (auto row = A.begin(); row != A.end(); ++row) {

	int d = row.index();
	if (overlap[d] == 0) {
	    
	    auto col = row->begin();
	    for (; col != row->end(); ++col) {
		int nab = col.index();

		newLoc[d][nab] = A[d][nab];
	    }
	}
	else {
	    newLoc[d][d] = A[d][d];
	}
    }

    //std::cout << "Nonzeros for before and after adjecency fix "<< A.nonzeroes() << " " << newLoc.nonzeroes() << " " << A.nonzeroes()-newLoc.nonzeroes() << " tot/ghost/inter"<< A.N()<< " "<< numg<< " " << A.N()-numg << std::endl;
}



template<class Mat, class Vec>
void invertAllBlocks(Mat A, Vec x)
{
    double min = 10;
    double max = 0;
    for (auto row = A.begin(); row != A.end(); ++row) {

	int d = row.index();
	auto col = row->begin();
	for (; col != row->end(); ++col) {
	    int nab = col.index();
	    
	    if (nab==d) {
		double det1 = A[d][nab].determinant ();
		A[d][nab].invert();
		double det =  A[d][nab].determinant();

		if (std::abs(det1)<min)
		    min = std::abs(det1);
		if (std::abs(det)>max)
		    max = std::abs(det);
		if (std::isnan(det))
		    std::cout<< d<<" "<<nab<<" "<<det1<<" "<<det<<std::endl;
		if (std::isnan(det1))
		    std::cout<< d<<" "<<nab<<" "<<det1<<" "<<det<<std::endl;
		if (std::isinf(det))
		    std::cout<< d<<" "<<nab<<" "<<det1<<" "<<det<<std::endl;
		if (std::isinf(det1))
		    std::cout<< d<<" "<<nab<<" "<<det1<<" "<<det<<std::endl;

		if (det==0)
		    std::cout<< d<<" "<<nab<<" "<<det1<<" "<<det<<std::endl;
		if (det1==0)
		    std::cout<< d<<" "<<nab<<" "<<det1<<" "<<det<<std::endl;
	    }
	}
    }

    std::cout << "Smallest det is "<< min<<" max is "<< max<<std::endl; 


    const int bs = Mat::block_type::rows;
    double vecMin = 10;
    for (size_t i = 0; i < x.size(); ++i) {

	for (size_t j = 0; j < bs; ++j) {
	    if (std::abs(x[i][j]) < vecMin) {

		vecMin = std::abs(x[i][j]);
	    }
	}
    }
    std::cout << "Smallest Vec is "<< vecMin<<std::endl;
}
template<class Mat>
void find_min_value(const Mat& A)
{
    const int bs = Mat::block_type::rows;

    bool firstnan = true;
    int firstnanid = -1;
    double small = 10;
    for (auto row = A.begin(); row != A.end(); ++row) {

	int d = row.index();
	auto col = row->begin();
	for (; col != row->end(); ++col) {

	    int nab = col.index();

	    if (d == 604143) {
		if (nab == 604143) {

		    for (size_t i = 0; i < bs; ++i) {
			for (size_t k = 0; k < bs; ++k) {
			    std::cout << A[d][nab][i][k]<<" ";
			}
			std::cout <<std::endl;
		    }
		}
	    }
	    for (size_t i = 0; i < bs; ++i) {
		for (size_t k = 0; k < bs; ++k) {

		    double val = std::abs(A[d][nab][i][k]);


		    
		    
		    if (std::isnan(val)) {
			if (firstnan) {
			    
			    std::cout << "nan id: "<< d<<" "<< nab <<" "<< i << " "<< k<<std::endl;
			    if (firstnanid<0)
				firstnanid = d;
			    else {
				if (firstnanid!=d)
				    firstnan = false;
			    }
			}
		    }
		    if (val>0) {
			if (val < small)
			    small = val;
		    }
		}
	    }
	}
    }

    std::cout << "Smallest mat val is "<< small <<std::endl;
}

template<class Mat, class Vec>
void findZeroDiag(const Mat& A, const Vec& v)
{

    const auto block_size = Vec::block_type::dimension;
    
    for (auto row = A.begin(); row != A.end(); ++row) {

	int d = row.index();
	auto col = row->begin();
	
	for (; col != row->end(); ++col) {
	    int nab = col.index();

	    if (d == nab) {

		for (int eq = 0; eq < block_size; eq++) {

		    if (A[d][nab][eq][eq] == 0) {
			std::cout << "Zero detexted: "<<d << " " << eq<< " " << A[d][nab][eq][eq] << std::endl;
			for (int er = 0; er < block_size; er++) {
			    for (int ec = 0; ec < block_size; ec++) {
				std::cout << A[d][nab][er][ec] << " ";
			    }
			    std::cout <<std::endl;
			}
			
		    }
		}
	    }
	}
    }
}

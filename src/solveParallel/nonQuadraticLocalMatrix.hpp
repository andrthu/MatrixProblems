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

#ifndef OPM_NONQUADRATICLOCALMATRIX_HEADER_INCLUDED
#define OPM_NONQUADRATICLOCALMATRIX_HEADER_INCLUDED

#endif // OPM_NONQUADRATICLOCALMATRIX_HEADER_INCLUDED

///Same functions with reordering 
///Construct adjecency pattern of local matrix based on global adjecency pattern.
template<class Mat>
Dune::MatrixIndexSet getAdjecencyRectangular(Mat& A, std::vector<int>& overlap, 
					     std::vector<int>& local2global, 
					     std::vector<int>& global2local, int rank,
					     int size, std::vector<int>& l2r, 
					     std::vector<int>& r2l)
{
    Dune::MatrixIndexSet op;
    op.resize(size, local2global.size());
    
    for (auto row = A.begin(); row != A.end(); ++row)
    {
	int d = row.index();
	
	if (overlap[d] == 0)
	{
	    int ld = global2local[d];
	    int rld = l2r[ld];
	    op.add(rld, ld);
	    auto col = row->begin();
	    if (overlap[d] == 0)
	    {
		for (; col != row->end(); ++col)
		{
		    int nab = col.index();
		    if (overlap[nab] >  -1)		    
		    {
			op.add(rld, global2local[nab]);
		    }
		}	
	    }
	}
    }
    return op;
}

template<class Mat>
Dune::MatrixIndexSet getAdjecencyRectangularRe(Mat& A, std::vector<int>& overlap, 
					       std::vector<int>& local2global, 
					       std::vector<int>& global2local,
					       int rank,
					       std::vector<int>& reorder, int size,
					       std::vector<int>& l2r, 
					       std::vector<int>& r2l)
{
    Dune::MatrixIndexSet op;
    op.resize(r2l.size(), local2global.size());
    
    for (auto row = A.begin(); row != A.end(); ++row)
    {
	int d = row.index();
	
	if (overlap[d] == 0)
	{
	    int rd = reorder[global2local[d]];
	    int rrd = l2r[rd];
	    op.add(rrd, rd);
	    auto col = row->begin();
	    if (overlap[d] == 0)
	    {
		for (; col != row->end(); ++col)
		{
		    int nab = col.index();
		    if (overlap[nab] >  -1)		    
		    {
			op.add(rrd, reorder[global2local[nab]]);
		    }
		}	
	    }
	}
    }
    return op;
}

///Build local matrix fromm global matrix.
template<class Mat> 
void buildLocalMatrixRectangular(Mat& A, Mat& loc, std::vector<int>& overlap, 
				 std::vector<int>& local2global, 
				 std::vector<int>& global2local, int rank,
				 std::vector<int>& reorder, int intSize,
				 std::vector<int>& l2r, std::vector<int>& r2l)
{
    /*    
    auto op = getAdjecencyRectangularRe(A, overlap, local2global, global2local, 
					rank, reorder, intSize, l2r, r2l);
    */
    auto op_d = getAdjecencyRectangular(A, overlap, local2global, 
					global2local, rank, intSize ,l2r, r2l);    
    op_d.exportIdx(loc);
    loc = 0;
    /*
    Mat dummy;
    op_d.exportIdx(dummy);
    */

    typename Mat::block_type diag(0.0);
    for (int el = 0; el < diag.size(); el++)
	diag[el][el] = 1.0e100;
        
    for (auto row = loc.begin(); row != loc.end(); ++row)
    {
	int d = row.index();
	//int rd = reorder[d];
	int ld = r2l[d];
	int gd = local2global[ld];
	auto block = A[gd][gd]; 
	if (overlap[gd] == 1)
	{
	    block = diag;
	}
	
	loc[d][ld] = block;
	
	auto col = row->begin();
	for (; col != row->end(); ++col)
	{
	    int nab = col.index();
	    //int rnab = reorder[nab];
	    int gnab = local2global[nab];
	    if (overlap[gd] == 0)
		loc[d][nab] = A[gd][gnab];
	}
    }
}

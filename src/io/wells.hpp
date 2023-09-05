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

#ifndef OPM_WELLS_HEADER_INCLUDED
#define OPM_WELLS_HEADER_INCLUDED

#endif // OPM_WELLS_HEADER_INCLUDED

#include <sstream>

Dune::MatrixIndexSet getOffAdj(std::vector<int> col, int N)
{
    Dune::MatrixIndexSet op;
    op.resize(1,N);    
    for (int i = 0; i < col.size(); ++i)
    {	
	op.add(0, col[i]);
    }
    return op;
}
template<class Mat>
void createSingleEntryMatrix(Mat& m)
{
    Dune::MatrixIndexSet op;
    op.resize(1, 1);    	
    op.add(0, 0);    
    op.exportIdx(m);
}

template<class Mat, class OffBlockList>
void createWellMatrix(Mat& m, OffBlockList& val, 
		      std::vector<int>& col, int N)
{
    auto op = getOffAdj(col, N);
    op.exportIdx(m);
    
    for (int i = 0; i < col.size(); ++i)
    {
	m[0][col[i]] = val[i];
    }
}

template<class MatList, class OffBlockListList>
void createOffMatrices(MatList& off, OffBlockListList& vals, 
		       std::vector<std::vector<int>>& cols, int N)
{
    typedef typename MatList::value_type Mat;
    auto val = vals.begin();    
    for (auto col = cols.begin(); col != cols.end(); ++col, ++val)
    {	
	Mat m;	
	createWellMatrix(m, *val, *col, N);
	off.push_back(m);	
    }
}

template<class Block, class Mat, class OffBlock, class WellVec>
class Wells
{
private:
    int numWells_;
    std::vector<Dune::BCRSMatrix<Block>> invD_;
    std::vector<Mat> B_;
    std::vector<Mat> C_;
    std::vector<std::vector<int>> wellCells_;
    std::vector<Dune::BlockVector<WellVec>> rhs_;
    
public:

    typedef Mat off_mat;
    
    Wells(){ numWells_ = 0; }

    void addWell(Dune::BCRSMatrix<Block> D, Mat C, Mat B, std::vector<int> wc,
		 Dune::BlockVector<WellVec> r)
    {
	invD_.push_back(D);
	C_.push_back(C);
	B_.push_back(B);
	wellCells_.push_back(wc);
	rhs_.push_back(r);
	numWells_++;
    }
    std::vector<Dune::BCRSMatrix<Block>> invD() {return invD_;}
    std::vector<Mat> B() {return B_;}
    std::vector<Mat> C() {return C_;}
    std::vector<std::vector<int>> wellCells() {return wellCells_;}
    std::vector<Dune::BlockVector<WellVec>> rhs() {return rhs_;}

    template<class Vec>
    void apply(const Vec& x, Vec& Ax) const
    {
	Dune::BlockVector<WellVec> Bx(1);
	Dune::BlockVector<WellVec> iDBx(1);
	for (int i = 0; i < numWells_; ++i)
	{	    
	    B_[i].mv(x, Bx);
	    invD_[i].mv(Bx, iDBx);
	    C_[i].mmtv(iDBx, Ax);
	}	
    }
    
    template<class Vec>
    void apply(Vec& rhs)
    {
	Dune::BlockVector<WellVec> iDr(1);
	for (int i = 0; i < numWells_; ++i)
	{
	    invD_[i].mv(rhs_[i],iDr);
	    C_[i].mmtv(iDr, rhs);
	}
    }

    void readWellsFromFiles(char* dirName, int N)
    {
	std::string dn = std::string(dirName);
	std::string B_name = dn + std::string("/offB.txt");
	std::string C_name = dn + std::string("/offC.txt");
	std::string D_name = dn + std::string("/InvD.txt");
	std::string R_name = dn + std::string("/rhsWell.txt");

	std::string line;
	int blockNum = 0;
	bool first_l = true;
	int blockDim = 0;
	int wellNum = 0;
	int blockCount = 0;
	
	//std::vector<std::vector<int>> Bcols;
	std::vector<std::vector<OffBlock>> Bvals;

	std::ifstream B_file(B_name);
	while (getline(B_file, line))
	{
	    if (!line.empty())
	    {
		if (line.compare(0,2,std::string("%%")) !=0 )
		{
		    if (line.compare(0,1,std::string("%")) !=0 )
		    {
			if (first_l)
			{
			    int col, row, nnz;
			    std::istringstream sline(line);
			    sline >> row;
			    sline >> col;
			    sline >> nnz;
			    
			    blockDim = col/3; 
			    first_l = false;			    
			}
			else
			{
			    int col, row;
			    double val;
			    std::istringstream sline(line);
			    sline >> row;
			    sline >> col;
			    sline >> val;
			    
			    if (blockCount % 12 == 0)
			    {
				Bvals[wellNum - 1].push_back(OffBlock(0.0));
				int blockCol = (col)/3;
				wellCells_[wellNum - 1].push_back(blockCol);

			    }
			    Bvals[wellNum - 1][blockCount/12][row - 1][blockCount % 3] = val;

			    blockCount++;
			}
		    }
		    else
		    {
			std::vector<OffBlock> new_well;
			std::vector<int> new_cols;
			Bvals.push_back(new_well);
			wellCells_.push_back(new_cols);
			
			first_l = true;
			wellNum++;
			blockCount = 0;
		    }
		}
	    }
	}
	
	B_file.close();
		
	createOffMatrices(B_, Bvals, wellCells_, N);
	
	//std::string line;
	blockNum = 0;
	first_l = true;
	blockDim = 0;
	wellNum = 0;
	blockCount = 0;
	
	std::vector<std::vector<int>> Ccols;
	std::vector<std::vector<OffBlock>> Cvals;

	std::ifstream C_file(C_name);
	while (getline(C_file, line))
	{
	    if (!line.empty())
	    {
		if (line.compare(0,2,std::string("%%")) !=0 )
		{
		    if (line.compare(0,1,std::string("%")) !=0 )
		    {
			if (first_l)
			{
			    int col, row, nnz;
			    std::istringstream sline(line);
			    sline >> row;
			    sline >> col;
			    sline >> nnz;
			    
			    blockDim = col/3; 
			    first_l = false;
			    
			}
			else
			{
			    int col, row;
			    double val;
			    std::istringstream sline(line);
			    sline >> row;
			    sline >> col;
			    sline >> val;
			    
			    if (blockCount % 12 == 0)
			    {
				Cvals[wellNum - 1].push_back(OffBlock(0.0));
				int blockCol = (col)/3;
				Ccols[wellNum - 1].push_back(blockCol);
			    }
			    Cvals[wellNum - 1][blockCount/12][row - 1][blockCount % 3] = val;

			    blockCount++;
			}
		    }
		    else
		    {
			std::vector<OffBlock> new_well;
			std::vector<int> new_cols;
			Cvals.push_back(new_well);
			Ccols.push_back(new_cols);
			first_l = true;
			wellNum++;
			blockCount = 0;
		    }
		}
	    }
	}
	C_file.close();
		
	createOffMatrices(C_, Cvals, Ccols, N);

	std::ifstream D_file(D_name);
	bool skip = true;
	while (getline(D_file, line))
	{
	    if (!line.empty())
	    {
		if (line.compare(0,2,std::string("%%")) !=0 )
		{
		    if (line.compare(0,1,std::string("%")) !=0 )
		    {
			if (skip)
			    skip = false;
			else
			{
			    int col, row;
			    double val;
			    std::istringstream sline(line);
			    sline >> row;
			    sline >> col;
			    sline >> val;			    
			    
			    invD_[blockNum - 1][0][0][row - 1][col - 1] = val;
			}
		    }
		    else
		    {
			blockNum++;
			Dune::BCRSMatrix<Block> inv_m;
			createSingleEntryMatrix(inv_m);
			inv_m[0][0] = Block(0.0);
			invD_.push_back(inv_m);
			skip = true;
		    }
		}
	    }
	}
	D_file.close();
	numWells_ = blockNum;
	
	blockNum = 0;
	int vecNum = 0;
	std::ifstream R_file(R_name);
	while (getline(R_file, line))
	{
	    if (!line.empty())
	    {
		if (line.compare(0,2,std::string("%%")) !=0 )
		{
		    if (line.compare(0,1,std::string("%")) !=0 )
		    {
			if (skip)
			    skip = false;
			else
			{
			    double val;
			    std::istringstream sline(line);			    
			    sline >> val;			    
			    
			    rhs_[blockNum - 1][0][vecNum] = val;
			    vecNum++;
			}
		    }
		    else
		    {
			Dune::BlockVector<WellVec> rhs_v(1);
			rhs_v = 0.0;
			rhs_.push_back(	Dune::BlockVector<WellVec>(rhs_v));
			skip = true;
			blockNum++;
			vecNum = 0;			
		    }
		}
	    }
	}
	R_file.close();
    }

    void printOff()
    {

	for (int i = 0; i < numWells_; ++i)
	{
	    auto c = C_[i][0].begin();
	    for (auto b = B_[i][0].begin(); b!=B_[i][0].end(); ++b, ++c)
	    {
		std::cout << b.index()<< " "<< c.index()<< " "<< i <<std::endl;
	    }
	    std::cout << std::endl;
	}
    }

    void printD()
    {
	for (int i = 0; i < numWells_; ++i)
	{
	    for (int k = 0; k < 4; ++k)
	    {
		for(int j = 0; j < 4; ++j)
		{
		    std::cout << invD_[i][0][0][k][j] << " ";
		}
		std::cout << std::endl;
	    }
	    std::cout << std::endl;
	}
	
    }

    void printRHS()
    {
	for (int i = 0; i < numWells_; ++i)
	{
	    for (int j = 0; j < 4; ++j)
	    {
		std::cout << rhs_[i][0][j] << std::endl; 
	    }
	    std::cout << std::endl;
	}
    }
    int numWells() {return numWells_;}
};

template<class Mat>
Dune::MatrixIndexSet getBCAdjecency(Mat& A, std::vector<int>& local2global, 
				    std::vector<int>& global2local, int rank)
{
    Dune::MatrixIndexSet op;
    op.resize(1, local2global.size());
    
    for (auto col = A[0].begin(); col != A[0].end(); ++col)
    {
	int idx = col.index();
	op.add( 0, global2local[idx] );
    }

    return op;    
}
template<class Mat>
void buildBCMatrix(Mat& A, Mat& loc, std::vector<int>& local2global, 
		   std::vector<int>& global2local, int rank)
{
    
    auto op = getBCAdjecency(A, local2global, global2local, rank);
    op.exportIdx(loc);
    loc = 0;
    
    for (auto row = loc.begin(); row != loc.end(); ++row)
    {
	int d = row.index();
	int gd = local2global[d];
	
	auto col = row->begin();
	for (; col != row->end(); ++col)
	{
	    int nab = col.index();
	    int gnab = local2global[nab];	    
	    loc[0][nab] = A[0][gnab];
	}
    }
}
template<class W>
void buildLocalWells(W& w, W& l, std::vector<int>& overlap,
		     std::vector<int>& local2global,
		     std::vector<int>& global2local, int rank)
{
    auto wc = w.wellCells();
    int numWells = w.numWells();
    for (int i = 0; i < numWells; ++i)
    {	
	if (overlap[wc[i][0]] == rank)
	{
	    typename W::off_mat c;
	    typename W::off_mat b;
	    
	    buildBCMatrix(w.C()[i], c, local2global, global2local, rank);
	    buildBCMatrix(w.B()[i], b, local2global, global2local, rank);
	    
	    l.addWell(w.invD()[i], c, b, wc[i], w.rhs()[i]);
	}
    }
}

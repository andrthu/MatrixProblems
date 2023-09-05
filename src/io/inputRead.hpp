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

#ifndef OPM_INPUTREAD_HEADER_INCLUDED
#define OPM_INPUTREAD_HEADER_INCLUDED

#endif // OPM_INPUTREAD_HEADER_INCLUDED

#include <boost/algorithm/string/predicate.hpp>
#include <boost/program_options.hpp>

template<class O>
void readMatMarketObject(O& o,const char* fileName)
{
    std::ifstream fileIn(fileName);
    Dune::readMatrixMarket(o, fileIn);
    fileIn.close();
}

template<class Mat3, class Mat1, class Vec>
void readFromDir(Mat3& A, Mat1& trans, Mat1& wells, Vec& rhs, char* dirName, int rank)
{
    std::string d = std::string(dirName);
    std::string A_name = d + std::string("/BlackoilMatrix.mtx");
    std::string t_name = d + std::string("/transAdj.mtx");
    std::string w_name = d + std::string("/wellAdj.mtx");
    std::string r_name = d + std::string("/BlackoilRHS.vec");

    if (rank == 0) {std::cout<<"Read Matrix"<<std::endl;}
    readMatMarketObject(A, A_name.data());
    if (rank == 0) {
	readMatMarketObject(trans, t_name.data());
	readMatMarketObject(wells, w_name.data());
    }
    if (rank == 0) {std::cout<<"Read RHS"<<std::endl;}
    readMatMarketObject(rhs, r_name.data());
}

template<class Mat3, class Mat1, class Vec, class D>
void handleMatrixSystemInput(int argc, char** argv, Mat3& A, Mat1& trans, 
			     Mat1& wells, Vec& rhs, D& DR, int rank)
{
    if (argc>3)
    {
	if (argc>4)
	    DR.read_file_and_update(argv[5]);
	if (rank==0)
	    DR.write_param();
	
	readMatMarketObject(A, argv[1]);
	readMatMarketObject(trans, argv[2]);
	readMatMarketObject(rhs, argv[3]);
	readMatMarketObject(wells, argv[4]);
    }
    else
    {
	if (argc>1)
	    DR.read_file_and_update(argv[2]);
	if (rank==0)
	    DR.write_param();
	
	readFromDir(A,trans,wells,rhs,argv[1], rank);
    }
}

template<class Mat3, class Mat1, class Vec, class D, class C>
void handleMatrixSystemInputSomeRanks(int argc, char** argv, Mat3& A, Mat1& trans, 
				      Mat1& wells, Vec& rhs, D& DR, const C& cc, bool readMat)
{
    int rank = cc.rank();
    
    if (argc>3)
    {
	if (argc>4)
	    DR.read_file_and_update(argv[5]);
	if (rank==0)
	    DR.write_param();
	if (readMat) {
	    readMatMarketObject(A, argv[1]);
	    readMatMarketObject(trans, argv[2]);
	    readMatMarketObject(rhs, argv[3]);
	    readMatMarketObject(wells, argv[4]);
	}
    }
    else
    {
	if (argc>2) {
	    if ( boost::algorithm::ends_with( argv[2], ".json") ) {
		if (rank==0)
		    std::cout << "json file: " << argv[2] << std::endl;
		DR.dict[12] = std::string(argv[2]);
	    } else {
		DR.read_file_and_update(argv[2]);
	    }
	} 
	if (rank==0)
	    DR.write_param();
	if (readMat) {
	    readFromDir(A,trans,wells,rhs,argv[1], rank);
	}
    }
    //cc.barrier();
}

bool help_detect(int argc, char** argv)
{
    namespace po = boost::program_options;

    po::options_description desc("Program solving linear systems generated from reservoir simulations. Runs in parallel and uses transmissibilty for partitioning. jsonSolve executable solves 3-phase black-oil systems and jsonSolveCO2 solves 2-phase CO2-store systems. \nHow to run: \n .\\jsonSolve matrix/dir/path [optional]jsonFile\n\nMatrix directory must contain four files: \n  *BlackoilMatrix.mtx: Matrix file.\n  *BlackoilRHS.vec: Vector file. \n  *transAdj.mtx: Transmissibility file for partitioning. \n  *wellAdj.mtx: Well connections for partitioning. \n\nBlack-oil system matrices can be found in the directory: /global/D1/homes/andreast/linear_systems/3-phase/ \nCO2 stsyems in: /global/D1/homes/andreast/linear_systems/2-phase/ \n\n Options");
    desc.add_options()
	("help", "produce help message");
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help")) {
	std::cout << desc << "\n";
	return true;
    }
    return false;
    
}

template<class MatB3, class MatB1>
void printMatrixAndInfo(MatB3& A, MatB1 trans)
{
    double nnz = 1.0;
    int max_nnz = 0;
    int min_nnz = 1000;
    
    for (auto row = A.begin(); row != A.end(); ++row)
    {
	auto col = row->begin();
	int row_nnz =0;
	for (; col != row->end(); ++col)
	{
	    std::cout << col.index() << " ";
	    row_nnz++;
	    nnz+=1.0;
	}
	if (row_nnz>max_nnz)
	    max_nnz=row_nnz;
	if (row_nnz<min_nnz)
	    min_nnz=row_nnz;

	std::cout<<std::endl;	 
    }
    for (auto row = trans.begin(); row != trans.end(); ++row)
    {
	auto col = row->begin();
	
	for (; col != row->end(); ++col)
	{
	    std::cout << *col << " ";
	}	
	std::cout<<std::endl;
    }

    std::cout << "Rows: "<<A.N()<<" Columns: "<<A.M()<<std::endl;
    double avg_nnz= nnz/A.N();
    std::cout << "avg: "<< avg_nnz<<" min: " <<min_nnz<< " max: "<<max_nnz << std::endl;
}


template<class MatB3>
void writeMatrixRowSizeHist(MatB3& A)
{
    std::ofstream histFile;
    histFile.open("NorneHist.txt");

    for (auto row = A.begin(); row != A.end(); ++row)
    {
	auto col = row->begin();
	int row_nnz = 0;
	for (; col != row->end(); ++col)
	{	   
	    row_nnz++;
	}
	histFile << row_nnz;
	histFile << "\n";
	//std::cout << row_nnz;
	//std::cout << std::endl;	 
    }
    histFile.close();
}

template<class Mat>
void storeRowSize(Mat& A, std::vector<int>& row_size)
{
    for (auto row = A.begin(); row != A.end(); ++row)
    {
	auto col = row->begin();
	int row_nnz = 0;
	for (; col != row->end(); ++col)
	{	   
	    row_nnz++;
	}
	row_size[row.index()] = row_nnz;
    }
}

template<class Mat, class C>
void storeRowSizeFromRoot(Mat& A, std::vector<int>& row_size, const C& cc)
{
    int rank = cc.rank();
    if (rank == 0) {
	for (auto row = A.begin(); row != A.end(); ++row)
	{
	    auto col = row->begin();
	    int row_nnz = 0;
	    for (; col != row->end(); ++col)
	    {
		row_nnz++;
	    }
	    row_size[row.index()] = row_nnz;
	}
    }
    cc.broadcast(row_size.data(), row_size.size(), 0);
}

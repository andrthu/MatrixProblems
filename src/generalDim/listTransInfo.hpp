/*
  Copyright 2025 Andreas Thune.

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
#include <algorithm>

template<class Mat>
double infoNNZ(Mat& t)
{
    std::vector<double> vals;
    double max = 0;
    double s = 0;
    for (auto row=t.begin();row!=t.end();++row) {
	auto col = row->begin();
	for (; col!=row->end(); ++col) {
	    auto val = *col;
	    if (val > max)
		max = val;
	    vals.push_back(val);
	    s += val;
	}
    }

    std::cout << max << " "<< vals.size()<<" "<<s/vals.size() <<std::endl;

    int num90 = 0;
    int num50 = 0;
    int num10 = 0;
    int num1 = 0;
    int num01 = 0;
    int num001 = 0;
    int num0 = 0;
    int numA = 0;
    for (int i = 0; i < vals.size(); ++i) {

	auto v = vals[i];
	if (v>0.9*max) {
	    num90++;
	    
	}
	if (v>0.5*max) {
	    num50++;
	    
	}
	if (v>0.1*max) {
	    num10++;
	    
	}
	if (v>0.01*max) {
	    num1++;
	    
	}
	if (v>0.001*max) {
	    num01++;
	    
	}
	if (v>0.0001*max) {
	    num001++;
	}
	if (v == 0)
	    num0++;
	if (v>s/vals.size())
	    numA++;
    }

    std::cout << num90 << " "<< num50 <<" "<< num10 <<" "<< num1 <<" "<< num01<<" "<<num001 <<
	" "<< num0<< " " << numA << std::endl;

    return max;
}

template<class Mat>
std::vector<double> infoNNZ2(Mat& t)
{
    std::vector<double> vals;
    double max = 0;
    double s = 0;
    for (auto row=t.begin();row!=t.end();++row) {
	auto col = row->begin();
	for (; col!=row->end(); ++col) {
	    auto val = *col;
	    if (val > max)
		max = val;
	    vals.push_back(val);
	    s += val;
	}
    }

    std::sort(vals.begin(),vals.end());

    return vals;
}

template<class Mat, class T>
void removeSmallTransNNZ(const Mat& A, Mat& M, T trans, double w)
{
    Dune::MatrixIndexSet op;
    op.resize( A.N(), A.N() );

    

    for (auto row = A.begin(); row != A.end(); ++row) {

	int d = row.index();
	op.add(d,d);
	auto col = row->begin();
	for (; col != row->end(); ++col) {
	    int nab = col.index();
	    if (trans.exists(d,nab)) {
		if (trans[d][nab] > w)
		    op.add(d,nab);
	    } else {
		op.add(d,nab);
	    }
	}
	
    }

    op.exportIdx(M);
    

    for (auto row = M.begin(); row != M.end(); ++row) {
	int d = row.index();
	M[d][d] = A[d][d];
	auto col = row->begin();
	for (; col != row->end(); ++col) {
	    int nab = col.index();
	    M[d][nab] = A[d][nab];
	}
    }
}


template<class Mat3, class Vec>
void gen_dim_list_trans_info(int argc, char** argv)
{

  auto systemDirs = parse_multiple_systems(argc, argv);
    

  typedef Dune::FieldMatrix<double,1,1> BlockMat1;
  typedef Dune::BCRSMatrix<BlockMat1> Mat;
    
  typedef Dune::MPIHelper::MPICommunicator MPICommunicator;
  typedef Dune::Communication<MPICommunicator> CollectiveCommunication;    

  typedef Dune::BiCGSTABSolver<Vec> Solver;
  typedef Dune::InverseOperatorResult Stat;
    
  typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
  typedef Dune::OverlappingSchwarzScalarProduct<Vec,Comm> ScalarProduct;
  typedef GhostLastMatrixAdapter<Mat3,Vec,Vec,Comm> GLO;                 // solveParallel/ghostLastOperations.hpp
  typedef Dune::OverlappingSchwarzOperator<Mat3,Vec,Vec,Comm> Operator;
  typedef Opm::ParallelOverlappingILU0<Mat3,Vec,Vec,Comm> ILU;
  
  const auto block_size = Vec::block_type::dimension;
    
  CollectiveCommunication cc(MPI_COMM_WORLD);
  int rank = cc.rank();

    DictRead DR;
    Comm comm(cc);
    std::shared_ptr<Comm> parComm(new(Comm));
    std::vector<int> mpiVec;
    
    std::vector<Mat> mats;
    std::vector<Mat3> systems;
    std::vector<Vec> rhs;
    std::vector<ScalarProduct> sps;
    
    for (int i = 0; i < systemDirs.size(); ++i) {

	if ( boost::algorithm::ends_with( systemDirs[i], ".json") ) {
	    DR.dict[12] = systemDirs[i];
	} else if ( boost::algorithm::ends_with( systemDirs[i], ".ini") ) {
	    DR.read_file_and_update(systemDirs[i].data());
	}    
	else {
	    Mat trans, wells;

	    readTransMatOnly(trans, systemDirs[i], rank);
	    readWellMatOnly(wells, systemDirs[i], rank);

	    std::vector<int> part;
	    std::vector<int> rs(trans.N()); 
	    storeRowSizeFromRoot(trans, rs, cc);
	    part.resize(trans.N(), rank);

	    zoltanPartitionFunction(part, trans, wells, cc, DR, rs, 10);

	    mats.push_back(trans);

	    Mat3 A_loc;
	    Vec rhs_loc;
	    mpiVec = readMatOnRootAndDist(systemDirs[i], A_loc, rhs_loc, DR, comm, parComm, cc, mpiVec, true, i!=0);
	    systems.push_back(A_loc);
	    rhs.push_back(rhs_loc);
	    sps.push_back(ScalarProduct(*parComm));
	}
    }

    for (int i = 0; i < mats.size(); ++i) {

	auto t = mats[i];
	std::cout << t.N() << " "<<t.nonzeroes()  <<std::endl;
	double max = infoNNZ(t);
	std::vector<double> sortTrans = infoNNZ2(t);
	Mat3 M;
	auto A_ = systems[i];
	GLO op(A_, *parComm);
	auto sp = sps[i];
	auto rhs_ = rhs[i];

	std::vector<double> W = {-1, 0.00001, 0.00002, 0.00005, 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.05, 0.1};
	std::vector<double> ratio = {-1, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.8};

	if (block_size == 2) {
	    W = {-1, 0.000001, 0.000002, 0.000005, 0.00001, 0.00002, 0.00005, 0.0001, 0.0002, 0.0005, 0.001, 0.005};

	} else {

	    if (t.N() > 74431){

		W = {-1, 0.00001, 0.00002, 0.00003, 0.00004, 0.00005, 0.0001, 0.0002, 0.0005, 0.001, 0.01};
	    }
	}

	for (double w : ratio) {

	    double threshold = -1;
	    if (w > 0)
		threshold = sortTrans[(int) (w*sortTrans.size())];
	    //removeSmallTransNNZ(A_, M, t, max, w);
	    removeSmallTransNNZ(A_, M, t, threshold);

	    double rr = ((double)M.nonzeroes())/A_.nonzeroes();
	    double rr2 = ((double)M.nonzeroes() - M.N() )/(A_.nonzeroes()-M.N());
	    //std::cout << M.N() << " " << M.nonzeroes() << " " << A_.nonzeroes() << " " << rr <<std::endl;

	    std::string use_ilu("ILU");
	    auto ilu_help = Opm::convertString2Milu(use_ilu);

	    double tol = std::stod(DR.dict[1]);
	    
	    ILU ilu(M, *parComm, 1, ilu_help, M.N(), false, false );
	    Solver bicg(op, sp, ilu, tol, 2000, 0);
	    Stat statistics;

	    Vec x(A_.N());
	    x = 0;
	    Vec rhsC(rhs_);
	    bicg.apply(x, rhsC, statistics);

	    std::cout << "Solve "<< w<< " " << statistics.iterations << " " << statistics.elapsed << " "
		      << statistics.elapsed/statistics.iterations<< " "<< rr << " "<< rr2 << std::endl;
	}
    }

}

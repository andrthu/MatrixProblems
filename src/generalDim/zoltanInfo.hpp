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


template<class Mat, class T>
void removeNonPart(const Mat& A, Mat& M, T trans, std::vector<int> part, std::vector<double> inter)
{
    Dune::MatrixIndexSet op;
    op.resize( A.N(), A.N() );

    int numOff = 0;
    int numOff2 = 0;
    std::vector<int> counts(inter.size() - 1, 0);
    
    for (auto row = A.begin(); row != A.end(); ++row) {

	int d = row.index();

	int dp = part[d];
	op.add(d,d);
	auto col = row->begin();
	for (; col != row->end(); ++col) {
	    int nab = col.index();
	    int np = part[nab];
	    if (dp == np) {
		op.add(d,nab);
	    }
	    else {
		
		numOff++;
		if (trans.exists(d,nab)) {
		    numOff2++;
		    double num = trans[d][nab];
		    for (size_t i = 0; i < inter.size() - 1; ++i) {
			if (num >= inter[i] && num < inter[i+1]) {
			    counts[i]++;
			}
		    }
		}
	    }
	}
    }

    std::cout << std::endl;
    std::cout << numOff<<" " << numOff2<< std::endl;
    for (size_t i = 0; i < inter.size() - 1; ++i) {
	std::cout << ((double) counts[i])/numOff << " ";
    }
    
    std::cout << std::endl;

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

template<class Mat>
std::vector<double> sortTrans(Mat t, std::vector<double> & st)
{
    for (auto row = t.begin(); row != t.end(); ++row) {

	auto col = row->begin();
	for (; col != row->end(); ++col) {
	    st.push_back(*col);
	}
    }

    std::sort(st.begin(), st.end());

    std::cout << st[0] << " " << st[st.size()/10] << " "<< st[st.size()/2] << " " << st[9*st.size()/10] << " " << st[st.size() -1 ] << std::endl;

    auto N = st.size();
    std::vector<double> p = {st[0], st[N/10], st[N/5], st[3*N/10], st[4*N/10], st[N/2], st[6*N/10], st[7*N/10], st[8*N/10], st[9*N/10], st[N-1]};

    return p;
    
}


template<class Mat3, class Vec>
void gen_dim_zoltan_info(int argc, char** argv)
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
    std::vector<std::vector<std::vector<int>>> sysParts;
    std::vector<Vec> rhs;
    std::vector<ScalarProduct> sps;
    std::vector<std::vector<double>> inters;


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
	    
	    std::vector<std::vector<int>> parts;
	    std::vector<int> rs(trans.N()); 
	    storeRowSizeFromRoot(trans, rs, cc);

	    std::vector<double> st;
	    auto inter = sortTrans(trans,st);
	    inters.push_back(inter);

	    std::vector<int> numP = {1, 2, 4, 8 , 16, 32, 64, 128, 256, 512};
	    for (int np : numP ) {
		std::vector<int> part;
		part.resize(trans.N(), rank);

		std::cout << "Part for " << np << " parrtitions" <<std::endl;
		if ( std::stoi(DR.dict[13]) > 0) {
		    //zoltanPartitionFunction(part, trans, wells, cc, DR, rs, np, inter[std::min(teller,9)]);
		    int precent =std::stoi(DR.dict[14]);
		    zoltanPartitionFunction(part, trans, wells, cc, DR, rs, np, inter[precent]);
		} else {
		    zoltanPartitionFunction(part, trans, wells, cc, DR, rs, np, -1);
		}

		
		serialPartEval(trans,wells,part,np);
		parts.push_back(part);
	    }

	    sysParts.push_back(parts);
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
	std::vector<double> st;
	auto inter = inters[i];//sortTrans(t,st);
	
	
	auto parts = sysParts[i];
	std::cout << t.N() << " "<<t.nonzeroes()  <<std::endl;
	double max = 0;
	Mat3 M;
	auto A_ = systems[i];
	GLO op(A_, *parComm);
	auto sp = sps[i];
	auto rhs_ = rhs[i];
	


	for (std::vector<int> part : parts) {
	    removeNonPart(A_, M, t, part, inter);

	    double rr = ((double)M.nonzeroes())/A_.nonzeroes();
	    double rr2 = ((double)M.nonzeroes() - M.N() )/(A_.nonzeroes()-M.N());
	    //std::cout << M.N() << " "<<M.nonzeroes()<< " " << A_.nonzeroes() << " "<< rr <<std::endl;

	    std::string use_ilu("ILU");
	    auto ilu_help = Opm::convertString2Milu(use_ilu);

	    double tol = std::stod(DR.dict[1]); //0.005;
	    if (block_size == 2) {
		tol = std::stod(DR.dict[1]); //0.005;
	    }
	    
	    ILU ilu(M, *parComm, 1, ilu_help, M.N(), false, false );
	    Solver bicg(op, sp, ilu, tol, 200, 0);
	    Stat statistics;

	    Vec x(A_.N());
	    x = 0;
	    Vec rhsC(rhs_);
	    bicg.apply(x, rhsC, statistics);

	    std::cout << "Solve "<< 0 << " " << statistics.iterations << " " << statistics.elapsed << " "
	    << statistics.elapsed/statistics.iterations<< " "<< rr << " "<< rr2 << std::endl;
	}
    }

}

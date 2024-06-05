/*
  Copyright 2024 Andreas Thune.

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

#ifndef OPM_SINGLEHIERARCHYMULTJSONSOLVEGENERAL_HEADER_INCLUDED
#define OPM_SINGLEHIERARCHYMULTJSONSOLVEGENERAL_HEADER_INCLUDED

#endif // OPM_SINGLEHIERARCHYMULTJSONSOLVEGENERAL_HEADER_INCLUDED


template<class Mat, class Vec>
void gen_dim_jsonSolve_mult_sys_same_hir(std::vector<std::string> systemDirs)
{
    typedef Dune::MPIHelper::MPICommunicator MPICommunicator;
    typedef Dune::CollectiveCommunication<MPICommunicator> CollectiveCommunication;
    typedef Dune::BiCGSTABSolver<Vec> Solver;
    typedef Dune::InverseOperatorResult Stat;
    
    typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
    typedef Dune::OverlappingSchwarzScalarProduct<Vec,Comm> ScalarProduct;
    typedef GhostLastMatrixAdapter<Mat,Vec,Vec,Comm> GLO;                 // solveParallel/ghostLastOperations.hpp
    typedef Dune::OverlappingSchwarzOperator<Mat,Vec,Vec,Comm> Operator;
    typedef Opm::ParallelOverlappingILU0<Mat,Vec,Vec,Comm> ILU;
    typedef Dune::FlexibleSolver<GLO> FlexibleSolverType;

    const auto block_size = Vec::block_type::dimension;
    
    CollectiveCommunication cc(MPI_COMM_WORLD);
    int rank = cc.rank();

    std::vector<Mat> systems;
    std::vector<Vec> rhs;

    std::vector<ScalarProduct> sps;
    
    DictRead DR;
    Comm comm(cc);
    std::shared_ptr<Comm> parComm(new(Comm));
    std::vector<int> mpiVec;
    for (int i = 0; i < systemDirs.size(); ++i) {
	
	if ( boost::algorithm::ends_with( systemDirs[i], ".json") ) {
	    DR.dict[12] = systemDirs[i];
	}
	else {
	    Mat A_loc;
	    Vec rhs_loc;

	    mpiVec = readMatOnRootAndDist(systemDirs[i], A_loc, rhs_loc, DR, comm, parComm, cc, mpiVec, true, i!=0);

	    systems.push_back(A_loc);
	    rhs.push_back(rhs_loc);
	    sps.push_back(ScalarProduct(*parComm));
	}
    }

    if (rank == 0) {std::cout << std::setprecision (15) << std::endl;}

    Opm::FlowLinearSolverParameters flsp_json;
    flsp_json.linsolver_ = DR.dict[12];

    Opm::PropertyTree prm_json(flsp_json.linsolver_);
    std::string pc_Type = prm_json.get<std::string>("preconditioner.type");
    if (pc_Type == "cpr") {
	prm_json.put("preconditioner.coarsesolver.preconditioner.verbosity", 10);
    }
    prm_json.put("preconditioner.verbosity", 10);
    prm_json.put("verbosity", 2);
    
    Mat* matrix = &systems[0];
    auto glo = std::make_unique<GLO>(*matrix, *parComm);

    // Create QuasiImpesWeights function used for CPR
    int pidx = 1;
    if (block_size == 2) { pidx = 0; }
    std::function<Vec()> quasi;
    quasi = [matrix, pidx]() {
	    return Opm::Amg::getQuasiImpesWeights<Mat, Vec>(*matrix, pidx, false);
	};

    //Create linear solver. AMG hierarchy is set-up here based on the systems[0] matrix 
    auto fs_json = std::make_unique<FlexibleSolverType>(*glo, *parComm, prm_json, quasi, pidx);
    
    for (int i = 0; i < systems.size(); ++i) {

	Dune::InverseOperatorResult stat;
	
	Vec crhs(rhs[i]);
	Vec x(crhs.size());
	x=0;

	matrix = &systems[i];
	fs_json->preconditioner().update();
	fs_json->apply(x, crhs, prm_json.get<double>("tol", 0.001), stat);

    }
}

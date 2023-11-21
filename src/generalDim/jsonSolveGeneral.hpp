/*
  Copyright 2023 Andreas Thune.

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

#ifndef OPM_JSONSOLVEGENERAL_HEADER_INCLUDED
#define OPM_JSONSOLVEGENERAL_HEADER_INCLUDED

#endif // OPM_JSONSOLVEGENERAL_HEADER_INCLUDED

template<class Mat, class Vec>
void gen_dim_jsonSolve(int argc, char** argv)
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
    typedef Dune::FlexibleSolver<Mat, Vec> FlexibleSolverType;

    const auto block_size = Vec::block_type::dimension;
    
    CollectiveCommunication cc(MPI_COMM_WORLD);
    int rank = cc.rank();

    Mat A_loc, A_loc_;
    Vec rhs, rhs_loc;

    DictRead DR;
    Comm comm(cc);
    std::shared_ptr<Comm> parComm(new(Comm));
    readMatOnRootAndDist(argc, argv, A_loc, rhs_loc, DR, comm, parComm, cc); // in partition/overlapCreation.hpp

    //findZeroDiag(A_loc, rhs_loc);
    ScalarProduct sp(*parComm);
    GLO linOp(A_loc, *parComm);
    
    Opm::FlowLinearSolverParameters flsp_json;
    flsp_json.linsolver_ = DR.dict[12];

    Opm::PropertyTree prm_json(flsp_json.linsolver_);

    std::string pc_Type = prm_json.get<std::string>("preconditioner.type");
    if (pc_Type == "cpr") {
	prm_json.put("preconditioner.coarsesolver.preconditioner.verbosity", 10);
    }
    prm_json.put("preconditioner.verbosity", 10);
    prm_json.put("verbosity", 2);
    //prm_json.put("tol", 0.001);
    //prm_json.put("maxiter", 30);
    std::function<Vec()> quasi;
    int pidx = 1;
    if (block_size == 2)
	pidx = 0;
    
    quasi = [A_loc, pidx]() {
	return Opm::Amg::getQuasiImpesWeights<Mat, Vec>(A_loc, pidx, false);
    };
    
    auto fs_json = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_json, quasi, pidx);

    Vec x(rhs_loc.size());
    x=0;
    Vec crhs(rhs_loc);
    Dune::InverseOperatorResult stat;
    fs_json->apply(x, crhs, prm_json.get<double>("tol", 0.001), stat);
}


template<class Mat, class Vec>
void gen_dim_jsonSolve_extended(int argc, char** argv)
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
    typedef Dune::FlexibleSolver<Mat, Vec> FlexibleSolverType;

    CollectiveCommunication cc(MPI_COMM_WORLD);
    int rank = cc.rank();

    Mat A_loc, A_loc_;
    Vec rhs, rhs_loc;

    DictRead DR;
    Comm comm(cc);
    std::shared_ptr<Comm> parComm(new(Comm));
    readMatOnRootAndDist(argc, argv, A_loc, rhs_loc, DR, comm, parComm, cc); // in partition/overlapCreation.hpp

    findZeroDiag(A_loc, rhs_loc);
    ScalarProduct sp(*parComm);
    GLO linOp(A_loc, *parComm);
    
    Opm::FlowLinearSolverParameters flsp_json;
    flsp_json.linsolver_ = DR.dict[12];

    Opm::PropertyTree prm_json(flsp_json.linsolver_);

    std::string pc_Type = prm_json.get<std::string>("preconditioner.type");
    if (pc_Type == "cpr") {
	prm_json.put("preconditioner.coarsesolver.preconditioner.verbosity", 10);
    }
    prm_json.put("preconditioner.verbosity", 10);
    prm_json.put("verbosity", 2);
    prm_json.put("tol", 0.001);
    prm_json.put("maxiter", 30);
    std::function<Vec()> quasi;
    int pidx = 1;

    quasi = [A_loc, pidx]() {
	return Opm::Amg::getQuasiImpesWeights<Mat, Vec>(A_loc, pidx, false);
    };
    
    auto fs_json = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_json, quasi, pidx);

    Vec x(rhs_loc.size());
    x=0;
    Vec crhs(rhs_loc);
    Dune::InverseOperatorResult stat;
    fs_json->apply(x, crhs, 0.001, stat);

    //flexibleSolverAndPreconTimer(*fs_json, rhs_loc, cc, std::string("CPR_SLEIPNER"), 0.0001, true);

    using CriterionBase
	= Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricDependency<Mat, Dune::Amg::FirstDiagonal>>;
    using Criterion = Dune::Amg::CoarsenCriterion<CriterionBase>;

    using Smoother = ILU;
    using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
    typedef Dune::Amg::AMGCPR<GLO, Vec, ILU, Comm> AMGCPR;
    
    //Criterion criterion(15, 1200,1.2,1.6,static_cast<Dune::Amg::AccumulationMode>(1));//, true);

    Criterion criterion(15, prm_json.get<int>("coarsenTarget", 1200));
    auto pc_child = prm_json.get_child_optional("preconditioner");
    setCrit(criterion, *pc_child );
    criterion.setRandomParallelGhostIndexOrder(false);
    SmootherArgs smootherArgs;
    setOpmILU0args(smootherArgs, *pc_child);
    auto amg = std::make_shared<AMGCPR>(linOp, criterion, smootherArgs, *parComm);

    int verb = 0;
    if (rank == 0) {verb=2;}
    Solver solver(linOp, sp, *amg, 0.01, 30, verb);
    Vec x2(rhs_loc.size());
    Vec crhs2(rhs_loc);
    x2=0;
    Dune::InverseOperatorResult stat2;
    solver.apply(x2, crhs2, stat2);    
}

template<class Mat, class Vec>
void gen_dim_jsonSolve_compare_parallel(int argc, char** argv)
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
    typedef Dune::FlexibleSolver<Mat, Vec> FlexibleSolverType;

    CollectiveCommunication cc(MPI_COMM_WORLD);
    int rank = cc.rank();

    Mat A_loc, A_loc_;
    Vec rhs, rhs_loc;

    DictRead DR;
    Comm comm(cc);
    std::shared_ptr<Comm> parComm(new(Comm));
    readMatOnRootAndDist(argc, argv, A_loc, rhs_loc, DR, comm, parComm, cc); // in partition/overlapCreation.hpp

    ScalarProduct sp(*parComm);
    GLO linOp(A_loc, *parComm);
    
    Opm::FlowLinearSolverParameters flsp_json;
    flsp_json.linsolver_ = DR.dict[12];

    Opm::PropertyTree prm_json(flsp_json.linsolver_);

    std::string pc_Type = prm_json.get<std::string>("preconditioner.type");
    if (pc_Type == "cpr") {
	prm_json.put("preconditioner.coarsesolver.preconditioner.verbosity", 10);
    }
    prm_json.put("preconditioner.verbosity", 10);
    prm_json.put("verbosity", 2);
    prm_json.put("tol", 0.000001);
    prm_json.put("maxiter", 300);
    std::function<Vec()> quasi;
    int pidx = 1;

    quasi = [A_loc, pidx]() {
	return Opm::Amg::getQuasiImpesWeights<Mat, Vec>(A_loc, pidx, false);
    };
    
    auto fs_json = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_json, quasi, pidx);

    Vec x(rhs_loc.size());
    x=0;
    Vec crhs(rhs_loc);
    Dune::InverseOperatorResult stat;
    fs_json->apply(x, crhs, prm_json.get<double>("tol", 0.001), stat);

    bool writeSolution = true;
    if (writeSolution) {

	if (cc.size() == 1) {

	    std::ofstream xfile("Solution.vec");
	    xfile.precision (std::numeric_limits<double>::digits10 + 1);
	    Dune::writeMatrixMarket(x, xfile);
	    xfile.close();
	}
	else {
	    typedef Dune::OwnerOverlapCopyAttributeSet::AttributeSet AttributeSet;
	    Vec serialSol;
	    readMatMarketObject(serialSol ,"Solution.vec");

	    auto indexSet = parComm->indexSet();

	    double sumDiff = 0;
	    double infNorm = 0;
	    int ownerCells = 0;

	    double pDiff = 0;
	    double satW = 0;
	    double satG = 0;
	    for (auto idx = indexSet.begin(); idx!=indexSet.end(); ++idx) {

		if (idx->local().attribute()==AttributeSet::owner) {

		    auto diff = x[idx->local().local()] - serialSol[idx->global()-1];
		    //auto diff = rhs_loc[idx->local().local()] - rhs[idx->global()];
		    
		    sumDiff += diff.two_norm ();
		    //infNorm = std::max(infNorm, diff.infinity_norm ());
		    ownerCells += 1;

		    pDiff += std::abs(diff[1]);
		    satW += std::abs(diff[0]);
		    satG += std::abs(diff[2]);

		    if (idx->global() == 0) {
			std::cout << idx->local().local()<< " " << idx->global()<< " "<< x[idx->local().local()][1]<< " " << serialSol[idx->global()-1][1] << std::endl;
		    }
		    if (false) {
			std::cout << idx->local().local()<< " " << idx->global()<< " "<< x[idx->local().local()][1]<< " " << serialSol[idx->global()-1][1] << std::endl;
		    }
		    
		}
	    }

	    //std::cout<< rank << ": "<< sumDiff<< " " << sumDiff/ownerCells<< " "<< infNorm<< std::endl;

	    std::cout<< rank << ": "<< pDiff<< " " << satW<< " " << satG<< " " << pDiff/ownerCells<< " " << satW/ownerCells<< " " << satG/ownerCells<< std::endl;
	}
    }
}


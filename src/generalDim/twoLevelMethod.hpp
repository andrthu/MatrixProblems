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

#ifndef OPM_TWOLEVELMETHOD_HEADER_INCLUDED
#define OPM_TWOLEVELMETHOD_HEADER_INCLUDED

#endif // OPM_TWOLEVELMETHOD_HEADER_INCLUDED

template<class Mat, class Vec>
void gen_dim_twoLevel(int argc, char** argv)
{
    typedef Dune::MPIHelper::MPICommunicator MPICommunicator;
    typedef Dune::CollectiveCommunication<MPICommunicator> CollectiveCommunication;
    typedef Dune::BiCGSTABSolver<Vec> Solver;
    typedef Dune::InverseOperatorResult Stat;

    typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
    typedef Dune::OverlappingSchwarzScalarProduct<Vec,Comm> ScalarProduct; 
    typedef Dune::OverlappingSchwarzOperator<Mat,Vec,Vec,Comm> Operator;
    typedef GhostLastMatrixAdapter<Mat,Vec,Vec,Comm> GLO;
    typedef Opm::ParallelOverlappingILU0<Mat,Vec,Vec,Comm> ILU;
    typedef Dune::Amg::AMGCPR<GLO, Vec, ILU, Comm> AMGCPR;
    typedef Dune::FlexibleSolver<Mat, Vec> FlexibleSolverType;

    typedef ILU Smoother;
    typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
    typedef Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricDependency<Mat, Dune::Amg::FirstDiagonal>> CriterionBase;
    typedef Dune::Amg::CoarsenCriterion<CriterionBase> Criterion;
   
    typedef Dune::SeqILU<Mat,Vec,Vec> DuneSmoother;
    typedef Dune::BlockPreconditioner<Vec, Vec, Comm, DuneSmoother> ParSmoother;
    typedef typename Dune::Amg::SmootherTraits<ParSmoother>::Arguments Smoother1Args;
    typedef Dune::Amg::AMGCPR<Operator, Vec, ParSmoother, Comm> AMGCPR_DUNE;

    //twoleveltypedef
    typedef Dune::FieldMatrix<double,1,1> CoarseBlockMat;
    typedef Dune::BCRSMatrix<CoarseBlockMat> CoarseMat;
    typedef Dune::FieldVector<double,1> CoarseBlockVec;
    typedef Dune::BlockVector<CoarseBlockVec> CoarseVec;
    typedef GhostLastMatrixAdapter<CoarseMat,CoarseVec,CoarseVec,Comm> CGLO;
    typedef Opm::ParallelOverlappingILU0<CoarseMat,CoarseVec,CoarseVec,Comm> CILU;
    typedef Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricDependency<CoarseMat, Dune::Amg::FirstDiagonal>> CoarseCriterionBase;
    typedef Dune::Amg::CoarsenCriterion<CoarseCriterionBase> CoarseCriterion;
    
    typedef Opm::PressureTransferPolicy<GLO, CGLO, Comm, false> LevelTransferPolicy;
    typedef OneStepAMGCoarseSolverPolicyCpr<CGLO, CILU, CoarseCriterion, Comm, LevelTransferPolicy> AMGSolver;
    typedef Dune::Amg::PressureSolverPolicy<CGLO,AMGSolver,LevelTransferPolicy> CoarseSolverPolicy;
    typedef Dune::Amg::TwoLevelMethodCpr<GLO, AMGSolver, ILU> TwoLevelMethod;

    const auto block_size = Vec::block_type::dimension;
    // ------------------------------------------------ 
    // --- start reading input 
    CollectiveCommunication cc(MPI_COMM_WORLD);
    int rank = cc.rank();

    Mat A_loc;
    Vec rhs_loc;

    DictRead DR;
    Comm comm(cc);
    std::shared_ptr<Comm> parComm(new(Comm));
    readMatOnRootAndDist(argc, argv, A_loc, rhs_loc, DR, comm, parComm, cc);
    // --- Complete reading input

    // ------------------------------------------------
    // --- Set up operators, preconditioners and solvers 
    ScalarProduct sp(*parComm);
    Operator linOp(A_loc, *parComm);
    GLO glLinOp(A_loc, *parComm);

    Opm::FlowLinearSolverParameters flsp_json;
    flsp_json.linsolver_ = std::string("/home/andreast/fork_opm/test/linear_solver_config_files/2021/cpr/sleipner_cpr_verbose_quasiimpes_tol3.json");
    Opm::PropertyTree prm_json(flsp_json.linsolver_);
    std::function<Vec()> quasi;
    int pidx = 1;
    if (block_size == 2)
	pidx = 0;
    quasi = [A_loc, pidx]() {
	return Opm::Amg::getQuasiImpesWeights<Mat, Vec>(A_loc, pidx, false);
    };

    //auto cpr = std::make_shared<Dune::OwningTwoLevelPreconditioner<GLO, Vec, false, Comm>>(glLinOp, prm_json, quasi, pidx, *parComm);
    auto fs_json = std::make_unique<FlexibleSolverType>(glLinOp, *parComm, prm_json, quasi, pidx);
    
    //flexibleSolverAndPreconTimer(*fs_json, rhs_loc, cc, std::string("CPR_SLEIPNER"), 0.01, true);
    
    std::string use_ilu("ILU");
    auto ilu_help = Opm::convertString2Milu(use_ilu);
    std::shared_ptr<ILU> ilu(new ILU(A_loc, *parComm, 1, ilu_help, A_loc.N(), false, false) );

    //Set up coarse solver and transfer policy
    CoarseCriterion critCPR(15, prm_json.get<int>("preconditioner.coarsesolver.preconditioner.coarsenTarget", 1200));
    setCritCPR(critCPR, prm_json);
    SmootherArgs smootherArgsCPR;
    setOpmILU0Defargs(smootherArgsCPR);
    std::shared_ptr<AMGSolver> csp(new AMGSolver(smootherArgsCPR, critCPR));
    auto qWgt = quasi();
    LevelTransferPolicy ltp(*parComm, qWgt, pidx);

    if (rank == 0) {std::cout << "Start tlm" << std::endl;}
    auto tlm = std::make_shared<TwoLevelMethod>(glLinOp, ilu, ltp, *csp, 1, 1);
    if (rank == 0) {std::cout << "Made tlm" << std::endl;}

    Vec x(rhs_loc.size());
    x=0;
    tlm->pre(x,rhs_loc);
    cc.barrier();
    Dune::Timer timer;
    timer.start();
    tlm->apply(x,rhs_loc);
    cc.barrier();
    double cprT = timer.stop();
    if (rank == 0)
	std::cout << "TLM-timer "<< cprT <<std::endl;
    auto amgS = tlm->getCoarseSolver();
    auto amgC = amgS->getCoarseAMG();
    
    timeFineLevel(glLinOp, rhs_loc, cc, parComm, ltp);
    
    CoarseVec rhs_coarse(rhs_loc.size());
    rhs_coarse=0;

    for (auto b = rhs_loc.begin(); b!=rhs_loc.end(); ++b) {

	auto bw = qWgt[b.index()];
	for (size_t i = 0; i < b->size(); ++i) {
	    rhs_coarse[b.index()] += (*b)[i] * bw[i];
	}
    }
       
    amgHiaInfo(amgC, rhs_coarse, cc);
    
}

template<class Mat, class Vec>
void gen_two_dim_old_co2(int argc, char** argv)
{
    typedef Dune::MPIHelper::MPICommunicator MPICommunicator;
    typedef Dune::CollectiveCommunication<MPICommunicator> CollectiveCommunication;
    typedef Dune::BiCGSTABSolver<Vec> Solver;
    typedef Dune::InverseOperatorResult Stat;

    typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
    typedef Dune::OverlappingSchwarzScalarProduct<Vec,Comm> ScalarProduct; 
    typedef Dune::OverlappingSchwarzOperator<Mat,Vec,Vec,Comm> Operator;
    typedef GhostLastMatrixAdapter<Mat,Vec,Vec,Comm> GLO;
    typedef Opm::ParallelOverlappingILU0<Mat,Vec,Vec,Comm> ILU;
    typedef Dune::Amg::AMGCPR<GLO, Vec, ILU, Comm> AMGCPR;
    typedef Dune::FlexibleSolver<Mat, Vec> FlexibleSolverType;

    typedef ILU Smoother;
    typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
    typedef Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricDependency<Mat, Dune::Amg::FirstDiagonal>> CriterionBase;
    typedef Dune::Amg::CoarsenCriterion<CriterionBase> Criterion;
   
    typedef Dune::SeqILU<Mat,Vec,Vec> DuneSmoother;
    typedef Dune::BlockPreconditioner<Vec, Vec, Comm, DuneSmoother> ParSmoother;
    typedef typename Dune::Amg::SmootherTraits<ParSmoother>::Arguments Smoother1Args;
    typedef Dune::Amg::AMGCPR<Operator, Vec, ParSmoother, Comm> AMGCPR_DUNE;

    //twoleveltypedef
    typedef Dune::FieldMatrix<double,1,1> CoarseBlockMat;
    typedef Dune::BCRSMatrix<CoarseBlockMat> CoarseMat;
    typedef Dune::FieldVector<double,1> CoarseBlockVec;
    typedef Dune::BlockVector<CoarseBlockVec> CoarseVec;
    typedef GhostLastMatrixAdapter<CoarseMat,CoarseVec,CoarseVec,Comm> CGLO;
    typedef Opm::ParallelOverlappingILU0<CoarseMat,CoarseVec,CoarseVec,Comm> CILU;
    typedef Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricDependency<CoarseMat, Dune::Amg::FirstDiagonal>> CoarseCriterionBase;
    typedef Dune::Amg::CoarsenCriterion<CoarseCriterionBase> CoarseCriterion;
    
    typedef Opm::PressureTransferPolicy<GLO, CGLO, Comm, false> LevelTransferPolicy;
    typedef OneStepAMGCoarseSolverPolicyCpr<CGLO, CILU, CoarseCriterion, Comm, LevelTransferPolicy> AMGSolver;
    typedef Dune::Amg::PressureSolverPolicy<CGLO,AMGSolver,LevelTransferPolicy> CoarseSolverPolicy;
    typedef Dune::Amg::TwoLevelMethodCpr<GLO, AMGSolver, ILU> TwoLevelMethod;
    
    // ------------------------------------------------ 
    // --- start reading input 
    CollectiveCommunication cc(MPI_COMM_WORLD);
    int rank = cc.rank();

    Mat A_loc;
    Vec rhs_loc;

    DictRead DR;
    Comm comm(cc);
    std::shared_ptr<Comm> parComm(new(Comm));
    readMatOnRootAndDist(argc, argv, A_loc, rhs_loc, DR, comm, parComm, cc);
    // --- Complete reading input

    // ------------------------------------------------
    // --- Set up operators, preconditioners and solvers 
    ScalarProduct sp(*parComm);
    Operator linOp(A_loc, *parComm);
    GLO glLinOp(A_loc, *parComm);

    Opm::FlowLinearSolverParameters flsp_json;
    flsp_json.linsolver_ = std::string("/home/andreast/fork_opm/test/linear_solver_config_files/2021/cpr/sleipner_cpr_verbose_quasiimpes_tol3.json");
    Opm::PropertyTree prm_json(flsp_json.linsolver_);
    std::function<Vec()> quasi;
    int pidx = 0;
    quasi = [A_loc, pidx]() {
	return Opm::Amg::getQuasiImpesWeights<Mat, Vec>(A_loc, pidx, false);
    };

    auto cpr = std::make_shared<Dune::OwningTwoLevelPreconditioner<GLO, Vec, false, Comm>>(glLinOp, prm_json, quasi, pidx, *parComm);

    std::string use_ilu("ILU");
    auto ilu_help = Opm::convertString2Milu(use_ilu);
    std::shared_ptr<ILU> ilu(new ILU(A_loc, *parComm, 1, ilu_help, A_loc.N(), false, false) );
    
    CoarseCriterion critCPR(15, prm_json.get<int>("preconditioner.coarsesolver.preconditioner.coarsenTarget", 1200));
    setCritCPR(critCPR, prm_json);
    SmootherArgs smootherArgsCPR;
    setOpmILU0Defargs(smootherArgsCPR);
    std::shared_ptr<AMGSolver> csp(new AMGSolver(smootherArgsCPR, critCPR));
    auto qWgt = quasi();
    LevelTransferPolicy ltp(*parComm, qWgt, pidx);

    if (rank == 0) {std::cout << "Start tlm" << std::endl;}
    auto tlm = std::make_shared<TwoLevelMethod>(glLinOp, ilu, ltp, *csp, 1, 1);
    if (rank == 0) {std::cout << "Made tlm" << std::endl;}

    auto amgS = tlm->getCoarseSolver();
    auto amgC = amgS->getCoarseAMG();

    timeFineLevel(glLinOp, rhs_loc, cc, parComm);
    
    CoarseVec rhs_coarse(rhs_loc.size());
    rhs_coarse=0;

    for (auto b = rhs_loc.begin(); b!=rhs_loc.end(); ++b) {

	auto bw = qWgt[b.index()];
	for (size_t i = 0; i < b->size(); ++i) {
	    rhs_coarse[b.index()] += (*b)[i] * bw[i];
	}
    }
    
    amgHiaInfo(amgC, rhs_coarse, cc);
}

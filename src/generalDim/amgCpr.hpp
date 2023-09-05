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

#ifndef OPM_AMGCPR_HEADER_INCLUDED
#define OPM_AMGCPR_HEADER_INCLUDED

#endif // OPM_AMGCPR_HEADER_INCLUDED

template<class Mat, class Vec>
void gen_dim_amgCprTest(int argc, char** argv)
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

    using Smoother = ILU;
    using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
    using CriterionBase
	= Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricDependency<Mat, Dune::Amg::FirstDiagonal>>;
    using Criterion = Dune::Amg::CoarsenCriterion<CriterionBase>;


    
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
    
    Opm::FlowLinearSolverParameters flsp_amg;
    Opm::PropertyTree prm_amg = setupAMG(std::string("amg"), flsp_amg);
    prm_amg.put("preconditioner.verbosity", 10);
    
    Criterion criterion(15, prm_amg.get<int>("coarsenTarget", 1200));
    setCrit(criterion, prm_amg);

    SmootherArgs smootherArgs;
    setOpmILU0args(smootherArgs, prm_amg);

    Smoother1Args smoother1Args;
    smoother1Args.iterations = prm_amg.get<int>("iterations", 1);
    
    //typedef std::numeric_limits< double > dbl;
    //std::cout.precision(dbl::max_digits10);
    
    auto amg = std::make_shared<AMGCPR>(glLinOp, criterion, smootherArgs, *parComm);
    if (rank == 0) {std::cout << std::endl;}
    if (rank == 0) {std::cout << std::endl;}
    //auto amg2 = std::make_shared<AMGCPR>(linOp, criterion, smootherArgs, *parComm);
    // --- Complete set up

    int verb=0;
    if (rank==0) {verb=2;}
    Solver solver(linOp, sp, *amg, 0.01, 30, verb);
    Solver solverCpr(glLinOp, sp, *cpr, 0.01, 30, verb);
    
    Vec x(rhs_loc.size());
    Vec x_cpr(rhs_loc.size());
    Vec crhs1(rhs_loc);
    Vec crhs2(rhs_loc);
    x=0;x_cpr=0;
    Dune::InverseOperatorResult stat,stat_cp3;
    if (rank==0) {std::cout << std::endl;}
    if (rank==0) {std::cout << "CPR_solve" <<std::endl;}
    solverCpr.apply(x_cpr,crhs2,stat_cp3);
    solver.apply(x, crhs1,stat);
    
    auto amg_dune = std::make_shared<AMGCPR_DUNE>(linOp, criterion, smoother1Args, *parComm);

    Solver solver2(linOp, sp, *amg_dune, 0.01, 30, verb);
    solverAndPreconTimer(solver, *amg, rhs_loc, cc, std::string("AMG"), 0.01);
    solverAndPreconTimer(solver2, *amg_dune, rhs_loc, cc, std::string("AMG_DUNE"), 0.01);
    /*
    Vec x2(rhs_loc.size());
    x2=0;
    Dune::InverseOperatorResult stat2;
    solver2.apply(x2, crhs2,stat);
    */
    std::vector<std::size_t> ag, ag2;
    amg->getCoarsestAggregateNumbers(ag);
    amgHiaInfo(*amg, rhs_loc, cc); // in partition/flexibleSolverTimer.hpp
}

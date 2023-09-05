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
void gen_dim_gl(int argc, char** argv)
{
    typedef Dune::MPIHelper::MPICommunicator MPICommunicator;
    typedef Dune::CollectiveCommunication<MPICommunicator> CollectiveCommunication;
    typedef Dune::BiCGSTABSolver<Vec> Solver;
    typedef Dune::InverseOperatorResult Stat;

    typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
    typedef Dune::OverlappingSchwarzScalarProduct<Vec,Comm> ScalarProduct; 
    typedef Dune::OverlappingSchwarzOperator<Mat,Vec,Vec,Comm> Operator;
    typedef GhostLastMatrixAdapter<Mat,Vec,Vec,Comm> GLO;
    typedef OverlappingSchwarzOperatorCopy<Mat,Vec,Vec,Comm> GLO2;
    typedef Opm::ParallelOverlappingILU0<Mat,Vec,Vec,Comm> ILU;
    typedef GhostLastILU0<Mat,Vec,Vec,Comm> GLILU;
    typedef Dune::Amg::AMGCPR<GLO, Vec, ILU, Comm> AMGCPR;
    typedef Dune::FlexibleSolver<Mat, Vec> FlexibleSolverType;

    typedef ILU Smoother;
    typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
    typedef Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricDependency<Mat, Dune::Amg::FirstDiagonal>> CriterionBase;
    typedef Dune::Amg::CoarsenCriterion<CriterionBase> Criterion;
   
    typedef ParallelOverlappingILU0<Mat,Vec,Vec,Comm> CopyILU;//
    typedef Dune::SeqILU<Mat,Vec,Vec> DuneSmoother;
    typedef Dune::BlockPreconditioner<Vec, Vec, Comm, DuneSmoother> ParSmoother;
    typedef typename Dune::Amg::SmootherTraits<CopyILU>::Arguments Smoother1Args;
    typedef Dune::Amg::AMGCPR<Operator, Vec, CopyILU, Comm> AMGCPR_CILU;
    
    typedef typename Dune::Amg::SmootherTraits<GLILU>::Arguments Smoother2Args;
    typedef Dune::Amg::AMGCPR<GLO, Vec, GLILU, Comm> AMGCPR_GLILU;

    //twoleveltypedef
    typedef Dune::FieldMatrix<double,1,1> CoarseBlockMat;
    typedef Dune::BCRSMatrix<CoarseBlockMat> CoarseMat;
    typedef Dune::FieldVector<double,1> CoarseBlockVec;
    typedef Dune::BlockVector<CoarseBlockVec> CoarseVec;
    typedef GhostLastMatrixAdapter<CoarseMat,CoarseVec,CoarseVec,Comm> CGLO;
    typedef Opm::ParallelOverlappingILU0<CoarseMat,CoarseVec,CoarseVec,Comm> CILU;
    typedef GhostLastILU0<CoarseMat,CoarseVec,CoarseVec,Comm> CGLILU;
    typedef Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricDependency<CoarseMat, Dune::Amg::FirstDiagonal>> CoarseCriterionBase;
    typedef Dune::Amg::CoarsenCriterion<CoarseCriterionBase> CoarseCriterion;
    
    typedef Opm::PressureTransferPolicy<GLO, CGLO, Comm, false> LevelTransferPolicy;
    typedef OneStepAMGCoarseSolverPolicyCpr<CGLO, CGLILU, CoarseCriterion, Comm, LevelTransferPolicy> AMGSolver;
    typedef Dune::Amg::PressureSolverPolicy<CGLO, AMGSolver,LevelTransferPolicy> CoarseSolverPolicy;
    typedef Dune::Amg::TwoLevelMethodCpr<GLO, AMGSolver, GLILU> TwoLevelMethod;

    //twoleveltypedef no GL

    typedef Dune::OverlappingSchwarzOperator<CoarseMat,CoarseVec,CoarseVec,Comm> COperator;
    typedef ParallelOverlappingILU0<CoarseMat,CoarseVec,CoarseVec,Comm> CCopyILU;
    typedef Opm::PressureTransferPolicy<Operator, COperator, Comm, false> NoGLLevelTransferPolicy;
    typedef OneStepAMGCoarseSolverPolicyCpr<COperator, CCopyILU, CoarseCriterion, Comm, NoGLLevelTransferPolicy> NoGLAMGSolver;
    typedef Dune::Amg::PressureSolverPolicy<COperator,NoGLAMGSolver,NoGLLevelTransferPolicy> NoGLCoarseSolverPolicy;
    typedef Dune::Amg::TwoLevelMethodCpr<Operator, NoGLAMGSolver, CopyILU> NoGLTwoLevelMethod;


    //twoleveltypedef no dune_GL

    typedef Dune::SeqILU<CoarseMat,CoarseVec,CoarseVec> CDuneSmoother;
    typedef Dune::BlockPreconditioner<CoarseVec, CoarseVec, Comm, CDuneSmoother> CParSmoother;
    typedef Opm::PressureTransferPolicy<GLO, CGLO, Comm, false> DuneLevelTransferPolicy;
    typedef OneStepAMGCoarseSolverPolicyCpr<CGLO, CParSmoother, CoarseCriterion, Comm, DuneLevelTransferPolicy> DuneAMGSolver;
    typedef Dune::Amg::PressureSolverPolicy<CGLO,DuneAMGSolver, DuneLevelTransferPolicy> DuneoarseSolverPolicy;
    typedef Dune::Amg::TwoLevelMethodCpr<GLO, DuneAMGSolver, ParSmoother> DuneTwoLevelMethod;
    
    const auto block_size = Vec::block_type::dimension;
    // ------------------------------------------------ 
    // --- start reading input 
    CollectiveCommunication cc(MPI_COMM_WORLD);
    int rank = cc.rank();

    Mat A_loc, A_loc_with_ghost, A_2;
    Vec rhs_loc, rhs_loc2;

    DictRead DR;
    Comm comm(cc);
    std::shared_ptr<Comm> parComm(new(Comm));
    std::shared_ptr<Comm> parComm2(new(Comm));
    
    readMatOnRootAndDist(argc, argv, A_loc, rhs_loc, DR, comm, parComm, cc, false);
    //readMatOnRootAndDist(argc, argv, A_loc, rhs_loc2, DR, comm, parComm2, cc, false);

    buildLocalMatrixFromLocAndComm(A_loc, A_2, *parComm);
    
    // --- Complete reading input

    // ------------------------------------------------
    // --- Set up operators, preconditioners and solvers 
    ScalarProduct sp(*parComm);
    Operator linOp(A_loc, *parComm);
    GLO glLinOp(A_2, *parComm);
    GLO2 glLinOp2(A_loc, *parComm);
    
    multipleMinLoopTimeSpMVAS(cc,linOp, rhs_loc, 3);
    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimeSpMVAS(cc,glLinOp, rhs_loc, 3, true);
    //if (rank == 0) {std::cout << std::endl;}
    //multipleMinLoopTimeSpMVAS(cc,glLinOp2, rhs_loc, 3);
    
    Opm::FlowLinearSolverParameters flsp_json;
    flsp_json.linsolver_ = std::string("/home/andreast/fork_opm/test/linear_solver_config_files/2021/cpr/sleipner_cpr_verbose_quasiimpes_tol3.json");
    Opm::PropertyTree prm_json(flsp_json.linsolver_);
    std::function<Vec()> quasi, quasi2;
    int pidx = 1;
    if (block_size == 2)
	pidx = 0;
    quasi = [A_loc, pidx]() {
	return Opm::Amg::getQuasiImpesWeights<Mat, Vec>(A_loc, pidx, false);
    };
    quasi = [A_2, pidx]() {
	return Opm::Amg::getQuasiImpesWeights<Mat, Vec>(A_2, pidx, false);
    };

    //auto cpr = std::make_shared<Dune::OwningTwoLevelPreconditioner<GLO, Vec, false, Comm>>(glLinOp, prm_json, quasi, pidx, *parComm);
    //auto fs_json = std::make_unique<FlexibleSolverType>(glLinOp, *parComm, prm_json, quasi, pidx);
    //auto fs_json2 = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_json, quasi, pidx);
    
    //flexibleSolverAndPreconTimer(*fs_json, rhs_loc, cc, std::string("CPR_SLEIPNER"), 0.01, true);
    //flexibleSolverAndPreconTimer(*fs_json2, rhs_loc, cc, std::string("CPR_SLEIPNER"), 0.01, true);
    
    std::string use_ilu("ILU");
    auto ilu_help = Opm::convertString2Milu(use_ilu);
    std::shared_ptr<ILU> ilu(new ILU(A_2, *parComm, 1, ilu_help, A_loc.N(), false, false) ); //opm ilu + gl
    std::shared_ptr<CopyILU> cilu(new CopyILU(A_loc, *parComm, 1, A_loc.N()) );                //copy of opm::ilu - gl
    std::shared_ptr<DuneSmoother> seq_ilu(new DuneSmoother(A_loc,0,1));                        //dune seq-ilu
    std::shared_ptr<ParSmoother> dune_ilu(new ParSmoother(seq_ilu, *parComm));                 //dune block-jac+seq-ilu
    std::shared_ptr<GLILU> glilu(new GLILU(A_2, *parComm, 1, A_loc.N()) );                   //copy of opm::ilu

    Opm::FlowLinearSolverParameters flsp_amg;
    Opm::PropertyTree prm_amg = setupAMG(std::string("amg"), flsp_amg);
    prm_amg.put("preconditioner.verbosity", 10);
    Criterion AMGcriterion(15, prm_amg.get<int>("coarsenTarget", 1200));
    setCrit(AMGcriterion, prm_amg);
    Smoother1Args AMGsmootherArgs1;
    Smoother2Args AMGsmootherArgs2;
    setOpmILU0noMILU(AMGsmootherArgs1, prm_amg);
    setOpmILU0noMILU(AMGsmootherArgs2, prm_amg);

    std::shared_ptr<AMGCPR_GLILU> amgGL (new AMGCPR_GLILU (glLinOp, AMGcriterion, AMGsmootherArgs2, *parComm) );
    std::shared_ptr<AMGCPR_CILU> amgNoGL (new AMGCPR_CILU(linOp, AMGcriterion, AMGsmootherArgs1, *parComm) );
    
    if (rank == 0) {std::cout << std::endl;}
    Dune::Timer updateTimer;
    cc.barrier();
    updateTimer.start();
    ilu->update();
    cc.barrier();
    double utime = updateTimer.stop();
    if (rank == 0) {std::cout << "Update OPM-ILU time "<< utime << std::endl;}

    cc.barrier();
    updateTimer.reset();
    updateTimer.start();
    glilu->updateTest();
    cc.barrier();
    double utime2 = updateTimer.stop();
    if (rank == 0) {std::cout << "UpdateTest GL-ILU time "<< utime2 << std::endl;}

    cc.barrier();
    updateTimer.reset();
    updateTimer.start();
    glilu->update();
    cc.barrier();
    double utime22 = updateTimer.stop();
    if (rank == 0) {std::cout << "Update GL-ILU time "<< utime22 << std::endl;}

    cc.barrier();
    updateTimer.reset();
    updateTimer.start();
    cilu->update();
    cc.barrier();
    double utime3 = updateTimer.stop();
    if (rank == 0) {std::cout << "Update C-ILU time "<< utime3 << std::endl;}
    
    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimePre(cc, *ilu, rhs_loc, std::string("ILU"), 3);
    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimePre(cc, *cilu, rhs_loc, std::string("CILU"), 3);
    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimePre(cc, *dune_ilu, rhs_loc, std::string("DILU"), 3);
    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimePre(cc, *glilu, rhs_loc, std::string("GLILU"), 3);

    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimePre(cc, *amgGL, rhs_loc, std::string("GLAMG"), 3);
    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimePre(cc, *amgNoGL, rhs_loc, std::string("CAMG"), 3);
    
    //Set up coarse solver and transfer policy
    CoarseCriterion critCPR(15, prm_json.get<int>("preconditioner.coarsesolver.preconditioner.coarsenTarget", 1200));
    setCritCPR(critCPR, prm_json);
    Smoother2Args smootherArgsCPR;
    //setOpmILU0Defargs(smootherArgsCPR);
    smootherArgsCPR.iterations = 1;
    smootherArgsCPR.relaxationFactor = 1.0;
    const int iluwitdh_ = 0;
    smootherArgsCPR.setN(iluwitdh_);
    std::shared_ptr<AMGSolver> csp(new AMGSolver(smootherArgsCPR, critCPR));
    auto qWgt = quasi();
    LevelTransferPolicy ltp(*parComm, qWgt, pidx);

    if (rank == 0) {std::cout << "Start tlm" << std::endl;}
    auto tlm = std::make_shared<TwoLevelMethod>(glLinOp, glilu, ltp, *csp, 1, 1);
    if (rank == 0) {std::cout << "Made tlm" << std::endl;}

    //Set up coarse solver and transfer policy NO GLO
    CoarseCriterion critCPR_NOGL(15, prm_json.get<int>("preconditioner.coarsesolver.preconditioner.coarsenTarget", 1200));
    setCritCPR(critCPR_NOGL, prm_json);
    Smoother1Args smootherArgsCPR_NOGL;
    smootherArgsCPR_NOGL.iterations = 1;
    const int iluwitdh = 0;
    smootherArgsCPR_NOGL.setN(iluwitdh);
    smootherArgsCPR_NOGL.relaxationFactor = 1.0;
    std::shared_ptr<NoGLAMGSolver> csp_NOGL(new NoGLAMGSolver(smootherArgsCPR_NOGL, critCPR_NOGL));
    auto qWgt_NOGL = quasi();
    NoGLLevelTransferPolicy ltp_NOGL(*parComm, qWgt_NOGL, pidx);


    if (rank == 0) {std::cout << "Start tlm" << std::endl;}
    auto tlm_NOGL = std::make_shared<NoGLTwoLevelMethod>(linOp, cilu, ltp_NOGL, *csp_NOGL, 1, 1);
    if (rank == 0) {std::cout << "Made tlm" << std::endl;}


    //Set up coarse solver and transfer policy Dune smoother
    CoarseCriterion critCPR_dune(15, prm_json.get<int>("preconditioner.coarsesolver.preconditioner.coarsenTarget", 1200));
    setCritCPR(critCPR_dune, prm_json);
    Smoother1Args smootherArgsCPR_dune;
    smootherArgsCPR_dune.iterations = 1;
    smootherArgsCPR_dune.setN(iluwitdh);
    smootherArgsCPR_dune.relaxationFactor = 1.0;
    std::shared_ptr<DuneAMGSolver> csp_dune(new DuneAMGSolver(smootherArgsCPR_dune, critCPR_dune));
    auto qWgt_dune = quasi();
    DuneLevelTransferPolicy ltp_dune(*parComm, qWgt_dune, pidx);

    if (rank == 0) {std::cout << "Start tlm" << std::endl;}
    auto tlm_dune = std::make_shared<DuneTwoLevelMethod>(glLinOp, dune_ilu, ltp_dune, *csp_dune, 1, 1);
    if (rank == 0) {std::cout << "Made tlm" << std::endl;}

    
    Dune::Timer timer;
    
    Vec xx(rhs_loc.size());
    xx=0;
    tlm_NOGL->pre(xx,rhs_loc);
    cc.barrier();
    timer.reset();
    timer.start();
    tlm_NOGL->apply(xx,rhs_loc);
    cc.barrier();
    double cprT_gl = timer.stop();
    if (rank == 0)
	std::cout << "TLM_NoGL-timer "<< cprT_gl <<std::endl;
    
    Vec x(rhs_loc.size());
    x=0;
    tlm->pre(x,rhs_loc);
    cc.barrier();
    timer.reset();
    timer.start();
    tlm->apply(x,rhs_loc);
    cc.barrier();
    double cprT = timer.stop();
    if (rank == 0)
	std::cout << "TLM-timer "<< cprT <<std::endl;

    Vec xxx(rhs_loc.size());
    xxx=0;
    tlm_dune->pre(xxx,rhs_loc);
    cc.barrier();
    timer.reset();
    timer.start();
    tlm_dune->apply(xxx,rhs_loc);
    cc.barrier();
    double cprT_dune = timer.stop();
    if (rank == 0)
	std::cout << "TLM_dune-timer "<< cprT_dune <<std::endl;

    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimePre(cc, *tlm, rhs_loc, std::string("TLM"), 3);
    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimePre(cc, *tlm_NOGL, rhs_loc, std::string("NoGL-TLM"), 3);
    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimePre(cc, *tlm_dune, rhs_loc, std::string("DTLM"), 3);
    
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

    auto amgS2 = tlm_NOGL->getCoarseSolver();
    auto amgC2 = amgS->getCoarseAMG();

    amgHiaInfo(amgC2, rhs_coarse, cc);
}

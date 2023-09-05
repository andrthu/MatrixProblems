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

#ifndef OPM_UPDATETEST_HEADER_INCLUDED
#define OPM_UPDATETEST_HEADER_INCLUDED

#endif // OPM_UPDATETEST_HEADER_INCLUDED
template<class Mat, class Vec>
void update_test(int argc, char** argv)
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
    typedef Dune::Amg::AMGCPR<Operator, Vec, GLILU, Comm> AMGCPR_GLILU;

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
    GLO2 glLinOp2(A_loc, *parComm);
    
    multipleMinLoopTimeSpMVAS(cc,linOp, rhs_loc, 3);
    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimeSpMVAS(cc,glLinOp, rhs_loc, 3);
    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimeSpMVAS(cc,glLinOp2, rhs_loc, 3);
    
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
    //auto fs_json = std::make_unique<FlexibleSolverType>(glLinOp, *parComm, prm_json, quasi, pidx);
    //auto fs_json2 = std::make_unique<FlexibleSolverType>(linOp, *parComm, prm_json, quasi, pidx);
    
    //flexibleSolverAndPreconTimer(*fs_json, rhs_loc, cc, std::string("CPR_SLEIPNER"), 0.01, true);
    //flexibleSolverAndPreconTimer(*fs_json2, rhs_loc, cc, std::string("CPR_SLEIPNER"), 0.01, true);
    
    std::string use_ilu("ILU");
    auto ilu_help = Opm::convertString2Milu(use_ilu);
    std::shared_ptr<ILU> ilu(new ILU(A_loc, *parComm, 1, ilu_help, A_loc.N(), false, false) ); //opm ilu + gl
    std::shared_ptr<CopyILU> cilu(new CopyILU(A_loc, *parComm, 1, A_loc.N()) );                //copy of opm::ilu - gl
    std::shared_ptr<DuneSmoother> seq_ilu(new DuneSmoother(A_loc,0,1));                        //dune seq-ilu
    std::shared_ptr<ParSmoother> dune_ilu(new ParSmoother(seq_ilu, *parComm));                 //dune block-jac+seq-ilu
    std::shared_ptr<GLILU> glilu(new GLILU(A_loc, *parComm, 1, A_loc.N()) );                   //copy of opm::ilu

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

    glilu->printLUOnRoot(rank);
}

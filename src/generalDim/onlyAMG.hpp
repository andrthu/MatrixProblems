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

#ifndef OPM_ONLYAMG_HEADER_INCLUDED
#define OPM_ONLYAMG_HEADER_INCLUDED

#endif // OPM_ONLYAMG_HEADER_INCLUDED


template<class C, class O, class Comm, class Mat, class Vec>
void sleipnerAMG(const C& cc, O glLinOp, const Comm& parComm, Mat A_loc, Vec rhs_loc)
{
    typedef Dune::FlexibleSolver<Mat, Vec> FlexibleSolverType;

    
    int rank = cc.rank();
    Opm::FlowLinearSolverParameters flsp_sleipner;
    flsp_sleipner.linsolver_ = std::string("/home/andreast/fork_opm/test/linear_solver_config_files/2021/amg/sleipner.json");
    Opm::PropertyTree prm_sleipner (flsp_sleipner.linsolver_);
    prm_sleipner.put("maxiter", 60);
    
    std::function<Vec()> weightsCalculator;
    Dune::Timer timer;
    if (rank == 0) {std::cout << std::endl;}
    if (rank == 0) {std::cout << std::endl;}
    cc.barrier();
    timer.reset();
    timer.start();
    auto fs_sleipner = std::make_unique<FlexibleSolverType>(glLinOp, parComm, prm_sleipner, weightsCalculator, 0);
    double sleipnerTime = timer.stop();
    if (rank == 0) {std::cout << "Setup AMG_SLEIPNER time: "<< sleipnerTime << std::endl;}

    flexibleSolverAndPreconTimer(*fs_sleipner, rhs_loc, cc, std::string("AMG_SLEIPNER"), 0.01, true, false);
}

template<class C, class O, class Comm, class Mat, class Vec>
void sleipnerCPR(const C& cc, O glLinOp, const Comm& parComm, Mat A_loc, Vec rhs_loc)
{
    typedef Dune::FlexibleSolver<Mat, Vec> FlexibleSolverType;

    const auto block_size = Vec::block_type::dimension;
    int rank = cc.rank();
    Opm::FlowLinearSolverParameters flsp_sleipner;
    flsp_sleipner.linsolver_ = std::string("/home/andreast/fork_opm/test/linear_solver_config_files/2021/cpr/sleipner_cpr_verbose_quasiimpes_tol3.json");

    Opm::PropertyTree prm_sleipner (flsp_sleipner.linsolver_);
    prm_sleipner.put("maxiter", 60);

    std::function<Vec()> quasi;

    int pidx = 1;
    if (block_size == 2)
	pidx = 0;
    quasi = [A_loc, pidx]() {
	return Opm::Amg::getQuasiImpesWeights<Mat, Vec>(A_loc, pidx, false);
    };

    auto fs_sleipner = std::make_unique<FlexibleSolverType>(glLinOp, parComm, prm_sleipner, quasi, pidx);

    flexibleSolverAndPreconTimer(*fs_sleipner, rhs_loc, cc, std::string("CPR_SLEIPNER"), 0.01, true, false);
}

template<class C, class O, class Comm, class Mat, class Vec>
void flexibleILU(const C& cc, O glLinOp, const Comm& parComm, Mat A_loc, Vec rhs_loc)
{

    typedef Dune::FlexibleSolver<Mat, Vec> FlexibleSolverType;

    int rank = cc.rank();
    Opm::FlowLinearSolverParameters flsp;
    Opm::PropertyTree prm_ilu = setupILU(std::string("ILU"), flsp);
    prm_ilu.put("verbosity", 2);
    prm_ilu.put("maxiter", 200);
    std::function<Vec()> weightsCalculator;
    
    auto fs = std::make_unique<FlexibleSolverType>(glLinOp, parComm, prm_ilu, weightsCalculator, 0);

    flexibleSolverAndPreconTimer(*fs, rhs_loc, cc, std::string("ILU"), 0.01, true, false);
}

template<class Mat, class Vec>
void gen_dim_only_AMG(int argc, char** argv)
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

    Opm::FlowLinearSolverParameters flsp_default;
    

    Opm::PropertyTree prm_def = setupAMG(std::string("amg"), flsp_default);

    /*
    Opm::FlowLinearSolverParameters flsp_sleipner;
    flsp_sleipner.linsolver_ = std::string("/home/andreast/fork_opm/test/linear_solver_config_files/2021/amg/sleipner.json");
    Opm::PropertyTree prm_sleipner (flsp_sleipner.linsolver_);
    */
    std::function<Vec()> weightsCalculator;

    Dune::Timer timer;
    cc.barrier();
    timer.start();
    auto fs_default = std::make_unique<FlexibleSolverType>(glLinOp, *parComm, prm_def, weightsCalculator, 0);
    double amgTime = timer.stop();
    if (rank == 0) {std::cout << "Setup AMG time: "<< amgTime<< std::endl;}

    /*
    if (rank == 0) {std::cout << std::endl;}
    if (rank == 0) {std::cout << std::endl;}
    cc.barrier();
    timer.reset();
    timer.start();
    auto fs_sleipner = std::make_unique<FlexibleSolverType>(glLinOp, *parComm, prm_sleipner, weightsCalculator, 0);
    double sleipnerTime = timer.stop();
    if (rank == 0) {std::cout << "Setup AMG_SLEIPNER time: "<< amgTime<< std::endl;}
    */
   
    flexibleSolverAndPreconTimer(*fs_default, rhs_loc, cc, std::string("AMG_DEFAULT"), 0.01, true, false);
    //flexibleSolverAndPreconTimer(*fs_sleipner, rhs_loc, cc, std::string("AMG_SLEIPNER"), 0.01, true, false);

    flexibleILU(cc, glLinOp, *parComm, A_loc, rhs_loc);
    sleipnerAMG(cc, glLinOp, *parComm, A_loc, rhs_loc);
    sleipnerCPR(cc, glLinOp, *parComm, A_loc, rhs_loc);
}

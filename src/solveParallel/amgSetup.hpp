/*
  Copyright 2019 Andreas Thune.

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

#ifndef OPM_AMGSETUP_HEADER_INCLUDED
#define OPM_AMGSETUP_HEADER_INCLUDED

#endif // OPM_AMGSETUP_HEADER_INCLUDED


namespace Dune
{
  namespace Amg
  {



      
  } //AMG
} //Dune

/**
 * @brief A policy class for solving the coarse level system using one step of AMG.
 * @tparam O The type of the linear operator used.
 * @tparam S The type of the smoother used in AMG.
 * @tparam C The type of the crition used for the aggregation within AMG.
 */
template<class O, class S, class C, class Comm, class L>
class OneStepAMGCoarseSolverPolicyCpr
{
public:
  /** @brief The type of the linear operator used. */
  typedef O Operator;
  /** @brief The type of the range and domain of the operator. */
  typedef typename O::range_type X;
  /** @brief The type of the crition used for the aggregation within AMG.*/
  typedef C Criterion;
  /** @brief The type of the smoother used in AMG. */
  typedef S Smoother;
  /** @brief The type of the arguments used for constructing the smoother. */
  typedef typename Dune::Amg::SmootherTraits<S>::Arguments SmootherArgs;
  /** @brief The type of the AMG construct on the coarse level.*/
  typedef Dune::Amg::AMGCPR<Operator,X,Smoother,Comm> AMGType;
  typedef L LTP;
  typedef Comm Communication;
  /**
   * @brief Constructs the coarse solver policy.
   * @param args The arguments used for constructing the smoother.
   * @param c The crition used for the aggregation within AMG.
   */
  OneStepAMGCoarseSolverPolicyCpr(const SmootherArgs& args, const Criterion& c)
    : smootherArgs_(args), criterion_(c)
  {}
  /** @brief Copy constructor. */
  OneStepAMGCoarseSolverPolicyCpr(const OneStepAMGCoarseSolverPolicyCpr& other)
  : coarseOperator_(other.coarseOperator_), smootherArgs_(other.smootherArgs_),
    criterion_(other.criterion_)
  {}
private:
  /**
   * @brief A wrapper that makes an inverse operator out of AMG.
   *
   * The operator will use one step of AMG to approximately solve
   * the coarse level system.
   */
    struct AMGInverseOperator : public Dune::InverseOperator<X,X>
  {
    AMGInverseOperator(const typename AMGType::Operator& op,
                       const Criterion& crit,
                       const typename AMGType::SmootherArgs& args)
      : first_(true)
    {
	amg_ = std::make_unique<AMGType> (op, crit, args);
    }

    AMGInverseOperator(const typename AMGType::Operator& op,
                       const Criterion& crit,
                       const typename AMGType::SmootherArgs& args,
		       const Communication& comm)
      : first_(true)
    {
	amg_ = std::make_unique<AMGType> (op, crit, args, comm);
    }

    void apply(X& x, X& b, double reduction, Dune::InverseOperatorResult& res)
    {
      DUNE_UNUSED_PARAMETER(reduction);
      DUNE_UNUSED_PARAMETER(res);
      if(first_)
      {
        amg_->pre(x,b);
        first_=false;
        x_=x;
      }
      amg_->apply(x,b);
    }

    void apply(X& x, X& b, Dune::InverseOperatorResult& res)
    {
      return apply(x,b,1e-8,res);
    }

    virtual Dune::SolverCategory::Category category() const
    {
      return amg_->category();
    }
    
    ~AMGInverseOperator()
    {
      if(!first_)
        amg_->post(x_);
    }
    AMGInverseOperator(const AMGInverseOperator& other)
    : x_(other.x_), amg_(other.amg_), first_(other.first_)
    {
    }
    AMGType getCoarseAMG() const
    {
	return *amg_;
    }
  private:
    X x_;
    std::unique_ptr<AMGType> amg_;
    bool first_;
  };
public:
  /** @brief The type of solver constructed for the coarse level. */
  typedef AMGInverseOperator CoarseLevelSolver;

  /**
    * @brief Constructs a coarse level solver.
    *
    * @param transferPolicy The policy describing the transfer between levels.
    * @return A pointer to the constructed coarse level solver.
    * @tparam P The type of the level transfer policy.
    */
  template<class P>
  CoarseLevelSolver* createCoarseLevelSolver(P& transferPolicy)
  {
    coarseOperator_=transferPolicy.getCoarseLevelOperator();
    auto& tp = dynamic_cast<LTP&>(transferPolicy);
    AMGInverseOperator* inv = new AMGInverseOperator(*coarseOperator_,
                                                     criterion_,
                                                     smootherArgs_,
						     tp.getCoarseLevelCommunication());

    return inv; //std::shared_ptr<InverseOperator<X,X> >(inv);

  }

private:
  /** @brief The coarse level operator. */
  std::shared_ptr<Operator> coarseOperator_;
  /** @brief The arguments used to construct the smoother. */
  SmootherArgs smootherArgs_;
  /** @brief The coarsening criterion. */
  Criterion criterion_;
};

template<class Crit>
void setCrit(Crit& criterion, Opm::PropertyTree prm_amg)
{
    criterion.setDefaultValuesIsotropic(2);
    criterion.setAlpha(prm_amg.get<double>("alpha", 0.33));
    criterion.setBeta(prm_amg.get<double>("beta", 1e-5));
    criterion.setMaxLevel(prm_amg.get<int>("maxlevel", 15));
    criterion.setSkipIsolated(prm_amg.get<bool>("skip_isolated", false));
    criterion.setNoPreSmoothSteps(prm_amg.get<int>("pre_smooth", 1));
    criterion.setNoPostSmoothSteps(prm_amg.get<int>("post_smooth", 1));
    criterion.setDebugLevel(10);//prm_amg.get<int>("verbosity", 10));
    criterion.setAccumulate(static_cast<Dune::Amg::AccumulationMode>(prm_amg.get<int>("accumulate", 1)));
    criterion.setProlongationDampingFactor(prm_amg.get<double>("prolongationdamping", 1.6));
    criterion.setMaxDistance(prm_amg.get<int>("maxdistance", 2));
    criterion.setMaxConnectivity(prm_amg.get<int>("maxconnectivity", 15));
    criterion.setMaxAggregateSize(prm_amg.get<int>("maxaggsize", 6));
    criterion.setMinAggregateSize(prm_amg.get<int>("minaggsize", 4));
    criterion.setRandomParallelGhostIndexOrder(prm_amg.get<bool>("random_coarse_ghost_index", true));
}

template<class SA>
void setOpmILU0args(SA& smootherArgs, Opm::PropertyTree prm_amg)
{
    smootherArgs.iterations = prm_amg.get<int>("iterations", 1);
    const int iluwitdh = prm_amg.get<int>("iluwidth", 0);
    smootherArgs.setN(iluwitdh);
    const Opm::MILU_VARIANT milu = Opm::convertString2Milu(prm_amg.get<std::string>("milutype", std::string("ilu")));
    smootherArgs.setMilu(milu);
    smootherArgs.relaxationFactor = prm_amg.get<double>("relaxation", 1.0);
}

template<class SA>
void setOpmILU0Defargs(SA& smootherArgs)
{
    smootherArgs.iterations = 1;
    const int iluwitdh = 0;
    smootherArgs.setN(iluwitdh);
    const Opm::MILU_VARIANT milu = Opm::convertString2Milu(std::string("ilu"));
    smootherArgs.setMilu(milu);
    smootherArgs.relaxationFactor = 1.0;
}

template<class SA>
void setOpmILU0noMILU(SA& smootherArgs, Opm::PropertyTree prm_amg)
{
    smootherArgs.iterations = prm_amg.get<int>("iterations", 1);
    const int iluwitdh = prm_amg.get<int>("iluwidth", 0);
    smootherArgs.setN(iluwitdh);
    smootherArgs.relaxationFactor = prm_amg.get<double>("relaxation", 1.0);
}

/*
    "maxiter": "10",
    "tol": "1e-3",
    "verbosity": "3",
    "solver": "bicgstab",
    "preconditioner": {
        "type": "cpr",
        "weight_type": "quasiimpes",
        "finesmoother": {
            "type": "ParOverILU0",
            "relaxation": "1"
        },
        "pre_smooth": "1",
        "post_smooth": "1",
        "pressure_var_index": "0",
        "verbosity": "10",
        "coarsesolver": {
            "maxiter": "1",
            "tol": "0.10000000000000001",
            "solver": "loopsolver",
            "verbosity": "0",
            "preconditioner": {
                "type": "amg",
                "alpha": "0.33333333333300003",
                "relaxation": "1",
                "iterations": "1",
                "coarsenTarget": "1200",
                "pre_smooth": "1",
                "post_smooth": "1",
                "beta": "0",
                "smoother": "ILU0",
                "verbosity": "10",
                "maxlevel": "15",
                "skip_isolated": "0",
                "accumulate": "1",
                "prolongationdamping": "1.0",
                "maxdistance": "2",
                "maxconnectivity": "15",
                "maxaggsize": "6",
                "minaggsize": "4"
            }
        }
    }

 */
template<class Crit>
void setCritCPR(Crit& criterion, Opm::PropertyTree prm_amg)
{
    criterion.setDefaultValuesIsotropic(2);
    criterion.setAlpha(prm_amg.get<double>("preconditioner.coarsesolver.preconditioner.alpha", 0.33));
    criterion.setBeta(prm_amg.get<double>("preconditioner.coarsesolver.preconditioner.beta", 1e-5));
    criterion.setMaxLevel(prm_amg.get<int>("preconditioner.coarsesolver.preconditioner.maxlevel", 15));
    criterion.setSkipIsolated(prm_amg.get<bool>("preconditioner.coarsesolver.preconditioner.skip_isolated", false));
    criterion.setNoPreSmoothSteps(prm_amg.get<int>("preconditioner.coarsesolver.preconditioner.pre_smooth", 1));
    criterion.setNoPostSmoothSteps(prm_amg.get<int>("preconditioner.coarsesolver.preconditioner.post_smooth", 1));
    criterion.setDebugLevel(10);//prm_amg.get<int>("verbosity", 10));
    criterion.setAccumulate(static_cast<Dune::Amg::AccumulationMode>(prm_amg.get<int>("preconditioner.coarsesolver.preconditioner.accumulate", 1)));
    criterion.setProlongationDampingFactor(prm_amg.get<double>("preconditioner.coarsesolver.preconditioner.prolongationdamping", 1.6));
    criterion.setMaxDistance(prm_amg.get<int>("preconditioner.coarsesolver.preconditioner.maxdistance", 2));
    criterion.setMaxConnectivity(prm_amg.get<int>("preconditioner.coarsesolver.preconditioner.maxconnectivity", 15));
    criterion.setMaxAggregateSize(prm_amg.get<int>("preconditioner.coarsesolver.preconditioner.maxaggsize", 6));
    criterion.setMinAggregateSize(prm_amg.get<int>("preconditioner.coarsesolver.preconditioner.minaggsize", 4));
}

template<class Operator, class Vec, class Pre, class Comm, class C>
void timeLevel(const Operator& linOp, Vec& rhs_large, Pre& pre, const Comm& comm, const C& cc)
{
    int rank = cc.rank();
    size_t N = linOp.getmat().N();
    Vec rhs(N);
    for (size_t i = 0; i < N; ++i) {rhs[i]=rhs_large[i];}

    if (rank == 0) {std::cout << std::endl;}
    cc.barrier();
    multipleMinLoopTimeComm(cc, comm, rhs);

    if (rank == 0) {std::cout << std::endl;}
    cc.barrier();
    multipleMinLoopTimeProject(cc, comm, rhs, 5);

    if (rank == 0) {std::cout << std::endl;}
    cc.barrier();
    multipleMinLoopTimeSpMVAS(cc,linOp, rhs, 3);

    if (rank == 0) {std::cout << std::endl;}
    cc.barrier();
    multipleMinLoopTimePre(cc, pre, rhs, 3);
}

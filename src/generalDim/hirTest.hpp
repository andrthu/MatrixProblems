
template<class Mat, class Agg>
void buildPro(Mat& P, Mat& PT, const Agg& aggM, int N, int M) {
    Dune::MatrixIndexSet op, opT;
    op.resize(M, N);
    opT.resize(N, M);

    int aggCount=0;
    for (auto aggI = aggM.begin(); aggI!= aggM.end(); ++aggI, aggCount++) {
	op.add(*aggI,aggCount);
	opT.add(aggCount,*aggI);
    }

    op.exportIdx(P);
    opT.exportIdx(PT);

    typename Mat::block_type diag(0.0);
    for (int el = 0; el < diag.size(); el++)
	diag[el][el] = 1.0;
    aggCount=0;
    for (auto aggI = aggM.begin(); aggI!= aggM.end(); ++aggI, aggCount++) {

	P[*aggI][aggCount]=diag;
	PT[aggCount][*aggI]=diag;
    }
    
}

template<class Mat, class Agg>
void buildAggRestrict(Mat& R, const Mat& A, const Agg& aggM, int N) {
    Dune::MatrixIndexSet op;
    op.resize(N, N);

    for (auto row= A.begin(); row != A.end(); ++row) {
	
	
	int d = row.index();
	op.add(d,d);

	auto aggO = aggM[d];
	auto col = row->begin();
	for (; col != row->end(); ++col) {
	    int nab = col.index();
	    if (aggO==aggM[nab])
		if (nab!=d)
		    op.add(d,nab);
	}
    }
    op.exportIdx(R);

    for (auto row= R.begin(); row != R.end(); ++row) {

	int d = row.index();
	

	auto aggO = aggM[d];
	auto col = row->begin();
	for (; col != row->end(); ++col) {
	    int nab = col.index();
	    R[d][nab] = A[d][nab];
	}
    }
    
}

template<class Mat>
void jacobiS(Mat& S, double w)
{

    typename Mat::block_type diag(0.0);
    for (int el = 0; el < diag.size(); el++)
	diag[el][el] = 1 - w;
    
    for (auto row = S.begin(); row != S.end(); ++row) {

	int d = row.index();
	auto inv = S[d][d];
	inv.invert();
	inv*=-w;
	//= -w * S[d][d].invert();
	auto col = row->begin();
	S[d][d]=diag;
	for (; col != row->end(); ++col) {

	    int nab = col.index();
	    if (nab!=d) {
		S[d][nab] =  S[d][nab].leftmultiply(inv);
	    }
	    
		
	}
    }
    
}

template<class Mat, class Vec>
void gen_dim_hirTest(int argc, char** argv)
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

    Mat A_loc, A_loc_;
    Vec rhs, rhs_loc;

    DictRead DR;
    Comm comm(cc);
    std::shared_ptr<Comm> parComm(new(Comm));
    readMatOnRootAndDist(argc, argv, A_loc, rhs_loc, DR, comm, parComm, cc); // in partition/overlapCreation.hpp

    //findZeroDiag(A_loc, rhs_loc);
    ScalarProduct sp(*parComm);
    GLO linOp(A_loc, *parComm);
    Operator op(A_loc, *parComm);

    Opm::FlowLinearSolverParameters flsp_json;
    flsp_json.linsolver_ = DR.dict[12];

    Opm::PropertyTree prm_json(flsp_json.linsolver_);
    
    using CriterionBase
	= Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricDependency<Mat, Dune::Amg::FirstDiagonal>>;
    using Criterion = Dune::Amg::CoarsenCriterion<CriterionBase>;

    Criterion criterion(15, prm_json.get<int>("coarsenTarget", 1200));
    auto pc_child = prm_json.get_child_optional("preconditioner");
    setCrit(criterion, *pc_child );


    typedef std::allocator<Operator> Allocator;
    typedef Dune::Amg::Hierarchy<Comm,Allocator> ParallelInformationHierarchy;
    typedef typename ParallelInformationHierarchy::Iterator PInfoIterator;
    
    ParallelInformationHierarchy parallelInformation(parComm);
    PInfoIterator infoLevel = parallelInformation.finest();

    infoLevel->buildGlobalLookup(op.getmat().N());
    
    typedef Dune::Amg::PropertiesGraphCreator<Operator,Comm> GraphCreator;
    typedef typename GraphCreator::PropertiesGraph PropertiesGraph;
    typedef typename GraphCreator::GraphTuple GraphTuple;
    typedef typename PropertiesGraph::VertexDescriptor Vertex;

    std::vector<bool> excluded(op.getmat().N(), false);

    typedef Dune::NegateSet<typename Comm::OwnerSet> OverlapFlags;
    
    GraphTuple graphs = GraphCreator::create(op, excluded, *parComm, OverlapFlags());

    typedef Dune::Amg::AggregatesMap<Vertex> AggMap;
    AggMap* aggregatesMap=new AggMap(std::get<1>(graphs)->maxVertex()+1);

    
    auto [noAggregates, isoAggregates, oneAggregates, skippedAggregates] =
          aggregatesMap->buildAggregates(op.getmat(), *(std::get<1>(graphs)), criterion, false);

    Mat P,PT;
    buildPro(P, PT, *aggregatesMap, op.getmat().N(), noAggregates);
    std::cout<<" A NxM "<< A_loc.N() << " " << A_loc.M()<<" "<< A_loc.nonzeroes()<<std::endl;
    std::cout<<" P NxM "<< P.N() << " " << P.M()<<std::endl;
    std::cout<<" PT NxM "<< PT.N() << " " << PT.M()<<std::endl;

    Vec x(op.getmat().N());
    Vec y(noAggregates);
    x=1;

    P.mv(x,y);
    //PT.mv(y,x);

    Mat AR;
    buildAggRestrict(AR,A_loc,*aggregatesMap,op.getmat().N());
    std::cout<<" AR NxM "<< AR.N() << " " << AR.M()<<" "<< AR.nonzeroes()<< " nnz(RART)/nnz(A): "<< (float)AR.nonzeroes()/A_loc.nonzeroes() <<std::endl;
    
    Mat PA,PAPT;
    Dune::matMultMat(PA,P,A_loc);
    std::cout<<" PA NxM "<< PA.N() << " " << PA.M()<<" "<< PA.nonzeroes()<<std::endl;
    Dune::matMultMat(PAPT,PA,PT);

    std::cout<<" PAPT NxM "<< PAPT.N() << " " << PAPT.M()<<" "<<PAPT.nonzeroes()<<std::endl;
    std::cout<<std::endl;
    //Mat S(A_loc);
    Mat S(AR);
    jacobiS(S, 0.9);
    std::cout<<" S NxM "<< S.N() << " " << S.M()<<" "<< S.nonzeroes()<<std::endl;

    Mat PS,PST;
    Dune::matMultMat(PS,P,S);
    Dune::matMultMat(PST,S,PT);

    std::cout<<" PS NxM "<< PS.N() << " " << PS.M()<<" "<< PS.nonzeroes()<<std::endl;
    std::cout<<" PST NxM "<< PST.N() << " " << PST.M()<<" "<< PST.nonzeroes()<<std::endl;

    Mat RA, RART;
    Dune::matMultMat(RA,PS,A_loc);
    std::cout<<" RA NxM "<< RA.N() << " " << RA.M()<<" "<< RA.nonzeroes()<<std::endl;
    Dune::matMultMat(RART,RA,PST);
    std::cout<<" RART NxM "<< RART.N() << " " << RART.M()<<" "<< RART.nonzeroes()<< " nnz(RART)/nnz(A): "<< (float)RART.nonzeroes()/A_loc.nonzeroes() <<std::endl;
    
    int aggCount = 0;
    //for (auto aggI = aggregatesMap->begin(); aggI!= aggregatesMap->end(); ++aggI, aggCount++)
    //if(rank==0 )
    //std::cout<<" Agg value "<< aggCount << " " << *aggI<< " "<< y[*aggI][0]<< std::endl;

    
    if(rank==0 )
	std::cout<<" Have built "<<noAggregates
		 <<" aggregates totally ("<<isoAggregates<<" isolated aggregates, "
		 << oneAggregates<<" aggregates of one vertex,  and skipped "
		 << skippedAggregates<<" aggregates)."<<std::endl;

    typedef typename Dune::Amg::ConstructionTraits<Comm>::Arguments CommunicationArgs;
    CommunicationArgs commargs(parComm->communicator(),parComm->category());
    parallelInformation.addCoarser(commargs);
    ++infoLevel;
    
    typedef Dune::Amg::VertexVisitedTag VertVetTag;
    typename Dune::PropertyMapTypeSelector<VertVetTag,PropertiesGraph>::Type visitedMap =
	get(VertVetTag(), *(std::get<1>(graphs)));

    
    int aggregates = Dune::Amg::IndicesCoarsener<Comm,OverlapFlags>
	::coarsen(*parComm,
		  *(std::get<1>(graphs)),
		  visitedMap,
		  *aggregatesMap,
		  *infoLevel,
		  noAggregates);

    std::cout<<"Aggregates: " << aggregates <<std::endl;

    infoLevel->buildGlobalLookup(aggregates);

    Dune::Amg::AggregatesPublisher<Vertex,OverlapFlags,Comm>::publish(*aggregatesMap,
								      *parComm,
								      infoLevel->globalLookup());

    std::vector<bool>& visited=excluded;

    typedef std::vector<bool>::iterator Iterator;
    typedef Dune::IteratorPropertyMap<Iterator, Dune::IdentityMap> VisitedMap2;
    Iterator end = visited.end();
    for(Iterator iter= visited.begin(); iter != end; ++iter)
	*iter=false;

    VisitedMap2 visitedMap2(visited.begin(), Dune::IdentityMap());
    
    Dune::Amg::GalerkinProduct<Comm> productBuilder;
    std::shared_ptr<typename Operator::matrix_type>
	coarseMatrix(productBuilder.build(*(std::get<0>(graphs)), visitedMap2,
					  *parComm,
					  *aggregatesMap,
					  aggregates,
					  OverlapFlags()));

    productBuilder.calculate(op.getmat(),*aggregatesMap, *coarseMatrix, *infoLevel, OverlapFlags());

    
    std::cout<<" CM NxM "<< coarseMatrix->N() << " " << coarseMatrix->M()<<" "<< coarseMatrix->nonzeroes()<<std::endl;

}

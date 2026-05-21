/*
  Copyright 2026 Andreas Thune.

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

template <typename MatrixType>
void createNegativeAbsRowSumMatrix(const MatrixType& mat1, 
                                   MatrixType& mat2) 
{
    const std::size_t rows = mat1.N();
    const std::size_t cols = mat1.M();

    // 1. Build the Sparsity Pattern
    Dune::MatrixIndexSet indices(rows, cols);
    for (auto i = mat1.begin(); i != mat1.end(); ++i) {
        for (auto j = i->begin(); j != i->end(); ++j) {
            indices.add(i.index(), j.index());
        }
    }

    // 2. Add the diagonal indices (ensures they exist in mat2)
    for (std::size_t i = 0; i < rows; ++i) {
        indices.add(i, i);
    }

    // 3. Apply the pattern to mat2 and initialize values to zero
    indices.exportIdx(mat2);
    mat2 = 0.0;

    // 4. Populate mat2 using mat1's data
    for (auto i = mat1.begin(); i != mat1.end(); ++i) {
        auto rowIdx = i.index();
        double rowSum = 0.0;

        for (auto j = i->begin(); j != i->end(); ++j) {
            auto colIdx = j.index();
            
            // Extract scalar value from 1x1 block
            double val = (*j)[0][0];
            
            // Set mat2 off-diagonal to -mat1
            // Note: If colIdx == rowIdx (unlikely per your description), 
            // this will be overwritten by the rowSum later anyway.
            mat2[rowIdx][colIdx][0][0] = -val;
            
            // Accumulate absolute values for the diagonal
            rowSum += std::abs(val);
        }

        // 5. Assign the calculated sum to the diagonal of mat2
        mat2[rowIdx][rowIdx][0][0] = rowSum;
    }
}

template<class Mat, class Hir>
void writeMergedTrans(Mat mat, const Hir& hir)
{
    //auto aggMap = hir->aggregatesMaps();
    //auto mats = hir->matrices();
    auto aggMap = hir.aggregatesMaps();
    auto mats = hir.matrices();
    auto aggA = aggMap.begin();//opHir1->aggregatesMaps().begin();
    int t = 0;
    std::vector<std::vector<int>> maps;
    std::vector<std::vector<int>> cum;
    std::vector<int> aggSize;
    auto cmat = mats.finest();

    typedef typename Hir::AggregatesMap AMap;
    
    for (auto aggA = aggMap.begin(); aggA!=aggMap.end(); aggA++) {
        std::vector<int> f2c;
	
	//aggSize.push_back((*aggA)->end()-(*aggA)->begin());
	int numIso = 0;
	int y =0;
	//std::cout << (*aggA)->end()-(*aggA)->begin() << " " << cmat->getmat().N() <<std::endl;
	for (auto aa= (*aggA)->begin(); aa!=(*aggA)->end(); aa++,y++) {
	    //std::cout << *aa << " " << t<<std::endl;
	    if (*aa == AMap::ISOLATED) {
		//std::cout << *aa << " " << t<<std::endl;
		f2c.push_back(-1);
		numIso++;
	    } else {
		f2c.push_back(*aa);
	    }
	}
	maps.push_back(f2c);
	if (t==0) {
	    cmat++;
	    std::cout << f2c.size() << " " << cmat->getmat().N() <<std::endl;
	    aggSize.push_back(cmat->getmat().N());
	    cum.push_back(f2c);
	} else {
	    std::vector<int> pre = cum[t-1];
	    std::vector<int> nc(mat.N(), 0);

	    if (f2c.size() > 0) {
		cmat++;
	        std::cout << f2c.size() << " " << cmat->getmat().N() <<std::endl;
		aggSize.push_back(cmat->getmat().N());
	        for (int i = 0; i < mat.N(); ++i) {
		    //std::cout << t<< " " <<i << " " <<  pre[i] << " "<< f2c[pre[i]] <<std::endl;
		    if (pre[i] == -1) 
			nc[i] = -1;
		    else
			nc[i] = f2c[pre[i]];
		}
		
		cum.push_back(nc);
	    }
	}
	std::cout << "Number of Isolated: "<<numIso << std::endl;
	t++;
    }
    std::cout << "End create map" << std::endl;
    /*
    for (int f = 0; f < cum.size(); f++) {
        std::string filename = "mergedTrans_" + std::to_string(f) + ".txt";
	std::ofstream outputFileC(filename);
	std::vector<int> f2c = cum[f];
	for (auto i = mat.begin(); i != mat.end(); ++i) {

	    auto rowIdx = i.index();

	    for (auto j = i->begin(); j != i->end(); ++j) {
	        auto colIdx = j.index();

		if (rowIdx != colIdx) {

		    if (f2c[rowIdx] == f2c[colIdx]) {
		        outputFileC << rowIdx << " " << colIdx << " " << std::abs(*j) <<std::endl; 
		    }
		}
	    }
	}
	
	outputFileC.close();
    }
    */

    for (int f = 0; f < cum.size(); f++) {

	std::vector<std::map<int, double> > edges(aggSize[f]);
	std::vector<std::vector<int>> nodes(aggSize[f]);
	
	for (auto i = mat.begin(); i != mat.end(); ++i) {

	    auto rowIdx = i.index();

	    auto cidx = cum[f][rowIdx];
	    if (cidx != -1) {
		//std::cout <<rowIdx<<" "<<cidx<<std::endl;
		nodes[cidx].push_back(rowIdx);
	    
		for (auto j = i->begin(); j != i->end(); ++j) {
		    auto colIdx = j.index();
		    auto cnab = cum[f][colIdx];

		    if ( edges[cidx].count(cnab) == 1 ) {
			edges[cidx][cnab] += 1;
		    } else {
			edges[cidx].insert({cnab,1});
		    }
		}
	    }
	}
	int nsize = 0;
	for (int x = 0; x < aggSize[f]; x++) {

	    if (nodes[x].size() > nsize)
		nsize = nodes[x].size();
	}

	std::cout << "Max node size: "<< nsize<< std::endl;
    }
}


template<class Mat, class Vec>
void gen_dim_transmissibility_hir_part(std::vector<std::string> systemDirs)
{
    typedef Dune::FieldVector<double,1> BlockVec1;
    typedef Dune::BlockVector<BlockVec1> VecT;
    typedef Dune::FieldMatrix<double,1,1> BlockMat1;
    typedef Dune::BCRSMatrix<BlockMat1> MatT;
    
    typedef Dune::MPIHelper::MPICommunicator MPICommunicator;
    typedef Dune::Communication<MPICommunicator> CollectiveCommunication;
    typedef Dune::BiCGSTABSolver<Vec> Solver;
    typedef Dune::InverseOperatorResult Stat;
    
    typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
    typedef Dune::OverlappingSchwarzScalarProduct<Vec,Comm> ScalarProduct;
    typedef GhostLastMatrixAdapter<Mat,Vec,Vec,Comm> GLO;                 // solveParallel/ghostLastOperations.hpp
    typedef GhostLastMatrixAdapter<MatT,VecT,VecT,Comm> GLOT;
    typedef Dune::OverlappingSchwarzOperator<Mat,Vec,Vec,Comm> Operator;
    typedef Opm::ParallelOverlappingILU0<Mat,Vec,Vec,Comm> ILU;
    typedef Opm::ParallelOverlappingILU0<MatT,VecT,VecT,Comm> ILUT;
    typedef Dune::Amg::AMGCPR<GLOT, VecT, ILUT, Comm> AMGCPR;
    typedef Dune::FlexibleSolver<GLO> FlexibleSolverType;

    using Smoother = ILUT;
    using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
    using CriterionBase
	= Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricDependency<MatT,Dune::Amg::FirstDiagonal>>;
    using Criterion = Dune::Amg::CoarsenCriterion<CriterionBase>;


    const auto block_size = Vec::block_type::dimension;
    
    CollectiveCommunication cc(MPI_COMM_WORLD);
    int rank = cc.rank();

    std::vector<Mat> systems;
    std::vector<Vec> rhs;
    std::vector<MatT> transMats;

    std::vector<ScalarProduct> sps;
    
    DictRead DR;
    Comm comm(cc);
    std::shared_ptr<Comm> parComm(new(Comm));
    std::vector<int> mpiVec;

    for (int i = 0; i < systemDirs.size(); ++i) {
	if ( boost::algorithm::ends_with( systemDirs[i], ".ini") ) {
	    
	    DR.read_file_and_update(systemDirs[i].data());
	    if (rank == 0) {
		DR.write_param();
	    }
	}
    }
    
    for (int i = 0; i < systemDirs.size(); ++i) {

	if ( boost::algorithm::ends_with( systemDirs[i], ".json") ) {
	    DR.dict[12] = systemDirs[i];
	} else if ( boost::algorithm::ends_with( systemDirs[i], ".ini") ) {}
	else {

	    Mat A_loc;
	    MatT trans;
	    Vec rhs_loc;
	    readTransMatOnly(trans, systemDirs[i], rank);

	    MatT poisson;
	    createNegativeAbsRowSumMatrix(trans, poisson);
	    
	    
	    mpiVec = readMatOnRootAndDist(systemDirs[i], A_loc, rhs_loc, DR, comm, parComm, cc, mpiVec, true, i!=0);

	    transMats.push_back(poisson);
	    systems.push_back(A_loc);
	    rhs.push_back(rhs_loc);
	    sps.push_back(ScalarProduct(*parComm));
	}
    }

    if (rank == 0) {std::cout << std::endl;}

    GLOT glLinOp(transMats[0], *parComm);
	    
    Opm::FlowLinearSolverParameters flsp_amg;
    
    Opm::PropertyTree prm_amg = setupAMG(std::string("amg"), flsp_amg);
    prm_amg.put("preconditioner.verbosity", 10);
    prm_amg.put("skip_isolated",false);
    Criterion criterion(15, prm_amg.get<int>("coarsenTarget", 1200));
    setCrit(criterion, prm_amg);
    SmootherArgs smootherArgs;
    setOpmILU0args(smootherArgs, prm_amg);

    auto amg = std::make_shared<AMGCPR>(glLinOp, criterion, smootherArgs, *parComm);

    std::vector<std::size_t> ag;
    amg->getCoarsestAggregateNumbers(ag);
    auto opHir1 = amg->operatorHirarchyList();
    auto mats = opHir1->matrices();
    auto cmat = mats.coarsest()->getmat();

    writeMergedTrans(transMats[0], *opHir1);
}

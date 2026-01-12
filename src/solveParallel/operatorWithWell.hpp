/*
  Copyright 2025 Andreas Thune.

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

#ifndef OPM_OPERATORWITHWELL_HEADER_INCLUDED
#define OPM_OPERATORWITHWELL_HEADER_INCLUDED

#endif

template<class M, class X, class Y, class W, class MSW, class C>
class GhostLastMatrixWellAdapter : public Dune::AssembledLinearOperator<M,X,Y>
{
public:
    typedef M matrix_type;
    typedef X domain_type;
    typedef Y range_type;
    typedef W WellModel;
    typedef MSW MultiSegmentWellModel;
    typedef typename X::field_type field_type;

    using PressureMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<field_type, 1, 1>>;
    typedef C communication_type;

    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::overlapping;
    }

    //! constructor: just store a reference to a matrix
    GhostLastMatrixWellAdapter (const M& A,
				const std::vector<W>& wells,
				const std::vector<MSW>& mswells,
				const communication_type& comm)
        : A_( Dune::stackobject_to_shared_ptr(A) ),wells_(wells),mswells_(mswells),comm_(comm)
    {
        interiorSize_ = setInteriorSize(comm_);
    }

    GhostLastMatrixWellAdapter (const std::shared_ptr<M> A,
				std::vector<W>& wells,
				const std::vector<MSW>& mswells,
				const communication_type& comm)
        : A_( A ), wells_(wells), mswells_(mswells), comm_(comm)
    {
        interiorSize_ = setInteriorSize(comm_);	
    }

    virtual void apply( const X& x, Y& y ) const override
    {
        for (auto row = A_->begin(); row.index() < interiorSize_; ++row)
        {
            y[row.index()]=0;
            auto endc = (*row).end();
            for (auto col = (*row).begin(); col != endc; ++col)
                (*col).umv(x[col.index()], y[row.index()]);
        }
	for (auto & well : wells_) {
	    well.apply(x,y);
	}
	for (auto & mwell : mswells_) {
	    mwell.apply(x,y);
	}

        ghostLastProject( y );
    }

    // y += \alpha * A * x
    virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const override
    {
	//if (comm_.communicator().rank() == 0) {std::cout << "Do AMG/CPR SpMV"<< std::endl;}
        for (auto row = A_->begin(); row.index() < interiorSize_; ++row)
        {
            auto endc = (*row).end();
            for (auto col = (*row).begin(); col != endc; ++col)
                (*col).usmv(alpha, x[col.index()], y[row.index()]);
        }
	for (auto & well : wells_) {
	    well.applyscaleadd(alpha, x, y);
	}
	for (auto & mwell : mswells_) {
	    mwell.applyscaleadd(alpha, x, y);
	}
        ghostLastProject( y );
    }

    virtual const matrix_type& getmat() const override { return *A_; }

    const communication_type& comm() { return comm_; }

    int getNumberOfExtraEquations() const
    {
	return wells_.size() + mswells_.size();
    }

    void addWellPressureEquationsStruct(PressureMatrix& jacobian) const
    {
	int nw = getNumberOfExtraEquations();
	int rdofs = A_->N();
	

	for(int i=0; i < nw; i++){
	    int wdof = rdofs + i;
	    jacobian.entry(wdof,wdof) = 1.0;
	}

	int wnum = 0;
	for (auto & well : wells_) {

	    auto wc = well.wellConnections();
	    for(int perfcell : wc) {
		int wdof = rdofs + wnum;
		jacobian.entry(wdof,perfcell) = 0.0;
		jacobian.entry(perfcell, wdof) = 0.0;
	    }
	    wnum++;
	}

	for (auto & mwell : mswells_) {

	    auto wc = mwell.wellConnections();
	    for(int perfcell : wc) {
		int wdof = rdofs + wnum;
		jacobian.entry(wdof,perfcell) = 0.0;
		jacobian.entry(perfcell, wdof) = 0.0;
	    }
	    wnum++;
	}
    }
    
    void addWellPressureEquations(PressureMatrix& jacobian,
				  const X& weights,
				  const bool use_well_weights) const
    {
	int nw = getNumberOfExtraEquations();
	int rdofs = A_->N();

	for(int i=0; i < nw; i++){
	    int wdof = rdofs + i;
	    jacobian[wdof][wdof] = 1.0;
	}
	int widx = 0;
	for (auto & well : wells_) {

	    well.addWellPressureEquations(jacobian, weights, use_well_weights, widx, 1);
	    widx++;
	}
    }
    
private:
    void ghostLastProject(Y& y) const
    {
        size_t end = y.size();
        for (size_t i = interiorSize_; i < end; ++i)
            y[i] = 0;
    }

    size_t setInteriorSize(const communication_type& comm) const
    {
        auto indexSet = comm.indexSet();
	auto rowIt = A_->end();
        size_t is = 0;
        for (auto idx = indexSet.begin(); idx!=indexSet.end(); ++idx) {

            if (idx->local().attribute()==1) {
                auto loc = idx->local().local();
                if (loc > is) {
                    is = loc;
                }
            }
        }
        return is + 1;
    }
    const std::shared_ptr<const matrix_type> A_ ;
    const communication_type&  comm_;
    const std::vector<WellModel>& wells_;
    const std::vector<MultiSegmentWellModel> mswells_;
    size_t interiorSize_;
    typename matrix_type::RowIterator endRow_;
};


namespace Dune {
    namespace Amg {

        template<class M, class X, class Y, class W, class MSW, class C>
	class ConstructionTraits<GhostLastMatrixWellAdapter<M,X,Y,W,MSW,C> >
	{
	public:
	    typedef ParallelOperatorArgs<M,C> Arguments;

   	    static inline std::shared_ptr<GhostLastMatrixWellAdapter<M,X,Y,W,MSW,C>> construct(const Arguments& args)
	    {
	        return std::make_shared<GhostLastMatrixWellAdapter<M,X,Y,W,MSW,C>>
		    (args.matrix_, args.comm_);
	    }
	};

    } // end namespace Amg
} // end namespace Dune

template <class X, class Y, class W, class MSW>
class WellModelsFromFileOperator : public Opm::LinearOperatorExtra<X, Y>
{
public:
    using Base = Opm::LinearOperatorExtra<X, Y>;
    using field_type = typename Base::field_type;
    using PressureMatrix = typename Base::PressureMatrix;

    typedef W WellModel;
    typedef MSW MultiSegmentWellModel;
    
    explicit WellModelsFromFileOperator(const std::vector<WellModel>& wells,
					const std::vector<MultiSegmentWellModel>& mswells)
	: wells_(wells), mswells_(mswells)
    {}

    void apply(const X& x, Y& y) const override
    {
	for (auto & well : wells_) {
	    well.apply(x,y);
	}
	for (auto & mwell : mswells_) {
	    mwell.apply(x,y);
	}
    }

    void applyscaleadd(field_type alpha, const X& x, Y& y) const override
    {
	for (auto & well : wells_) {
	    well.applyscaleadd(alpha, x, y);
	}
	for (auto & mwell : mswells_) {
	    mwell.applyscaleadd(alpha, x, y);
	}
    }

    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

    int getNumberOfExtraEquations() const override
    {
	return wells_.size() + mswells_.size();
    }

    void addWellPressureEquationsStruct(PressureMatrix& jacobian) const override
    {
	int nw = getNumberOfExtraEquations();
	int rdofs = jacobian.N() - nw;

	for(int i=0; i < nw; i++){

	    int wdof = rdofs + i;
	    //std::cout << "entry loop " << i << " "<< wdof << std::endl; 
	    jacobian.entry(wdof,wdof) = 1.0;
	}

	
	int wnum = 0;
	for (auto & well : wells_) {

	    auto wc = well.wellConnections();
	    for(int perfcell : wc) {
		int wdof = rdofs + wnum;
		jacobian.entry(wdof,perfcell) = 0.0;
		jacobian.entry(perfcell, wdof) = 0.0;
	    }
	    wnum++;
	}
	

	for (int mwi = 0; mwi < mswells_.size(); ++mwi) {


	    const std::vector<int>& mwc = mswells_[mwi].wellConnections();

	    for(int perfcell : mswells_[mwi].wellConnections()) {
		int wdof = rdofs + wnum;

		jacobian.entry(wdof,perfcell) = 0.0;
		jacobian.entry(perfcell, wdof) = 0.0;
	    }
	    wnum++;
	}

	//std::cout << "Dofs in jacobian: " << jacobian.N() << std::endl;
    }
    
    void addWellPressureEquations(PressureMatrix& jacobian,
				  const X& weights,
				  const bool use_well_weights) const override
    {
	int nw = getNumberOfExtraEquations();
	int rdofs =jacobian.N()-nw;

	for(int i=0; i < nw; i++){
	    int wdof = rdofs + i;
	    jacobian[wdof][wdof] = 1.0;
	}

	int widx = 0;
	for (auto & well : wells_) {

	    well.addWellPressureEquations(jacobian, weights, use_well_weights, widx, 1);
	    widx++;
	}

	for (auto & mwell : mswells_) {

	    mwell.addWellPressureEquations(jacobian, weights, use_well_weights, widx, 1);
	    widx++;
	}
    }
    
private:
    const std::vector<WellModel>& wells_;
    const std::vector<MultiSegmentWellModel>& mswells_;
};

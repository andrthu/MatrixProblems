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

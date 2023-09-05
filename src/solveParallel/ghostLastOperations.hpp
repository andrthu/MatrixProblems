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

#ifndef OPM_GHOSTLASTOPERATIONS_HEADER_INCLUDED
#define OPM_GHOSTLASTOPERATIONS_HEADER_INCLUDED

#endif // OPM_GHOSTLASTOPERATIONS_HEADER_INCLUDED

template<class X, class C>
class GhostLastScalarProduct : public Dune::ScalarProduct<X>
{
public:
    typedef X domain_type;
    typedef typename X::field_type field_type;
    typedef typename Dune::FieldTraits<field_type>::real_type real_type;
    
    typedef C communication_type;
    
    enum {category=Dune::SolverCategory::overlapping};

    GhostLastScalarProduct (const communication_type& com, unsigned interiorSize)
	: cc(com), interiorSize_(interiorSize)
    {}
    virtual field_type dot (const X& x, const X& y)
    {
	field_type result = 0;
	for (unsigned i = 0; i < interiorSize_; ++i)
	    result += x[i]*(y[i]);
	return cc.sum(result);
    }
    
    virtual real_type norm (const X& x)
    {
	return std::sqrt(dot(x,x));
    }
    
private:
    const communication_type& cc;
    unsigned interiorSize_;
};

/*
template<class B>
class GhostLastBCRSMatrix : public Dune::BCRSMatrix<B>
{
public:
    typedef typename Dune::BCRSMatrix<B>::size_type size_type; 
    GhostLastBCRSMatrix() : Dune::BCRSMatrix<B>() {}
    GhostLastBCRSMatrix() : Dune::BCRSMatrix<B>()
}
*/


template<class M, class X, class Y>
class GhostLastOperator : public Dune::AssembledLinearOperator<M,X,Y>
{
public:

    typedef M matrix_type;
    typedef X domain_type;
    typedef Y range_type;
    typedef typename X::field_type field_type;
    typedef typename M::size_type size_type;
    typedef Dune::OwnerOverlapCopyCommunication<int,int> communication_type;
    
    enum {category=Dune::SolverCategory::overlapping};
    
    explicit GhostLastOperator (const M& A, size_type interiorSize, const communication_type &comm) : _A_(A), interiorSize_(interiorSize), comm_(comm){}
    
    //! apply operator to x:  \f$ y = A(x) \f$
    virtual void apply (const X& x, Y& y) const
    {
	y=0;

	auto first_row = _A_.begin();
	for (auto row = first_row; std::distance(first_row, row) < interiorSize_; ++row)
	{
	    auto endc = (*row).end();
	    for (auto col = (*row).begin(); col != endc; ++col)
		(*col).umv(x[col.index()], y[row.index()]);
	}
    }
    
    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
    {
	auto first_row = _A_.begin();
	for (auto row = first_row; std::distance(first_row, row) < interiorSize_; ++row)
	{
	    auto endc = (*row).end();
	    for (auto col = (*row).begin(); col != endc; ++col)
		(*col).usmv(alpha, x[col.index()], y[row.index()]);
	}
    }
    
    //! get matrix via *
    virtual const M& getmat () const
    {
	return _A_;
    }

    communication_type comm()
    {
	return comm_;
    }
    
private:
    const M& _A_;
    size_type interiorSize_;
    const communication_type& comm_;
};

template <class B, class A=std::allocator<B> >
class EasyCopyMatrix : public Dune::BCRSMatrix<B,A>
{
public:
    //===== type definitions and constants

    //! export the type representing the field
    typedef typename B::field_type field_type;

    //! export the type representing the components
    typedef B block_type;

    //! export the allocator type
    typedef A allocator_type;

    //! implement row_type with compressed vector
    //typedef BCRSMatrix<B,A>::row_type row_type;

    //! The type for the index access and the size
    typedef typename A::size_type size_type;

    //! increment block level counter
    enum {blocklevel = B::blocklevel+1};

    /** \brief Default constructor */
    EasyCopyMatrix() : Dune::BCRSMatrix<B,A>() {}

    explicit EasyCopyMatrix(int size)
      : Dune::BCRSMatrix<B,A>(size, size, Dune::BCRSMatrix<B,A>::random) {

      for (int i=0; i<size; i++)
        this->Dune::BCRSMatrix<B,A>::setrowsize(i, 1);

      this->Dune::BCRSMatrix<B,A>::endrowsizes();

      for (int i=0; i<size; i++)
        this->Dune::BCRSMatrix<B,A>::addindex(i, i);

      this->Dune::BCRSMatrix<B,A>::endindices();

    }

    //! assignment
    EasyCopyMatrix& operator= (const EasyCopyMatrix& other) {
      this->Dune::BCRSMatrix<B,A>::operator=(other);
      return *this;
    }

    //! assignment from scalar
    EasyCopyMatrix& operator= (const field_type& k) {
      this->Dune::BCRSMatrix<B,A>::operator=(k);
      return *this;
    }

    /** \brief Inverts the matrix */
    void invert() {
      for (int i=0; i<this->N(); i++)
        (*this)[i][i].invert();
    }
};


template<class M, class X, class Y, class C>
class GhostLastMatrixAdapter : public Dune::AssembledLinearOperator<M,X,Y>
{
public:
    typedef M matrix_type;
    typedef X domain_type;
    typedef Y range_type;
    typedef typename X::field_type field_type;


    typedef C communication_type;

    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::overlapping;
    }

    //! constructor: just store a reference to a matrix
    GhostLastMatrixAdapter (const M& A,
                            const communication_type& comm)
        : A_( Dune::stackobject_to_shared_ptr(A) ), comm_(comm)
    {
        interiorSize_ = setInteriorSize(comm_);
    }

    GhostLastMatrixAdapter (const std::shared_ptr<M> A,
                            const communication_type& comm)
        : A_( A ), comm_(comm)
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
    size_t interiorSize_;
    typename matrix_type::RowIterator endRow_;
};


namespace Dune {
    namespace Amg {

	template<class M, class X, class Y, class C>
	class ConstructionTraits<GhostLastMatrixAdapter<M,X,Y,C> >
	{
	public:
	    typedef ParallelOperatorArgs<M,C> Arguments;

	    static inline std::shared_ptr<GhostLastMatrixAdapter<M,X,Y,C>> construct(const Arguments& args)
	    {
		return std::make_shared<GhostLastMatrixAdapter<M,X,Y,C>>
		    (args.matrix_, args.comm_);
	    }
	};

    } // end namespace Amg
} // end namespace Dune


template<class M, class X, class Y, class C>
class OverlappingSchwarzOperatorCopy : public Dune::AssembledLinearOperator<M,X,Y>
{
public:

    typedef M matrix_type;
    typedef X domain_type;
    //! \brief The type of the range.
    //!
    //! E.g. BlockVector or another type fulfilling the ISTL
    //! vector interface.
    typedef Y range_type;
    //! \brief The field type of the range
    typedef typename X::field_type field_type;
    //! \brief The type of the communication object.
    //!
    //! This must either be OwnerOverlapCopyCommunication or a type
    //! implementing the same interface.
    typedef C communication_type;

    /**
     * @brief constructor: just store a reference to a matrix.
     *
     * @param A The assembled matrix.
     * @param com The communication object for syncing overlap and copy
     * data points. (E.~g. OwnerOverlapCopyCommunication )
     */
    OverlappingSchwarzOperatorCopy (const matrix_type& A, const communication_type& com)
      : _A_(stackobject_to_shared_ptr(A)), communication(com)
    {}

    OverlappingSchwarzOperatorCopy (const std::shared_ptr<matrix_type> A, const communication_type& com)
      : _A_(A), communication(com)
    {}

    //! apply operator to x:  \f$ y = A(x) \f$
    virtual void apply (const X& x, Y& y) const
    {
      y = 0;
      _A_->umv(x,y);     // result is consistent on interior+border
      communication.project(y);     // we want this here to avoid it before the preconditioner
                                    // since there d is const!
    }

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
    {
      _A_->usmv(alpha,x,y);     // result is consistent on interior+border
      communication.project(y);     // we want this here to avoid it before the preconditioner
                                    // since there d is const!
    }

    //! get the sequential assembled linear operator.
    virtual const matrix_type& getmat () const
    {
      return *_A_;
    }

    //! Category of the linear operator (see SolverCategory::Category)
    virtual Dune::SolverCategory::Category category() const
    {
	return Dune::SolverCategory::overlapping;
    }

private:
    const std::shared_ptr<const matrix_type>_A_;
    const communication_type& communication;
};

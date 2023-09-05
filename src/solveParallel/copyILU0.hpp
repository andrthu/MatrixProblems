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

#ifndef OPM_COPYILU0_HEADER_INCLUDED
#define OPM_COPYILU0_HEADER_INCLUDED

#endif // OPM_COPYILU0_HEADER_INCLUDED

#include <dune/common/version.hh>
#include <dune/istl/preconditioner.hh>
#include <dune/istl/ilu.hh>

template<class Matrix, class Domain, class Range, class ParallelInfo = Dune::Amg::SequentialInformation>
class ParallelOverlappingILU0;

template<class Matrix, class Domain, class Range, class ParallelInfo = Dune::Amg::SequentialInformation>
class GhostLastILU0;

template<class F>
class ParallelOverlappingILU0Args
    : public Dune::Amg::DefaultSmootherArgs<F>
{
 public:
    ParallelOverlappingILU0Args( )
        :  n_(0)
    {}
    
    void setN(int n)
    {
        n_ = n;
    }
    int getN() const
    {
        return n_;
    }
 private:
    int n_;
};

template<class M, class CRS, class InvVector>
void convertToCRS(const M& A, CRS& lower, CRS& upper, InvVector& inv , size_t interiorSize)
{
    // No need to do anything for 0 rows. Return to prevent indexing a
    // a zero sized array.
    if ( A.N() == 0 )
    {
	return;
    }

    typedef typename M :: size_type size_type;

    lower.clear();
    upper.clear();
    inv.clear();
    lower.resize( A.N() );
    upper.resize( A.N() );
    inv.resize( A.N() );

    // Count the lower and upper matrix entries.
    size_type numLower = 0;
    size_type numUpper = 0;
    const auto endi = A.end();
    for (auto i = A.begin(); i != endi; ++i) {
	const size_type iIndex = i.index();
	size_type numLowerRow = 0;
	for (auto j = (*i).begin(); j.index() < iIndex; ++j) {
	    ++numLowerRow;
	}
	numLower += numLowerRow;
	numUpper += (*i).size() - numLowerRow - 1;
    }
    assert(numLower + numUpper + A.N() == A.nonzeroes());

    lower.reserveAdditional( numLower );

    // implement left looking variant with stored inverse
    size_type row = 0;
    size_type colcount = 0;
    lower.rows_[ 0 ] = colcount;
    for (auto i=A.begin(); i.index() < interiorSize; ++i, ++row)
    {
	const size_type iIndex  = i.index();

	// eliminate entries left of diagonal; store L factor
	for (auto j=(*i).begin(); j.index() < iIndex; ++j )
	{
            lower.push_back( (*j), j.index() );
            ++colcount;
	}
	lower.rows_[ iIndex+1 ] = colcount;
    }
    for ( ; row < A.N(); ++row ) {
	lower.rows_[ row+1 ] = colcount;
    }

    const auto rendi = A.beforeBegin();
    row = 0;
    colcount = 0;
    upper.rows_[ 0 ] = colcount ;

    upper.reserveAdditional( numUpper );

    const int numEq = M :: block_type :: rows;
    typename M :: block_type diag_block(0.0);
    for (int eq = 0; eq < numEq; ++eq)
	diag_block[eq][eq] = 1.0;
    
    // NOTE: upper and inv store entries in reverse order, reverse here
    // relative to ILU
    for (auto i=A.beforeEnd(); i!=rendi; --i, ++ row )
    {
	const size_type iIndex = i.index();

	if (iIndex < interiorSize) {
	    // store in reverse row order
	    // eliminate entries left of diagonal; store L factor
	    for (auto j=(*i).beforeEnd(); j.index()>=iIndex; --j )
	    {
		const size_type jIndex = j.index();
		if( j.index() == iIndex )
	        {
		    inv[ row ] = (*j);
		    break;
		}
		else if ( j.index() >= i.index() )
                {
		    upper.push_back( (*j), jIndex );
		    ++colcount ;
		}
	    }
	    upper.rows_[ row+1 ] = colcount;
	}
	else {
	    upper.rows_[ row+1 ] = colcount;
	    inv[ row ] = diag_block;
	}
    }
}

template<class M, class CRS, class InvVector>
void updateConvertToCRS(const M& A, CRS& lower, CRS& upper, InvVector& inv, size_t interiorSize)
{
    // No need to do anything for 0 rows. Return to prevent indexing a
    // a zero sized array.
    if ( A.N() == 0 )
    {
	return;
    }

    typedef typename M :: size_type size_type;

    // implement left looking variant with stored inverse

    const auto endi = A.end();
    size_type nnz_count = 0;
    for (auto i=A.begin(); i.index() < interiorSize; ++i)
    {
	const size_type iIndex  = i.index();

	// eliminate entries left of diagonal; store L factor
	for (auto j=(*i).begin(); j.index() < iIndex; ++j )
	{
            lower.values_[nnz_count] = (*j);
            ++nnz_count;
	}
    }

    const auto rendi = A.beforeBegin();
    size_type row = 0;
    nnz_count = 0;

    // NOTE: upper and inv store entries in reverse order, reverse here
    // relative to ILU
    for (auto i=A.beforeEnd(); i!=rendi; --i, ++ row )
    {
	const size_type iIndex = i.index();

	if (iIndex < interiorSize) {
	    // store in reverse row order
	    // eliminate entries left of diagonal; store L factor
	    for (auto j=(*i).beforeEnd(); j.index()>=iIndex; --j )
	    {
		const size_type jIndex = j.index();
		if( j.index() == iIndex )
	        {
		    inv[ row ] = (*j);
		    break;
		}
		else if ( j.index() >= i.index() )
                {
		    upper.values_[nnz_count] = (*j);
		    ++nnz_count ;
		}
	    }
	}
    }
}
template<class M, class CRS, class InvVector>
void update_ghost_last_bilu0_decomposition (const M& A, CRS& lower, CRS& upper, InvVector& inv, size_t interiorSize)
{
    typedef typename M::RowIterator rowiterator;
    typedef typename M::ColIterator coliterator;
    typedef typename M::block_type block;

    for (rowiterator i = A.begin(); i.index() < interiorSize; ++i) {
	// coliterator is diagonal after the following loop
	coliterator endij=(*i).end();           // end of row i
	coliterator ij;

	auto id = i.index();
	auto invId = interiorSize - i.index(); 
	
    }
    
}

namespace Dune
{

namespace Amg
{

template<class M, class X, class Y, class C>
struct SmootherTraits<ParallelOverlappingILU0<M,X,Y,C> >
{
    using Arguments = ParallelOverlappingILU0Args<typename M::field_type>;
};

/// \brief Tells AMG how to construct the Opm::ParallelOverlappingILU0 smoother
/// \tparam Matrix The type of the Matrix.
/// \tparam Domain The type of the Vector representing the domain.
/// \tparam Range The type of the Vector representing the range.
/// \tparam ParallelInfo The type of the parallel information object
///         used, e.g. Dune::OwnerOverlapCommunication
template<class Matrix, class Domain, class Range, class ParallelInfo>
struct ConstructionTraits<ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfo> >
{
    typedef ParallelOverlappingILU0<Matrix,Domain,Range,ParallelInfo> T;
    typedef DefaultParallelConstructionArgs<T,ParallelInfo> Arguments;

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
    typedef std::shared_ptr< T > ParallelOverlappingILU0Pointer;
#else
    typedef T*                   ParallelOverlappingILU0Pointer;
#endif

    static inline ParallelOverlappingILU0Pointer construct(Arguments& args)
    {
        return ParallelOverlappingILU0Pointer(
                new T(args.getMatrix(),
                      args.getComm(),
                      args.getArgs().getN(),
                      args.getArgs().relaxationFactor) );
    }

#if ! DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
    // this method is not needed anymore in 2.7 since std::shared_ptr is used
    static inline void deconstruct(T* bp)
    {
        delete bp;
    }
#endif

};

} // end namespace Amg
} // end namespace Dune

template<class F>
class GhostLastILU0Args
    : public Dune::Amg::DefaultSmootherArgs<F>
{
 public:
    GhostLastILU0Args( )
        :  n_(0)
    {}
    
    void setN(int n)
    {
        n_ = n;
    }
    int getN() const
    {
        return n_;
    }
 private:
    int n_;
};

namespace Dune
{

namespace Amg
{

template<class M, class X, class Y, class C>
struct SmootherTraits<GhostLastILU0<M,X,Y,C> >
{
    using Arguments = GhostLastILU0Args<typename M::field_type>;
};

/// \brief Tells AMG how to construct the Opm::ParallelOverlappingILU0 smoother
/// \tparam Matrix The type of the Matrix.
/// \tparam Domain The type of the Vector representing the domain.
/// \tparam Range The type of the Vector representing the range.
/// \tparam ParallelInfo The type of the parallel information object
///         used, e.g. Dune::OwnerOverlapCommunication
template<class Matrix, class Domain, class Range, class ParallelInfo>
struct ConstructionTraits<GhostLastILU0<Matrix,Domain,Range,ParallelInfo> >
{
    typedef GhostLastILU0<Matrix,Domain,Range,ParallelInfo> T;
    typedef DefaultParallelConstructionArgs<T,ParallelInfo> Arguments;

#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
    typedef std::shared_ptr< T > GhostLastILU0Pointer;
#else
    typedef T*                   GhostLastILU0Pointer;
#endif

    static inline GhostLastILU0Pointer construct(Arguments& args)
    {
        return GhostLastILU0Pointer(
                new T(args.getMatrix(),
                      args.getComm(),
                      args.getArgs().getN(),
                      args.getArgs().relaxationFactor) );
    }

#if ! DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
    // this method is not needed anymore in 2.7 since std::shared_ptr is used
    static inline void deconstruct(T* bp)
    {
        delete bp;
    }
#endif

};

} // end namespace Amg
} // end namespace Dune

template <class PI>
class InteriorSizeSetter
{
public:
    size_t set_interiorSize(size_t N, size_t interiorSize, const PI& comm)
    {
	DUNE_UNUSED_PARAMETER(N);
	DUNE_UNUSED_PARAMETER(comm);
	return interiorSize;
    }
};

template<>
class InteriorSizeSetter<Dune::OwnerOverlapCopyCommunication<int,int>>
{
public:
    size_t set_interiorSize(size_t N, size_t interiorSize, const Dune::OwnerOverlapCopyCommunication<int,int>& comm)
    {
	if (interiorSize<N)
	    return interiorSize;
	auto indexSet = comm.indexSet();

	size_t new_is = 0;
	for (auto idx = indexSet.begin(); idx!=indexSet.end(); ++idx) {

	    if (idx->local().attribute()==1) {
		auto loc = idx->local().local();
		if (loc > new_is) {
		    new_is = loc;
		}
	    }
	}
	return new_is + 1;
    }
};

template<class Matrix, class Domain, class Range, class ParallelInfoT>
class ParallelOverlappingILU0
: public Dune::PreconditionerWithUpdate<Domain,Range>
{
    typedef ParallelInfoT ParallelInfo;


public:
    //! \brief The matrix type the preconditioner is for.
    typedef typename std::remove_const<Matrix>::type matrix_type;
    //! \brief The domain type of the preconditioner.
    typedef Domain domain_type;
    //! \brief The range type of the preconditioner.
    typedef Range range_type;
    //! \brief The field type of the preconditioner.
    typedef typename Domain::field_type field_type;

    typedef typename matrix_type::block_type  block_type;
    typedef typename matrix_type::size_type   size_type;

protected:
    struct CRS
    {
      CRS() : nRows_( 0 ) {}

      size_type rows() const { return nRows_; }

      size_type nonZeros() const
      {
        assert( rows_[ rows() ] != size_type(-1) );
        return rows_[ rows() ];
      }

      void resize( const size_type nRows )
      {
          if( nRows_ != nRows )
          {
            nRows_ = nRows ;
            rows_.resize( nRows_+1, size_type(-1) );
          }
      }

      void reserveAdditional( const size_type nonZeros )
      {
          const size_type needed = values_.size() + nonZeros ;
          if( values_.capacity() < needed )
          {
              const size_type estimate = needed * 1.1;
              values_.reserve( estimate );
              cols_.reserve( estimate );
          }
      }

      void push_back( const block_type& value, const size_type index )
      {
          values_.push_back( value );
          cols_.push_back( index );
      }

      void clear()
      {
          rows_.clear();
          values_.clear();
          cols_.clear();
          nRows_= 0;
      }

      std::vector< size_type  > rows_;
      std::vector< block_type > values_;
      std::vector< size_type  > cols_;
      size_type nRows_;
    };
public:
    Dune::SolverCategory::Category category() const override
    {
      return std::is_same<ParallelInfoT, Dune::Amg::SequentialInformation>::value ?
              Dune::SolverCategory::sequential : Dune::SolverCategory::overlapping;
    }

    /*! \brief Constructor.

      Constructor gets all parameters to operate the prec.
      \param A The matrix to operate on.
      \param n ILU fill in level (for testing). This does not work in parallel.
      \param w The relaxation factor.
      \param redblack Whether to use a red-black ordering.
      \param reorder_sphere If true, we start the reordering at a root node.
                            The vertices on each layer aound it (same distance) are
                            ordered consecutivly. If false, we preserver the order of
                            the vertices with the same color.
    */
    template<class BlockType, class Alloc>
    ParallelOverlappingILU0 (const Dune::BCRSMatrix<BlockType,Alloc>& A,
                             const int n, const field_type w)
        : lower_(),
          upper_(),
          inv_(),
          comm_(nullptr), w_(w),
          relaxation_( std::abs( w - 1.0 ) > 1e-15 ),
          A_(&reinterpret_cast<const Matrix&>(A)), iluIteration_(n)
    {
        interiorSize_ = A.N();
        // BlockMatrix is a Subclass of FieldMatrix that just adds
        // methods. Therefore this cast should be safe.
        update();
    }
    
    template<class BlockType, class Alloc>
    ParallelOverlappingILU0 (const Dune::BCRSMatrix<BlockType,Alloc>& A,
                             const ParallelInfo& comm,
                             const field_type w,
                             size_type interiorSize)
        : lower_(),
          upper_(),
          inv_(),
          comm_(&comm), w_(w),
          relaxation_( std::abs( w - 1.0 ) > 1e-15 ),
          interiorSize_(interiorSize),
          A_(&reinterpret_cast<const Matrix&>(A)), iluIteration_(0)
    {
        // BlockMatrix is a Subclass of FieldMatrix that just adds
        // methods. Therefore this cast should be safe.
        update( );
    }

    /*!
      \brief Prepare the preconditioner.

      \copydoc Preconditioner::pre(X&,Y&)
    */
    virtual void pre (Domain& x, Range& b) override
    {
        DUNE_UNUSED_PARAMETER(x);
        DUNE_UNUSED_PARAMETER(b);
    }

    /*!
      \brief Apply the preconditoner.

      \copydoc Preconditioner::apply(X&,const Y&)
    */
    virtual void apply (Domain& v, const Range& d) override
    {
	//Range& md = d;
	//Domain& mv = v;
        // iterator types
        typedef typename Range ::block_type  dblock;
        typedef typename Domain::block_type  vblock;

        const size_type iEnd = lower_.rows();
        const size_type lastRow = iEnd - 1;
        size_type upperLoppStart = iEnd - interiorSize_;
        size_type lowerLoopEnd = interiorSize_;
        if( iEnd != upper_.rows() )
        {
            OPM_THROW(std::logic_error,"ILU: number of lower and upper rows must be the same");
        }

        // lower triangular solve
        for( size_type i=0; i<lowerLoopEnd; ++ i )
        {
          dblock rhs( d[ i ] );
          const size_type rowI     = lower_.rows_[ i ];
          const size_type rowINext = lower_.rows_[ i+1 ];

          for( size_type col = rowI; col < rowINext; ++ col )
          {
            lower_.values_[ col ].mmv( v[ lower_.cols_[ col ] ], rhs );
          }

          v[ i ] = rhs;  // Lii = I
        }

        for( size_type i=upperLoppStart; i<iEnd; ++ i )
        {
            vblock& vBlock = v[ lastRow - i ];
            vblock rhs ( vBlock );
            const size_type rowI     = upper_.rows_[ i ];
            const size_type rowINext = upper_.rows_[ i+1 ];

            for( size_type col = rowI; col < rowINext; ++ col )
            {
                upper_.values_[ col ].mmv( v[ upper_.cols_[ col ] ], rhs );
            }

            // apply inverse and store result
            inv_[ i ].mv( rhs, vBlock);
        }

        copyOwnerToAll( v );

        if( relaxation_ ) {
            v *= w_;
        }
    }

    template <class V>
    void copyOwnerToAll( V& v ) const
    {
        if( comm_ ) {
            comm_->copyOwnerToAll(v, v);
	}
    }

    /*!
      \brief Clean up.

      \copydoc Preconditioner::post(X&)
    */
    virtual void post (Range& x) override
    {
        DUNE_UNUSED_PARAMETER(x);
    }

    virtual void update() override
    {
        
        int ilu_setup_successful = 1;
        std::string message;
        const int rank = ( comm_ ) ? comm_->communicator().rank() : 0;

        std::unique_ptr< Matrix > ILU;

	if( iluIteration_ == 0 ) {
	    /*
	      if (comm_) {
	      detail::InteriorSizeSetter<ParallelInfoT> iss;
	      interiorSize_ = iss.set_interiorSize(A_->N(), interiorSize_, *comm_);
	      }
	    */
	    // create ILU-0 decomposition
                
	    ILU.reset( new Matrix( *A_ ) );
	    Opm::detail::ghost_last_bilu0_decomposition(*ILU, interiorSize_);
	}
        // store ILU in simple CRS format
	Opm::detail::convertToCRS( *ILU, lower_, upper_, inv_ );
    }
protected:
    //! \brief The ILU0 decomposition of the matrix.
    CRS lower_;
    CRS upper_;
    std::vector< block_type > inv_;

    const ParallelInfo* comm_;
    //! \brief The relaxation factor to use.
    const field_type w_;
    const bool relaxation_;
    size_type interiorSize_;
    const Matrix* A_;
    int iluIteration_;
};


template<class Matrix, class Domain, class Range, class ParallelInfoT>
class GhostLastILU0
: public Dune::PreconditionerWithUpdate<Domain,Range>
{
    typedef ParallelInfoT ParallelInfo;


public:
    //! \brief The matrix type the preconditioner is for.
    typedef typename std::remove_const<Matrix>::type matrix_type;
    //! \brief The domain type of the preconditioner.
    typedef Domain domain_type;
    //! \brief The range type of the preconditioner.
    typedef Range range_type;
    //! \brief The field type of the preconditioner.
    typedef typename Domain::field_type field_type;

    typedef typename matrix_type::block_type  block_type;
    typedef typename matrix_type::size_type   size_type;

protected:
    struct CRS
    {
      CRS() : nRows_( 0 ) {}

      size_type rows() const { return nRows_; }

      size_type nonZeros() const
      {
        assert( rows_[ rows() ] != size_type(-1) );
        return rows_[ rows() ];
      }

      void resize( const size_type nRows )
      {
          if( nRows_ != nRows )
          {
            nRows_ = nRows ;
            rows_.resize( nRows_+1, size_type(-1) );
          }
      }

      void reserveAdditional( const size_type nonZeros )
      {
          const size_type needed = values_.size() + nonZeros ;
          if( values_.capacity() < needed )
          {
              const size_type estimate = needed * 1.1;
              values_.reserve( estimate );
              cols_.reserve( estimate );
          }
      }

      void push_back( const block_type& value, const size_type index )
      {
          values_.push_back( value );
          cols_.push_back( index );
      }

      void clear()
      {
          rows_.clear();
          values_.clear();
          cols_.clear();
          nRows_= 0;
      }

      std::vector< size_type  > rows_;
      std::vector< block_type > values_;
      std::vector< size_type  > cols_;
      size_type nRows_;
    };
public:
    Dune::SolverCategory::Category category() const override
    {
      return std::is_same<ParallelInfoT, Dune::Amg::SequentialInformation>::value ?
              Dune::SolverCategory::sequential : Dune::SolverCategory::overlapping;
    }

    /*! \brief Constructor.

      Constructor gets all parameters to operate the prec.
      \param A The matrix to operate on.
      \param n ILU fill in level (for testing). This does not work in parallel.
      \param w The relaxation factor.
      \param redblack Whether to use a red-black ordering.
      \param reorder_sphere If true, we start the reordering at a root node.
                            The vertices on each layer aound it (same distance) are
                            ordered consecutivly. If false, we preserver the order of
                            the vertices with the same color.
    */
    template<class BlockType, class Alloc>
    GhostLastILU0 (const Dune::BCRSMatrix<BlockType,Alloc>& A,
                             const int n, const field_type w)
        : lower_(),
          upper_(),
          inv_(),
          comm_(nullptr), w_(w),
          relaxation_( std::abs( w - 1.0 ) > 1e-15 ),
          A_(&reinterpret_cast<const Matrix&>(A)), iluIteration_(n)
    {
        interiorSize_ = A.N();
        // BlockMatrix is a Subclass of FieldMatrix that just adds
        // methods. Therefore this cast should be safe.
        update();
    }


    
    template<class BlockType, class Alloc>
    GhostLastILU0 (const Dune::BCRSMatrix<BlockType,Alloc>& A,
                             const ParallelInfo& comm,
                             const field_type w,
                             size_type interiorSize)
        : lower_(),
          upper_(),
          inv_(),
          comm_(&comm), w_(w),
          relaxation_( std::abs( w - 1.0 ) > 1e-15 ),
          interiorSize_(interiorSize),
          A_(&reinterpret_cast<const Matrix&>(A)), iluIteration_(0)
    {
        // BlockMatrix is a Subclass of FieldMatrix that just adds
        // methods. Therefore this cast should be safe.
        update( );
    }

    /*!
      \brief Prepare the preconditioner.

      \copydoc Preconditioner::pre(X&,Y&)
    */
    virtual void pre (Domain& x, Range& b) override
    {
        DUNE_UNUSED_PARAMETER(x);
        DUNE_UNUSED_PARAMETER(b);
    }

    /*!
      \brief Apply the preconditoner.

      \copydoc Preconditioner::apply(X&,const Y&)
    */
    virtual void apply (Domain& v, const Range& d) override
    {
	//Range& md = d;
	//Domain& mv = v;
        // iterator types
        typedef typename Range ::block_type  dblock;
        typedef typename Domain::block_type  vblock;

        const size_type iEnd = lower_.rows();
        const size_type lastRow = iEnd - 1;
        size_type upperLoppStart = iEnd - interiorSize_;
        size_type lowerLoopEnd = interiorSize_;
        if( iEnd != upper_.rows() )
        {
            OPM_THROW(std::logic_error,"ILU: number of lower and upper rows must be the same");
        }

        // lower triangular solve
        for( size_type i=0; i<lowerLoopEnd; ++ i )
        {
          dblock rhs( d[ i ] );
          const size_type rowI     = lower_.rows_[ i ];
          const size_type rowINext = lower_.rows_[ i+1 ];

          for( size_type col = rowI; col < rowINext; ++ col )
          {
            lower_.values_[ col ].mmv( v[ lower_.cols_[ col ] ], rhs );
          }

          v[ i ] = rhs;  // Lii = I
        }

        for( size_type i=upperLoppStart; i<iEnd; ++ i )
        {
            vblock& vBlock = v[ lastRow - i ];
            vblock rhs ( vBlock );
            const size_type rowI     = upper_.rows_[ i ];
            const size_type rowINext = upper_.rows_[ i+1 ];

            for( size_type col = rowI; col < rowINext; ++ col )
            {
                upper_.values_[ col ].mmv( v[ upper_.cols_[ col ] ], rhs );
            }

            // apply inverse and store result
            inv_[ i ].mv( rhs, vBlock);
        }

        copyOwnerToAll( v );

        if( relaxation_ ) {
            v *= w_;
        }
    }

    template <class V>
    void copyOwnerToAll( V& v ) const
    {
        if( comm_ ) {
            comm_->copyOwnerToAll(v, v);
	}
    }

    /*!
      \brief Clean up.

      \copydoc Preconditioner::post(X&)
    */
    virtual void post (Range& x) override
    {
        DUNE_UNUSED_PARAMETER(x);
    }

    virtual void update() override
    {
        
        int ilu_setup_successful = 1;
        std::string message;
        const int rank = ( comm_ ) ? comm_->communicator().rank() : 0;

        std::unique_ptr< Matrix > ILU;

	if( iluIteration_ == 0 ) {
	    
	      if (comm_) {
		  InteriorSizeSetter<ParallelInfoT> iss;
		  interiorSize_ = iss.set_interiorSize(A_->N(), interiorSize_, *comm_);
	      }
	    
	    // create ILU-0 decomposition
                
	    ILU.reset( new Matrix( *A_ ) );
	    Opm::detail::ghost_last_bilu0_decomposition(*ILU, interiorSize_);
	}
        // store ILU in simple CRS format
	convertToCRS( *ILU, lower_, upper_, inv_ , interiorSize_ );
    }

    void updateTest()
    {
        
        int ilu_setup_successful = 1;
        std::string message;
        const int rank = ( comm_ ) ? comm_->communicator().rank() : 0;

        std::unique_ptr< Matrix > ILU;

	if( iluIteration_ == 0 ) {
	    
	    // create ILU-0 decomposition
                
	    ILU.reset( new Matrix( *A_ ) );
	    Opm::detail::ghost_last_bilu0_decomposition(*ILU, interiorSize_);
	}
        // store ILU in simple CRS format
	updateConvertToCRS( *ILU, lower_, upper_, inv_, interiorSize_ );
    }

    void printLUOnRoot(int rank)
    {
	if (rank == 0) {

	    std::cout << "Lower" << std::endl;

	    auto start = lower_.rows_[0];
	    for (auto r = 1; r < lower_.rows_.size(); ++r) {

		std::cout << "row "<< r-1 << ": ";
		for (auto c = start; c < lower_.rows_[r]; ++c) {
		    std::cout << lower_.cols_[c] << " ";
		}
		start = lower_.rows_[r];;
		std::cout <<std::endl;
	    }
	    
	    std::cout <<std::endl;
	    std::cout << "Upper" << std::endl;
	    start = upper_.rows_[0];
	    for (auto r = 1; r < upper_.rows_.size(); ++r) {

		std::cout << "row "<< r-1 << " "<< upper_.rows_.size() - r - 1 <<": ";
		for (auto c = start; c < upper_.rows_[r]; ++c) {
		    std::cout << upper_.cols_[c] << " ";
		}
		start = upper_.rows_[r];;
		std::cout <<std::endl;
	    }

	    std::cout << lower_.rows_.size() << " " << upper_.rows_.size() <<std::endl;
	}
	
    }

    
protected:
    //! \brief The ILU0 decomposition of the matrix.
    CRS lower_;
    CRS upper_;
    std::vector< block_type > inv_;

    const ParallelInfo* comm_;
    //! \brief The relaxation factor to use.
    const field_type w_;
    const bool relaxation_;
    size_type interiorSize_;
    const Matrix* A_;
    int iluIteration_;
};

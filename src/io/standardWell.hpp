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

#ifndef OPM_STANDARDWELL_HEADER_INCLUDED
#define OPM_STANDARDWELL_HEADER_INCLUDED
#endif


class StandardWell
{
public:
    typedef Dune::FieldMatrix<double,4,3> BlockOff;
    typedef Dune::FieldMatrix<double,4,4> BlockD;
    typedef Dune::BCRSMatrix<BlockOff> MatOff;
    typedef Dune::BCRSMatrix<BlockD> MatD;

    typedef Dune::BlockVector<Dune::FieldVector<double,3>> Vec;
    typedef Dune::BlockVector<Dune::FieldVector<double,4>> WellVec;
    
    StandardWell(MatOff B, MatOff C, MatD D, std::vector<int> well_cells)
	: B_(B),C_(C),D_(D),well_cells_(well_cells)
    {
	Bx_.resize(D_.N());
	invD_.resize(D_.N());
    }

    void apply(const Vec& x, Vec& Ax) const
    {
	x_loc_.resize(well_cells_.size());
	Ax_loc_.resize(well_cells_.size());

	for (size_t i = 0; i < well_cells_.size(); ++i) {
	    x_loc_[i] = x[well_cells_[i]];
	    Ax_loc_[i] = Ax[well_cells_[i]];
	}
	
	B_.mv(x_loc_, Bx_);

	auto& invDBx = invD_;
	D_.mv(Bx_, invDBx);

	C_.mmtv(invDBx, Ax_loc_);
	
	for (size_t i = 0; i < well_cells_.size(); ++i) {
	    Ax[well_cells_[i]] =Ax_loc_[i];
	}
    }

    void applyscaleadd(const double alpha, const Vec& x, Vec& Ax) const
    {
	if( scaleAddRes_.size() != Ax.size() ) {
	    scaleAddRes_.resize( Ax.size() );
	}
	scaleAddRes_ = 0.0;
	apply( x, scaleAddRes_ );
	Ax.axpy( alpha, scaleAddRes_ );
    }
    
private:
    MatOff B_;
    MatOff C_;
    MatD D_;
    std::vector<int> well_cells_;

    mutable Vec x_loc_;
    mutable Vec Ax_loc_;
    mutable Vec scaleAddRes_;
    mutable WellVec Bx_;
    mutable WellVec invD_;
};

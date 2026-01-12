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
    typedef Dune::FieldMatrix<double,1,1> BlockP;
    typedef Dune::BCRSMatrix<BlockOff> MatOff;
    typedef Dune::BCRSMatrix<BlockD> MatD;
    typedef Dune::BCRSMatrix<Opm::MatrixBlock<double, 1, 1>> PressureMatrix;

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

    std::vector<int> wellConnections() const {return well_cells_;}

    void addWellPressureEquations(PressureMatrix& jacobian,
				  const Vec& weights,
				  const bool use_well_weights,
				  int well_index,
				  const int pressureVarIndex) const
    {
	int nperf = 0;
	auto cell_weights = weights[0];
	cell_weights = 0.0;
	const int number_cells = weights.size();
	const int welldof_ind = number_cells + well_index;
	
	bool bhp_control = false;

	int bhp_var_index = 3;
	if (!bhp_control || use_well_weights) {

	    for (auto colC = C_[0].begin(),
                  endC = C_[0].end(); colC != endC; ++colC) {

		const auto row_index = well_cells_[colC.index()];
		const auto& bw = weights[row_index];
		double matel = 0;
		assert((*colC).M() == bw.size());
		for (std::size_t i = 0; i < bw.size(); ++i) {
		    matel += (*colC)[bhp_var_index][i] * bw[i];
		}
		jacobian[row_index][welldof_ind] = matel;
		cell_weights += bw;
		nperf += 1;
	    }
	}
	cell_weights /= nperf;
	std::size_t blockSz = D_[0][0].size();
	WellVec bweights(1);
	
	bweights[0] = 0.0;
	double diagElem = 0;
	if (use_well_weights ) {
	    double abs_max = 0;
	    WellVec rhs(1);
	    rhs[0][bhp_var_index] = 1.0;
	    BlockD inv_diag_block = D_[0][0];
	    BlockD inv_diag_block_transpose = D_[0][0];

	    for (std::size_t i = 0; i < blockSz; ++i) {
		bweights[0][i] = 0;
		for (std::size_t j = 0; j < blockSz; ++j) {
		    bweights[0][i] += inv_diag_block_transpose[i][j] * rhs[0][j];
		}
		abs_max = std::max(abs_max, std::fabs(bweights[0][i]));
	    }
	    assert(abs_max > 0.0);
	    for (std::size_t i = 0; i < blockSz; ++i) {
		bweights[0][i] /= abs_max;
	    }
	    diagElem = 1.0 / abs_max;
	} else {
	    if ( bhp_control) {
		bweights[0][blockSz-1] = 1.0;
		diagElem = 1.0;
	    } else {
		for (std::size_t i = 0; i < cell_weights.size(); ++i) {
		    bweights[0][i] = cell_weights[i];
		}
		bweights[0][blockSz-1] = 0.0;
		diagElem = 0.0;
		const auto& locmat =     D_[0][0];
		for (std::size_t i = 0; i < cell_weights.size(); ++i) {
		    diagElem += locmat[i][bhp_var_index] * cell_weights[i];
		}
	    }
	}
	jacobian[welldof_ind][welldof_ind] = diagElem;
	if (!bhp_control || use_well_weights) {
	    for (auto colB = B_[0].begin(),
                  endB = B_[0].end(); colB != endB; ++colB) {

		const auto col_index = well_cells_[colB.index()];
		const auto& bw = bweights[0];
		double matel = 0;
		for (std::size_t i = 0; i < bw.size(); ++i) {
		    matel += (*colB)[i][pressureVarIndex] * bw[i];
		}
		jacobian[welldof_ind][col_index] = matel;
	    }
	}
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

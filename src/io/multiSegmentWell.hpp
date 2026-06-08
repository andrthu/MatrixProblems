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
/*
namespace Dune {
template<class M> class UMFPack;
}
#include <dune/istl/umfpack.hh>
*/
class MultiSegmentWell
{
public:
    typedef Dune::FieldMatrix<double,4,3> BlockOff;
    typedef Dune::FieldMatrix<double,4,4> BlockD;
    typedef Dune::BCRSMatrix<BlockOff> MatOff;
    typedef Dune::BCRSMatrix<BlockD> MatD;
    typedef Dune::BCRSMatrix<Opm::MatrixBlock<double, 1, 1>> PressureMatrix;
    
    typedef Dune::BlockVector<Dune::FieldVector<double,3>> Vec;
    typedef Dune::BlockVector<Dune::FieldVector<double,4>> WellVec;
    
    MultiSegmentWell(const MatOff& B, const MatOff& C, const MatD&  D, const std::vector<int>& well_cells)
	: B_(B),C_(C),D_(D),well_cells_(well_cells)
    {
	Bx_.resize(D_.N());
	invD_.resize(D_.N());
	duneDSolver_ = std::make_shared<Dune::UMFPack<MatD>>(D_, 0);
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

	WellVec bx(Bx_);
	WellVec invDBx(bx.size());
	
	invDBx=0.0;
	Dune::InverseOperatorResult res;
	duneDSolver_->apply(invDBx, bx, res);
	
	C_.mmtv(invDBx, Ax_loc_);
	
	for (size_t i = 0; i < well_cells_.size(); ++i) {
	    Ax[well_cells_[i]] = Ax_loc_[i];
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

    const std::vector<int>& wellConnections() const {return well_cells_;}

    void addWellPressureEquations(PressureMatrix& jacobian,
				  const Vec& weights,
				  const bool use_well_weights,
				  int well_index,
				  const int pressureVarIndex) const
    {
	const int number_cells = weights.size();
	const int welldof_ind = number_cells + well_index;
	
	bool bhp_control = false;
	const int seg_pressure_var_ind = 3;
	if (!bhp_control) {
	    for (std::size_t rowC = 0; rowC < C_.N(); ++rowC) {
		for (auto colC = C_[rowC].begin(),
			 endC = C_[rowC].end(); colC != endC; ++colC) {

		    const auto row_index = well_cells_[colC.index()];
		    const auto& bw = weights[row_index];
		    double matel = 0.0;

		    for (std::size_t i = 0; i< bw.size(); ++i) {
			matel += bw[i]*(*colC)[seg_pressure_var_ind][i];
		    }
		    jacobian[row_index][welldof_ind] += matel;
		}
	    }
	}
	if (!bhp_control) {
	    auto well_weight = weights[0];
	    well_weight = 0.0;
	    int num_perfs = 0;
	    for (std::size_t rowB = 0; rowB < B_.N(); ++rowB) {
		for (auto colB = B_[rowB].begin(),
			 endB = B_[rowB].end(); colB != endB; ++colB) {
		    const auto col_index = well_cells_[colB.index()];
		    const auto& bw = weights[col_index];
		    well_weight += bw;
		    num_perfs += 1;
		}
	    }
	    well_weight /= num_perfs;
	    assert(num_perfs > 0);
	

	    double diag_ell = 0.0;
	    for (std::size_t rowB = 0; rowB < B_.N(); ++rowB) {
		const auto& bw = well_weight;
		for (auto colB = B_[rowB].begin(),
			 endB = B_[rowB].end(); colB != endB; ++colB) {
		    const auto col_index = well_cells_[colB.index()];
		    double matel = 0.0;
		    for (std::size_t i = 0; i< bw.size(); ++i) {
			matel += bw[i] *(*colB)[i][pressureVarIndex];
		    }
		    jacobian[welldof_ind][col_index] += matel;
		    diag_ell -= matel;
		}
	    }
	    jacobian[welldof_ind][welldof_ind] = diag_ell;
	} else {
	    jacobian[welldof_ind][welldof_ind] = 1.0;
	}
    }

private:
    MatOff B_;
    MatOff C_;
    MatD D_;
    std::vector<int> well_cells_;

    mutable std::shared_ptr<Dune::UMFPack<MatD>> duneDSolver_;
    
    mutable Vec x_loc_;
    mutable Vec Ax_loc_;
    mutable Vec scaleAddRes_;
    mutable WellVec Bx_;
    mutable WellVec invD_;
};

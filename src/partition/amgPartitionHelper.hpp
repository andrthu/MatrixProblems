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
void createTransDiagonal(const MatrixType& mat1,
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
            // Note: If colIdx == rowIdx  
            // this will be overwritten by the rowSum later anyway.
            mat2[rowIdx][colIdx][0][0] = -val;
            
            // Accumulate absolute values for the diagonal
            rowSum += std::abs(val);
        }

        // 5. Assign the calculated sum to the diagonal of mat2
        mat2[rowIdx][rowIdx][0][0] = rowSum;
    }
}

template<class Crit>
void setCritForPart(Crit& criterion, Opm::PropertyTree prm_amg)
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
    //criterion.setRandomParallelGhostIndexOrder(prm_amg.get<bool>("random_coarse_ghost_index", true));
}

template<class SA>
void setILU0argsForPart(SA& smootherArgs, Opm::PropertyTree prm_amg)
{
    smootherArgs.iterations = prm_amg.get<int>("iterations", 1);
    const int iluwitdh = prm_amg.get<int>("iluwidth", 0);
    smootherArgs.setN(iluwitdh);
    const Opm::MILU_VARIANT milu = Opm::convertString2Milu(prm_amg.get<std::string>("milutype", std::string("ilu")));
    smootherArgs.setMilu(milu);
    smootherArgs.relaxationFactor = prm_amg.get<double>("relaxation", 1.0);
}

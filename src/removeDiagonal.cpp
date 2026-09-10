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

#include "config.h"
#include <iostream>
#include <fstream>
#include <string>

#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/fmatrix.hh>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmarket.hh>

// Assumes A is a non-block matrix, i.e. block size 1x1.
typedef Dune::FieldMatrix<double,1,1> BlockType;
typedef Dune::BCRSMatrix<BlockType> Matrix;

// Returns a copy of A with all diagonal entries removed from the
// sparsity pattern (not just zeroed out).
Matrix removeDiagonal(const Matrix& A)
{
    size_t nnz = 0;
    for (auto row = A.begin(); row != A.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            if (col.index() != row.index()) {
                ++nnz;
            }
        }
    }

    Matrix B;
    B.setBuildMode(Matrix::row_wise);
    B.setSize(A.N(), A.M(), nnz);

    for (auto row = B.createbegin(); row != B.createend(); ++row) {
        const auto i = row.index();
        for (auto col = A[i].begin(); col != A[i].end(); ++col) {
            if (col.index() != i) {
                row.insert(col.index());
            }
        }
    }

    for (auto row = A.begin(); row != A.end(); ++row) {
        const auto i = row.index();
        for (auto col = row->begin(); col != row->end(); ++col) {
            const auto j = col.index();
            if (j != i) {
                B[i][j] = *col;
            }
        }
    }

    return B;
}

int main(int argc, char** argv)
{
    try {
        Dune::MPIHelper::instance(argc, argv);

        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " <input.mtx> <output.mtx>" << std::endl;
            std::cerr << "Reads a (non-block) MatrixMarket matrix, removes its diagonal "
                       << "entries and writes the result to a new MatrixMarket file." << std::endl;
            return 1;
        }

        const std::string inputFile = argv[1];
        const std::string outputFile = argv[2];

        Matrix A;
        {
            std::ifstream file(inputFile);
            if (!file) {
                std::cerr << "Could not open input file " << inputFile << std::endl;
                return 1;
            }
            Dune::readMatrixMarket(A, file);
        }

        std::cout << "Read matrix with " << A.N() << " rows, " << A.M()
                   << " columns and " << A.nonzeroes() << " nonzeros from "
                   << inputFile << std::endl;

        const Matrix B = removeDiagonal(A);

        std::cout << "Removed diagonal: matrix now has " << B.nonzeroes()
                   << " nonzeros." << std::endl;

        {
            std::ofstream file(outputFile);
            Dune::writeMatrixMarket(B, file);
        }

        std::cout << "Wrote matrix to " << outputFile << std::endl;

        return 0;
    }
    catch (Dune::Exception& e) {
        std::cerr << "Dune reported error: " << e << std::endl;
        return 1;
    }
    catch (...) {
        std::cerr << "Unknown exception thrown!" << std::endl;
        return 1;
    }
}

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

// Identical to transMatrixCompare.cpp, except it assumes a 2x2-block
// system matrix (2-phase CO2-store systems, matching jsonSolveCO2.cpp)
// instead of a 3x3-block one (3-phase black-oil). See transMatrixCompare.cpp
// for the full explanation of what this tool does.

#include "config.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <limits>

#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/fmatrix.hh>

#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmarket.hh>

// Matrix block sizes matching jsonSolveCO2.cpp: a 2x2-block system matrix
// (2-phase CO2-store) and a scalar (1x1) transmissibility/connectivity
// matrix indexed by the same cell numbering.
typedef Dune::FieldMatrix<double,2,2> BlockMat2;
typedef Dune::BCRSMatrix<BlockMat2> Mat2;
typedef Dune::FieldMatrix<double,1,1> BlockMat1;
typedef Dune::BCRSMatrix<BlockMat1> Mat1;

void printUsage(const char* prog)
{
    std::cerr << "Compares transmissibility values against the corresponding off-diagonal\n"
	      << "block entries of a CO2-store system matrix and plots the result.\n\n"
	      << "Usage: " << prog << " matrix/dir/path [outputPrefix]\n\n"
	      << "matrix/dir/path must contain the same two files jsonSolveCO2 reads its\n"
	      << "system from:\n"
	      << "  *BlackoilMatrix.mtx: 2x2-block system matrix.\n"
	      << "  *transAdj.mtx: scalar transmissibility matrix for partitioning.\n\n"
	      << "If matrix/dir/path/porosity.csv is also present (columns\n"
	      << "cell,i,j,k,poro -- see the writePorosity tool), a second comparison\n"
	      << "of transmissibility/porosity vs. matrix entries is produced as well.\n\n"
	      << "outputPrefix (default \"trans_vs_matrix\") names the .csv and .png\n"
	      << "files written to the current directory." << std::endl;
}

// Reads a porosity.csv as written by writePorosity (header "cell,...,poro",
// one row per active cell, "cell" numbered 0..N-1 matching the matrix
// row/column numbering). Returns poro[cell] for cell in [0, N), or an empty
// vector if the file can't be read.
std::vector<double> readPorosity(const std::string& path)
{
    std::ifstream f(path);
    if (!f) {
        return {};
    }

    std::string header;
    std::getline(f, header);
    std::vector<std::string> cols;
    {
        std::stringstream ss(header);
        std::string col;
        while (std::getline(ss, col, ',')) {
            cols.push_back(col);
        }
    }
    int cellCol = -1, poroCol = -1;
    for (size_t k = 0; k < cols.size(); ++k) {
        if (cols[k] == "cell") cellCol = static_cast<int>(k);
        if (cols[k] == "poro") poroCol = static_cast<int>(k);
    }
    if (cellCol < 0 || poroCol < 0) {
        std::cerr << "Warning: " << path
                   << " does not have both a 'cell' and a 'poro' column, ignoring it"
                   << std::endl;
        return {};
    }

    std::vector<double> poro;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::vector<std::string> fields;
        std::stringstream ss(line);
        std::string field;
        while (std::getline(ss, field, ',')) {
            fields.push_back(field);
        }
        if (static_cast<int>(fields.size()) <= std::max(cellCol, poroCol)) {
            continue;
        }
        const int cell = std::stoi(fields[cellCol]);
        const double p = std::stod(fields[poroCol]);
        if (cell < 0) continue;
        if (static_cast<size_t>(cell) >= poro.size()) {
            poro.resize(cell + 1, std::numeric_limits<double>::quiet_NaN());
        }
        poro[cell] = p;
    }
    return poro;
}

// Runs one of the bundled plotting scripts on a CSV, labelled with xLabel
// on the x-axis. Returns the script's exit code (0 on success).
int runPlot(const std::string& scriptName, const std::string& csvName,
            const std::string& pngName, const std::string& xLabel)
{
    const std::string scriptPath = std::string(TRANSMATRIXCOMPARE_SCRIPT_DIR) + "/" + scriptName;
    const std::string cmd = "python3 \"" + scriptPath + "\" \"" + csvName + "\" \"" + pngName
        + "\" \"" + xLabel + "\"";

    std::cout << "Running: " << cmd << std::endl;
    const int ret = std::system(cmd.c_str());
    if (ret != 0) {
        std::cerr << "Plotting script failed (exit code " << ret << ") for " << pngName << std::endl;
    } else {
        std::cout << "Wrote plot to " << pngName << std::endl;
    }
    return ret;
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    if (argc < 2 || std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h") {
        printUsage(argv[0]);
        return (argc < 2) ? 1 : 0;
    }

    const std::string dir = argv[1];
    const std::string outPrefix = (argc > 2) ? argv[2] : "trans_vs_matrix";

    Mat2 A;
    Mat1 trans;

    try {
        const std::string aName = dir + "/BlackoilMatrix.mtx";
        const std::string tName = dir + "/transAdj.mtx";

        std::cout << "Reading " << aName << std::endl;
        {
            std::ifstream f(aName);
            if (!f) {
                std::cerr << "Could not open " << aName << std::endl;
                return 1;
            }
            Dune::readMatrixMarket(A, f);
        }

        std::cout << "Reading " << tName << std::endl;
        {
            std::ifstream f(tName);
            if (!f) {
                std::cerr << "Could not open " << tName << std::endl;
                return 1;
            }
            Dune::readMatrixMarket(trans, f);
        }
    }
    catch (Dune::Exception& e) {
        std::cerr << "Dune reported error while reading matrices: " << e << std::endl;
        return 1;
    }

    if (A.N() != trans.N() || A.M() != trans.M()) {
        std::cerr << "Warning: system matrix is " << A.N() << "x" << A.M()
                   << " blocks but transmissibility matrix is " << trans.N()
                   << "x" << trans.M() << " -- expected matching cell counts."
                   << std::endl;
    }

    const std::string poroName = dir + "/porosity.csv";
    const std::vector<double> poro = readPorosity(poroName);
    const bool havePoro = !poro.empty();
    if (havePoro) {
        std::cout << "Found " << poroName << ", also comparing transmissibility/porosity" << std::endl;
        if (poro.size() < A.N()) {
            std::cerr << "Warning: " << poroName << " only covers " << poro.size()
                       << " cells but the matrix has " << A.N()
                       << " rows -- rows beyond that will be skipped in the /porosity comparison"
                       << std::endl;
        }
    } else {
        std::cout << "No usable " << poroName
                   << " found, skipping the transmissibility/porosity comparison" << std::endl;
    }

    const std::string csvName = outPrefix + ".csv";
    const std::string poroCsvName = outPrefix + "_over_poro.csv";

    std::ofstream csv(csvName);
    if (!csv) {
        std::cerr << "Could not open " << csvName << " for writing" << std::endl;
        return 1;
    }
    csv.precision(std::numeric_limits<double>::digits10 + 1);
    csv << "i,j,trans,m00,m01,m10,m11\n";

    std::ofstream poroCsv;
    if (havePoro) {
        poroCsv.open(poroCsvName);
        if (!poroCsv) {
            std::cerr << "Could not open " << poroCsvName << " for writing" << std::endl;
            return 1;
        }
        poroCsv.precision(std::numeric_limits<double>::digits10 + 1);
        poroCsv << "i,j,trans,m00,m01,m10,m11\n";
    }

    size_t nPairs = 0;
    size_t nMissingTrans = 0;
    size_t nMissingPoro = 0;

    for (auto row = A.begin(); row != A.end(); ++row) {
        const auto i = row.index();
        for (auto col = row->begin(); col != row->end(); ++col) {
            const auto j = col.index();
            if (j == i) {
                continue; // diagonal blocks have no transmissibility counterpart
            }
            if (!trans.exists(i, j)) {
                ++nMissingTrans;
                continue;
            }
            const double t = trans[i][j];
            const BlockMat2& block = *col;

            csv << i << ',' << j << ',' << t;
            for (int r = 0; r < 2; ++r) {
                for (int c = 0; c < 2; ++c) {
                    csv << ',' << block[r][c];
                }
            }
            csv << '\n';
            ++nPairs;

            if (havePoro) {
                // Divide row i's transmissibility by row i's (i.e. cell i's)
                // porosity -- the equation this off-diagonal block belongs
                // to is equation i.
                const bool ok = i < poro.size() && !std::isnan(poro[i]) && poro[i] > 0;
                if (!ok) {
                    ++nMissingPoro;
                } else {
                    const double tOverPhi = t / poro[i];
                    poroCsv << i << ',' << j << ',' << tOverPhi;
                    for (int r = 0; r < 2; ++r) {
                        for (int c = 0; c < 2; ++c) {
                            poroCsv << ',' << block[r][c];
                        }
                    }
                    poroCsv << '\n';
                }
            }
        }
    }
    csv.close();
    if (havePoro) {
        poroCsv.close();
    }

    std::cout << "Wrote " << nPairs << " off-diagonal block pairs to " << csvName << std::endl;
    if (nMissingTrans > 0) {
        std::cout << "Skipped " << nMissingTrans
                   << " off-diagonal matrix blocks with no matching transmissibility entry"
                   << std::endl;
    }
    if (havePoro && nMissingPoro > 0) {
        std::cout << "Skipped " << nMissingPoro
                   << " pairs with no usable porosity for the row's cell (writing " << poroCsvName << ")"
                   << std::endl;
    }

    int failures = 0;
    failures += (runPlot("plot_trans_vs_matrix.py", csvName, outPrefix + ".png",
                          "transmissibility") != 0);
    failures += (runPlot("plot_trans_vs_matrix_abs.py", csvName, outPrefix + "_abs.png",
                          "transmissibility") != 0);
    if (havePoro) {
        failures += (runPlot("plot_trans_vs_matrix.py", poroCsvName, outPrefix + "_over_poro.png",
                              "transmissibility / porosity") != 0);
        failures += (runPlot("plot_trans_vs_matrix_abs.py", poroCsvName, outPrefix + "_over_poro_abs.png",
                              "transmissibility / porosity") != 0);
    }

    if (failures > 0) {
        std::cerr << failures << " plotting script(s) failed. "
                   << "The comparison data is still available in " << csvName
                   << (havePoro ? (" and " + poroCsvName) : std::string()) << std::endl;
        return 1;
    }

    return 0;
}

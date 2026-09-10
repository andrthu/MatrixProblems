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

// Reads the grid (EGRID) and static properties (INIT) of an Eclipse-format
// simulation case and writes the porosity of every active cell to a CSV
// file: one row per active cell, numbered 0..nactive-1, together with its
// (I,J,K) grid location (1-based, the usual Eclipse reporting convention)
// since active-cell numbering alone does not reveal the cell's position in
// the structured grid.
//
// Usage: writePorosity <rootName> <output.csv>
//   <rootName> is the case name without extension, e.g.
//   ../reservoir_out/example_norne_out/NORNE_ATW2013
//   (both NORNE_ATW2013.EGRID and NORNE_ATW2013.INIT must exist next to it)

#include "config.h"

#include <opm/io/eclipse/EGrid.hpp>
#include <opm/io/eclipse/EInit.hpp>

#include <array>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " <rootName> <output.csv>" << std::endl;
        std::cerr << "  reads <rootName>.EGRID and <rootName>.INIT" << std::endl;
        return EXIT_FAILURE;
    }

    const std::string rootName = argv[1];
    const std::string outName = argv[2];

    Opm::EclIO::EGrid grid(rootName + ".EGRID");
    Opm::EclIO::EInit init(rootName + ".INIT");

    const int nactive = grid.activeCells();

    const std::vector<float>& poro = init.getInitData<float>("PORO");
    if (static_cast<int>(poro.size()) != nactive) {
        std::cerr << "warning: PORO array size (" << poro.size()
                  << ") does not match number of active cells (" << nactive
                  << "); writing min(size, nactive) rows" << std::endl;
    }
    const int nrows = std::min(nactive, static_cast<int>(poro.size()));

    std::ofstream out(outName);
    if (!out) {
        std::cerr << "error: could not open " << outName << " for writing" << std::endl;
        return EXIT_FAILURE;
    }

    out << "cell,i,j,k,poro\n";
    for (int cell = 0; cell < nrows; ++cell) {
        const std::array<int, 3> ijk = grid.ijk_from_active_index(cell);
        out << cell << ','
            << (ijk[0] + 1) << ',' << (ijk[1] + 1) << ',' << (ijk[2] + 1) << ','
            << poro[cell] << '\n';
    }

    std::cout << "wrote " << nrows << " active cells to " << outName << std::endl;
    return EXIT_SUCCESS;
}

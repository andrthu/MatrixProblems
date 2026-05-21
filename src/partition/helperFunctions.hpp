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

template<class Mat>
void getWeightedEdgeCut(const Mat& W, const std::vector<int>& mpivec, int rank)
{

    if (rank == 0) {

	int edgeCut=0;
	double wgtEdgeCut=0;

	for (auto i = W.begin(); i != W.end(); ++i) {

	    auto rowIdx = i.index();
	    
	    for (auto j = i->begin(); j != i->end(); ++j) {
		auto colIdx = j.index();

		if (mpivec[rowIdx] != mpivec[colIdx]) {
		    edgeCut += 1;
		    wgtEdgeCut += (*j)[0][0];

		}
	    }
	}

	std::cout << std::endl;
	std::cout << "EdgeCutNormal: "<< edgeCut << std::endl;
	std::cout << "EdgeCutWeight: "<< wgtEdgeCut << std::endl;
	std::cout << std::endl;
    }
}

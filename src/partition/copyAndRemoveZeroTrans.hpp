/*
  Copyright 2024 Andreas Thune.

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


template<class Mat, class T>
void copyAndRemoveZeroTrans(const Mat& A, Mat& AN, const T& trans, const T& wells)
{
    Dune::MatrixIndexSet op;
    op.resize(A.N(),A.N());

    int Annz = 0;
    int nnnz = 0;
    for (auto row = A.begin(); row != A.end(); ++row) {
        auto d = row.index();
	op.add(d,d);

	for (auto col = row->begin(); col!= row->end(); ++col) {

	    auto n = col.index();
	    Annz++;
	    if (trans.exists(d,n)) {

	        if (trans[row.index()][col.index()] > 0) {
		    op.add(d,n);
		} else {
		    if (wells.exists(d,n)) {
		        op.add(d,n);
		    }
		}
	    } else {
	        op.add(d,n);
	    }
	}
    }

    op.exportIdx(AN);

    for (auto row = AN.begin(); row != AN.end(); ++row) {
        auto d = row.index();


	for (auto col = row->begin(); col!= row->end(); ++col) {

	    auto n = col.index();
	    nnnz++;
	    AN[d][n] = A[d][n];
	}
    }

    std::cout << "Old nnz: "<< Annz << " new nnz: " << nnnz << std::endl;
}

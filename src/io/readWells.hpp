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

#ifndef OPM_READWELLS_HEADER_INCLUDED
#define OPM_READWELLS_HEADER_INCLUDED

#endif // OPM_READWELLS_HEADER_INCLUDED

#include <boost/algorithm/string/predicate.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

class StoreRow
{
public:
    std::vector<int> cols;
    std::vector<double> vals;

    void insert(int c, double v) {
	cols.push_back(c);
	vals.push_back(v);
    }
    
};

template<class Mat>
void readWellMat(Mat& mat, std::string path)
{
    std::ifstream file(path);

    std::string line;
    // 1. Ignore the header/banner and comments (lines starting with '%')
    while (std::getline(file, line) && line[0] == '%') {
        if (line.substr(0, 15) == "%%MatrixMarket ") {
            // Optional: Parse the banner line to check format, storage, etc.
        }
    }

    std::stringstream ss(line);
    int M, N, L;
    if (!(ss >> M >> N >> L)) {
        std::cerr << "Error reading matrix dimensions." << std::endl;
        return;
    }

    //std::cout << M<< " " << N << " " << L << std::endl;

    typedef typename Mat::block_type Block; 
    int brow = Block::rows;
    int bcol = Block::cols;
    int bsize = brow*bcol;

    //std::cout << M/brow<< " " << N/bcol << " " << L/bsize << std::endl;

    int m = M/brow;
    int n = N/bcol;
    int nnz = L/bsize;

    mat.setSize(m,n,nnz);
    

    std::vector<StoreRow> rows(M);
    for (int v = 0; v < L; ++v) {

	std::getline(file, line);
	std::stringstream data_ss(line);
	int mrow, mcol;
	double val;
	data_ss >> mrow >> mcol >> val;

	//std::cout << mrow << " " << mcol << " " << val << std::endl;
	rows[mrow-1].insert(mcol-1, val);
    }

    mat.setBuildMode(Mat::row_wise);

    for(typename Mat::CreateIterator iter=mat.createbegin();  iter!= mat.createend(); ++iter) {

	for(std::size_t b=iter.index()*brow; b<iter.index()*brow+brow;++b) {
	    for (auto c = rows[b].cols.begin(); c!=rows[b].cols.end(); ++c) {
		//std::cout << iter.index() << " " << *c << " " << std::endl;
		iter.insert((*c)/bcol);
	    }
	}
    }
    mat=0;

    for (int r = 0; r < M; ++r) {

	for (int c=0; c < rows[r].cols.size(); ++c) {

	    int cc = rows[r].cols[c];
	    mat [r/brow][cc/bcol][r%brow][cc%bcol] = rows[r].vals[c];
	}
    }
}

void read_perf(std::vector<int>& perf, std::string path)
{
    std::ifstream file(path);
    int idx;
    while (file >> idx) {
	perf.push_back(idx);
    }
}
template<class WellMod, class MSW>
void readWellDir(std::string dirName, std::vector<WellMod>& wellMods, std::vector<MSW>& msWellMods)
{
    namespace fs = boost::filesystem;
    std::string wellDir = dirName + "/well";
    int file_count = 0;
    for (const auto& entry : fs::directory_iterator(wellDir)) {

        //std::cout << entry.path().filename().string() << std::endl;
	file_count++;
    }

    typedef Dune::FieldMatrix<double,4,3> BlockOff;
    typedef Dune::FieldMatrix<double,4,4> BlockD;

    typedef Dune::BCRSMatrix<BlockOff> MatOff;
    typedef Dune::BCRSMatrix<BlockD> MatD;

    if (file_count % 4 == 0) {

	int numWells = file_count/4;

	//numWells = 1;
	for (int id = 0; id < numWells; id++) {
	    
	    MatOff B;
	    MatOff C;
	    MatD D;
	    std::vector<int> well_cells;
	    bool standard = true;
	    
	    std::string base = wellDir + "/well" + std::to_string(id) + "_standard_";
	    if (! fs::exists(base+ "B.mtx")) {
		base = wellDir + "/well" + std::to_string(id) + "_multiSegment_";
		standard = false;
		//std::cout << "Well" << id << " is multi-segment" << std::endl;
	    }
	    std::string bname = base + "B.mtx";
	    std::string cname = base + "C.mtx";
	    std::string dname = base + "D.mtx";
	    std::string pname = base + "perf.txt";

	    read_perf(well_cells, pname);
	    readWellMat(B, bname);
	    readWellMat(C, cname);
	    readMatMarketObject(D,dname.data());

	    if (standard) {
		WellMod wm(B,C,D,well_cells);
		wellMods.push_back(wm);
	    } else {
		MSW msw(B,C,D,well_cells);
		msWellMods.push_back(msw);
	    }
	}
    }
}

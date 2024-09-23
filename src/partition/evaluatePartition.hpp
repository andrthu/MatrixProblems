/*
  Copyright 2018 Andreas Thune.

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

#ifndef OPM_EVALUATEPARTITION_HEADER_INCLUDED
#define OPM_EVALUATEPARTITION_HEADER_INCLUDED

#endif // OPM_EVALUATEPARTITION_HEADER_INCLUDED


template<class Mat, class Vec, class Comm>
void findCommunication(Mat& trans, Mat& wells, std::vector<int>& partition,
		       Vec& comTab, Vec& receiveTab, Comm cc)
{
    std::vector<int> visited;
    visited.resize(partition.size(), -1);
    
    int rank = cc.rank();
    
    double weightCut = 0.0;
    int cut = 0;
    for (auto row = trans.begin(); row != trans.end(); ++row)
    {
	if (partition[row.index()] == rank)
	{
	    auto col = row->begin();
	    for (; col != row->end(); ++col)
	    {
		int nabId = col.index();
		int nabPart = partition[nabId];
		if (nabPart != rank)
		{
		    comTab[nabPart]++;
		    cut++;
		    weightCut += *col;
		    if (visited[nabId]==-1)
		    {
			visited[nabId] = rank;
			receiveTab[nabPart]++;
		    }
		}
	    }
	}
    }

    double weightCut_ = cc.sum(weightCut);
    int cut_ = cc.sum(cut);
    
    if (rank == 0)
    {
	std::cout << "Edge-cut: "<< cut_/2 << std::endl;
	std::cout << "Weighted edge-cut: "<< weightCut_/2 << std::endl;
    }
}

template<class Vec>
void printComTab(Vec& comTab, int rank)
{
    std::cout << "ComTab rank " << rank << ": ";
    for (int i = 0; i < comTab.size(); ++i)
    {
	std::cout << comTab[i] << " ";  
    }
    std::cout << std::endl;
}

template<class C>
void printComTabOnRoot(C cc, std::vector<int>& comTab)
{
    //--- Send ---
    int size = cc.size();
    int rank = cc.rank();
    
    int tabSize = size*size;

    int* sendTab = comTab.data();
    int entireTab[tabSize];

    int start[size];
    int length[size];
    for (int i = 0; i < size; ++i) {
	start[i] = i*size;
	length[i] = size;
    }

    cc.gatherv(sendTab, size, entireTab, length, start, 0);

    //--- Print ---

    if (rank == 0) 
    {
	std::cout << std::endl;
	
	for (int prc = 0; prc < size; ++prc) 
	{
	    std::cout << "ComTab rank " << prc <<": ";
	    for ( int nab = 0; nab < size; ++nab)
	    {
		int idx = prc*size + nab;
		std::cout << entireTab[idx] << " ";
	    }
	    std::cout << std::endl;
	}
    }
}

template<class Vec>
int comNabCalc(Vec& comTab)
{
    int nabs = 0;
    for (int i = 0; i < comTab.size(); ++i)
    {
	if (comTab[i] > 0)
	{
	    nabs++;
	}
    }
    return nabs;
}

template<class Mat, class Comm>
void evaluatePartition(Mat& trans, Mat& wells, std::vector<int>& partition, Comm& cc)
{
    typedef Dune::BlockVector<Dune::FieldVector<double,1>> TabVec;
    
    int rank = cc.rank();
    int mpi_size = cc.size();
    
    std::vector<int> comTab(mpi_size, 0);
    std::vector<int> receiveTab(mpi_size, 0);
    //comTab = 0;
    //receiveTab = 0;
    
    findCommunication(trans, wells, partition, comTab, receiveTab, cc);
    int comNab = comNabCalc(comTab);
    
    //printComTab(comTab, rank);
    //printComTab(receiveTab, rank);
    printComTabOnRoot(cc, receiveTab);
    //std::cout << "ComNab rank " << rank << ": " << comNab << std::endl;
} 

template<class CC>
void printNumCells(const CC& cc, int outValue, int outType)
{
    int numCells = outValue;
    int size = cc.size();
    int rank = cc.rank();
    
    int cellVec[size];
    
    cc.gather(&numCells, cellVec, 1, 0);

    if (rank == 0)
    {
	std::cout << std::endl;
	for (int i = 0; i < size; ++i)
	{
	    int ic = cellVec[i];
	    if (outType == 0)
		std::cout << "Rank " << i << " numCells: "<< ic <<std::endl;
	    else if (outType == 1)
		std::cout << "Rank " << i << " InteriorSize: "<< ic <<std::endl;
	    else if (outType == 2)
		std::cout << "Rank "<< i << " nnz: " << ic << std::endl;
	    else if (outType == 3)
		std::cout << "Rank "<< i << " IntNnz: " << ic << std::endl;

	}
	std::cout << std::endl;
    }
}

template<class C, class M>
void evalWellCommOnRoot(const std::vector<int>& part, const M& wells, const C& cc)
{
    int numWells = -1;

    for (auto row = wells.begin(); row != wells.end(); ++row) {

	for (auto col = row->begin(); col != row->end(); ++col) {
	    int val = *col;
	    if (val >numWells)
		numWells = val;
	}
    }
    numWells++;
    if (cc.rank() == 0) {std::cout<< "Number of wells: "<< numWells<< std::endl;}

    std::vector<int> wellCells(cc.size(), 0);    
    std::vector<std::vector<int>> a2aCT(cc.size(), std::vector<int>(cc.size(), 0));
    if (cc.rank() == 0) {
	

	for (int wellIdx = 0; wellIdx < numWells; wellIdx++) {
	    std::vector<int> well_procs(cc.size(), 0);
	    for (auto row = wells.begin(); row != wells.end(); ++row) {
		//std::cout << row.index() << " " << part.size() << std::endl;
		int rank = part[row.index()];
		auto col = row->begin();
		if (col!= row->end()) {
		    if (*row->begin() == wellIdx) {
			//std::cout << rank << " " << *row->begin() << std::endl;
			wellCells[rank]++;
			well_procs[rank]++;
		    }
		}
	    }
	    std::cout << "Well " << wellIdx<<" present on rank: ";
	    for (int r = 0; r<cc.size(); ++r) {
		if (well_procs[r] > 0) {
		    std::cout << r <<" ";

		    for (int nab = 0; nab<cc.size();++nab){
			if (well_procs[nab] > 0) {
			    if (r!=nab)
				a2aCT[r][nab]++;
			}
		    }
		}
	    }
	    std::cout <<std::endl;
	}
    }
    if (cc.rank() == 0) {
	int sumwc = 0;
	std::cout <<std::endl;
	for (int r = 0; r < wellCells.size(); ++r ) {
	    sumwc += wellCells[r];
	    std::cout << "WellCells on " << r << " " << wellCells[r] << std::endl;
	}
	
	std::cout << "Total WellCells " << sumwc<< std::endl;
	std::cout <<std::endl;
	for (int prc = 0; prc < cc.size(); ++prc) {
	    std::cout << "a2aCT rank " << prc <<": ";
	    for ( int nab = 0; nab < cc.size(); ++nab) {
		std::cout <<a2aCT[prc][nab] << " ";
	    }
	    std::cout << std::endl;
	}
    }
}

template<class C, class M, class G>
void evalGhostOnRoot(const std::vector<int>& part, const M& A, const G& trans, const C& cc)
{
    
    if (cc.rank() == 0) {

	for (int rank = 0; rank < cc.size(); ++rank ) {

	    std::map<int,int> ghost;
	    std::map<int,int> ghost2;
	    int bs = 0;
	    for (auto row = A.begin(); row != A.end(); ++row) {

		if (rank == part[row.index()]) {

		    for (auto col = row->begin(); col!= row->end(); ++col) {

			if (rank != part[col.index()]) {

			    ghost[col.index()] = 1;
			    if (trans.exists(row.index(),col.index())) {
				if (trans[row.index()][col.index()] > 0) {
				    ghost2[col.index()] = 1;
				}
			    } else {
				bs++;
			    }
			}
		    }
		}
		
	    }
	    std::cout << "Pre-dist ghost Rank " << rank << ": " << ghost.size() << " remove fault: "<< ghost2.size()<< " mystery: "<< bs << std::endl; 
	    
	}
    }
}

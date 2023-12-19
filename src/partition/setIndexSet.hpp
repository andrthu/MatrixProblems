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

#ifndef OPM_SETINDEXSET_HEADER_INCLUDED
#define OPM_SETINDEXSET_HEADER_INCLUDED

#endif // OPM_SETINDEXSET_HEADER_INCLUDED

template<class PI>
void setParallelLocalIndex(PI& indexSet, std::vector<int>& local2global,
			   std::vector<int>& overlapMap)
{
    typedef Dune::OwnerOverlapCopyAttributeSet::AttributeSet AttributeSet;
    typedef Dune::ParallelLocalIndex<AttributeSet> PLI;
    indexSet.beginResize();
    
    for (int loc = 0; loc < local2global.size(); ++loc)
    {
	int glob = local2global[loc];
	
	if (overlapMap[glob]==0)
	{
	    indexSet.add(glob,PLI(loc,AttributeSet::owner));
	}
	else 
	{
	    indexSet.add(glob,PLI(loc,AttributeSet::copy));
        }
    }
    indexSet.endResize();
}
template<class P, class R>
Dune::IndexInfoFromGrid<int,int> getParallelInfo(P& pid, const R& rid)
{
    Dune::IndexInfoFromGrid<int,int> info;
    
    for (auto && i : pid)
    {	
	info.addLocalIndex(std::tuple<int,int,int>(i.global(),i.local().local(),i.local().attribute()));
    }
    
    for (auto && r : rid)
    {
	int rank = r.first;
	for(auto && rr : *r.second.first)
	{
	    auto glob = rr.localIndexPair().global();
	    auto ra = rr.attribute();
	    auto la = rr.localIndexPair().local().attribute();
	    info.addRemoteIndex(std::tuple<int,int,int>(rank,glob,ra));
	}
    }
    return info;
}

template<class P, class R>
Dune::IndexInfoFromGrid<int,int> getParallelInfoReorder(P& pid, const R& rid,
							std::vector<int> reorder)
{
    Dune::IndexInfoFromGrid<int,int> info;
    
    for (auto && i : pid)
    {	
	info.addLocalIndex(std::tuple<int,int,int>(i.global(),reorder[i.local().local()],i.local().attribute()));
    }
    
    for (auto && r : rid)
    {
	int rank = r.first;
	for(auto && rr : *r.second.first)
	{
	    auto glob = rr.localIndexPair().global();
	    auto ra = rr.attribute();
	    auto la = rr.localIndexPair().local().attribute();
	    info.addRemoteIndex(std::tuple<int,int,int>(rank,glob,ra));
	}
    }
    return info;
}


template<class Comm, class C>
void getIndexSetInfo(Comm comm, const C& cc, std::vector<int>& comTab, std::vector<int>& rowType, bool cellTogether=false, bool printInfo=true)
{
    auto indexSet = comm->indexSet();

    int rank = cc.rank();
    int size = cc.size();

    int ghost = 0;
    int total = 0;
    for (auto && i : indexSet) {

	int global = i.global();
	int local = i.local().local();
	int attr = i.local().attribute();
	total += 1;
	if (attr == 3) {
	    ghost += 1;
	    rowType[local] = 1;
	}
    }

    int ghostVec[size];
    int totalVec[size];

    cc.gather(&ghost, ghostVec, 1, 0);
    cc.gather(&total, totalVec, 1, 0);

    if (rank == 0) {
	if (printInfo) {
	    if (cellTogether) {
		std::cout << std::endl;
		for (int i = 0; i < size; ++i) {
		    int ic = totalVec[i];
		    int gc = ghostVec[i];
		    std::cout << "Rank " << i << " numCells,ghostCells,InteriorSize: "<< ic << " "<<gc<< " "<< ic-gc <<std::endl;
		}
		std::cout << std::endl;
	    }
	    else {
		std::cout << std::endl;
		for (int i = 0; i < size; ++i) {
		    int ic = totalVec[i];
		    std::cout << "Rank " << i << " numCells: "<< ic <<std::endl;
		}
		std::cout << std::endl;

		std::cout << std::endl;
		for (int i = 0; i < size; ++i) {
		    int ic = totalVec[i];
		    int gc = ghostVec[i];
		    std::cout << "Rank " << i << " ghostCells: "<< gc <<std::endl;
		}
		std::cout << std::endl;

		std::cout << std::endl;
		for (int i = 0; i < size; ++i) {
		    int ic = totalVec[i];
		    int gc = ghostVec[i];
		    std::cout << "Rank " << i << " InteriorSize: "<< ic-gc <<std::endl;
		}
		std::cout << std::endl;
	    }
	}
    }

    comTab.resize(size, 0);
    
    for (auto ri = comm->remoteIndices().begin(); ri!=comm->remoteIndices().end(); ++ri) {
	int numRemote = 0;
	int num3aa = 0;
	int num1aa = 0;

	for(auto && rr : *ri->second.first) {

	    auto glob = rr.localIndexPair().global();
	    //remote index attribute. If ra==1, then rank ri->first owns vertex glob.
	    auto ra = rr.attribute();
	    auto la = rr.localIndexPair().local().attribute();
	    if (ra == 1)
		num1aa++;
	    if (ra == 3)
		num3aa++;
	    numRemote++;
	}
	comTab[ri->first] = num1aa;
    }
}

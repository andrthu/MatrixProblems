/*
  Copyright 2019 Andreas Thune.

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

#ifndef OPM_REORDER_HEADER_INCLUDED
#define OPM_REORDER_HEADER_INCLUDED

#endif // OPM_REORDER_HEADER_INCLUDED

template<class Mat>
void CMK_reordering(Mat& A, std::vector<int>& reorder,
		    std::vector<int>& local2global,
		    std::vector<int>& global2local,
		    std::vector<int>& overlap)
{
    // Size of graph
    int V = local2global.size();

    std::vector<int> visited;
    visited.resize(V, 0);    
    reorder.resize(V, -1);
    
    int id = 0;
    for (int start = 0; start < V; ++start)
    {
	if (visited[start] == 0)
	{
	    std::list<int> queue;
	    queue.push_back(start);
	    visited[start] = 1;

	    while(!queue.empty())
	    {
		int v = queue.front();
		queue.pop_front();
		
		reorder[v] = id;
		id++;
		
		auto gv = A[local2global[v]];
		for (auto gnab = gv.begin(); gnab != gv.end(); ++gnab)
		{
		    int gidx = gnab.index();
		    int nab = global2local[gidx];
		    
		    if (visited[nab] == 0)
		    {
			if (overlap[gidx] > -1)
			{
			    queue.push_back(nab);
			    visited[nab] = 2;
			}
		    }
		}
	    }
	}
    }
}

template<class Mat>
void RCM_reordering(Mat& A, std::vector<int>& reorder,
		    std::vector<int>& local2global,
		    std::vector<int>& global2local,
		    std::vector<int>& overlap)
{
    CMK_reordering(A, reorder, local2global, global2local, overlap);
    
    int V = reorder.size();
    for (int v = 0; v < V; ++v)
    {
	int id = reorder[v];
	reorder[v] = V - id - 1;
    }
}

template<class Mat>
void trivial_reordering(Mat& A, std::vector<int>& reorder,
			std::vector<int>& local2global,
			std::vector<int>& global2local,
			std::vector<int>& overlap)
{
    int V = local2global.size();
    reorder.resize(V, 0);
    for (int i = 0; i < V; ++i)
    {
	reorder[i] = i;
    }
}

template<class Mat>
void overlap_reordering(Mat& A, std::vector<int>& reorder,
			std::vector<int>& local2global,
			std::vector<int>& global2local,
			std::vector<int>& overlap)
{
    int V = local2global.size();
    reorder.resize(V, 0);

    int idx = 0;
    for (int i = 0; i < V; ++i)
    {
	int gid = local2global[i];
	if (overlap[gid] == 0)
	{
	    reorder[i] = idx;
	    idx++;
	}
    }
    for (int i = 0; i < V; ++i)
    {
	int gid = local2global[i];
	if (overlap[gid] == 1)
	{
	    reorder[i] = idx;
	    idx++;	    
	}
    }
}
template<class Mat, class D>
void do_reorder(Mat& A, std::vector<int>& reorder,
		std::vector<int>& local2global,
		std::vector<int>& global2local,
		std::vector<int>& overlap, D dr)
{
    int method = std::stoi(dr.dict[8]);

    if (method == 0)
	trivial_reordering(A,reorder,local2global,global2local,overlap);
    else if (method == 1)
	CMK_reordering(A,reorder,local2global,global2local,overlap);
    else if (method == 2)
	RCM_reordering(A,reorder,local2global,global2local,overlap);
    else if (method == 3)
	overlap_reordering(A,reorder,local2global,global2local,overlap);
}

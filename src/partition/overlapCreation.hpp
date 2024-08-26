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

#ifndef OPM_OVERLAPCREATION_HEADER_INCLUDED
#define OPM_OVERLAPCREATION_HEADER_INCLUDED

#endif // OPM_OVERLAPCREATION_HEADER_INCLUDED

template<class Comm>
void myIds(std::vector<int>& local2global, std::vector<int>& global2local, 
	   std::vector<int>& overlap, Comm comm)
{
    global2local.resize(overlap.size(), -1);
    int lid = 0;
    for (int i = 0; i < overlap.size(); ++i)
    {
	if (overlap[i]>-1)
	{
	    local2global.push_back(i);
	    global2local[i] = lid;
	    lid++;
	}	
    }
}

void myRecIds(std::vector<int>& local2global, std::vector<int>& global2local, 
	      std::vector<int>& overlap, std::vector<int>& l2r, std::vector<int>& r2l)
{
    l2r.resize(local2global.size(), -1);
    int rid = 0;
    for (int i = 0; i < local2global.size(); ++i)
    {
	if (overlap[local2global[i]] == 0)
	{
	    r2l.push_back(i);
	    l2r[i] = rid;
	    rid++;
	}	
    }
}

template<class Mat, class Comm>
int addOverlap(std::vector<std::set<int>>& overlap, 
	       std::vector<int>& overlap2, 
	       const std::vector<int>& part, 
	       Mat& A, Comm comm)
{
    overlap2.resize(part.size(), -1);
    
    int rank = comm.rank();
    int partSize = 0;
    for (auto row=A.begin();row!=A.end();++row)
    {
	int id = row.index();
	if (part[id]==rank)
	{
	    partSize++;
	    overlap2[id] = 0;
	    
	    auto col = row->begin();
	    for (; col!=row->end(); ++col)
	    {
		int nab = col.index();
		if (part[nab]!=rank)
		{
		    overlap[nab].insert(rank);
		    overlap[id].insert(part[nab]);
		    
		    overlap2[nab] = 1;
		}
	    }
	}
    }
    return partSize;
}


template<class Mat, class Comm>
int addOverlapFromRoot(std::vector<std::set<int>>& overlap, 
		       std::vector<int>& overlap2, 
		       const std::vector<int>& part, 
		       Mat& A, Comm comm)
{
    int rank = comm.rank();
    
    overlap2.resize(part.size(), -1);
    if (rank == 0) {

	for (size_t proc = 0; proc < comm.size(); ++proc) {
	    std::vector<int> overlap2P;
	    overlap2P.resize(part.size(), -1);
	    int partSize = 0;
	    for (auto row=A.begin();row!=A.end();++row) {
		int id = row.index();
		if (part[id]==proc) {
		    partSize++;
		    overlap2P[id] = 0;
		    auto col = row->begin();
		    for (; col!=row->end(); ++col) {
			int nab = col.index();
			if (part[nab]!=proc) {
			    //overlap[nab].insert(rank);
			    //overlap[id].insert(part[nab]);
			    overlap2P[nab] = 1;
			}
		    }
		}
	    }
	}
    }
    
    
    
    int partSize = 0;
    for (auto row=A.begin();row!=A.end();++row)
    {
	int id = row.index();
	if (part[id]==rank)
	{
	    partSize++;
	    overlap2[id] = 0;
	    
	    auto col = row->begin();
	    for (; col!=row->end(); ++col)
	    {
		int nab = col.index();
		if (part[nab]!=rank)
		{
		    overlap[nab].insert(rank);
		    overlap[id].insert(part[nab]);
		    
		    overlap2[nab] = 1;
		}
	    }
	}
    }
    return partSize;
}

template<class Mat, class Comm>
void addNoOverlap(std::vector<std::set<int>>& overlap, 
		  std::vector<int>& overlap2, 
		  const std::vector<int>& part, 
		  Mat& A, Comm comm)
{
    overlap2.resize(part.size(), -1);
    
    int rank = comm.rank();
    for (auto row=A.begin();row!=A.end();++row)
    {
	int id = row.index();
	if (part[id]==rank)
	{
	    overlap2[id] = 0;	   	    
	}
    }
}

template<class Mat, class Vec, class Comm, class C>
void constructLocalFromRoot(Mat A, Mat& A_loc, int N, Vec& rhs, Vec& rhs_loc, const std::vector<int>& mpivec,
			    Comm& comm, std::shared_ptr<Comm>& parComm, const C& cc)
{
    int rank = cc.rank();
    typedef Dune::Amg::MatrixGraph<Mat> MatrixGraph;

    MatrixGraph graph = MatrixGraph(A);
    
    comm.buildGlobalLookup(graph.noVertices());
    Dune::fillIndexSetHoles(graph, comm);

    int nparts = cc.size();
    int myDomain = -1;
    std::vector<int> domainMapping(nparts);
    Dune::getDomain(cc, mpivec.data(), N, nparts, &myDomain, domainMapping);

    std::vector<int> setPartition(comm.indexSet().size(), -1);
    typedef typename  Comm::OwnerSet OwnerSet;
    typedef typename  Comm::ParallelIndexSet::const_iterator Iterator;
    std::size_t i=0;
    for(Iterator index = comm.indexSet().begin(); index != comm.indexSet().end(); ++index) {
	if(OwnerSet::contains(index->local().attribute())) {
	    setPartition[index->local()]=domainMapping[mpivec[i++]];
	}
    }
    comm.copyOwnerToAll(setPartition, setPartition);
    if (rank == 0) {std::cout << "Building coms etc. complete"<< std::endl;}
    cc.barrier();
    Dune::RedistributeInformation<Comm> redistInf;
    bool ret = Dune::buildCommunication(graph, setPartition, comm, parComm, redistInf.getInterface(), 0);
    redistInf.setSetup();
    cc.barrier();
    
    if (rank == 0) {std::cout << "Building coms2 etc. complete"<< std::endl;}
    if (cc.size() > 1)
	Dune::redistributeMatrix(A, A_loc, comm, *parComm, redistInf);
    else
	A_loc = A;
    //std::cout<< rank <<": myDoFs: "<< A_loc.N() <<std::endl;
    cc.barrier();
    if (rank == 0) {std::cout << "Redistribute Matrices complete"<< std::endl;}
    rhs_loc.resize(A_loc.N());
    redistInf.redistribute(rhs, rhs_loc);
}

template<class Mat, class Vec, class D, class Comm, class C>
std::vector<int> readMatOnRootAndDist(int argc, char** argv, Mat& A_loc, Vec& rhs_loc,
				      D& DR, Comm& comm, std::shared_ptr<Comm>& parComm,
				      const C& cc, bool gro=true)
{
    std::vector<int> dummy;
    return readMatOnRootAndDist(argc,argv,A_loc,rhs_loc,DR,comm,parComm,cc,dummy,gro,false);
}

template<class Mat, class Vec, class D, class Comm, class C>
std::vector<int> readMatOnRootAndDist(int argc, char** argv, Mat& A_loc, Vec& rhs_loc, D& DR, Comm& comm, std::shared_ptr<Comm>& parComm, const C& cc, std::vector<int>& partVec, bool gro, bool usePartVec)
{
    typedef Dune::FieldMatrix<double,1,1> Block;
    typedef Dune::BCRSMatrix<Block> Graph;

    int rank = cc.rank();

    Graph trans, wells;
    
    std::vector<int> mpivec;
    
    if (cc.size() > 0) {
	Vec rhs;
	Mat A, A_loc_;
    
	handleMatrixSystemInputSomeRanks(argc, argv, A, trans, wells, rhs, DR, cc, rank==0);
	DR.dict[5] = "2"; //Zoltan debug level
	if (rank == 0) {std::cout << "Reading Matrices complete"<< std::endl;}
	
	int N;
	if (rank == 0)
	    N = A.N();
	cc.broadcast(&N, 1, 0);

	std::vector<int> row_size(N);
	storeRowSizeFromRoot(A, row_size, cc);

	//partition matrix
	mpivec.resize(N, rank);
	if (cc.size() > 1) {
	    zoltanPartitionFunction(mpivec, trans, wells, cc, DR, row_size);
	    //evalWellCommOnRoot(mpivec,wells,cc);
	    //evalGhostOnRoot(mpivec,A,trans,cc);
	}
	if (rank == 0) {std::cout << "Zoltan partition complete"<< std::endl;}
	cc.barrier();
	constructLocalFromRoot(A, A_loc_, N, rhs, rhs_loc, mpivec, comm, parComm, cc);
	if (rank == 0) {std::cout << "Local Matrix construction complete"<< std::endl;}
	parComm->remoteIndices().template rebuild<false>();
	cc.barrier();
	std::vector<int> rowType(A_loc_.N(), 0);
	std::vector<int> comTab;
	getIndexSetInfo(parComm, cc, comTab, rowType);
	printComTabOnRoot(cc, comTab);

	//remove off-diagonal NNZ on ghost rows.
	if (cc.size() > 1)
	    if (gro)
		buildLocalMatrixFromLoc(A_loc_, A_loc, rowType);
	    else
		A_loc = A_loc_;
	else
	    A_loc = A_loc_;
	printNumCells(cc, A_loc.nonzeroes(), 2);
    }

    else {
	handleMatrixSystemInputSomeRanks(argc, argv, A_loc, trans, wells, rhs_loc, DR, cc, rank==0);
	//parComm = std::shared_ptr<Comm> (new Comm(cc, comm.category(), false));

	Dune::Amg::MatrixGraph<Mat> g = Dune::Amg::MatrixGraph<Mat>(A_loc);

	int nparts = cc.size();
	Dune::RedistributeInformation<Comm> redistInf;
	std::vector<int> setPartition(comm.indexSet().size(), -1);
	std::vector<int> mpivec(comm.indexSet().size(), 0);
	std::vector<int> domainMapping(nparts);
	int myDomain = -1;
	Dune::getDomain(cc, mpivec.data(), A_loc.N(), nparts, &myDomain, domainMapping);
	typedef typename  Comm::OwnerSet OwnerSet;
	typedef typename  Comm::ParallelIndexSet::const_iterator Iterator;
	std::size_t i=0;
	for(Iterator index = comm.indexSet().begin(); index != comm.indexSet().end(); ++index) {
	    if(OwnerSet::contains(index->local().attribute())) {
		setPartition[index->local()]=domainMapping[mpivec[i++]];
	    }
	}
	
	bool ret = Dune::buildCommunication(g, setPartition, comm, parComm, redistInf.getInterface(), 0);
	parComm->remoteIndices().template rebuild<false>();
    }
    return mpivec;
}

template<class Mat, class Vec, class D, class Comm, class C>
std::vector<int> readMatOnRootAndDist(std::string systemDir, Mat& A_loc, Vec& rhs_loc, D& DR, Comm& comm, std::shared_ptr<Comm>& parComm, const C& cc, bool gro=true)
{
    std::vector<int> dummy;
    return readMatOnRootAndDist(systemDir,A_loc,rhs_loc,DR,comm,parComm,cc,dummy,gro,false);
}

template<class Mat, class Vec, class D, class Comm, class C>
std::vector<int> readMatOnRootAndDist(std::string systemDir, Mat& A_loc, Vec& rhs_loc, D& DR, Comm& comm, std::shared_ptr<Comm>& parComm, const C& cc, std::vector<int> partVec, bool gro=true, bool usePartVec=false)
{
    typedef Dune::FieldMatrix<double,1,1> Block;
    typedef Dune::BCRSMatrix<Block> Graph;

    int rank = cc.rank();

    Graph trans, wells;

    std::vector<int> mpivec;
    if (cc.size() > 0) {
	Vec rhs;
	Mat A, A_loc_;

	if (rank == 0)
	    readFromDir(A,trans,wells,rhs,systemDir,rank, !usePartVec);

	DR.dict[5] = "2"; //Zoltan debug level
	if (rank == 0) {std::cout << "Reading Matrices complete"<< std::endl;}
	
	int N;
	if (rank == 0)
	    N = A.N();
	cc.broadcast(&N, 1, 0);

	std::vector<int> row_size(N);
	storeRowSizeFromRoot(A, row_size, cc);

	//partition matrix
	mpivec.resize(N, rank);
	if (cc.size() > 1) {
	    if (usePartVec)
		mpivec = partVec;
	    else
		zoltanPartitionFunction(mpivec, trans, wells, cc, DR, row_size);
	    //evalWellCommOnRoot(mpivec,wells,cc);
	    //evalGhostOnRoot(mpivec,A,trans,cc);
	}
	if (!usePartVec)
	    if (rank == 0) {std::cout << "Zoltan partition complete"<< std::endl;}
	cc.barrier();
	constructLocalFromRoot(A, A_loc_, N, rhs, rhs_loc, mpivec, comm, parComm, cc);
	if (rank == 0) {std::cout << "Local Matrix construction complete"<< std::endl;}
	parComm->remoteIndices().template rebuild<false>();
	cc.barrier();
	std::vector<int> rowType(A_loc_.N(), 0);
	std::vector<int> comTab;
	getIndexSetInfo(parComm, cc, comTab, rowType, false, !usePartVec);
	if (!usePartVec) {

	    printComTabOnRoot(cc, comTab);
	}
	//remove off-diagonal NNZ on ghost rows.
	if (cc.size() > 1)
	    if (gro)
		buildLocalMatrixFromLoc(A_loc_, A_loc, rowType);
	    else
		A_loc = A_loc_;
	else
	    A_loc = A_loc_;
	if (!usePartVec)
	    printNumCells(cc, A_loc.nonzeroes(), 2);
	else
	    if (rank == 0) {std::cout << std::endl;}
    }

    else {
	if (rank == 0)
	    readFromDir(A_loc,trans,wells,rhs_loc,systemDir,rank);
	//handleMatrixSystemInputSomeRanks(argc, argv, A_loc, trans, wells, rhs_loc, DR, cc, rank==0);
	//parComm = std::shared_ptr<Comm> (new Comm(cc, comm.category(), false));

	Dune::Amg::MatrixGraph<Mat> g = Dune::Amg::MatrixGraph<Mat>(A_loc);

	int nparts = cc.size();
	Dune::RedistributeInformation<Comm> redistInf;
	std::vector<int> setPartition(comm.indexSet().size(), -1);
	std::vector<int> mpivec(comm.indexSet().size(), 0);
	std::vector<int> domainMapping(nparts);
	int myDomain = -1;
	Dune::getDomain(cc, mpivec.data(), A_loc.N(), nparts, &myDomain, domainMapping);
	typedef typename  Comm::OwnerSet OwnerSet;
	typedef typename  Comm::ParallelIndexSet::const_iterator Iterator;
	std::size_t i=0;
	for(Iterator index = comm.indexSet().begin(); index != comm.indexSet().end(); ++index) {
	    if(OwnerSet::contains(index->local().attribute())) {
		setPartition[index->local()]=domainMapping[mpivec[i++]];
	    }
	}
	
	bool ret = Dune::buildCommunication(g, setPartition, comm, parComm, redistInf.getInterface(), 0);
	parComm->remoteIndices().template rebuild<false>();
    }
    return mpivec;
}

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

#include "amgPartitionHelper.hpp"

struct WgtIdx {

    double wgt;
    int idx;

    bool operator<(const WgtIdx& other) const {
	return wgt < other.wgt;
    }
};

template <class Mat>
class TransWellGraph
{
public:

    typedef Dune::BlockVector<Dune::FieldVector<double,1>> Vec;
    
    
    TransWellGraph(Mat T, Mat W, std::vector<int> rs, int useNormal, double logBase, bool isRoot)
    {
	if (isRoot) {
	    wgtType_ = useNormal;
	    trans = T;
	    wells = W;
	    row_size = rs;
	    base_ = std::exp(logBase);
	    if (useNormal < 2)
	    {
		edgeWgt_ = T;
		scaler_ = 1.0e18;
		sortTrans();
	    }
	    else if (useNormal == 2)
	    {
		findMin();
		createLogWeights();
		scaler_ = 1.0;
	    }
	    else if (useNormal == 4)
	    {
		sortTrans();
		createCatWeights();
		scaler_ = 1.0;
	    }
	}
    }

    void findMin()
    {
	double minVal = 1.0e18;
	for (auto trow = trans.begin(); trow != trans.end(); trow++)
	{	    		    	
	    auto tcol = trow->begin();
	    for (; tcol != trow->end(); ++tcol)
	    {
		if (*tcol != 0.0)
		{
		    if (minVal > *tcol)
		    {
			minVal = *tcol;
		    }
		}
	    }
	}
	minLogWgt = std::log(minVal);
	//std::cout << minLogWgt << " "<< minVal << std::endl; 
    }
    
    void createLogWeights()
    {
	Mat logWgt(trans);
	for (auto trow = trans.begin(); trow != trans.end(); trow++)
	{	    		    	
	    auto tcol = trow->begin();
	    int rid = trow.index();
	    for (; tcol != trow->end(); ++tcol)
	    {
		int cid = tcol.index();
		if (*tcol != 0.0)
		{
		    logWgt[rid][cid] = 1.0 + (std::log(*tcol) - minLogWgt)/std::log(base_);
		}
		else
		{
		    logWgt[rid][cid] = 0.0;
		}
		//if (logWgt[rid][cid] < 0)
		//std::cout << rid << " " << cid<< " "<< logWgt[rid][cid] << " " << *tcol<<" "<< trans[rid][cid]<<std::endl; 
	    }
	}	
	edgeWgt_ = logWgt;
    }

    void sortTrans()
    {
	unsigned nnz = trans.nonzeroes();
	std::vector<double>trans_list(nnz, 0.0);

	int idx = 0;
	for (auto row = trans.begin(); row != trans.end(); ++row)
	{
	    auto col = row->begin();
	    for (; col!=row->end(); ++col)
	    {
		trans_list[idx] = *col;	
	    }
	}
	
	std::sort(trans_list.begin(), trans_list.end());

	trans_bound_.push_back(trans_list[nnz/4]);
	trans_bound_.push_back(trans_list[3*(nnz/4)]);

    }

    double sortTransFindThreshold(double w)
    {
	std::vector<double>trans_list;

	for (auto row = trans.begin(); row != trans.end(); ++row) {
	    auto col = row->begin();
	    for (; col != row->end(); ++col) {
		trans_list.push_back(*col);
	    }
	}

	int nnz = trans_list.size();
	std::sort(trans_list.begin(), trans_list.end());

	//std::cout << "sort threshold "<< w << " " <<static_cast<int>( w*nnz)<< " " <<nnz << " " << trans_list[static_cast<int> (w*nnz)] <<std::endl;
	return trans_list[static_cast<int> (w*nnz)];
    }

    void createCatWeights()
    {
	Mat catWgt(trans);
	for (auto trow = trans.begin(); trow != trans.end(); trow++)
	{	    		    	
	    auto tcol = trow->begin();
	    int rid = trow.index();
	    for (; tcol != trow->end(); ++tcol)
	    {
		int cid = tcol.index();
		double t = *tcol;
		if (t == 0.0)
		{
		    catWgt[rid][cid] = 0.1;
		}
		else if (t < trans_bound_[0])
		{
		    catWgt[rid][cid] = 1.0;
		}
		else if ( t < trans_bound_[1] && t > trans_bound_[0])
		{
		    catWgt[rid][cid] = 10.0;
		}
		else
		{
		    catWgt[rid][cid] = 100.0;
		}
	    }
	}	
	edgeWgt_ = catWgt;
    }

    template<class R>
    void dps(R row, int v, int master, double w, std::vector<bool>& visited,
	     std::vector<int>& f2c, std::vector<int>& cnode, std::vector<std::tuple<int,int,double> >& edges) {

	visited[v] = true;
	f2c[v] = master;
	cnode.push_back(v);
	
	auto col = row.begin();
	for (; col != row.end(); ++col) {
	    int nab = col.index();

	    if (trans[v][nab] > w) {
		if (!visited[nab]) {
		    dps(trans[nab],nab,master,w,visited,f2c,cnode,edges);
		} else {
		    if (f2c[v]!=f2c[nab]) {
			std::cout << "Problem " << nab << " " << v <<
			    " " << f2c[v] << " " << f2c[nab] <<std::endl; 
		    }
		}
	    }
	}
	col = row.begin();
	for (; col != row.end(); ++col) {
	    int nab = col.index();
	    if (f2c[v]!=f2c[nab]) {
		edges.push_back({v,nab,trans[v][nab]});
	    }
	}
    }

    void createCoarseEdges(std::vector<std::vector<std::tuple<int,int,double> >> gEdges,
			   std::vector<int> f2c) {

	//std::vector<std::map<int, double> > cedges;

	for (std::vector<std::tuple<int,int,double> > es : gEdges ) {

	    std::map<int, double> ce;

	    for (std::tuple<int,int,double> fe : es) {

		int coarseNab = f2c[std::get<1>(fe)];
		double weight = wgtType_ == 0 ? 1.0 : std::get<2>(fe);
		if (std::get<2>(fe) > 0) {
		    if ( ce.count(coarseNab) == 1 ) {
			ce[coarseNab] += weight;
		    } else {
			ce.insert({coarseNab,weight});
		    }
		}
	    }
	    cedges.push_back(ce);
	}
    }
    
    void coarsenGraph(double w) {

	int N = trans.N();
	std::vector<bool> visited(N, false);

	f2c.resize(N, 0);
	std::vector<int> c2f;

	int biggest = 0;
	int single = 0;

	std::vector<std::vector<std::tuple<int,int,double> >> gEdges;
	int newV = 0;
	for (int v = 0; v < N; ++v) {

	    if (!visited[v]) {

		//f2c[v] = v;
		c2f.push_back(v);
		std::vector<int> cnode;
		std::vector<std::tuple<int,int,double> > edges;		
		dps(trans[v],v,newV,w,visited,f2c,cnode,edges);
		newV++;
		if (cnode.size() > biggest)
		    biggest = cnode.size();
		if (cnode.size() == 1)
		    single++;
		gEdges.push_back(edges);
		courseNodes_.push_back(cnode);
	    }
	}
	createCoarseEdges(gEdges, f2c);
	std::cout << "Coarse graph size: " << c2f.size()<< " "<< biggest<< " "<< single<< " " << N << std::endl;
    }

    template<class R, class Q>
    void dps(R row, Q q, int v, int master, double w, int mns, std::vector<bool>& visited,
	     std::vector<int>& f2c, std::vector<int>& cnode,
	     std::vector<std::tuple<int,int,double> >& edges) {

	visited[v] = true;
	f2c[v] = master;
	cnode.push_back(v);

	//int totSize = cnode.size() + q.size();
	auto col = row.begin();
	for (; col != row.end(); ++col) {
	    int nab = col.index();

	    if (trans[v][nab] > w) {
		if (!visited[nab]) {
		    q.push({trans[v][nab], nab});
		} 
	    }
	}


	if ( cnode.size() < mns ) {
	    if (!q.empty()) {
		auto strongCon = q.top();
		int nab = strongCon.idx;
		q.pop();
		while (visited[nab] && !q.empty()) {
		    strongCon = q.top();
		    nab = strongCon.idx;
		    q.pop();
		}
		if (!visited[nab])
		    dps(trans[nab],q,nab,master,w,mns,visited,f2c,cnode,edges);
		
	    }
	}
	
	col = row.begin();
	for (; col != row.end(); ++col) {
	    int nab = col.index();
	    if (f2c[v]!=f2c[nab]) {
		edges.push_back({v,nab,trans[v][nab]});
	    }
	}
    }

    void coarsenGraphMaxNodeSize(double w, int maxNodeSize, int rank) {

	int N = trans.N();
	std::vector<bool> visited(N, false);

	f2c.resize(N, 0);
	std::vector<int> c2f;

	int biggest = 0;
	int single = 0;
	int msns = 0;

	std::vector<std::vector<std::tuple<int,int,double> >> gEdges;
	int newV = 0;
	for (int v = 0; v < N; ++v) {

	    if (!visited[v]) {

		//std::cout << "start " << v<< std::endl;
		//f2c[v] = v;
		std::priority_queue<WgtIdx> q;
		c2f.push_back(v);
		std::vector<int> cnode;
		std::vector<std::tuple<int,int,double> > edges;		
		dps(trans[v],q,v,newV,w,maxNodeSize,visited,f2c,cnode,edges);
		newV++;
		if (cnode.size() > biggest)
		    biggest = cnode.size();
		if (cnode.size() == 1)
		    single++;
		if (cnode.size() == maxNodeSize)
		    msns++;
		gEdges.push_back(edges);
		courseNodes_.push_back(cnode);
		//std::cout << v <<" fin "<< cnode.size() << std::endl;
	    }
	}
	createCoarseEdges(gEdges, f2c);
	if (rank == 0)
	    std::cout << "Coarse graph size(csize,biggest,numSingle,fsize,numBig): " << c2f.size()<< " "<< biggest<< " "<< single<< " " << N << " "<<  msns << std::endl;
    }

    template<class Comm>
    void createAmgGraph(const Comm& comm, int level)
    {
	int rank = comm.communicator().rank();
	typedef Dune::MatrixAdapter<Mat, Vec, Vec> GLO;
	typedef Dune::SeqSSOR<Mat, Vec, Vec> ILU;
	typedef Dune::Amg::AMGCPR<GLO, Vec, ILU> AMGCPR;
	typedef typename AMGCPR::OperatorHierarchy Hir;
	using Smoother = ILU;
	using SmootherArgs = typename Dune::Amg::SmootherTraits<Smoother>::Arguments;
	using CriterionBase
	= Dune::Amg::AggregationCriterion<Dune::Amg::SymmetricDependency<Mat,Dune::Amg::FirstDiagonal>>;
	using Criterion = Dune::Amg::CoarsenCriterion<CriterionBase>;

	if (rank == 0) {
	    Mat diagT;
	    createTransDiagonal(trans, diagT);
	    //GLO glLinOp(diagT, comm);
	    GLO glLinOp(diagT);

	    Opm::FlowLinearSolverParameters flsp_amg;
	    Opm::PropertyTree prm_amg = setupAMG(std::string("amg"), flsp_amg);
	    prm_amg.put("skip_isolated",false);
	    Criterion criterion(15, prm_amg.get<int>("coarsenTarget", 1200));
	    setCritForPart(criterion, prm_amg);
	    SmootherArgs smootherArgs;
	    //setILU0argsForPart(smootherArgs, prm_amg);

	    //auto amg = std::make_shared<AMGCPR>(glLinOp, criterion, smootherArgs, comm);
	    auto amg = std::make_shared<AMGCPR>(glLinOp, criterion, smootherArgs);

	    
	    auto hir = amg->operatorHirarchyList();

	    auto aggMap = hir->aggregatesMaps();
	    auto mats = hir->matrices();
	    auto aggA = aggMap.begin();//opHir1->aggregatesMaps().begin();
	    int t = 0;
	    std::vector<std::vector<int>> maps;
	    std::vector<std::vector<int>> cum;
	    std::vector<int> aggSize;
	    auto cmat = mats.finest();

	    typedef typename Hir::AggregatesMap AMap;
    
	    for (auto aggA = aggMap.begin(); aggA!=aggMap.end(); aggA++) {
		std::vector<int> f2c_;

		//aggSize.push_back((*aggA)->end()-(*aggA)->begin());
		int numIso = 0;
		int y =0;
		//std::cout << (*aggA)->end()-(*aggA)->begin() << " " << cmat->getmat().N() <<std::endl;
		for (auto aa= (*aggA)->begin(); aa!=(*aggA)->end(); aa++,y++) {
		    //std::cout << *aa << " " << t<<std::endl;
		    if (*aa == AMap::ISOLATED) {
			//std::cout << *aa << " " << t<<std::endl;
			f2c_.push_back(-1);
			numIso++;
		    } else {
			f2c_.push_back(*aa);
		    }
		}
		maps.push_back(f2c_);
		if (t==0) {
		    cmat++;
		    //std::cout << f2c_.size() << " " << cmat->getmat().N() <<std::endl;
		    aggSize.push_back(cmat->getmat().N());
		    cum.push_back(f2c_);
		} else {
		    std::vector<int> pre = cum[t-1];
		    std::vector<int> nc(diagT.N(), 0);

		    if (f2c_.size() > 0) {
			cmat++;
			//std::cout << f2c_.size() << " " << cmat->getmat().N() <<std::endl;
			aggSize.push_back(cmat->getmat().N());
			for (int i = 0; i < diagT.N(); ++i) {
			    if (pre[i] == -1) 
				nc[i] = -1;
			    else
				nc[i] = f2c_[pre[i]];
			}
			cum.push_back(nc);
		    }
		}
		//std::cout << "Number of Isolated: "<<numIso << std::endl;
		t++;
	    }
	    std::cout << "End create map" << std::endl;

	    int f = level;

	    courseNodes_.resize(aggSize[f]);
	    cedges.resize(aggSize[f]);
	    f2c = cum[f];
	    for (auto i = diagT.begin(); i != diagT.end(); ++i) {
	    
		auto rowIdx = i.index();

		auto cidx = cum[f][rowIdx];
		if (cidx != -1) {
		    //std::cout <<rowIdx<<" "<<cidx<<std::endl;
		    courseNodes_[cidx].push_back(rowIdx);
	    
		    for (auto j = i->begin(); j != i->end(); ++j) {
			auto colIdx = j.index();
			auto cnab = cum[f][colIdx];

			if (std::abs(*j) > 0) {
			    if ( cedges[cidx].count(cnab) == 1 ) {
				cedges[cidx][cnab] += 1;
			    } else {
				cedges[cidx].insert({cnab,1});
			    }
			}
		    }
		}
	    }
	    int nsize = 0;
	    for (int x = 0; x < aggSize[f]; x++) {

		if (courseNodes_[x].size() > nsize)
		    nsize = courseNodes_[x].size();
	    }

	    std::cout << "Max node size: "<< nsize<< std::endl;
	}
    }

    Mat trans;
    Mat wells;
    Mat edgeWgt_;
    double minLogWgt;
    double scaler_;
    double base_;
    int wgtType_;

    
    std::vector<std::vector<int>> courseNodes_;
    std::vector<std::map<int, double> > cedges;
    std::vector<int> f2c;

    std::vector<double> trans_bound_;
    std::vector<int> row_size;
};

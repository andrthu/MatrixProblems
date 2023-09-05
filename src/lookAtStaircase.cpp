/*
  Copyright 2022 Andreas Thune.

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

#include "config.h"
#include <array>
#include <vector>
#include <cstdlib>
#include <memory>
#include <cmath>
#include <tuple>
#include <limits>
#include <chrono>
#include <unordered_set>
#include <fstream>
#include <string>
#include <map>
#include <numeric>
#include <algorithm>
#include <chrono>

#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/parallel/indexset.hh>
#include <dune/common/timer.hh>

#include <dune/common/function.hh>

#include <mpi.h>

#include "hardwarePaperFiles/nodeInfo.hpp"
#include "hardwarePaperFiles/cotaIO.hpp"

void cavium_map(int P, int N, bool bySocket, std::vector< std::vector<int> >& map) {

    int perNode = P/N;

    int perSoc = 32;

    for (int p = 0; p<P; ++p) {

	int ppn = p % perNode;
	int node = p/perNode;
	int core = ppn;
	int soc = ppn/perSoc;
	if (bySocket) {
	    soc = ppn%2;
	    core = ppn/2 + soc*perSoc;
	}
	std::vector<int> pMap = {core, soc, node};
	map.push_back(pMap);
    }
}

void fugaku_map(int P, int N, bool bySocket, std::vector< std::vector<int> >& map) {

    int perNode = P/N;

    int perSoc = 1;

    //std::cout << perNode << " " << P << " " << N << std::endl;
    for (int p = 0; p<P; ++p) {

	int ppn = p % perNode;
	int node = p/perNode;
	int core = ppn;
	int soc = ppn;
	if (bySocket) {
	    soc = ppn%2;
	    core = ppn/2 + soc*perSoc;
	}
	//std::cout << core << " " << soc << " " << node << std::endl;
	std::vector<int> pMap = {core, soc, node};
	map.push_back(pMap);
    }
    
}

void rome_map(int P, int N, bool bySocket, std::vector< std::vector<int> >& map, bool roundRobin) {

    int perNode = P/N;

    int perSoc = 64;

    for (int p = 0; p<P; ++p) {

	int ppn = p % perNode;
	int node = p/perNode;
	int core = ppn;
	int soc = ppn/perSoc;
	if (roundRobin) {

	    node = p % N;
	    core = p/N;
	    soc = core%2;
	}
	//std::cout << core << " "<< soc<< " " << node <<" "<<std::endl;
	std::vector<int> pMap = {core, soc, node};
	map.push_back(pMap);
    }
}

/*
    intraS = {0:0, 0.5:7.5, 1:14.6, 2:25.5, 3:26.8, 4:32.5, 5:37.5, 6:41.7, 7:43.3, 8:45, 9:45, 10:46.1,
              11:47.1, 12:47.1, 13:47.8, 14:49.9, 15:52.0, 16:54.0 }

    interS = {0:0, 1:13.0, 2:26.3, 3:30.95, 4:35.6, 5:35.6, 6:35.6, 7:36.1, 8:36.6, 9:36.7, 10:36.8, 11:37.1,
              12:37.4, 13:37.8, 14:38.2, 15:38.8, 16:39.4, 17:39.4, 18:39.4, 19:39.4, 20:39.4, 21:39.4, 22:39.4,
              23:39.4, 24:39.4, 25:39.4, 26:39.4, 27:39.4, 28:39.4, 29:42.4, 30:45.4, 31:45.4,
              32:45.4}

 */

void cavium_bw(std::vector< double>& on, std::vector< double>& off)
{
    on  = {0.00, 7.50, 14.6, 20.05, 25.5, 26.15, 26.8, 29.65, 32.5, 35.0, 37.5, 39.6, 41.7, 42.5, 43.3, 44.15, 45.0, 45.0, 45.0,
    	   45.55, 46.1, 46.6, 47.1, 47.1, 47.1, 47.45, 47.8, 48.85, 49.9, 50.95, 52.0, 53.0, 54.0 };

    off = {0.00, 13.0, 27.4, 31.5, 35.6, 35.6, 35.6, 36.1, 36.6, 36.7, 36.8, 37.1, 37.4, 37.8, 38.2, 38.8, 39.4, 39.4, 39.4,
	   39.4, 39.4, 39.4, 39.4, 39.4, 39.4, 39.4, 39.4, 39.4, 39.4, 42.4, 45.4, 45.4, 45.4};

    //on  = {0.00, 7.50, 14.6, 20.05, 25.5, 25.5, 26.8, 26.8, 32.5, 32.5, 37.6, 37.6, 41.7, 41.7, 43.3, 43.3, 45.0, 45.0, 45.0,
    //	   45.0, 46.1, 46.1, 47.1, 47.1, 47.1, 47.1, 47.8, 47.8, 49.9, 49.9, 52.0, 52.0, 54.0 };
    //off = {0.00, 10.5, 13.3, 13.3, 23.3, 23.3, 27.4, 27.4, 34.4, 34.4, 34.4, 34.4, 34.4, 34.4, 34.4, 34.4, 34.4, 34.4, 36.9,
    //36.9, 37.4, 37.4, 38.3, 38.3, 38.3, 38.3, 38.3, 40.1, 40.1, 40.1, 40.1, 45.0, 45.0};
}

void rome_bw(std::vector< double>& on, std::vector< double>& off)
{
    /*
    intraS = {0:0, 0.5: 10.2, 1:16.8, 2:17.6, 3:17.6, 4:19.2, 5:19.7, 
              6:19.7, 7:19.7, 8:23.4, 9:23.4, 10:23.4,
              11:23.4, 12:23.4, 13:23.4, 14:23.4, 15:23.4, 
              16:34.7, 17:34.7, 18:44.8, 19:44.8, 20:48.3,
              21:48.3, 22:48.3, 23:48.3, 24:48.3, 25:51.0, 
              26:51.0, 27:51.0, 28:51.0, 29:51.0, 30:51.0, 
              31:51, 32:51}

    interS = {0:0, 1:10.6, 2:17.4, 3:19.8, 4:22.2, 5:22.2, 6:22.2, 7:23.3, 8:24.4, 9:24.7, 10:25, 11:25, 12:25,
              13:25.5, 14:26.0, 15:26.0, 16:26.0, 17:31.8, 18:37.6, 19:37.6, 20:37.6, 21:39.3, 22:41.0, 23:41.0,
              24:41.0, 25:42.9, 26:44.8, 27:44.8, 28:44.8, 29:44.8, 30:44.8, 31:44.8, 32:44.8, 33:51.2, 34:57.6,
              35:57.6, 36:57.6, 37:57.6, 38:57.6, 39:58.0, 40:58.4, 41:58.4, 42:58.4, 43:58.4, 44:58.4, 45:58.4,
              46:58.4, 47:58.4, 48:58.4, 49:59.5, 50:60.6, 51:60.6, 52:60.6, 53:60.6, 54:60.6,55:60.6, 56:60.6,
              57:60.7, 58:60.6, 59:60.6, 60:60.6, 61:60.6, 62:60.6, 63:60.6, 64:60.6 }

     */
    on = {0, 10.2, 16.8, 17.2, 17.6, 17.6, 17.6, 18.4, 19.2, 19.45, 19.7,
	  19.7, 19.7, 19.7, 19.7, 21.55, 23.4, 23.4, 23.4, 23.4, 23.4,
	  23.4, 23.4, 23.4, 23.4, 23.4, 23.4, 23.4, 23.4, 23.4, 23.4,
	  29.05,34.7, 34.7, 34.7, 39.75, 44.8, 44.8, 44.8, 46.55, 48.3,
	  48.3, 48.3, 48.3, 48.3, 48.3, 48.3, 48.3, 48.3, 49.65, 51.0,
	  51.0, 51.0, 51.0, 51.0, 51.0, 51.0, 51.0, 51.0, 51.0, 51.0,
	  51.0, 51.0, 51.0, 51.0, 51.0};

    off = {0, 10.6, 17.4, 19.8, 22.2, 22.2, 22.2, 23.3, 24.4, 24.7, 25.0,
	   25.0, 25.0, 25.5, 26.0, 26.0, 26.0, 31.8, 37.6, 37.6, 37.6,
	   39.3, 41.0, 41.0, 41.0, 42.9, 44.8, 44.8, 44.8, 44.8, 44.8,
	   44.8, 44.8, 51.2, 57.6, 57.6, 57.6, 57.6, 57.6, 58.0, 58.4,
	   58.4, 58.4, 58.4, 58.4, 58.4, 58.4, 58.4, 58.4, 59.5, 60.6,
	   60.6, 60.6, 60.6, 60.6, 60.6, 60.6, 60.6, 60.6, 60.6, 60.6,
	   60.6, 60.6, 60.6, 60.6 };
}

void create_communication_type(std::vector< std::vector<int> >& cota, std::vector< std::vector<int> >& map,
			       std::vector<NodeInfo>& nodeInfo)
{
    auto P = cota.size();

    int eagerMes = 0;
    double eagerVol = 0;
    double totVol = 0;
    
    for (size_t p = 0; p<P; ++p) {
	int soc = map[p][1];
	int node = map[p][2];

	//std::cout << p << " " << soc << " " << node << " " << map.size() <<std::endl;
	int volOn = 0;
	int volOff = 0;
	int volInter = 0;

	int onMes = 0;
	int offMes = 0;
	int interMes = 0;
	
	std::vector< std::pair<int,int> > onNab;
	std::vector< std::pair<int,int> > offNab;
	std::vector< std::pair<int,int> > intraNab;
	std::vector< std::pair<int,int> > intraNabOn;
	std::vector< std::pair<int,int> > intraNabOff;
	std::vector< std::pair<int,int> > interNab;

	for (size_t nab = 0; nab < P; ++nab) {
	    
	    //int vol = cota[nab][p];
	    int vol = cota[p][nab];
	    int sendVol = vol;//cota[p][nab];
	    int nabNode = map[nab][2];
	    int nabSoc = map[nab][1];

	    // uncomment for comTab print 1/2
	    //std::cout << vol << " ";

	    if (vol < 8*1024) {
		eagerVol += vol;
		eagerMes += 1;
	    }
	    totVol += vol;
	    
	    if (vol > 0 ) {
		//std::cout <<soc<<" "<< p << " " << nab<< " " << vol<< std::endl;
		if (node == nabNode) {

		    if (soc == nabSoc) {
			onMes += 1;
			volOn += vol;
			onNab.push_back( std::make_pair(vol, nab) );
			intraNabOn.push_back( std::make_pair(sendVol, nab) );
		    } else {
			offMes+=1;
			volOff += vol;
			offNab.push_back( std::make_pair(vol, nab) );
			intraNabOff.push_back( std::make_pair(sendVol, nab) );
		    }
		    intraNab.push_back( std::make_pair(sendVol, nab) );
		} else {
		    interMes+=1;
		    volInter += vol;
		    interNab.push_back( std::make_pair(sendVol, nab) );
		}
	    }
	}
	// uncomment for comTab print 2/2
	//std::cout <<std::endl;
	nodeInfo[node].sockets_[soc].onVol_.push_back(volOn);
	nodeInfo[node].sockets_[soc].offVol_.push_back(volOff);
	if (volOn+volOff > 0)
	    nodeInfo[node].sockets_[soc].theta_.push_back( (double) volOn/((double)volOn+volOff));
	else
	    nodeInfo[node].sockets_[soc].theta_.push_back( 1.0 );
	nodeInfo[node].interVol_.push_back(volInter);

	nodeInfo[node].onMes_.push_back(onMes);
	nodeInfo[node].offMes_.push_back(offMes);
	nodeInfo[node].interMes_.push_back(interMes);
	
	nodeInfo[node].sockets_[soc].nabOn_.push_back(onNab);
	nodeInfo[node].sockets_[soc].nabOff_.push_back(offNab);
	nodeInfo[node].nabIntra_.push_back(intraNab);
	nodeInfo[node].nabIntraOn_.push_back(intraNabOn);
	nodeInfo[node].nabIntraOff_.push_back(intraNabOff);
	nodeInfo[node].nabInter_.push_back(interNab);
	nodeInfo[node].time_.push_back(0.0);
	nodeInfo[node].interTime_.push_back(0.0);
	nodeInfo[node].ranks_.push_back(p);
	nodeInfo[node].sockets_[soc].ranks_.push_back( nodeInfo[node].ranks_.size() - 1 );
    }

    //std::cout << "eager and totvol " <<eagerVol << " " << totVol << " " << eagerVol/totVol << std::endl;
}


void setEstTime(NodeInfo& Node, std::vector< std::vector<int> >& map, double onLat, double offLat)
{
    for (size_t p = 0; p < Node.time_.size(); ++ p) {
	size_t rank = Node.ranks_[p];
	int soc = map[p][1];
	int core = map[p][0];

	double recv = Node.time_[p];
	for ( const auto &nabS : Node.sendIntra_[p] ) {
	    int nab = nabS.first;
	    double rt = Node.sendIntra_[nab][p];
	    if (rt > recv)
		recv = rt;
	}
	double lat = Node.onMes_[p]*onLat + Node.offMes_[p]*offLat;
	Node.est_.push_back(recv + lat);
    }
}

void singleTypeEst(std::vector<NodeInfo>& nodeInfoVec, std::vector< std::vector<int> >& map)
{
    double BW = 10.0;
    double GB = 1024*1024*1024;
    std::vector<double> caviumOnBW;
    std::vector<double> caviumOffBW;
    cavium_bw(caviumOnBW, caviumOffBW);

    double onLat = 2.3e-6;
    double offLat = 4.4e-6;

    for (size_t N = 0; N < nodeInfoVec.size(); ++N) {
	//std::cout << "Node "<<N << " "<<nodeInfoVec.size()<<std::endl;
	auto Node = nodeInfoVec[N];
	for (int soc=0; soc < Node.numS_; ++soc) {

	    std::vector< int > onVol = Node.sockets_[soc].onVol_;
	    std::vector<size_t> idx(onVol.size());
	    std::iota(idx.begin(), idx.end(), 0);

	    std::sort(idx.begin(),idx.end(),
		      [&onVol](size_t i1, size_t i2) {return onVol[i1]<onVol[i2];});

	    int sentData = 0;
	    double last_t = 0.0;
	    for (size_t p = 0; p < idx.size(); ++ p) {

		size_t i = idx[p];
		size_t proc_i = i + soc * Node.sockets_[soc].P_;
		double t = (double)( idx.size() - p ) * ((double) (onVol[i]-sentData) ) / ( GB* caviumOnBW[idx.size() - p]);
		//std::cout << i <<" " << p <<" "<< t<< " "<< last_t + t << " " <<onVol[i] <<" "<< sentData << " "<< caviumOnBW[idx.size() - p]<<std::endl;
		Node.time_[proc_i] = last_t + t;
		last_t = last_t + t;
		sentData = onVol[i];

		//incoming message completion time
		std::vector< std::pair<int,int> > nabs = Node.nabIntra_[proc_i];
		std::sort(nabs.begin(), nabs.end(),
			  [](std::pair<int,int> i1, std::pair<int,int> i2 ) {return std::get<0>(i1) < std::get<0>(i2);} );

		std::map< int, double > sendTime;
		double st = 0.0;
		int prevMes = 0;
		for (size_t ni = 0; ni < nabs.size(); ++ni) {
		    int sij = std::get<0>(nabs[ni]);
		    st = st + t * (nabs.size()-ni)*(sij-prevMes)/(double)(sentData);
		    prevMes = sij;
		    sendTime[std::get<1>(nabs[ni])] = st;
		}
		nodeInfoVec[N].sendIntra_[proc_i] = sendTime;
	    }
	    //std::cout << "SentData: " << sentData << std::endl;
	}
	//setEstTime(Node, map, onLat, offLat);
	for (size_t p = 0; p < Node.time_.size(); ++ p) {
	    size_t rank = Node.ranks_[p];
	    int soc = map[p][1];
	    int core = map[p][0];

	    double recv = Node.time_[p];
	    for ( const auto &nabS : Node.sendIntra_[p] ) {
		int nab = nabS.first;
		double rt = Node.sendIntra_[nab][p];
		if (rt > recv)
		    recv = rt;
	    }
	    double lat = Node.onMes_[p]*onLat + Node.offMes_[p]*offLat;
	    //std::cout << recv <<" "<< lat << " "<< Node.time_[p] << std::endl;
	    nodeInfoVec[N].est_.push_back(recv + lat);
	}
    }
}

void interTypeEst(std::vector<NodeInfo>& nodeInfoVec, std::vector< std::vector<int> >& map, bool wait)
{

    double GB = 1024*1024*1024;
    
    //Rome vals
    //double BW = 11.5;
    //double rome_BW = 11.5;
    //double interLat = 1e-6;
    //double singleBW = 9.5;
    
    //Cavium vals
    //double BW = 10.5;
    //double cavium_BW = 10.5;
    //double singleBW = 9.1;
    //double interLat = 1.5e-6;

    //Fugaku vals rend
    //double BW = 6.2;
    //double rome_BW = 6.2;
    //double interLat = 2.1e-6;
    //double singleBW = 6;

    //Fugaku vals eager
    double BW = 0.52;
    double rome_BW = 0.52;
    double interLat = 0.9e-6;
    double singleBW = 0.5;
    
    auto begin = std::chrono::high_resolution_clock::now();
    for (size_t N = 0; N < nodeInfoVec.size(); ++N) {
	//std::cout << "Node "<<N << " "<<nodeInfoVec.size()<<std::endl;

	std::vector< int > interVol = nodeInfoVec[N].interVol_;
	std::vector<size_t> idx(interVol.size());
	std::iota(idx.begin(), idx.end(), 0);

	std::sort(idx.begin(),idx.end(),
		  [&interVol](size_t i1, size_t i2) {return interVol[i1]<interVol[i2];});

	BW = rome_BW;
	//BW = cavium_BW;
	int sentData = 0;
	double last_t = 0.0;
	for (size_t p = 0; p < idx.size(); ++ p) {

	    size_t i = idx[p];
	    //std::cout << "p: " << p << " " << i << std::endl;

	    size_t proc_i = i;
	    if (idx.size() - p == 1)
		BW = singleBW;
	    double t = (double)( idx.size() - p ) * ((double) (interVol[i]-sentData) ) / ( GB * BW );
	    //std::cout << i <<" "<< p <<" "<< t<< " "<< last_t + t << " " << interVol[i]-sentData  <<" "<< sentData << " "<< BW <<std::endl;
	    nodeInfoVec[N].interTime_[proc_i] = last_t + t;
	    last_t = last_t + t;
	    sentData = interVol[i];

	    //incoming message completion time
	    std::vector< std::pair<int,int> > nabs = nodeInfoVec[N].nabInter_[proc_i];
	    if (nabs.size() > 0) {
		std::sort(nabs.begin(), nabs.end(),
			  [](std::pair<int,int> i1, std::pair<int,int> i2 ) {return std::get<0>(i1) < std::get<0>(i2);} );
	    }
	    std::map< int, double > sendTime;
	    int trueSentData = 0;
	    for (size_t ni = 0; ni < nabs.size(); ++ni) {
		trueSentData += std::get<0>(nabs[ni]);
	    }
	    double st = 0.0;
	    int prevMes = 0;
	    for (size_t ni = 0; ni < nabs.size(); ++ni) {
		int sij = std::get<0>(nabs[ni]);
		st = st + last_t * (nabs.size()-ni)*(sij-prevMes)/((double) trueSentData);
		prevMes = sij;
		sendTime[std::get<1>(nabs[ni])] = st;
		//std::cout << last_t << " "<< st << std::endl;
	    }
	    //if (proc_i == 0) {
	    //		std::cout << "inter proc0 "<< last_t << " "<< st << std::endl;
	    //}
	    nodeInfoVec[N].sendInter_[proc_i] = sendTime;
	}
	//std::cout << "SentData: " << sentData << std::endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsedC = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf("Time measured before wait: %.6f seconds.\n", elapsedC.count() * 1e-9);
    
    for (size_t N = 0; N < nodeInfoVec.size(); ++N) {
	for (size_t p = 0; p < nodeInfoVec[N].interTime_.size(); ++ p) {
	    size_t rank = nodeInfoVec[N].ranks_[p];

	    double recv = nodeInfoVec[N].interTime_[p];
	    for ( const auto &nabS : nodeInfoVec[N].sendInter_[p] ) {

		int nabRank = nabS.first;
		int nabNodeId = map[nabRank][2];
		int nabCore = map[nabRank][0];

		double rt = nodeInfoVec[nabNodeId].sendInter_[nabCore][rank];
		//if (rank == 0) {
		//    std::cout <<"times for r0 inter: "<< recv << " " << rt<<std::endl; 
		//}
		if (rt > recv && wait) {
		    recv = rt;
		}
	    }
	    double lat = nodeInfoVec[N].interMes_[p]*interLat;
	    //std::cout << recv <<" "<< lat << " "<< Node.time_[p] << std::endl;
	    nodeInfoVec[N].interEst_.push_back(recv + lat);
	}
    }
    
}

double bwTet(double on, double off, double tet) {return tet*on + (1-tet)*off;}

void mixTypeEst(std::vector<NodeInfo>& nodeInfoVec, std::vector< std::vector<int> >& map, bool include_wait)
{
    std::vector<double> caviumOnBW(12.5,64);
    std::vector<double> caviumOffBW(10.5,64);

    double GB = 1024.0 * 1024.0 * 1024.0;
    std::vector<double> archOnBW;
    std::vector<double> archOffBW;
    rome_bw(archOnBW, archOffBW);
    //cavium_bw(archOnBW, archOffBW);
    
    //double onLat = 2.3e-6;
    //double offLat = 4.4e-6;

    double onLat = 1.7e-6;
    double offLat = 2.9e-6;

    for (size_t N = 0; N < nodeInfoVec.size(); ++N)
    {
	for (int soc=0; soc < nodeInfoVec[N].numS_; ++soc)
	{
	    // List of on-volume, off-volume and thetas on socket (soc) on Node (N)
	    std::vector< int > onVol = nodeInfoVec[N].sockets_[soc].onVol_;
	    std::vector< int > offVol = nodeInfoVec[N].sockets_[soc].offVol_;
	    std::vector< double > thetas = nodeInfoVec[N].sockets_[soc].theta_;

	    // Vector that stores core index of processs on socket (soc), range=0,...N_s-1 
	    std::vector<size_t> idx(onVol.size());
	    std::iota(idx.begin(), idx.end(), 0);

	    // Keep track of data recieved on each process
	    std::vector< double > recvTotal(offVol.size(), 0.0);

	    //running time estimate
	    double last_t = 0.0;
	    //std::cout << "Node socket step proc_q rank qIt" <<std::endl; 
	    // loop over processes on socket
	    for (size_t k = 0; k < onVol.size(); ++k)
	    {
		//set number of current active processes
		size_t active = onVol.size() - k;

		//on- and off-socket bandwidth given N=active active processes
		double bwOn  = archOnBW[active];
		double bwOff = archOffBW[active]/2;

		//Estimate complete time of first remaining process in idx 
		auto qIt = idx.begin();
		double firstTet = thetas[*qIt]; 
		double fbw = GB*bwTet(bwOn, bwOff, firstTet)/active;
		double tq = ( std::max( onVol[*qIt] + offVol[*qIt] - recvTotal[*qIt], 0.0 ) )/fbw;
		
		//std::cout << N<<" "<< soc<<" "<< k << " " << idx.size() << " "<< *qIt<< " "<< tq<<" "<< tet<< " "<<bwTet(bwOn, bwOff, tet) <<" "<< bwK<< std::endl;

		for (auto it = idx.begin() + 1 ; it!=idx.end(); ++it)
		{
		    double oTet = thetas[*it];
		    double obw = GB * bwTet(bwOn, bwOff, oTet)/active;
		    double tj = ( std::max( onVol[*it] + offVol[*it] - recvTotal[*it], 0.0 ) )/obw;
		    
		    if (tj<tq) {
			tq = tj;
			qIt = it;
			fbw = obw;
			firstTet = oTet;
		    }
		}

		//proc_q is process id (on node N) of fastest process
		int proc_q = *qIt + soc * nodeInfoVec[N].sockets_[soc].P_;
		int rank_fp = proc_q + 2*N*onVol.size();
		//if (N == 1 && soc == 1)
		//    std::cout << N<<" "<< soc<<" "<< k << " " << proc_q << " "<< rank_fp<<" "<< *qIt<< " "<< last_t + tq <<   " " << fbw/GB << " " << firstTet << std::endl;

		//Set recieve time estimate for proc_q
		last_t = last_t + tq;
		nodeInfoVec[N].time_[proc_q] = last_t;

		//total data recieved on proc_q is everything. 
		recvTotal[*qIt] = onVol[*qIt] + offVol[*qIt];
		int rem_q = *qIt;

		//Remove proc_q from active process list
		idx.erase(qIt,qIt+1);

		//Loop throgh all remaining active processes and update recieved data
		for (auto it = idx.begin(); it!=idx.end(); ++it) {
		    
		    double oTet = thetas[*it];
		    recvTotal[*it] += tq * GB*bwTet(bwOn, bwOff, oTet)/active;
		}

		//-- Start of incoming message completion time part

		//List of intra-node neighbours
		std::vector< std::pair<int,int> > nabsOn = nodeInfoVec[N].nabIntraOn_[proc_q];
		std::vector< std::pair<int,int> > nabsOff = nodeInfoVec[N].nabIntraOff_[proc_q];
		//Total amoun of data sent/recieved
		int totSentIntraOn = 0;
		for (size_t ni = 0; ni < nabsOn.size(); ++ni) {
		    totSentIntraOn += std::get<0>(nabsOn[ni]);
		}
		int totSentIntraOff = 0;
		for (size_t ni = 0; ni < nabsOff.size(); ++ni) {
		    totSentIntraOff += std::get<0>(nabsOff[ni]);
		}
		//Sort incoming messages
		std::sort(nabsOn.begin(), nabsOn.end(),
			  [](std::pair<int,int> i1, std::pair<int,int> i2 ) {return std::get<0>(i1) < std::get<0>(i2);} );
		std::sort(nabsOff.begin(), nabsOff.end(),
			  [](std::pair<int,int> i1, std::pair<int,int> i2 ) {return std::get<0>(i1) < std::get<0>(i2);} );

		//Track when each message is completed
		std::map< int, double > sendTime;
		double st = 0.0;
		int prevMes = 0;
		for (size_t ni = 0; ni < nabsOn.size(); ++ni) {
		    int sij = std::get<0>(nabsOn[ni]);
		    st = st + last_t * (nabsOn.size()-ni)*(sij-prevMes)/((double) totSentIntraOn );
		    prevMes = sij;
		    sendTime[std::get<1>(nabsOn[ni])] = st;
		}

		st = 0.0;
		prevMes = 0;
		for (size_t ni = 0; ni < nabsOff.size(); ++ni) {
		    int sij = std::get<0>(nabsOff[ni]);
		    st = st + last_t * (nabsOff.size()-ni)*(sij-prevMes)/((double) totSentIntraOff );
		    prevMes = sij;
		    sendTime[std::get<1>(nabsOff[ni])] = st;
		}
		nodeInfoVec[N].sendIntra_[proc_q] = sendTime;
	    }
	}
	//setEstTime(nodeInfoVec[N], map, onLat, offLat);
	for (size_t p = 0; p < nodeInfoVec[N].time_.size(); ++ p) {
	    size_t rank = nodeInfoVec[N].ranks_[p];
	    int soc = map[p][1];
	    int core = map[p][0];

	    double recv = nodeInfoVec[N].time_[p];
	    for ( const auto &nabS : nodeInfoVec[N].sendIntra_[p] ) {
		int nab = nabS.first;
		int nabCore = map[nab][0];
		double rt = nodeInfoVec[N].sendIntra_[nabCore][rank];
		//if (rank == 0 && N==0) {
		//    std::cout <<"times for r0: "<< recv << " " << rt<<" "<< nabCore <<" "<<   std::endl; 
		//}
		if (rt > recv && include_wait) {
		    recv = rt;
		}
	    }
	    double lat = nodeInfoVec[N].onMes_[p]*onLat + nodeInfoVec[N].offMes_[p]*offLat;
	    //std::cout << recv <<" "<< lat << " "<< Node.time_[p] << std::endl;
	    nodeInfoVec[N].est_.push_back(recv + lat);
	}
    }
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    std::vector< std::vector<int> > cota;
    std::vector<int> incoming;

    std::ifstream cotaFile;
    auto begin = std::chrono::high_resolution_clock::now();
    cotaFile.open(argv[1]);

    read_comTab2(cotaFile, cota, incoming);
    cotaFile.close();
    auto endF = std::chrono::high_resolution_clock::now();
    auto elapsedF = std::chrono::duration_cast<std::chrono::nanoseconds>(endF - begin);
    printf("Time measured for read: %.5f seconds.\n", elapsedF.count() * 1e-9);

    size_t ns = 4;
    size_t nsa = ns + 1;
    
    //int numNodes = cota.size()/nsa + 1;
    int numNodes = cota.size()/ns;
    std::vector< std::vector<int> > map;
    //cavium_map(cota.size(), cota.size()/nsa + 1 , false, map);
    //rome_map(cota.size(), cota.size()/nsa + 1 , false, map, true);
    //rome_map(cota.size(), cota.size()/nsa + 1 , false, map, false);
    //rome_map(cota.size(), 4 , false, map, false);
    fugaku_map(cota.size(), numNodes, false, map);
    std::vector<NodeInfo> nodeInfoVec;
    for (int n = 0; n < numNodes; ++n)
	//nodeInfoVec.push_back(NodeInfo(ns, 2));
	nodeInfoVec.push_back(NodeInfo(ns, 4));
    
    create_communication_type(cota, map, nodeInfoVec);

    //printComTabForRoundRobinMap(cota, 4);
    //printComTabForRoundRobinSocketMap(cota, 2);
    auto endMap = std::chrono::high_resolution_clock::now();
    auto elapsedM = std::chrono::duration_cast<std::chrono::nanoseconds>(endMap - endF);
    printf("Time measured map gen: %.6f seconds.\n", elapsedM.count() * 1e-9);

    //singleTypeEst(nodeInfoVec, map);
    mixTypeEst(nodeInfoVec, map, true);
    interTypeEst(nodeInfoVec, map, false);

    double fugaku_mean=0.0;
    double fugaku_max=0.0;
    
    for (int N = 0; N < nodeInfoVec.size(); ++N) {
	//for (auto& i : nodeInfoVec[N].interEst_) {
	for (size_t i = 0; i < nodeInfoVec[N].interEst_.size(); ++i) {
	    std::cout.precision(8);
	    std::cout << nodeInfoVec[N].interEst_[i]+ nodeInfoVec[N].est_[i] <<std::endl;
	    //std::cout << nodeInfoVec[N].est_[i] << std::endl;
	    //std::cout << nodeInfoVec[N].interEst_[i] <<std::endl;
	    fugaku_mean += (nodeInfoVec[N].interEst_[i] + nodeInfoVec[N].est_[i]);
	    if (nodeInfoVec[N].interEst_[i]+ nodeInfoVec[N].est_[i] > fugaku_max) {fugaku_max=nodeInfoVec[N].interEst_[i]+ nodeInfoVec[N].est_[i];}
	    //std::cout << nodeInfoVec[N].interEst_[i] + nodeInfoVec[N].est_[i] << " "<< nodeInfoVec[N].interEst_[i]<< std::endl;
	}
    }
    /*
    for (int N = 0; N < nodeInfoVec.size(); ++N) {

	for (int soc=0; soc < nodeInfoVec[N].numS_; ++soc) {

	    for (int i = 0; i < nodeInfoVec[N].sockets_[soc].onVol_.size();++i)
		std::cout << (nodeInfoVec[N].sockets_[soc].onVol_[i] + nodeInfoVec[N].sockets_[soc].offVol_[i])<< std::endl;
	}
    }

    std::cout << std::endl;

    for (int N = 0; N < nodeInfoVec.size(); ++N) {

	for (int i = 0; i < nodeInfoVec[N].interVol_.size();++i)
	    std::cout <<  nodeInfoVec[N].interVol_[i] << std::endl;
    }
    */
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsedC = std::chrono::duration_cast<std::chrono::nanoseconds>(end - endMap);
    //printf("Time measured calc: %.6f seconds.\n", elapsedC.count() * 1e-9);
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    //printf("Time measured: %.5f seconds.\n", elapsed.count() * 1e-9);

    //printf("mean: %.5f and max %.5f \n", fugaku_mean/512, fugaku_max);
    return 0;
}

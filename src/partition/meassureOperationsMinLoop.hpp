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

#ifndef OPM_MEASSUREOPERATIONSMINLOOP_HEADER_INCLUDED
#define OPM_MEASSUREOPERATIONSMINLOOP_HEADER_INCLUDED

#endif // OPM_MEASSUREOPERATIONSMINLOOP_HEADER_INCLUDED

template<class Op, class Vec, class Comm>
double timeMinLoopSpMV(Op& o, Vec& x, Comm& cc, int loopSize)
{
    Dune::Timer timer;
    Vec y(x.size());

#pragma omp parallel for
    for (size_t j = 0; j < x.size(); ++j) 
	y[j] = 0;

    double times[loopSize];
    

    for (int i = 0; i < loopSize; ++i) {
	cc.barrier();
	timer.reset();
	timer.start();
	o.apply(x, y);
	times[i] = timer.stop();
    }

    double tm = times[0];
    for (int i = 1; i < loopSize; ++i) 
	if (times[i]<tm) {tm = times[i];}

    return tm;
}

template<class Comm, class O, class Vec>
void multipleMinLoopTimeSpMV(Comm cc, O& linOp, Vec& x, int I=10, bool gl=false, int well=0)
{
    int rank = cc.rank();
    
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t1 = timeMinLoopSpMV(linOp, x, cc, 100);
	double times1[cc.size()];
	
	cc.gather(&t1, times1, 1, 0);

	if (rank==0)
	{
	    if (gl) {
		if (well == 0)
		    std::cout << "MeasGLapply+1 " << cc.size() << ": ";
		else if (well == 1)
		    std::cout << "MeasGLcompAdjApply+1 " << cc.size() << ": ";
		else if (well == 2)
		    std::cout << "MeasGLsepBCDApply+1 " << cc.size() << ": ";
	    }
	    else
		std::cout << "MeasSpMV+1 " << cc.size() << ": ";
	    for (int i = 0; i < cc.size(); i++)
		std::cout << times1[i] << " ";
	    std::cout << std::endl;
	}
    }
}


template<class Op, class Vec, class Comm>
double timeMinLoopSpMVAS(Op& o, Vec& x, Comm& cc, int loopSize)
{
    Dune::Timer timer;
    Vec y(x.size());

#pragma omp parallel for
    for (size_t j = 0; j < x.size(); ++j) 
	y[j] = 0;

    double times[loopSize];
    

    for (int i = 0; i < loopSize; ++i) {
	cc.barrier();
	timer.reset();
	timer.start();
	o.applyscaleadd(-1, x, y);
	times[i] = timer.stop();
    }

    double tm = times[0];
    for (int i = 1; i < loopSize; ++i) 
	if (times[i]<tm) {tm = times[i];}

    return tm;
}

template<class Comm, class O, class Vec>
void multipleMinLoopTimeSpMVAS(Comm cc, O& linOp, Vec& x, int I=10, bool gl=false, int well=0)
{
    int rank = cc.rank();
    
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t1 = timeMinLoopSpMVAS(linOp, x, cc, 100);
	double times1[cc.size()];
	
	cc.gather(&t1, times1, 1, 0);

	if (rank==0)
	{
	    if (gl)
		std::cout << "MeasSpMV-GL+1 " << cc.size() << ": ";
	    else
		std::cout << "MeasSpMV+1 " << cc.size() << ": ";
	    for (int i = 0; i < cc.size(); i++)
		std::cout << times1[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template<class Op, class Vec, class Comm>
double timeMinLoopTLP(Op& o, Vec& x, Comm& cc, int loopSize)
{
    Dune::Timer timer;
    
    double times[loopSize];
    
    for (int i = 0; i < loopSize; ++i) {
	cc.barrier();
	timer.reset();
	timer.start();
	o.moveToCoarseLevel(x);
	times[i] = timer.stop();
    }

    double tm = times[0];
    for (int i = 1; i < loopSize; ++i) 
	if (times[i]<tm) {tm = times[i];}

    return tm;
}

template<class Comm, class O, class Vec>
void multipleMinLoopTimeTLP(Comm cc, O& linOp, Vec& x, int I=10, bool gl=false, int well=0)
{
    int rank = cc.rank();
    
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t1 = timeMinLoopTLP(linOp, x, cc, 100);
	double times1[cc.size()];
	
	cc.gather(&t1, times1, 1, 0);

	if (rank==0)
	{
	    std::cout << "MeasTLP+1 " << cc.size() << ": ";
	    for (int i = 0; i < cc.size(); i++)
		std::cout << times1[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template<class Op, class Vec, class Comm>
double timeMinLoopTLPC(Op& o, Vec& x, Comm& cc, int loopSize)
{
    Dune::Timer timer;
    
    double times[loopSize];
    
    for (int i = 0; i < loopSize; ++i) {
	cc.barrier();
	timer.reset();
	timer.start();
	o.moveToFineLevel(x);
	times[i] = timer.stop();
    }

    double tm = times[0];
    for (int i = 1; i < loopSize; ++i) 
	if (times[i]<tm) {tm = times[i];}

    return tm;
}

template<class Comm, class O, class Vec>
void multipleMinLoopTimeTLPC(Comm cc, O& linOp, Vec& x, int I=10, bool gl=false, int well=0)
{
    int rank = cc.rank();
    
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t1 = timeMinLoopTLPC(linOp, x, cc, 100);
	double times1[cc.size()];
	
	cc.gather(&t1, times1, 1, 0);

	if (rank==0)
	{
	    std::cout << "MeasTLPC+1 " << cc.size() << ": ";
	    for (int i = 0; i < cc.size(); i++)
		std::cout << times1[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template<class Vec, class Comm, class C>
double timeMinLoopComm(Vec& x, const C& comm, Comm& cc, int loopSize)
{
    Dune::Timer timer;
    Vec y=x;

    double times[loopSize];
      
    for (int i = 0; i < loopSize; ++i) {
	cc.barrier();
	timer.reset();
	timer.start();
	comm.copyOwnerToAll(y, y);
	times[i] = timer.stop();

	y.axpy(0.01, x);
    }

    double tm = times[0];
    for (int i = 1; i < loopSize; ++i) { 
	//std::cout << times[i] << " "<< cc.rank()<<std::endl;
	if (times[i]<tm) {
	    tm = times[i];
	    
	}
    }

    return tm;
}

template<class Comm, class C, class Vec>
void multipleMinLoopTimeComm(Comm cc, const C& comm, Vec& x, int I=10)
{
    int rank = cc.rank();
    
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t1 = timeMinLoopComm(x, comm, cc, 100);
	double times1[cc.size()];
	
	cc.gather(&t1, times1, 1, 0);

	if (rank==0)
	{
	    std::cout << "MeasComm+1 " << cc.size() << ": ";
	    for (int i = 0; i < cc.size(); i++)
		std::cout << times1[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template<class Vec, class Comm, class C>
double timeMinLoopProject(Vec& x, const C& comm, Comm& cc, int loopSize)
{
    Dune::Timer timer;
    Vec y=x;

    double times[loopSize];
      
    for (int i = 0; i < loopSize; ++i) {
	cc.barrier();
	timer.reset();
	timer.start();
	comm.project(y);
	times[i] = timer.stop();

	y.axpy(0.01, x);
    }

    double tm = times[0];
    for (int i = 1; i < loopSize; ++i) { 
	//std::cout << times[i] << " "<< cc.rank()<<std::endl;
	if (times[i]<tm) {
	    tm = times[i];
	    
	}
    }

    return tm;
}

template<class Comm, class C, class Vec>
void multipleMinLoopTimeProject(Comm cc, const C& comm, Vec& x, int I=10)
{
    int rank = cc.rank();
    
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t1 = timeMinLoopProject(x, comm, cc, 100);
	double times1[cc.size()];
	
	cc.gather(&t1, times1, 1, 0);

	if (rank==0)
	{
	    std::cout << "MeasProject+1 " << cc.size() << ": ";
	    for (int i = 0; i < cc.size(); i++)
		std::cout << times1[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template <class S, class Vec, class Comm>
double timeMinLoopIP(S& sp, Vec& x, Comm cc, int loopSize)
{
    Vec y(x.size());
    y.axpy(0.000178,x);
    //y *= 0.178;

    double val = 0;
    double val2 = 0;
    Dune::Timer timer;

    double times[loopSize];
      
    for (int i = 0; i < loopSize; ++i) {
	cc.barrier();
	timer.reset();
	timer.start();
	val = sp.dot(x, y);
	times[i] = timer.stop();
	y.axpy(0.001*val,x);
	//std::cout << val << " " << times[i] <<std::endl;
	val2+=val;
    }
    
    double tm = times[0];
    for (int i = 1; i < loopSize; ++i) { 
	//std::cout << times[i] << " "<< cc.rank()<<std::endl;
	if (times[i]<tm) {
	    tm = times[i];
	    
	}
    }

    return tm;
}

template<class Comm, class O, class Vec>
void multipleMinLoopTimeSP(Comm cc, O sp, Vec& x, int I=10, bool gl=false)
{
    int rank = cc.rank();
    
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t1 = timeMinLoopIP(sp, x, cc, 100);
	double times1[cc.size()];
	
	cc.gather(&t1, times1, 1, 0);

	if (rank==0)
	{
	    if (gl)
		std::cout << "MeasGLsp+1 " << cc.size() << ": ";
	    else
		std::cout << "MeasSP+1 " << cc.size() << ": ";
	    for (int i = 0; i < cc.size(); i++)
		std::cout << times1[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template <class S, class Vec, class Comm>
double timeMinLoopSolver(S& solver, Vec& rhs, Comm cc, int loopSize)
{
    Vec x(rhs.size());
    Vec rhs_ = rhs;
    x=0;
    Dune::InverseOperatorResult stat;
    
    Dune::Timer timer;
    double times[loopSize];
    
    for (int i = 0; i < loopSize; ++i) {
	cc.barrier();
	timer.reset();
	timer.start();
	solver.apply(x, rhs_, stat);
	times[i] = timer.stop();
	x=0;
	rhs_=rhs;
	//std::cout << val << " " << times[i] <<std::endl;
    }
    
    double tm = times[0];
    for (int i = 1; i < loopSize; ++i) { 
	//std::cout << times[i] << " "<< cc.rank()<<std::endl;
	if (times[i]<tm) {
	    tm = times[i];
	    
	}
    }

    return tm;
}

template<class Comm, class O, class Vec>
void multipleMinLoopTimeSolver(Comm cc, O solver, Vec& rhs, int I=10, bool gl=false, bool hundred=false)
{
    int rank = cc.rank();
    
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t1 = timeMinLoopSolver(solver, rhs, cc, 10);
	double times1[cc.size()];
	
	cc.gather(&t1, times1, 1, 0);

	if (rank==0)
	{
	    if (hundred) {
		if (gl)
		    std::cout << "MeasGLhundred+1 " << cc.size() << ": ";
		else
		    std::cout << "MeasHundred+1 " << cc.size() << ": ";
		
	    }
	    else {
		if (gl)
		    std::cout << "MeasGLoneIter+1 " << cc.size() << ": ";
		else
		    std::cout << "MeasOneIter+1 " << cc.size() << ": ";
		
	    }
	    for (int i = 0; i < cc.size(); i++)
		std::cout << times1[i] << " ";
	    std::cout << std::endl;
	}
    }
}


template<class Solver, class Vec>
double timeCoarseSolverOnRoot(Solver solver, Vec& x)
{
    Dune::Timer timer;
    Vec y (x.size());
    y=0;
    int loopSize = 10;
    double times[loopSize];
    for (int i = 0; i < loopSize; ++i) {
	Dune::InverseOperatorResult res;
	timer.reset();
	timer.start();
	solver.apply(y, x, res);
	times[i] = timer.stop();
    }
    double tm = times[0];
    for (int i = 1; i < loopSize; ++i) {
	if (times[i]<tm) {
	    tm = times[i];
	}
    }
    return tm;
}

template <class Pre, class Vec, class Comm>
double timeMinLoopPre(Pre& pre, Vec& x, Comm cc, int loopSize)
{
    Dune::Timer timer;
    Vec y(x.size());
    y = 0;
    double times[loopSize];
    for (int i = 0; i < loopSize; ++i) {
	cc.barrier();
	timer.reset();
	timer.start();
	pre.apply(y, x);
	times[i] = timer.stop();
	
    }

    double tm = times[0];
    for (int i = 1; i < loopSize; ++i) { 
	if (times[i]<tm) {
	    tm = times[i];	    
	}
    }

    return tm;
}

template <class Pre, class Vec, class Comm>
double timeMinLoopFS(Pre& fs, Vec& x, Comm cc, int loopSize)
{
    Dune::Timer timer;
    Vec y(x.size());
    y = 0;
    double times[loopSize];
    for (int i = 0; i < loopSize; ++i) {
	cc.barrier();
	timer.reset();
	timer.start();
	fs.preconditioner().apply(y, x);
	times[i] = timer.stop();
	
    }

    double tm = times[0];
    for (int i = 1; i < loopSize; ++i) { 
	if (times[i]<tm) {
	    tm = times[i];	    
	}
    }

    return tm;
}

template<class Comm, class O, class Vec>
void multipleMinLoopTimePre(Comm cc, O& pre, Vec& x, int I=10, bool gl=false, bool dun=false, bool cpr=false)
{
    int rank = cc.rank();

    Vec y(x.size());
    y=0;
    pre.pre(y,x);
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t1 = timeMinLoopPre(pre, x, cc, 100);
	double times1[cc.size()];
	
	cc.gather(&t1, times1, 1, 0);
	if (rank==0)
	{
	    if (dun)
		std::cout << "MeasDUNEpre+1 " << cc.size() << ": ";
	    else if (cpr)
		std::cout << "MeasCPRpre+1 " << cc.size() << ": ";
	    else
		std::cout << "MeasPre+1 " << cc.size() << ": ";
	    for (int i = 0; i < cc.size(); i++)
		std::cout << times1[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template<class Comm, class O, class Vec>
void multipleMinLoopTimePre(Comm cc, O& pre, Vec& x, std::string pcname, int I=10)
{
    int rank = cc.rank();

    Vec y(x.size());
    y=0;
    pre.pre(y,x);
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t1 = timeMinLoopPre(pre, x, cc, 100);
	double times1[cc.size()];
	
	cc.gather(&t1, times1, 1, 0);
	if (rank==0)
	{
	    
	    std::cout << "Meas"<< pcname.c_str()<<"+1 " << cc.size() << ": ";
	    for (int i = 0; i < cc.size(); i++)
		std::cout << times1[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template<class Comm, class O, class Vec>
void multipleMinLoopTimeFS(Comm cc, O fs, Vec& x, std::string pcname, int I=10, bool gl=false, bool dun=false)
{
    int rank = cc.rank();

    Vec y(x.size());
    y=0;
    //auto pre = fs.preconditioner();
    fs.preconditioner().pre(y,x);
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t1 = timeMinLoopFS(fs, x, cc, 100);
	double times1[cc.size()];
	
	cc.gather(&t1, times1, 1, 0);
	if (rank==0)
	{
	    std::cout << "Meas" << pcname.c_str()<<"+1 "<< cc.size() << ": ";
	    for (int i = 0; i < cc.size(); i++)
		std::cout << times1[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template<class Vec, class Comm>
double timeMinLoopAxpy(Vec& x, Comm& cc, int loopSize)
{
    Dune::Timer timer;
    Vec y(x.size());
    Vec z=x;
    y = 1.0;
    
    double alpha = 0.23;

    double times[loopSize];
      
    for (int i = 0; i < loopSize; ++i) {
	cc.barrier();
	timer.reset();
	timer.start();
	z.axpy(alpha, y);
	times[i] = timer.stop();
    }

    double tm = times[0];
    for (int i = 1; i < loopSize; ++i) { 
	if (times[i]<tm) {
	    tm = times[i];
	    
	}
    }
    return tm;
}

template<class Comm, class Vec>
void multipleMinLoopTimeAxpy(Comm cc, Vec& x, int I=10, bool gl=false)
{
    int rank = cc.rank();
    
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t1 = timeMinLoopAxpy(x, cc, 100);
	double times1[cc.size()];
	
	cc.gather(&t1, times1, 1, 0);

	if (rank==0)
	{
	    if (gl)
		std::cout << "MeasGLaxpy+1 " << cc.size() << ": ";
	    else
		std::cout << "Measaxpy+1 " << cc.size() << ": ";
	    for (int i = 0; i < cc.size(); i++)
		std::cout << times1[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template<class Comm>
double timeMinLoopSimpleAdd(double* x, double* y, size_t size, Comm& cc, int loopSize)
{
    Dune::Timer timer;
    double alpha = 0.23;

    double times[loopSize];
      
    for (int i = 0; i < loopSize; ++i) {
	cc.barrier();
	timer.reset();
	timer.start();

	for (size_t i = 0; i< size; ++i)
	    x[i] += alpha*y[i];

	times[i] = timer.stop();
    }

    double tm = times[0];
    for (int i = 1; i < loopSize; ++i) { 
	if (times[i] < tm) {
	    tm = times[i];
	    
	}
    }
    return tm;
}

template<class Comm, class Vec>
void multipleMinLoopSimpleAdd(Comm cc, Vec& x, size_t size, int I=10, bool gl=false)
{
    int rank = cc.rank();

    const auto block_size = Vec::block_type::dimension;
    
    std::vector<double> x_(size*block_size);
    std::vector<double> y(size*block_size, 1.0);


    size_t bb = 0;
    for (size_t ii = 0; ii < size; ++ii) {
	for (size_t jj = 0; jj < block_size; ++jj) {
	    x_[bb] = x[ii][jj];
	    bb++;
	}
    }
    
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t1 = timeMinLoopSimpleAdd(x_.data(), y.data(), block_size*size, cc, 100);
	double times1[cc.size()];
	
	cc.gather(&t1, times1, 1, 0);

	if (rank==0)
	{
	    if (gl)
		std::cout << "MeasGLvecAdd " << cc.size() << ": ";
	    else
		std::cout << "MeasvectorSimpleAdd+1 " << cc.size() << ": ";
	    for (int i = 0; i < cc.size(); i++)
		std::cout << times1[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template<class Vec, class Comm, class C>
double timeMinLoopDummyCommBuild(Vec& x, const C& comm, Comm& cc, int loopSize)
{
    Dune::Timer timer;

    double times[loopSize];
      
    for (int i = 0; i < loopSize; ++i) {
	cc.barrier();
	timer.reset();
	timer.start();
	comm.buildCommunicator(x);
	times[i] = timer.stop();
    }

    double tm = times[0];
    for (int i = 1; i < loopSize; ++i) { 
	//std::cout << times[i] << " "<< cc.rank()<<std::endl;
	if (times[i]<tm) {
	    tm = times[i];
	    
	}
    }

    return tm;
}

template<class Comm, class C, class Vec>
void multipleMinLoopTimeDummyCommBuild(Comm cc, const C& comm, Vec& x, int I=10)
{
    int rank = cc.rank();
    
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t1 = timeMinLoopDummyCommBuild(x, comm, cc, 100);
	double times1[cc.size()];
	
	cc.gather(&t1, times1, 1, 0);

	if (rank==0)
	{
	    std::cout << "MeasDummyBuildC+1 " << cc.size() << ": ";
	    for (int i = 0; i < cc.size(); i++)
		std::cout << times1[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template<class CC, class Comm, class Vec, class SP, class GSP>
void measureStandardOperations(const CC& cc, const Comm& comm, Vec rhs, SP sp, GSP gsp,
			       int numCells, int iSize)
{
    int rank = cc.rank();

    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimeAxpy(cc, rhs);

    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopSimpleAdd(cc, rhs, numCells);

    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopSimpleAdd(cc, rhs, iSize, 10 , true);

    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimeComm(cc, comm, rhs);

    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimeSP(cc, sp, rhs);

    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimeSP(cc, gsp, rhs, 10 , true);
}

template<class CC, class Vec, class WellOp, class WellAdj, class NoWell>
void measureWellSpMV(const CC& cc, Vec rhs, WellOp& wop, WellAdj& aop, NoWell& nwop)
{
    int rank = cc.rank();

    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimeSpMV(cc, wop, rhs, 10, true, 2);

    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimeSpMV(cc, aop, rhs, 10, true, 0);

    if (rank == 0) {std::cout << std::endl;}
    multipleMinLoopTimeSpMV(cc, nwop, rhs, 10, true, 1);

}

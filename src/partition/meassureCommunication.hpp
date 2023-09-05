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

#ifndef OPM_MEASSURECOMMUNICATION_HEADER_INCLUDED
#define OPM_MEASSURECOMMUNICATION_HEADER_INCLUDED

#endif // OPM_MEASSURECOMMUNICATION_HEADER_INCLUDED

template<class C, class PI, class Comm>
double meassureComm(C cc, Comm& comm, PI info, int vecSize, bool barrier1, bool barrier2)
{
    //typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
    typedef Dune::FieldVector<double,3>                  BlockVec;
    typedef Dune::BlockVector<BlockVec>                  Vec;
        
    int rank = cc.rank();
    //Comm comm(info, cc);

    Vec x(vecSize);
    double value = rank + 0.112;
    x = value;
    
    Dune::Timer timer;
    if (barrier1)
	cc.barrier();
    
    timer.start();
    comm.copyOwnerToAll(x, x);
    if(barrier2)
	cc.barrier();
    double t = timer.stop();
    
    //std::cout <<"Rank " << rank << " time: " << t << std::endl;
    return t;
}

template<class C, class O, class Vec>
double timeSpMV(O& linOp, Vec& x, C cc, bool barrier1, bool barrier2)
{        
    Dune::Timer timer;    
        
    Vec y(x.size());
    y = 0;

    if (barrier1)
	cc.barrier();
    
    timer.start();
    linOp.apply(x, y);
    if(barrier2)
	cc.barrier();
    double t = timer.stop();
    
    return t;
}

template<class C, class Vec, class Pre>
double timePre(Pre& pre, Vec& x, C cc, bool barrier1, bool barrier2)
{        
    Dune::Timer timer;    
        
    Vec y(x.size());
    y = 0;

    if (barrier1)
	cc.barrier();
    
    timer.start();
    pre.apply(y, x);
    if(barrier2)
	cc.barrier();
    double t = timer.stop();
    
    return t;
}

template<class C, class Pre, class Vec>
void multipleTimePre(C cc, Pre& pre, Vec& x, int I=10)
{
    int rank = cc.rank();

    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = timePre(pre, x, cc, false, false);
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "MeasPre+1 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
    
    if (rank==0)
	std::cout << std::endl;
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = timePre(pre, x, cc, true, false);
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "MeasPre+2 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template<class C, class Vec, class Solver, class Stat>
double timeSolver(Solver& s, Vec& x, C cc, Stat ss,bool barrier1, bool barrier2)
{        
    Dune::Timer timer;    
        
    Vec y(x.size());
    y = 0;

    if (barrier1)
	cc.barrier();
    
    timer.start();
    s.apply(y, x, ss);
    if(barrier2)
	cc.barrier();
    double t = timer.stop();
    
    return t;
}

template<class C, class S, class Vec>
void multipleTimeSolver(C cc, S& s, Vec& x, int I=10)
{
    int rank = cc.rank();
    Dune::InverseOperatorResult stat;

    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = timeSolver(s, x, cc, stat, false, false);
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "Measbicg+1 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
    
    if (rank==0)
	std::cout << std::endl;
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = timeSolver(s, x, cc, stat, true, false);
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "Measbicg+2 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
}


template<class C, class W, class Vec>
double timeWell(W& well, Vec& x, C cc, bool barrier1, bool barrier2)
{        
    Dune::Timer timer;    
        
    Vec y(x.size());
    y = 0;

    if (barrier1)
	cc.barrier();
    
    timer.start();
    well.apply(x, y);
    if(barrier2)
	cc.barrier();
    double t = timer.stop();
    
    return t;
}

template<class C, class W, class Vec>
void multipleTimeWell(C cc, W& well, Vec& x, int I=10)
{
    int rank = cc.rank();        

    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = timeWell(well, x, cc, false, false);    
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "MeasWtime+1 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
    
    if (rank==0)
	std::cout << std::endl;
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = timeWell(well, x, cc, true, false);
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "MeasWtime+2 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template<class C, class PI, class Mat, class Vec>
double timeUMV(Mat& A, Vec& x, C cc, PI info, bool barrier1, bool barrier2)
{        
    Dune::Timer timer;    
    
    Vec y(A.N());
    y = 0;

    if (barrier1)
	cc.barrier();
    
    timer.start();    
    A.umv(x, y);
    if(barrier2)
	cc.barrier();
    double t = timer.stop();
        
    return t;
}

template<class C, class Mat, class Vec>
double timeRecUMV(Mat& A, Vec& x, C cc, bool barrier1, bool barrier2)
{        
    Dune::Timer timer;
            
    Vec y(A.N());
    y = 0;

    if (barrier1)
	cc.barrier();
    
    timer.start();    
    A.umv(x, y);
    if(barrier2)
	cc.barrier();
    double t = timer.stop();
    
    return t;
}


template<class C, class Vec, class SP>
double timeScalarProduct(Vec& x, SP& sp, C cc, bool barrier1, bool barrier2)
{
    Vec y(x.size());
    y = 1.0;
    
    Dune::Timer timer;
    if (barrier1)
	cc.barrier();
    
    timer.start();
    double val = sp.dot(x, y);
    if(barrier2)
	cc.barrier();
    double t = timer.stop();
    
    //std::cout <<"Rank " << rank << " time: " << t << std::endl;
    return t;
}

template<class SP, class C, class Vec>
void multipleTimeScalarProduct2(C cc, Vec& x, SP sp, int I=10)
{
    int rank = cc.rank();
    
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = timeScalarProduct(x, sp, cc, false, false);    
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "MeasGLsp+1 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
    
    if (rank==0)
	std::cout << std::endl;
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = timeScalarProduct(x, sp, cc, true, false);    
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "MeasGLsp+2 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template<class C, class PI, class Vec>
void multipleTimeScalarProduct(C cc, PI info, Vec& x, int I=10)
{
    int rank = cc.rank();
    typedef Dune::OwnerOverlapCopyCommunication<int,int>    Comm;
    typedef Dune::OverlappingSchwarzScalarProduct<Vec,Comm> ScalarProduct; 
    
    Comm comm(info, cc);
    ScalarProduct sp(comm);
    
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = timeScalarProduct(x, sp, cc, false, false);    
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "MeasSP1 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
    
    if (rank==0)
	std::cout << std::endl;
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = timeScalarProduct(x, sp, cc, true, false);    
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "MeasSP2 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template<class C, class PI, class Mat, class Vec>
void multipleTimeSpMV(C cc, PI info, Mat& A, Vec& x, int I=10)
{
    int rank = cc.rank();
    
    typedef Dune::OwnerOverlapCopyCommunication<int,int>       Comm;
    typedef Dune::OverlappingSchwarzOperator<Mat,Vec,Vec,Comm> Operator;
    
    Comm comm(info, cc);
    Operator linOp(A, comm);

    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = timeSpMV(linOp, x, cc, false, false);    
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "MeasSpMV+1 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
    
    if (rank==0)
	std::cout << std::endl;
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = timeSpMV(linOp, x, cc, true, false);    
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "MeasSpMV+2 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template<class C, class O, class Vec>
void multipleTimeSpMV2(C cc, O linOp, Vec& x, int I=10)
{
    int rank = cc.rank();
    
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = timeSpMV(linOp, x, cc, false, false);    
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "MeasGLapply+1 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
    
    if (rank==0)
	std::cout << std::endl;
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = timeSpMV(linOp, x, cc, true, false);    
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "MeasGLapply+2 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template<class C, class PI, class Vec, class O>
void multipleTimeWellSpMV(C cc, PI info, O& linOp, Vec& x, int I=10)
{
    int rank = cc.rank();        

    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = timeSpMV(linOp, x, cc, false, false);    
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "MeasWellOp+1 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
    
    if (rank==0)
	std::cout << std::endl;
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = timeSpMV(linOp, x, cc, true, false);    
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "MeasWellOp+2 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template<class C, class PI, class Mat, class Vec>
void multipleTimeUMV(C cc, PI info, Mat& A, Vec& x, int I=10)
{
    int rank = cc.rank();
    
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = timeUMV(A, x, cc, info, false, false);    
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "MeasUMV+1 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
    
    if (rank==0)
	std::cout << std::endl;
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = timeUMV(A, x, cc, info, true, false);    
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "MeasUMV+2 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template<class C, class PI, class Mat, class Vec>
void multipleTimeRec(C cc, PI info, Mat& A, Vec& x, int I=10)
{
    int rank = cc.rank();
    
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = timeRecUMV(A, x, cc, false, false);    
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "MeasREC+1 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
    
    if (rank==0)
	std::cout << std::endl;
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = timeRecUMV(A, x, cc, true, false);    
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "MeasREC+2 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template<class C, class PI>
void multipleMeassureComm(C cc, PI info, int vecSize, int I=10)
{
    int rank = cc.rank();
    
    typedef Dune::OwnerOverlapCopyCommunication<int,int> Comm;
    Comm comm(info, cc);

    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = meassureComm(cc, comm, info, vecSize, false, false);    
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "MeasComm1 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
    
    if (rank==0)
	std::cout << std::endl;
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t = meassureComm(cc, comm, info, vecSize, true, false);    
	double times[cc.size()];
	cc.gather(&t, times, 1, 0);
	
	if (rank==0)
	{
	    std::cout << "MeasComm2 " << cc.size() << ": ";
	    for (int i = 0; i<cc.size(); i++)
		std::cout << times[i] << " ";
	    std::cout << std::endl;
	}
    }
}

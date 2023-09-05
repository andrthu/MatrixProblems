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

#ifndef OPM_TIMEMULTIPLETIMES_HEADER_INCLUDED
#define OPM_TIMEMULTIPLETIMES_HEADER_INCLUDED

#endif // OPM_TIMEMULTIPLETIMES_HEADER_INCLUDED

template<class C, class O, class Vec>
double timeSpMVTimes10(O& linOp, Vec& x, C cc, bool barrier1)
{        
    Dune::Timer timer;
    Vec y(x.size());
    y = 0;

    if (barrier1)
	cc.barrier();
    
    timer.start();

    linOp.apply(x, y);
    linOp.apply(x, y);
    linOp.apply(x, y);
    linOp.apply(x, y);
    linOp.apply(x, y);
    linOp.apply(x, y);
    linOp.apply(x, y);
    linOp.apply(x, y);
    linOp.apply(x, y);
    linOp.apply(x, y);

    double t = timer.stop();
    return t/10.0;
}


template<class C, class Vec, class SP>
double timeScalarProductTimes10(Vec& x, SP& sp, C cc, bool barrier1)
{
    Vec y(x.size());
    y = 1.0;
    
    Dune::Timer timer;
    if (barrier1)
	cc.barrier();
    
    double val = 0;
    timer.start();
    
    val = sp.dot(x, y);
    val = sp.dot(x, y);
    val = sp.dot(x, y);
    val = sp.dot(x, y);
    val = sp.dot(x, y);
    val = sp.dot(x, y);
    val = sp.dot(x, y);
    val = sp.dot(x, y);
    val = sp.dot(x, y);
    val = sp.dot(x, y);

    double t = timer.stop();    
    return t/10.0;
}

template<class C, class Vec, class Pre>
double timePreTimes10(Pre& pre, Vec& x, C cc, bool barrier1)
{        
    Dune::Timer timer;    
        
    Vec y(x.size());
    y = 0;

    if (barrier1)
	cc.barrier();
    
    timer.start();

    pre.apply(y, x);
    pre.apply(y, x);
    pre.apply(y, x);
    pre.apply(y, x);
    pre.apply(y, x);
    pre.apply(y, x);
    pre.apply(y, x);
    pre.apply(y, x);
    pre.apply(y, x);
    pre.apply(y, x);

    double t = timer.stop();
    return t/10.0;
}

template<class C, class Pre, class Vec>
void multipleAverageTimePre(C cc, Pre& pre, Vec& x, int I=10)
{
    int rank = cc.rank();

    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t1 = timePreTimes10(pre, x, cc, false);
	double t2 = timePreTimes10(pre, x, cc, true);
	double times1[cc.size()];
	double times2[cc.size()];
	cc.gather(&t1, times1, 1, 0);
	cc.gather(&t2, times2, 1, 0);

	if (rank==0)
	{
	    std::cout << "MeasPre+1 " << cc.size() << ": ";
	    for (int i = 0; i < cc.size(); i++)
		std::cout << times1[i] << " ";
	    std::cout << std::endl;

	    std::cout << "MeasPre+2 " << cc.size() << ": ";
	    for (int i = 0; i < cc.size(); i++)
		std::cout << times2[i] << " ";
	    std::cout << std::endl;
	}
    }
}


template<class C, class O, class Vec>
void multipleAverageTimeSpMV(C cc, O linOp, Vec& x, int I=10, bool gl=false)
{
    int rank = cc.rank();
    
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t1 = timeSpMVTimes10(linOp, x, cc, false);
	double t2 = timeSpMVTimes10(linOp, x, cc, true);

	double times1[cc.size()];
	double times2[cc.size()];
	cc.gather(&t1, times1, 1, 0);
	cc.gather(&t2, times2, 1, 0);

	if (rank==0)
	{
	    if (gl)
		std::cout << "MeasGLapply+1 " << cc.size() << ": ";
	    else
		std::cout << "MeasSpMV+1 " << cc.size() << ": ";
	    for (int i = 0; i < cc.size(); i++)
		std::cout << times1[i] << " ";
	    std::cout << std::endl;

	    if (gl)
		std::cout << "MeasGLapply+2 " << cc.size() << ": ";
	    else
		std::cout << "MeasSpMV+2 " << cc.size() << ": ";
	    for (int i = 0; i < cc.size(); i++)
		std::cout << times2[i] << " ";
	    std::cout << std::endl;
	}
    }
}

template<class C, class SP, class Vec>
void multipleAverageTimeScalarProduct(C cc, SP sp, Vec& x, int I=10, bool gl=false)
{
    int rank = cc.rank();
    
    for (int j = 0; j < I; ++j)
    {
	cc.barrier();
	double t1 = timeScalarProductTimes10(x, sp, cc, false);
	double t2 = timeScalarProductTimes10(x, sp, cc, true);

	double times1[cc.size()];
	double times2[cc.size()];
	cc.gather(&t1, times1, 1, 0);
	cc.gather(&t2, times2, 1, 0);

	if (rank==0)
	{
	    if (gl)
		std::cout << "MeasGLsp+1 " << cc.size() << ": ";
	    else
		std::cout << "MeasSP+1 " << cc.size() << ": ";
	    for (int i = 0; i < cc.size(); i++)
		std::cout << times1[i] << " ";
	    std::cout << std::endl;

	    if (gl)
		std::cout << "MeasGLsp+2 " << cc.size() << ": ";
	    else
		std::cout << "MeasSP+2 " << cc.size() << ": ";
	    for (int i = 0; i < cc.size(); i++)
		std::cout << times2[i] << " ";
	    std::cout << std::endl;
	}
    }
}

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
        double t1 = timeMinLoopSpMVAS(linOp, x, cc, 50);
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

// Same as timeMinLoopSpMVAS, but without the cc.barrier() call that
// re-synchronizes all ranks before every single call to applyscaleadd. See
// the comment above timeMinLoopCommNoBarrier for the rationale: the barrier
// inside the loop hides waiting/imbalance cost under the "communication"
// label instead of attributing it to the operation itself. Removing it
// gives a timing closer to what the operation costs "in the wild", at the
// price of being noisier.
template<class Op, class Vec, class Comm>
double timeMinLoopSpMVASNoBarrier(Op& o, Vec& x, Comm& cc, int loopSize)
{
    Dune::Timer timer;
    Vec y(x.size());

#pragma omp parallel for
    for (size_t j = 0; j < x.size(); ++j)
        y[j] = 0;

    double times[loopSize];

    for (int i = 0; i < loopSize; ++i) {
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
void multipleMinLoopTimeSpMVASNoBarrier(Comm cc, O& linOp, Vec& x, int I=10, bool gl=false, int well=0)
{
    int rank = cc.rank();

    for (int j = 0; j < I; ++j)
    {
        double t1 = timeMinLoopSpMVASNoBarrier(linOp, x, cc, 50);
        double times1[cc.size()];

        // cc.gather is itself collective, so it still bounds when rank 0
        // can print, but unlike cc.barrier() it does not force all ranks
        // to rendezvous before *every* timed call above.
        cc.gather(&t1, times1, 1, 0);

        if (rank==0)
        {
            if (gl)
                std::cout << "MeasSpMV-GLNoBarrier+1 " << cc.size() << ": ";
            else
                std::cout << "MeasSpMVNoBarrier+1 " << cc.size() << ": ";
            for (int i = 0; i < cc.size(); i++)
                std::cout << times1[i] << " ";
            std::cout << std::endl;
        }
    }
}

// A third way of measuring SpMVAS (applyscaleadd), aimed at giving a more
// trustworthy number than either of the two variants above. See the
// comment above timeMinLoopCommGood for the full rationale:
//  - An untimed warm-up call is made before the loop starts.
//  - Ranks are synchronized with a single cc.barrier() once, right before
//    the timed loop, and there is no barrier inside the loop.
//  - Both the minimum and the mean of the loopSize samples are recorded
//    per rank, and both are reduced across ranks with cc.max() since the
//    cost felt by the application is bounded by the slowest rank.
template<class Op, class Vec, class Comm>
double timeMinLoopSpMVASGood(Op& o, Vec& x, Comm& cc, int loopSize, double& meanOut)
{
    Dune::Timer timer;
    Vec y(x.size());

#pragma omp parallel for
    for (size_t j = 0; j < x.size(); ++j)
        y[j] = 0;

    // Untimed warm-up call.
    o.applyscaleadd(-1, x, y);

    // Align all ranks once, right before timing starts.
    cc.barrier();

    double times[loopSize];

    for (int i = 0; i < loopSize; ++i) {
        timer.reset();
        timer.start();
        o.applyscaleadd(-1, x, y);
        times[i] = timer.stop();
    }

    double tm = times[0];
    double sum = 0.0;
    for (int i = 0; i < loopSize; ++i) {
        sum += times[i];
        if (times[i]<tm) {tm = times[i];}
    }
    meanOut = sum / loopSize;

    return tm;
}

template<class Comm, class O, class Vec>
void multipleMinLoopTimeSpMVASGood(Comm cc, O& linOp, Vec& x, int I=10, bool gl=false, int well=0)
{
    int rank = cc.rank();

    for (int j = 0; j < I; ++j)
    {
        cc.barrier();
        double mean1 = 0.0;
        double min1 = timeMinLoopSpMVASGood(linOp, x, cc, 50, mean1);

        double minTimes1[cc.size()];
        double meanTimes1[cc.size()];
        cc.gather(&min1, minTimes1, 1, 0);
        cc.gather(&mean1, meanTimes1, 1, 0);

        if (rank==0)
        {
            std::string tag = gl ? "MeasSpMV-GL" : "MeasSpMV";

            std::cout << tag << "GoodMin+1 " << cc.size() << ": ";
            for (int i = 0; i < cc.size(); i++)
                std::cout << minTimes1[i] << " ";
            std::cout << std::endl;

            std::cout << tag << "GoodMean+1 " << cc.size() << ": ";
            for (int i = 0; i < cc.size(); i++)
                std::cout << meanTimes1[i] << " ";
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

// Same as timeMinLoopComm, but without the cc.barrier() call that
// re-synchronizes all ranks before every single call to copyOwnerToAll.
// Since copyOwnerToAll is itself a (point-to-point) communication
// primitive that synchronizes the ranks it exchanges data with, adding a
// barrier before each call forces *all* ranks - including ones that do
// not directly exchange data with each other - to rendezvous first. That
// hides any waiting/imbalance cost that would normally be attributed to
// the communication step, and it also means the timer partially measures
// how long the *other* ranks needed to catch up to the barrier rather
// than the cost of the operation itself. Removing the barrier gives a
// timing that is closer to what the operation costs "in the wild", at
// the price of being noisier and more sensitive to what ranks were doing
// just before the call.
template<class Vec, class Comm, class C>
double timeMinLoopCommNoBarrier(Vec& x, const C& comm, Comm& cc, int loopSize)
{
    Dune::Timer timer;
    Vec y=x;

    double times[loopSize];

    for (int i = 0; i < loopSize; ++i) {
        timer.reset();
        timer.start();
        comm.copyOwnerToAll(y, y);
        times[i] = timer.stop();

        y.axpy(0.01, x);
    }

    double tm = times[0];
    for (int i = 1; i < loopSize; ++i) {
        if (times[i]<tm) {
            tm = times[i];
        }
    }

    return tm;
}

template<class Comm, class C, class Vec>
void multipleMinLoopTimeCommNoBarrier(Comm cc, const C& comm, Vec& x, int I=10)
{
    int rank = cc.rank();

    for (int j = 0; j < I; ++j)
    {
        double t1 = timeMinLoopCommNoBarrier(x, comm, cc, 100);
        double times1[cc.size()];

        // cc.gather is itself collective, so it still bounds when rank 0
        // can print, but unlike cc.barrier() it does not force all ranks
        // to rendezvous before *every* timed call above.
        cc.gather(&t1, times1, 1, 0);

        if (rank==0)
        {
            std::cout << "MeasCommNoBarrier+1 " << cc.size() << ": ";
            for (int i = 0; i < cc.size(); i++)
                std::cout << times1[i] << " ";
            std::cout << std::endl;
        }
    }
}

// A third way of measuring copyOwnerToAll, aimed at giving a more
// trustworthy number than either of the two variants above:
//
//  - A warm-up call (untimed) is made before the loop starts, so that
//    one-off costs (e.g. lazy allocation of communication buffers inside
//    DUNE/MPI on the first call) do not pollute the measurements.
//  - Ranks are synchronized with a single cc.barrier() *once*, right
//    before the timed loop starts, so all ranks begin the loop on an
//    equal footing. There is no barrier *inside* the loop: copyOwnerToAll
//    already synchronizes the ranks that exchange data, so repeatedly
//    barrier-syncing everyone in between calls (as multipleMinLoopTimeComm
//    does) is both unnecessary and actively distorts the measurement, as
//    explained above.
//  - Per rank we record both the minimum (best-case, cache-warm cost) and
//    the mean (typical cost, including any waiting caused by load
//    imbalance) of the loopSize samples, since the minimum alone can
//    hide systematic imbalance between ranks.
//  - Because the operation is collective in effect (every rank has to
//    wait for its slowest communication partner), the cost actually felt
//    by the application is bounded by the *slowest* rank, not by any
//    single rank's local timing. We therefore also reduce both the min-
//    and mean-times across ranks with cc.max(), in addition to gathering
//    the per-rank numbers for inspection.
template<class Vec, class Comm, class C>
double timeMinLoopCommGood(Vec& x, const C& comm, Comm& cc, int loopSize, double& meanOut)
{
    Dune::Timer timer;
    Vec y=x;

    // Untimed warm-up call.
    comm.copyOwnerToAll(y, y);

    // Align all ranks once, right before timing starts.
    cc.barrier();

    double times[loopSize];

    for (int i = 0; i < loopSize; ++i) {
        timer.reset();
        timer.start();
        comm.copyOwnerToAll(y, y);
        times[i] = timer.stop();

        y.axpy(0.01, x);
    }

    double tm = times[0];
    double sum = 0.0;
    for (int i = 0; i < loopSize; ++i) {
        sum += times[i];
        if (times[i]<tm) {
            tm = times[i];
        }
    }
    meanOut = sum / loopSize;

    return tm;
}

template<class Comm, class C, class Vec>
void multipleMinLoopTimeCommGood(Comm cc, const C& comm, Vec& x, int I=10)
{
    int rank = cc.rank();

    for (int j = 0; j < I; ++j)
    {
        cc.barrier();
        double mean1 = 0.0;
        double min1 = timeMinLoopCommGood(x, comm, cc, 100, mean1);

        double minTimes1[cc.size()];
        double meanTimes1[cc.size()];
        cc.gather(&min1, minTimes1, 1, 0);
        cc.gather(&mean1, meanTimes1, 1, 0);

        if (rank==0)
        {
            std::cout << "MeasCommGoodMin+1 " << cc.size() << ": ";
            for (int i = 0; i < cc.size(); i++)
                std::cout << minTimes1[i] << " ";
            std::cout << std::endl;

            std::cout << "MeasCommGoodMean+1 " << cc.size() << ": ";
            for (int i = 0; i < cc.size(); i++)
                std::cout << meanTimes1[i] << " ";
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
        double t1 = timeMinLoopFS(fs, x, cc, 50);
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

// Same as timeMinLoopFS, but without the cc.barrier() call that
// re-synchronizes all ranks before every single call to
// preconditioner().apply(). See the comment above timeMinLoopCommNoBarrier
// for the rationale: the barrier inside the loop hides waiting/imbalance
// cost under the "operation" label instead of attributing it correctly.
// Removing it gives a timing closer to what the operation costs "in the
// wild", at the price of being noisier.
template <class Pre, class Vec, class Comm>
double timeMinLoopFSNoBarrier(Pre& fs, Vec& x, Comm cc, int loopSize)
{
    Dune::Timer timer;
    Vec y(x.size());
    y = 0;
    double times[loopSize];
    for (int i = 0; i < loopSize; ++i) {
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
void multipleMinLoopTimeFSNoBarrier(Comm cc, O fs, Vec& x, std::string pcname, int I=10, bool gl=false, bool dun=false)
{
    int rank = cc.rank();

    Vec y(x.size());
    y=0;
    fs.preconditioner().pre(y,x);
    for (int j = 0; j < I; ++j)
    {
        double t1 = timeMinLoopFSNoBarrier(fs, x, cc, 50);
        double times1[cc.size()];

        // cc.gather is itself collective, so it still bounds when rank 0
        // can print, but unlike cc.barrier() it does not force all ranks
        // to rendezvous before *every* timed call above.
        cc.gather(&t1, times1, 1, 0);
        if (rank==0)
        {
            std::cout << "Meas" << pcname.c_str()<<"NoBarrier+1 "<< cc.size() << ": ";
            for (int i = 0; i < cc.size(); i++)
                std::cout << times1[i] << " ";
            std::cout << std::endl;
        }
    }
}

// A third way of measuring the preconditioner apply, aimed at giving a
// more trustworthy number than either of the two variants above. See the
// comment above timeMinLoopCommGood for the full rationale:
//  - An untimed warm-up call is made before the loop starts.
//  - Ranks are synchronized with a single cc.barrier() once, right before
//    the timed loop, and there is no barrier inside the loop.
//  - Both the minimum and the mean of the loopSize samples are recorded
//    per rank, and both are reduced across ranks with cc.max() since the
//    cost felt by the application is bounded by the slowest rank.
template <class Pre, class Vec, class Comm>
double timeMinLoopFSGood(Pre& fs, Vec& x, Comm cc, int loopSize, double& meanOut)
{
    Dune::Timer timer;
    Vec y(x.size());
    y = 0;

    // Untimed warm-up call.
    fs.preconditioner().apply(y, x);

    // Align all ranks once, right before timing starts.
    cc.barrier();

    double times[loopSize];
    for (int i = 0; i < loopSize; ++i) {
        timer.reset();
        timer.start();
        fs.preconditioner().apply(y, x);
        times[i] = timer.stop();
    }

    double tm = times[0];
    double sum = 0.0;
    for (int i = 0; i < loopSize; ++i) {
        sum += times[i];
        if (times[i]<tm) {
            tm = times[i];
        }
    }
    meanOut = sum / loopSize;

    return tm;
}

template<class Comm, class O, class Vec>
void multipleMinLoopTimeFSGood(Comm cc, O fs, Vec& x, std::string pcname, int I=10)
{
    int rank = cc.rank();

    Vec y(x.size());
    y=0;
    fs.preconditioner().pre(y,x);
    for (int j = 0; j < I; ++j)
    {
        cc.barrier();
        double mean1 = 0.0;
        double min1 = timeMinLoopFSGood(fs, x, cc, 50, mean1);

        double minTimes1[cc.size()];
        double meanTimes1[cc.size()];
        cc.gather(&min1, minTimes1, 1, 0);
        cc.gather(&mean1, meanTimes1, 1, 0);

        if (rank==0)
        {
            std::cout << "Meas" << pcname.c_str() << "GoodMin+1 " << cc.size() << ": ";
            for (int i = 0; i < cc.size(); i++)
                std::cout << minTimes1[i] << " ";
            std::cout << std::endl;

            std::cout << "Meas" << pcname.c_str() << "GoodMean+1 " << cc.size() << ": ";
            for (int i = 0; i < cc.size(); i++)
                std::cout << meanTimes1[i] << " ";
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

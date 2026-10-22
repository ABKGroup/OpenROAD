/**************************************************************************
***    
*** Copyright (c) 1995-2001 Regents of the University of California,
***               Ken Boese, Sun Cho, Andrew E. Caldwell, Andrew B. Kahng,
***               Igor L. Markov and Mike Oliver
***
***  Contact author(s): abk@cs.ucsd.edu, oliver@cs.ucla.edu
***  Affiliation:   UCLA, Computer Science Department,
***                 Los Angeles, CA 90095-1596 USA
***
***  Permission is hereby granted, free of charge, to any person obtaining 
***  a copy of this software and associated documentation files (the
***  "Software"), to deal in the Software without restriction, including
***  without limitation 
***  the rights to use, copy, modify, merge, publish, distribute, sublicense, 
***  and/or sell copies of the Software, and to permit persons to whom the 
***  Software is furnished to do so, subject to the following conditions:
***
***  The above copyright notice and this permission notice shall be included
***  in all copies or substantial portions of the Software.
***
*** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
*** EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
*** OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
*** IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
*** CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
*** OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
*** THE USE OR OTHER DEALINGS IN THE SOFTWARE.
***
***
***************************************************************************/




//////////////////////////////////////////////////////////////////////
//
// optimizer.h: interface for the Optimizer class.
//              This is the base class for objects which optimize
//              ScanChain objects.
//////////////////////////////////////////////////////////////////////

#if !defined(_ABKSCANOPT_OPTMIZER_H)
#define _ABKSCANOPT_OPTMIZER_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <ScanChain/scanChain.h>
#include "path.h"

namespace abkscanopt
{
class Optimizer  
{
public:
    class Stats
    {
    public:
        double Cost;
        double Time;
        unsigned Count;
        unsigned MeasurementUnits;
        Stats():Cost(0.0),Time(0.0),Count(0),MeasurementUnits(0){}
    };

    enum OptimizerKind
    {
        OPTIMIZER1
    };

    class Params
    {
    public:
        bool          maintainStats;
        unsigned      statsInterval;
        OptimizerKind subKind;
        Params        *subParams; //makes a linked list; if you reach
                                  //a NULL, use default values.  Params
                                  //object does *not* own the memory (no
                                  //dtor) so use automatic variables

        //block for Optimizer1
        unsigned majorLoops;
        unsigned nDescents;
        unsigned kickMove;
        unsigned nnear; //number of nearest cells to consider
        bool     only2Opt;
        bool     temp_control;
        bool     zeroTemp;
        //end Optimizer1 params

        Params():maintainStats(false),statsInterval(1),
            subKind(OPTIMIZER1),subParams(NULL),majorLoops(100),
            kickMove(15),nnear(20),only2Opt(false),temp_control(false),
            zeroTemp(true){}

    };
protected:
    ScanChain         &_chain;   //scan chain to optimize
    Path               _bsfPath; //best path so far
    RandomRawUnsigned &_randuns; //RNG
    Params   _params;
    vector<Stats>     _stats,_sumStats,_sumSquareStats;
    void              _updateStats(unsigned idx,unsigned measurementUnits,
                                   double cost,double time);

    double     _cost;  //cost at this level and below

    //since Optimizer is a friend to ScanChain, all classes
    //derived from Optimizer can modify path using this method
    vector<unsigned> &_getChainPath(){return _chain._path;}
public:


    //named ctor; returns new object of derived type specified by "kind"
    static Optimizer *newOptimizer(ScanChain &chain,
        RandomRawUnsigned &randuns,Params params,
        OptimizerKind kind);
	Optimizer(ScanChain &chain,RandomRawUnsigned &randuns,
             Params params);
	virtual ~Optimizer();

    double getCost() const {return _cost;}

    virtual void rawStatsOut(::ostream &os) const=0;

};
};

#endif // !defined(_ABKSCANOPT_OPTMIZER_H)

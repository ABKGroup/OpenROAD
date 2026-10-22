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
// optimizer.cxx: implementation of the Optimizer class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include "optimizer.h"
#include "optimizer1.h"
using namespace abkscanopt;
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Optimizer::Optimizer(ScanChain &chain,RandomRawUnsigned &randuns,
                     Params params):
_chain(chain),_bsfPath(chain),_randuns(randuns),_params(params)
{
    if (!_chain.isPathValid())
    {
        _chain.initialPath(_randuns);
        _chain._fillInCellInfo();
        _bsfPath.save();
    }

}

Optimizer::~Optimizer()
{

}

Optimizer *Optimizer::newOptimizer(ScanChain &chain,RandomRawUnsigned &randuns,
                                   Params params,OptimizerKind kind)
{
    switch(kind)
    {
    case Optimizer::OPTIMIZER1:
        {
            return new Optimizer1(chain,randuns,params);
            break;
        }
    default:
        abkfatal(false,"Unknown optimizer kind");
    }
}

void  Optimizer:: _updateStats(unsigned idx,unsigned measurementUnits,
                                   double cost,double time)
{
    const unsigned numStats=_stats.size();
    abkfatal(numStats==_sumStats.size() && numStats==_sumSquareStats.size(),
        "Inconsistent stats sizes");
    if (idx+1>numStats)
    {
        unsigned i;
        Stats emptyStats;
        for (i=numStats;i<idx+1;i++)
        {
            _stats.push_back(emptyStats);
            _sumStats.push_back(emptyStats);
            _sumSquareStats.push_back(emptyStats);

        }
        Stats &stat=_stats[idx];
        stat.Count=1;
        stat.MeasurementUnits=measurementUnits;
        stat.Cost=cost;
        stat.Time=time;

        Stats &sumStat=_sumStats[idx];
        sumStat.Count += 1;
        sumStat.MeasurementUnits+=measurementUnits;
        sumStat.Cost+=cost;
        sumStat.Time +=time;

        Stats &sumSquareStat=_sumSquareStats[idx];
        sumSquareStat.Count += 1;
        sumSquareStat.MeasurementUnits+=measurementUnits*measurementUnits;
        sumSquareStat.Cost+=cost*cost;
        sumSquareStat.Time +=time*time;
    }
}



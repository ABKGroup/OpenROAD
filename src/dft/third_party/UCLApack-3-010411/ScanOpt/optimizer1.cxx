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
// optimizer1.cxx: implementation of the Optimizer1 class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include "optimizer1.h"
#include "scanTourDZ.h"
#include <float.h>
using namespace abkscanopt;
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Optimizer1::Optimizer1(ScanChain &chain,RandomRawUnsigned &randuns,
                       Params params):
Optimizer(chain,randuns,params),_lastPath(chain)
{
    unsigned i;
    unsigned g;
    unsigned loop;
    unsigned NewIn,NewOut;
    unsigned OldIn,OldOut;
    double Cost;
    double LevelCost;
    double bsfCost = DBL_MAX;
    double lastCost = DBL_MAX;
    const unsigned nSubgroups=_chain.getNumSubgroups();
    vector<double> SubGroupCost(nSubgroups,DBL_MAX);
    double ElapsedTime;
    double OldSubGroupCost;
    Params defaultParams;

    //Parameters for subgroups
    Params &subParams=(_params.subParams)?*_params.subParams:defaultParams;
    
    //#ifdef _DEBUG
    _chain.checkPath();
    //#endif
    
    Timer m_Timer;
    if (_params.maintainStats)
    {
        bsfCost = _bsfPath.cost(); //_bsfPath initialized in Optimizer ctor
        _updateStats(0,-1,bsfCost,0.0);
    }
    
    LevelCost = this->_optimizeLevel();
    
    //#ifdef _DEBUG
    _chain.checkPath();
    //#endif
    
    for (i=0;i<nSubgroups;i++)
    {
        Optimizer *pOpt=newOptimizer(*_chain.getCells()[i].subgroup,_randuns,
            subParams,_params.subKind);
        SubGroupCost[i]=pOpt->getCost();
        delete pOpt;
    }
    
    Cost = LevelCost;
    for (i=0;i<nSubgroups;i++)
        Cost += SubGroupCost[i];
    
    if (Cost<bsfCost)
    {
        _bsfPath.save();
        bsfCost = Cost;
    }
    
    if (Cost<lastCost && _params.zeroTemp)
    {
        _lastPath.save();
        lastCost = _lastPath.cost();
    }
    
    if (_params.maintainStats)
    {
        m_Timer.split();
        ElapsedTime = m_Timer.getCombTime();
        m_Timer.resume();
        _updateStats(1,0,bsfCost,ElapsedTime);
    }
    
    for (loop=1;loop<_params.majorLoops;loop++)
    {
        if (nSubgroups > 0)
        {
            // choose a group g and change its in, out cells
            
            g = (loop-1)%nSubgroups;
            OldSubGroupCost = SubGroupCost[g];
            ScanChain &subgroup=*_chain.getCells()[g].subgroup;
            unsigned nSubcells=subgroup.getNumCells();
            const vector<unsigned> &subpath=subgroup.getPath();
            OldIn = subpath[0];
            OldOut = subpath[nSubcells-1];
            
            // If group g has only one legal in cell and one
            // legal out cell, we can even skip it.
            
            // LOOK!  We haven't taken care of the situation where,
            //        say, we have one legal in cell and two legal
            //        out cells, but one of the legal out cells is also
            //        the only legal in cell.  In that case the program
            //        will hang.

            const vector<unsigned> &sublegalins=subgroup.getLegalins();
            unsigned nSublegalins=sublegalins.size();
            const vector<unsigned> &sublegalouts=subgroup.getLegalouts();
            unsigned nSublegalouts=sublegalouts.size();
            
            if (nSublegalins > 1 || nSublegalouts > 1)
            {
                do
                {
                    NewIn = sublegalins[_randuns%nSublegalins];
                    NewOut = sublegalouts[_randuns%nSublegalouts];
                } while 
                    (NewIn==NewOut || (NewIn==OldIn&&NewOut==OldOut));
                
                //takes care of changing inX,inY,outX,outY of
                //the corresponding cell, as well as changing the path
                _chain.subModifyPath(g,NewIn,NewOut);

                Optimizer *pOpt=newOptimizer(subgroup,_randuns,
                                             subParams,_params.subKind);
                SubGroupCost[g]=pOpt->getCost();
                delete pOpt;
            }
        }   //end "if more that one subgrp
        
        LevelCost=_optimizeLevel();
        Cost = LevelCost;
        
        for (i=0;i<nSubgroups;i++)
            Cost += SubGroupCost[i];
        
        if (Cost<bsfCost)
        {
            _bsfPath.save();
            bsfCost=Cost;
        }
        
        if (_params.zeroTemp)
        {
            if (Cost<lastCost)
            {
                _lastPath.save();
                lastCost = Cost;
            }
            else
            {
                _lastPath.restore();
#ifdef ABKDEBUG
                _chain.checkPath();
#endif
                Cost = _lastPath.cost();
                if (nSubgroups > 0)
                    SubGroupCost[g]=OldSubGroupCost;
            }
        }
                
        if (_params.maintainStats)
        {
            if (loop%_params.statsInterval == 0)
            {
                m_Timer.split();
                ElapsedTime = m_Timer.getCombTime();
                m_Timer.resume();
                
                _updateStats(loop/_params.statsInterval+1,loop,
                             bsfCost,ElapsedTime);
            }
            
        }
        
    } // end loop through number of "major loops"
    
    
    _chain.checkPath();
    _cost = Cost;

}

//***************************************************************************

Optimizer1::~Optimizer1()
{

}
//***************************************************************************
double Optimizer1::_optimizeLevel()
{

    ScanTourDZ::OptParams optParams;
    optParams.nDescents = _params.nDescents;
    optParams.kickMove = _params.kickMove;
    optParams.only2Opt = _params.only2Opt;
    optParams.temp_control = _params.temp_control;
    optParams.bZeroFix=true;
    unsigned collapseSize=_chain.getNumCells()-1;
    unsigned inIdx=_chain.getPath()[0];
    unsigned outIdx=_chain.getPath()[collapseSize];

    unsigned nnear=_params.nnear;
    if (nnear>collapseSize-1) nnear=collapseSize-1;

    ScanTourDZ tour(_chain,nnear,_randuns,optParams);

    //Really, if the chain is a subpath (i.e. fully specified
    //by POs) then we shouldn't even be creating the ScanTourDZ
    //object, which will take some time.  However the waste is
    //probably small, and we do need to get the cost of the level.
    //mro 04apr01
    if (!_chain.isSubpath())
    {

        tour.scoOptAll();


        unsigned lesser=inIdx,greater=outIdx;
        if (lesser>greater) std::swap(lesser,greater);

        const vector<unsigned> &collapsedPath=tour.getPath();
        unsigned i;
        for (i=1;i<collapseSize;i++)
        {
            _getChainPath()[i]=
                ScanCells::uncollapseIndex(collapsedPath[i],lesser,greater);
        }

        _chain.computeInverse();

    #ifdef ABKDEBUG
        _chain.checkPath();
    #endif
    }


    return tour.scoTourCost();

}
//***************************************************************************
void Optimizer1::rawStatsOut(::ostream &os) const
{
    unsigned i;
    const unsigned iMax=_stats.size();
    os << "Cost\tTime\tMajor Loops"<<endl;
    for (i=0;i<iMax;i++)
    {
        const Stats &rawstat=_stats[i];
        if (rawstat.Count!=0)
        {
            abkfatal(rawstat.Count==1,"Internal error in stats tracking");
            os << rawstat.Cost << "\t" << rawstat.Time << "\t"
                << rawstat.MeasurementUnits << endl;
        }
    }
}


//***************************************************************************
//***************************************************************************

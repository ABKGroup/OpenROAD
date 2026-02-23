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
// scanTourDistinguishedZero.cxx: implementation of the ScanTourDZ class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _MSC_VER
#pragma warning(disable:4786 4355)
#endif

#include "scanTourDZ.h"
#include "kenKDTree.h"
#include <float.h>
#define FUZZ      0.00001

#include <algorithm>

using namespace abkscanopt;
using namespace std;


#ifdef COUNT
long nTwoOptCount;
long nThreeOptCount;
long nFindTwoOptCount;
long nFindThreeOptCount;
long nRunTwoOptCount;
long nRunThreeOptCount;
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


ScanTourDZ::ScanTourDZ(ScanChain &origChain,
                     unsigned nnear,
                     RandomRawUnsigned &randuns,
                     OptParams optParams):
ScanChain(origChain,collapseInOut()),
_kdTree(*this,nnear),
_randuns(randuns),
_n(getNumCells()),
_optParams(optParams),
_ZeroShift(0),
//_bZero2Optable(_n+1,true),
_bZero2Optable(_n,true),
_zeroConstraintsValid(false),
_2optLimsValid(false),
_zero2optableValid(false),
_RightLim2Opt(_n,_n-1),
_scoRevGain(2*_n+2,DBL_MAX),
_temp(_n)

{
    _ComputeConstraintExistence();
}


ScanTourDZ::~ScanTourDZ()
{
    
}

void ScanTourDZ:: _scoInitRevGain()
{
    unsigned i;
    unsigned tt1, tt2;
    abkfatal(_pathValid,"Tour invalid");
    
    _scoRevGain[0] = 0.0;
    for (i=1; i<=_n; i++) 
    {
        tt1 = _path[i-1];
        tt2 = _path[i % _n];
        _scoRevGain[i] = _scoRevGain[i-1] + realDist(tt1,tt2)
            - realDist(tt2,tt1);
    }
}
void ScanTourDZ::_scoReversePath(unsigned t1,unsigned t3)
{
    vector<unsigned> tint(_n);
    unsigned i;
    int diff;
    abkfatal(_pathValid,"Invalid tour");
    abkfatal(_inverseValid,"Invalid inverse tour");
    vector<unsigned>::iterator tt = _path.begin();
    vector<unsigned>::iterator ti = _pathinv.begin();
    
    if ( t1 < t3) // note: if t1 was 0 and _optParams.bZeroFix is true,
        // now t1 is _n, so this cond will be
        //       false.
    {
        /* reverse the subpath between t1 & t3 */
        for (i=t1; i<=t3; i++)                       
            tint[i-t1] = tt[i];                      
        
        for (i=t1; i<=t3; i++)                       
        {                                        
            tt[i] = tint[t3-i];                      
            ti[tt[i]] = i;                           
        }  /* scoReversePath() */
    }
    else
    { // wraparound case.  This depends on the fact that
        // if t1 was originally 0, it has been reset to _n.
        
        //copy_n(tt,_n,tint.begin());
        copy(tt,tt+_n,tint.begin());
        
        for (i=1;i<=_n-t1;i++)
        {
            tt[i] = tint[_n-i];
            ti[tt[i]] = i;
        }
        
        //copy_n(tint.begin()+(t3+1),(t1-t3-1),tt+(_n-t1+1));
        copy(tint.begin()+(t3+1),tint.begin()+t1,tt+(_n-t1+1));
        diff = t1+t3-_n;
        
        for (i=t3+1;i<t1;i++) //LOOK! check this carefully
            ti[tint[i]] -= diff;
        
        for (i=_n-t3;i<_n;i++)
        {
            tt[i] = tint[_n-i];
            ti[tt[i]] = i;
        }
        
        
    }
#ifdef _CHECKINV
    for (i=0;i<_n;i++)
    {
        abkfatal(tt[ti[i]]==i,"Inverse doesn't check!\_n");
    }
#endif
    
}
bool ScanTourDZ::_scoFind2Opt(unsigned t1 )
{
    unsigned i;
    bool flag;
    unsigned t2, t3, t4;
    unsigned tt1, tt2, tt3, tt4;
    unsigned bestt3;
    double g1, gain, bestg;
    vector<vector<int> > const &nnd = _kdTree.getDistMatrix();
    vector<vector<unsigned> > const &nni = _kdTree.getIndexMatrix();
    unsigned nnear = _kdTree.getNNear();
    abkfatal(_pathValid,"Invalid tour!");
    abkfatal(_inverseValid,"Invalid inverse!");
    vector<unsigned>::iterator tt = _path.begin();
    vector<unsigned>::iterator ti = _pathinv.begin();
#ifdef COUNT
    nFindTwoOptCount++;
#endif
    
    tt1 = tt[t1];
    
    if (t1==0 && _optParams.bZeroFix)
        t1 = _n;
    
    t2 = t1-1;  // LOOK! this is a real problem in the !_optParams.bZeroFix
    //       situation.  How do we calculate gain?
    
    if (t2<0)    // LOOK! this is not necessary if _optParams.bZeroFix is true.
        t2=_n-1;  //       if _optParams.bZeroFix is false we have other problems.
    
    tt2 = tt[t2];
    
    /* 2-opt group rule: if one edge cuts a group, the other edge
    must have 1 end-point in the group */
    
    /* try replacing edge tt1-tt2 with each of the nnear edges from
    tt1 and from tt2, use the one that gives the best result */
    
    bestg = 0.0;
    g1 = realDist(tt2,tt1);
    
    /* first try all nnear edges from t1 to t4 */
    for (i=0; i<nnear; i++) 
    {
        tt4 = nni[tt1][i];
        t4 = ti[tt4];
        
        /*
        t3 = t4-1;
        
          if (t3 < 0)
          t3 = _n-1;
        */
        if (t4>0) t3=t4-1;
        else      t3=_n-1;
        
        tt3 = tt[t3];
        if (t3==t1 ||
            t3==t2-1  //LOOK!  What's wrong with having t3 just before t2?
            || (t3==0 && t1==_n)
            ) continue;        /* scho 3/96 */
        
        if (t1<t3)
            gain = g1 - nnd[tt1][i] + realDist(tt3,tt4) - realDist(tt2,tt3) + 
            _scoRevGain[t3] - _scoRevGain[t1];
        else
        { // wraparound case
            abkfatal(_optParams.bZeroFix,"need _optParams.bZeroFix set");
            gain = g1 - nnd[tt1][i] + realDist(tt3,tt4) - realDist(tt2,tt3) + 
                _scoRevGain[t3] + _scoRevGain[_n] - _scoRevGain[t1];
        }
        
        if (gain > bestg) 
        {
            if (!_bNoConstraints) //NCellPO > 0)
                flag = _legal2OptMove(t1, t3);
            else
                flag = true;
            
            if (flag) 
            {
                bestt3 = t3;
                bestg = gain;
            }
        }
    } /*  for (i=0;i<nnear;i++) */
    
    /* next try all nnear edges from t2 to t3 */
    for (i=0;i<nnear;i++) 
    {
        tt3 = nni[tt2][i];
        t3 = ti[tt3];
        
        if (t3==t1 ||
            t3==t2-1 || //LOOK!  what's wrong with this?
            (t3==0 && t1==_n)) continue;        /* scho 3/96 */
        
        t4 = t3+1;
        
        if (t4==_n)
            t4 = 0;
        
        tt4 = tt[t4];
        
        if (t1<t3)
            gain = g1 - nnd[tt2][i] + realDist(tt3,tt4) - realDist(tt1,tt4) + 
            _scoRevGain[t3] - _scoRevGain[t1];
        else
        { // wraparound case
            abkfatal(_optParams.bZeroFix,"need _optParams.bZeroFix set");
            gain = g1 - nnd[tt2][i] + realDist(tt3,tt4) - realDist(tt1,tt4) + 
                _scoRevGain[t3] + _scoRevGain[_n] - _scoRevGain[t1];
        }
        
        if (gain > bestg) 
        {
            if (!_bNoConstraints) //nCellPO > 0)
                flag = _legal2OptMove( t1, t3);
            else
                flag = true;
            
            if (flag) 
            {
                bestt3 = t3;
                bestg = gain;
            }
        }
    }  /*  for (i=0;i<tnnear[tt2];i++)  */
    
    if (bestg > FUZZ) 
    {
#ifdef COUNT
        nTwoOptCount++;
#endif
        
        _scoReversePath(t1,bestt3);
        //scoResetRevGain(tt,t1,bestt3);
        _scoInitRevGain();   //LOOK! REALLY! this is clearly a perf. loss
        //      wrt scoResetRevGain()
        
        if (_optParams.bZeroFix)
            _ComputeZeroConstraints();
        _Compute2OptLims();
        
        //printf("\_n *** scoReversePath in 2opt, t1 %d, bestt3 %d \_n\_n", t1, bestt3);
        //scoPrintCities(_n);
        
#ifdef _DEBUG
        abkfatal(_checkPO(),"Partial order failure in _scoFind2Opt");
        //{
        //fprintf(stderr,"\_n PO Fail in scoFind2Opt \_n");
        //scoPrintCities(_n); exit(1);
        //}
        abkassert(_isatour(),"Not a tour");
#endif
        
        return( true );
    } 
    else 
        return( false );
    
}
    
bool ScanTourDZ::_scoFind3Opt(unsigned t1 )
{
    abkfatal(_pathValid,"Invalid tour!");
    abkfatal(_inverseValid,"Invalid inverse!");
    vector<unsigned>::iterator tt = _path.begin();
    vector<unsigned>::iterator ti = _pathinv.begin();
    unsigned t2, t3, t4, t5, t6;
    unsigned tt1, tt2, tt3, tt4, tt5, tt6;
    unsigned greatest,least,middle;
    unsigned j, k;
    bool     flag;
    double g1, g2, g3;
    double d1, d2;
    vector<vector<int> > const &nnd = _kdTree.getDistMatrix();
    vector<vector<unsigned> > const &nni = _kdTree.getIndexMatrix();
    unsigned nnear = _kdTree.getNNear();
    vector<POType> const &CellPO=ScanCells::getPOs();
    unsigned NCellPO = CellPO.size();
#ifdef COUNT
    nFindThreeOptCount++;
#endif
    
    tt1 = tt[t1];
    
    if (t1==0)
        t1 = _n;
    
    t2 = t1-1;   tt2 = tt[t2];
    
    g1 = realDist(tt2,tt1);
    
    for (j=0; j < nnear; j++) 
    {
        
        if (nnd[tt2][j] > g1) 
            break;
        
        tt3 = nni[tt2][j];   t3 = ti[tt3];
        
        if (t3 == 0)
            t3 = _n;
        
        if ((t3==t1) || (t3==t2-1)) 
            continue;  /* scho */
        
        t4 = t3-1;   tt4 = tt[t4];
        
        d1 = realDist(tt4,tt3);
        g2 = g1 - nnd[tt2][j] + d1;
        
        /* newc = scoTourCost( &tt[i], _n ); */
        for (k=0; k < nnear; k++) 
        {
            
            if (nnd[tt4][k] > g2) 
                break;
            
            tt5 = nni[tt4][k];
            t5 = ti[tt5];
            
            if (t5 == 0)
                t5 = _n;
            
            //if (t5 <= t3 && (t1>t4 || t5>t2))
            //    continue;   /* city t5 must occur after t3 */
            
            if (t1>t3)
            {
                greatest = t1;
                least    = t3;
                middle   = t5;
            }
            else if (t3>t5)
            {
                greatest = t3;
                least    = t5;
                middle   = t1;
            }
            else
            {
                greatest = t5;
                least    = t1;
                middle   = t3;
            }
            
            if (middle <= least || greatest <= middle)
                continue;
            
            t6 = t5-1;
            tt6 = tt[t6];
            
            if (NCellPO > 0)
                flag = _legal3OptMove(greatest, least, middle);
            else
                flag = true;
            
            if (flag) 
            {     /* Legal Move?  3/96 */
                d1 = realDist(tt6,tt5);   d2 = realDist(tt6,tt1);
                g3 = g2 - nnd[tt4][k] + d1 - d2;
                
                if (g3 > FUZZ) 
                {
#ifdef COUNT
                    nThreeOptCount++;
#endif
                    
                    _scoMakeChange(greatest,least,middle);
                    
                    abkfatal(_checkPO(),"PO violation");
                    abkassert(_isatour(),"Not a tour");
                    return true;
                    
                    //if (_checkPO())
                    //    return( true );
                    
                    //else 
                    //{
                    //scoPrintCities(nCities);
                    //fprintf(stderr,"EXIt in a 3Opt; Not a tour"); exit(1);
                    //}
                    
                    //if (isatour(tt, _n)) 
                    //    return( TRUE );
                    //else 
                    //    {
                    //    fprintf(stderr,"EXIt in a 3Opt; Not a tour"); exit(1);
                    //   }
                    
                    //return( TRUE );
                } /*  if (g3 > FUZZ) */
            } /* if (flag) */
        } /* for (k=0; k < nnear; k++) */ 
        
    } /* for (j=0; j < nnear; j++) */
    
    return(false);  /* no improving move found */
    
}


void ScanTourDZ::_ComputeZeroConstraints()
{
    abkfatal(_zero2optableValid,"Invalid zero2optable flags");
    _Z2OptRtLim = 0;
    _zeroConstraintsValid=true;
    
    while (_Z2OptRtLim<_n-1 && _bZero2Optable[_path[_Z2OptRtLim+1]] )
        _Z2OptRtLim++;
    
    if (_Z2OptRtLim == _n-1)
    {
        _Z2OptLtLim = 1; // no constraints
        return;
    }
    
    else
    {
        _Z2OptLtLim = _n;
        
        while (_bZero2Optable[_path[_Z2OptLtLim-1]])
            _Z2OptLtLim--;
        
    }
    
    return;
}

void ScanTourDZ::_Compute2OptLims()
{
    abkfatal(_inverseValid,"Inverse invalid");
    unsigned i,j;
    unsigned From,To,PosnTo,PosnFrom;
    unsigned NCellPO = getPOs().size();
    
    // Set least restrictive possible initial values
    for (i=0;i<_n;i++)
        _RightLim2Opt[i] = _n-1;
    
    for (i=0;i<NCellPO;i++)
    {
        From=getPOs()[i].from;
        To=getPOs()[i].to;
        
        // If _optParams.bZeroFix is TRUE then PO constraints from and to "Cell 0"
        // are taken care of automatically by the fact that it's the
        // input and output cell
        
        if ((From && To) || !_optParams.bZeroFix)
        {
            PosnTo = _pathinv[To];     
            PosnFrom = _pathinv[From]; 
            
            for (j=PosnFrom;j!=unsigned(-1);j--)
            {
                if (_RightLim2Opt[j] < PosnTo)
                    break;  // RightLim2Opt[j] is nondecreasing in j,
                // so we can leave the loop now
                else
                    _RightLim2Opt[j] = PosnTo-1;
            }
            
        }
    }
    _2optLimsValid=true;
    return;
}

void ScanTourDZ::_ComputeZero2Optables()
{
    unsigned i;
    int From,To;
    unsigned NCellPO = getPOs().size();
    
    _bZero2Optable[0]=false;

    for (i=1;i<_n;i++)
        _bZero2Optable[i] = true;
    
    //_bZero2Optable[_n] = false; //sentinel value
    
    for (i=0;i<NCellPO;i++)
    {
        if ((From=getPOs()[i].from) && (To=getPOs()[i].to))
        {
            _bZero2Optable[From] = false;
            _bZero2Optable[To]   = false;
        }
    }
    _zero2optableValid=true;
}

void ScanTourDZ::_ComputeConstraintExistence()
{
    unsigned i;
    int From,To;
    unsigned NCellPO = getPOs().size();
    
    _bNoConstraints = true;
    
    for (i=0;i<NCellPO;i++)
    {
        From=getPOs()[i].from;
        To=getPOs()[i].to;
        
        if ((From && To)|| !_optParams.bZeroFix)
        {
            _bNoConstraints = false;
            return;
        }
    }
}

bool ScanTourDZ::_checkPO()
{
    int from, to, floc, tloc, p;
    unsigned NCellPO = getPOs().size();
    abkfatal(_inverseValid,"Invalid inverse tour");
    vector<unsigned>::iterator ti = _pathinv.begin();
    
    if (_bNoConstraints)
        return true;
    
    /* Check Direct PO relation based on PS & PE */
    for (p=NCellPO-1; p>=0; p--) 
    {
        from = getPOs()[p].from;
        to = getPOs()[p].to;
        floc = ti[from];
        tloc = ti[to];
        if ((floc % _n) > (tloc % _n))
        {
            // In the bZeroFix case, we can ignore PO constraints
            // to and from the zero cell.  However ones *from* the
            // zero cell will be satisfied anyway, since it has index
            // zero.
            
            if (!_optParams.bZeroFix || (to != 0) )
                return false;
        }
    }  /* for */
    
    return true;
}

#ifdef ABKDEBUG
//The reason we make this compilation conditional on ABKDEBUG
//is that there's a weird optimizer bug in gcc2.95.2 and 2.95.3
//which prevents it from compiling.  Therefore all refs to
//it must be abkassert, which is also conditional on ABKDEBUG
bool ScanTourDZ::_isatour()
{
    unsigned i;
    vector<unsigned> tt(_path);
    
    sort(tt.begin(),tt.end());
  
    for (i=0;i<tt.size();i++)
    {
        if (tt[i]!=i) 
        {
            return( false );
        }
    }
 
    return( true );
    
}
#endif

bool ScanTourDZ::_legal2OptMove(unsigned t1,unsigned t3)
{
    if (_optParams.bZeroFix)
    {
        if (t1 < t3)
            return (t3 <= _RightLim2Opt[t1]);
        else //wraparound case
            return ((t1 >= _Z2OptLtLim) && (t3 <= _Z2OptRtLim) );
    }
    else
        return ( (t1<t3) && (t3 <= _RightLim2Opt[t1]) ); 
}

double ScanTourDZ::scoTourCost() const
{
    unsigned i;
    double cost;
    
    cost = realDist(_path[_n-1],_path[0]);
    
    for (i=1; i<_n; i++) 
        cost += realDist(_path[i-1],_path[i]);
    
    return cost;
    
}

bool ScanTourDZ::_scoRun2Opt(double oldc )
{
    /* run a 2opt heuristic (i.e., try reversing subpaths in the tour */
#ifdef COUNT
    nRunTwoOptCount++;
#endif
    
    unsigned i;
    bool tchanged, changed;
    double newc;
    
    _scoInitRevGain();  /* initialize the store cost for reversing any
    subpath */
    
    if (_optParams.bZeroFix)
        _ComputeZeroConstraints();
    _Compute2OptLims();
    
    _ZeroShift = 0;
    changed = false;
    for (i=0; i<_n; i++) 
    {
        tchanged = _scoFind2Opt(i+_ZeroShift );
        if (tchanged) 
        { 
            newc = scoTourCost();
            abkfatal(newc <= oldc,"Increasing cost!");
            //if (newc > oldc)
            //    {
            //    fprintf(stderr,"increasing cost: %lf  to  %lf\_n",oldc,newc);
            //    exit (-507);
            //    }
            oldc = newc;
#ifdef DEBUG
            printf("new cost=%lf  t1=%d\_n",oldc,i);
#endif
            changed = true;
        }
        //tt[i+_n] = tt[i];
        //ti[tt[i]] = i+_n;
        //scoRevGain[i+_n] = scoRevGain[i+_n-1] + realDist(tt[i+_n-1],tt[i]) 
        //                    - realDist(tt[i],tt[i+_n-1]);
    }
    
    /* move tour back to start at 0 */
    //for (j=0; j<i; j++) 
    //    {
    //    tt[j] = tt[j+_n];
    //     ti[tt[j]] = j;
    //     }
    
    return( changed );
    
}

bool ScanTourDZ::_scoRun3Opt(double oldc )
{
    unsigned i;
    bool tchanged, changed;
    double newc;
    
#ifdef COUNT
    nRunThreeOptCount++;
#endif
    
    changed = false;
    for (i=0; i<_n; i++) 
    {
        tchanged = _scoFind3Opt( i );
        
        if (tchanged) 
        { 
            changed = true;
            /* newc = scoTourCost( &tt[i], n ); */
            newc = scoTourCost(  );
            abkfatal(newc<=oldc,"Increasing cost in 3-opt");
            
            //if (newc > oldc) 
            //    {
            //    fprintf(stderr,"increasing cost: %lf  to  %lf\n",oldc,newc);
            //    fprintf(stderr,"new cost=%lf  t1=%d\n",oldc,i);
            //    exit(-508);
            //    }
#ifdef DEBUG
            printf("new cost=%lf  t1=%d\n",oldc,i);
#endif
            oldc = newc;
        }
        //tt[i+n] = tt[i];
        //ti[tt[i]] = i+n;
    }
    
    /* move tour back to start at 0 */
    //for (i=0; i<n; i++) 
    //    {
    //    tt[i] = tt[i+n];
    //    ti[tt[i]] = i;
    //    }
    
    return(  changed );
    
}

/* Can [t3,t6] & [t5,t2} be swapped? */
bool ScanTourDZ::_legal3OptMove(int t1,int t3,int t5) 
{
    unsigned NCellPO = getPOs().size();
    abkfatal(_inverseValid,"Invalid inverse tour");
    vector<unsigned>::iterator ti = _pathinv.begin();
    int t2, t6, from, to, floc, tloc, p;
    
    t2 = t1 - 1;
    t6 = t5 - 1;
    
    for (p=NCellPO-1; p>=0; p--) 
    {
        from = getPOs()[p].from;
        to = getPOs()[p].to;
        
        // In the bZeroFix case, we don't worry about PO constraints
        // to or from the 0 cell.
        
        if ((!from || !to) && _optParams.bZeroFix)
            continue;
        
        floc = ti[from];
        tloc = ti[to];
        
        if ((t3 <= floc && floc <= t6) && (t5 <= tloc && tloc <= t2))
            return (false);
        
        if ((t3 <= tloc && tloc <= t6) && (t5 <= floc && floc <= t2))
            return (false);
        
        //printf("\_n legal3 (from %d to %d) (t1 %d t3 %d %t5 %d) floc %d,"
        //       " tloc %d",
        //        from, to, t1,t3, t3, floc, tloc);
    }  /* for */
    return (true);
}

void ScanTourDZ::_scoMakeChange(int t1,int t3,int t5)
{  /* make the 3-opt move on tour t (updating tour inverse array ti) */
    unsigned t2, t6;
    unsigned i,n1,n2;
    abkfatal(_pathValid,"Invalid tour");
    abkfatal(_inverseValid,"Invalid inverse tour");
    vector<unsigned>::iterator tt = _path.begin();
    vector<unsigned>::iterator ti = _pathinv.begin();
    
    t2 = t1 - 1;
    t6 = t5 - 1;
    
    n1 = t5-t3;
    n2 = t2+1-t5;
    
    /* figure out which of the two swapped subtours is longer; save
    that one to array temp.  Was not correct. 3/96 */
    if ( n2 > n1 ) {
        for (i=0; i<n2; i++)
            _temp[i] = tt[t5+i];
        for (i=0; i<n1; i++)
            tt[t2-i] = tt[t6-i];
        for (i=0; i<n2; i++)
            tt[t3+i] = _temp[i];
    }
    else {
        for (i=0; i<n1; i++)
            _temp[i] = tt[t3+i];
        for (i=0; i<n2; i++)
            tt[t3+i] = tt[t5+i];
        for (i=0; i<n1; i++)
            tt[t3+n2+i] = _temp[i];
    }
    /* update the inverse tour ti */
    for (i=t3; i<=t2; i++) {
        ti[tt[i]] = i;
        //if (i >= _n) 
        //    tt[i-_n] = tt[i];
    }
    
} /* scoMakeChange() */

void ScanTourDZ::_scoRunOpts()
{
    double oldc;
    bool changed, tchanged;
    
    oldc = scoTourCost();
#ifdef DEBUG
    printf("init cost=%lf\n",oldc);
#endif
    
    do 
    { /* while (changed) */
        if (!_optParams.only2Opt) 
        {
#ifdef DEBUG
            printf("    starting 3-opt\n");
#endif
            changed = _scoRun3Opt(oldc );
            oldc = scoTourCost( );
        }
        else 
            changed = false;
#ifdef DEBUG
        printf("    starting 2-opt\n");
#endif
        
        tchanged = _scoRun2Opt(oldc );
        oldc = scoTourCost();
        
        changed = changed || tchanged;
        
    } while (changed);
    
}

double ScanTourDZ::scoOptAll()
{
    unsigned j;
    double cost, lastCost, bsfCost, dCost,  tRand, tVal;
    int num_cons,  num_no_improvement=0;
    bool in_temp_control;
    long TriesCount;
    //int ti[MAX_CITIES];
    //int newTi[MAX_CITIES];
    vector<unsigned> newTi;
    //int bsfTour[MAX_CITIES];
    vector<unsigned> bsfTour;
    vector<unsigned> savedTour;
    
    /* initialize the t-inverse array again */
    //for (i=0;i<ncities;i++)
    //    ti[t[i]] = i;
    computeInverse();
    
    _ComputeConstraintExistence();
    
    if (_optParams.bZeroFix)
        _ComputeZero2Optables();
    
    //#ifdef _DEBUG
    //if (!isatour(t,ncities))
    //exit (-237);
    //#endif
    abkassert(_isatour(),"Not a tour");
    
    _scoRunOpts();
    bsfCost = scoTourCost();
    //copyTour(t,bsfTour,ti,ncities); //09/17/96 mro
    //_path.copyTo(bsfTour);
    bsfTour=_path;
    lastCost = bsfCost;
    
    
    in_temp_control = false;
    
    TriesCount = 0L;
    for (j=2; j<=_optParams.nDescents; j++) 
    {
        //_path.copyTo(savedTour);
        savedTour=_path;
        do 
        {
            switch (_optParams.kickMove) 
            {
            case 15:
                //_path.copyFrom(savedTour);
                _path=savedTour;_inverseValid=false;
                computeInverse();
                _randCorrect2Opt( );
                _randCorrect2Opt( );
                _randCorrect2Opt( );
                _randCorrect2Opt( );
                _randCorrect2Opt( );
                break;
                
            case 14:
                //_path.copyFrom(savedTour);
                _path=savedTour;_inverseValid=false;
                computeInverse();
                _randCorrect2Opt( );
                _randCorrect2Opt( );
                _randCorrect2Opt( );
                _randCorrect2Opt( );
                break;
                
            case 13:
                //_path.copyFrom(savedTour);
                _path=savedTour;_inverseValid=false;
                computeInverse();
                _randCorrect2Opt( );
                _randCorrect2Opt( );
                _randCorrect2Opt( );
                break;
                
            case 12:
                //_path.copyFrom(savedTour);
                _path=savedTour;
                computeInverse();_inverseValid=false;
                _randCorrect2Opt( );
                _randCorrect2Opt( );
                break;
                
            case 11:
                //_path.copyFrom(savedTour);
                _path=savedTour;_inverseValid=false;
                computeInverse();
                _randCorrect2Opt( );
                break;
                
            case 10:
                _random10Opt(savedTour);
                break;
            case 6:
                _random6Opt(savedTour);
                break;
            case 5:
                _random5Opt(savedTour);
                break;
            case 4:
                _random4Opt(savedTour);
                break;
            case 3:
                _random3Opt(savedTour);
                break;
            default:
                abkfatal3(false,
                    "Kick move ",_optParams.kickMove," is not implemented");
            } /* switch (kickMove) */
            
            if (++TriesCount > MAX_TRIES_FINDING_KICKMOVE)
                break;
        } while (!_checkPO());
        
        //#ifdef _DEBUG
        //        if (!isatour(newTour,ncities)) 
        //            exit (-238);
        //#endif
        abkassert(_isatour(),"Not a tour");
        
        _scoRunOpts ();
        cost = scoTourCost();
        
        if ((_optParams.temp_control) && 
            (!in_temp_control && num_no_improvement >= 3)) 
        {
            in_temp_control = true;
            _optParams.temperature = cost/_optParams.t_div;
            num_cons = 3;
        }
        
        if (in_temp_control) 
        {
            --num_cons;
            if (num_cons < 0) 
            {
                num_cons = 0;
                _optParams.temperature = 0.0;
                num_no_improvement = 0;
                in_temp_control = false;
            }
        }
        
        dCost = cost - lastCost;
        
        if (dCost < 0.0) 
        {
            in_temp_control = false;
            num_no_improvement = 0;
            _optParams.temperature = 0.0;
            //copyTour(newTour, t,ti, ncities);
            lastCost = cost;
        }
        else 
        {
            num_no_improvement++;
            if (_optParams.temperature > 0.0) 
            {
                tRand = double(_randuns.operator unsigned())/UINT_MAX;
                //UNIFORM_RAND;
                tVal = exp( -dCost / _optParams.temperature );
                if (tRand < tVal) 
                {
                    //copyTour(newTour, t,ti, ncities);
                    lastCost = cost;
                    printf("dCost %f temp %f tRand %f tVal %f\n",
                        dCost, _optParams.temperature, tRand, tVal);
                }
                else
                    _path=savedTour;_inverseValid=false;
            }
            else
                _path=savedTour;_inverseValid=false;
        }
        
        if (cost < bsfCost) 
        {
            num_no_improvement = 0;
            in_temp_control = false;
            num_cons = 0;
            _optParams.temperature = 0.0;
            bsfCost = cost;
            //copyTour(newTour,bsfTour,ti,ncities); //mro 9/17/96
            //_path.copyTo(bsfTour);
            bsfTour=_path;
        }
        
        
        if (TriesCount > MAX_TRIES_FINDING_KICKMOVE)
        {
            fprintf(stdout,"Timeout in scoOptAll()\n");
            break;
        }
        
        if (TriesCount > WARN_TRIES_FINDING_KICKMOVE)
        {
            fprintf(stdout,"Near timeout in scoOptAll()\n");
        }
    }  /* for (j=2; j<=_optParams.nDescents; j++) */
    
    //copyTour(bsfTour,t,ti,ncities);  //mro 9/17/96
    //_path.copyFrom(bsfTour);
    _path=bsfTour;_inverseValid=false;
    computeInverse();
    return bsfCost;            //mro 9/6/96
}



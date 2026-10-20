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




// randOpts.cxx:  Source file for those methods of ScanTourDZ which
//                generate a random n-opt move.

#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include "scanTourDZ.h"
using namespace abkscanopt;

const unsigned ScanTourDZ::maxTriesFindingRandopt=1000;

void ScanTourDZ::_randCorrect2Opt()
{
    bool found;
    unsigned start,end;
    int MaxDif;
    unsigned triesCount=0;
    
    //for (i=0;i<_n;i++)
    //    ti[tt[i]] = i;
    
    this->computeInverse();
    
    if (_optParams.bZeroFix)
        _ComputeZeroConstraints();
    _Compute2OptLims();
    
    found = false;
    
    while (!found)
    {
        if (++triesCount>maxTriesFindingRandopt)
        {
            abkfatal(false,"Timeout in ScanTourDZ::_randCorrect2Opt()");
        }
        
        start = _randuns%_n;
        //(int)(drand48()*_n);
        
        if (_bNoConstraints)
        {
            MaxDif = (_optParams.bZeroFix) ? (_n-1) : (_n-1-start);
        }
        
        else
        {
            if (_optParams.bZeroFix)
            {
                if (start == 0)
                {
                    start = _n;
                }
                
                if (start < _Z2OptLtLim)
                {
                    MaxDif = _RightLim2Opt[start] - start;
                }
                else
                {
                    MaxDif = _Z2OptRtLim - start + _n;
                }
                
            }
            else
            {
                MaxDif = _RightLim2Opt[start] - start;
            }
        }
        
        if (MaxDif<=1)
            continue;
        else
            found = true;
    }
    
    end = (start + _randuns%(MaxDif-1) + 1) % _n;
    //end = (start + (int)(drand48()*(MaxDif-1)) + 1) % _n;
    _scoReversePath(start,end);
    _scoInitRevGain();  
    
    //#ifdef _DEBUG
    //if (!checkPO(ti,_n))
    //    {
    //    fprintf(stderr,"PO violation in correct random 2-opts\_n");
    //    exit(-921);
    //    }
    //#endif
    abkassert(_checkPO(),"PO violation in correct random 2-opts");
    
}
void
ScanTourDZ::_random3Opt(vector<unsigned> &oldPath)
{
    vector<unsigned> perm(_n,UINT_MAX);
    unsigned temp;
    unsigned i, j;
    unsigned t1, t2, t3, t4, t5, t6;
    
    /* first choose random three edges to cut */
    for (i=0; i<_n; i++)
        perm[i] = i;
    for (i=0; i<2; i++) {
        j = i + _randuns%(_n-i);
        // (lrand48() % (_n-i));
        temp = perm[i]; perm[i] = perm[j]; perm[j] = perm[i];
    }
    /* we're interested in first 3 elements in perm, but make sure
    they're in order (bubble sort) */
    for (j=2; j>0; j--)
        for (i=0; i<j; i++)
            if (perm[i]>perm[i+1]) {
                temp = perm[i];  perm[i] = perm[i+1];  perm[i+1] = temp;
            }
            
            t1 = perm[0];  t2 = t1+1;
            t3 = perm[1];  t4 = t3+1;
            t5 = perm[2];  t6 = (t5+1)%_n;
            
            for (i=0; i<=t1; i++)
                _path[i] = oldPath[i];
            j = t1+1;
            for (i=t4; i<=t5; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            for (i=t2; i<=t3; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            if (t6>0)
                for (i=t6; i<_n; i++) {
                    _path[j] = oldPath[i];
                    j++;
                }
                
                
                computeInverse();
                
                abkfatal2(j==_n,"wrong # nodes in 3-opt tour = ",j);
                
} /* random3Opt() */
/*************************************************************/
void
ScanTourDZ::_random4Opt(vector<unsigned> &oldPath)		/* double bridge */
{
    vector<unsigned> perm(_n,UINT_MAX);
    unsigned temp;
    unsigned i, j;
    unsigned t1, t2, t3, t4, t5, t6, t7, t8;
    
    /* first choose random four edges to cut */
    for (i=0; i<_n; i++)
        perm[i] = i;
    for (i=0; i<4; i++) {
        j = i + _randuns%(_n-i);
        //(lrand48() % (_n-i));
        temp = perm[i]; perm[i] = perm[j]; perm[j] = perm[i];
    }
    
    /* we're interested in first 4 elements in perm, but make sure
    they're in order (bubble sort) */
    for (j=3; j>0; j--)
        for (i=0; i<j; i++)
            if (perm[i]>perm[i+1]) {
                temp = perm[i];  perm[i] = perm[i+1];  perm[i+1] = temp;
            }
            
            t1 = perm[0];  t2 = t1+1;
            t3 = perm[1];  t4 = t3+1;
            t5 = perm[2];  t6 = t5+1;
            t7 = perm[3];  t8 = (t7+1)%_n;
            
            for (i=0; i<=t1; i++)
                _path[i] = oldPath[i];
            j = t1+1;
            for (i=t6; i<=t7; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            for (i=t4; i<=t5; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            for (i=t2; i<=t3; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            if (t8>0) 
                for (i=t8; i<_n; i++) {
                    _path[j] = oldPath[i];
                    j++;
                }
                
                
                computeInverse();
                
                abkfatal2(j==_n,"wrong # nodes in 4-opt tour = ",j);
                
} /* random4Opt() */
/*************************************************************/
void
ScanTourDZ::_random5Opt(vector<unsigned> &oldPath)
{
    vector<unsigned> perm(_n,UINT_MAX);
    unsigned temp;
    unsigned i, j;
    unsigned t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    
    /* first choose random five edges to cut */
    for (i=0; i<_n; i++)
        perm[i] = i;
    for (i=0; i<5; i++) {
        j = i + _randuns%(_n-i);
        //(lrand48() % (_n-i));
        temp = perm[i]; perm[i] = perm[j]; perm[j] = perm[i];
    }
    
    /* we're interested in first 5 elements in perm, but make sure
    they're in order (bubble sort) */
    for (j=4; j>0; j--)
        for (i=0; i<j; i++)
            if (perm[i]>perm[i+1]) {
                temp = perm[i];  perm[i] = perm[i+1];  perm[i+1] = temp;
            }
            
            t1 = perm[0];  t2 = t1+1;
            t3 = perm[1];  t4 = t3+1;
            t5 = perm[2];  t6 = t5+1;
            t7 = perm[3];  t8 = t7+1;
            t9 = perm[4];  t10 = (t9 + 1)%_n;
            
            for (i=0; i<=t1; i++)
                _path[i] = oldPath[i];
            j = t1+1;
            for (i=t6; i<=t7; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            for (i=t2; i<=t3; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            for (i=t8; i<=t9; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            for (i=t4; i<=t5; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            if (t10>0) 
                for (i=t10; i<_n; i++) {
                    _path[j] = oldPath[i];
                    j++;
                }
                
                computeInverse();
                
                abkfatal2(j==_n,"wrong # nodes in 5-opt tour = ",j);
                
} /* random5Opt() */
/*************************************************************/
void
ScanTourDZ::_random6Opt(vector<unsigned> &oldPath)
{
    vector<unsigned> perm(_n,UINT_MAX);
    unsigned temp;
    unsigned i, j;
    unsigned t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
    
    /* first choose random six edges to cut */
    for (i=0; i<_n; i++)
        perm[i] = i;
    for (i=0; i<6; i++) {
        j = i + _randuns%(_n-i);
        //(lrand48() % (_n-i));
        temp = perm[i]; perm[i] = perm[j]; perm[j] = perm[i];
    }
    
    /* we're interested in first 6 elements in perm, but make sure
    they're in order (bubble sort) */
    for (j=5; j>0; j--)
        for (i=0; i<j; i++)
            if (perm[i]>perm[i+1]) {
                temp = perm[i];  perm[i] = perm[i+1];  perm[i+1] = temp;
            }
            
            t1 = perm[0];  t2 = t1+1;
            t3 = perm[1];  t4 = t3+1;
            t5 = perm[2];  t6 = t5+1;
            t7 = perm[3];  t8 = t7+1;
            t9 = perm[4];  t10 = t9 + 1;
            t11 = perm[5]; t12 = (t11 + 1)%_n; 
            
            for (i=0; i<=t1; i++)
                _path[i] = oldPath[i];
            j = t1+1;
            for (i=t8; i<=t9; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            for (i=t4; i<=t5; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            for (i=t10; i<=t11; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            for (i=t2; i<=t3; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            
            for (i=t6; i<=t7; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            if (t12>0)
                for (i=t12; i<_n; i++) {
                    _path[j] = oldPath[i];
                    j++;
                }
                
                
                computeInverse();
                
                abkfatal2(j==_n,"wrong # nodes in 3-opt tour = ",j);
                
} /* random6Opt() */
/*************************************************************/
void
ScanTourDZ::_random10Opt(vector<unsigned> &oldPath)
{
    vector<unsigned> perm(_n,UINT_MAX);
    unsigned temp;
    unsigned i, j;
    unsigned t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20;
    
    /* first choose random ten edges to cut */
    for (i=0; i<_n; i++)
        perm[i] = i;
    for (i=0; i<10; i++) {
        j = i + _randuns%(_n-i);
        //(lrand48() % (_n-i));
        temp = perm[i]; perm[i] = perm[j]; perm[j] = perm[i];
    }
    
    /* we're interested in first 10 elements in perm, but make sure
    they're in order (bubble sort) */
    for (j=9; j>0; j--)
        for (i=0; i<j; i++)
            if (perm[i]>perm[i+1]) {
                temp = perm[i];  perm[i] = perm[i+1];  perm[i+1] = temp;
            }
            
            t1 = perm[0];  t2 = t1+1;
            t3 = perm[1];  t4 = t3+1;
            t5 = perm[2];  t6 = t5+1;
            t7 = perm[3];  t8 = t7+1;
            t9 = perm[4];  t10 = t9 + 1;
            t11 = perm[5]; t12 = t11 + 1;
            t13 = perm[6]; t14 = t13 + 1;
            t15 = perm[7]; t16 = t15 + 1;
            t17 = perm[8]; t18 = t17 + 1;
            t19 = perm[9]; t20 = (t19 + 1)%_n;
            
            for (i=0; i<=t1; i++)
                _path[i] = oldPath[i];
            j = t1+1;
            for (i=t10; i<=t11; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            for (i=t6; i<=t7; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            for (i=t14; i<=t15; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            for (i=t2; i<=t3; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            
            for (i=t18; i<=t19; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            for (i=t8; i<=t9; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            for (i=t16; i<=t17; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            for (i=t4; i<=t5; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            for (i=t12; i<=t13; i++) {
                _path[j] = oldPath[i];
                j++;
            }
            if (t20>0)
                for (i=t20; i<_n; i++) {
                    _path[j] = oldPath[i];
                    j++;
                }
                
                
                computeInverse();
                
                abkfatal2(j==_n,"wrong # nodes in 3-opt tour = ",j);
                
} /* random10Opt() */

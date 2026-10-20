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
// scanTourDZ.h: interface for the ScanTourDZ class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_SCANTOURDZ_H__INCLUDED_)
#define _SCANTOURDZ_H__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <ScanChain/scanChain.h>
#include "kenKDTree.h"
namespace abkscanopt
{

#define MAX_TRIES_FINDING_KICKMOVE 1000000
#define WARN_TRIES_FINDING_KICKMOVE 10000

class KenKDTree;

// The "DZ" stands for "distinguished zero".  It's a "tour"
// rather than a "chain" because it represents a closed
// loop -- the zero cell is the input *and* output.

class ScanTourDZ : public ScanChain  
{
public:
    const static unsigned maxTriesFindingRandopt;
    class OptParams
    {
    public:
        bool     only2Opt; //don't use 3-opts in greedy descent
        unsigned nDescents;
        unsigned kickMove;
        bool     temp_control;
        double   temperature;
        double   t_div;
        bool     bZeroFix;
        OptParams():
        only2Opt(false),
            nDescents(5),
            kickMove(4),
            temp_control(false),
            temperature(0.0),
            t_div(100.0),
            bZeroFix(true){}
    };
    protected:
        const KenKDTree _kdTree;
        RandomRawUnsigned &_randuns;
        unsigned _n; //number of scan cells
        
        OptParams _optParams;
        
        //**************************
        // This next block of variables deals with the capability of treating
        // cell "0" (as an index to scanCell[]) as a "special" cell representing
        // both the input and the output of the path we're searching for, so that
        // we look for a tour and the length of the tour is the same as the
        // length of the corresponding path.
        //
        // When we do this, we want to keep cell "0" fixed to facilitate the special
        // calculations involving the partial order constraints.  We keep it
        // fixed in position "0" of the tour[] array.
        //
        //**************************
        
        //bool _bZeroFix;   // treat the tour in the manner described above.
        // moved to OptParams
        
        bool _bNoConstraints;    // If this is TRUE it means that there are
        // no partial order constraints (except, in
        // the bZeroFix case, those going to or
        // from "0").
        
        int  _ZeroShift;     // In the loop in ScoRun2Opt(), we repeatedly
        // call ScoFind2Opt().  If we do a 2-opt that
        // moves the position of cell "0", we have
        // to move it back.  This variable keeps track
        // of the cumulative effect of such moves.
        // NOTE: seems to be unused at the moment
        
        vector<bool> _bZero2Optable; // This will be TRUE for a city
        // just in case there are no PO constraints
        // involving the city, except possibly]
        // those going to or from city "0"
        // A 2-opt move that reverses a subpath
        // containing city "0" is possible just
        // in case this value is TRUE for all cities
        // in the subpath.
        
        unsigned _Z2OptRtLim;                 // these are the positions of the rightmost and
        unsigned _Z2OptLtLim;                 // leftmost
        // cells that can be included in a 2-opt
        // reversing a subpath including "0"
        // Since we're considering wrapping around
        // "0", Z2OptLtLim will actually be greater
        // than Z2OptRtLim, except in the case
        // that all cities have bZero2Optable TRUE.
        // Note that these are positions, not
        // cell indices.
        
        bool _zeroConstraintsValid;
        bool _2optLimsValid;
        bool _zero2optableValid;
        
        //**************************
        // This next variable is for the 2-opt speedup.
        //**************************
        
        vector<unsigned> _RightLim2Opt;    // For a city in position i in the tour,
        // RightLim2Opt[i] is the rightmost j such
        // that you can reverse the subpath [i,j].
        // Note that both i and j are positions,
        // not city indices.  Note also that j
        // will never be less than i; subpaths that
        // wrap around through position 0 will be
        // considered separately.
        
        vector<double> _scoRevGain;
        
        mutable vector<unsigned> _temp; // scratch space for doing 3-opts
        
        void _scoInitRevGain();
        void _scoReversePath(unsigned t1,unsigned t3);
        bool _scoFind2Opt(unsigned t1 );
        bool _scoFind3Opt(unsigned t1 );
        bool _legal2OptMove(unsigned t1, unsigned t3);
        bool _checkPO();
        bool _isatour();
        
        void _ComputeZeroConstraints();
        void _Compute2OptLims();
        void _ComputeZero2Optables();
        void _ComputeConstraintExistence();
        bool _scoRun2Opt(double oldc );
        bool _scoRun3Opt(double oldc );
        bool _legal3OptMove(int t1,int t3,int t5);
        void _scoMakeChange(int t1,int t3,int t5);
                
        void _randCorrect2Opt();
        void _random3Opt(vector<unsigned> &oldTour);
        void _random4Opt(vector<unsigned> &oldTour);
        void _random5Opt(vector<unsigned> &oldTour);
        void _random6Opt(vector<unsigned> &oldTour);
        void _random10Opt(vector<unsigned> &oldTour);
        void _scoRunOpts();
        
        
    public:
        ScanTourDZ(ScanChain &origChain,
                     unsigned nnear,
                     RandomRawUnsigned &randuns,
                     OptParams optParams);
        
        double scoOptAll();
        double scoTourCost() const;
        
        virtual ~ScanTourDZ();
        
        
    };
};    
#endif // !defined(_SCANTOURDZ_H__INCLUDED_)

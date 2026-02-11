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
// optimizer1.h: interface for the Optimizer1 class.
//               This is the only existing concrete class derived from
//               Optimizer
//////////////////////////////////////////////////////////////////////

#if !defined(_ABKSCANOPT_OPTMIZER1_H)
#define _ABKSCANOPT_OPTMIZER1_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "optimizer.h"

namespace abkscanopt
{
class Optimizer1 : public Optimizer 
{
protected:
    Path       _lastPath; //path to restore if optimization failed to improve
                          //cost

    double     _optimizeLevel(); //return value is cost

public:
	Optimizer1(ScanChain &chain,RandomRawUnsigned &randuns,
               Params params);
	virtual ~Optimizer1();

    virtual void rawStatsOut(::ostream &os) const;
};
};

#endif // !defined(_ABKSCANOPT_OPTMIZER_H)

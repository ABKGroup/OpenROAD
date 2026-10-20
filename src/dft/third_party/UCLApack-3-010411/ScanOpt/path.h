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
// path.h: interface for the Path class.
//         This class stores possible solutions to a ScanChain
//         optimization problem, even if the path in subgroups
//         differs from the current one
//////////////////////////////////////////////////////////////////////

#if !defined(_SCANOPT_PATH_H)
#define _SCANOPT_PATH_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include <ScanChain/scanChain.h>

namespace abkscanopt
{

    class Cell;
//This class describes a path through the scan chain hierarchically,
//going into all the subgroups.
class Path  
{
    ScanChain &_chain;
    vector<unsigned> _path;
    vector<Path*> _subpaths;
    bool          _pathValid;
public:
	Path(ScanChain &chain);
	virtual ~Path();

    //save path of _chain in this Path object
    void save();

    //restore path from this Path object to _chain
    void restore();

    double cost() const;

    //gets the "real" cell in position i of *this*
    //path, i.e. with in/out locs determined by the
    //Path object rather than the current path in the chain.
    void getCellInLoc(unsigned i,int &inX,int &inY) const;
    void getCellOutLoc(unsigned i,int &outX,int &outY) const;

};
}; //end namespace abkscanopt

#endif // !defined(_SCANOPT_PATH_H)


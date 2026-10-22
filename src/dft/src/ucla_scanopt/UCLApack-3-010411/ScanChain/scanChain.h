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
// ScanChain.h: interface for the ScanChain class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_SCANCHAIN_H__INCLUDED_)
#define _SCANCHAIN_H__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <ABKCommon/abkcommon.h>
#include <vector>
//#include <Combi/permut.h>
#include "scanCells.h"

namespace abkscanopt
{

class ScanChain : public ScanCells  
{
public:
    class collapseInOut{}; //used to identify "collapsing ctor"
    friend class ScanCells; //needs to call _modifyPath()
    friend class Path;
    friend class Optimizer;
protected:
    
    // _path tells you the order in which the cells are
    // visited; _pathinv is its inverse.
    vector<unsigned>          _path;
    vector<unsigned>          _pathinv;
    bool              _inverseValid; // if false, then _pathinv meaningless
    bool              _pathValid; // if false, then _path meaningless

    
    void _pathReversePOViolations();

    //protected ctor.  To create a ScanChain object, derive
    //a class with a public ctor, from which you can call
    //ScanCells::_addCell(), ScanCells::_addSubgroups(),
    //ScanCells::_transAddPOElt().
    ScanChain():
        _path(0),_pathinv(0),_inverseValid(false),_pathValid(false){}

    //"collapsing ctor" used by ScanTourDZ
    ScanChain(ScanChain &orig,collapseInOut);
    
private:
    //Change the input and output of the path
    void _modifyPath(unsigned newIn,unsigned newOut);

public:
    
    virtual ~ScanChain();

    bool isPathValid() const {return _pathValid;}
    
    void computeInverse();
    
    const vector<unsigned> &getPath() const {return _path;}
    void printPathTopDown(ostream &os,unsigned nTabs) const;
    void printPathXgraph(ostream &os) const;

    //Choose an initial path at random
    void initialPath(RandomRawUnsigned &randuns);

    void checkPath() const;
    
    
};


}; // end namespace abkscanopt
#endif // !defined(_SCANCHAIN_H__INCLUDED_)

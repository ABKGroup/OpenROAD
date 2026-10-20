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
// path.cxx: implementation of the Path class.
//
//////////////////////////////////////////////////////////////////////

#include "path.h"
using namespace abkscanopt;
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Path::Path(ScanChain &chain)
:_chain(chain),_path(_chain.getPath()),_subpaths(_chain.getNumSubgroups()),
_pathValid(_chain.isPathValid())
{
    const unsigned numSubgroups=_chain.getNumSubgroups();
    unsigned i;
    for (i=0;i<numSubgroups;i++)
    {
        _subpaths[i] = new Path(*_chain.getCells()[i].subgroup);
    }

}
//***************************************************************************

Path::~Path()
{
    const unsigned numSubpaths=_subpaths.size();
    unsigned i;
    for (i=0;i<numSubpaths;i++)
        delete _subpaths[i];

}
//***************************************************************************

void Path::save()
{
    _pathValid=_chain.isPathValid();
    _path=_chain.getPath();
    const unsigned numSubpaths=_subpaths.size();
    unsigned i;
    for (i=0;i<numSubpaths;i++)
        _subpaths[i]->save();

}
//***************************************************************************

void Path::restore()
{
    _chain._path=_path;
    _chain._pathValid=_pathValid;
    _chain.computeInverse();
    const unsigned numSubpaths=_subpaths.size();
    unsigned i;
    for (i=0;i<numSubpaths;i++)
        _subpaths[i]->restore();

}

//***************************************************************************
double Path::cost() const
{
    abkwarn(_pathValid,"computing cost of invalid path");
    double retval=0;
    const unsigned numSubpaths=_subpaths.size();
    unsigned i;
    for (i=0;i<numSubpaths;i++)
        retval += _subpaths[i]->cost();
    const unsigned numCells=_path.size();


    for (i=0;i<numCells-1;i++)
    {
        int X1,X2,Y1,Y2;
        getCellOutLoc(i,X1,Y1);
        getCellInLoc(i+1,X2,Y2);
        retval +=ScanChain::realDistStatic(X1,Y1,X2,Y2);
    }
    return retval;
}

//***************************************************************************
void Path::getCellInLoc(unsigned i,int &inX,int &inY) const
{
    const unsigned cellIdx=_path[i];
    if (cellIdx>=_subpaths.size())
    {
        const Cell &cell=_chain.getCells()[cellIdx];
        inX=cell.inX;
        inY=cell.inY;
    }
    else
    {
        const Path &subpath=*_subpaths[cellIdx];
        subpath.getCellInLoc(0,inX,inY);
    }
}
//***************************************************************************
void Path::getCellOutLoc(unsigned i,int &outX,int &outY) const
{
    const unsigned cellIdx=_path[i];
    if (cellIdx>=_subpaths.size())
    {
        const Cell &cell=_chain.getCells()[cellIdx];
        outX=cell.outX;
        outY=cell.outY;
    }
    else
    {
        const Path &subpath=*_subpaths[cellIdx];
        subpath.getCellOutLoc(subpath._path.size()-1,outX,outY);
    }
}


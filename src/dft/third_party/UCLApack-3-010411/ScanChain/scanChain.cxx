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
// ScanChain.cxx: implementation of the ScanChain class.
//
//////////////////////////////////////////////////////////////////////
#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include "scanChain.h"
#include "compare.h"
#include <algorithm>
using namespace abkscanopt;
using std::copy;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//"collapsing ctor" used by ScanTourDZ
ScanChain::ScanChain(ScanChain &orig,collapseInOut)
:ScanCells(orig,orig.getPath()[0],orig.getPath()[orig.getNumCells()-1]),
_path(orig.getNumCells()-1,UINT_MAX)
{
    unsigned inIdx=orig.getPath()[0];
    unsigned outIdx=orig.getPath()[orig.getNumCells()-1];
    //This is a specialized ctor used only by ScanTourDZ
    unsigned lesser=inIdx,greater=outIdx;
    if (lesser>greater) std::swap(lesser,greater);
    unsigned i,iMax=orig.getNumCells()-1;
    const vector<unsigned> &opath=orig.getPath();
    _path[0]=0;
    for (i=0;i<iMax;i++)
    {
        _path[i] = ScanCells::collapseIndex(opath[i],lesser,greater);
    }
    abkfatal(_path[0]==0,"Internal error");
    _pathValid=true;
    computeInverse();

}

//***************************************************************************

void ScanChain::computeInverse()
{
    abkfatal(_pathValid,"Tried to invert invalid path");
    const unsigned iMax=_path.size();
    _pathinv.clear();
    _pathinv.insert(_pathinv.end(),iMax,UINT_MAX);
    unsigned i;
    for (i=0;i<iMax;i++)
        _pathinv[_path[i]]=i;
    _inverseValid=true;
}
//***************************************************************************

void ScanChain::printPathTopDown(ostream &os,unsigned nTabs) const
{
    abkfatal(_pathValid,"Attempt to print an invalid path");
    unsigned i,j;
    for (i=0;i<getNumCells();i++)
    {
        unsigned idx=_path[i];
        for (j=0;j<nTabs;j++) os << "    ";
        os << idx ;
        os << " (" <<  this->identifyIndex(idx).c_str() << ")" << endl;
        
        if (idx<getNumSubgroups())
        {
            const Cell &cell=getCells()[idx];
            abkfatal(cell.subgroup,"Indexing error!");
            cell.subgroup->printPathTopDown(os,nTabs+1);
        }
    }
    os << endl;
}
//***************************************************************************
void ScanChain::printPathXgraph(ostream &os) const
{
    unsigned i;
    const vector<Cell> &cells=getCells();

    for (i=0;i<getNumCells()-1;i++)
    {
        const Cell &from = cells[_path[i]];
        const Cell &to   = cells[_path[i+1]];
        os << "move " << from.outX << " " << from.outY << endl;
        os << "draw " << to.inX << " " << to.inY << endl;
    }

    os << "\"" << getName().c_str() << "\"" << endl << endl;

    for (i=0;i<getNumSubgroups();i++)
    {
        cells[i].subgroup->printPathXgraph(os);
    }
}
//***************************************************************************

ScanChain::~ScanChain()
{
    
}

//***************************************************************************

void ScanChain::initialPath(RandomRawUnsigned &randuns)
{
    abkfatal(!_pathValid,"Error:  called initialPath() "
        "when path already existed");
    //Adapted from Sun Cho's function randomPath()
    unsigned i;
    unsigned NewIn,NewOut;
    
    ScanCells::_getRandomPermutation(_path,getNumCells(),randuns);
    _pathValid=true;_inverseValid=false;

    computeInverse();
    
    
    const unsigned nLegalins=getLegalins().size();
    const unsigned nLegalouts=getLegalouts().size();
    do
    {
        NewIn = getLegalins()[randuns%nLegalins];
        NewOut = getLegalouts()[randuns%nLegalouts];
    } while  (NewIn==NewOut);
    
    _modifyPath(NewIn,NewOut);
    
    //m_nInCellIdx = m_pnPathIdxArray[0];
    //m_nOutCellIdx = m_pnPathIdxArray[m_nNumberSubGroups+m_nNumberLooseCells-1];
    
    // Now we get an initial path for each subgroup, and fill in the
    // "shadow cells" in m_ppcSubGroupsArray
    for (i=0;i<getNumSubgroups();i++)
    {
        getCells()[i].subgroup->initialPath(randuns);
    }
    
    checkPath();
}

//***************************************************************************

void ScanChain::_modifyPath(unsigned newIn,unsigned newOut)
{
    vector<unsigned> DoubledArray;
    int gap;
    int NewInPosn,NewOutPosn;
    
    NewInPosn = _pathinv[newIn];
    NewOutPosn = _pathinv[newOut];
    
    gap = NewOutPosn-NewInPosn;
    
    if (gap<0) gap += getNumCells();
    
    DoubledArray.reserve(_path.size()*2);
    // Double up array
    DoubledArray.insert(DoubledArray.end(),_path.begin(),_path.end());
    DoubledArray.insert(DoubledArray.end(),_path.begin(),_path.end());
    
    //Now copy into new order
    copy(DoubledArray.begin()+NewInPosn,DoubledArray.begin()+NewInPosn+gap,
        _path.begin());
    
    copy(DoubledArray.begin()+NewOutPosn+1,DoubledArray.begin()
        +NewOutPosn+getNumCells()-gap,_path.begin()+gap);
    
    _path[getNumCells()-1] = newOut;_inverseValid=false;
    computeInverse();
    
    
    //swap PO violations
    _pathReversePOViolations();

#ifdef ABKDEBUG
    checkPath(); //TODO: remove
#endif
    

}

//***************************************************************************
// The proper functioning of this next function depends *essentially* on
// the fact that the partial order elements are maintained in such a way
// that if you have A<B,B<C and A<C, the last one deducible by transitivity
// from the first two, then A<C always appears *after* the other two.
// Then we have to loop through the PO constraints *backwards*.
//***************************************************************************

void ScanChain::_pathReversePOViolations()
{
    unsigned from,to;
    unsigned i;
    
    //LOOK! it's not clear this is any faster than
    //      the brute-force sort below, given the
    //      inefficiency of the comparison operator
    if (isSubpath())
    {
        CompareBasedOnTotalOrder c(*this);
        std::sort(_path.begin(),_path.end(),c);_inverseValid=false;
        computeInverse();
    }
    
    //Now for the usual case
    
    for (i=getPOs().size()-1;i!=static_cast<unsigned>(-1);i--)
    {
        const POType &po=getPOs()[i];
        from = po.from;
        to   = po.to;
        if  (_pathinv[to]<_pathinv[from])
        {
            _path[_pathinv[to]] = from;
            _path[_pathinv[from]] = to;
            std::swap(_pathinv[to],_pathinv[from]);
        }
    }
}

//***************************************************************************

void ScanChain::checkPath() const
{
    unsigned i;
    unsigned from,to;
    
    for (i=0;i<getNumCells();i++)
        if (_path[_pathinv[i]] != i)
            abkfatal(false,"_path and _pathinv are not inverses");
        
    bool found=false;
    const unsigned nLegalins=getLegalins().size();
    const unsigned inIdx=_path[0];
    for (i=0;i<nLegalins;i++)
    {
        if (inIdx == getLegalins()[i])
        {
            found = true;
            break;
        }
    }
    
    abkfatal(found,"Illegal in cell in InitialPath()\n");
    
    found = false;
    
    const unsigned nLegalouts=getLegalouts().size();
    const unsigned outIdx=_path[getNumCells()-1];
    for (i=0;i<nLegalouts;i++)
    {
        if (outIdx == getLegalouts()[i])
        {
            found = true;
            break;
        }
    }
    
    abkfatal(found,"Illegal out cell in InitialPath()\n");
    
    const unsigned nPOs=getPOs().size();
    for (i=0;i<nPOs;i++)
    {
        const POType &po=getPOs()[i];
        from = po.from;
        to   = po.to;
        abkfatal(_pathinv[to]>=_pathinv[from],"PO violation");
    }
    
}

//***************************************************************************
//***************************************************************************

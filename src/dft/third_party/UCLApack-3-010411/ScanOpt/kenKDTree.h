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
// KenKDTree.h: interface for the KenKDTree class.
// k-d tree implementation adapted from work of Ken Boese
//////////////////////////////////////////////////////////////////////

#if !defined(_KENKDTREE_H____INCLUDED_)
#define _KENKDTREE_H____INCLUDED_

#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif


#include <vector>
#include <ScanChain/scanChain.h>
namespace abkscanopt
{

struct kdnode /* nodes for K-d trees */
{
    bool  bucket;
    int  cutdim;
    int cutval;
    kdnode* lowson;
    kdnode* hison;
    int lowpt, hipt;
    bool nodeEmpty;
    int size;
    kdnode* father;
};

typedef kdnode* kdIter;

const unsigned kdCutOff=5;


class ScanCells;

class KenKDTree  
{
    class sortByX
    {
    private:
        std::vector<Cell> const &_scanCell;
    public:
        sortByX(std::vector<Cell> const &scanCellPassed):
          _scanCell(scanCellPassed){}
          virtual bool operator()(const unsigned i1, const unsigned i2);
    };
    
    class sortByY
    {
    private:
        std::vector<Cell> const &_scanCell;
    public:
        sortByY(std::vector<Cell> const &scanCellPassed):
          _scanCell(scanCellPassed){}
          virtual bool operator()(const unsigned i1, const unsigned i2);
    };
    
    sortByX _sortX;
    sortByY _sortY;
    
    std::vector<Cell> const &_scanCell;
    std::vector<POType> const &_CellPO;
    std::vector<unsigned> const &_POEnd;
    
    unsigned _n,_nnear;

    //distances to nearest cells
    std::vector<std::vector<int> >        _nnd;

    //indices of nearest cells
    std::vector<std::vector<unsigned> >    _nni;
    
    std::vector<int> _D;
    std::vector<unsigned> _A;
    
    vector<unsigned> _perm;
    vector<unsigned> _yperm;
    std::vector<kdnode> _nodes;
    
    kdIter _kdMemIndex;
    
    unsigned _nntarget;
    int _nndist;
    unsigned _currentNum;
    unsigned _heapSize;
    
    kdIter _kdRoot;
    
    std::vector<kdIter> _bucketptr;
    int _optXIn(unsigned i,unsigned d);
    int _optXOut(unsigned i,unsigned d);
    
    int _findmaxspread(unsigned low,unsigned u);
    void _newSelect(unsigned low,unsigned u,unsigned m,unsigned cutdim );
    unsigned _runPartition(unsigned low,unsigned hi,unsigned m,
        unsigned cutdim );
    
    void _initHeap();
    void _treeBuild();
    kdIter _newNode();
    kdIter _kdTreeBuild(unsigned low,unsigned u,kdIter father);
    void _deleteNode(unsigned node);
    void _undeleteNode(unsigned node);
    void _nns(unsigned node,std::vector<unsigned> &neighbors);
    void _rnnMultiple(kdIter p);
    unsigned _extractTopHeap();
    void _insertHeap(unsigned index,int d );
    void _insertNode(kdIter p);
    void _heapify(unsigned i );
    
    
    public:
        KenKDTree(ScanCells const &scanCells,unsigned nnear);
        
        inline int realDist(unsigned i,unsigned j) const
        {return ScanCells::realDistStatic(_scanCell[i].outX,_scanCell[i].outY,
                               _scanCell[j].inX,_scanCell[j].inY);}
        
        unsigned getNNear() const {return _nnear;}
        
        std::vector<std::vector<int> > const &getDistMatrix() const
        {return _nnd;}
        std::vector<std::vector<unsigned> > const &getIndexMatrix() const
        {return _nni;}
        
};
};

#endif // _KENKDTREE_H____INCLUDED_

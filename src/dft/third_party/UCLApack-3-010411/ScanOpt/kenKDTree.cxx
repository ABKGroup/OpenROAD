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
// KenKDTree.cxx: implementation of the KenKDTree class.
// k-d tree implementation adapted from work of Ken Boese
//////////////////////////////////////////////////////////////////////

#include "kenKDTree.h"
#include <algorithm>
#include <float.h>
#define parent(i) ((i-1)>>1)
#define left(i) (((i+1)<<1)-1)
#define right(i) ((i+1)<<1)

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

using namespace abkscanopt;
using std::sort;
KenKDTree::KenKDTree(ScanCells const &scanCells,unsigned nnear):
_sortX(scanCells.getCells()),
_sortY(scanCells.getCells()),
_scanCell(scanCells.getCells()),
_CellPO(scanCells.getPOs()),
_POEnd(scanCells.getPOHeadTos()),
_n(scanCells.getNumCells()),_nnear(nnear),
_nnd(_n,vector<int>(_nnear)),
_nni(_n,vector<unsigned>(_nnear)),
_D(_n),_A(_n),
_perm(_n,UINT_MAX),
_yperm(_n,UINT_MAX),
_nodes(_n),
_bucketptr(_n)

{
    unsigned  i, j, k;
    for (i=0;i<_n;i++)
        _perm[i]=_yperm[i]=i;
    sort(_perm.begin(),_perm.end(),_sortX);
    unsigned t;
    bool flag;
    vector<unsigned> neighbors(_nnear,UINT_MAX);
    
    
    _treeBuild();
    
    for (i=0;i<_n;i++) 
    {
        _deleteNode( i );
        _nns( i, neighbors);
        
        for (j=0;j<_nnear;j++) 
        {
            k = neighbors[j];
            t = _POEnd[i];           /* check a PO constraint 3/96  */
            _nnd[i][j] = realDist(i,k);
            _nni[i][j] = k;
            flag = true;
            
            while (flag && (t!=UINT_MAX)) 
            {     /* check a PO constraint 3/96  */
                //if (GCellPO[nGroupNo][t].from == k) 
                if (_CellPO[t].from == k) 
                {
                    flag = false;
                    _nnd[i][j] = INT_MAX/2-1;
                } 
                else 
                    //t = GCellPO[nGroupNo][t].tnext;
                    t = _CellPO[t].tnext;
            }
        }   /* for (j=0;j<nnear;j++) */
        _undeleteNode( i );
    } /* for (i=0;i<n;i++) */
    
}

void KenKDTree::_initHeap()
{
    _heapSize = 0;
}

void KenKDTree::_nns(unsigned node,vector<unsigned> &neighbors)
{
    int i;
    
    _nntarget = node;
    _nndist = INT_MAX/2; //without the "/2", get critical  comparison error
    _initHeap();
    //numNear = k;
    _currentNum = 0;
    _rnnMultiple( _kdRoot );
    for (i=_nnear-1; i>=0; i--)
        neighbors[i] = _extractTopHeap();
    
}

kdIter KenKDTree::_newNode()
{
    return _kdMemIndex++;
}

void KenKDTree::_treeBuild()
{
    
    _kdMemIndex = _nodes.data();
    _kdRoot = _kdTreeBuild(0, _n-1, kdIter(NULL));
    
}

kdIter KenKDTree::_kdTreeBuild(unsigned low,unsigned u,kdIter father)
{
    kdIter p;
    unsigned i, m;
    
    /*p = (kdnode *) malloc( sizeof(kdnode) );*/
    p = _newNode();
    p->lowpt = low;
    p->hipt = u;
    p->father = father;
    p->nodeEmpty = false;
    p->size = u - low + 1;
    if (u-low+1 <= kdCutOff) 
    {
        p->bucket = true;
        p->lowson = NULL;
        p->hison  = NULL;
        
        for (i=low; i<=u; i++)
            _bucketptr[_perm[i]] = p;
        
    }
    else 
    {
        p->bucket = false;
        p->cutdim = _findmaxspread(low,u);
        m = (low + u) / 2;
        _newSelect(low, u, m, p->cutdim);
        p->cutval = _optXIn(_perm[m],p->cutdim);
        p->lowson = _kdTreeBuild(low,m,p);
        p->hison  = _kdTreeBuild(m+1,u,p);
    }
    return p;
}

void KenKDTree::_rnnMultiple(kdIter p)
{
    int thisdist;
    int val;
    int thisx;
    bool flag;
    int i;
    unsigned k;
    unsigned low, hi, size;
    
    low = p->lowpt;
    hi = p->hipt;
    size = p->size;
    if  (p->bucket) 
    {
        for (i = p->lowpt; i<= p->hipt; i++) 
        {
            flag = true;
            k=_POEnd[_nntarget];	  	/* check a PO constraint 3/96  */
            
            while (flag && (k!=UINT_MAX)) 
            { 	/* check a PO constraint 3/96  */
                if (_CellPO[k].from == _perm[i])
                    flag = false; 	
                else
                    k = _CellPO[k].tnext;
	           }
            
	           if (flag)
                   thisdist = realDist(_nntarget,_perm[i]);	//was geomDist.  Difference?
               else
                   thisdist = INT_MAX/2-1; //subtract one to make it less than initial _nndist
                              
               if (thisdist < _nndist) 
               {
                   if (_currentNum >= _nnear)
                       _extractTopHeap();
                   else
                       _currentNum++;
                   
                   _insertHeap(_perm[i],-thisdist);
                   
                   if (_currentNum >= _nnear)
                       _nndist = -_D[_A[0]];
               }
        }
    }   
    else if (_currentNum+size <= _nnear) 
    {
        _insertNode( p->lowson );
        _insertNode( p->hison );
        _currentNum += size;
        if (_currentNum >= _nnear)
            _nndist = -_D[_A[0]];
    }
    else 
    {
        val = p->cutval;
        thisx = _optXOut(_nntarget, p->cutdim);
        /*thisx = xOut(nntarget, p->cutdim);*/
        if (thisx < val) 
        {
            _rnnMultiple( p->lowson);
            if (thisx + _nndist > val)
                _rnnMultiple(p->hison);
        }
        else 
        {
            _rnnMultiple( p->hison);
            if (thisx - _nndist < val)
                _rnnMultiple( p->lowson);
        }
    }
    
}  /* rnnMultiple() */

void KenKDTree::_insertHeap(unsigned  index,int d )
{
    int i, p;
    
    _D[index] = d;
    i = _heapSize;
    _heapSize++;
    p = parent(i);
    while (i>0 && _D[_A[p]] > d) {
        _A[i] = _A[p];
        i = p;
        p = parent(i);
    }
    _A[i] = index;
    
} /* insertHeap() */


unsigned KenKDTree::_extractTopHeap()
{
    int minIndex;
    
    abkfatal(_heapSize>=1,"Heap underflow");
    
    minIndex = _A[0];
    _heapSize--;
    _A[0] = _A[_heapSize];
    _heapify( 0 );
    
    return minIndex;
    
} /*  extractTopHeap() */

void KenKDTree::_heapify(unsigned i )
{
    unsigned L, R, smallest, tInt;
    int d;
    
    L = left( i );
    R = right( i );
    d = _D[_A[i]];
    
    if (L < _heapSize && _D[_A[L]] < d)  /* note that want MAX value at top of heap
        */
        smallest = L;
    else
        smallest = i;
    
    if (R < _heapSize && _D[_A[R]] < _D[_A[smallest]])
        smallest = R;
    
    if (smallest != i) {
        tInt = _A[i];
        
        _A[i] = _A[smallest];
        _A[smallest] = tInt;
        _heapify(smallest);
    }
    
} /* heapify() */

int KenKDTree::_optXIn(unsigned i,unsigned d)
{
    if (d==0)
        return  _scanCell[i].inX;
    else
        return  _scanCell[i].inY;
}

int KenKDTree::_optXOut(unsigned i,unsigned d)
{
    if (d==0)
        return  _scanCell[i].outX;
    else
        return  _scanCell[i].outY;
}

void KenKDTree::_newSelect(unsigned low,unsigned u,unsigned m,unsigned cutdim )
{
    unsigned q;
    
    if (low==u)
        return;
    
    q = _runPartition( low, u, m, cutdim );
    if (m <= q)
        _newSelect( low, q, m, cutdim );
    else
        _newSelect( q+1, u, m, cutdim );
    
} /* newSelect() */

unsigned KenKDTree::_runPartition(unsigned low,unsigned hi,unsigned m,unsigned cutdim )
{
    unsigned  i, j;
    unsigned pVal = _perm[m];
    //_perm.swapValues(m,low);
    std::swap(_perm[m],_perm[low]);
    
    //perm[m] = perm[low]; perm[low] = pVal;
    i = low-1;
    j = hi+1;
    if (cutdim==0)
    {
        while (true) 
        {
            do 
            {j-=1;} while (_sortX(pVal,_perm[j]));
            do 
            {i+=1;} while (_sortX(_perm[i],pVal));
            if (i<j) 
            {
                //_perm.swapValues(i,j);
                std::swap(_perm[i],_perm[j]);
            }
            else
                return j;
        }
    }
    else
    {
        while (true) 
        {
            do 
            {j-=1;} while (_sortY(pVal,_perm[j]));
            do 
            {i+=1;} while (_sortY(_perm[i],pVal));
            if (i<j) 
            {
                //_perm.swapValues(i,j);
                std::swap(_perm[i],_perm[j]);
            }
            else
                return j;
        }
    }
    
} /* runPartition() */
bool KenKDTree::sortByX::operator ()(const unsigned i1, const unsigned i2)
{
    if (_scanCell[i1].inX < _scanCell[i2].inX)
        return true;
    else if (_scanCell[i1].inX > _scanCell[i2].inX)
        return false;
    else if (_scanCell[i1].inY < _scanCell[i2].inY)
        return true;
    else if (_scanCell[i1].inY > _scanCell[i2].inY)
        return false;
    else if (i1<i2)
        return true;
    return false;
}

bool KenKDTree::sortByY::operator ()(const unsigned i1, const unsigned i2)
{
    if (_scanCell[i1].inY < _scanCell[i2].inY)
        return true;
    else if (_scanCell[i1].inY > _scanCell[i2].inY)
        return false;
    else if (_scanCell[i1].inX < _scanCell[i2].inX)
        return true;
    else if (_scanCell[i1].inX > _scanCell[i2].inX)
        return false;
    else if (i1<i2)
        return true;
    return false;
}

int KenKDTree::_findmaxspread(unsigned low,unsigned u)
{  /* returns the dimension with largest difference between min and max
    among the points in perm[low...u] */
    unsigned i, d, maxD;
    int val, min, max, maxDiff;
    
    maxDiff = 0;
    for (d=0;d<2;d++) 
    {
        val = _optXIn(_perm[u],d);
        min = val;
        max = val;
        for (i=low; i<u; i++) 
        {
            val = _optXIn(_perm[i],d);
            if (val > max) max = val;
            else if (val < min) min = val;
        }
        val = max - min;
        /*
        if (d==1)
            val *= _vweight;
            */
        if (val > maxDiff)
        {
            maxDiff = val;
            maxD = d;
        }
    } /*  for (d=0;d<2;d++) */
    
    return maxD;
    
} /* findmaxspread() */

/*************************************************************/
void KenKDTree::_undeleteNode(unsigned  pointnum )
{
    kdIter p;
    unsigned j;
    
    p = _bucketptr[pointnum];
    (p->hipt)++;
    j = p->hipt;
    while (_perm[j] != pointnum)
        j++;
    //_perm.swapValues(j, p->hipt);
    std::swap(_perm[j],_perm[p->hipt]);
    p->nodeEmpty = false;
    (p->size)++;
    while ((p=p->father)) {
        p->nodeEmpty=false;
        (p->size)++;
    }
    
} /* undelete() */
/*************************************************************/
void KenKDTree::_deleteNode(unsigned  pointnum )
{
    kdnode *p;
    unsigned j;
    
    p = _bucketptr[pointnum];
    j = p->lowpt;
    while (_perm[j] != pointnum)
        j++;
    //_perm.swapValues(j, p->hipt);
    std::swap(_perm[j],_perm[p->hipt]);
    (p->hipt)--;
    (p->size)--;
    if (p->lowpt > p->hipt) p->nodeEmpty = true;
    /*while ((p=p->father) && p->lowson->empty && p->hison->empty)
    p->empty=TRUE;*/
    while ((p=p->father))
    {
        if (p->lowson->nodeEmpty && p->hison->nodeEmpty)
            p->nodeEmpty=true;
        (p->size)--;	
    }
    
} /* delete() */

void KenKDTree::_insertNode(kdIter p)
{
    int i,flag;
    unsigned k;
    int thisdist;
    
    if (p->nodeEmpty) return;
    
    if (!(p->bucket))
    {
        _insertNode( p->lowson);
        _insertNode( p->hison);
    }
    else
    {
        for (i=p->lowpt; i<=p->hipt; i++) 
        {
            flag = true;
            k=_POEnd[_nntarget];           /* check a PO constraint 3/96  */
            while (flag && (k!=UINT_MAX))
            {     /* check a PO constraint 3/96  */
                if (_CellPO[k].from == _perm[i])
                    flag = false;
                else
                    k = _CellPO[k].tnext;
            }
            
            if (flag)
                thisdist = realDist(_nntarget,_perm[i]);//was geomDist.  Diff?
            else
                thisdist = INT_MAX/2-1;
            
            
            _insertHeap(_perm[i],-thisdist);
        }
    }
        
} /* insertNode() */

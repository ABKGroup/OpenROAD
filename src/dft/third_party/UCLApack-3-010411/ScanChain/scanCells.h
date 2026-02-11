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



// scanCells.h: interface for the ScanCells class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_SCANCELLS_H____INCLUDED_)
#define _SCANCELLS_H____INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include <ABKCommon/abkcommon.h>
#include <string>
namespace abkscanopt
{
using std::vector;
using std::string;
class ScanChain;
class Cell
{
public:

    // "Named constructors" -- look up in Stroustroup to see
    // if this is done right.
    void r1Cell(const string &name,RandomRawUnsigned &randuns);
    void normCell(double xCenter,double yCenter,
        double sigma, double offset, const string &name,
        RandomNormal &randnorm);
    void uniformCell(double xCenter,double yCenter,
        double sigma, double offset,const string &name,
        RandomRawDouble &randdoub);
    void externalCell(double xCenter,double yCenter,double sigma,
        double offset,const string &name,RandomRawDouble &randdoub);

    string _name;
        //LOOK! we'll probably privatize these once more code converted
    int inX, inY;
    int outX, outY;
    ScanChain *subgroup;  //This pointer is non-NULL just in case
                          //this is a fictitious cell representing
                          //a subgroup.  In that case it's the responsibilty
                          //of other code to set inX,inY,outX,outY
                          //to the in/outs of cells in the subgroup
    Cell():inX(INT_MAX),inY(INT_MAX),outX(INT_MAX),outY(INT_MAX),
        subgroup(NULL){}

    const string &getName() const {return _name;}

    void print(ostream &os,const char TabString[]) const;
} ;

struct POType
{
    unsigned from, to;
    unsigned fnext, tnext; //index of next po elt from or to same city
    bool     basisElt;  //no combination of other pos implies this one
    POType():from(UINT_MAX),
             to(UINT_MAX),
             fnext(UINT_MAX),
             tnext(UINT_MAX),
             basisElt(false){}
             
};


class ScanCells  
{
private:
    
    // _scanCell and _poEnd are parallel arrays.
    // _poHeadTos[i] is the index (in _pos) of the
    // first partial order element pointing to
    // _scanCell[i], in the linked list of all these.
    // If _poEnd[i] is UINT_MAX, then there are no p.o. elements
    // pointing to _scanCell[i].
    // _poHeadFroms[i] is similar, but begins the linked list
    // of pos leading *from* _scanCell[i].
    // _poTailTos and _poTailFroms are the *tails* of
    // these linked lists (needed to be able to add to the end).

    unsigned _numSubgroups;
    unsigned _numCells;
    
    vector<Cell> _scanCell;
    vector<unsigned>     _poHeadTos;
    vector<unsigned>     _poHeadFroms;
    vector<unsigned>     _poTailTos;
    vector<unsigned>     _poTailFroms;
    
    vector<POType>       _pos;
    
    unsigned _numBasisPOElts;  //number of po elts that can't be deduced
                               //by transitivity from other po elts

    bool                 _isSubpath; //true if partial order completely
                                     //determines the path

    vector<unsigned>     _legalins;  // indices of cities allowed to be
                                     // the input of the path

    vector<unsigned>     _legalouts; // same as above, but outputs

    string _name;

    //must be called every time _scanCell increases in size
    void _addSpace();

    //low-level adding of a partial-order element (does
    //not compute transitive consequences.
    void _addPOElt(unsigned from,unsigned to,bool basisElt);

    //This private version is called only by _transAddPOElt(),
    //which is why it can figure out whether a PO is a basis
    //element or not.
    bool _poConstraintExists(unsigned from,unsigned to,
                             bool maintainBasisStatus);

protected:


    //the subgroup will be owned by this ScanCells object
    void _addSubgroup(ScanChain *subgroup);

    //for adding real cells, not subgroups
    void _addCell(int inX, int inY, int outX, int outY,const string &name);
    void _addCell(const Cell &cell);

    void _addLegalin(int legalin);
    void _addLegalout(int leagalout);


    bool _transAddPOElt(unsigned from,unsigned to);
    //LOOK! we can probably replace this with random_shuffle.
    //      I want to leave as is for now to get the same
    //      results as old code, for comparison purposes.
    static void _getRandomPermutation(vector<unsigned> &randarray,
                               unsigned numElts,
                               RandomRawUnsigned &randuns);

    void _clearPO();

    void _chooseLegalIns();
    void _chooseLegalOuts();

    void _setName(const string &name);
    // puts the partial order into normal form, with basis
    // elements at the beginning.  Make sure you call this
    // after adding po elements!
    void _normalizePO();

    //This ctor creates an empty ScanCells object.  Fill it
    //in by calling other protected methods of ScanCells,
    //in the following order:  First call _addSubgroup()
    //to add as many subgroups as you're going to (the ScanCells
    //dtor will delete these, so create them with "new").  Then
    //use _addCell() to add loose cells.  Then fill in
    //any partial-order constraints with _transAddPOElt(),
    //and finalize the POs with _normalizePO().  Then
    //fix the legal input and output cells, either by
    //calling _addLegalin() and _addLegalout() for each
    //of them, or by calling _chooseLegalIns() and _chooseLegalOuts()
    //to make all possible cells legal inputs and outputs
    //(subject to the partial order).  Finally, call
    // _checkIndexRanges() and _checkInOutLegality()
    //to make sure everything is consistent.
    ScanCells();

    //"Collapsing ctor"
    //This ctor creates a ScanCells object derived from
    //the one passed to it, except that the cells
    //indexed by inIdx and outIdx will be identified
    //with each other and assigned to index 0 (all other
    //indices adjusted accordingly).  Also, the new
    //ScanCells object will be flat (all cells will have
    //NULL for their "subgroup" member).
    //
    //The intended use is by ScanTourDZ
    //The only reason "orig" is not const is that we
    //have to call _fillInCellInfo() on it (thus its
    //subgroups had better have valid paths)
    ScanCells(ScanCells &orig,unsigned inIdx,unsigned outIdx);

    //This looks at all cells which are subgroups, and
    //makes sure their inX,inY,outX,outY are set correctly.
    //Don't use unless all subgroups have valid paths.
    void _fillInCellInfo();

    void _checkIndexRanges();    
    void _checkInOutLegality();

    
public:


    //use this method to determine if "from" must come before "to",
    //without modifying any basisElt flags
    bool poConstraintExists(unsigned from,unsigned to) const;

    string identifyIndex(unsigned city) const;
    void printTopDown(ostream &os,unsigned nTabs) const;
    
    virtual ~ScanCells();

    inline static int realDistStatic(int X1,int Y1, int X2, int Y2)
    {return abs(X2-X1)+abs(Y2-Y1);}
    
    inline int realDist(unsigned i,unsigned j) const
    {return realDistStatic(_scanCell[i].outX,_scanCell[i].outY,
                           _scanCell[j].inX,_scanCell[j].inY);}
    //{return abs(_scanCell[i].outX - _scanCell[j].inX) + 
    //        abs(_scanCell[i].outY - _scanCell[j].inY) ; }
    
    inline unsigned getNumCells() const {return _numCells;}
    
    inline unsigned getNumSubgroups() const {return _numSubgroups;}
    inline vector<POType> const &getPOs() const {return _pos;}
    inline vector<unsigned> const &getPOHeadTos() const {return _poHeadTos;}
    inline vector<unsigned> const &getPOHeadFroms() const {return _poHeadFroms;}
    inline vector<Cell> const &getCells() const {return _scanCell;}
    inline vector<unsigned> const &getLegalins() const {return _legalins;}
    inline vector<unsigned> const &getLegalouts() const {return _legalouts;}

    bool isSubpath() const {return _isSubpath;}
    bool hasOptimizableSubgroups() const;

    const string &getName() const {return _name;}



    //utility for determining the new index when input/output
    //index are identified as 0 (used by "collapsing ctor"
    //of ScanCells.
    static unsigned collapseIndex(int idx,int lesser,int greater);
    static unsigned uncollapseIndex(int idx,int lesser,int greater);
    
    //modify path of subgroup, changing in/out of corresp cell
    void subModifyPath(unsigned grpIdx,unsigned newIn,unsigned newOut);
};
};
ostream &operator<<(ostream &os, const abkscanopt::ScanCells &cells);
#endif // !defined(_SCANCELLS_H____INCLUDED_)

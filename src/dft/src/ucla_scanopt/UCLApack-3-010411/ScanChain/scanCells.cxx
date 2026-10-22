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
// scanCells.cxx: implementation of the ScanCells class.
//
//////////////////////////////////////////////////////////////////////

#include "scanCells.h"
#include "scanChain.h"
using namespace abkscanopt;
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//"collapsing ctor"
ScanCells::ScanCells(ScanCells &orig,unsigned inIdx,unsigned outIdx):
_numSubgroups(0),
                     _numCells(0),
                     _scanCell(0),
                     _poHeadTos(0),
                     _poHeadFroms(0),
                     _poTailTos(0),
                     _poTailFroms(0),
                     _pos(0),_numBasisPOElts(0),
                     _isSubpath(false),
                     _legalins(0),_legalouts(0)
{
    orig._fillInCellInfo();
    _scanCell.reserve(orig.getNumCells()-1);
    const vector<Cell> &ocells=orig._scanCell;
    unsigned lesser=inIdx,greater=outIdx;
    if (lesser>greater) std::swap(lesser,greater);

    //the "zero cell" made by identifying the input and output
    //cells of the path
    _addCell(ocells[outIdx].inX,ocells[outIdx].inY,
             ocells[inIdx].outX,ocells[inIdx].outY,"");

    //Now we copy the other cells.  Note that subgroup pointers
    //do *not* get copied.
    unsigned i;
    for (i=0;i<lesser;i++)
    {
        const Cell &cell=ocells[i];
        _addCell(cell);
        _poHeadTos.back()=orig._poHeadTos[i];
        _poHeadFroms.back()=orig._poHeadFroms[i];
        _poTailTos.back()=orig._poTailTos[i];
        _poTailFroms.back()=orig._poTailFroms[i];
    }
    for (i=lesser+1;i<greater;i++)
    {
        const Cell &cell=ocells[i];
        _addCell(cell);
        _poHeadTos.back()=orig._poHeadTos[i];
        _poHeadFroms.back()=orig._poHeadFroms[i];
        _poTailTos.back()=orig._poTailTos[i];
        _poTailFroms.back()=orig._poTailFroms[i];
    }
    for (i=greater+1;i<orig.getNumCells();i++)
    {
        const Cell &cell=ocells[i];
        _addCell(cell);
        _poHeadTos.back()=orig._poHeadTos[i];
        _poHeadFroms.back()=orig._poHeadFroms[i];
        _poTailTos.back()=orig._poTailTos[i];
        _poTailFroms.back()=orig._poTailFroms[i];
    }

    //Now we translate partial-order constraints

    const vector<POType> &opos=orig.getPOs();
    const unsigned nPOs=opos.size();
    POType blank;
    _pos.insert(_pos.end(),nPOs,blank);
    for (i=0;i<nPOs;i++)
    {
        POType &newPo=_pos[i];
        const POType &oldPo=opos[i];
        newPo.from=ScanCells::collapseIndex(oldPo.from,lesser,greater);
        newPo.to=ScanCells::collapseIndex(oldPo.to,lesser,greater);
        newPo.fnext=oldPo.fnext;
        newPo.tnext=oldPo.tnext;
    }

    const vector<unsigned> &olos=orig.getLegalouts();
    const unsigned nolos=olos.size();
    _legalouts.reserve(nolos);
    for (i=0;i<nolos;i++)
    {
        unsigned lo=ScanCells::collapseIndex(olos[i],lesser,greater);
        if (lo!=0)
            _legalouts.push_back(lo);
    }

    _legalins.push_back(0); // 0 must be the first cell in path
    _checkIndexRanges();
    _checkInOutLegality();
}



ScanCells::~ScanCells()
{
    unsigned i;
    for (i=0;i<_numSubgroups;i++)
    {
        Cell &cell=_scanCell[i];
        abkfatal(cell.subgroup !=NULL,"Unexpected error on exit: "
            "wrong num subgroups");
        delete cell.subgroup;
    }
}

void Cell::normCell(double xCenter,double yCenter,double sigma,double offset,
               const string &name,RandomNormal &randnorm)
{
    double dinX  = xCenter + randnorm*sigma;
    double dinY  = yCenter + randnorm*sigma;
    inX =static_cast<int>(dinX);
    inY =static_cast<int>(dinY);
    outX = static_cast<int>(dinX + offset);
    outY = inY;
    _name=name;

}

void Cell::r1Cell(const string &name,RandomRawUnsigned &randuns)
{
    int Pos = randuns%10000;
    inX = Pos;
    Pos = randuns%10000;
    inY = Pos;
    int Dif = randuns%10-5;
    outX = inX + Dif;
    Dif = randuns%10-5;
    outY = inY + Dif;
    _name=name;
}

void Cell::uniformCell(double xCenter,double yCenter,double sigma,
                       double offset,const string &name,
                       RandomRawDouble &randdoub)
{
    double dinX = xCenter + (randdoub-0.5)*sigma*pow(12.0,0.5);
    inX = static_cast<int>(dinX);
    double dinY = yCenter + (randdoub-0.5)*sigma*pow(12.0,0.5);
    inY = static_cast<int>(dinY);
    outX = static_cast<int>(dinX + offset);
    outY = inY;
    _name=name;
}

void Cell::externalCell(double xCenter,double yCenter,double sigma,
        double offset,const string &name,RandomRawDouble &randdoub)
{
    int sign;
    int XorY;
    double dinX;
    
    XorY = (int)(randdoub*2);
    
    if (XorY==0)
    {
        sign = ((int)(randdoub*2))*2 - 1; // +1 or -1
        dinX = xCenter + sign*0.5*sigma*pow(12.0,0.5);
        inY = static_cast<int>(yCenter + (randdoub-0.5)
                                       *sigma*pow(12.0,0.5));
    }
    else
    {
        sign = ((int)(randdoub*2))*2 - 1; // +1 or -1
        inY = static_cast<int>(yCenter + sign*0.5*sigma*pow(12.0,0.5));
        dinX = xCenter + (randdoub-0.5)
            *sigma*pow(12.0,0.5);
    }

    inX=static_cast<int>(dinX);
    
    outX = static_cast<int>(dinX + offset);
    outY = inY;
    _name=name;
}



void Cell::print(ostream &os,const char TabString[]) const
{
    
    //We've removed the possibility of multiple ins or
    //multiple outs for a cell
    os << TabString << "    CELL \""<< _name.c_str()
        <<"\" (INs 1 OUTs 1)\n";
    
    os << TabString <<"    INs:\n";
    
    os << TabString << "    (0)    IN "
        << inX << " " << inY <<"\n";

    os << TabString << "    ENDINs\n";
    os << TabString <<"    OUTs:\n";
    
    os << TabString << "    (0)    OUT "
        << outX << " " << outY <<"\n";

    os << TabString << "    ENDOUTs\n";
    
    os << TabString << "    ENDCELL\n";
}
ScanCells::ScanCells()                     
                     :_numSubgroups(0),
                     _numCells(0),
                     _scanCell(0),
                     _poHeadTos(0),
                     _poHeadFroms(0),
                     _poTailTos(0),
                     _poTailFroms(0),
                     _pos(0),_numBasisPOElts(0),
                     _isSubpath(false),
                     _legalins(0),_legalouts(0)

{
}
bool ScanCells::_transAddPOElt(unsigned from,unsigned to)
{

    
    // this array will contain "from" and the indices
    // of all cities constrained to come before "from"
    vector<unsigned> FromCityIdxArray(_numCells,UINT_MAX);
    
    // this array will contain "to" and the indices
    // of all cities constrained to come after "to"
    vector<unsigned> ToCityIdxArray(_numCells,UINT_MAX);
    
    // This var counts the number of cities constrained to come
    // before the one given by "from".
    unsigned NumberFromCities=0;
    
    // This var counts the number of cities constrained to come
    // after the one given by "to"
    unsigned NumberToCities=0;
    
    unsigned iPOElt;
    unsigned i,j;
    
    FromCityIdxArray[0] = from ;
    
    iPOElt=_poHeadTos[from];
    if (iPOElt!=UINT_MAX)
    {
        do
        {
            const POType &po=_pos[iPOElt];
            FromCityIdxArray[++NumberFromCities] = po.from;
            iPOElt=po.tnext;
            
        } while (iPOElt != UINT_MAX);
    }
    
    ToCityIdxArray[0] = to ;
    
    iPOElt=_poHeadFroms[to];
    if (iPOElt!=UINT_MAX)
    {
        do
        {
            const POType &po=_pos[iPOElt];
            ToCityIdxArray[++NumberToCities] = po.to;
            iPOElt=po.fnext;
            
        } while (iPOElt != UINT_MAX);
    }
    
    // the " <= " is because we didn't count FromCityIdx itself.
    // similarly in the "for j" loop.
    for (i=0;i<=NumberFromCities;i++)
    {
        for (j=0;j<=NumberToCities;j++)
        {
            // The function POConstraintExists() will take care of
            // marking transitive consequences as non-basis elts
            
            if (!_poConstraintExists(FromCityIdxArray[i],
                ToCityIdxArray[j],
                (i!=0 || j!=0)))
                
                // We mark the direct PO element as a basis elt;
                // this may change later.
                _addPOElt(FromCityIdxArray[i],
                ToCityIdxArray[j],
                (i==0 && j==0));
            
            //#ifdef _DEBUG
            if (_poConstraintExists(ToCityIdxArray[j],
                FromCityIdxArray[i],false))
            {
                abkwarn(false,"Cycle in partial order");
                return false;
            }
            //#endif
        }
    }
    return true;
}


bool ScanCells::_poConstraintExists(unsigned from,unsigned to,
                             bool maintainBasisStatus)
{
    unsigned  iPOElt;
    iPOElt=_poHeadFroms[from];
    
    if (iPOElt!=UINT_MAX)
    {
        do
        {
            POType &po=_pos[iPOElt];
            if (po.to == to)
            {
                if (maintainBasisStatus)
                    po.basisElt = false;
                return true;
            }
            iPOElt=po.fnext;
        } while (iPOElt!=UINT_MAX);
    }
    
    return false;
}

bool ScanCells::poConstraintExists(unsigned from,unsigned to) const
{
    unsigned  iPOElt;
    iPOElt=_poHeadFroms[from];
    
    if (iPOElt!=UINT_MAX)
    {
        do
        {
            const POType &po=_pos[iPOElt];
            if (po.to == to) return true;
            iPOElt=po.fnext;
        } while (iPOElt!=UINT_MAX);
    }
    
    return false;
}

void ScanCells::_addPOElt(unsigned from,unsigned to,bool basisElt)
{
    POType p;
    const unsigned index=_pos.size();
    _pos.push_back(p);
    POType &po=_pos[index];
    
    po.from=from;
    po.to=to;
    po.basisElt=basisElt;

    if (_poHeadFroms[from]==UINT_MAX)
        _poHeadFroms[from]=index;
        
    if (_poHeadTos[to]==UINT_MAX)
        _poHeadTos[to]=index;
        
    const unsigned tailfromidx=_poTailFroms[from];
    if (tailfromidx!=UINT_MAX)
        _pos[tailfromidx].fnext=index;
    _poTailFroms[from]=index;

    const unsigned tailtoidx=_poTailTos[to];
    if (tailtoidx!=UINT_MAX)
        _pos[tailtoidx].tnext=index;
        
    _poTailTos[to]=index;
        
}


//LOOK! we can probably replace this with random_shuffle.
//      I want to leave as is for now to get the same
//      results as old code, for comparison purposes.
void ScanCells::_getRandomPermutation(vector<unsigned> &randarray,
                               unsigned numElts,
                               RandomRawUnsigned &randuns)
{
    randarray.clear();
    randarray.insert(randarray.end(),numElts,UINT_MAX);
    unsigned i;
    unsigned randIdx;
    for (i=0;i<numElts;i++)
        randarray[i] = i;
    
    
    for (i=numElts-1;i!=static_cast<unsigned>(-1);i--)
    {
        randIdx = randuns%(i+1);
        std::swap(randarray[i],randarray[randIdx]);
    }
    
}

void ScanCells::_clearPO()
{
    _numBasisPOElts=0;
    _pos.clear();
    _poHeadFroms.clear();
    _poHeadFroms.insert(_poHeadFroms.end(),_numCells,UINT_MAX);
    _poHeadTos.clear();
    _poHeadTos.insert(_poHeadTos.end(),_numCells,UINT_MAX);
    _poTailFroms.clear();
    _poTailFroms.insert(_poTailFroms.end(),_numCells,UINT_MAX);
    _poTailTos.clear();
    _poTailTos.insert(_poTailTos.end(),_numCells,UINT_MAX);

}



void ScanCells::_normalizePO()
{
    const unsigned nTransPOElts=_pos.size();
    unsigned i;
    unsigned n;
    
    n=_numCells;
    if (2*nTransPOElts == n*(n-1))
    {
        _isSubpath = true;    // We have the max number of POs, so this
        return;               // group is completely constrained.
    }

    vector<POType> TempPOArray=_pos;
    _clearPO();
    
    //first pass:  _transAddPOElt() will make sure that, as you
    //             add basis elements (decided previously by
    //             _transAddPOElt() ), consequences always follow
    //             premisses.
    
    for (i=0;i<nTransPOElts;i++)
    {
        const POType &temppo=TempPOArray[i];
        if (TempPOArray[i].basisElt)
        {
            _transAddPOElt(temppo.from,temppo.to);
        }
    }

    TempPOArray=_pos; //copy back again, after adding basis elements
    _clearPO();
    
    //second pass:  Now we make sure basis elements come at the beginning
    
    for (i=0;i<nTransPOElts;i++)
    {
        const POType &temppo=TempPOArray[i];
        if (temppo.basisElt)
        {
            _addPOElt(temppo.from,temppo.to,true);
            _numBasisPOElts++;
        }
        
    }
    
    for (i=0;i<nTransPOElts;i++)
    {
        const POType &temppo=TempPOArray[i];
        if (!temppo.basisElt)
        {
            _addPOElt(temppo.from,temppo.to,false);
        }
        
    }
    
}

// This determines which cities are legal inputs
void ScanCells::_chooseLegalIns()
{
    _legalins.clear();
    unsigned i;
    for (i=0;i<_numCells;i++)
    {
        if (_poHeadTos[i]==UINT_MAX)
            _legalins.push_back(i);
    }
}
void ScanCells::_chooseLegalOuts()
{
    _legalouts.clear();
    unsigned i;
    for (i=0;i<_numCells;i++)
    {
        if (_poHeadFroms[i]==UINT_MAX)
            _legalouts.push_back(i);
    }
}

string ScanCells::identifyIndex(unsigned city) const
{
    if (city < _numSubgroups)
        return string("group \"")+_scanCell[city].subgroup->getName()+"\"";
    else
        return string("cell \"")+_scanCell[city].getName()+"\"";
}

void ScanCells::printTopDown(ostream &os,unsigned nTabs) const
{
    char TabString[255]; // = "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t";
    unsigned i;
    const unsigned nPOs=_pos.size();
    
    TabString[4*nTabs] = '\0';
    memset(TabString,' ',4*nTabs);
    
    os << TabString << "GROUP \""<< _name.c_str()
        <<"\" (POs "<< nPOs
        <<" SUBGROUPs "<<_numSubgroups
        <<" LOOSECELLs " << _numCells-_numSubgroups 
        <<" LEGALINs " << _legalins.size()
        << " LEGALOUTs " << _legalouts.size() <<")\n\n";
    os << "\n";
    
    os << TabString << "SUBGROUPs (group \"" << _name.c_str() <<"\"):\n";
    
    for (i=0;i<_numSubgroups;i++)
    {
        os << TabString << "(" << i << ")\n";
        _scanCell[i].subgroup->printTopDown(os,nTabs+1);
    }
    os << TabString << "ENDSUBGROUPs (group \""<< _name.c_str() <<"\")\n\n";
    
    os << TabString << "LOOSECELLs (group \"" << _name.c_str() <<"\"):\n";
    
    for (i=_numSubgroups;i<_numCells;i++)
    {
        os << TabString <<  "("<< i<< ")\n";
        _scanCell[i].print(os,TabString);
    }
    
    os << TabString << "ENDLOOSECELLs (group \""<<_name.c_str()<<"\")\n\n";
    
    os << TabString << "POs (group \"" << _name.c_str() << "\"):\n";

    
    for (i=0;i<nPOs;i++)
    {
        os << TabString << "(" << i <<")    FROM " << _pos[i].from <<
            " (" << identifyIndex(_pos[i].from).c_str() << ")\n";
        
        os << TabString <<"                                 TO "
            << _pos[i].to <<" ("<< identifyIndex(_pos[i].to).c_str()
            <<") ";
        if (_pos[i].basisElt)
            os <<"(basis element)";
        os << "\n";
        
    }
    
    os << TabString << "ENDPOs (group \""
        << _name.c_str() << "\")\n";
    
    os << "\n";
    
    os << TabString << "LEGALINs (group \"" << _name.c_str() <<"\"):\n";
    
    const unsigned nLegalins=_legalins.size();
    for (i=0;i<nLegalins;i++)
    {
        os << TabString << "("<<i<<")    LEGALIN "<<_legalins[i]
            <<" ("<<identifyIndex(_legalins[i]).c_str()<<")\n";
    }
    
    os << TabString << "ENDLEGALINs (group \"" << _name.c_str()
        << "\")\n";
    
    os <<"\n";
    
    os << TabString << "LEGALOUTs (group \"" << _name.c_str() <<"\"):\n";
    
    const unsigned nLegalouts=_legalouts.size();
    for (i=0;i<nLegalouts;i++)
    {
        os << TabString << "("<< i << ")    LEGALOUT " << _legalouts[i] <<
            " (" << identifyIndex(_legalouts[i]).c_str() <<")\n";
    }
    
    
    os << TabString << "ENDLEGALOUTs (group \"" << _name.c_str() <<"\")\n";
    
    os << TabString << 
    "ENDGROUP (group \"" << _name.c_str() <<"\")\n\n";
    os.flush();
}

void ScanCells::_addSpace()
{
    _numCells++;
    _poHeadFroms.push_back(UINT_MAX);
    _poHeadTos.push_back(UINT_MAX);
    _poTailFroms.push_back(UINT_MAX);
    _poTailTos.push_back(UINT_MAX);
}

void ScanCells::_addCell(const Cell &cell)
{
    _addCell(cell.inX,cell.inY,cell.outX,cell.outY,cell._name);
}

void ScanCells::_addCell(int inX,int inY, int outX,int outY,const string &name)
{
    Cell tempCell;
    tempCell.inX=inX;tempCell.inY=inY;tempCell.outX=outX;tempCell.outY=outY;
    tempCell._name=name;
    tempCell.subgroup=NULL;
    _scanCell.push_back(tempCell);
    _addSpace();
}

void ScanCells::_addSubgroup(ScanChain *subgroup)
{
    abkfatal(_numCells==_numSubgroups,
        "Can''t add a subgroup after adding loose cells");
    Cell tempCell;
    tempCell.subgroup=subgroup;
    tempCell._name=subgroup->getName();
    _scanCell.push_back(tempCell);
    _addSpace();
    _numSubgroups++;
}

void ScanCells::_addLegalin(int legalin)
{
    _legalins.push_back(legalin);
}

void ScanCells::_addLegalout(int legalout)
{
    _legalouts.push_back(legalout);
}

void ScanCells::_setName(const string &name)
{
    _name=name;
}

ostream &operator<<(ostream& os, const ScanCells &cells)
{
    cells.printTopDown(os,0);
    return os;
}
//************************************************************************
void ScanCells::_checkInOutLegality()
{
    unsigned i;
    unsigned CellIdx;
    unsigned iPO;
    const vector<unsigned> &legalins=getLegalins();
    const unsigned nLegalins=legalins.size();
    const vector<unsigned> &poHeadTos=getPOHeadTos();
    const vector<POType> &pos=getPOs();
    
    for (i=0;i<nLegalins;i++)
    {
        CellIdx = legalins[i];
        iPO=poHeadTos[CellIdx];
        if (iPO!=UINT_MAX)
        {
            char errtxt[1023];
            sprintf(errtxt,"Group %s:\n"
                "Index %d has been specified as a legal in cell.\n"
                "This is illegal because there is a partial order\n"
                "constraint requiring that index %d come before it",
                getName().c_str(),
                CellIdx,
                pos[iPO].from);
            abkfatal(false,errtxt);
        }
    }

    const vector<unsigned> &legalouts=getLegalouts();
    const unsigned nLegalouts=legalouts.size();
    const vector<unsigned> &poHeadFroms=getPOHeadFroms();
    
    for (i=0;i<nLegalouts;i++)
    {
        CellIdx = legalouts[i];
        iPO=poHeadFroms[CellIdx];
        if (iPO!=UINT_MAX)
        {
            char errtxt[1023];
            sprintf(errtxt,"Group %s:\n"
                "Index %d has been specified as a legal out cell.\n"
                "This is illegal because there is a partial order\n"
                "constraint requiring that index %d come after it",
                getName().c_str(),
                CellIdx,
                pos[iPO].to);
            abkfatal(false,errtxt);
        }
    }
}

//************************************************************************
void ScanCells::_checkIndexRanges()
{
    unsigned TopOfRange;
    unsigned i;
    unsigned From,To;
    
    TopOfRange = getNumCells();
    const vector<unsigned> &legalins=getLegalins();
    const unsigned nLegalins=legalins.size();
    
    for (i=0;i<nLegalins;i++)
    {
        if ( legalins[i] > TopOfRange)
        {
            char errtxt[1023];
            sprintf(errtxt,"Group %s:\n"
                "Index %d has been specified as a legal in cell.\n"
                "The range of available indices for this group is"
                " 0 to %d\n",getName().c_str(),legalins[i],TopOfRange);
            abkfatal(false,errtxt);
        }
    }
    const vector<unsigned> &legalouts=getLegalouts();
    const unsigned nLegalouts=legalouts.size();
    for (i=0;i<nLegalouts;i++)
    {
        if (legalouts[i] > TopOfRange)
        {
            char errtxt[1023];
            sprintf(errtxt,"Group %s:\n"
                "Index %d has been specified as a legal out cell.\n"
                "The range of available indices for this group is"
                " 0 to %d\n",getName().c_str(),legalouts[i],TopOfRange);
            abkfatal(false,errtxt);
        }
    }
    
    const vector<POType> &pos=getPOs();
    const unsigned nTrans=pos.size();
    for (i=0;i<nTrans;i++)
    {
        const POType &po=pos[i];
        From = po.from;
        To   = po.to;
        
        if (From > TopOfRange)
        {
            char errtxt[1023];
            sprintf(errtxt,"Group %s:\n"
                "Index %d has been specified in"
                " a partial order constraint (as a FROM cell).\n"
                "The range of available indices for this group is"
                " 0 to %d\n",getName().c_str(),From,TopOfRange);
            abkfatal(false,errtxt);
        }
        if (To > TopOfRange)
        {
            char errtxt[1023];
            sprintf(errtxt,"Group %s:\n"
                "Index %d has been specified in"
                " a partial order constraint (as a TO cell).\n"
                "The range of available indices for this group is"
                " 0 to %d\n",getName().c_str(),To,TopOfRange);
            abkfatal(false,errtxt);
        }
    }
    
}
//************************************************************************
void ScanCells::_fillInCellInfo()
{
    unsigned i;
    for (i=0;i<_numSubgroups;i++)
    {
        Cell &cell=_scanCell[i];
        ScanChain *pSubgroup=cell.subgroup;
        abkfatal(pSubgroup!=NULL,"Indexing error");
        abkfatal(pSubgroup->isPathValid(),
            "Can''t call _fillInCellInfo when there''s invalid path");
        pSubgroup->_fillInCellInfo(); //recurse
        const vector<unsigned> &subpath=pSubgroup->getPath();
        unsigned subNumCells=pSubgroup->getNumCells();
        unsigned subInIdx=subpath[0];
        unsigned subOutIdx=subpath[subNumCells-1];
        const vector<Cell> &subcells=pSubgroup->getCells();
        const Cell &subInCell=subcells[subInIdx];
        const Cell &subOutCell=subcells[subOutIdx];
        cell.inX=subInCell.inX;
        cell.inY=subInCell.inY;
        cell.outX=subOutCell.outX;
        cell.outY=subOutCell.outY;
    }
}
//************************************************************************
unsigned ScanCells::collapseIndex(int idx,int lesser,int greater)
{
    abkfatal(lesser!=greater,"In/out idx the same");

    if (idx == lesser || idx == greater)
        return 0;
    else if (idx < lesser)
        return idx + 1;
    else if (idx < greater)
        return idx;
    else
        return idx-1;
}
//************************************************************************
unsigned ScanCells::uncollapseIndex(int idx,int lesser,int greater)
{
    abkfatal(lesser!=greater,"In/out idx the same");
    abkfatal(idx!=0,"Can''t uncollapse zero index");
    if (idx < lesser+1)
        return idx-1;
    else if (idx < greater)
        return idx;
    else
        return idx+1;
}
//************************************************************************

bool ScanCells::hasOptimizableSubgroups() const
{
    unsigned i;
    
    for (i=0;i<_numSubgroups;i++)
    {
        if (!(_scanCell[i].subgroup->isSubpath()))
            return true;
    }
    return false;
}

void ScanCells::subModifyPath(unsigned grpIdx,unsigned newIn,unsigned newOut)
{
    abkfatal(grpIdx<_numSubgroups,"subgroup index out of range");
    Cell &cell=_scanCell[grpIdx];
    cell.subgroup->_modifyPath(newIn,newOut);
    const vector<unsigned> &subpath=cell.subgroup->getPath();
    unsigned inCellIdx=subpath[0];
    unsigned outCellIdx=subpath[subpath.size()-1];
    const vector<Cell> &subcells=cell.subgroup->getCells();
    const Cell &inCell=subcells[inCellIdx], &outCell=subcells[outCellIdx];
    cell.inX=inCell.inX;
    cell.inY=inCell.inY;
    cell.outX=outCell.outX;
    cell.outY=outCell.outY;
}


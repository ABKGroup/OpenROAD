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




// scanChainRandom.cxx:  Implements constructor for a random ScanChain
// object
#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include "scanChainRandom.h"

using namespace abkscanopt;

ScanChainRandom::ScanChainRandom(RandGenInfo info[],const string &name)                     
                     :ScanChain()

{
    _setName(name);
    unsigned i;
    unsigned from,to;
    string SubGrpName;
    string CellName;
    char NumberString[10];
    int GridSize,XGridCoord,YGridCoord;
    double SideLength;
    
    RandomRawUnsigned randuns("ScanChainRandom::ScanChainRandom,randuns");
    RandomRawDouble   randdoub("ScanChainRandom::ScanChainRandom,randdoub");
    RandomNormal      randnorm(0,1,"ScanChainRandom::ScanChainRandom,randnorm");
    
    for (i=0;i<info[0].numberSubGroups;i++)
    {
        SubGrpName=name;
        sprintf(NumberString,"%3d\0",i);
        SubGrpName+=NumberString;
        switch (info[0].subgroupsCenter)
        {
        case SAME_CENTER:
            {
                info[1].xCenter = info[0].xCenter;
                info[1].yCenter = info[0].yCenter;
                break;
            }
        case SAME_SCHEME:
            {
                switch (info[0].randCellType)
                {
                case 0:
                    break;
                case 1:
                    info[1].xCenter = randnorm*info[0].sigma
                        +info[0].xCenter;
                    
                    info[1].yCenter = randnorm*info[0].sigma
                        +info[0].yCenter;
                    break;
                case 2:
                    {
                        info[1].xCenter = (randdoub-0.5)*
                            info[0].sigma*pow(12.0,0.5)
                            +info[0].xCenter;
                        
                        info[1].yCenter = (randdoub-0.5)*
                            info[0].sigma*pow(12.0,0.5) 
                            +info[0].yCenter;
                        break;
                    }
                default:
                    abkfatal(false,"bad RandCellType\n");
                }
                break;
            }
        case GRIDDED:
            {
                GridSize = (int)ceil(sqrt(info[0].numberSubGroups));
                XGridCoord = i % GridSize;
                YGridCoord = i / GridSize;
                
                SideLength = info[0].sigma*pow(12.0,0.5);
                
                info[1].xCenter = info[0].xCenter +
                    ((2*XGridCoord+1)/(2.0*GridSize) - 0.5) *
                    SideLength;
                info[1].yCenter = info[0].yCenter +
                    ((2*YGridCoord+1)/(2.0*GridSize) - 0.5) *
                    SideLength;
            }
        } //switch (info[0].SubgroupCenter)

        _addSubgroup(new ScanChainRandom(info+1,SubGrpName));
        
    };

    int StartIndex=(info[0].externalInOut) ? 2 : 0;

    for (i=0;i<info[0].numberLooseCells-StartIndex;i++)
    {
        Cell cell;
        CellName=name;
        sprintf(NumberString,",%d",i);
        CellName += NumberString;
        switch (info[0].randCellType)
        {
        case 0:
            cell.r1Cell(CellName,randuns);
            break;
            
        case 1:
            cell.normCell(info[0].xCenter,
                info[0].yCenter,
                info[0].sigma,
                info[0].offset,
                CellName,randnorm);
            break;
            
        case 2:
            cell.uniformCell(info[0].xCenter,
                info[0].yCenter,
                info[0].sigma,
                info[0].offset,
                CellName,randdoub);
            break;
            
        default:
            abkfatal(false,"bad RandCellType\n");
        }
        _addCell(cell);
    }
      
    if (info[0].externalInOut)
    {
        Cell cell1;
        cell1.externalCell(info[0].xCenter,
            info[0].yCenter,
            info[0].sigma,
            info[0].offset,
            "INPUT",randdoub);

        _addCell(cell1);

        Cell cell2;

        cell2.externalCell(info[0].xCenter,
            info[0].yCenter,
            info[0].sigma,
            info[0].offset,
            "OUTPUT",randdoub);

        _addCell(cell2);
    }
    
    
    // We need some way to keep from getting cycles in the PO.
    // Sun Cho's solution was to add only PO elements whose "to" indices
    // were numerically greater than their "from" indices.  That won't
    // work for us because we want to be able to have PO constraints
    // between loose cells and subgroups, in either direction.  Therefore
    // we pick a random permutation of the available indices.  The
    // eventual PO will be a suborder of the inverse of this permutation.
    
    vector<unsigned> PermArray;
    _getRandomPermutation(PermArray,getNumCells()-StartIndex,randuns);
    
    const unsigned numCells=getNumCells();

    while (getPOs().size() < info[0].desiredNumberPOElts)
    {
        from  = randuns%(numCells-StartIndex);
        do
            to = randuns%(numCells-StartIndex);
        while (from == to);
        
        if (PermArray[from] > PermArray[to]) 
        {
            std::swap(from,to);
        }
        
        bool ok=_transAddPOElt(from,to);
        abkfatal(ok,"Unexpected error: Cycle in PO added while "
            "generating ScanChainRandom object randomly");
    } ;
    
    _normalizePO();
    
    
  
    if (info[0].externalInOut)
    {
        _addLegalin(numCells-2);
        _addLegalout(numCells-1);
    }
    else
        
    {
        _chooseLegalIns();
        _chooseLegalOuts();
    }
    

    _checkIndexRanges();    
    _checkInOutLegality();
}


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




// This file is simple sample code showing how to generate a random problem

#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include <ScanChain/scanChainRandom.h>

using namespace std;
using namespace abkscanopt;

//********************************************************************
int main(int argc,char **argv)
{
    SeedHandler::turnOffLogging();
    SeedHandler::overrideExternalSeed(12);
    RandomRawUnsigned randForPath("main,randForPath");
    ScanChainRandom::RandGenInfo rinfo[2];
    rinfo[0].desiredNumberPOElts=10;
    rinfo[0].externalInOut=false;
    rinfo[0].numberLooseCells=40;
    rinfo[0].numberSubGroups=2;
    rinfo[0].offset=10;
    rinfo[0].randCellType=1;
    rinfo[0].sigma=5000;
    rinfo[0].subgroupsCenter=ScanChainRandom::SAME_SCHEME;
    rinfo[0].xCenter=0;
    rinfo[0].yCenter=0;
    rinfo[1].desiredNumberPOElts=5;
    rinfo[1].externalInOut=false;
    rinfo[1].numberLooseCells=10;
    rinfo[1].numberSubGroups=0;
    rinfo[1].offset=10;
    rinfo[1].randCellType=1;
    rinfo[1].sigma=1000;
    rinfo[1].subgroupsCenter=ScanChainRandom::SAME_SCHEME;
    rinfo[1].xCenter=0;
    rinfo[1].yCenter=0;


    ScanChainRandom chain(rinfo,"Sim 0");
    chain.initialPath(randForPath);
    
    ::cout << chain ;
    chain.printPathTopDown(::cout,0);

    return 0;
}


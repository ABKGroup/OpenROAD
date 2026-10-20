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




// This file for testing random generation of problems.

#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif


#include <ScanChain/scanChainRandom.h>

using namespace std;
using namespace abkscanopt;
struct RunParams;
void parseArgs(int argc,char **argv,RunParams &runParams);
void SetRandGenParms(ScanChainRandom::RandGenInfo  randGenInfo[],
                     const RunParams &runParams);
void PrintOutHelp(RunParams &runParams);


struct RunParams
{
    std::string commandName;
    bool PrintHelp;
    long RandSeed;
    int LooseCells;
    int GroupLooseCells;
    int Subgroups;
    int POs;
    int GroupPOs;
    double Sigma;
    double GroupSigma;
    double Offset;
    bool ExternalInOut;
    bool SubgroupsSameCenter;
    bool SubgroupsGridded;
    bool UniformCellDist;
    bool GroupUniformCellDist;
    RunParams():
    commandName(""),
        PrintHelp(false),
        RandSeed(12),
        LooseCells(40),
        GroupLooseCells(10),
        Subgroups(2),
        POs(10),
        GroupPOs(5),
        Sigma(5000),
        GroupSigma(1000),
        Offset(10),
        ExternalInOut ( false),
        SubgroupsSameCenter ( false),
        SubgroupsGridded ( false),
        UniformCellDist ( false),
        GroupUniformCellDist ( false){}
    
};


//********************************************************************
int main(int argc,char **argv)
{
    ScanChainRandom::RandGenInfo  randGenInfo[2];
    //long                        PathSeed[2];
    
    RunParams runParams;
    parseArgs(argc,argv,runParams);
    if (runParams.PrintHelp)
    {
        PrintOutHelp(runParams);
        return 0;
    }
    SeedHandler::overrideExternalSeed(runParams.RandSeed);
    SeedHandler::turnOffLogging();
    SetRandGenParms(randGenInfo,runParams);
    
    

    ScanChainRandom cells(randGenInfo,"Sim 0");
    
    ::cout << cells;

    return 0;
}
//*********************************************************************
void parseArgs(int argc,char **argv,RunParams &runParams)
{
    int i;
    runParams.commandName=argv[0];
    
    for (i=1;i<argc;i++) 
    {
        if (strcmp("-h",argv[i]) == 0)  
        {
            runParams.PrintHelp = true;
        }
        else if (strcmp("-rs",argv[i]) == 0)  
        {
            i++;
            sscanf(argv[i],"%ld",&runParams.RandSeed);
        }
        else if (strcmp("-lc",argv[i]) == 0)  
        {
            i++;
            sscanf(argv[i],"%d",&runParams.LooseCells);
        }
        else if (strcmp("-glc",argv[i]) == 0)  
        {
            i++;
            sscanf(argv[i],"%d",&runParams.GroupLooseCells);
        }
        else if (strcmp("-sg",argv[i]) == 0)  
        {
            i++;
            sscanf(argv[i],"%d",&runParams.Subgroups);
        }
        else if (strcmp("-po",argv[i]) == 0)  
        {
            i++;
            sscanf(argv[i],"%d",&runParams.POs);
        }
        else if (strcmp("-gpo",argv[i]) == 0)  
        {
            i++;
            sscanf(argv[i],"%d",&runParams.GroupPOs);
        }
        else if (strcmp("-sigma",argv[i]) == 0)  
        {
            i++;
            sscanf(argv[i],"%lf",&runParams.Sigma);
        }
        else if (strcmp("-gsigma",argv[i]) == 0)  
        {
            i++;
            sscanf(argv[i],"%lf",&runParams.GroupSigma);
        }
        else if (strcmp("-offset",argv[i]) == 0)  
        {
            i++;
            sscanf(argv[i],"%lf",&runParams.Offset);
        }
        else if (strcmp("-eio",argv[i]) == 0)
        {
            runParams.ExternalInOut = true;
        }
        else if (strcmp("-ssc",argv[i]) == 0)
        {
            runParams.SubgroupsSameCenter = true;
        }
        else if (strcmp("-sgg",argv[i]) == 0)
        {
            runParams.SubgroupsGridded = true;
        }
        else if (strcmp("-unif",argv[i]) == 0)
        {
            runParams.UniformCellDist = true;
        }
        else if (strcmp("-gunif",argv[i]) == 0)
        {
            runParams.GroupUniformCellDist = true;
        }
        else 
        {
            char errtxt[1023];
            sprintf(errtxt,"command argument %d %s unrecognized; program exiting\n",
                i,argv[i]);
            abkfatal(false,errtxt);
        }
    }
}
//*********************************************************************
void SetRandGenParms(ScanChainRandom::RandGenInfo randGenInfo[],
    const RunParams &runParams)
{
    randGenInfo[0].numberLooseCells = runParams.LooseCells;
    randGenInfo[0].numberSubGroups  = runParams.Subgroups;
    randGenInfo[0].desiredNumberPOElts = runParams.POs;
    if (runParams.UniformCellDist)
        randGenInfo[0].randCellType = 2;  // use uniform distribution
    else
        randGenInfo[0].randCellType = 1;  // use normal distribution
    randGenInfo[0].xCenter = 0;
    randGenInfo[0].yCenter = 0;
    randGenInfo[0].sigma = runParams.Sigma;
    randGenInfo[0].offset = runParams.Offset;
    randGenInfo[0].externalInOut = runParams.ExternalInOut;
    
    if (runParams.SubgroupsSameCenter)
        randGenInfo[0].subgroupsCenter = ScanChainRandom::SAME_CENTER;
    else if (runParams.SubgroupsGridded)
        randGenInfo[0].subgroupsCenter = ScanChainRandom::GRIDDED;
    else
        randGenInfo[0].subgroupsCenter = ScanChainRandom::SAME_SCHEME;
    
    
    randGenInfo[1].numberLooseCells = runParams.GroupLooseCells;
    randGenInfo[1].numberSubGroups  = 0;
    randGenInfo[1].desiredNumberPOElts = runParams.GroupPOs;
    if (runParams.GroupUniformCellDist)
        randGenInfo[1].randCellType = 2;  // use uniform distribution
    else
        randGenInfo[1].randCellType = 1;  // use normal distribution
    randGenInfo[1].xCenter = 0;
    randGenInfo[1].yCenter = 0;
    randGenInfo[1].sigma = runParams.GroupSigma;
    randGenInfo[1].offset = runParams.Offset; // same for subgroups as for
    // top-level group
    
    randGenInfo[1].externalInOut = false;
    randGenInfo[1].subgroupsCenter = ScanChainRandom::SAME_SCHEME;
}

//*****************************************************************************
void PrintOutHelp(RunParams &runParams)
{
    cout<<"Command line syntax:"<<endl;
    cout<<""<<endl;
    cout<< runParams.commandName.c_str();
    cout<<"    [-sg SubGroups] [-po PartialOrderConstraints]"<<endl;
    cout<<"    [-gpo SubgroupPartialOrderConstraints]"<<endl;
    cout<<"    [-sigma Sigma] [-gsigma SubgroupSigma]"<<endl;
    cout<<"    [-offset Offset]"<<endl;
    cout<<"    [-eio] [-ssc | -sgg] [-unif] [-gunif]"<<endl;
    cout<<"    [-h]"<<endl;
    cout<<""<<endl;
    cout<<""<<endl;
    cout<<"-lc LooseCells      The number of loose cells (i.e. not part of any"<<endl;
    cout<<"            subgroup) in the top-level group.  Default value is 40,"<<endl;
    cout<<"            but no importance should be ascribed to this fact."<<endl;
    cout<<""<<endl;
    cout<<"-glc SubgroupLooseCells Same, for subgroups of the top-level group."<<endl;
    cout<<"                        Default is 10, again not for any reason"<<endl;
    cout<<"                        likely to be important to the user."<<endl;
    cout<<""<<endl;
    cout<<"-sg Subgroups       Number of subgroups of the top-level group."<<endl;
    cout<<"                    Default (unimportant) is 2."<<endl;
    cout<<""<<endl;
    cout<<"-po PartialOrderConstraints"<<endl;
    cout<<"            The desired number of constraints in the transitive"<<endl;
    cout<<"            partial order.  The program may sometimes overshoot"<<endl;
    cout<<"            this value a little; it keeps adding constraints"<<endl;
    cout<<"            until there are at least this many in the transitive"<<endl;
    cout<<"            closure.  This value is for the top-level group.  Default"<<endl;
    cout<<"            (unimportant) is 10."<<endl;
    cout<<""<<endl;
    cout<<"-gpo SubgroupPartialOrderConstraints"<<endl;
    cout<<"            Same, for subgroups of the top-level group.  Default"<<endl;
    cout<<"            (unimportant) is 5."<<endl;
    cout<<""<<endl;
    cout<<"-sigma Sigma        Loose cells (and centers of subgroups) for the"<<endl;
    cout<<"            top-level group are chosen with X and Y coordinates"<<endl;
    cout<<"            independently normally distributed with mean 0"<<endl;
    cout<<"            and standard deviation Sigma.  Default (unimportant) is 5000."<<endl;
    cout<<""<<endl;
    cout<<"-gsigma SubgroupSigma   Same as above, for subgroups of the top-level group"<<endl;
    cout<<"            except that the mean X and Y coordinates are chosen"<<endl;
    cout<<"            in the same way as the coordinates of a loose cell"<<endl;
    cout<<"            in the top-level group.  Default (unimportant) is 1000."<<endl;
    cout<<""<<endl;
    cout<<"-offset Offset      The output pin of a loose cell will be chosen"<<endl;
    cout<<"            this many units to the right of its input pin.  Default"<<endl;
    cout<<"            (unimportant) is 10."<<endl;
    cout<<""<<endl;
    cout<<"-eio                 Force the input and output cell of the top-"<<endl;
    cout<<"                     level group to be near the outer edge"<<endl;
    cout<<""<<endl;
    cout<<"-ssc                 Cause subgroups to be centered at the same point"<<endl;
    cout<<"                     as the top-level group.  Takes precedence over"<<endl;
    cout<<"                     -sgg option"<<endl;
    cout<<""<<endl;
    cout<<"-sgg                 Lay out subgroups in a grid pattern"<<endl;
    cout<<""<<endl;
    cout<<"-unif                Use a uniform distribution (rather than normal) for"<<endl;
    cout<<"                     the cities in the top-level group"<<endl;
    cout<<""<<endl;
    cout<<"-gunif               Same as -unif, but for subgroups"<<endl;
    cout<<""<<endl;
    cout<<"-h                   Print this help message and exit"<<endl;
}


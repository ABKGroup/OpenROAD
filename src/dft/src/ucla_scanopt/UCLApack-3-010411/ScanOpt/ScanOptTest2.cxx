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




//This source file is for the original application (as in strip2.tar)
//that could either parse problems or generate random ones, and then
//would do a number of simulations, to get a statistical universe for
//performance analysis.
#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include <ScanChain/scanChainRandom.h>
#include <ScanChain/scanChainParsed.h>
#include <ScanOpt/optimizer1.h>

using namespace abkscanopt;
struct RunParams;
void parseArgs(int argc,char **argv,RunParams &runParams);
void SetRandGenParms(ScanChainRandom::RandGenInfo RandGenInfo[],
                     const RunParams &runParams);
void SetOptimParms(Optimizer::Params OptInfo[],
                   const RunParams &runParams);
void WriteParms(FILE *fpParms,const RunParams &runParams);
void PrintOutHelp(RunParams &runParams);


struct RunParams
{
    std::string commandName;
    bool PrintHelp;
    unsigned RandSeed;
    int LooseCells;
    int GroupLooseCells;
    int Subgroups;
    int POs;
    int GroupPOs;
    double Sigma;
    double GroupSigma;
    double Offset;
    int MajorLoops;
    int Descents;
    int GroupDescents;
    int KickMove;
    int GroupKickMove;
    int PrintFreq;
    bool Only2Opt;
    bool GroupOnly2Opt;
    bool ProblemPrint;
    bool bReadFromFile;
    bool ZeroTemp;
    bool GroupZeroTemp;
    char szInFileName[255];
    char szInFileFullPath[512];
    int  NumSimulations;
    bool ExternalInOut;
    bool SubgroupsSameCenter;
    bool SubgroupsGridded;
    bool UniformCellDist;
    bool GroupUniformCellDist;
    bool Euclidean;
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
        MajorLoops(100),
        Descents(-1),
        GroupDescents(-1),
        KickMove(-1),
        GroupKickMove(-1),
        PrintFreq(1),
        Only2Opt(false),
        GroupOnly2Opt(false),
        ProblemPrint ( false),
        bReadFromFile ( false),
        ZeroTemp(true),
        GroupZeroTemp(true),
        NumSimulations(1),
        ExternalInOut ( false),
        SubgroupsSameCenter ( false),
        SubgroupsGridded ( false),
        UniformCellDist ( false),
        GroupUniformCellDist ( false),
        Euclidean ( false) {}
    
};

#ifdef COUNT
extern long nTwoOptCount;
extern long nThreeOptCount;
extern long nFindTwoOptCount;
extern long nFindThreeOptCount;
extern long nRunTwoOptCount;
extern long nRunThreeOptCount;
#endif

//********************************************************************
int main(int argc,char **argv)
{
    int SimulationIndex;
#ifdef _EUCLIDEAN
    char Name[] = "SiE             ";
    char TourName[] = "touE                                    ";
    char GraphName[] = "xE                                     ";
    char StatsName[] = "StatE                                  ";
    char ProblemName[] = "ProbleE                              ";
    char NumString[20];
    char DotText[]=".txt";
#else
    char Name[] = "Sim             ";
    char TourName[] = "tour                                    ";
    char GraphName[] = "xg                                     ";
    char StatsName[] = "Stats                                  ";
    char ProblemName[] = "Problem                              ";
    char NumString[20];
    char DotText[]=".txt";
#endif
    
    
    
    
    
    
    ScanChainRandom::RandGenInfo RandGenInfo[2];
    Optimizer::Params OptInfo[3];
    //long                        PathSeed[2];
    
    RunParams runParams;
    parseArgs(argc,argv,runParams);
    SeedHandler::overrideExternalSeed(runParams.RandSeed);
    SeedHandler::turnOffLogging();
    SetRandGenParms(RandGenInfo,runParams);
    SetOptimParms(OptInfo,runParams);
    RandomRawUnsigned randunsSim("main,randunsSim");
    //RandomRawUnsigned randunsFile("main,randunsFile");
    
    
    ScanChain *pProblem; //the reason we use a pointer here
                         //is so we can change the problem for
                         //each simulation if using random generation,
                         //but use the same one if we're reading from a file

    if (runParams.PrintHelp)
    {
        PrintOutHelp(runParams);
        return 0;
    }
    
    if (runParams.bReadFromFile)
    {
        ::ifstream is(runParams.szInFileName);
        pProblem = new ScanChainParsed(is);

        bool hasOpt=pProblem->hasOptimizableSubgroups();
        
        if (runParams.KickMove==-1)
            runParams.KickMove = (hasOpt) ?
            15:12;
        
        if (runParams.GroupKickMove==-1 && hasOpt)
            runParams.GroupKickMove = 15;
        
        if (runParams.Descents==-1)
            runParams.Descents = (hasOpt) ?
            5:15;
        
        if (runParams.GroupDescents==-1 && hasOpt)
            runParams.GroupDescents=15;
        
        SetOptimParms(OptInfo,runParams);
        
        
    }
    

    FILE *fpParms = fopen("Parms.txt","w");
    WriteParms(fpParms,runParams);
    fclose(fpParms);
    
    for (SimulationIndex=0;SimulationIndex<runParams.NumSimulations;SimulationIndex++)
    {
        //RandGenInfo[0].RandSeed = RandSeed+SimulationIndex;
        //RandGenInfo[1].RandSeed = 0;
        sprintf(Name+3,"%2d",SimulationIndex);
        
        TourName[4]= '\0';
        GraphName[2] = '\0';
        StatsName[5] = '\0';
        ProblemName[7] = '\0';
        
        sprintf(NumString,"%03d",SimulationIndex);
        strcat(NumString,DotText);
        strcat(TourName,NumString);
        strcat(GraphName,NumString);
        strcat(StatsName,NumString);
        strcat(ProblemName,NumString);
        
        if (!runParams.bReadFromFile)
        {
            pProblem = new ScanChainRandom(RandGenInfo,Name);
            bool hasOpt=pProblem->hasOptimizableSubgroups();
            if (runParams.KickMove==-1)
                runParams.KickMove = (hasOpt) ?
                15:12;
            
            if (runParams.GroupKickMove==-1 && hasOpt)
                runParams.GroupKickMove = 15;
            
            if (runParams.Descents==-1)
                runParams.Descents = (hasOpt) ?
                5:15;
            
            if (runParams.GroupDescents==-1 && hasOpt)
                runParams.GroupDescents=15;
            
            SetOptimParms(OptInfo,runParams);
        }
        
        if (runParams.ProblemPrint && (!runParams.bReadFromFile ||
            !SimulationIndex))
        {
            ::ofstream os(ProblemName);
            pProblem->printTopDown(os,0);
        }
        
        Optimizer::Params optParams=OptInfo[0];
        optParams.subParams=&OptInfo[1];
        
        Timer totalTime;
        Optimizer1 opt(*pProblem,randunsSim,optParams);
        totalTime.stop();
        cout << "Total time to optimize: " << totalTime << endl;

        ::ofstream fStats(StatsName);
        opt.rawStatsOut(fStats);

        ::ofstream fTour(TourName);
        pProblem->printPathTopDown(fTour,0);

        ::ofstream fGraph(GraphName);
        pProblem->printPathXgraph(fGraph);
        
        
        if (!runParams.bReadFromFile)
        {
            delete pProblem;pProblem=NULL;
        }
    }

    if (pProblem) delete pProblem;
    
#ifdef COUNT
    printf("scoRun2Opt called %ld times\n",nRunTwoOptCount);
    printf("scoRun3Opt called %ld times\n",nRunThreeOptCount);
    printf("scoFind2Opt called %ld times\n",nFindTwoOptCount);
    printf("scoFind3Opt called %ld times\n",nFindThreeOptCount);
    printf("actually performed %ld 2-opts\n",nTwoOptCount);
    printf("actually performed %ld 3-opts\n",nThreeOptCount);
#endif
    
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
            sscanf(argv[i],"%d",&runParams.RandSeed);
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
        else if (strcmp("-ml",argv[i]) == 0)  
        {
            i++;
            sscanf(argv[i],"%d",&runParams.MajorLoops);
        }
        else if (strcmp("-des",argv[i]) == 0) 
        {
            i++;
            sscanf(argv[i],"%d",&runParams.Descents);
        }
        else if (strcmp("-gdes",argv[i]) == 0)  
        {
            i++;
            sscanf(argv[i],"%d",&runParams.GroupDescents);
        }
        else if (strcmp("-km",argv[i]) == 0)  
        {
            i++;
            sscanf(argv[i],"%d",&runParams.KickMove);
        }
        else if (strcmp("-gkm",argv[i]) == 0)  
        {
            i++;
            sscanf(argv[i],"%d",&runParams.GroupKickMove);
        }
        else if (strcmp("-pfreq",argv[i]) == 0)  
        {
            i++;
            sscanf(argv[i],"%d",&runParams.PrintFreq);
        }
        else if (strcmp("-sims",argv[i]) == 0)  
        {
            i++;
            sscanf(argv[i],"%d",&runParams.NumSimulations);
        }
        else if (strcmp("-only2",argv[i]) == 0) 
        {
            runParams.Only2Opt = true;
        }
        else if (strcmp("-gonly2",argv[i]) == 0)  
        {
            runParams.GroupOnly2Opt = true;
        }
        else if (strcmp("-probprnt",argv[i]) == 0)  
        {
            runParams.ProblemPrint = true;
        }
        else if (strcmp("-rfile",argv[i]) == 0)
        {
            runParams.bReadFromFile = true;
            strcpy(runParams.szInFileName,argv[++i]);
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
        else if (strcmp("-zt",argv[i]) == 0)
        {
            runParams.ZeroTemp = true;
        }
        else if (strcmp("-gzt",argv[i]) == 0)
        {
            runParams.GroupZeroTemp = true;
        }
        else if (strcmp("-euclid",argv[i]) == 0)
        {
            runParams.Euclidean = true;
        }
        else 
        {
            char errtxt[1023];
            sprintf(errtxt,"command argument %d %s unrecognized; program exiting\n",
                i,argv[i]);
            abkfatal(false,errtxt);
        }
    }
#ifdef _EUCLIDEAN
    if (!runParams.Euclidean)
    {
        abkfatal(false,"This code was compiled to use "
            "Euclidean distance; you must specify "
            "-euclid option\n");
    }
#else
    if (runParams.Euclidean)
    {
        abkfatal(false,"This code was compiled to use "
            "Manhattan distance; you may not specify "
            "-euclid option\n");
    }
#endif
}
//*********************************************************************
void SetRandGenParms(ScanChainRandom::RandGenInfo RandGenInfo[],
    const RunParams &runParams)
{
    RandGenInfo[0].numberLooseCells = runParams.LooseCells;
    RandGenInfo[0].numberSubGroups  = runParams.Subgroups;
    RandGenInfo[0].desiredNumberPOElts = runParams.POs;
    if (runParams.UniformCellDist)
        RandGenInfo[0].randCellType = 2;  // use uniform distribution
    else
        RandGenInfo[0].randCellType = 1;  // use normal distribution
    RandGenInfo[0].xCenter = 0;
    RandGenInfo[0].yCenter = 0;
    RandGenInfo[0].sigma = runParams.Sigma;
    RandGenInfo[0].offset = runParams.Offset;
    RandGenInfo[0].externalInOut = runParams.ExternalInOut;
    
    if (runParams.SubgroupsSameCenter)
        RandGenInfo[0].subgroupsCenter = ScanChainRandom::SAME_CENTER;
    else if (runParams.SubgroupsGridded)
        RandGenInfo[0].subgroupsCenter = ScanChainRandom::GRIDDED;
    else
        RandGenInfo[0].subgroupsCenter = ScanChainRandom::SAME_SCHEME;
    
    
    RandGenInfo[1].numberLooseCells = runParams.GroupLooseCells;
    RandGenInfo[1].numberSubGroups  = 0;
    RandGenInfo[1].desiredNumberPOElts = runParams.GroupPOs;
    if (runParams.GroupUniformCellDist)
        RandGenInfo[1].randCellType = 2;  // use uniform distribution
    else
        RandGenInfo[1].randCellType = 1;  // use normal distribution
    RandGenInfo[1].xCenter = 0;
    RandGenInfo[1].yCenter = 0;
    RandGenInfo[1].sigma = runParams.GroupSigma;
    RandGenInfo[1].offset = runParams.Offset; // same for subgroups as for
    // top-level group
    
    RandGenInfo[1].externalInOut = false;
    RandGenInfo[1].subgroupsCenter = ScanChainRandom::SAME_SCHEME;
}
//*********************************************************************
void SetOptimParms(Optimizer::Params OptInfo[],
    const RunParams &runParams)
{
    OptInfo[0].majorLoops = runParams.MajorLoops;
    OptInfo[0].zeroTemp = runParams.ZeroTemp;
    OptInfo[0].nDescents = runParams.Descents;
    OptInfo[0].kickMove = runParams.KickMove;
    OptInfo[0].only2Opt = runParams.Only2Opt;
    OptInfo[0].maintainStats = true;
    OptInfo[0].statsInterval = runParams.PrintFreq;
    OptInfo[0].subKind=Optimizer::OPTIMIZER1;
    
    OptInfo[1].majorLoops = 1;    // do more by
                                    // increasing
                                    // descents
                                    // instead
    OptInfo[1].nDescents = runParams.GroupDescents;
    OptInfo[1].kickMove = runParams.GroupKickMove;
    OptInfo[1].only2Opt = runParams.GroupOnly2Opt;
    
    //Next line is actually irrelevant as long as subgroups have
    //no subgroups and MajorLoops are 1 at subgroup level
    OptInfo[1].zeroTemp = runParams.GroupZeroTemp;
    OptInfo[1].maintainStats = false;
    OptInfo[1].statsInterval = INT_MAX;
    /*
    OptInfo[2].m_eOptimizer = OLDCODE; // this is a quick fix to
    // let subpaths of subgroups work
    OptInfo[2].m_eOptGroup = IOFIX;
    OptInfo[2].m_bMaintainStats = false;
    OptInfo[2].m_nStatsInterval = INT_MAX;
    OptInfo[2].m_bPrintStatsLocally = false;
    */
    
}
//*********************************************************************
void WriteParms(FILE *fpParms, const RunParams &runParams)
{
    if (!runParams.bReadFromFile)
    {
        fprintf(fpParms,"Problem generation parameters:\n");
        fprintf(fpParms,"    Top level group:\n");
        fprintf(fpParms,"    Number Loose Cells = %d\n",runParams.LooseCells);
        fprintf(fpParms,"    Number Subgroups = %d\n",runParams.Subgroups);
        fprintf(fpParms,"    Desired transitive PO elements = %d\n",runParams.POs);
        fprintf(fpParms,"    Std Dev X and Y coords = %f\n",runParams.Sigma);
        if (runParams.ExternalInOut)
            fprintf(fpParms,"    In/Out cells constrained to be outside\n");
        
        if (runParams.UniformCellDist)
            fprintf(fpParms,"    Uniform cell distribution used\n");
        
        if (runParams.SubgroupsSameCenter)
            fprintf(fpParms,
            "    Subgroups have same center as top-level group\n");
        
        fprintf(fpParms,"    Offset from in pin to out pin =%f\n",runParams.Offset);
        fprintf(fpParms,"    Initial seed for random problem generation = %d\n",
            runParams.RandSeed);
        fprintf(fpParms,"\n");
        fprintf(fpParms,"    Subgroups:\n");
        fprintf(fpParms,"    Number Loose Cells = %d\n",runParams.GroupLooseCells);
        fprintf(fpParms,"    Number Subgroups = %d\n",0);
        fprintf(fpParms,"    Desired transitive PO elements = %d\n",runParams.GroupPOs);
        fprintf(fpParms,"    Std Dev X and Y coords = %f\n",runParams.GroupSigma);
        if (runParams.GroupUniformCellDist)
            fprintf(fpParms,"    Uniform cell distribution used\n");
        
        fprintf(fpParms,"\n");
    }
    else
    {
        fprintf(fpParms,"Input file is %s\n",runParams.szInFileFullPath);
        fprintf(fpParms,"\n");
    }
    fprintf(fpParms,"Optimization parameters:\n");
    fprintf(fpParms,"    Top-level group:\n");
    fprintf(fpParms,"    Number of major loops = %d\n",runParams.MajorLoops);
    fprintf(fpParms,"    Number of greedy descents per LSMC optimization ="
        " %d\n",runParams.Descents);
    fprintf(fpParms,"    Type of kickmove = %d\n",runParams.KickMove);
    if (runParams.Only2Opt)
        fprintf(fpParms,"    Greedy descent phase uses only 2-opts\n");
    else
        fprintf(fpParms,"    Greedy descent phase uses both 2-opts and"
        " 3-opts\n");
    
    if (runParams.ZeroTemp)
        fprintf(fpParms,"    Zero-temp annealing at pin-change level\n");
    else
        fprintf(fpParms,"    Infinite-temp annealing at pin-change level\n");
    
    
    fprintf(fpParms,"\n");
    fprintf(fpParms,"    Subgroups:\n");
    fprintf(fpParms,"    Number of major loops = %d\n",1);
    fprintf(fpParms,"    Number of greedy descents per LSMC optimization ="
        " %d\n",runParams.GroupDescents);
    fprintf(fpParms,"    Type of kickmove = %d\n",runParams.GroupKickMove);
    if (runParams.GroupOnly2Opt)
        fprintf(fpParms,"    Greedy descent phase uses only 2-opts\n");
    else
        fprintf(fpParms,"    Greedy descent phase uses both 2-opts and"
        " 3-opts\n");
#ifdef _EUCLIDEAN
    fprintf(fpParms,"EUCLIDEAN DISTANCE USED!!!\n");
#endif
}
//*****************************************************************************
void PrintOutHelp(RunParams &runParams)
{
    cout<<"Command line syntax:"<<endl;
    cout<<""<<endl;
    cout<< runParams.commandName.c_str() << " [-rfile inputfile] [-probprnt] [-ml MajorLoops]"<<endl;
cout<<"    [-des Descents] [-gdes SubgroupDescents]"<<endl;
cout<<"    [-km  KickMove] [-gkm SubgroupKickMove]"<<endl;
cout<<"    [-only2] [-gonly2]"<<endl;
cout<<"    [-pfreq PrintFrequency] [-rs RandSeed]"<<endl;
cout<<"    [-sims NumberSimulations]"<<endl;
cout<<""<<endl;
cout<<"    [-lc LooseCells] [-glc SubgroupLooseCells]"<<endl;
cout<<"    [-sg SubGroups] [-po PartialOrderConstraints]"<<endl;
cout<<"    [-gpo SubgroupPartialOrderConstraints]"<<endl;
cout<<"    [-sigma Sigma] [-gsigma SubgroupSigma]"<<endl;
cout<<"    [-offset Offset]"<<endl;
cout<<"    [-h]"<<endl;
cout<<""<<endl;
cout<<"Note:  The second block of command-line options is for use only if you"<<endl;
cout<<"want the program to generate the problem for you.  If you specify"<<endl;
cout<<"an input file, any options from the second block will simply be ignored."<<endl;
cout<<""<<endl;
cout<<"-rfile inputfile    Read in the problem from file inputfile"<<endl;
cout<<""<<endl;
cout<<"-probprnt       Print out the problem in standard format."<<endl;
cout<<"            If you read the problem in with -rfile,"<<endl;
cout<<"            it will be printed out to file Problem000.txt"<<endl;
cout<<"            along with machine-generated comments."<<endl;
cout<<"            If you had the software generate the problems,"<<endl;
cout<<"            they will be printed out (one per simulation)"<<endl;
cout<<"            to files called Problemnnn.txt where nnn is"<<endl;
cout<<"            the index of the simulation."<<endl;
cout<<""<<endl;
cout<<"-ml MajorLoops      This is the largest unit of measurement for"<<endl;
cout<<"            how long the simulation is allowed to run."<<endl;
cout<<"            In each major loop, the program chooses a subgroup"<<endl;
cout<<"            of the top-level group (i.e. of the whole"<<endl;
cout<<"            problem), changes its input/output cells,"<<endl;
cout<<"            and reoptimizes that group and the top-level"<<endl;
cout<<"            group.  Default: 100"<<endl;
cout<<""<<endl;
cout<<"-des Descents       Each time the top-level group is optimized,"<<endl;
cout<<"            the program will descend to a local minimum"<<endl;
cout<<"            in terms of 2-opts and 3-opts this many times,"<<endl;
cout<<"            doing kickmoves in between.  Default: 5 if"<<endl;
cout<<"            there are optimizable subgroups of the top-level"<<endl;
cout<<"            group, 15 otherwise."<<endl;
cout<<""<<endl;
cout<<"-gdes SubgroupDescents  Same thing, but for subgroups of the top-level"<<endl;
cout<<"            group.  Note that the data structures permit"<<endl;
cout<<"            subgroups to have subgroups themselves, but"<<endl;
cout<<"            the front end is not currently set up to"<<endl;
cout<<"            use them.  That should be an easy enhancement."<<endl;
cout<<"            The optimization code *is* set up to deal"<<endl;
cout<<"            with subgroups of subgroups (and so on) but"<<endl;
cout<<"            this feature is not yet tested.  Default: 15"<<endl;
cout<<""<<endl;
cout<<"-km KickMove        Specifies the type of kickmove to be used when"<<endl;
cout<<"            optimizing the top-level group.  The currently-"<<endl;
cout<<"            implemented kickmoves are:"<<endl;
cout<<""<<endl;
cout<<"            3   Do a random 3-opt."<<endl;
cout<<"            4   Do a random 4-opt."<<endl;
cout<<"            5   Do a random 5-opt."<<endl;
cout<<"            6   Do a random 6-opt."<<endl;
cout<<"            7   Do a random 7-opt."<<endl;
cout<<"            10  Do a random 10-opt."<<endl;
cout<<"            11  Do a random 2-opt."<<endl;
cout<<"            12  Do 2 random 2-opts."<<endl;
cout<<"            13  Do 3 random 2-opts."<<endl;
cout<<"            14  Do 4 random 2-opts."<<endl;
cout<<"            15  Do 5 random 2-opts."<<endl;
cout<<""<<endl;
cout<<"            Default: 15 if there are optimizable subgroups of"<<endl;
cout<<"            the top-level group, 12 otherwise."<<endl;
cout<<""<<endl;
cout<<"-gkm SubgroupKickMove   Specifies type of kickmove to be used when optimizing"<<endl;
cout<<"            subgroups of the top-level group.  Default: 15"<<endl;
cout<<""<<endl;
cout<<"-only2          Specifies that only 2-opts (i.e. not 3-opts)"<<endl;
cout<<"            shall be used during the greedy-descent"<<endl;
cout<<"            phase for the top-level group."<<endl;
cout<<""<<endl;
cout<<"-gonly2         Same thing, for subgroups of the top-level"<<endl;
cout<<"            group."<<endl;
cout<<""<<endl;
cout<<""<<endl;
cout<<"-pfreq PrintFrequency   Every PrintFrequency-many major loops, a"<<endl;
cout<<"            line will be added to the statistics giving"<<endl;
cout<<"            the cost of the best solution found so far,"<<endl;
cout<<"            the cpu time used, the count (which will always"<<endl;
cout<<"            be 1; it's there for other purposes), and the"<<endl;
cout<<"            number of major loops.  This will go into"<<endl;
cout<<"            a file called Statsnnn.txt where nnn is"<<endl;
cout<<"            the index of the simulation.  Default value is 1."<<endl;
cout<<""<<endl;
cout<<"-rs RandSeed    Seeds the random-number generator before"<<endl;
cout<<"            generating the intial random path and optimizing."<<endl;
cout<<"            Also controls the random number generator used"<<endl;
cout<<"            to generate random problems.  However those"<<endl;
cout<<"            random number generators are independent of"<<endl;
cout<<"            each other."<<endl;
cout<<"            Default: 12"<<endl;
cout<<""<<endl;
cout<<"-sims Simulations   Run this many simulations.  Default value is 1."<<endl;
cout<<""<<endl;
cout<<"-h                  Prints out help file"<<endl;
cout<<""<<endl;
cout<<"*********The following options are ignored if you specify an input file*****"<<endl;
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
cout<<"            is 10."<<endl;
cout<<""<<endl;
cout<<"-gpo SubgroupPartialOrderConstraints"<<endl;
cout<<""<<endl;
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
cout<<"Files:  This program produces output files called Statsnnn.txt, Parms.txt"<<endl;
cout<<"xgnnn.txt, tournnn.txt and (sometimes) Problemnnn.txt, where nnn varies from"<<endl;
cout<<"000 to 1 less than the number of simulations.  There is no checking"<<endl;
cout<<"as to whether these files already exist, so if you don't want them"<<endl;
cout<<"overwritten, run it in a different directory."<<endl;
cout<<""<<endl;
cout<<"The best solution found is written to tournnn.txt.  xgnnn.txt is an"<<endl;
cout<<"xgraph file showing this solution graphically.  Parms.txt shows the"<<endl;
cout<<"parameters used to run the program."<<endl;
cout<<""<<endl;
cout<<"Environment variables:  You can set an environment variable called"<<endl;
cout<<"ATSP_FILES.  If the file specified in the -rfile option is not found"<<endl;
cout<<"in the current working directory, the program will look in the"<<endl;
cout<<"directory given by ATSP_FILES."<<endl;
cout<<""<<endl;
cout<<"Bug :  If you use -rfile with an absolute pathname,"<<endl;
cout<<"the file name recorded in Parms.txt will look strange because it"<<endl;
cout<<"will have the current working directory prepended to the absolute"<<endl;
cout<<"pathname."<<endl;
cout<<""<<endl;
cout<<"********************Note on defaults******************************"<<endl;
cout<<""<<endl;
cout<<"The optimization defaults (i.e. number of descents, type of kickmove)"<<endl;
cout<<"have been chosen after a course of experimentation that was significant"<<endl;
cout<<"but by no means exhaustive.  Furthermore the parameters chosen as defaults"<<endl;
cout<<"were not necessarily the best in every case; rather they were ones that"<<endl;
cout<<"seemed to give good performance frequently.  It is entirely possible"<<endl;
cout<<"that some different combination of parameters will give better results"<<endl;
cout<<"for a particular problem."<<endl;
cout<<""<<endl;
cout<<"The defaults for problem generation (number of cells etc.) are included"<<endl;
cout<<"merely for completeness and are not intended to be taken as representative."<<endl;
cout<<"It is not anticipated that many users will use the program in this mode"<<endl;
cout<<"in any case; it was used primarily for the experimentation mentioned"<<endl;
cout<<"in the paragraph above."<<endl;
cout<<""<<endl;
cout<<"********************Note on timing*********************************"<<endl;
cout<<""<<endl;
cout<<"The length of time required to execute a given number of \"major loops\""<<endl;
cout<<"depends strongly on the number of cities involved, how they are"<<endl;
cout<<"divided into groups, and the number of partial-order constraints."<<endl;
cout<<"Some indicative times for 100 major loops follow, run on a 150-MHz"<<endl;
cout<<"Pentium machine:"<<endl;
cout<<""<<endl;
cout<<"400 cells, no subgroups, no partial-order constraints:  about 64 seconds"<<endl;
cout<<""<<endl;
cout<<"400 cells divided in top-level group + 4 subgroups, 10 partial-order"<<endl;
cout<<"constraints (in the transitive closure) in each group: 42 seconds."<<endl;
cout<<""<<endl;
cout<<"400 cells divided in top-level group + 9 subgroups, 10 partial-order"<<endl;
cout<<"constraints in each group:  19 seconds."<<endl;
cout<<""<<endl;
cout<<"1000 cells, no subgroups, no partial-order constraints: about 633 seconds"<<endl;
cout<<""<<endl;
cout<<"1000 cells divided in top-level group + 4 subgroups, no partial-order"<<endl;
cout<<"constraints:  about 142 seconds."<<endl;
cout<<""<<endl;
cout<<"1000 cells divided in top-level group + 4 subgroups, 25 partial-order"<<endl;
cout<<"constraints in each group: about 151 seconds."<<endl;
cout<<""<<endl;
cout<<""<<endl;
cout<<""<<endl;
cout<<"-------------------INPUT FILE FORMAT-------------------------------"<<endl;
cout<<""<<endl;
cout<<""<<endl;
cout<<"This describes the format of the input files to be used with the"<<endl;
cout<<"-rfile option."<<endl;
cout<<""<<endl;
cout<<"Machine-generated comments are enclosed in parentheses ()."<<endl;
cout<<"User comments are C++ style; // makes everything till the end of the line"<<endl;
cout<<"a comment."<<endl;
cout<<"Colons are ignored between tokens (not sure if I've tested this); they're"<<endl;
cout<<"for readability.  Keywords are not case-sensitive."<<endl;
cout<<""<<endl;
cout<<"For an example, see Problem000.txt in this directory.  It includes"<<endl;
cout<<"software-generated comments in parentheses."<<endl;
cout<<""<<endl;
cout<<"A group is declared as follows:"<<endl;
cout<<""<<endl;
cout<<"    GROUP \"group name\""<<endl;
cout<<""<<endl;
cout<<"        [Subgroups list]"<<endl;
cout<<"        [Loose cells list]"<<endl;
cout<<"        [Partial order constraints list]"<<endl;
cout<<"        [Legal inputs list]"<<endl;
cout<<"        [Legal outputs list]"<<endl;
cout<<""<<endl;
cout<<"    ENDGROUP"<<endl;
cout<<""<<endl;
cout<<"The subgroups list must come first, and the loose cells list second."<<endl;
cout<<""<<endl;
cout<<"[Partial order constraints list] has the following format:"<<endl;
cout<<""<<endl;
cout<<"    POs :"<<endl;
cout<<"    FROM fromcity TO tocity"<<endl;
cout<<"    FROM fromcity TO tocity"<<endl;
cout<<"    ..."<<endl;
cout<<"    ENDPOs"<<endl;
cout<<""<<endl;
cout<<"Here you may have 0 or more of the FROM ... TO constructions.  The"<<endl;
cout<<"indices for fromciy and tocity are zero-based and assume the subgroups"<<endl;
cout<<"come before the loose cells.  So for example if you have 5 subgroups"<<endl;
cout<<"for this group and you want to refer to the 2nd loose cell declared"<<endl;
cout<<"in the loose cells list, you would refer to it as 6."<<endl;
cout<<""<<endl;
cout<<"[Legal inputs list] has the following format:"<<endl;
cout<<""<<endl;
cout<<"    LEGALINs :"<<endl;
cout<<"    LEGALIN legalincity"<<endl;
cout<<"    LEGALIN legalincity"<<endl;
cout<<"    ..."<<endl;
cout<<"    ENDLEGALINs"<<endl;
cout<<""<<endl;
cout<<"Again the indices for legalincity are zero-based and assume the subgroups"<<endl;
cout<<"come before the loose cells."<<endl;
cout<<""<<endl;
cout<<"If you leave out the legal inputs list, the program will calculate for"<<endl;
cout<<"you all cities which are not excluded by partial order constraints from"<<endl;
cout<<"being legal input cities."<<endl;
cout<<""<<endl;
cout<<"[Legal outputs list] has the same format as the legal inputs list, except"<<endl;
cout<<"that the keywords are, respectively, LEGALOUTs, LEGALOUT, ENDLEGALOUTs."<<endl;
cout<<"Again you can omit the list and have the software calculate it maximally"<<endl;
cout<<"with respect to the partial order."<<endl;
cout<<""<<endl;
cout<<"[Subgroups list] has the following format:"<<endl;
cout<<""<<endl;
cout<<"    SUBGROUPs :"<<endl;
cout<<"    group"<<endl;
cout<<"    group"<<endl;
cout<<"    ..."<<endl;
cout<<"    ENDSUBGROUPs"<<endl;
cout<<""<<endl;
cout<<"where the format for group is the one described in this document (i.e."<<endl;
cout<<"it begins GROUP \"subgroup name\" and ends ENDGROUP)."<<endl;
cout<<""<<endl;
cout<<"[Loose cells list] has the following format:"<<endl;
cout<<""<<endl;
cout<<"    LOOSECELLs :"<<endl;
cout<<"    cell"<<endl;
cout<<"    cell"<<endl;
cout<<"    ..."<<endl;
cout<<"    ENDLOOSECELLs"<<endl;
cout<<""<<endl;
cout<<"where the format for cell is as follows:"<<endl;
cout<<""<<endl;
cout<<"    CELL \"cell name\""<<endl;
cout<<"    [ins list]"<<endl;
cout<<"    [outs list]"<<endl;
cout<<"    ENDCELL"<<endl;
cout<<""<<endl;
cout<<"[ins list] has the following format:"<<endl;
cout<<""<<endl;
cout<<"    INs :"<<endl;
cout<<"    IN xcoord ycoord"<<endl;
cout<<"    IN xcoord ycoord"<<endl;
cout<<"    ..."<<endl;
cout<<"    ENDINs"<<endl;
cout<<""<<endl;
cout<<"[outs list] has the same format except that the keywords are, respectively,"<<endl;
cout<<"OUTs, OUT and ENDOUTs."<<endl;
cout<<""<<endl;
cout<<"At this point it is important to note two things:"<<endl;
cout<<""<<endl;
cout<<"    1) the xcoord and ycoord may be specified as floating-point"<<endl;
cout<<"    values, but the current optimization code truncates them to"<<endl;
cout<<"    integers before using them."<<endl;
cout<<""<<endl;
cout<<"    2) while it is possible in the data structure to specify more than"<<endl;
cout<<"    one legal in-pin or out-pin for a cell, the current code only"<<endl;
cout<<"    considers the first coordinate-pair specified."<<endl;
cout<<""<<endl;
cout<<""<<endl;
cout<<"*****************Note on fully-specified subpaths*******************"<<endl;
cout<<""<<endl;
cout<<"If you know something about a desired solution, you can specify"<<endl;
cout<<"it as a subpath in the following somewhat brute-force manner:"<<endl;
cout<<"Make it a subgroup, then specify partial-order constraints that"<<endl;
cout<<"completely determine the path.  E.g. if you have three cells in"<<endl;
cout<<"the subpath, your POs list could look like"<<endl;
cout<<""<<endl;
cout<<"    POs :"<<endl;
cout<<"    FROM 2 TO 1"<<endl;
cout<<"    FROM 1 TO 0"<<endl;
cout<<"    ENDPOs"<<endl;
cout<<""<<endl;
cout<<"This would be placed inside a GROUP...ENDGROUP block as"<<endl;
cout<<"explained above.  Each cell (that is CELL...ENDCELL block) corresponding"<<endl;
cout<<"to these three cells would need to be moved inside the loose cells list"<<endl;
cout<<"(that is LOOSECELLs...ENDLOOSECELLs block) of this new subgroup."<<endl;
cout<<"The subgroups list should be empty (i.e. would consist of the"<<endl;
cout<<"keyword SUBGROUPs followed by ENDSUBGROUPs with no declarations in"<<endl;
cout<<"between) whereas the legal input and legal output lists should"<<endl;
cout<<"be omitted altogether (i.e. the keywords LEGALINs, LEGALOUTs,"<<endl;
cout<<"ENDLEGALINs and ENDLEGALOUTs would not appear in this three-cell"<<endl;
cout<<"subgroup).  So the whole declaration might look like:"<<endl;
cout<<""<<endl;
cout<<"    group \"subpath\""<<endl;
cout<<""<<endl;
cout<<""<<endl;
cout<<"    subgroups:"<<endl;
cout<<"    endsubgroups"<<endl;
cout<<"    loosecells:"<<endl;
cout<<"        cell \"subpath cell 0\""<<endl;
cout<<"            ins:"<<endl;
cout<<"                in 0 0"<<endl;
cout<<"            endins"<<endl;
cout<<"            outs:"<<endl;
cout<<"                out 0 10"<<endl;
cout<<"            endouts"<<endl;
cout<<"        endcell"<<endl;
cout<<"        cell \"subpath cell 1\""<<endl;
cout<<"            ins:"<<endl;
cout<<"                in 0 4000"<<endl;
cout<<"            endins"<<endl;
cout<<"            outs:"<<endl;
cout<<"                out 10 4000"<<endl;
cout<<"            endouts"<<endl;
cout<<"        endcell"<<endl;
cout<<"        cell \"subpath cell 2\""<<endl;
cout<<"            ins:"<<endl;
cout<<"                in 4000 0"<<endl;
cout<<"            endins"<<endl;
cout<<"            outs:"<<endl;
cout<<"                out 4000 10"<<endl;
cout<<"            endouts"<<endl;
cout<<"        endcell"<<endl;
cout<<"    endloosecells"<<endl;
cout<<"    pos:"<<endl;
cout<<"    from 2 to 1"<<endl;
cout<<"    from 1 to 0"<<endl;
cout<<"    endpos"<<endl;
cout<<""<<endl;
cout<<""<<endl;
cout<<"    endgroup"<<endl;
cout<<" "<<endl;
cout<<"This entire declaration would appear either in the subgroups list of "<<endl;
cout<<"the top-level group or that of one of its subgroups.  The declarations of"<<endl;
cout<<"the three cells, \"subpath cell 1\" etc., would have been moved from"<<endl;
cout<<"the loosecells list of that group."<<endl;
cout<<" "<<endl;
cout<<"Then you know that the final solution will have these cells in"<<endl;
cout<<"the order 2,1,0 with nothing in between."<<endl;
cout<<""<<endl;
cout<<"For this purpose you can make subgroups of subgroups.  The data"<<endl;
cout<<"format supports subgroups of subgroups, but the code is not"<<endl;
cout<<"yet tested for optimizing them.  However in this special case"<<endl;
cout<<"there is nothing to optimize."<<endl;
cout<<""<<endl;
cout<<"****************Note on large number of partial-order constraints*********"<<endl;
cout<<""<<endl;
cout<<"The code is tested primarily for the case that the partial-order"<<endl;
cout<<"constraints are rather sparse (not counting the 'subpath' case mentioned"<<endl;
cout<<"above).  If you put in a large number of constraints, performance"<<endl;
cout<<"will bog down.  In extreme cases you could cause an infinite loop."<<endl;
cout<<""<<endl;
cout<<"****************Note on legal input/output cells*************************"<<endl;
cout<<""<<endl;
cout<<"Ordinarily you will specify just one legal in cell and one legal out cell"<<endl;
cout<<"for the top-level group (i.e. the scan chain as a whole).  For the subgroups,"<<endl;
cout<<"you will ordinarily omit the LEGALINs and LEGALOUTs list, thereby allowing"<<endl;
cout<<"the program to choose any in/out cell for a subgroup that is consistent"<<endl;
cout<<"with the partial-order constraints."<<endl;
cout<<""<<endl;
cout<<"If you specify multiple LEGALINs and LEGALOUTs for the top-level group,"<<endl;
cout<<"the program will *not* attempt to choose the best one of each but will"<<endl;
cout<<"simply choose one of each at random.  This is a possible area for enhancement"<<endl;
cout<<"to the program."<<endl;
}


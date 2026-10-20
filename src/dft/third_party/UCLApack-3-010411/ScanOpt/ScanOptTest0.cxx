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




//This file is for the main ScanOpt optimization application, which
//parses a problem from an input file and optimizes it.

#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include <ScanChain/scanChainRandom.h>
#include <ScanChain/scanChainParsed.h>
#include <ScanOpt/optimizer1.h>

using namespace abkscanopt;
struct RunParams;
void parseArgs(int argc,char **argv,RunParams &runParams);
void SetOptimParms(Optimizer::Params OptInfo[],
                   const RunParams &runParams);
void PrintOutHelp(RunParams &runParams);


struct RunParams
{
    std::string commandName;
    bool PrintHelp;
    unsigned RandSeed;
    int MajorLoops;
    int Descents;
    int GroupDescents;
    int KickMove;
    int GroupKickMove;
    int PrintFreq;
    unsigned nnear;
    bool Only2Opt;
    bool GroupOnly2Opt;
    bool ZeroTemp;
    bool GroupZeroTemp;
    std::string InFileName;
    std::string StatsFile;
    std::string XgraphFile;
    std::string TourFile;
    RunParams():
    commandName(""),
    PrintHelp(false),
        RandSeed(12),
        MajorLoops(100),
        Descents(-1),
        GroupDescents(-1),
        KickMove(-1),
        GroupKickMove(-1),
        PrintFreq(1),
        nnear(20),
        Only2Opt(false),
        GroupOnly2Opt(false),
        ZeroTemp(true),
        GroupZeroTemp(true){}
    
};

//********************************************************************
int main(int argc,char **argv)
{
    
    Optimizer::Params OptInfo[3];
    
    RunParams runParams;
    parseArgs(argc,argv,runParams);
    SeedHandler::overrideExternalSeed(runParams.RandSeed);
    SeedHandler::turnOffLogging();
    SetOptimParms(OptInfo,runParams);
    RandomRawUnsigned randunsSim("main,randunsSim");
    //RandomRawUnsigned randunsFile("main,randunsFile");
    
    
    if (runParams.PrintHelp)
    {
        PrintOutHelp(runParams);
        return 0;
    }
    
    ::ifstream is(runParams.InFileName.c_str(),ios::in | ios::nocreate);
    abkfatal3(is.good(),"File ",runParams.InFileName.c_str()," not found");
    ScanChainParsed chain(is);

    bool hasOpt=chain.hasOptimizableSubgroups();
    
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

    Optimizer::Params optParams=OptInfo[0];
    optParams.subParams=&OptInfo[1];
    
    Timer totalTime;
    Optimizer1 opt(chain,randunsSim,optParams);
    totalTime.stop();
    cout << "Total time to optimize: " << totalTime << endl;
    
    if (!runParams.StatsFile.empty())
    {
        ::ofstream fStats(runParams.StatsFile.c_str());
        opt.rawStatsOut(fStats);
    }

    if (!runParams.TourFile.empty())
    {
        ::ofstream fTour(runParams.TourFile.c_str());
        chain.printPathTopDown(fTour,0);
    }

    if (!runParams.XgraphFile.empty())
    {
        ::ofstream fGraph(runParams.XgraphFile.c_str());
        chain.printPathXgraph(fGraph);
    }
    
    

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
        else if (strcmp("-nn",argv[i]) == 0)  
        {
            i++;
            sscanf(argv[i],"%d",&runParams.nnear);
        }
        else if (strcmp("-only2",argv[i]) == 0) 
        {
            runParams.Only2Opt = true;
        }
        else if (strcmp("-gonly2",argv[i]) == 0)  
        {
            runParams.GroupOnly2Opt = true;
        }
        else if (strcmp("-rfile",argv[i]) == 0)
        {
            runParams.InFileName=argv[++i];
        }
        else if (strcmp("-xfile",argv[i]) == 0)
        {
            runParams.XgraphFile=argv[++i];
        }
        else if (strcmp("-sfile",argv[i]) == 0)
        {
            runParams.StatsFile=argv[++i];
        }
        else if (strcmp("-tfile",argv[i]) == 0)
        {
            runParams.TourFile=argv[++i];
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
void SetOptimParms(Optimizer::Params OptInfo[],
    const RunParams &runParams)
{
    OptInfo[0].majorLoops = runParams.MajorLoops;
    OptInfo[0].zeroTemp = runParams.ZeroTemp;
    OptInfo[0].nDescents = runParams.Descents;
    OptInfo[0].kickMove = runParams.KickMove;
    OptInfo[0].nnear=runParams.nnear;
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
    OptInfo[1].nnear=runParams.nnear;
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
//*****************************************************************************
void PrintOutHelp(RunParams &runParams)
{
    cout<<"Command line syntax:"<<endl;
    cout<<""<<endl;
    cout<< runParams.commandName.c_str() << " -rfile inputfile  [-ml MajorLoops]"<<endl;
    cout<<"    [-des Descents] [-gdes SubgroupDescents]"<<endl;
    cout<<"    [-km  KickMove] [-gkm SubgroupKickMove]"<<endl;
    cout<<"    [-nn NumNear] [-only2] [-gonly2]"<<endl;
    cout<<"    [-pfreq PrintFrequency] [-rs RandSeed]"<<endl;
    cout<<"    [-tfile TourOutFile] [-sfile StatsOutFile] [-xfile XgraphOutFile]"<<endl;
    cout<<""<<endl;
    cout<<"    [-h]"<<endl;
    cout<<""<<endl;
    cout<<"-rfile inputfile    Read in the problem from file inputfile"<<endl;
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
    cout<<"-nn NumNear Specifies maximum number of nearest neighbors to consider"<<endl;
    cout<<"            when doing 2-opts and 3-opts"<<endl;
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
    cout<<"-tfile TourOutFile    If specified, the final optimized tour of"<<endl;
    cout<<"            the scan chain will be saved to this file."<<endl;
    cout<<""<<endl;
    cout<<"-sfile StatsOutFile   If specified, the cost achieved at various"<<endl;
    cout<<"            points in the optimization, and the times required, will"<<endl;
    cout<<"            be saved to this file."<<endl;
    cout<<""<<endl;
    cout<<"-xfile XgraphOutFile  If specified, the final optimized tour of the"<<endl;
    cout<<"             scan chain will be saved to this file in xgraph format."<<endl;
    cout<<""<<endl;
    cout<<"-h                  Prints out this help message"<<endl;
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


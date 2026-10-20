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




// scanChainParsed.cxx: implements parser for ScanChain
#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include "scanChainParsed.h"
using namespace abkscanopt;
const unsigned ScanChainParsed::bufLength=512;
//****************************************************************************
//************************************************************************
ScanChainParsed::ScanChainParsed(istream &is):ScanChain(),
m_error(ERROR_NONE),_in(is),
m_bPOsRead(false),
m_bLegalinsRead(false),
m_bLegaloutsRead(false),m_bSubgroupsRead(false),m_bLoosecellsRead(false),
_topLevelLeftoverChars(NULL),
m_pLastCharPos(_topLevelLeftoverChars),_topLevelCurrentLine(0),
m_nCurrentLine(_topLevelCurrentLine),
_topLevelLastChar('\0'),
m_cLastChar(_topLevelLastChar),_topLevelBuffer(new char[bufLength]),
m_szLineBuffer(_topLevelBuffer)
{
    ReadInGroup(true);
}

//************************************************************************

ScanChainParsed::ScanChainParsed(istream &is,
                                 char &LastChar,
                                 char *&pLeftoverChars,
                                 int  &CurrentLine,
                                 bufType &szLineBuffer)
:ScanChain(),m_error(ERROR_NONE),_in(is),m_bPOsRead(false),
m_bLegalinsRead(false),
m_bLegaloutsRead(false),m_bSubgroupsRead(false),m_bLoosecellsRead(false),
m_pLastCharPos(pLeftoverChars),m_nCurrentLine(CurrentLine),
m_cLastChar(LastChar),_topLevelBuffer(NULL),
m_szLineBuffer(szLineBuffer)
{
    ReadInGroup(false);
}
//************************************************************************
bool ScanChainParsed::ReadInGroup(bool bTopLevelGroup)
{
    char *pToken;
    TokenType token;
    
    if (!bTopLevelGroup)
    {
        string name;
        ReadQuotedName(name);
        _setName(name);
        m_state = STATE_READLISTS;
    }
    else
        m_state = STATE_INITIAL;
    
    // Parse the file
    //ResetScanner();
    while(!m_error && (pToken = NextToken()))
    {
        token = TranslateToken(pToken);
        switch (m_state)
        {
        case STATE_INITIAL: 
            {
                if (token==TOKEN_GROUP)
                {
                    string name;
                    ReadQuotedName(name);
                    _setName(name);
                    m_state = STATE_READLISTS;
                }
                else
                    m_error = ERROR_UNEXPECTED_TOKEN;
                break;
            }
        case STATE_READLISTS:
            {
                ReadListsHandler(token);
                break;
            } // case STATE_READLISTS
        case STATE_READPOS:
            {
                ReadInPOs(token);
                break;
            }
        case STATE_READLEGALINS:
            {
                ReadInLegalins(token);
                break;
            }
        case STATE_READLEGALOUTS:
            {
                ReadInLegalouts(token);
                break;
            }
        case STATE_READSUBGROUPS:
            {
                ReadInSubgroups(token);
                break;
            }
        case STATE_READLOOSECELLS:
            {
                ReadInLoosecells(token);
                break;
            }
            
        default:
            abkfatal(false,"Parse error");
        }   //switch (state)
        
        if (m_state == STATE_ENDGROUP)
            break;  //don't want to read another token
        
    }  // while
    
    if(!m_error && m_state != STATE_ENDGROUP)
        m_error = ERROR_INCOMPLETE;
    
    if (m_error)
        ErrorHandler();
    
    
    _checkIndexRanges();    
    _checkInOutLegality();
    
    _normalizePO();
    
    if (!m_bLegalinsRead)
        _chooseLegalIns();
    
    if (!m_bLegaloutsRead)
        _chooseLegalOuts();
    
    return true;
}
//************************************************************************
void ScanChainParsed::ReadListsHandler(TokenType token)
{
    switch (token)
    {
    case TOKEN_ENDGROUP:
        m_state = STATE_ENDGROUP;
        break;
    case TOKEN_POS:
        {
            if (m_bPOsRead)
                m_error = ERROR_DUPLICATE_POS;
            else
            {
                m_bPOsRead = true;
                m_state = STATE_READPOS;
            }
            break;
        }
    case TOKEN_LEGALINS:
        {
            if (m_bLegalinsRead)
                m_error = ERROR_DUPLICATE_LEGALINS;
            else
            {
                m_bLegalinsRead = true;
                m_state = STATE_READLEGALINS;
            }
            break;
        }
    case TOKEN_LEGALOUTS:
        {
            if (m_bLegaloutsRead)
                m_error = ERROR_DUPLICATE_LEGALOUTS;
            else
            {
                m_bLegaloutsRead = true;
                m_state = STATE_READLEGALOUTS;
            }
            break;
        }
    case TOKEN_SUBGROUPS:
        {
            if (m_bSubgroupsRead)
                m_error = ERROR_DUPLICATE_SUBGROUPS;
            else
            {
                m_bSubgroupsRead = true;
                m_state = STATE_READSUBGROUPS;
            }
            break;
        }
    case TOKEN_LOOSECELLS:
        {
            if (m_bLoosecellsRead)
                m_error = ERROR_DUPLICATE_LOOSECELLS;
            else
            {
                m_bLoosecellsRead = true;
                m_state = STATE_READLOOSECELLS;
            }
            break;
        }
    default:
        abkfatal(0,"Parse error");
    }
}
//************************************************************************
char *ScanChainParsed::NextToken()
{ //This code was stolen from Leonardo with minor modifications.
    char    *pToken,
        *pTemp;
    
    bool bInQuote = false;
    
    while(1)
    {
        if(m_pLastCharPos)
        {
            *m_pLastCharPos = m_cLastChar;
            pTemp = m_pLastCharPos;
        }
        else
        {
            m_szLineBuffer[0] = 0;
            _in.getline(m_szLineBuffer, bufLength);
            m_nCurrentLine++;
            
            if(!m_szLineBuffer[0] && _in.eof())
                return NULL;
            
            pTemp = m_szLineBuffer;
        }
        
        while(*pTemp == ' ' || *pTemp == '\t' || *pTemp == ':')
            pTemp++;
        
        if (*pTemp == '(')
        {
            while (*pTemp != ')' && *pTemp)
            {
                if (!*pTemp)
                {
                    m_error = ERROR_EOL_IN_COMMENT;
                    return NULL;
                }
                else
                    pTemp++;
            }
            pTemp++;
            m_cLastChar = *pTemp;
            
            if (m_cLastChar)
                m_pLastCharPos = pTemp;
            else
                m_pLastCharPos = NULL;
            continue;
        }
        
        if(!*pTemp)
        {
            m_pLastCharPos = NULL;
            continue;
        }
        
        if(*pTemp == '/' && *(pTemp + 1) == '/')
        {
            m_pLastCharPos = NULL; // the effect of this is to ignore
            // the rest of the line.
            continue;
        }
        
        
        if (*pTemp == '\"')
            bInQuote = true;
        
        
        pToken = pTemp;
        do
        {
            pTemp++;
        }
        while((*pTemp
            && *pTemp != ':'
            && *pTemp != '('
            && *pTemp != ')'
            && *pTemp != ' '
            && *pTemp != '\t'
            && *pTemp != '\"'
            && *pTemp != '/' && !bInQuote) 
            ||
            (*pTemp != '\"' && *pTemp && bInQuote));
        
        
        if (bInQuote)  
        {
            if (!pTemp)
            {
                m_error = ERROR_EOL_IN_QUOTE;
                return NULL;
            }
            pTemp++;    //want final quote mark included
        }
        
        m_cLastChar = *pTemp;
        m_pLastCharPos = pTemp;
        
        *pTemp = 0;
        
        
        return pToken;
    }
    
    // We'll never get here. This statement is just to satisfy the compiler.
    return NULL;
}
//************************************************************************

void ScanChainParsed::ReadQuotedName(string &Dest)
{
    char *pToken;
    int len;
    
    pToken = NextToken();
    
    if (*pToken != '\"')
    {
        m_error = ERROR_NAME_EXPECTED;
        return;
    }
    else
    {
        len = strlen(pToken);
        char *pc=new char[len-1];
        memcpy(pc,pToken+1,(len-2)*sizeof(char));
        pc[len-2] = '\0';
        Dest=pc;
        delete [] pc;
    }
    
    return;
}
//************************************************************************
ScanChainParsed::TokenType ScanChainParsed::TranslateToken(char *pToken)
{
    if (!strcasecmp(pToken,"GROUP"))
        return TOKEN_GROUP;
    else if (!strcasecmp(pToken,"ENDGROUP"))
        return TOKEN_ENDGROUP;
    else if (!strcasecmp(pToken,"POs"))
        return TOKEN_POS;
    else if (!strcasecmp(pToken,"ENDPOs"))
        return TOKEN_ENDPOS;
    else if (!strcasecmp(pToken,"LEGALINs"))
        return TOKEN_LEGALINS;
    else if (!strcasecmp(pToken,"ENDLEGALINs"))
        return TOKEN_ENDLEGALINS;
    else if (!strcasecmp(pToken,"LEGALOUTs"))
        return TOKEN_LEGALOUTS;
    else if (!strcasecmp(pToken,"ENDLEGALOUTs"))
        return TOKEN_ENDLEGALOUTS;
    else if (!strcasecmp(pToken,"SUBGROUPs"))
        return TOKEN_SUBGROUPS;
    else if (!strcasecmp(pToken,"ENDSUBGROUPs"))
        return TOKEN_ENDSUBGROUPS;
    else if (!strcasecmp(pToken,"LOOSECELLs"))
        return TOKEN_LOOSECELLS;
    else if (!strcasecmp(pToken,"ENDLOOSECELLs"))
        return TOKEN_ENDLOOSECELLS;
    else if (!strcasecmp(pToken,"FROM"))
        return TOKEN_FROM;
    else if (!strcasecmp(pToken,"TO"))
        return TOKEN_TO;
    else if (!strcasecmp(pToken,"LEGALIN"))
        return TOKEN_LEGALIN;
    else if (!strcasecmp(pToken,"LEGALOUT"))
        return TOKEN_LEGALOUT;
    else if (!strcasecmp(pToken,"CELL"))
        return TOKEN_CELL;
    else if (!strcasecmp(pToken,"ENDCELL"))
        return TOKEN_ENDCELL;
    else if (!strcasecmp(pToken,"INs"))
        return TOKEN_INS;
    else if (!strcasecmp(pToken,"ENDINs"))
        return TOKEN_ENDINS;
    else if (!strcasecmp(pToken,"OUTs"))
        return TOKEN_OUTS;
    else if (!strcasecmp(pToken,"ENDOUTs"))
        return TOKEN_ENDOUTS;
    else if (!strcasecmp(pToken,"IN"))
        return TOKEN_IN;
    else if (!strcasecmp(pToken,"OUT"))
        return TOKEN_OUT;
    
    return TOKEN_UNKNOWN;
}
//************************************************************************
void ScanChainParsed::ReadInPOs(TokenType token)
{
    char *pToken;
    bool done = false;
    unsigned From,To;
    
    m_state = STATE_EXPECT_FROM;
    
    do
    {
        switch (m_state)
        {
        case STATE_EXPECT_FROM:
            {
                switch (token)
                {
                case TOKEN_FROM:
                    {
                        pToken = NextToken();
                        From = atoi(pToken);
                        m_state = STATE_EXPECT_TO;
                        pToken = NextToken();
                        break;
                    }
                case TOKEN_ENDPOS:
                    {
                        done = true;
                        m_state = STATE_READLISTS;
                        break;
                    }
                default:
                    {
                        m_error = ERROR_FROM_EXPECTED;
                        break;
                    }
                }
                break;
            }
        case STATE_EXPECT_TO:
            {
                switch (token)
                {
                case TOKEN_TO:
                    {
                        pToken = NextToken();
                        To = atoi(pToken);
                        m_state = STATE_EXPECT_FROM;
                        pToken = NextToken();
                        if (!_transAddPOElt(From,To))
                        {
                            m_error = ERROR_CYCLE_IN_PO;
                        }
                        break;
                    }
                default:
                    {
                        m_error = ERROR_UNEXPECTED_TOKEN;
                        break;
                    }
                }
                break;
            }
        default:
            {
                m_error = ERROR_UNEXPECTED_STATE;
                return;
            }
        }
        token = TranslateToken(pToken);
    } while (!done && !m_error);
    
}
//************************************************************************
void ScanChainParsed::ReadInLegalins(TokenType token)
{
    char *pToken;
    bool done = false;
    int Legalin;
    
    do
    {
        switch (token)
        {
        case TOKEN_LEGALIN:
            {
                pToken = NextToken();
                Legalin = atoi(pToken);
                _addLegalin(Legalin);
                pToken = NextToken();
                break;
            }
        case TOKEN_ENDLEGALINS:
            {
                done = true;
                m_state = STATE_READLISTS;
                break;
            }
        default:
            {
                m_error = ERROR_LEGALIN_EXPECTED;
                break;
            }
        }
        token = TranslateToken(pToken);
    } while (!done && !m_error);
    
}
//************************************************************************
void ScanChainParsed::ReadInLegalouts(TokenType token)
{
    char *pToken;
    bool done = false;
    int Legalout;
    
    do
    {
        switch (token)
        {
        case TOKEN_LEGALOUT:
            {
                pToken = NextToken();
                Legalout = atoi(pToken);
                _addLegalout(Legalout);
                pToken = NextToken();
                break;
            }
        case TOKEN_ENDLEGALOUTS:
            {
                done = true;
                m_state = STATE_READLISTS;
                break;
            }
        default:
            {
                m_error = ERROR_LEGALOUT_EXPECTED;
                break;
            }
        }
        token = TranslateToken(pToken);
    } while (!done && !m_error);
    
}
//************************************************************************
void ScanChainParsed::ReadInSubgroups(TokenType token)
{
    char *pToken;
    bool done = false;
    m_bSubgroupsRead=true;
    string emptyString;


    if (m_bLoosecellsRead)
    {
        m_error=ERROR_LOOSECELLS_BEFORE_SUBGROUPS;
        return;
    }
    
    do
    {
        switch (token)
        {
        case TOKEN_GROUP:
            {
                ScanChain *pChain =
                    new ScanChainParsed(_in,m_cLastChar,m_pLastCharPos,
                    m_nCurrentLine,
                    m_szLineBuffer);

                _addSubgroup(pChain);

                pToken = NextToken();
                break;
            }
        case TOKEN_ENDSUBGROUPS:
            {
                done = true;
                m_state = STATE_READLISTS;
                break;
            }
        default:
            {
                m_error = ERROR_GROUP_EXPECTED;
                break;
            }
        }
        if (!done) token = TranslateToken(pToken);
    } while (!done && !m_error);
    
}
//************************************************************************
//************************************************************************
void ScanChainParsed::ReadInLoosecells(TokenType token)
{
    m_bLoosecellsRead=true;
    char *pToken;
    bool done = false;
    
    do
    {
        switch (token)
        {
        case TOKEN_CELL:
            {
                Cell cell;
                ReadQuotedName(cell._name);
                ReadInPorts(cell);
                _addCell(cell);
                pToken = NextToken();
                break;
            }
        case TOKEN_ENDLOOSECELLS:
            {
                done = true;
                m_state = STATE_READLISTS;
                break;
            }
        default:
            {
                m_error = ERROR_LOOSECELL_EXPECTED;
                break;
            }
        }
        token = TranslateToken(pToken);
    } while (!done && !m_error);
    
}
//************************************************************************
void ScanChainParsed::ReadInPorts(Cell &cell)
{
    char *pToken;
    TokenType token;
    bool done = false;
    //m_state = STATE_READINSOUTS;
    while(!m_error && (pToken = NextToken()))
    {
        token = TranslateToken(pToken);
        switch (token)
        {
        case TOKEN_INS:
            {
                ReadInIns(cell);
                break;
            }
        case TOKEN_OUTS:
            {
                ReadInOuts(cell);
                break;
            }
        case TOKEN_ENDCELL:
            {
                done = true;
                m_state = STATE_READLOOSECELLS;
                break;
            }
        default:
            {
                m_error = ERROR_UNEXPECTED_TOKEN;
            }
        }
        if (done)
            break;  //don't want to read in another token
    }// while
}
//************************************************************************
void ScanChainParsed::ReadInIns(Cell &cell)
{
    unsigned countIns=0;
    char *pToken;
    bool done = false;
    TokenType token;
    double Xcoord,Ycoord;
    pToken = NextToken();
    token = TranslateToken(pToken);
    do
    {
        switch (token)
        {
        case TOKEN_IN:
            {
                pToken = NextToken();
                Xcoord  = atof(pToken);
                pToken = NextToken();
                Ycoord = atof(pToken);
                cell.inX=static_cast<int>(Xcoord);
                cell.inY=static_cast<int>(Ycoord);
                countIns++;
                abkwarn(countIns!=2,
                    "Multiple INs for a cell, ignoring all but last)");
                pToken = NextToken();
                break;
            }
        case TOKEN_ENDINS:
            {
                done = true;
                break;
            }
        default:
            {
                m_error = ERROR_IN_EXPECTED;
                break;
            }
        }
        token = TranslateToken(pToken);
    } while (!done && !m_error);
    
}
//************************************************************************
void ScanChainParsed::ReadInOuts(Cell &cell)
{
    unsigned countOuts=0;
    char *pToken;
    bool done = false;
    TokenType token;
    double Xcoord,Ycoord;
    pToken = NextToken();
    token = TranslateToken(pToken);
    do
    {
        switch (token)
        {
        case TOKEN_OUT:
            {
                pToken = NextToken();
                Xcoord  = atof(pToken);
                pToken = NextToken();
                Ycoord = atof(pToken);
                cell.outX=static_cast<int>(Xcoord);
                cell.outY=static_cast<int>(Ycoord);
                countOuts++;
                abkwarn(countOuts!=2,
                    "Multiple OUTs for a cell, ignoring all but last)");
                pToken = NextToken();
                break;
            }
        case TOKEN_ENDOUTS:
            {
                done = true;
                break;
            }
        default:
            {
                m_error = ERROR_OUT_EXPECTED;
                break;
            }
        }
        token = TranslateToken(pToken);
    } while (!done && !m_error);
    
}
//************************************************************************
void ScanChainParsed::ErrorHandler()
{
    switch (m_error)
    {
    case ERROR_NONE:
        {
            return;
            break;
        }
    case ERROR_UNEXPECTED_TOKEN:
        {
            fprintf(stderr,"Unexpected token found\n");
            break;
        }
    case ERROR_DUPLICATE_POS:
        {
            fprintf(stderr,"More than one block of PO constraints\n");
            break;
        }
    case ERROR_DUPLICATE_LEGALINS:
        {
            fprintf(stderr,"More than one block of legal in cells\n");
            break;
        }
    case ERROR_DUPLICATE_LEGALOUTS:
        {
            fprintf(stderr,"More than one block of legal out cells\n");
            break;
        }
    case ERROR_DUPLICATE_SUBGROUPS:
        {
            fprintf(stderr,"More than one block of subgroups\n");
            break;
        }
    case ERROR_DUPLICATE_LOOSECELLS:
        {
            fprintf(stderr,"More than one block of loose cells\n");
            break;
        }
    case ERROR_EOL_IN_COMMENT:
        {
            fprintf(stderr,"Newline in comment\n");
            break;
        }
    case ERROR_EOL_IN_QUOTE:
        {
            fprintf(stderr,"Newline in quoted name\n");
            break;
        }
    case ERROR_NAME_EXPECTED:
        {
            fprintf(stderr,"Quoted name expected\n");
            break;
        }
    case ERROR_FROM_EXPECTED:
        {
            fprintf(stderr,"Expected token FROM or ENDPOs\n");
            break;
        }
    case ERROR_LEGALIN_EXPECTED:
        {
            fprintf(stderr,"Expected token LEGALIN or ENDLEGALINs\n");
            break;
        }
    case ERROR_LEGALOUT_EXPECTED:
        {
            fprintf(stderr,"Expected token LEGALOUT or ENDLEGALOUTs\n");
            break;
        }
    case ERROR_GROUP_EXPECTED:
        {
            fprintf(stderr,"Expected token GROUP or ENDSUBGROUPs\n");
            break;
        }
    case ERROR_LOOSECELL_EXPECTED:
        {
            fprintf(stderr,"Expected token LOOSECELL or ENDLOOSECELLs\n");
            break;
        }
    case ERROR_UNEXPECTED_STATE:
        {
            fprintf(stderr,"Internal error, unexpected state\n");
            break;
        }
    case ERROR_IN_EXPECTED:
        {
            fprintf(stderr,"Expected token IN or ENDINS\n");
            break;
        }
    case ERROR_OUT_EXPECTED:
        {
            fprintf(stderr,"Expected token OUT or ENDOUTS\n");
            break;
        }
    case ERROR_CYCLE_IN_PO:
        {
            fprintf(stderr,"Cycle detected in partial order\n");
            break;
        }
    case ERROR_INCOMPLETE:
        {
            fprintf(stderr,"Group declaration incomplete\n");
            break;
        }
    case ERROR_LOOSECELLS_BEFORE_SUBGROUPS:
        {
            fprintf(stderr,"Can''t parse SUBGROUPS after LOOSECELLS "
                "in a given group");
            break;
        }
    }
    char errtxt[1023];
    sprintf(errtxt,"Group is %s\n"
        "Line number %d reached in source file\n",
        getName().c_str(),
        m_nCurrentLine);
    abkfatal(false,errtxt);
    
}

ScanChainParsed::~ScanChainParsed()
{
    if (_topLevelBuffer) delete [] _topLevelBuffer;
}

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




//
// scanChainParsed.h: parser for ScanChain
//
#ifndef _SCANCHIANPARSED_H
#define _SCANCHIANPARSED_H
#include "scanChain.h"

namespace abkscanopt
{
typedef char *bufType;
class ScanChainParsed : public ScanChain
{
public:
    const static unsigned bufLength;
protected:
    enum TokenType
    {
        TOKEN_GROUP,
            TOKEN_ENDGROUP,
            TOKEN_POS,
            TOKEN_ENDPOS,
            TOKEN_LEGALINS,
            TOKEN_ENDLEGALINS,
            TOKEN_LEGALOUTS,
            TOKEN_ENDLEGALOUTS,
            TOKEN_SUBGROUPS,
            TOKEN_ENDSUBGROUPS,
            TOKEN_LOOSECELLS,
            TOKEN_ENDLOOSECELLS,
            TOKEN_FROM,
            TOKEN_TO,
            TOKEN_LEGALIN,
            TOKEN_LEGALOUT,
            TOKEN_CELL,
            TOKEN_ENDCELL,
            TOKEN_INS,
            TOKEN_ENDINS,
            TOKEN_OUTS,
            TOKEN_ENDOUTS,
            TOKEN_IN,
            TOKEN_OUT,
            TOKEN_UNKNOWN
    };
    enum StateEnum
    {
        STATE_INITIAL,
            STATE_ENDGROUP,
            STATE_READLISTS,
            STATE_READPOS,
            STATE_READLEGALINS,
            STATE_READLEGALOUTS,
            STATE_READSUBGROUPS,
            STATE_READLOOSECELLS,
            STATE_EXPECT_FROM,
            STATE_EXPECT_TO
            //STATE_READINSOUTS
    } m_state;
    
    enum ErrorEnum
    {
        ERROR_NONE = 0,
            ERROR_UNEXPECTED_TOKEN,
            ERROR_DUPLICATE_POS,
            ERROR_DUPLICATE_LEGALINS,
            ERROR_DUPLICATE_LEGALOUTS,
            ERROR_DUPLICATE_SUBGROUPS,
            ERROR_DUPLICATE_LOOSECELLS,
            ERROR_EOL_IN_COMMENT,
            ERROR_EOL_IN_QUOTE,
            ERROR_NAME_EXPECTED,
            ERROR_FROM_EXPECTED,
            ERROR_LEGALIN_EXPECTED,
            ERROR_LEGALOUT_EXPECTED,
            ERROR_GROUP_EXPECTED,
            ERROR_LOOSECELL_EXPECTED,
            ERROR_UNEXPECTED_STATE,
            ERROR_IN_EXPECTED,
            ERROR_OUT_EXPECTED,
            ERROR_CYCLE_IN_PO,
            ERROR_INCOMPLETE,
            ERROR_LOOSECELLS_BEFORE_SUBGROUPS
    } m_error;

    istream &_in;
    
    bool m_bPOsRead;
    bool m_bLegalinsRead;
    bool m_bLegaloutsRead;
    bool m_bSubgroupsRead;
    bool m_bLoosecellsRead;
    
    char *_topLevelLeftoverChars;
    char * &m_pLastCharPos;
    int _topLevelCurrentLine;
    int  &m_nCurrentLine;
    char _topLevelLastChar;
    char &m_cLastChar;
    bufType _topLevelBuffer; //NULL unless this is the top level
    bufType &m_szLineBuffer;

    bool ReadInGroup( bool bTopLevelGroup);
    
    //void ResetScanner();
    
    char *NextToken();
    
    void ReadQuotedName(string &Dest);
    
    TokenType TranslateToken(char*);
    
    void ReadListsHandler(TokenType token);
    
    void ReadInPOs(TokenType token);
    void ReadInLegalins(TokenType token);
    void ReadInLegalouts(TokenType token);
    void ReadInSubgroups(TokenType token);
    void ReadInLoosecells(TokenType token);
    void ReadInPorts(Cell &cell);
    void ReadInIns(Cell &cell);
    void ReadInOuts(Cell &cell);
    
    void ErrorHandler();

    //parsing ctor called hierarchically by the
    //simple parsing ctor, for subgroups.  Don't
    //call directly!
    ScanChainParsed(istream &is,char &LastChar,
                                 char *&pLeftoverChars,
                                 int  &CurrentLine,
                                 bufType &szLineBuffer);
    
public:
    ScanChainParsed(istream &is);
    virtual ~ScanChainParsed();
    
    
};
};
#endif //!defined(_SCANCHIANPARSED_H)

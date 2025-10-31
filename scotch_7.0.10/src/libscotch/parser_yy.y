%{
/* Copyright 2004,2007,2008,2011,2014,2018,2019,2021,2023-2025 IPB, Universite de Bordeaux, INRIA & CNRS
**
** This file is part of the Scotch software package for static mapping,
** graph partitioning and sparse matrix ordering.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
**
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
**
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
**
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
/************************************************************/
/**                                                        **/
/**   NAME       : parser_yy.y                             **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the syntactic parser     **/
/**                which processes strategy strings.       **/
/**                                                        **/
/**   DATES      : # Version 3.1  : from : 07 nov 1995     **/
/**                                 to   : 13 jun 1996     **/
/**                # Version 3.2  : from : 24 sep 1996     **/
/**                                 to   : 27 feb 1997     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to   : 01 oct 1998     **/
/**                # Version 4.0  : from : 20 dec 2001     **/
/**                                 to   : 11 jun 2004     **/
/**                # Version 5.1  : from : 30 oct 2007     **/
/**                                 to   : 24 jul 2011     **/
/**                # Version 6.0  : from : 30 sep 2014     **/
/**                                 to   : 27 apr 2018     **/
/**                # Version 7.0  : from : 02 mar 2018     **/
/**                                 to   : 11 sep 2024     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define SCOTCH_PARSER_YY

#include "module.h"
#include "common.h"
#include "parser.h"
#include "parser_yy.h"
#include "parser_ly.h"
#include "parser_lh.h"
#include "parser_ll.h"

#ifdef SCOTCH_DEBUG_PARSER3
#define YYDEBUG                     1
#endif /* SCOTCH_DEBUG_PARSER3 */

/*
**  The static and global definitions.
**  See also at the end of this file.
*/

/*+ Method token conversion array. +*/

extern unsigned int         parsermethtokentab[];  /* Pre-definition */

/*
**  The function prototypes.
*/

YY_DECL;                                          /* Definition of yylex(); Flex and Bison interaction is crap */
%}

%define api.pure full
%locations
%param {void * scanptr}
%parse-param {ParserEnv * penvptr}

%initial-action {
  @$.cobenum =
  @$.libenum =
  @$.coennum =
  @$.liennum = 1;
  @$.tebeptr =
  @$.teenptr = penvptr->textptr;
}

%union {
  char                      CASEVAL;              /* Case value          */
  StratTest *               TEST;                 /* Test type           */
  StratTestType             TESTOP;               /* Relational type     */
  double                    DOUBLE;               /* Double-precision    */
  INT                       INTEGER;              /* Integer             */
  char                      STRING[PARSERSTRINGLEN]; /* Character string */
  struct {
    const StratTab *        tabl;                 /* Current tables    */
    Strat *                 strat;                /* Current method    */
    const StratParamTab *   param;                /* Current parameter */
  } SAVE;                                         /* Parameter type    */
  Strat *                   STRAT;                /* Strategy tree     */
}

%token <STRING>      METHODNAME
%token <STRING>      PARAMNAME
%token <CASEVAL>     VALCASE
%token <DOUBLE>      VALDOUBLE
%token <INTEGER>     VALINT
%token <STRING>      VALSTRING
%token               VALSTRAT      VALPARAM      VALTEST

%type <TEST>         TEST          TESTOR        TESTAND       TESTNOT
%type <TEST>         TESTREL       TESTEXPR1     TESTEXPR2     TESTEXPR3
%type <TEST>         TESTEXPR4     TESTVAL       TESTVAR
%type <TESTOP>       TESTRELOP     TESTEXPR1OP   TESTEXPR2OP   TESTEXPR3OP
%type <DOUBLE>       VALSDOUBLE
%type <INTEGER>      VALSINT
%type <STRAT>        STRATCONCAT   STRATTEST     STRATTESTELSE STRATEMPTY
%type <STRAT>        STRATGROUP    STRATMETHOD   STRATSELECT

%start STRAT

%%

/*
**  These rules define the strategy grammar.
*/

STRAT         : STRATSELECT
              {
                penvptr->straptr = ($1);          /* Save pointer to root of tree */
              }
              ;

STRATSELECT   : STRATSELECT '|' STRATEMPTY
              {
                Strat *             straptr;

                if ((straptr = (Strat *) memAlloc (sizeof (Strat))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (1)");
                  stratExit  ($1);
                  stratExit  ($3);
                  YYABORT;
                }

                straptr->tablptr                 = penvptr->stratab;
                straptr->typeval                 = STRATNODESELECT;
                straptr->data.seledat.stratab[0] = ($1);
                straptr->data.seledat.stratab[1] = ($3);

                ($$) = straptr;
              }
              | STRATEMPTY
              ;

STRATEMPTY    : STRATCONCAT
              |
              {
                Strat *             straptr;

                if ((straptr = (Strat *) memAlloc (sizeof (Strat))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (2)");
                  YYABORT;
                }

                straptr->tablptr = penvptr->stratab;
                straptr->typeval = STRATNODEEMPTY;

                ($$) = straptr;
              }
              ;

STRATCONCAT   : STRATCONCAT STRATTEST
              {
                Strat *             straptr;

                if ((straptr = (Strat *) memAlloc (sizeof (Strat))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (3)");
                  stratExit  ($1);
                  stratExit  ($2);
                  YYABORT;
                }

                straptr->tablptr                 = penvptr->stratab;
                straptr->typeval                 = STRATNODECONCAT;
                straptr->data.concdat.stratab[0] = ($1);
                straptr->data.concdat.stratab[1] = ($2);

                ($$) = straptr;
              }
              | STRATTEST
              ;

STRATTEST     :
              {
                PARSERLLBEGIN (VALTEST);          /* Parse parameter tokens */
              }
                '/' TEST
              {
                PARSERLLBEGIN (VALSTRAT);         /* Parse strategy tokens */
              }
                '?' STRATSELECT STRATTESTELSE ';'
              {
                Strat *             straptr;

                if ((straptr = (Strat *) memAlloc (sizeof (Strat))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (4)");
                  stratExit  ($6);
                  if (($7) != NULL)
                    stratExit ($7);
                  stratTestExit ($3);
                  YYABORT;
                }

                straptr->tablptr                 = penvptr->stratab;
                straptr->typeval                 = STRATNODECOND;
                straptr->data.conddat.testptr    = ($3);
                straptr->data.conddat.stratab[0] = ($6);
                straptr->data.conddat.stratab[1] = ($7);

                ($$) = straptr;
              }
              | STRATGROUP
              ;

STRATTESTELSE : ':' STRATSELECT
              {
                ($$) = ($2);
              }
              |
              {
                ($$) = NULL;
              }
              ;

STRATGROUP    : '(' STRATSELECT ')'
              {
                ($$) = ($2);
              }
              | STRATMETHOD
              ;

STRATMETHOD   : METHODNAME
              {
                Strat *             straptr;
                int                 methnum;
                size_t              methlen;
                StratMethodTab *    methtab;
                int                 i;
                size_t              j;

                methnum = 0;
                methlen = 0;                      /* No method recognized yet     */
                methtab = penvptr->stratab->methtab; /* Point to the method table */
                for (i = 0; methtab[i].nameptr != NULL; i ++) {
                  if ((strncasecmp (($1),         /* Find longest matching code name */
                       methtab[i].nameptr,
                       j = strlen (methtab[i].nameptr)) == 0) &&
                      (j > methlen)) {
                    methnum = methtab[i].methnum;
                    methlen = j;
                  }
                }
                if (methlen == 0) {               /* If method name not known */
                  errorPrint ("stratParserParse: invalid method name \"%s\", line %d, column %d, at \"%s\"",
                              ($1), (@1).libenum, (@1).cobenum, (@1).tebeptr);
                  YYABORT;
                }
                if ((straptr = (Strat *) memAlloc (sizeof (Strat))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (5)");
                  YYABORT;
                }

                straptr->tablptr              = penvptr->stratab;
                straptr->typeval              = STRATNODEMETHOD;
                straptr->data.methdat.methnum = methnum; /* Set method type        */
                if (methtab[methnum].dataptr != NULL) /* If default values exist   */
                  memcpy (&straptr->data.methdat.datadat, /* Set values to default */
                          methtab[methnum].dataptr,
                          sizeof (StratNodeMethodData));

                penvptr->straptr = straptr;       /* Structure available for parameter processing */
              }
                METHODPARAM
              {
                StratParamTab *     paratab;
                int                 paraidx;

                paratab = penvptr->stratab->paratab; /* Point to the parameter table */
                for (paraidx = 0; paratab[paraidx].nameptr != NULL; paraidx ++) {
                  if ((paratab[paraidx].methnum == penvptr->straptr->data.methdat.methnum) && /* If a strategy parameter found for this method */
                      (paratab[paraidx].typeval == STRATPARAMSTRAT)) {
                    if (*((Strat **) ((byte *) &penvptr->straptr->data.methdat.datadat + /* And this parameter has not been set */
                        (paratab[paraidx].doffptr - paratab[paraidx].dbasptr))) == NULL)
                      errorPrintW ("stratParserParse: strategy parameter \"%s\" of method \"%s\" not set, line %d, column %d, before \"%s\"",
                                   paratab[paraidx].nameptr, penvptr->stratab->methtab[penvptr->straptr->data.methdat.methnum].nameptr,
                                   (@2).libenum, (@2).cobenum, (@2).tebeptr);
                  }
                }

                ($$) = penvptr->straptr;          /* Return current structure */
                penvptr->straptr = NULL;          /* No current structure     */
              }
              ;

METHODPARAM   :
              {
                PARSERLLBEGIN (VALPARAM);         /* Parse parameter tokens */
              }
                '{' PARAMLIST
              {
                PARSERLLBEGIN (VALSTRAT);         /* Parse strategy tokens */
              }
                '}'
              |                                   /* No parameters at all */
              ;

PARAMLIST     : PARAMLIST ',' PARAMPARAM
              | PARAMPARAM
              ;

PARAMPARAM    : PARAMNAME
              {
                int               paraidx;
                size_t            paralen;
                StratParamTab *   paratab;
                int               i;
                size_t            j;

                paraidx = 0;
                paralen = 0;                      /* No parameter recognized yet     */
                paratab = penvptr->stratab->paratab; /* Point to the parameter table */
                for (i = 0; paratab[i].nameptr != NULL; i ++) {
                  if ((paratab[i].methnum == penvptr->straptr->data.methdat.methnum) &&
                      (strncasecmp (($1),         /* Find longest matching parameter name */
                                    paratab[i].nameptr,
                                    j = strlen (paratab[i].nameptr)) == 0) &&
                      (j > paralen)) {
                    paraidx = i;
                    paralen = j;
                  }
                }
                if (paralen == 0) {
                  errorPrint ("stratParserParse: invalid method parameter name \"%s\", line %d, column %d, before \"%s\"",
                              ($1), (@1).libenum, (@1).cobenum, (@1).tebeptr);
                  YYABORT;
                }

                ($<SAVE>$).tabl = penvptr->stratab; /* Save current strategy tables   */
                penvptr->paraptr = &paratab[paraidx]; /* Save current parameter value */
                PARSERLLBEGIN (stratmethtokentab[penvptr->paraptr->typeval & ~STRATPARAMDEPRECATED]); /* Get non-deprecated type */
                if (penvptr->paraptr->typeval == STRATPARAMSTRAT) /* If parameter is a strategy                                  */
                  penvptr->stratab = (StratTab *) penvptr->paraptr->dselptr; /* Use new strategy tables                          */
              }
                '=' PARAMVAL
              {
                PARSERLLBEGIN (VALPARAM);         /* Go-on reading parameters          */
                penvptr->stratab = ($<SAVE>2).tabl; /* Restore current strategy tables */
              }
              ;

PARAMVAL      : VALCASE
              {
                char              c;              /* Character read             */
                char *            p;              /* Pointer to selector string */
                int               i;              /* Index in selector string   */

                if ((penvptr->paraptr->typeval & STRATPARAMDEPRECATED) == 0) { /* If parameter is not deprecated */
                  c = ($1);                       /* First, use char as is */
                  for (p = (char *) penvptr->paraptr->dselptr, i = 0;
                       (*p != '\0') && (*p != c);
                       p ++, i ++) ;
                  if (*p == '\0') {               /* Char was not found         */
                    c = tolower (c);              /* Convert char to lower case */
                    for (p = (char *) penvptr->paraptr->dselptr, i = 0;
                         (*p != '\0') && (*p != c);
                         p ++, i ++) ;
                    if (*p == '\0') {
                      errorPrint ("stratParserParse: invalid method parameter switch \"%s=%c\", line %d, column %d, before \"%s\"",
                                  penvptr->paraptr->nameptr, ($1), (@1).libenum, (@1).cobenum, (@1).tebeptr);
                      YYABORT;
                    }
                  }

#ifdef SCOTCH_DEBUG_PARSER2
                  if (((penvptr->paraptr->doffptr - penvptr->paraptr->dbasptr) + sizeof (int)) > sizeof (StratNodeMethodData)) {
                    errorPrint ("stratParserParse: internal error (1)");
                    YYABORT;
                  }
#endif /* SCOTCH_DEBUG_PARSER2 */

                  *((int *) ((byte *) &penvptr->straptr->data.methdat.datadat +
                             (penvptr->paraptr->doffptr - penvptr->paraptr->dbasptr))) = i;
                }
              }
              | VALSDOUBLE
              {
                if ((penvptr->paraptr->typeval & STRATPARAMDEPRECATED) == 0) { /* If parameter is not deprecated */
#ifdef SCOTCH_DEBUG_PARSER2
                  if (((penvptr->paraptr->doffptr - penvptr->paraptr->dbasptr) + sizeof (double)) > sizeof (StratNodeMethodData)) {
                    errorPrint ("stratParserParse: internal error (2)");
                    YYABORT;
                  }
#endif /* SCOTCH_DEBUG_PARSER2 */

                  *((double *) ((byte *) &penvptr->straptr->data.methdat.datadat +
                                (penvptr->paraptr->doffptr - penvptr->paraptr->dbasptr))) = ($1);
                }
              }
              | VALSINT
              {
                if ((penvptr->paraptr->typeval & STRATPARAMDEPRECATED) == 0) { /* If parameter is not deprecated */
#ifdef SCOTCH_DEBUG_PARSER2
                  if (((penvptr->paraptr->doffptr - penvptr->paraptr->dbasptr) + sizeof (INT)) > sizeof (StratNodeMethodData)) {
                    errorPrint ("stratParserParse: internal error (3)");
                    YYABORT;
                  }
#endif /* SCOTCH_DEBUG_PARSER2 */

                  *((INT *) ((byte *) &penvptr->straptr->data.methdat.datadat +
                             (penvptr->paraptr->doffptr - penvptr->paraptr->dbasptr))) = (INT) ($1);
                }
              }
              | VALSTRING
              {
                if ((penvptr->paraptr->typeval & STRATPARAMDEPRECATED) == 0) { /* If parameter is not deprecated */
#ifdef SCOTCH_DEBUG_PARSER2
                  if (((penvptr->paraptr->doffptr - penvptr->paraptr->dbasptr) + strlen ($1) + 1) > sizeof (StratNodeMethodData)) {
                    errorPrint ("stratParserParse: internal error (4)");
                    YYABORT;
                  }
#endif /* SCOTCH_DEBUG_PARSER2 */

                  strcpy ((char *) ((byte *) &penvptr->straptr->data.methdat.datadat +
                                    (penvptr->paraptr->doffptr - penvptr->paraptr->dbasptr)),
                          ($1));
                }
              }
              |
              {
                ($<SAVE>$).strat = penvptr->straptr;
                ($<SAVE>$).param = penvptr->paraptr;
                penvptr->straptr = NULL;
                penvptr->paraptr = NULL;
              }
                STRATSELECT
              {
                penvptr->straptr = ($<SAVE>1).strat; /* Restore current method    */
                penvptr->paraptr = ($<SAVE>1).param; /* Restore current parameter */

                if ((penvptr->paraptr->typeval & STRATPARAMDEPRECATED) == 0) { /* If parameter is not deprecated */
#ifdef SCOTCH_DEBUG_PARSER2
                  if (((penvptr->paraptr->doffptr - penvptr->paraptr->dbasptr) + sizeof (Strat *)) > sizeof (StratNodeMethodData)) {
                    errorPrint ("stratParserParse: internal error (5)");
                    YYABORT;
                  }
#endif /* SCOTCH_DEBUG_PARSER2 */

                  *((Strat **) ((byte *) &penvptr->straptr->data.methdat.datadat +
                                (penvptr->paraptr->doffptr - penvptr->paraptr->dbasptr))) = ($2);
                }
              }
              | error
              {
                errorPrint ("stratParserParse: invalid value for parameter \"%s\" of method \"%s\", line %d, column %d, before \"%s\"",
                            penvptr->paraptr->nameptr, penvptr->straptr->tablptr->methtab[penvptr->straptr->data.methdat.methnum].nameptr,
                            (@1).libenum, (@1).cobenum, (@1).tebeptr);
                YYABORT;
              }
              ;

TEST          : TESTOR
              ;

TESTOR        : TESTOR '|' TESTAND
              {
                StratTest *         testptr;

                if ((testptr = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (6)");
                  stratTestExit ($1);
                  stratTestExit ($3);
                  YYABORT;
                }

                testptr->testval         = STRATTESTOR;
                testptr->nodeval         = STRATPARAMLOG;
                testptr->data.testtab[0] = ($1);
                testptr->data.testtab[1] = ($3);

                ($$) = testptr;
              }
              | TESTAND
              ;

TESTAND       : TESTAND '&' TESTNOT
              {
                StratTest *         testptr;

                if ((testptr = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (7)");
                  stratTestExit ($1);
                  stratTestExit ($3);
                  YYABORT;
                }

                testptr->testval         = STRATTESTAND;
                testptr->nodeval         = STRATPARAMLOG;
                testptr->data.testtab[0] = ($1);
                testptr->data.testtab[1] = ($3);

                ($$) = testptr;
              }
              | TESTNOT
              ;

TESTNOT       : '!' TESTNOT
              {
                StratTest *         testptr;

                if ((testptr = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (8)");
                  stratTestExit ($2);
                  YYABORT;
                }

                testptr->testval         = STRATTESTNOT;
                testptr->nodeval         = STRATPARAMLOG;
                testptr->data.testtab[0] = ($2);

                ($$) = testptr;
              }
              | '(' TESTOR ')'
              {
                ($$) = ($2);
              }
              | TESTREL
              ;

TESTREL       : TESTEXPR1 TESTRELOP TESTEXPR1
              {
                StratTest *         testptr;

                if ((testptr = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (9)");
                  stratTestExit ($1);
                  stratTestExit ($3);
                  YYABORT;
                }
                testptr->testval         = ($2);
                testptr->nodeval         = STRATPARAMLOG;
                testptr->data.testtab[0] = ($1);
                testptr->data.testtab[1] = ($3);

                ($$) = testptr;
              }
              ;

TESTRELOP     : '<'
              {
                ($$) = STRATTESTLT;
              }
              | '='
              {
                ($$) = STRATTESTEQ;
              }
              | '>'
              {
                ($$) = STRATTESTGT;
              }
              ;

TESTEXPR1     : TESTEXPR1 TESTEXPR1OP TESTEXPR2
              {
                StratTest *         testptr;

                if ((testptr = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (10)");
                  stratTestExit ($1);
                  stratTestExit ($3);
                  YYABORT;
                }
                testptr->testval         = ($2);
                testptr->data.testtab[0] = ($1);
                testptr->data.testtab[1] = ($3);

                ($$) = testptr;
              }
              | TESTEXPR2
              ;

TESTEXPR1OP   : '+'
              {
                ($$) = STRATTESTADD;
              }
              | '-'
              {
                ($$) = STRATTESTSUB;
              }
              ;

TESTEXPR2     : TESTEXPR2 TESTEXPR2OP TESTEXPR3
              {
                StratTest *         testptr;

                if ((testptr = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  stratTestExit ($1);
                  stratTestExit ($3);
                  errorPrint    ("stratParserParse: out of memory (11)");
                  YYABORT;
                }
                testptr->testval         = ($2);
                testptr->data.testtab[0] = ($1);
                testptr->data.testtab[1] = ($3);

                ($$) = testptr;
              }
              | TESTEXPR3
              ;

TESTEXPR2OP   : '*'
              {
                ($$) = STRATTESTMUL;
              }
              ;

TESTEXPR3     : TESTEXPR3 TESTEXPR3OP TESTEXPR4
              {
                StratTest *         testptr;

                if ((testptr = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint    ("stratParserParse: out of memory (12)");
                  stratTestExit ($1);
                  stratTestExit ($3);
                  YYABORT;
                }
                testptr->testval         = ($2);
                testptr->data.testtab[0] = ($1);
                testptr->data.testtab[1] = ($3);

                ($$) = testptr;
              }
              | TESTEXPR4
              ;

TESTEXPR3OP   : '%'
              {
                ($$) = STRATTESTMOD;
              }
              ;

TESTEXPR4     : '(' TESTEXPR1 ')'
              {
                ($$) = ($2);
              }
              | TESTVAL
              | TESTVAR
              ;

TESTVAL       : VALSDOUBLE
              {
                StratTest *         testptr;

                if ((testptr = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (13)");
                  YYABORT;
                }

                testptr->testval         = STRATTESTVAL;
                testptr->nodeval         = STRATPARAMDOUBLE;
                testptr->data.val.valdbl = ($1);

                ($$) = testptr;
              }
              | VALSINT
              {
                StratTest *         testptr;

                if ((testptr = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (14)");
                  YYABORT;
                }

                testptr->testval         = STRATTESTVAL;
                testptr->nodeval         = STRATPARAMINT;
                testptr->data.val.valint = ($1);

                ($$) = testptr;
              }
              ;

TESTVAR       : PARAMNAME
              {
                StratTest *       testptr;
                StratParamTab *   condtab;
                int               para;
                size_t            paralen;
                int               i;
                size_t            j;

                para    = 0;
                paralen = 0;                      /* No parameter recognized yet */
                condtab = penvptr->stratab->condtab; /* Point to parameter table */
                for (i = 0; condtab[i].nameptr != NULL; i ++) {
                  if ((strncasecmp (($1),         /* Find longest matching parameter name */
                                    condtab[i].nameptr,
                                    j = strlen (condtab[i].nameptr)) == 0) &&
                      (j > paralen)) {
                    para    = i;
                    paralen = j;
                  }
                }
                if (paralen == 0) {
                  errorPrint ("stratParserParse: invalid graph parameter name \"%s\", line %d, column %d, before \"%s\"",
                              ($1), (@1).libenum, (@1).cobenum, (@1).tebeptr);
                  YYABORT;
                }

                if ((testptr = (StratTest *) memAlloc (sizeof (StratTest))) == NULL) {
                  errorPrint ("stratParserParse: out of memory (15)");
                  YYABORT;
                }

                testptr->testval          = STRATTESTVAR;
                testptr->nodeval          = condtab[para].typeval;
                testptr->data.var.datatab = penvptr->stratab;
                testptr->data.var.dataoft = condtab[para].doffptr - condtab[para].dbasptr;

                ($$) = testptr;
              }
              ;

VALSDOUBLE    : TESTEXPR1OP VALDOUBLE
              {
                ($$) = (($1) == STRATTESTSUB) ? - ($2) : ($2);
              }
              | VALDOUBLE
              ;

VALSINT       : TESTEXPR1OP VALINT
              {
                ($$) = (($1) == STRATTESTSUB) ? - ($2) : ($2);
              }
              | VALINT
              ;

%%

/*
**  The static and global definitions (bis).
**  These are put at the end of the file because
**  the token values that they use are not yet
**  defined in the first section of the file.
*/

unsigned int                stratmethtokentab[] = { /* Table for parameter/token type conversion */
                              VALCASE,
                              VALDOUBLE,
                              VALINT,
                              -1,                 /* No logical parameters */
                              VALSTRAT,
                              VALSTRING,
                              -1                  /* One more value to detect array overflow */
                            };

/************************************/
/*                                  */
/* These routines drive the parser. */
/*                                  */
/************************************/

/* This routine is the entry point for
** the strategy parser.
** It returns:
** - !NULL  : pointer to the strategy.
** - NULL   : on error.
*/

Strat *
stratParserParse (
const StratTab * const      stratab,              /*+ Pointer to parsing tables +*/
const char * const          textptr)              /*+ Strategy string to parse  +*/
{
  YY_BUFFER_STATE     buffdat;
  ParserEnv           penvdat;                    /* Parser environment    */
  yyscan_t            scandat;                    /* Pointer to lex memory */
  int                 o;

  penvdat.stratab = stratab;                      /* Point to the parsing tables             */
  penvdat.straptr = NULL;                         /* Clear up the temporary strategy pointer */
  penvdat.textptr = textptr;                      /* Initialize the lexical parser           */

  if (stratParserLexInit (&scandat) != 0) {
    errorPrint ("stratParserParse: cannot initialize reentrant parser");
    return     (NULL);
  }
  buffdat = stratParserLexScanString (textptr, scandat); /* Let's hope nothing breaks; error management in flex is just crap */
  stratParserLexBufSwitch (buffdat, scandat);

  o = yyparse (scandat, &penvdat);                /* Parse the strategy string */

  stratParserLexBufDelete (buffdat, scandat);
  stratParserLexDestroy   (scandat);

  if (o != 0) {
    if (penvdat.straptr != NULL)
      stratExit (penvdat.straptr);
    return (NULL);
  }

  return (penvdat.straptr);                       /* Return strategy pointer */
}

/* This routine displays the parser error message.
** It returns:
** - void  : in all cases.
*/

static
void
yyerror (
const ParserLocation * const  plocptr,            /*+ Scan location +*/
void * const                  scanptr,            /*+ Not used      +*/
const ParserEnv * const       penvptr,            /*+ Not used      +*/
const char * const            mesgptr)            /*+ Not used      +*/
{
  errorPrint ("stratParserParse: invalid strategy string, line %d, column %d, at \"%s\"",
              plocptr->libenum, plocptr->cobenum, plocptr->tebeptr);
}

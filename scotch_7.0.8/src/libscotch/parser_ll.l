%top{
/* Copyright 2004,2007,2018-2021,2023-2025 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : parser_ll.l                             **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Clement BARTHELEMY                      **/
/**                                                        **/
/**   FUNCTION   : This module is the lexical parser       **/
/**                which processes strategy strings.       **/
/**                                                        **/
/**   DATES      : # Version 3.1  : from : 07 nov 1995     **/
/**                                 to   : 23 aug 1996     **/
/**                # Version 3.2  : from : 24 sep 1996     **/
/**                                 to   : 05 jun 1997     **/
/**                # Version 3.3  : from : 01 oct 1998     **/
/**                                 to   : 11 sep 2001     **/
/**                # Version 4.0  : from : 20 dec 2001     **/
/**                                 to   : 23 dec 2001     **/
/**                # Version 5.1  : from : 09 jun 2009     **/
/**                                 to   : 24 jul 2011     **/
/**                # Version 7.0  : from : 02 mar 2018     **/
/**                                 to   : 07 jun 2025     **/
/**                                                        **/
/**   NOTES      : # In order for flex to read its input   **/
/**                  with getc() instead of fread, we set  **/
/**                  YY_ALWAYS_INTERACTIVE to 1. This may  **/
/**                  not always work with future releases. **/
/**                # Following common.h, we stop unistd.h  **/
/**                  from being included on Windows-like   **/
/**                  platforms by setting YY_NO_UNISTD_H.  **/
/**                                                        **/
/************************************************************/

#define YY_DECL                     int yylex (YYSTYPE * yylval_param, YYLTYPE * yylloc, yyscan_t yyscanner)
}

%{
/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"

#if (defined COMMON_OS_WINDOWS || defined HAVE_NOT_UNISTD_H)
#define YY_NO_UNISTD_H
#endif /* (defined COMMON_OS_WINDOWS || defined HAVE_NOT_UNISTD_H) */

#undef INTEGER                                    /* In case someone defined them */
#undef DOUBLE

#include "parser.h"
#include "parser_yy.h"
#include "parser_ly.h"                            /*+ Definitions produced by yacc +*/

#define YY_USER_ACTION              parserLocationUpdate (yylloc, yytext);
%}

%option bison-bridge
%option noinput
%option nounput
%option noyywrap
%option reentrant

IDENT                       [A-Za-z][0-9A-Za-z]*
INTEGER                     [0-9]+
FLOAT                       [0-9]+(\.[0-9]+)?([Ee][-+]?[0-9]+)?
STRING                      \"[^\"]*\"

%s lparam
%s lparamcase
%s lparamdouble
%s lparamint
%s lparamstring
%s ltest

%%
<INITIAL>[0-9A-Za-z]        {
                              strncpy (PARSERYYLVAL (STRING), yytext, PARSERSTRINGLEN);
                              PARSERYYLVAL (STRING[PARSERSTRINGLEN - 1]) = '\0';
                              return (METHODNAME);
                            }

<lparam>{IDENT}             {
                              strncpy (PARSERYYLVAL (STRING), yytext, PARSERSTRINGLEN);
                              PARSERYYLVAL (STRING[PARSERSTRINGLEN - 1]) = '\0';
                              return (PARAMNAME);
                            }

<lparamcase>[0-9A-Za-z]     {
                              PARSERYYLVAL (CASEVAL) = yytext[0];
                              return (VALCASE);
                            }

<lparamdouble>({FLOAT}|{INTEGER}) {
                              PARSERYYLVAL (DOUBLE) = atof (yytext);
                              return (VALDOUBLE);
                            }

<lparamint>{INTEGER}        {
                              PARSERYYLVAL (INTEGER) = (INT) atol (yytext);
                              return (VALINT);
                            }
<lparamint>{FLOAT}          {
                              PARSERYYLVAL (INTEGER) = (INT) atof (yytext); /* FLOAT is put after so that INTEGER can be matched */
                              return (VALINT);
                            }

<lparamstring>{STRING}      {
                              yytext[yyleng - 1] = '\0'; /* Remove the heading and trailing \" */
                              strncpy (PARSERYYLVAL (STRING), yytext + 1, PARSERSTRINGLEN);
                              PARSERYYLVAL (STRING[PARSERSTRINGLEN - 1]) = '\0';
                              return (VALSTRING);
                            }

<ltest>{INTEGER}            {
                              PARSERYYLVAL (INTEGER) = (INT) atol (yytext);
                              return (VALINT);
                            }
<ltest>{FLOAT}              {
                              PARSERYYLVAL (DOUBLE) = atof (yytext);
                              return (VALDOUBLE);
                            }
<ltest>{IDENT}              {
                              strncpy (PARSERYYLVAL (STRING), yytext, PARSERSTRINGLEN);
                              PARSERYYLVAL (STRING[PARSERSTRINGLEN - 1]) = '\0';
                              return (PARAMNAME);
                            }

[ \t\n]*                    {
                              parserLocationUpdate (yylloc, yytext); /* Update position without returning */
                            }
.                           return (yytext[0]);

%%
#include "parser_ll.h"

/*******************************************/
/*                                         */
/* These routines handle data input to the */
/* lexical analyzer.                       */
/*                                         */
/*******************************************/

/* This routine selects the sub-parser
** to parse the input.
** It returns:
** - VOID  : in all cases.
*/

void
stratParserLexBegin (
yyscan_t          yyscanner,
unsigned int      type)
{
  struct yyguts_t * yyg = (struct yyguts_t *) yyscanner; /* Awful hack to comply with the syntax of BEGIN */

  switch (type) {
    case VALCASE :
      BEGIN lparamcase;
      break;
    case VALDOUBLE :
      BEGIN lparamdouble;
      break;
    case VALINT :
      BEGIN lparamint;
      break;
    case VALSTRING :
      BEGIN lparamstring;
      break;
    case VALPARAM :
      BEGIN lparam;
      break;
    case VALSTRAT :
      BEGIN INITIAL;
      break;
    case VALTEST :
      BEGIN ltest;
      break;
  }
}

/* These routines are hooks for internal lexer
** routines which are called by the parser and
** cannot be properly accessed when Scotch suffixes
** are used.
*/

void
stratParserLexBufDelete (
YY_BUFFER_STATE             b,
yyscan_t                    yyscanner)
{
  yy_delete_buffer (b, yyscanner);
}

void
stratParserLexBufSwitch (
YY_BUFFER_STATE             new_buffer,
yyscan_t                    yyscanner)
{
  yy_switch_to_buffer (new_buffer, yyscanner);
}

int
stratParserLexDestroy (
yyscan_t                    yyscanner)
{
  return (yylex_destroy (yyscanner));
}

int
stratParserLexInit (
yyscan_t *                  ptr_yy_globals)
{
  return (yylex_init (ptr_yy_globals));
}

YY_BUFFER_STATE
stratParserLexScanString (
const char *                yystr,
yyscan_t                    yyscanner)
{
  return (yy_scan_string (yystr, yyscanner));
}

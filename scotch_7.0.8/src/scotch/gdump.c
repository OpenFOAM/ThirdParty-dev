/* Copyright 2019,2020,2023 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : gdump.c                                 **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This program creates a source code      **/
/**                that represents a source graph.         **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 23 nov 2019     **/
/**                                 to   : 22 jan 2020     **/
/**                # Version 7.0  : from : 21 jan 2023     **/
/**                                 to   : 21 jan 2023     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "gdump.h"

/*
**  The static definitions.
*/

static int                  C_fileNum = 0;        /* Number of file in arg list */
static File                 C_fileTab[C_FILENBR] = { /* The file array          */
                              { FILEMODER },
                              { FILEMODEW } };

static const char *         C_usageList[] = {
  "gdump [<input graph file> [<output code file>]] <options>",
  "  -b<base>    : Set base value",
  "  -h          : Display this help",
  "  -p<prefix>  : Use prefix name",
  "  -s<suffix>  : Use suffix name",
  "  -V          : Print program version and copyright",
  NULL };

/******************************/
/*                            */
/* This is the main function. */
/*                            */
/******************************/

int
main (
int                         argc,
char *                      argv[])
{
  SCOTCH_Graph        grafdat;                    /* Source graph      */
  SCOTCH_Num          baseval;                    /* Base value        */
  const char *        prefptr;                    /* Pointer to prefix */
  const char *        suffptr;                    /* Pointer to suffix */
  int                 i;

  errorProg ("gdump");

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    return     (EXIT_SUCCESS);
  }

  fileBlockInit (C_fileTab, C_FILENBR);           /* Set default stream pointers */
  baseval = -1;                                   /* Use graph base by default   */
  prefptr =
  suffptr = "";

  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes                        */
    if ((argv[i][0] != '-') || (argv[i][1] == '\0') || (argv[i][1] == '.')) { /* If found a file name */
      if (C_fileNum < C_FILEARGNBR)               /* File name has been given                         */
        fileBlockName (C_fileTab, C_fileNum ++) = argv[i];
      else
        errorPrint ("main: too many file names given");
    }
    else {                                        /* If found an option name */
      switch (argv[i][1]) {
        case 'H' :                                /* Give the usage message */
        case 'h' :
          usagePrint (stdout, C_usageList);
          return     (EXIT_SUCCESS);
        case 'P' :                                /* Prefix */
        case 'p' :
          prefptr = &argv[i][2];
          break;
        case 'S' :                                /* Suffix */
        case 's' :
          suffptr = &argv[i][2];
          break;
        case 'V' :
          fprintf (stderr, "gdump, version " SCOTCH_VERSION_STRING "\n");
          fprintf (stderr, SCOTCH_COPYRIGHT_STRING "\n");
          fprintf (stderr, SCOTCH_LICENSE_STRING "\n");
          return  (EXIT_SUCCESS);
        default :
          errorPrint ("main: unprocessed option '%s'", argv[i]);
      }
    }
  }

  fileBlockOpen (C_fileTab, C_FILENBR);           /* Open all files */

  SCOTCH_graphInit  (&grafdat);
  SCOTCH_graphLoad  (&grafdat, C_filepntrsrcinp, baseval, 0);
  SCOTCH_graphCheck (&grafdat);

  if (SCOTCH_graphDump (&grafdat, prefptr, suffptr, C_filepntrcodout) != 0)
    errorPrint ("main: cannot dump graph");

  SCOTCH_graphExit (&grafdat);

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (de)compression tasks */

  return (EXIT_SUCCESS);
}

/* Copyright 2004,2007,2008,2010-2012,2014,2018,2019,2023,2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : amk_m2.c                                **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Creates the distance map for            **/
/**                bidimensional mesh graphs, to be used   **/
/**                to build the architecture description   **/
/**                files for these graphs.                 **/
/**                                                        **/
/**   DATES      : # Version 1.3  : from : 21 apr 1994     **/
/**                                 to   : 21 apr 1994     **/
/**                # Version 2.0  : from : 12 jul 1994     **/
/**                                 to   : 12 nov 1994     **/
/**                # Version 3.0  : from : 17 jul 1995     **/
/**                                 to   : 19 sep 1995     **/
/**                # Version 3.2  : from : 31 may 1997     **/
/**                                 to   : 02 jun 1997     **/
/**                # Version 3.4  : from : 03 feb 2000     **/
/**                                 to   : 03 feb 2000     **/
/**                # Version 4.0  : from : 09 feb 2004     **/
/**                                 to   : 09 feb 2004     **/
/**                # Version 5.0  : from : 23 dec 2007     **/
/**                                 to   : 21 apr 2008     **/
/**                # Version 5.1  : from : 01 jul 2010     **/
/**                                 to   : 14 feb 2011     **/
/**                # Version 6.0  : from : 01 jan 2012     **/
/**                                 to   : 24 sep 2019     **/
/**                # Version 7.0  : from : 21 jan 2023     **/
/**                                 to   : 10 sep 2024     **/
/**                                                        **/
/**   NOTES      : # The vertices of the (dX,dY) mesh are  **/
/**                  numbered as terminals so that         **/
/**                  t (0,0)=0, t (1,0) = 1,               **/
/**                  t (dX - 1, 0) = dX - 1,               **/
/**                  t (0,1) = dX, and                     **/
/**                  t(x,y) = (y * dX) + x.                **/
/**                                                        **/
/**                # The nested dissection method should   **/
/**                  behave like the architecture built in **/
/**                  the mapper.                           **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "arch.h"
#include "arch_mesh.h"
#include "amk_m2.h"

/*
**  The static definitions.
*/

static int                  C_paraNum = 0;        /* Number of parameters       */
static int                  C_fileNum = 0;        /* Number of file in arg list */
static File                 C_fileTab[C_FILENBR] = { /* File array              */
                              { FILEMODEW } };

static const char *         C_usageList[] = {
  "amk_m2 <dimX> [<dimY> [<output target file>]] <options>",
  "  -h          : Display this help",
  "  -m<method>  : Decomposition method",
  "                  n  : Nested dissection (cut biggest dimension)",
  "                  o  : One-way dissection (y, then x)",
  "  -V          : Print program version and copyright",
  "",
  "Default option set is : '-Mn'",
  NULL };

/*************************************************/
/*                                               */
/* The main routine, which computes the distance */
/* triangular table.                             */
/*                                               */
/*************************************************/

int
main (
int                         argc,
char *                      argv[])
{
  ArchMesh2        arch;                          /* Mesh dimensions            */
  ArchMesh2Dom     dom;                           /* Initial domain             */
  C_MethType       typeval;                       /* Bipartitioning method      */
  Anum             termnbr;                       /* Number of terminal domains */
  Anum             termnum;
  Anum             termmax;                       /* Maximum terminal number    */
  Anum *           termtab;                       /* Terminal numbers table     */
  Anum             x0, y0, x1, y1;
  int              i;

  errorProg ("amk_m2");

  if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
    usagePrint (stdout, C_usageList);
    return     (EXIT_SUCCESS);
  }

  typeval   = C_METHNESTED;
  arch.c[0] =                                     /* Preset mesh dimensions */
  arch.c[1] = 1;

  fileBlockInit (C_fileTab, C_FILENBR);           /* Set default stream pointers */

  for (i = 1; i < argc; i ++) {                   /* Loop for all option codes                        */
    if ((argv[i][0] != '-') || (argv[i][1] == '\0') || (argv[i][1] == '.')) { /* If found a file name */
      if (C_paraNum < 2) {                        /* If number of parameters not reached              */
        if ((arch.c[C_paraNum ++] = atoi (argv[i])) < 1) /* Get the dimension                         */
          errorPrint ("main: invalid dimension '%s'", argv[i]);
        continue;                                 /* Process the other parameters */
      }
      if (C_fileNum < C_FILEARGNBR)               /* A file name has been given */
        fileBlockName (C_fileTab, C_fileNum ++) = argv[i];
      else
        errorPrint ("main: too many file names given");
    }
    else {                                        /* If found an option name */
      switch (argv[i][1]) {
        case 'M' :                                /* Use a built-in method */
        case 'm' :
          switch (argv[i][2]) {
            case 'N' :                            /* Nested dissection */
            case 'n' :
              typeval = C_METHNESTED;
              break;
            case 'O' :                            /* One-way dissection */
            case 'o' :
              typeval = C_METHONEWAY;
              break;
            default :
              errorPrint ("main: unprocessed option '%s'", argv[i]);
          }
          break;
        case 'H' :                               /* Give the usage message */
        case 'h' :
          usagePrint (stdout, C_usageList);
          return     (EXIT_SUCCESS);
        case 'V' :
          fprintf (stderr, "amk_m2, version " SCOTCH_VERSION_STRING "\n");
          fprintf (stderr, SCOTCH_COPYRIGHT_STRING "\n");
          fprintf (stderr, SCOTCH_LICENSE_STRING "\n");
          return  (EXIT_SUCCESS);
        default :
          errorPrint ("main: unprocessed option '%s'", argv[i]);
      }
    }
  }

  fileBlockOpen (C_fileTab, C_FILENBR);           /* Open all files */

  dom.c[0][0] = 0;                                /* Set the initial domain */
  dom.c[0][1] = arch.c[0] - 1;
  dom.c[1][0] = 0;
  dom.c[1][1] = arch.c[1] - 1;

  termnbr = arch.c[0] * arch.c[1];                /* Compute number of terminals                  */
  termmax = 0;                                    /* Maximum terminal value not known yet         */
  if ((termtab = (Anum *) memAlloc (termnbr * sizeof (Anum))) == NULL) /* Allocate terminal array */
    errorPrint ("main: out of memory");

  memset (termtab, -1, termnbr * sizeof (unsigned int)); /* Initilize mapping table */

  C_domnBipart (&arch, &dom, 1, termtab, &termmax, /* Compute terminal numbers */
                (typeval == C_METHNESTED) ? archMesh2DomBipart : archMesh2DomBipartO);

  fprintf (C_filepntrarcout, "deco\n0\n" ANUMSTRING "\t" ANUMSTRING "\n", /* Print file header */
           termnbr,
           termmax);
  for (termnum = 0; termnum < termnbr; termnum ++) /* For all terminals                   */
    fprintf (C_filepntrarcout, ANUMSTRING "\t1\t" ANUMSTRING "\n", /* Print terminal data */
             termnum, termtab[termnum]);

  for (y0 = 0; y0 < arch.c[1]; y0 ++) {           /* For all vertices */
    for (x0 = 0; x0 < arch.c[0]; x0 ++) {
      for (y1 = 0; y1 <= y0; y1 ++) {             /* Compute distance to smaller vertices */
        for (x1 = 0; (x1 < arch.c[0]) && ((y1 < y0) || (x1 < x0)); x1 ++)
          fprintf (C_filepntrarcout,
                   ((x1 == 0) && (y1 == 0)) ? ANUMSTRING : " " ANUMSTRING,
                   C_domnDist (x0, y0, x1, y1));
      }
      fprintf (C_filepntrarcout, "\n");
    }
  }

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (un)compression tasks */

  memFree (termtab);                              /* Free terminal number array */

  return (EXIT_SUCCESS);
}

/* This routine recursively determines the values
** of all the terminal vertices of the mesh domain,
** and puts them in table.
*/

void
C_domnBipart (
ArchMesh2 *                 archptr,
ArchMesh2Dom *              domptr,
Anum                        num,
Anum *                      termtab,
Anum *                      termmax,
DomnBipartFunc              funcptr)
{
  ArchMesh2Dom        dom0;
  ArchMesh2Dom        dom1;

  if (funcptr (archptr, domptr, &dom0, &dom1) == 0) { /* If we can bipartition                          */
    C_domnBipart (archptr, &dom0, num + num,     termtab, termmax, funcptr); /* Bipartition recursively */
    C_domnBipart (archptr, &dom1, num + num + 1, termtab, termmax, funcptr);
  }
  else {                                          /* If we have reached the end */
    termtab[domptr->c[1][0] * archptr->c[0] +     /* Set the terminal number    */
            domptr->c[0][0]] = num;
    if (*termmax < num)                           /* If we have reached a new maximum */
      *termmax = num;                             /* Record it                        */
  }
}

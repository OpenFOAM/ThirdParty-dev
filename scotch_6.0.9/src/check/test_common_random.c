/* Copyright 2012,2014,2016,2018 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_common_random.c                    **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the random number     **/
/**                generator module.                       **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 01 oct 2014     **/
/**                                 to     24 aug 2019     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE               600
#endif /* _XOPEN_SOURCE */
#ifndef __USE_XOPEN2K
#define __USE_XOPEN2K                             /* For POSIX pthread_barrier_t */
#endif /* __USE_XOPEN2K */

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../libscotch/module.h"
#include "../libscotch/common.h"
#include "scotch.h"

#define RANDNBR                     100

/*********************/
/*                   */
/* The main routine. */
/*                   */
/*********************/

int
main (
int                 argc,
char *              argv[])
{
  INT *               randtab;
  int                 randnum;
  FILE *              fileptr;
  int                 passnum;

  if (argc != 3) {
    SCOTCH_errorPrint ("usage: %s state_file passnum", argv[0]);
    exit (EXIT_FAILURE);
  }

  if ((randtab = malloc (RANDNBR * sizeof (INT))) == NULL) {
    SCOTCH_errorPrint ("main: out of memory");
    exit (EXIT_FAILURE);
  }

  intRandInit ();                                 /* Initialize random generator */

  for (randnum = 0; randnum < RANDNBR; randnum ++)
    randtab[randnum] = intRandVal (INTVALMAX);

  intRandReset ();

  passnum = (atoi (argv[2]) == 0);                /* First pass to write file; second pass to read it */

  if ((fileptr = fopen (argv[1], (passnum) ? "w+" : "r")) == NULL) {
    SCOTCH_errorPrint ("main: cannot open file");
    exit (EXIT_FAILURE);
  }

  if (passnum) {                                  /* If first pass */
    for (randnum = 0; randnum < RANDNBR; randnum ++) {
      if (randtab[randnum] != intRandVal (INTVALMAX)) {
        SCOTCH_errorPrint ("main: cannot replay random sequence");
        exit (EXIT_FAILURE);
      }
    }

    if (fwrite (randtab, sizeof (INT), RANDNBR, fileptr) < RANDNBR) {
      SCOTCH_errorPrint ("main: cannot write to file (1)");
      exit (EXIT_FAILURE);
    }

    switch (SCOTCH_randomSave (fileptr)) {        /* Try to save random state, if enabled */
      case 0 :
        if (fprintf (fileptr, "#") != 1) {        /* Write separator character */
          SCOTCH_errorPrint ("main: cannot write to file (2)");
          exit (EXIT_FAILURE);
        }

        for (randnum = 0; randnum < RANDNBR; randnum ++)
          randtab[randnum] = intRandVal (INTVALMAX);
        if (fwrite (randtab, sizeof (INT), RANDNBR, fileptr) < RANDNBR) {
          SCOTCH_errorPrint ("main: cannot write to file (3)");
          exit (EXIT_FAILURE);
        }
        break;
      case 1 :
        printf ("Random state cannot be saved\n");
        break;
      default :
        SCOTCH_errorPrint ("Could not save random state");
        break;
    }

    sleep (2);                                    /* Next run will not get the same time() value */
  }
  else {                                          /* Second pass */
    char *              charptr;
    int                 o;

    if (fread (randtab, sizeof (INT), RANDNBR, fileptr) < RANDNBR) {
      SCOTCH_errorPrint ("main: cannot read from file (1)");
      exit (EXIT_FAILURE);
    }

    for (randnum = 0; randnum < RANDNBR; randnum ++) {
      if (randtab[randnum] != intRandVal (INTVALMAX))
        break;
    }

    o = (randnum == RANDNBR);
    charptr = (o) ? "same" : "different";
#if ((defined COMMON_DEBUG) || (defined COMMON_RANDOM_FIXED_SEED) || (defined SCOTCH_DETERMINISTIC))
    o ^= 1;
#endif /* ((defined COMMON_DEBUG) || (defined COMMON_RANDOM_FIXED_SEED) || (defined SCOTCH_DETERMINISTIC)) */

    if (o) {
      SCOTCH_errorPrint ("main: two consecutive runs yield %s values.", charptr);
      exit (EXIT_FAILURE);
    }
    printf ("Two consecutive runs yield %s values.\n", charptr);

    intRandReset ();                              /* Reset random seed to be sure */

    switch (SCOTCH_randomLoad (fileptr)) {
      case 0 :
        while (1) {                               /* Discard all CR(LF) before separator character */
          int                 c;

          c = getc (fileptr);
          if (c == '#')                           /* If separator character found */
            break;
          if (c == EOF) {
            SCOTCH_errorPrint ("main: cannot read from file (2)");
            exit (EXIT_FAILURE);
          }
        }

        if (fread (randtab, sizeof (INT), RANDNBR, fileptr) < RANDNBR) {
          SCOTCH_errorPrint ("main: cannot read from file (3)");
          exit (EXIT_FAILURE);
        }

        for (randnum = 0; randnum < RANDNBR; randnum ++) {
          if (randtab[randnum] != intRandVal (INTVALMAX)) {
            SCOTCH_errorPrint ("main: state not properly saved/restored");
            exit (EXIT_FAILURE);
          }
        }
        break;
      case 1 :
        SCOTCH_errorPrint ("main: random state cannot be loaded");
        break;
      default :
        SCOTCH_errorPrint ("main: could not save random state");
        break;
    }
  }

  fclose (fileptr);
  free   (randtab);

  exit (EXIT_SUCCESS);
}

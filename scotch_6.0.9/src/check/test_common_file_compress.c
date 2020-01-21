/* Copyright 2018 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_common_file_compress.c             **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the random number     **/
/**                generator module.                       **/
/**                                                        **/
/**   DATES      : # Version 6.0  : from : 09 jul 2018     **/
/**                                 to     10 jul 2018     **/
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

#define C_FILENBR                   2            /* Number of files in list */

#define C_BUFFSIZ                   65536

#define C_filepntrinp               fileBlockFile (C_fileTab, 0) /* Input file  */
#define C_filepntrout               fileBlockFile (C_fileTab, 1) /* Output file */

/*
**  The static and global variables.
*/

static File                 C_fileTab[2] = {      /* File array */
                              { FILEMODER },
                              { FILEMODEW } };

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
  byte *              bufftab;

  if (argc != 3) {
    SCOTCH_errorPrint ("usage: %s infile outfile", argv[0]);
    exit (EXIT_FAILURE);
  }

  if ((bufftab = malloc (C_BUFFSIZ * sizeof (byte))) == NULL) {
    SCOTCH_errorPrint ("main: out of memory");
    exit (EXIT_FAILURE);
  }

  fileBlockInit (C_fileTab, C_FILENBR);           /* Set default stream pointers */
  fileBlockName (C_fileTab, 0) = argv[1];
  fileBlockName (C_fileTab, 1) = argv[2];

  switch (fileBlockOpen (C_fileTab, C_FILENBR)) { /* Open all files */
    case 2 :
      SCOTCH_errorPrint ("main: (de)compression method not implemented");
      free (bufftab);
      exit (EXIT_SUCCESS);
    case 1 :
      SCOTCH_errorPrint ("main: cannot open files");
      free (bufftab);
      exit (EXIT_FAILURE);
  }

  while (1) {
    size_t              bytenbr;

    if ((bytenbr = fread (bufftab, 1, C_BUFFSIZ, C_filepntrinp)) == 0) {
      if (ferror (C_filepntrinp))
        SCOTCH_errorPrint ("main: read error");
      break;
    }
    if (fwrite (bufftab, 1, bytenbr, C_filepntrout) != bytenbr) {
      SCOTCH_errorPrint ("main: write error");
      break;
    }
  }

  free (bufftab);

  fileBlockClose (C_fileTab, C_FILENBR);          /* Always close explicitely to end eventual (de)compression tasks */

#ifdef COMMON_PTHREAD
  pthread_exit ((void *) 0);                      /* Allow potential (de)compression tasks to complete */
#endif /* COMMON_PTHREAD */
  exit (EXIT_SUCCESS);
}

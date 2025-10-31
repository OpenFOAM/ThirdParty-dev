/* Copyright 2007-2009,2020,2023,2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgraph_fold_dup.c                       **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module folds a distributed graph   **/
/**                into two distinct copies, which may     **/
/**                be different when the number of         **/
/**                processes is odd.                       **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 10 aug 2006     **/
/**                                 to   : 20 jun 2007     **/
/**                # Version 5.1  : from : 14 nov 2008     **/
/**                                 to   : 28 oct 2009     **/
/**                # Version 6.0  : from : 28 sep 2014     **/
/**                                 to   : 28 sep 2014     **/
/**                # Version 7.0  : from : 03 sep 2020     **/
/**                                 to   : 12 sep 2024     **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#define SCOTCH_DGRAPH_FOLD_DUP

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "dgraph_fold_dup.h"

/******************************/
/*                            */
/* This routine handles       */
/* distributed source graphs. */
/*                            */
/******************************/

/* This routine builds two folded graphs of
** a given graph on each of the two halves
** of the processes. The number of processes
** does not need to be even. There is a
** multi-threaded version, as well as a
** sequential one.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

static
void
dgraphFoldDup2 (
Context * restrict const          contptr,        /*+ (Sub-)context (not used)               +*/
const int                         spltnum,        /*+ Rank of sub-context in initial context +*/
const DgraphFoldDupSplit * const  spltptr)
{
  int                 o;

  o = dgraphFold2 (spltptr->splttab[spltnum].orggrafptr, spltnum, spltptr->fldgrafptr,
                   spltptr->splttab[spltnum].fldproccomm, spltptr->orgdataptr, spltptr->flddataptr, spltptr->datatype);

  if (o != 0)
    *spltptr->revaptr = 1;                        /* No mutex protection */
}

int
dgraphFoldDup (
const Dgraph * restrict const orggrafptr,
Dgraph * restrict const       fldgrafptr,
void * restrict const         orgdataptr,         /*+ Un-based array of data which must be folded, e.g. coarmulttab +*/
void ** restrict const        flddataptr,         /*+ Un-based array of data which must be folded, e.g. coarmulttab +*/
MPI_Datatype                  datatype,           /*+ Data type for array of data which must be folded              +*/
Context * restrict const      contptr)            /*+ Context                                                       +*/
{
#ifdef SCOTCH_PTHREAD_MPI
  Dgraph              orggrafdat;
  int                 thrdprolvl;
  int                 thrdglbmin;
#endif /* SCOTCH_PTHREAD_MPI */
  int                 thrdval;                    /* Flag set if multithreaded process is possible */
  int                 fldprocnbr;
  int                 fldprocnum;
  int                 fldproccol;
  DgraphFoldDupSplit  fldspltdat;
  int                 o;

  fldprocnbr = (orggrafptr->procglbnbr + 1) / 2;  /* Median cut on number of processors     */
  if (orggrafptr->proclocnum < fldprocnbr) {      /* Compute color and rank in two subparts */
    fldproccol = 0;
    fldprocnum = orggrafptr->proclocnum;
  }
  else {
    fldproccol = 1;
    fldprocnum = orggrafptr->proclocnum - fldprocnbr;
  }
  if (MPI_Comm_split (orggrafptr->proccomm, fldproccol, fldprocnum, &fldspltdat.splttab[fldproccol].fldproccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphFoldDup: communication error (1)");
    return (1);
  }
  fldspltdat.splttab[fldproccol ^ 1].fldproccomm = MPI_COMM_NULL;

  fldspltdat.splttab[0].orggrafptr = orggrafptr;
  fldspltdat.orgdataptr = orgdataptr;
  fldspltdat.flddataptr = flddataptr;
  fldspltdat.fldgrafptr = fldgrafptr;
  fldspltdat.datatype   = datatype;
  fldspltdat.revaptr    = &o;

  o = 0;                                          /* Assume splitting will go well       */
  thrdval = 0;                                    /* Assume splitting will be sequential */

#ifdef SCOTCH_PTHREAD_MPI
  MPI_Query_thread (&thrdprolvl);                 /* Get thread level of MPI implementation */
  if (thrdprolvl >= MPI_THREAD_MULTIPLE) {        /* If multiple threads can be used        */
    int                 thrdlocmin;

    thrdlocmin = contextThreadNbr (contptr);
    if (MPI_Allreduce (&thrdlocmin, &thrdglbmin, 1, MPI_INT, MPI_MIN, orggrafptr->proccomm) != MPI_SUCCESS) {
      errorPrint ("dgraphFoldDup: communication error (2)");
      return (1);
    }

    if (thrdglbmin > 1) {                         /* If all processes have multiple threads available */
      fldspltdat.splttab[1].orggrafptr = &orggrafdat;
      orggrafdat = *orggrafptr;                   /* Create a separate graph structure to change its communicator */

      if (MPI_Comm_dup (orggrafptr->proccomm, &orggrafdat.proccomm) != MPI_SUCCESS) { /* Duplicate communicator to avoid interferences in communications */
        errorPrint ("dgraphFoldDup: communication error (3)");
        return (1);
      }

#ifndef DGRAPHFOLDDUPNOTHREAD
      if (contextThreadLaunchSplit (contptr, (ContextSplitFunc) dgraphFoldDup2, &fldspltdat) == 0) /* If context could be split to run concurrently */
        thrdval = 1;                                /* No need to go through sequantial run */
#endif /* DGRAPHFOLDDUPNOTHREAD */
    }
  }
#endif /* SCOTCH_PTHREAD_MPI */

  if (thrdval == 0) {                             /* If need to go through the sequential run               */
    fldspltdat.splttab[1].orggrafptr = orggrafptr; /* No need for separate graph with separate communicator */

    dgraphFoldDup2 (contptr, 0, &fldspltdat);     /* Run tasks in sequence */
    if (o == 0)
      dgraphFoldDup2 (contptr, 1, &fldspltdat);
  }
#ifdef SCOTCH_PTHREAD_MPI
  if (thrdprolvl >= MPI_THREAD_MULTIPLE) {        /* If duplicated communicator was created, free it  */
    if (thrdglbmin > 1)                           /* If all processes have multiple threads available */
      MPI_Comm_free (&orggrafdat.proccomm);
  }
#endif /* SCOTCH_PTHREAD_MPI */

  fldgrafptr->pkeyglbval = fldproccol;            /* Discriminate between folded communicators at same level */

  return (o);
}

/* Copyright 2007,2009-2011,2023 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : dgraph_fold_comm.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module computes the communication  **/
/**                pattern of the distributed graph        **/
/**                folding process.                        **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 23 may 2006     **/
/**                                 to   : 10 sep 2007     **/
/**                # Version 5.1  : from : 18 jan 2009     **/
/**                                 to   : 10 sep 2011     **/
/**                # Version 7.0  : from : 18 jun 2023     **/
/**                                 to   : 12 aug 2023     **/
/**                                                        **/
/************************************************************/

/*
** The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "dgraph.h"
#include "dgraph_coarsen.h"
#include "dgraph_fold_comm.h"

/* This routine computes an upper bound on the size of
** local vertex arrays that have to be allocated before
** coarsening (and possibly folding) take place, notably
** the multinode array.
** Since this value depends on the folding distribution
** algorithm, it is placed in the same source file.
** The extreme case is when vertglbmax receivers have
** (coarvertlocavg - 1) vertices, and a sender has
** vertglbmax vertices as well. Only DGRAPHFOLDCOMMNBR
** vertices will be transferred to achieve balance, and
** (vertglbmax - DGRAPHFOLDCOMMNBR) vertices will have
** to be spread across DGRAPHFOLDCOMMNBR processes.
** It returns:
** - >= 0  : size of vertex arrays
*/

Gnum
dgraphCoarsenVertLocMax (
const Dgraph * restrict const finegrafptr,        /*+ Graph to coarsen and possibly fold +*/
const int                     flagval)            /*+ Coarsening type                    +*/
{
  Gnum                coarvertlocavg;             /* Average load on all processes if perfect balance achieved */
  Gnum                coarvertlocmax;             /* Possible Maximum load on a specific process               */
  int                 foldval;                    /* Type of folding (do not consider merging)                 */

  const int                         procglbnbr = finegrafptr->procglbnbr;
  const Gnum                        vertlocnbr = finegrafptr->vertlocnbr;
  const Gnum                        vertglbnbr = finegrafptr->vertglbnbr;
  const Gnum                        vertglbmax = finegrafptr->vertglbmax;

  foldval = flagval & (DGRAPHCOARSENFOLD | DGRAPHCOARSENFOLDDUP); /* Only consider folding flags */

  if ((foldval == DGRAPHCOARSENNONE) ||           /* If plain coarsening                                                      */
      (procglbnbr == 1))                          /* Or user called the routine with a mono-process context                   */
    return (vertlocnbr);                          /* Local number of vertices may not change if all local matchings succeeded */

  if (foldval == DGRAPHCOARSENFOLD)               /* If simple folding, with coarsening ratio 1.0              */
    coarvertlocavg = ((vertglbnbr * 2) / procglbnbr) + 1; /* Max ratio for FOLD is 2 -> 1                      */
  else                                            /* Folding with duplication with coarsening ratio 1.0        */
    coarvertlocavg = ((vertglbnbr * 2) / (procglbnbr - (procglbnbr % 2))) + 1; /* Max ratio FOLD-DUP is 3 -> 1 */

  if (procglbnbr < (2 * DGRAPHFOLDCOMMNBR))       /* If there is always enough communications to spread the load of a single process */
    coarvertlocmax = coarvertlocavg;              /* No overload is possible in this case                                            */
  else
    coarvertlocmax = coarvertlocavg + ((vertglbmax - DGRAPHFOLDCOMMNBR + (DGRAPHFOLDCOMMNBR - 1)) / DGRAPHFOLDCOMMNBR); /* Add possible overload */

  return (coarvertlocmax);
}

/* This routine computes an optimized communication
** scheme for folding the data of a distributed graph.
** It is currently based on a maximum fixed number of
** communications per process. If this maximum is reached,
** the algorithm will fail.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
dgraphFoldComm (
const Dgraph * restrict const                   grafptr,
const int                                       partval, /* 0 for first half, 1 for second half                           */
int * restrict const                            commmaxptr, /* Pointer to maximum number of communications per process    */
int * restrict const                            commtypvalptr, /* Process will be sender or receiver                      */
DgraphFoldCommData * restrict * restrict const  commdattabptr, /* Slots for communication                                 */
Gnum * restrict * restrict const                commvrttabptr, /* Slots of starting global vertex send indices            */
Gnum * restrict const                           proccnttab, /* Receive count array, for receivers                         */
int * restrict const                            vertadjnbrptr, /* Number of adjustment ranges, for receivers              */
Gnum * restrict * restrict const                vertadjtabptr, /* Pointer to global index adjustment array, for receivers */
Gnum * restrict * restrict const                vertdlttabptr) /* Pointer to global delta adjustment array, for receivers */
{
  Gnum * restrict               fldproccnttab;    /* Sanitized access to proccnttab of folded graph (used as flag)  */
  DgraphFoldCommData * restrict commdattab;       /* Array of communication slots for current process               */
  Gnum * restrict               commvrttab;       /* Adday of starting global vertex indices of communication slots */
  DgraphFoldCommData * restrict procsrttab;       /* Sort array for sender and receiver processes                   */
  int                           procfldnbr;       /* Number of processes in folded part                             */
  int                           procrcvbas;
  int                           procrcvnnd;
  int                           procsndbas;
  int                           procsndnnd;
  int                           sortsndbas;       /* Start index of first non-empty sender in sort array            */
  int * restrict                slotsndtab;       /* Array of send slot indices for adjustment range computations   */
  Gnum * restrict               vertadjtab;       /* Arrays of global index ajustments for all receiver slots       */
  Gnum * restrict               vertdlttab;
  int                           commmax;          /* Current maximum number of communications per process           */
  int                           procnum;
#ifdef SCOTCH_DEBUG_DGRAPH2
  DgraphFoldCommData * restrict procchktab;
  int                           chekloctab[2];
  int                           chekglbtab[2];
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  const int                         orgprocglbnbr = grafptr->procglbnbr;
  const int                         orgproclocnum = grafptr->proclocnum;
  const Gnum * restrict const       orgproccnttab = grafptr->proccnttab;
  const Gnum * restrict const       orgprocvrttab = grafptr->procvrttab;

#ifdef SCOTCH_DEBUG_DGRAPH2
  if (orgprocglbnbr < 2) {                        /* Folding is useless when graph is on a single process */
    errorPrint ("dgraphFoldComm: invalid parameters (1)");
    return (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  if ((procsrttab = memAlloc (orgprocglbnbr * sizeof (DgraphFoldCommData))) == NULL) {
    errorPrint ("dgraphFoldComm: out of memory (1)");
    return (1);
  }

  procfldnbr = (orgprocglbnbr + 1) / 2;           /* Get number of processes in part 0 (always more than in part 1) */
  if (partval == 0) {                             /* If part 0 will receive the data                                */
    procrcvbas = 0;                               /* Receive by ascending weight order in first part                */
    procrcvnnd =
    procsndbas = procfldnbr;                      /* Send by descending weight order in second part */
    procsndnnd = orgprocglbnbr;
  }
  else {                                          /* Part 1 will receive the data                  */
    procsndbas = 0;                               /* Send by descending weight order in first part */
    procsndnnd =
    procrcvbas = procfldnbr;                      /* Receive by ascending weight order in second part */
    procrcvnnd = orgprocglbnbr;
    procfldnbr = orgprocglbnbr - procfldnbr;      /* Number of processes in folded graph is in fact the smaller half */
  }

  if ((orgproclocnum >= procrcvbas) && (orgproclocnum < procrcvnnd)) { /* If calling process is a receiver          */
    fldproccnttab  = proccnttab;                  /* Point to vertex count array of receiver process                */
    *commtypvalptr = DGRAPHFOLDCOMMRECV;          /* Set local process type as receiver (pure receiver, by default) */
  }
  else {
    fldproccnttab  = NULL;                        /* Flag process as a sender */
    *commtypvalptr = DGRAPHFOLDCOMMSEND;
  }
#ifdef SCOTCH_DEBUG_DGRAPH2
  if ((*commtypvalptr == DGRAPHFOLDCOMMRECV) &&   /* Check that a receiver has the proper receiver array locations */
      ((proccnttab == NULL) || (vertadjtabptr == NULL) || (vertdlttabptr == NULL))) {
    errorPrint ("dgraphFoldComm: invalid parameters (2)");
    memFree    (procsrttab);                      /* Free group leader */
    return (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  for (procnum = procsndbas; procnum < procsndnnd; procnum ++) { /* Senders will not be changed in loop  */
    procsrttab[procnum].vertnbr = orgproccnttab[procnum]; /* All vertices of senders must be disposed of */
    procsrttab[procnum].procnum = (Gnum) procnum;
  }
  intSort2asc1 (procsrttab + procsndbas, procsndnnd - procsndbas); /* Sort sender part once for good */

  commdattab = NULL;                              /* Vertex redistribution arrays not allocated yet                        */
  for (commmax = DGRAPHFOLDCOMMNBR; ; commmax ++) { /* Start with small number of communications per process               */
    int * restrict      slotsndtax;               /* Based access to slot array, with respect to (procrcvbas * commsiz)    */
    int                 sortsndnnd;               /* Index of current end of send sub-array in sort array                  */
    int                 sortrcvnnd;               /* Index of current end of receive sub-array in sort array               */
    int                 sortrcvnum;               /* Index of current receiver process in sort array                       */
    Gnum                vertglbdlt;               /* Global overload accepted for receiver vertices                        */
    int                 commrcvnum;               /* Number of communications already performed by current receiver        */
    int                 commsiz;                  /* Size of receiver communication arrays for local process (commmax + 1) */
    int                 commnum;

    for (procnum = procrcvbas; procnum < procrcvnnd; procnum ++) { /* Compute difference with respect to average for receivers */
      procsrttab[procnum].vertnbr = orgproccnttab[procnum] - DATASIZE (grafptr->vertglbnbr, procfldnbr, procnum - procrcvbas);
      procsrttab[procnum].procnum = (Gnum) procnum;
    }
    intSort2asc1 (procsrttab + procrcvbas, procrcvnnd - procrcvbas); /* (Re)sort receiver part as may have been modified by previous loop */

    for (sortsndbas = procsndbas; sortsndbas < procsndnnd; sortsndbas ++) { /* Discard empty senders */
      if (procsrttab[sortsndbas].vertnbr != 0)    /* Stop at first non-empty sender                  */
        break;
    }

    if (commdattab != NULL)                       /* If retry, delete old arrays before rebuilding them */
      memFree (commdattab);

    commsiz = (fldproccnttab != NULL) ? (commmax + 1) : 0; /* Receiver-side arrays only necessary if process is a receiver */
    if (memAllocGroup ((void **) (void *)
                       &commdattab, (size_t) (commmax * sizeof (DgraphFoldCommData)),
                       &commvrttab, (size_t) (commmax * sizeof (Gnum)),
                       &vertadjtab, (size_t) (commsiz * orgprocglbnbr * sizeof (Gnum)),
                       &vertdlttab, (size_t) (commsiz * orgprocglbnbr * sizeof (Gnum)),
#ifdef SCOTCH_DEBUG_DGRAPH2
                       &procchktab, (size_t) (commmax * orgprocglbnbr * sizeof (DgraphFoldCommData)),
#endif /* SCOTCH_DEBUG_DGRAPH2 */
                       &slotsndtab, (size_t) (commsiz * procfldnbr * sizeof (int)), NULL) == NULL) {
      errorPrint ("dgraphFoldComm: out of memory (2)");
      memFree    (procsrttab);
      return (1);
    }

    for (commnum = 0; commnum < commmax; commnum ++) { /* Assume we will perform no communication at all */
      commdattab[commnum].vertnbr = 0;
      commdattab[commnum].procnum = -1;
    }
    if (fldproccnttab != NULL) {                  /* If we are a receiver process, initialize count and adjustment arrays */
      int                 procrcvnum;

      memSet (vertadjtab, ~0, commsiz * orgprocglbnbr * sizeof (Gnum)); /* Initialize index adjustment arrays */
      memSet (slotsndtab, ~0, commsiz * procfldnbr * sizeof (int));
      slotsndtax = slotsndtab - (procrcvbas * commsiz); /* Base access to slot array for receivers */

      for (procrcvnum = procrcvbas; procrcvnum < procrcvnnd; procrcvnum ++) { /* For all slots of receiver processes */
        int                 slotrcvnum;           /* First index in slot array for receiver process                  */

        slotrcvnum = procrcvnum * commsiz;
        vertadjtab[slotrcvnum] = orgprocvrttab[procrcvnum]; /* Set global start index of local chunk            */
        vertdlttab[slotrcvnum] = orgproccnttab[procrcvnum]; /* Set number of vertices to be kept in local chunk */
        slotsndtax[slotrcvnum] = slotrcvnum;      /* First slot of receiver is active                           */
      }
    }

    sortsndnnd = procsndnnd - 1;                  /* TRICK: point to end processes in arrays rather than just after them */
    sortrcvnnd = procrcvnnd - 1;
    sortrcvnum = procrcvbas;                      /* First receiver to start with is receiver of smallest load               */
    commrcvnum = 0;                               /* Current receiver has not started any remote communication               */
    vertglbdlt = 0;                               /* No overload yet among receiver processes                                */
    while (1) {                                   /* Loop on candidate processes for setting up communications               */
      int                 flagsndval;             /* Type of current sender: pure sender or sender receiver                  */
      int                 procsndnum;             /* Index of current sender process                                         */
      int                 slotsndnum;             /* Index of first slot of current sender                                   */
      int                 sortrcvtmp;             /* Temporary receive index for aggregating receiver capacities             */
      int                 sortsndnum;             /* Index of sender process in sort array (either send or receive)          */
      Gnum                vertsndnbr;             /* Number of vertices to send from current sender process                  */
      Gnum                vertadjval;             /* Global index of first vertex to send in a chunk                         */
      Gnum                vertrcvnbr;             /* Number of vertices a receiver can/will take at a time                   */
      Gnum                vertrcvdlt;             /* Room available to interleave pure senders with receiver sender messages */
      int                 commsndnum;             /* Number of communications already performed by current sender            */
      int                 commnbr;                /* Number of elementary communcations to achive vertex transfer            */

      vertsndnbr = 0;                             /* Nothing sent yet                  */
      if (sortsndnnd >= sortsndbas) {             /* If a sender process remains       */
        sortsndnum = sortsndnnd --;               /* Select and remove it              */
        vertsndnbr = procsrttab[sortsndnum].vertnbr; /* Get number of vertices to send */
        flagsndval = DGRAPHFOLDCOMMSEND;          /* Sender is a true sender           */
      }
      if (sortrcvnnd >= sortrcvnum) {             /* If a receiver process remains */
        Gnum                vertsndtmp;

        vertsndtmp = procsrttab[sortrcvnnd].vertnbr; /* Get number of vertices possibly to send                       */
        if ((vertsndtmp > vertsndnbr) &&          /* If receiver has something to send (> 0) and would be prioritary  */
            (vertsndtmp > vertglbdlt)) {          /* And sending it would contribute not to degrade imbalance further */
          if (vertsndnbr > 0)                     /* If a sender vertex was chosen before                             */
            sortsndnnd ++;                        /* Restore it                                                       */

          sortsndnum = sortrcvnnd --;             /* Select and remove sender receiver              */
          vertsndnbr = vertsndtmp - vertglbdlt;   /* Number of vertices to send, keeping local part */
          flagsndval = DGRAPHFOLDCOMMSEND | DGRAPHFOLDCOMMRECV; /* Sender is a sender receiver      */
        }
      }
      else {                                      /* No more receiver processes available               */
        if (vertsndnbr > 0)                       /* If a sender was willing to communicate             */
          goto redo;                              /* We could not compute a valid communication pattern */
      }
      if (vertsndnbr <= 0)                        /* If no more communications to perform          */
        goto done;                                /* We have computed a valid communcation pattern */

      vertrcvdlt = 0;                             /* Assume we will not be able to interleave pure senders      */
      for (commnbr = 0, vertrcvnbr = 0, sortrcvtmp = sortrcvnum; /* Compute maximum available room on receivers */
           (commnbr < commmax) && (sortrcvtmp <= sortrcvnnd); commnbr ++, sortrcvtmp ++)
        vertrcvnbr += vertglbdlt - procsrttab[sortrcvtmp].vertnbr;
      if (vertrcvnbr < vertsndnbr) {              /* If receiver processes cannot hold the data without additional overload */
        if (flagsndval == DGRAPHFOLDCOMMSEND)     /* If sender is a true sender                                             */
          vertglbdlt += (vertsndnbr - vertrcvnbr + (commnbr - 1)) / commnbr; /* Spread extra overload across receivers      */
        else {                                    /* Sender is a sender receiver                                            */
          Gnum                vertglbtmp;

          vertglbtmp  = (vertsndnbr - vertrcvnbr + commnbr) / (commnbr + 1); /* Compute overload including sender process */
          vertglbdlt += vertglbtmp;               /* Add extra overload to current overload                               */
          vertsndnbr -= vertglbtmp;               /* Extra overload vertices will not be sent by sender receiver process  */
        }
      }
      else {                                      /* There is room for adding messages       */
        if (flagsndval != DGRAPHFOLDCOMMSEND)     /* If selected sender is a sender receiver */
          vertrcvdlt = vertrcvnbr - vertsndnbr;   /* Compute remaining space among receivers */
      }

      procsndnum = (int) procsrttab[sortsndnum].procnum; /* Get rank of sender process                                             */
      slotsndnum = procsndnum * commsiz;          /* Get first slot of sender process                                              */
      commsndnum = 0;                             /* Sender has not communicated yet                                               */
      vertadjval = orgprocvrttab[procsndnum] + orgproccnttab[procsndnum] - vertsndnbr; /* Get global index of first vertex to send */

      do {                                        /* Loop on communications to exchange chunks of vertices */
        int                 procrcvnum;
        int                 slotrcvnum;           /* Index of first slot in receiver part */
        int                 flagrcvval;           /* Flag set if receiver will be full    */

#ifdef SCOTCH_DEBUG_DGRAPH2
        if (commsndnum >= commmax) {              /* Overload should have been computed so that this never happens */
          errorPrint ("dgraphFoldComm: internal error (1)");
          return (1);
        }
        if (sortrcvnum > sortrcvnnd) {            /* Overload should have been computed so that this never happens */
          errorPrint ("dgraphFoldComm: internal error (2)");
          return (1);
        }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

        procrcvnum = (int) procsrttab[sortrcvnum].procnum; /* Get receiver index                                    */
        slotrcvnum = procrcvnum * commsiz + 1;        /* Start slot of receiver is after local slot (even if empty) */
        vertrcvnbr = vertglbdlt - procsrttab[sortrcvnum].vertnbr; /* Compute capacity of receiver process           */

        if ((vertrcvdlt >  0) &&                  /* If there is a possibility to interleave a pure sender */
            (commrcvnum == 0) &&                  /* And we are considering a new receiver process         */
            (sortsndbas <= sortsndnnd) &&         /* And there is a pure sender available                  */
            (vertrcvdlt >= procsrttab[sortsndbas].vertnbr) && /* And it can fit in the available space     */
            (vertrcvnbr >= procsrttab[sortsndbas].vertnbr)) { /* And it will not overload the receiver     */
          int                 procsndtmp;
          int                 slotsndtmp;
          Gnum                vertrcvtmp;
          Gnum                vertadjtmp;

          procsndtmp = (int) procsrttab[sortsndbas].procnum; /* Sender process for communication is interleaved sender */
          slotsndtmp = procsndtmp * commsiz;      /* Interleaved sender will use only its first communication          */
          vertrcvtmp = procsrttab[sortsndbas].vertnbr; /* Receive all of its vertices in one chunk                     */
          vertadjtmp = orgprocvrttab[procsndtmp]; /* Send all of its vertices in one chunk                             */
          sortsndbas ++;                          /* Discard interleaved sender                                        */
          vertrcvdlt -= vertrcvtmp;               /* Reduce interleaving capacity                                      */

          if (fldproccnttab != NULL) {            /* If we are a receiver process, fill count and adjustment arrays */
            if (procrcvnum == orgproclocnum) {    /* If we are the receiver process                                 */
              commdattab[commrcvnum].vertnbr = vertrcvtmp; /* Record receive communication                          */
              commdattab[commrcvnum].procnum = (Gnum) procsndtmp; /* Sender is interleaved sender                   */
              commvrttab[commrcvnum] = vertadjtmp;
            }

            vertadjtab[slotsndtmp] = vertadjtmp;  /* Record communication in slot */
            vertdlttab[slotsndtmp] = vertrcvtmp;
            slotsndtax[slotrcvnum + commrcvnum] = slotsndtmp; /* Receiver slot points to sender slot */
          }
          else {                                  /* We are a sender process      */
            if (procsndtmp == orgproclocnum) {    /* If we are the sender process */
              commdattab[0].vertnbr = vertrcvtmp;   /* Record communication         */
              commdattab[0].procnum = (Gnum) procrcvnum;
              commvrttab[0] = vertadjtmp;
            }
          }

          commrcvnum ++;                          /* One more communication performed for receiver processes */
          vertrcvnbr -= vertrcvtmp;               /* Receiver's capacity has been reduced                    */
          if (vertrcvnbr <= 0) {                  /* If receiver's capacity has been exhausted               */
            sortrcvnum ++;                        /* Switch to next receiver                                 */
            commrcvnum = 0;

            procrcvnum = (int) procsrttab[sortrcvnum].procnum; /* Get receiver index                                    */
            slotrcvnum = procrcvnum * commsiz + 1;        /* Start slot of receiver is after local slot (even if empty) */
            vertrcvnbr = vertglbdlt - procsrttab[sortrcvnum].vertnbr; /* Compute capacity of receiver process           */
          }
        }

        flagrcvval = 1;                           /* Assume receiver will be full                       */
        if (vertsndnbr < vertrcvnbr) {            /* If less vertices (remaining) to send than capacity */
          vertrcvnbr = vertsndnbr;                /* Only send and receive the given amount             */
          flagrcvval = 0;                         /* Receiver will not be considered full               */
        }
        vertsndnbr -= vertrcvnbr;                 /* Remove number of received vertices from amount to send */

        if (fldproccnttab != NULL) {              /* If we are a receiver process, fill count and adjustment arrays */
          int                 slotsndtmp;         /* Index of current slot for communication send                   */

          slotsndtmp = 0;                         /* Assume sender is a pure sender                               */
          if (flagsndval != DGRAPHFOLDCOMMSEND) { /* If sender is a sender receiver                               */
            vertdlttab[slotsndnum] -= vertrcvnbr; /* Remove vertices from the initial slot of the sender receiver */
            slotsndtmp = 1;                       /* Send slots will be created after local receiver slot         */

            if (procsndnum == orgproclocnum) {    /* If we are the sender receiver process */
              *commtypvalptr = flagsndval;        /* Indicate it (maybe multiple times!)   */
              commdattab[commsndnum].vertnbr = vertrcvnbr; /* Record send communication    */
              commdattab[commsndnum].procnum = (Gnum) procrcvnum;
              commvrttab[commsndnum] = vertadjval;
            }
          }
          if (procrcvnum == orgproclocnum) {      /* If we are the receiver process      */
            commdattab[commrcvnum].vertnbr = vertrcvnbr; /* Record receive communication */
            commdattab[commrcvnum].procnum = (Gnum) procsndnum;
            commvrttab[commrcvnum] = vertadjval;
          }

          slotsndtmp += slotsndnum + commsndnum;  /* Compute actual send slot index */
          vertadjtab[slotsndtmp] = vertadjval;    /* Record communication in slot   */
          vertdlttab[slotsndtmp] = vertrcvnbr;
          slotsndtax[slotrcvnum + commrcvnum] = slotsndtmp; /* Receiver slot points to sender slot */
        }
        else {                                    /* We are a sender process      */
          if (procsndnum == orgproclocnum) {      /* If we are the sender process */
            commdattab[commsndnum].vertnbr = vertrcvnbr; /* Record communication  */
            commdattab[commsndnum].procnum = (Gnum) procrcvnum;
            commvrttab[commsndnum] = vertadjval;
          }
        }

        commrcvnum ++;                            /* One more communication performed for receiver and sender processes */
        commsndnum ++;
        if ((flagrcvval != 0) ||                  /* If receiver has been filled by previous communication    */
            (commrcvnum == commmax)) {            /* Or maximum number of communications reached for receiver */
          sortrcvnum ++;                          /* Switch to next receiver                                  */
          commrcvnum = 0;
        }
        vertadjval += vertrcvnbr;                 /* Update start global vertex index of next chunk  */
      } while (vertsndnbr > 0);                   /* As long as there are vertices to send           */
      if (commrcvnum > 0) {                       /* If current receiver may receive more messages   */
        Gnum                vertrcvtmp;

        vertrcvtmp = procsrttab[sortrcvnum].vertnbr + vertrcvnbr; /* Update its capacity from last communication */
        if (vertrcvtmp < 0)                       /* If vertex has not been overloaded by previous communication */
          procsrttab[sortrcvnum].vertnbr = vertrcvtmp; /* Update capacity to be taken into account in next round */
        else {                                    /* Else switch to next receiver                                */
          sortrcvnum ++;
          commrcvnum = 0;
        }
      }
    }
redo: ;                                           /* Increase maximum number of communications and retry */
  }
done:                                             /* A valid communication pattern has been produced */

  memFree (procsrttab);                           /* Sort array is no longer necessary */

#ifdef SCOTCH_DEBUG_DGRAPH2
  chekloctab[0] = - commmax;
  chekloctab[1] =   commmax;
  if (MPI_Allreduce (chekloctab, chekglbtab, 2, MPI_INT, MPI_MAX, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphFoldComm: communication error (1)");
    memFree    (commdattab);                      /* Free group leader */
    return (1);
  }
  if ((chekglbtab[0] != chekloctab[0]) ||
      (chekglbtab[1] != chekloctab[1])) {
    errorPrint ("dgraphFoldComm: internal error (3)");
    memFree    (commdattab);                      /* Free group leader */
    return (1);
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  if (fldproccnttab != NULL) {                    /* If we are a receiver process, finalize arrays  */
    Gnum *              fldproccntptr;            /* Pointer to current cell of proccnttab to fill  */
    Gnum                fldvertlocnbr;            /* Number of vertices in folded receiver process  */
    int                 vertadjnbr;               /* Number of adjustment slots                     */
    Gnum                vertadjsum;               /* Current new starting position of current slot  */
    int                 slotglbnbr;               /* Number of slots in the full slot array         */
    int                 slotglbnum;               /* Index of current slot in the full slot array   */
    int                 slotrcvnbr;               /* Number of slots in the receiver part of arrays */
    int                 slotrcvnum;               /* Index of current slot in receiver part         */
    int                 commsiz;                  /* Number of communication slots per process      */
    int                 commnum;                  /* Index of communication for current receiver    */

    commsiz = commmax + 1;
    for (slotrcvnum = 0, commnum = 0, fldvertlocnbr = 0, fldproccntptr = fldproccnttab, /* For all receiver processes and slots */
         slotrcvnbr = commsiz * procfldnbr, vertadjsum = grafptr->baseval; /* Compute slot adjustments and number of vertices   */
         slotrcvnum < slotrcvnbr; slotrcvnum ++, commnum ++) {
      Gnum                vertadjtmp;
      int                 slotsndnum;             /* Index of sender slot matching the receiver slot */

      if (commnum == commsiz) {                   /* If reached the slots of a new receiver process */
        *(fldproccntptr ++) = fldvertlocnbr;      /* Record number of vertices in proccnttab        */
        fldvertlocnbr = 0;                        /* Reset number of vertices and of communications */
        commnum = 0;
      }

      slotsndnum = slotsndtab[slotrcvnum];        /* Get index of matching sender slot */
      if (slotsndnum == -1)                       /* Skip empty slots                  */
        continue;

      vertadjtmp = vertdlttab[slotsndnum];        /* Accumulate sender slot indices and compute adjustments */
      vertdlttab[slotsndnum] = vertadjsum - vertadjtab[slotsndnum];
      vertadjsum += vertadjtmp;
      fldvertlocnbr += vertadjtmp;
    }
    *fldproccntptr = fldvertlocnbr;               /* Record last cell in proccnttab */

    for (slotglbnum = vertadjnbr = 0, slotglbnbr = (commmax + 1) * orgprocglbnbr; /* Compact vertex adjustment arrays */
         slotglbnum < slotglbnbr; slotglbnum ++) {
      if (vertadjtab[slotglbnum] != -1) {         /* If slot is not empty, compact it */
        vertadjtab[vertadjnbr] = vertadjtab[slotglbnum];
        vertdlttab[vertadjnbr] = vertdlttab[slotglbnum];
        vertadjnbr ++;
      }
    }
    vertadjtab[vertadjnbr] = orgprocvrttab[orgprocglbnbr]; /* Set upper bound on global vertex indices */

    *vertadjnbrptr = vertadjnbr;
    *vertadjtabptr = vertadjtab;
    *vertdlttabptr = vertdlttab;
  }
  *commdattabptr = commdattab;                    /* Set group leader */
  *commvrttabptr = commvrttab;
  *commmaxptr    = commmax;

#ifdef SCOTCH_DEBUG_DGRAPH2
  if (fldproccnttab == NULL) {                    /* If we are a sender process       */
    Gnum                vertsndnbr;               /* Number of vertices actually sent */
    int                 commnum;

    for (commnum = 0, vertsndnbr = 0;             /* Accumulate number of vertices sent */
         (commnum < commmax) && (commdattab[commnum].procnum >= 0); commnum ++)
      vertsndnbr += commdattab[commnum].vertnbr;
    if (vertsndnbr != grafptr->vertlocnbr) {      /* All vertices should have been sent */
      errorPrint ("dgraphFoldComm: internal error (4)");
      commdattab[0].procnum = -2;                 /* Set collective error flag */
    }
  }
  if (MPI_Allgather (commdattab, 2 * commmax, GNUM_MPI,
                     procchktab, 2 * commmax, GNUM_MPI, grafptr->proccomm) != MPI_SUCCESS) {
    errorPrint ("dgraphFoldComm: communication error (2)");
    goto fail;
  }

  for (procnum = 0; procnum < orgprocglbnbr; procnum ++) { /* For all processes in initial graph */
    int                 commnum;

    if (commdattab[0].procnum == -2)              /* Delayed error handling */
      goto fail;

    for (commnum = 0; (commnum < commmax) && (procchktab[commmax * procnum + commnum].procnum != -1); commnum ++) { /* For all communications */
      Gnum                procend;
      Gnum                vertnbr;
      int                 commend;

      procend = procchktab[commmax * procnum + commnum].procnum; /* Get communication parameters */
      vertnbr = procchktab[commmax * procnum + commnum].vertnbr;

      for (commend = 0; commend < commmax; commend ++) { /* Search for matching communication */
        if ((procchktab[procend * commmax + commend].procnum == procnum) &&
            (procchktab[procend * commmax + commend].vertnbr == vertnbr))
          break;                                  /* If matching communication found, terminate search */
      }
      if (commend >= commmax) {                   /* If no matching communication found */
        errorPrint ("dgraphFoldComm: internal error (5)");
fail:
        memFree (commdattab);                     /* Free group leader */
        return  (1);
      }
    }
  }
#endif /* SCOTCH_DEBUG_DGRAPH2 */

  return (0);
}

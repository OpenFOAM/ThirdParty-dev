/* Copyright 2004,2007,2009,2020,2022,2023 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : esmumps.c                               **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains a MUMPS interface  **/
/**                for the ordering routines of the        **/
/**                libSCOTCH + Emilio libfax libraries.    **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 16 may 2001     **/
/**                                 to   : 04 jun 2001     **/
/**                # Version 0.1  : from : 13 feb 2002     **/
/**                                 to   : 13 feb 2002     **/
/**                # Version 1.0  : from : 06 dec 2004     **/
/**                                 to   : 06 dec 2004     **/
/**                # Version 5.1  : from : 22 jan 2009     **/
/**                                 to   : 22 jan 2009     **/
/**                # Version 6.0  : from : 22 jan 2020     **/
/**                                 to   : 22 jan 2020     **/
/**                # Version 6.1  : from : 28 aug 2020     **/
/**                                 to   : 05 sep 2020     **/
/**                # Version 7.0  : from : 01 dec 2022     **/
/**                                 to   : 21 jan 2023     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "graph.h"
#include "dof.h"
#include "symbol.h"
#include "order.h"
#include "fax.h"
#include "esmumps.h"

/**************************************/
/*                                    */
/* This routine acts as an interface  */
/* between ordering software such as  */
/* MUMPS and Scotch+Emilio.           */
/*                                    */
/**************************************/

/* Meaning of the parameters :
** - n : order of the system (that is, number of columns).
** - iwlen : not used. Here for compatibility.
** - pe : on input, position in array iw of the extra-diagonal
**   terms for the considered column.
**   on output, -pe(i) is the father of node i in the elimination
**   tree if i is a principal variable, or it is the index of
**   the principal variable if i is a secondary variable.
** - pfree : number of extra-diagonal terms for the considered
**   node (that is, the number of arcs dans le graph for this
**   vertex).
** - len : array holding the number of extra-diagonal terms for
**   each column.
** - iw : array of extra-diagonal terms (preserved).
** - nv : on output, nv(i) = 0 if variable i is a secondary
**   variable, else nv(i) is the number of columns that are
**   merged into principal variable i.
** - elen : on output, direct permutation (for MUMPS; the
**   meaning of the "direct" and "inverse" permutations is
**   just the opposite for Scotch) :
**   k=elen(i) <==> column i is the k-th pivot.
** - last : on output, inverse permutation (for MUMPS) :
**   i=last(k) <==> column i est le k-th pivot
*/

int
esmumps2 (
const INT                   n,
const INT                   pfree,                /*+ Number of arcs, plus 1                   +*/
INT * restrict const        petab,
INT * restrict const        lentab,
INT * restrict const        iwtab,
INT * restrict const        velotab,              /*+ NULL or nvtab                            +*/
INT * restrict const        nvtab,
INT * restrict const        elentab,              /*+ Permutations computed for debugging only +*/
INT * restrict const        lasttab)              /*+ Permutations computed for debugging only +*/
{
  INT                         baseval;            /* Base value                          */
  Graph                       grafdat;            /* Graph                               */
  Order                       ordedat;            /* Graph ordering                      */
  SymbolMatrix                symbdat;            /* Block factored matrix               */
  INT                         vertnum;
  INT * restrict              vendtab;            /* Vertex end array                    */
  const INT * restrict        velotax;            /* Based access to inverse permutation */
  const INT * restrict        peritax;            /* Based access to inverse permutation */
  const SymbolCblk * restrict cblktax;            /* Based access to column block array  */
  const SymbolBlok * restrict bloktax;            /* Based access to block array         */
  INT * restrict              nvtax;
  INT * restrict              petax;
  INT                         cblknum;

  if ((vendtab = memAlloc (n * sizeof (INT))) == NULL) {
    errorPrint ("esmumps2: out of memory");
    return (1);
  }
  for (vertnum = 0; vertnum < n; vertnum ++)
    vendtab[vertnum] = petab[vertnum] + lentab[vertnum];

  graphInit        (&grafdat);
  graphBuildGraph2 (&grafdat, 1, n, pfree - 1, petab, vendtab, velotab, NULL, iwtab, NULL); /* Assume Fortran-based indexing */

  orderInit  (&ordedat);
  orderGraph (&ordedat, &grafdat);                /* Compute ordering with Scotch */

#ifdef ESMUMPS_ORDER_DUMP
  {
    const char *      filestr;
    FILE *            fileptr;

    filestr = envGetStr ("ESMUMPS_ORDER_DUMP_FILE", "/tmp/esmumps_ordering.ord");

    if (filestr[0] != '\0') {                     /* If environment variable is not empty, perform dump of block ordering */
      if ((fileptr = fopen (filestr, "w+")) == NULL) {
        errorPrint ("esmumps2: cannot open ordering dump file");
        return (1);
      }

      orderSave (&ordedat, fileptr);

      fclose (fileptr);
    }
  }
#endif /* ESMUMPS_ORDER_DUMP */

#ifdef ESMUMPS_DEBUG                              /* Permutations are output for debugging only */
  memCpy (elentab, ordedat.permtab, n * sizeof (INT)); /* Copy permutations                     */
  memCpy (lasttab, ordedat.peritab, n * sizeof (INT));
#endif /* ESMUMPS_DEBUG */

  symbolInit     (&symbdat);
  symbolFaxGraph (&symbdat, &grafdat, &ordedat);  /* Compute block symbolic factorizaion */

#ifdef ESMUMPS_DEBUG_OUTPUT
  {
    Dof                 deofdat;                  /* Matrix DOF structure */
    double              fnnzval;
    double              fopcval;

    dofInit    (&deofdat);
    dofGraph   (&deofdat, &grafdat, 1, ordedat.peritab); /* Base on graph weights or constant load of 1 */
    symbolCost (&symbdat, &deofdat, SYMBOLCOSTLDLT, &fnnzval, &fopcval); /* Compute factorization cost  */

    fprintf (stderr, "ESMUMPS CblkNbr=" INTSTRING ", BlokNbr=" INTSTRING ", NNZ=%lg, OPC=%lg\n",
             symbdat.cblknbr, symbdat.bloknbr, fnnzval, fopcval);

    dofExit (&deofdat);
  }
#endif /* ESMUMPS_DEBUG_OUTPUT */

  baseval = 1;                                    /* Assume Fortran-based indexing */
  peritax = ordedat.peritab - baseval;            /* Based accesses                */
  cblktax = symbdat.cblktab - baseval;
  bloktax = symbdat.bloktab - baseval;
  nvtax   = nvtab - baseval;
  petax   = petab - baseval;
  velotax = (velotab != NULL) ? velotab - baseval : velotab;

  for (cblknum = baseval; cblknum < (symbdat.cblknbr + baseval); cblknum ++) { /* For all column blocks */
    INT                 degrval;                  /* True degree of column block                        */
    INT                 bloknum;
    INT                 colunum;

    if (velotax == NULL) {                        /* If graph has no vertex weights */
      for (bloknum = cblktax[cblknum].bloknum, degrval = 0;
           bloknum < cblktax[cblknum + 1].bloknum; bloknum ++)
        degrval += bloktax[bloknum].lrownum - bloktax[bloknum].frownum + 1;
    }
    else {                                        /* Graph has vertex weights */
      for (bloknum = cblktax[cblknum].bloknum, degrval = 0;
           bloknum < cblktax[cblknum + 1].bloknum; bloknum ++) {
        INT                 brownum;              /* Block row index in permuted matrix */

        for (brownum  = bloktax[bloknum].frownum;
             brownum <= bloktax[bloknum].lrownum; brownum ++)
          degrval += velotax[peritax[brownum]];
      }
    }
    nvtax[peritax[cblktax[cblknum].fcolnum]] = degrval; /* Set true block degree */
    petax[peritax[cblktax[cblknum].fcolnum]] = (cblktax[cblknum].bloknum == cblktax[cblknum + 1].bloknum - 1) /* If column block has no extra-diagonals */
      ? 0                                         /* Then mark block as root of subtree */
      : - peritax[cblktax[bloktax[cblktax[cblknum].bloknum + 1].cblknum].fcolnum];

    for (colunum  = cblktax[cblknum].fcolnum + 1;  /* For all secondary variables */
         colunum <= cblktax[cblknum].lcolnum; colunum ++) {
      nvtax[peritax[colunum]] = 0;                 /* Set nv = 0 and pe = - principal variable */
      petax[peritax[colunum]] = - peritax[cblktax[cblknum].fcolnum];
    }
  }

  symbolExit (&symbdat);
  orderExit  (&ordedat);
  graphExit  (&grafdat);
  memFree    (vendtab);

  return (0);
}

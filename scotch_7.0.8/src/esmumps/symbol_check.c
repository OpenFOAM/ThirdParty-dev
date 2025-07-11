/* Copyright 2004,2020,2023 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : symbol_check.c                          **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                This module checks the consistency of   **/
/**                symbolic matrices.                      **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 29 sep 1998     **/
/**                                 to   : 07 oct 1998     **/
/**                # Version 1.0  : from : 03 jun 2002     **/
/**                                 to   : 03 jun 2002     **/
/**                # Version 3.0  : from : 29 feb 2004     **/
/**                                 to   : 29 feb 2004     **/
/**                # Version 6.0  : from : 06 feb 2020     **/
/**                                 to   : 06 feb 2020     **/
/**                # Version 6.1  : from : 24 feb 2020     **/
/**                                 to   : 24 feb 2020     **/
/**                # Version 7.0  : from : 21 jan 2023     **/
/**                                 to   : 21 jan 2023     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "symbol.h"

/******************************************/
/*                                        */
/* The symbolic matrix handling routines. */
/*                                        */
/******************************************/

/*+ This routine checks the consistency
*** of the given symbolic block matrix.
*** Because of incomplete factorization,
*** from version 1.0, no check is performed
*** regarding the existence of facing blocks
*** in facing columns.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
symbolCheck (
const SymbolMatrix * const  symbptr)
{
  INT                         baseval;            /* Base value                           */
  const SymbolCblk * restrict cblktax;            /* Based access to cblktab              */
  INT                         cblkmax;            /* Maximum column block index           */
  INT                         cblknum;            /* Based number of current column block */
  const SymbolBlok * restrict bloktax;            /* Based access to bloktab              */
  INT                         blokmax;            /* Maximum block index                  */
  INT                         bloknum;            /* Based number of current block        */
  INT                         nodemax;            /* Maximum node index                   */

  baseval = symbptr->baseval;
  cblktax = symbptr->cblktab - baseval;
  cblkmax = symbptr->cblknbr + (baseval - 1);
  bloktax = symbptr->bloktab - baseval;
  blokmax = symbptr->bloknbr + baseval;
  nodemax = symbptr->nodenbr + (baseval - 1);

  for (cblknum = bloknum = baseval;
       cblknum <= cblkmax; cblknum ++) {
    if ((cblktax[cblknum].fcolnum     <  baseval)                  ||
        (cblktax[cblknum].lcolnum     >  nodemax)                  ||
        (cblktax[cblknum].bloknum     >  blokmax)                  ||
        (cblktax[cblknum].fcolnum     >  cblktax[cblknum].lcolnum) ||
        (cblktax[cblknum + 1].fcolnum <= cblktax[cblknum].lcolnum) ||
        (cblktax[cblknum + 1].bloknum <= cblktax[cblknum].bloknum)) {
      errorPrint ("symbolCheck: invalid column block array");
      return     (1);
    }

    if ((bloktax[bloknum].frownum != cblktax[cblknum].fcolnum) ||
        (bloktax[bloknum].lrownum != cblktax[cblknum].lcolnum) ||
#ifdef SYMBOL_HAS_LEVFVAL
        (bloktax[bloknum].levfval != 0) ||
#endif /* SYMBOL_HAS_LEVFVAL */
        (bloktax[bloknum].cblknum != cblknum)) {
      errorPrint ("symbolCheck: invalid diagonal block");
      return     (1);
    }

    for (bloknum ++; bloknum < cblktax[cblknum + 1].bloknum; bloknum ++) {
      if ((bloktax[bloknum].cblknum <  baseval)                      ||
          (bloktax[bloknum].cblknum >  cblkmax)                      ||
          (bloktax[bloknum].frownum <= bloktax[bloknum - 1].lrownum) ||
          (bloktax[bloknum].cblknum <  bloktax[bloknum - 1].cblknum)) {
        errorPrint ("symbolCheck: invalid block array");
        return     (1);
      }
    }
  }

  return (0);
}

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
/**   NAME       : order_check.c                           **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                This module checks the consistency of   **/
/**                orderings.                              **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 06 oct 1998     **/
/**                                 to   : 06 oct 1998     **/
/**                # Version 1.0  : from : 19 nov 2003     **/
/**                                 to   : 20 nov 2003     **/
/**                # Version 2.0  : from : 28 feb 2004     **/
/**                                 to   : 28 feb 2004     **/
/**                # Version 6.0  : from : 06 feb 2020     **/
/**                                 to   : 06 feb 2020     **/
/**                # Version 7.0  : from : 21 jan 2023     **/
/**                                 to   : 21 jan 2023     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "order.h"

/***********************************/
/*                                 */
/* The ordering handling routines. */
/*                                 */
/***********************************/

/*+ This routine checks the consistency
*** of the given ordering.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
orderCheck (
const Order * restrict const  ordeptr)
{
  INT                   baseval;                  /* Node base value            */
  INT                   vnodnnd;                  /* Based number of nodes      */
  INT                   vnodnum;                  /* Number of current node     */
  INT                   rangnum;                  /* Current column block index */
  const INT * restrict  peritax;                  /* Based access to peritab    */
  const INT * restrict  permtax;                  /* Based access to permtab    */

  if (ordeptr->cblknbr < 0) {
    errorPrint ("orderCheck: invalid nunber of column blocks");
    return (1);
  }

  baseval = ordeptr->rangtab[0];                  /* Get base value */
  if (baseval < 0) {
    errorPrint ("orderCheck: invalid vertex node base number");
    return (1);
  }

  peritax = ordeptr->peritab - baseval;           /* Set based accesses */
  vnodnnd = ordeptr->rangtab[ordeptr->cblknbr];

  for (rangnum = 0; rangnum < ordeptr->cblknbr; rangnum ++) {
    if ((ordeptr->rangtab[rangnum] <  baseval) ||
        (ordeptr->rangtab[rangnum] >= vnodnnd) ||
        (ordeptr->rangtab[rangnum] >= ordeptr->rangtab[rangnum + 1])) {
      errorPrint ("orderCheck: invalid range array");
      return (1);
    }
  }

  permtax = ordeptr->permtab - baseval;

  for (vnodnum = baseval; vnodnum < vnodnnd; vnodnum ++) {
    INT                   vnodold;

    vnodold = peritax[vnodnum];
    if ((vnodold <  baseval) ||
        (vnodold >= vnodnnd) ||
        (permtax[vnodold] != vnodnum)) {
      errorPrint ("orderCheck: invalid permutation arrays");
      return (1);
    }
  }

  return (0);
}

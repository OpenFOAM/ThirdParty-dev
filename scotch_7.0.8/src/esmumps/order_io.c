/* Copyright 2004,2007,2020,2023 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : order_io.c                              **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the input/output        **/
/**                routines for the ordering structure.    **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 04 jan 1999     **/
/**                                 to     05 jan 1999     **/
/**                # Version 2.0  : from : 28 feb 2004     **/
/**                                 to     28 feb 2004     **/
/**                # Version 7.0  : from : 21 apr 2022     **/
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

/*+ This routine reads the given ordering
*** structure from the given stream.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
orderLoad (
Order * const               ordeptr,
FILE * const                stream)
{
  INT                 versval;                    /* Version number */
  INT                 cblknbr;
  INT                 cblknum;
  INT                 vertnbr;
  INT                 vertnum;
  INT                 vertnnd;
  INT *               permtax;
  INT *               peritax;
  int                 i;

  if ((intLoad (stream, &versval) +
       intLoad (stream, &cblknbr) +
       intLoad (stream, &vertnbr) != 3) ||
      (versval != 0)                    ||
      (cblknbr > vertnbr)) {
    errorPrint ("orderLoad: bad input (1)");
    return (1);
  }

  if (((ordeptr->rangtab = (INT *) memAlloc ((cblknbr + 1) * sizeof (INT))) == NULL) ||
      ((ordeptr->permtab = (INT *) memAlloc (vertnbr       * sizeof (INT))) == NULL) ||
      ((ordeptr->peritab = (INT *) memAlloc (vertnbr       * sizeof (INT))) == NULL)) {
    errorPrint ("orderLoad: out of memory");
    orderExit  (ordeptr);
    orderInit  (ordeptr);
    return (1);
  }
  ordeptr->cblknbr = cblknbr;

  for (cblknum = 0, i = 1; (i == 1) && (cblknum <= cblknbr); cblknum ++) /* Read column-block data */
    i = intLoad (stream, &ordeptr->rangtab[cblknum]);

  for (vertnum = 0; (i == 1) && (vertnum < vertnbr); vertnum ++) /* Read direct permutation */
    i = intLoad (stream, &ordeptr->permtab[vertnum]);

  if (i != 1) {
    errorPrint ("orderLoad: bad input (2)");
    orderExit  (ordeptr);
    orderInit  (ordeptr);
    return (1);
  }

  permtax = ordeptr->permtab - ordeptr->rangtab[0];
  peritax = ordeptr->peritab - ordeptr->rangtab[0];
  for (vertnum = ordeptr->rangtab[0], vertnnd = vertnum + vertnbr; /* Build inverse permutation */
       vertnum < vertnnd; vertnum ++)
    peritax[permtax[vertnum]] = vertnum;

  return (0);
}

/*+ This routine saves the given
*** ordering structure to the
*** given stream.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
orderSave (
const Order * const         ordeptr,
FILE * const                stream)
{
  INT                 vertnbr;
  INT                 vertnum;
  INT                 cblknum;
  int                 o;

  if (ordeptr->rangtab == NULL) {
    errorPrint ("orderSave: cannot save ordering without column block data");
    return (1);
  }
  if (ordeptr->permtab == NULL) {
    errorPrint ("orderSave: cannot save ordering without direct permutation data");
    return (1);
  }

  vertnbr = ordeptr->rangtab[ordeptr->cblknbr] -  /* Get number of nodes */
            ordeptr->rangtab[0];

  if (fprintf (stream, "0\n%ld\t%ld\n",
               (long) ordeptr->cblknbr,
               (long) vertnbr) == EOF) {
    errorPrint ("orderSave: bad output (1)");
    return (1);
  }

  for (cblknum = 0, o = 1; (o == 1) && (cblknum < ordeptr->cblknbr); cblknum ++) { /* Save column-block range array */
    o = intSave (stream, ordeptr->rangtab[cblknum]);
    putc (((cblknum & 7) == 7) ? '\n' : '\t', stream);
  }
  o = intSave (stream, ordeptr->rangtab[cblknum]);
  putc ('\n', stream);

  for (vertnum = 0; (o == 1) && (vertnum < (vertnbr - 1)); vertnum ++) { /* Save direct permutation */
    o = intSave (stream, ordeptr->permtab[vertnum]);
    putc (((vertnum & 7) == 7) ? '\n' : '\t', stream);
  }
  o = intSave (stream, ordeptr->permtab[vertnum]);
  putc ('\n', stream);

  if (o != 1)
    errorPrint ("orderSave: bad output (2)");

  return (1 - o);
}

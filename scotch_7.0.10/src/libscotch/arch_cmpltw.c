/* Copyright 2007,2008,2010,2011,2014,2015,2018,2023,2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : arch_cmpltw.c                           **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module handles the weighted        **/
/**                complete graph target architecture.     **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 11 dec 2007     **/
/**                                 to   : 11 aug 2010     **/
/**                # Version 6.0  : from : 14 feb 2011     **/
/**                                 to   : 12 apr 2015     **/
/**                # Version 7.0  : from : 18 feb 2018     **/
/**                                 to   : 20 sep 2024     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "graph.h"
#include "arch.h"
#include "arch_cmpltw.h"

/******************************************/
/*                                        */
/* These are the complete graph routines. */
/*                                        */
/******************************************/

/* This routine builds a complete weighted
** graph architecture from the given load array.
** The load array is sorted recursively in such
** a way that recusive bipartitions of subdomains
** yield parts that are not too much imbalanced.
** It returns:
** - 0   : if the architecture has been successfully built.
** - !0  : on error.
*/

static
void
archCmpltwArchBuild3 (
ArchCmpltwLoad * restrict const velotab,
ArchCmpltwLoad * restrict const vesotab,
Anum                            vertnbr,
Anum                            velosum)
{
  Anum                vel0sum;
  Anum                vel1sum;
  Anum                ver0nbr;
  Anum                ver1nbr;
  Anum                ver0num;
  Anum                ver1num;
  Anum                vertnum;

  ver0num =
  ver1num = vertnbr - 1;
  vel0sum = velotab[ver0num --].veloval;
  vel1sum = 0;
  for (vertnum = ver0num; vertnum >= 0; vertnum --) {
    if (vel1sum < vel0sum) {
      vel1sum            += velotab[vertnum].veloval;
      vesotab[ver1num --] = velotab[vertnum];
    }
    else {
      vel0sum            += velotab[vertnum].veloval;
      velotab[ver0num --] = velotab[vertnum];
    }
  }

  if (vel0sum >= vel1sum) {
    ver0nbr = vertnbr - ver0num - 1;
    ver1nbr = vertnbr - ver0nbr;
    memMov (velotab,           velotab + ver1nbr, ver0nbr * sizeof (ArchCmpltwLoad));
    memCpy (velotab + ver0nbr, vesotab + ver0nbr, ver1nbr * sizeof (ArchCmpltwLoad));
  }
  else {
    Anum                velotmp;

    ver0nbr = vertnbr - ver1num - 1;
    ver1nbr = vertnbr - ver0nbr;
    memCpy (velotab, vesotab + ver1nbr, ver0nbr * sizeof (ArchCmpltwLoad));
    velotmp = vel0sum;
    vel0sum = vel1sum;
    vel1sum = velotmp;
  }

  if (ver0nbr > 2)
    archCmpltwArchBuild3 (velotab, vesotab, ver0nbr, vel0sum);
  if (ver1nbr > 2)
    archCmpltwArchBuild3 (velotab + ver0nbr, vesotab + ver0nbr, ver1nbr, vel1sum);
}

static
int
archCmpltwArchBuild2 (
ArchCmpltw * restrict const archptr)
{
  ArchCmpltwLoad * restrict vesotab;              /* Auxiliary sort array for weighted vertices */

  if (archptr->vertnbr < 3)                       /* No need to sort if less than 3 vertices */
    return (0);

  if ((vesotab = (ArchCmpltwLoad *) memAlloc (archptr->vertnbr * sizeof (ArchCmpltwLoad))) == NULL) {
    errorPrint ("archCmpltwArchBuild2: out of memory");
    memFree (archptr->velotab);
    archptr->velotab = NULL;
    return (1);
  }

  intSort2asc2 (archptr->velotab, archptr->vertnbr); /* Sort load array by both keys to be portable across sorting implementations */

  archCmpltwArchBuild3 (archptr->velotab, vesotab, archptr->vertnbr, archptr->velosum);

  memFree (vesotab);

  return (0);
}

int
archCmpltwArchBuild (
ArchCmpltw * restrict const archptr,
const Anum                  vertnbr,
const Anum * restrict const velotab)
{
  ArchCmpltwLoad *    vecwtab;
  Anum                velosum;
  Anum                vertnum;

#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchCmpltw)    > sizeof (ArchDummy)) ||
      (sizeof (ArchCmpltwDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archCmpltwArchBuild: invalid type specification");
    return (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (vertnbr <= 0) {
    errorPrint ("archCmpltwArchBuild: invalid parameters (1)");
    return (1);
  }

  if ((vecwtab = (ArchCmpltwLoad *) memAlloc (vertnbr * sizeof (ArchCmpltwLoad))) == NULL) {
    errorPrint ("archCmpltwArchBuild: out of memory");
    return (1);
  }

  for (vertnum = 0, velosum = 0; vertnum < vertnbr; vertnum ++) { /* Fill vertex load array */
    Anum                veloval;

    veloval = velotab[vertnum];
    if (veloval <= 0) {                           /* Target weights cannot be negative nor null */
      errorPrint ("archCmpltwArchBuild: invalid parameters (2)");
      memFree    (vecwtab);
      return (1);
    }

    velosum += veloval;
    vecwtab[vertnum].veloval = veloval;
    vecwtab[vertnum].vertnum = vertnum;
  }

  archptr->vertnbr = (Anum) vertnbr;
  archptr->velotab = vecwtab;
  archptr->velosum = (Anum) velosum;

  return (archCmpltwArchBuild2 (archptr));
}

/* This routine loads the weighted complete
** graph architecture.
** It returns:
** - 0   : if the architecture has been successfully read.
** - !0  : on error.
*/

int
archCmpltwArchLoad (
ArchCmpltw * restrict const archptr,
FILE * restrict const       stream)
{
  ArchCmpltwLoad *    vecwtab;
  Anum                velosum;
  long                vertnbr;
  long                vertnum;

#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchCmpltw)    > sizeof (ArchDummy)) ||
      (sizeof (ArchCmpltwDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archCmpltwArchLoad: invalid type specification");
    return (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if ((fscanf (stream, "%ld", &vertnbr) != 1) ||
      (vertnbr < 1)) {
    errorPrint ("archCmpltwArchLoad: bad input (1)");
    return (1);
  }

  if ((vecwtab = (ArchCmpltwLoad *) memAlloc (vertnbr * sizeof (ArchCmpltwLoad))) == NULL) {
    errorPrint ("archCmpltwArchLoad: out of memory");
    return (1);
  }

  for (vertnum = 0, velosum = 0; vertnum < vertnbr; vertnum ++) {
    long                velotmp;
    Anum                veloval;

    if ((fscanf (stream, "%ld", &velotmp) != 1) ||
        (velotmp < 1)) {
      errorPrint ("archCmpltwArchLoad: bad input (2)");
      return (1);
    }
    veloval = (Anum) velotmp;
    if (veloval <= 0) {                           /* Target weights cannot be negative nor null */
      errorPrint ("archCmpltwArchLoad: bad input (3)");
      memFree    (vecwtab);
      return (1);
    }

    velosum += veloval;
    vecwtab[vertnum].veloval = veloval;
    vecwtab[vertnum].vertnum = vertnum;
  }

  archptr->vertnbr = (Anum) vertnbr;
  archptr->velotab = vecwtab;
  archptr->velosum = (Anum) velosum;

  return (archCmpltwArchBuild2 (archptr));
}

/* This routine saves the weighted
** complete graph architecture.
** It returns:
** - 0   : if the architecture has been successfully written.
** - !0  : on error.
*/

int
archCmpltwArchSave (
const ArchCmpltw * const    archptr,
FILE * restrict const       stream)
{
  Anum                vertnum;

#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchCmpltw)    > sizeof (ArchDummy)) ||
      (sizeof (ArchCmpltwDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archCmpltwArchSave: invalid type specification");
    return (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (fprintf (stream, ANUMSTRING, (Anum) archptr->vertnbr) == EOF) {
    errorPrint ("archCmpltwArchSave: bad output (1)");
    return (1);
  }

  for (vertnum = 0; vertnum < archptr->vertnbr; vertnum ++) { /* For all weights to output */
    Anum                verttmp;

    for (verttmp = 0; verttmp < archptr->vertnbr; verttmp ++) { /* For all vertex indices: O(n^2) loop but we don't really care */
      if (archptr->velotab[verttmp].vertnum == vertnum) {
        if (fprintf (stream, " " ANUMSTRING, (Anum) archptr->velotab[verttmp].veloval) == EOF) {
          errorPrint ("archCmpltwArchSave: bad output (2)");
          return (1);
        }
        break;
      }
      if (verttmp == archptr->vertnbr) {
        errorPrint ("archCmpltwArchSave: internal error");
        return (1);
      }
    }
  }

  if (fprintf (stream, "\n") == EOF) {
    errorPrint ("archCmpltwArchSave: bad output (3)");
    return (1);
  }

  return (0);
}

/* This routine frees the weighted complete
** graph architecture data structures.
** It returns:
** - 0   : if the architecture has been successfully freed.
** - !0  : on error.
*/

int
archCmpltwArchFree (
ArchCmpltw * const          archptr)
{
#ifdef SCOTCH_DEBUG_ARCH1
  if ((sizeof (ArchCmpltw)    > sizeof (ArchDummy)) ||
      (sizeof (ArchCmpltwDom) > sizeof (ArchDomDummy))) {
    errorPrint ("archCmpltwArchFree: invalid type specification");
    return (1);
  }
#endif /* SCOTCH_DEBUG_ARCH1 */

  if (archptr->velotab != NULL) {
    memFree (archptr->velotab);
    archptr->velotab = NULL;
  }

  return (0);
}

/* This function returns the smallest number
** of terminal domain included in the given
** domain.
*/

ArchDomNum
archCmpltwDomNum (
const ArchCmpltw * const    archptr,
const ArchCmpltwDom * const domnptr)
{
  return (archptr->velotab[domnptr->vertmin].vertnum); /* Return vertex number */
}

/* This function returns the terminal domain associated
** with the given terminal number in the architecture.
**
** It returns:
** - 0  : if label is valid and domain has been updated.
** - 1  : if label is invalid.
** - 2  : on error.
*/

int
archCmpltwDomTerm (
const ArchCmpltw * const    archptr,
ArchCmpltwDom * const       domnptr,
const ArchDomNum            domnnum)
{
  if (domnnum < archptr->vertnbr) {               /* If valid label */
    Anum                vertnum;

    for (vertnum = 0; vertnum < archptr->vertnbr; vertnum ++) { /* Search for terminal domain index matching vertex label */
      if (archptr->velotab[vertnum].vertnum == domnnum)
        break;
    }
#ifdef SCOTCH_DEBUG_ARCH2
    if (vertnum == archptr->vertnbr) {            /* If index not found */
      errorPrint ("archCmpltwDomTerm: internal error");
      return (2);
    }
#endif /* SCOTCH_DEBUG_ARCH2 */

    domnptr->vertmin = vertnum;                   /* Set the domain */
    domnptr->vertnbr = 1;
    domnptr->veloval = archptr->velotab[vertnum].veloval;

    return (0);
  }

  return (1);                                     /* Cannot set domain */
}

/* This function returns the number of
** elements in the complete domain.
*/

Anum
archCmpltwDomSize (
const ArchCmpltw * const    archptr,
const ArchCmpltwDom * const domnptr)
{
  return (domnptr->vertnbr);
}

/* This function returns the weight of
** the complete domain.
*/

Anum
archCmpltwDomWght (
const ArchCmpltw * const    archptr,
const ArchCmpltwDom * const domnptr)
{
  return (domnptr->veloval);
}

/* This function returns the average
** distance between two complete
** subdomains.
*/

Anum
archCmpltwDomDist (
const ArchCmpltw * const    archptr,
const ArchCmpltwDom * const dom0ptr,
const ArchCmpltwDom * const dom1ptr)
{
  return (((dom0ptr->vertmin == dom1ptr->vertmin) && /* All domains are at distance 1 */
           (dom0ptr->vertnbr == dom1ptr->vertnbr)) ? 0 : 1); /* If they are different */
}

/* This function sets the biggest
** domain available for this
** architecture.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archCmpltwDomFrst (
const ArchCmpltw * const        archptr,
ArchCmpltwDom * restrict const  domnptr)
{
  domnptr->vertmin = 0;
  domnptr->vertnbr = archptr->vertnbr;
  domnptr->veloval = archptr->velosum;

  return (0);
}

/* This routine reads domain information
** from the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archCmpltwDomLoad (
const ArchCmpltw * const        archptr,
ArchCmpltwDom * restrict const  domnptr,
FILE * const                    stream)
{
  long                vertmin;
  long                vertnbr;
  Anum                vertnum;
  Anum                vertnnd;
  Anum                velosum;

  if ((fscanf (stream, "%ld%ld",
               &vertmin,
               &vertnbr) != 2) ||
      (vertnbr < 1)            ||
      (vertnbr + vertmin > (long) archptr->vertnbr)) {
    errorPrint ("archCmpltwDomLoad: bad input");
    return (1);
  }
  domnptr->vertmin = (Anum) vertmin;
  domnptr->vertnbr = (Anum) vertnbr;

  for (vertnum = domnptr->vertmin, vertnnd = vertnum + domnptr->vertnbr, velosum = 0;
       vertnum < vertnnd; vertnum ++)
    velosum += archptr->velotab[vertnum].veloval;

  domnptr->veloval += velosum;

  return (0);
}

/* This routine saves domain information
** to the given stream.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
archCmpltwDomSave (
const ArchCmpltw * const    archptr,
const ArchCmpltwDom * const domnptr,
FILE * const                stream)
{
  if (fprintf (stream, ANUMSTRING " " ANUMSTRING " ",
               (Anum) domnptr->vertmin,
               (Anum) domnptr->vertnbr) == EOF) {
    errorPrint ("archCmpltwDomSave: bad output");
    return (1);
  }

  return (0);
}

/* This function tries to split a complete
** graph domain into two subdomains.
** It returns:
** - 0  : if bipartitioning succeeded.
** - 1  : if bipartitioning could not be performed.
** - 2  : on error.
*/

int
archCmpltwDomBipart (
const ArchCmpltw * const        archptr,
const ArchCmpltwDom * const     domnptr,
ArchCmpltwDom * restrict const  dom0ptr,
ArchCmpltwDom * restrict const  dom1ptr)
{
  Anum                vertnum;
  Anum                velosum1;
  Anum                velosum2;                   /* Half of overall load sum */

  if (domnptr->vertnbr <= 1)                      /* Return if cannot bipartition more */
    return (1);

  vertnum  = domnptr->vertmin + domnptr->vertnbr - 1;
  velosum1 = (Anum) archptr->velotab[vertnum].veloval;
  velosum2 = domnptr->veloval / 2;
  for (vertnum --; vertnum > domnptr->vertmin; vertnum --) {
    Anum                velotmp;

    velotmp = velosum1 + (Anum) archptr->velotab[vertnum].veloval;
    if (velotmp > velosum2)                       /* Domain 1 is always the least loaded */
      break;
    velosum1 = velotmp;
  }

  dom0ptr->vertmin = domnptr->vertmin;            /* Bipartition vertices */
  dom1ptr->vertmin = vertnum + 1;
  dom0ptr->vertnbr = dom1ptr->vertmin - domnptr->vertmin;
  dom1ptr->vertnbr = domnptr->vertnbr - dom0ptr->vertnbr;
  dom0ptr->veloval = domnptr->veloval - velosum1;
  dom1ptr->veloval = velosum1;

  return (0);
}

/* This function checks if dom1 is
** included in dom0.
** It returns:
** - 0  : if dom1 is not included in dom0.
** - 1  : if dom1 is included in dom0.
** - 2  : on error.
*/

int
archCmpltwDomIncl (
const ArchCmpltw * const    archptr,
const ArchCmpltwDom * const dom0ptr,
const ArchCmpltwDom * const dom1ptr)
{
  if ((dom1ptr->vertmin >= dom0ptr->vertmin) &&
      ((dom1ptr->vertmin + dom1ptr->vertnbr) <= (dom0ptr->vertmin + dom0ptr->vertnbr)))
    return (1);

  return (0);
}

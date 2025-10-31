/* Copyright 2004,2007-2012,2014,2018,2019,2021,2023,2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : library_graph_map.c                     **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                Sebastien FOURESTIER (v6.0)             **/
/**                                                        **/
/**   FUNCTION   : This module is the API for the mapping  **/
/**                routines of the libSCOTCH library.      **/
/**                                                        **/
/**   DATES      : # Version 3.2  : from : 19 aug 1998     **/
/**                                 to   : 20 aug 1998     **/
/**                # Version 3.3  : from : 19 oct 1998     **/
/**                                 to   : 30 mar 1999     **/
/**                # Version 3.4  : from : 01 nov 2001     **/
/**                                 to   : 01 nov 2001     **/
/**                # Version 4.0  : from : 13 jan 2004     **/
/**                                 to   : 13 nov 2005     **/
/**                # Version 5.1  : from : 29 oct 2007     **/
/**                                 to   : 24 jul 2011     **/
/**                # Version 6.0  : from : 03 mar 2011     **/
/**                                 to   : 15 may 2018     **/
/**                # Version 7.0  : from : 07 may 2019     **/
/**                                 to   : 11 sep 2024     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "context.h"
#include "parser.h"
#include "graph.h"
#include "arch.h"
#include "arch_dist.h"
#include "mapping.h"
#include "kgraph.h"
#include "kgraph_map_st.h"
#include "library_mapping.h"
#include "scotch.h"

/************************************/
/*                                  */
/* These routines are the C API for */
/* the mapping routines.            */
/*                                  */
/************************************/

/*+ This routine initializes an API opaque
*** mapping with respect to the given source
*** graph and the locations of output parameters.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphMapInit (
const SCOTCH_Graph * const  grafptr,              /*+ Graph to map                    +*/
SCOTCH_Mapping * const      mappptr,              /*+ Mapping structure to initialize +*/
const SCOTCH_Arch * const   archptr,              /*+ Target architecture used to map +*/
SCOTCH_Num * const          parttab)              /*+ Mapping array                   +*/
{
  LibMapping * restrict const lmapptr = (LibMapping *) mappptr;

  lmapptr->flagval = LIBMAPPINGNONE;              /* No options set */
  lmapptr->grafptr = (Graph *) CONTEXTOBJECT (grafptr);
  lmapptr->archptr = (Arch *)  archptr;
  if (parttab == NULL) {
    if ((lmapptr->parttab = (Gnum *) memAlloc (lmapptr->grafptr->vertnbr * sizeof (Gnum))) == NULL) {
      errorPrint (STRINGIFY (SCOTCH_graphMapInit) ": out of memory");
      return (1);
    }
    memSet (lmapptr->parttab, 0, lmapptr->grafptr->vertnbr * sizeof (Anum)); /* All vertices mapped to first domain    */
    lmapptr->flagval |= LIBMAPPINGFREEPART;       /* The user did not provided the partition array, so we will free it */
  }
  else
    lmapptr->parttab = (Gnum *) parttab;

  return (0);
}

/*+ This routine frees an API mapping.
*** It returns:
*** - VOID  : in all cases.
+*/

void
SCOTCH_graphMapExit (
const SCOTCH_Graph * const  grafptr,
SCOTCH_Mapping * const      mappptr)
{
  LibMapping * restrict const lmapptr = (LibMapping *) mappptr;

  if (((lmapptr->flagval & LIBMAPPINGFREEPART) != 0) && /* If parttab must be freed */
      (lmapptr->parttab != NULL))                 /* And if exists                  */
    memFree (lmapptr->parttab);                   /* Free it                        */

  memSet (lmapptr, 0, sizeof (LibMapping));
}

/*+ This routine computes a mapping or a
*** remapping, with or without fixed
*** vertices, of the API mapping
*** structures given in input, with
*** respect to the given strategy.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

static
int
graphMapCompute2 (
SCOTCH_Graph * const        actgrafptr,           /*+ Graph to order                         +*/
SCOTCH_Mapping * const      mappptr,              /*+ Mapping to compute                     +*/
SCOTCH_Mapping * const      mapoptr,              /*+ Old mapping                            +*/
const double                emraval,              /*+ Edge migration ratio                   +*/
const SCOTCH_Num *          vmlotab,              /*+ Vertex migration cost array            +*/
const SCOTCH_Num            vfixnbr,              /*+ Number of fixed vertices in part array +*/
SCOTCH_Strat * const        straptr)              /*+ Mapping strategy                       +*/
{
  ArchDom             domnorg;                    /* Initial domain                       */
  Kgraph              mapgrafdat;                 /* Effective mapping graph              */
  const Strat *       mapstraptr;                 /* Pointer to mapping strategy          */
  CONTEXTDECL        (actgrafptr);
  Anum *              pfixtax;
  Gnum                crloval;                    /* Coefficient load for regular edges   */
  Gnum                cmloval;                    /* Coefficient load for migration edges */
  const Gnum *        vmlotax;                    /* Vertex migration cost array          */
  int                 o;

  o = 1;                                          /* Assume an error */

  if (CONTEXTINIT (actgrafptr)) {
    errorPrint (STRINGIFY (graphMapCompute2) ": cannot initialize context");
    return (o);
  }

  const Graph * restrict const      grafptr = (Graph *) CONTEXTGETOBJECT (actgrafptr);
  LibMapping * restrict const       lmapptr = (LibMapping *) mappptr;
  const LibMapping * restrict const lmaoptr = (LibMapping *) mapoptr;
  const Gnum                        baseval = grafptr->baseval;

#ifdef SCOTCH_DEBUG_LIBRARY1
  if (grafptr != lmapptr->grafptr) {
    errorPrint ("graphMapCompute2: output mapping does not correspond to input graph");
    goto abort;
  }
  if (lmaoptr != NULL) {
    Gnum                vertnbr;
    Gnum                vertnum;

    if (grafptr != lmaoptr->grafptr) {
      errorPrint ("graphMapCompute2: old mapping does not correspond to input graph");
      goto abort;
    }
    if (lmapptr->archptr != lmaoptr->archptr) {
      errorPrint ("graphMapCompute2: output and old mappings do not correspond to same architecture");
      goto abort;
    }

    for (vertnum = 0, vertnbr = lmaoptr->grafptr->vertnbr; vertnum < vertnbr; vertnum ++) {
      ArchDom                     domndat;

      if ((lmaoptr->parttab[vertnum] >= 0) &&
          (archDomTerm (lmapptr->archptr, &domndat, lmaoptr->parttab[vertnum]) != 0)) {
        errorPrint ("graphMapCompute2: invalid old partition");
        goto abort;
      }
    }
  }
#endif /* SCOTCH_DEBUG_LIBRARY1 */
#ifdef SCOTCH_DEBUG_LIBRARY2
  if (graphCheck (grafptr) != 0) {                /* Vertex loads can be 0 if we have fixed vertices */
    errorPrint ("graphMapCompute2: invalid input graph");
    goto abort;
  }
#endif /* SCOTCH_DEBUG_LIBRARY2 */

  archDomFrst (lmapptr->archptr, &domnorg);       /* Compute initial domain to map to */

  if (*((Strat **) straptr) == NULL)              /* Set default mapping strategy if necessary */
    SCOTCH_stratGraphMapBuild (straptr, SCOTCH_STRATDEFAULT, archDomSize (lmapptr->archptr, &domnorg), 0.01);

  mapstraptr = *((Strat **) straptr);
#ifdef SCOTCH_DEBUG_LIBRARY1
  if (mapstraptr->tablptr != &kgraphmapststratab) {
    errorPrint ("graphMapCompute2: not a graph mapping strategy");
    goto abort;
  }
#endif /* SCOTCH_DEBUG_LIBRARY1 */

  if (vfixnbr > 0) {                              /* We have fixed vertices */
#ifdef SCOTCH_DEBUG_LIBRARY1
    ArchDom             domndat;
    Gnum                vertnbr;
    Gnum                vertnum;

    if (lmapptr->parttab == NULL) {               /* We must have fixed vertex information */
      errorPrint ("graphMapCompute2: missing output mapping part array");
      goto abort;
    }
    for (vertnum = 0, vertnbr = lmapptr->grafptr->vertnbr; vertnum < vertnbr; vertnum ++) {
      if ((lmapptr->parttab[vertnum] >= 0) &&
          (archDomTerm (lmapptr->archptr, &domndat, lmapptr->parttab[vertnum]) != 0)) {
        errorPrint ("graphMapCompute2: invalid fixed partition");
        goto abort;
      }
    }
#endif /* SCOTCH_DEBUG_LIBRARY1 */
    if (vfixnbr >= grafptr->vertnbr) {            /* If all fixed vertices, do nothing */
      o = 0;
      goto abort;
    }

    pfixtax = lmapptr->parttab - baseval;
  }
  else
    pfixtax = NULL;

  if (lmaoptr != NULL) {                          /* We are doing a repartitioning */
    Gnum                        numeval;
    Gnum                        denoval;

    vmlotax = (vmlotab != NULL) ? vmlotab - baseval : NULL;
    numeval = (INT) ((emraval * 100.0) + 0.5);
    denoval = intGcd (numeval, 100);
    cmloval = numeval / denoval;
    crloval = 100     / denoval;
  }
  else {
    vmlotax = NULL;
    cmloval =
    crloval = 1;
  }

  if (kgraphInit (&mapgrafdat, grafptr, lmapptr->archptr, &domnorg, vfixnbr, pfixtax, crloval, cmloval, vmlotax) != 0)
    goto abort;

  if (lmaoptr != NULL) {                          /* If we are doing a repartitioning, fill old mapping structure */
    if (mapAlloc (&mapgrafdat.r.m) != 0) {        /* Allocate part and domain arrays                              */
      errorPrint ("kgraphInit: cannot initialize remapping (1)");
      kgraphExit (&mapgrafdat);
      goto abort;
    }

    memSet (mapgrafdat.r.m.parttax + baseval, ~0, grafptr->vertnbr * sizeof (Anum)); /* Pre-set unknown vertex domains */
    mapgrafdat.r.m.flagval |= MAPPINGINCOMPLETE;  /* Mapping may contain incomplete information                        */

    if (mapBuild (&mapgrafdat.r.m, lmaoptr->parttab - baseval) != 0) { /* Merge old part information to incomplete mapping */
      errorPrint ("kgraphInit: cannot initialize remapping (2)");
      kgraphExit (&mapgrafdat);
      goto abort;
    }
  }

  mapgrafdat.contptr = CONTEXTGETDATA (actgrafptr);

  o = kgraphMapSt (&mapgrafdat, mapstraptr);
  mapTerm (&mapgrafdat.m, lmapptr->parttab - baseval); /* Propagate mapping result to part array */

  kgraphExit (&mapgrafdat);

abort:
  CONTEXTEXIT (actgrafptr);
  return (o);
}

/*+ This routine computes a mapping
*** of the API mapping structure with
*** respect to the given strategy.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphMapCompute (
SCOTCH_Graph * const        grafptr,              /*+ Graph to order     +*/
SCOTCH_Mapping * const      mappptr,              /*+ Mapping to compute +*/
SCOTCH_Strat * const        straptr)              /*+ Mapping strategy   +*/
{
  return (graphMapCompute2 (grafptr, mappptr, NULL, 1, NULL, 0, straptr));
}

/*+ This routine computes a mapping
*** with fixed vertices of the API
*** mapping structure with respect
*** to the given strategy.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphMapFixedCompute (
SCOTCH_Graph * const        grafptr,              /*+ Graph to order     +*/
SCOTCH_Mapping * const      mappptr,              /*+ Mapping to compute +*/
SCOTCH_Strat * const        straptr)              /*+ Mapping strategy   +*/
{
  return (SCOTCH_graphRemapFixedCompute (grafptr, mappptr, NULL, 1, NULL, straptr));
}

/*+ This routine computes a remapping
*** of the API mapping structure with
*** respect to the given strategy.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphRemapCompute (
SCOTCH_Graph * const        grafptr,              /*+ Graph to order              +*/
SCOTCH_Mapping * const      mappptr,              /*+ Mapping to compute          +*/
SCOTCH_Mapping * const      mapoptr,              /*+ Old mapping                 +*/
const double                emraval,              /*+ Edge migration ratio        +*/
const SCOTCH_Num *          vmlotab,              /*+ Vertex migration cost array +*/
SCOTCH_Strat * const        straptr)              /*+ Mapping strategy            +*/
{
  return (graphMapCompute2 (grafptr, mappptr, mapoptr, emraval, vmlotab, 0, straptr));
}

/*+ This routine computes a remapping
*** with fixed vertices of the API
*** mapping structure with respect
*** to the given strategy.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphRemapFixedCompute (
SCOTCH_Graph * const        grafptr,              /*+ Graph to order              +*/
SCOTCH_Mapping * const      mappptr,              /*+ Mapping to compute          +*/
SCOTCH_Mapping * const      mapoptr,              /*+ Old mapping                 +*/
const double                emraval,              /*+ Edge migration ratio        +*/
const SCOTCH_Num *          vmlotab,              /*+ Vertex migration cost array +*/
SCOTCH_Strat * const        straptr)              /*+ Mapping strategy            +*/
{
  Gnum                vfixnbr;
  Gnum                vertnbr;
  Gnum                vertnum;

  const Anum * restrict const pfixtab = ((LibMapping *) mappptr)->parttab;

  for (vertnum = 0, vertnbr = ((Graph *) CONTEXTOBJECT (grafptr))->vertnbr, vfixnbr = 0; /* Compute number of fixed vertices */
       vertnum < vertnbr; vertnum ++) {
    if (pfixtab[vertnum] != ~0)
      vfixnbr ++;
  }

  return (graphMapCompute2 (grafptr, mappptr, mapoptr, emraval, vmlotab, vfixnbr, straptr));
}

/*+ This routine computes a mapping of the
*** given graph structure onto the given
*** target architecture with respect to the
*** given strategy.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphMap (
SCOTCH_Graph * const        grafptr,              /*+ Graph to map        +*/
const SCOTCH_Arch * const   archptr,              /*+ Target architecture +*/
SCOTCH_Strat * const        straptr,              /*+ Mapping strategy    +*/
SCOTCH_Num * const          parttab)              /*+ Partition array     +*/
{
  SCOTCH_Mapping      mappdat;
  int                 o;

  SCOTCH_graphMapInit (grafptr, &mappdat, archptr, parttab);
  o = SCOTCH_graphMapCompute (grafptr, &mappdat, straptr);
  SCOTCH_graphMapExit (grafptr, &mappdat);

  return (o);
}

/*+ This routine computes a mapping of the
*** given graph structure onto the given
*** target architecture with respect to the
*** given strategy and the fixed vertices in
*** maptab.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphMapFixed (
SCOTCH_Graph * const        grafptr,              /*+ Graph to map        +*/
const SCOTCH_Arch * const   archptr,              /*+ Target architecture +*/
SCOTCH_Strat * const        straptr,              /*+ Mapping strategy    +*/
SCOTCH_Num * const          parttab)              /*+ Partition array     +*/
{
  SCOTCH_Mapping      mappdat;
  int                 o;

  SCOTCH_graphMapInit (grafptr, &mappdat, archptr, parttab);
  o = SCOTCH_graphMapFixedCompute (grafptr, &mappdat, straptr);
  SCOTCH_graphMapExit (grafptr, &mappdat);

  return (o);
}

/*+ This routine computes a remapping of the
*** given graph structure onto the given
*** target architecture with respect to the
*** given strategy.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphRemap (
SCOTCH_Graph * const        grafptr,              /*+ Graph to map                +*/
const SCOTCH_Arch * const   archptr,              /*+ Target architecture         +*/
SCOTCH_Num * const          parotab,              /*+ Old partition array         +*/
const double                emraval,              /*+ Edge migration ratio        +*/
const SCOTCH_Num *          vmlotab,              /*+ Vertex migration cost array +*/
SCOTCH_Strat * const        straptr,              /*+ Mapping strategy            +*/
SCOTCH_Num * const          parttab)              /*+ Partition array             +*/
{
  SCOTCH_Mapping      mappdat;
  SCOTCH_Mapping      mapodat;
  int                 o;

  SCOTCH_graphMapInit (grafptr, &mappdat, archptr, parttab);
  SCOTCH_graphMapInit (grafptr, &mapodat, archptr, parotab);
  o = SCOTCH_graphRemapCompute (grafptr, &mappdat, &mapodat, emraval, vmlotab, straptr);
  SCOTCH_graphMapExit (grafptr, &mapodat);
  SCOTCH_graphMapExit (grafptr, &mappdat);

  return (o);
}

/*+ This routine computes a remapping of the
*** given graph structure onto the given
*** target architecture with respect to the
*** given strategy and the fixed vertices in
*** maptab.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphRemapFixed (
SCOTCH_Graph * const        grafptr,              /*+ Graph to map                +*/
const SCOTCH_Arch * const   archptr,              /*+ Target architecture         +*/
SCOTCH_Num * const          parotab,              /*+ Old partition array         +*/
const double                emraval,              /*+ Edge migration ratio        +*/
const SCOTCH_Num *          vmlotab,              /*+ Vertex migration cost array +*/
SCOTCH_Strat * const        straptr,              /*+ Mapping strategy            +*/
SCOTCH_Num * const          parttab)              /*+ Partition array             +*/
{
  SCOTCH_Mapping      mappdat;
  SCOTCH_Mapping      mapodat;
  int                 o;

  SCOTCH_graphMapInit (grafptr, &mappdat, archptr, parttab);
  SCOTCH_graphMapInit (grafptr, &mapodat, archptr, parotab);
  o = SCOTCH_graphRemapFixedCompute (grafptr, &mappdat, &mapodat, emraval, vmlotab, straptr);
  SCOTCH_graphMapExit (grafptr, &mapodat);
  SCOTCH_graphMapExit (grafptr, &mappdat);

  return (o);
}

/*+ This routine computes a partition of
*** the given graph structure with respect
*** to the given strategy.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphPart (
SCOTCH_Graph * const        grafptr,              /*+ Graph to map     +*/
const SCOTCH_Num            partnbr,              /*+ Number of parts  +*/
SCOTCH_Strat * const        straptr,              /*+ Mapping strategy +*/
SCOTCH_Num * const          parttab)              /*+ Partition array  +*/
{
  SCOTCH_Arch         archdat;
  int                 o;

  SCOTCH_archInit  (&archdat);
  SCOTCH_archCmplt (&archdat, partnbr);
  o = SCOTCH_graphMap (grafptr, &archdat, straptr, parttab);
  SCOTCH_archExit (&archdat);

  return (o);
}

/*+ This routine computes a partition of
*** the given graph structure with respect
*** to the given strategy and the fixed
*** vertices in maptab.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphPartFixed (
SCOTCH_Graph * const        grafptr,              /*+ Graph to map     +*/
const SCOTCH_Num            partnbr,              /*+ Number of parts  +*/
SCOTCH_Strat * const        straptr,              /*+ Mapping strategy +*/
SCOTCH_Num * const          parttab)              /*+ Partition array  +*/
{
  SCOTCH_Arch         archdat;
  int                 o;

  SCOTCH_archInit  (&archdat);
  SCOTCH_archCmplt (&archdat, partnbr);
  o = SCOTCH_graphMapFixed (grafptr, &archdat, straptr, parttab);
  SCOTCH_archExit (&archdat);

  return (o);
}

/*+ This routine computes a repartitionning
*** of the given graph structure with
*** respect to the given strategy.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphRepart (
SCOTCH_Graph * const        grafptr,              /*+ Graph to map                +*/
const SCOTCH_Num            partnbr,              /*+ Number of parts             +*/
SCOTCH_Num * const          parotab,              /*+ Old partition array         +*/
const double                emraval,              /*+ Edge migration ratio        +*/
const SCOTCH_Num * const    vmlotab,              /*+ Vertex migration cost array +*/
SCOTCH_Strat * const        straptr,              /*+ Mapping strategy            +*/
SCOTCH_Num * const          parttab)              /*+ Partition array             +*/
{
  SCOTCH_Arch         archdat;
  int                 o;

  SCOTCH_archInit  (&archdat);
  SCOTCH_archCmplt (&archdat, partnbr);
  o = SCOTCH_graphRemap (grafptr, &archdat, parotab, emraval, vmlotab, straptr, parttab);
  SCOTCH_archExit (&archdat);

  return (o);
}

/*+ This routine computes a repartitionning
*** of the given graph structure with
*** respect to the given strategy and the
*** fixed vertices in maptab.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
SCOTCH_graphRepartFixed (
SCOTCH_Graph * const        grafptr,              /*+ Graph to map                +*/
const SCOTCH_Num            partnbr,              /*+ Number of parts             +*/
SCOTCH_Num * const          parotab,              /*+ Old partition array         +*/
const double                emraval,              /*+ Edge migration ratio        +*/
const SCOTCH_Num *          vmlotab,              /*+ Vertex migration cost array +*/
SCOTCH_Strat * const        straptr,              /*+ Mapping strategy            +*/
SCOTCH_Num * const          parttab)              /*+ Partition array             +*/
{
  SCOTCH_Arch         archdat;
  int                 o;

  SCOTCH_archInit  (&archdat);
  SCOTCH_archCmplt (&archdat, partnbr);
  o = SCOTCH_graphRemapFixed (grafptr, &archdat, parotab, emraval, vmlotab, straptr, parttab);
  SCOTCH_archExit (&archdat);

  return (o);
}

/*+ This routine parses the given
*** mapping strategy.
*** It returns:
*** - 0   : if string successfully scanned.
*** - !0  : on error.
+*/

int
SCOTCH_stratGraphMap (
SCOTCH_Strat * const        straptr,
const char * const          string)
{
  if (*((Strat **) straptr) != NULL)
    stratExit (*((Strat **) straptr));

  if ((*((Strat **) straptr) = stratInit (&kgraphmapststratab, string)) == NULL) {
    errorPrint (STRINGIFY (SCOTCH_stratGraphMap) ": error in mapping strategy");
    return (1);
  }

  return (0);
}

/*+ This routine provides predefined
*** mapping strategies.
*** It returns:
*** - 0   : if string successfully initialized.
*** - !0  : on error.
+*/

int
SCOTCH_stratGraphMapBuild (
SCOTCH_Strat * const        straptr,              /*+ Strategy to create              +*/
const SCOTCH_Num            flagval,              /*+ Desired characteristics         +*/
const SCOTCH_Num            partnbr,              /*+ Number of expected parts/size   +*/
const double                kbalval)              /*+ Desired imbalance ratio         +*/
{
  char                bufftab[8192];              /* Should be enough */
  char                bbaltab[64];
  char                kbaltab[64];
  char                kmovtab[64];
  char                mvrttab[64];
  const char *        difkptr;
  const char *        difsptr;
  const char *        exasptr;
  const char *        exaxptr;

  sprintf (bbaltab, "%lf", kbalval);
  sprintf (kbaltab, "%lf", kbalval);
  sprintf (kmovtab, GNUMSTRING, (Gnum) (((flagval & SCOTCH_STRATQUALITY) != 0) ? 200 : 80));
  sprintf (mvrttab, GNUMSTRING, (Gnum) (MAX ((20 * partnbr), 10000)));

  strcpy (bufftab, ((flagval & SCOTCH_STRATRECURSIVE) != 0)
          ? "<RECU>"                              /* Use only the recursive bipartitioning framework */
          : "m{vert=<MVRT>,low=<RECU>,asc=b{bnd=<DIFK>f{bal=<KBAL>,move=<KMOV>},org=f{bal=<KBAL>,move=<KMOV>}}}<EXAX>");
  stringSubst (bufftab, "<RECU>", "r{job=t,map=t,poli=S,bal=<KBAL>,sep=<BSEP><EXAS>}");
  stringSubst (bufftab, "<BSEP>", ((flagval & SCOTCH_STRATQUALITY) != 0) ?  "<BSEQ>|<BSEQ>|<BSEQ>" :  "<BSEQ>|<BSEQ>");
  stringSubst (bufftab, "<BSEQ>", "m{vert=120,low=h{pass=10}f{bal=<BBAL>,move=120},asc=b{bnd=<DIFS>f{bal=<BBAL>,move=120},org=f{bal=<BBAL>,move=120}}}");

  if ((flagval & SCOTCH_STRATSAFETY) != 0)
    difsptr = "";
  else
    difsptr = "(d{pass=40}|)";
  difkptr = "d{pass=40}";

  if ((flagval & SCOTCH_STRATBALANCE) != 0) {
    exasptr = "f{bal=<KBAL>}";
    exaxptr = "x{bal=<KBAL>}f{bal=<KBAL>,move=<KMOV>}";
  }
  else {
    exasptr = "";
    exaxptr = "";
  }

  stringSubst (bufftab, "<MVRT>", mvrttab);
  stringSubst (bufftab, "<EXAX>", exaxptr);
  stringSubst (bufftab, "<EXAS>", exasptr);
  stringSubst (bufftab, "<DIFS>", difsptr);
  stringSubst (bufftab, "<DIFK>", difkptr);
  stringSubst (bufftab, "<KMOV>", kmovtab);
  stringSubst (bufftab, "<KBAL>", kbaltab);
  stringSubst (bufftab, "<BBAL>", bbaltab);

  if (SCOTCH_stratGraphMap (straptr, bufftab) != 0) {
    errorPrint (STRINGIFY (SCOTCH_stratGraphMapBuild) ": error in sequential mapping strategy");
    return (1);
  }

  return (0);
}

/*+ This routine provides predefined
*** clustering strategies.
*** It returns:
*** - 0   : if string successfully initialized.
*** - !0  : on error.
+*/

int
SCOTCH_stratGraphClusterBuild (
SCOTCH_Strat * const        straptr,              /*+ Strategy to create      +*/
const SCOTCH_Num            flagval,              /*+ Desired characteristics +*/
const SCOTCH_Num            pwgtval,              /*+ Threshold part weight   +*/
const double                densval,              /*+ Threshold density value +*/
const double                bbalval)              /*+ Maximum imbalance ratio +*/
{
  char                bufftab[8192];              /* Should be enough */
  char                bbaltab[32];
  char                pwgttab[32];
  char                denstab[32];
  char *              difsptr;
  char *              exasptr;

  sprintf (bbaltab, "%lf", bbalval);
  sprintf (denstab, "%lf", densval);
  sprintf (pwgttab, GNUMSTRING, pwgtval);

  strcpy (bufftab, "r{job=u,map=t,poli=L,sep=/((load><PWGT>)&!(edge>vert*<DENS>*(vert-1)))?(<BIPA>m{vert=80,low=h{pass=10}f{bal=<BBAL>,move=80},asc=b{bnd=<DIFS>f{bal=<BBAL>,move=80},org=f{bal=<BBAL>,move=80}}})<EXAS>;}");
  stringSubst (bufftab, "<BIPA>", ((flagval & SCOTCH_STRATSPEED) != 0) ? ""
               : "m{vert=80,low=h{pass=10}f{bal=<BBAL>,move=80},asc=b{bnd=<DIFS>f{bal=<BBAL>,move=80},org=f{bal=<BBAL>,move=80}}}|");

  if ((flagval & SCOTCH_STRATBALANCE) != 0)
    exasptr = "f{bal=0}";
  else
    exasptr = "";

  if ((flagval & SCOTCH_STRATSAFETY) != 0)
    difsptr = "";
  else
    difsptr = "(d{pass=40}|)";

  stringSubst (bufftab, "<EXAS>", exasptr);
  stringSubst (bufftab, "<DIFS>", difsptr);
  stringSubst (bufftab, "<BBAL>", bbaltab);
  stringSubst (bufftab, "<DENS>", denstab);
  stringSubst (bufftab, "<PWGT>", pwgttab);

  if (SCOTCH_stratGraphMap (straptr, bufftab) != 0) {
    errorPrint (STRINGIFY (SCOTCH_stratGraphClusterBuild) ": error in sequential mapping strategy");
    return (1);
  }

  return (0);
}

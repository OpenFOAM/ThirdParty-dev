/* Copyright 2019,2020 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_hamf.c                             **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module tests the original Fortran  **/
/**                HAMF code provided by the MUMPS team.   **/
/**                                                        **/
/**   DATES      : # Version 6.1  : from : 24 nov 2019     **/
/**                                 to   : 10 feb 2020     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "common.h"
#include "module.h"
#include "graph.h"
#include "arch.h"
#include "bgraph.h"
#include "bgraph_bipart_gp.h"
#include "hgraph.h"
#include "scotch.h"

/*
**  The function prototypes.
*/

void                        mumps_hamf4_ (int *, int *, int *, int *, int64_t *, int64_t *, int64_t *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);

/*********************/
/*                   */
/* The test routine. */
/*                   */
/*********************/

int
test_hamf_graph1 (
Hgraph *                    grafptr)
{
  Gnum                vertnum;
  Gnum                vertadj;
  Gnum                vertnew;
  int64_t             edgenew;
  int                 norig;
  int                 n;
  int                 nbelts;
  int                 nbbuck;
  int                 ncmpa;
  int64_t             iwlen;
  int64_t             pfree;
  int *               degtab;
  int *               lentab;
  int *               iwtab;
  int *               nvtab;
  int *               elentab;
  int *               headtab;
  int *               nexttab;
  int *               lasttab;
  int *               wtab;
  int *               wftab;
  int64_t *           petab;

  const Gnum * restrict const verttax = grafptr->s.verttax;
  const Gnum * restrict const vendtax = grafptr->s.vendtax;
  const Gnum * restrict const velotax = grafptr->s.velotax;
  const Gnum * restrict const edgetax = grafptr->s.edgetax;

  n      = grafptr->s.vertnbr;
  norig  = grafptr->s.velosum;
  nbelts = 0;
  nbbuck = norig * 2;
  iwlen  = (int64_t) ((double) grafptr->s.edgenbr * 1.2) + 32;

  petab   = malloc (n * sizeof (int64_t));
  lentab  = malloc (n * sizeof (int));
  nvtab   = malloc (n * sizeof (int));
  elentab = malloc (n * sizeof (int));
  lasttab = malloc (n * sizeof (int));
  degtab  = malloc (n * sizeof (int));
  wftab   = malloc (n * sizeof (int));
  nexttab = malloc (n * sizeof (int));
  wtab    = malloc (n * sizeof (int));
  headtab = malloc ((nbbuck + 2) * sizeof (int));
  iwtab   = malloc (iwlen * sizeof (int));

  int64_t * restrict const  petax   = petab   - 1; /* Base HAMF arrays at base 1 */
  int * restrict const      iwtax   = iwtab   - 1;
  int * restrict const      lentax  = lentab  - 1;
  int * restrict const      nvtax   = nvtab   - 1;
  int * restrict const      elentax = elentab - 1;

  vertadj = 1 - grafptr->s.baseval;
  for (vertnum = grafptr->s.baseval, vertnew = edgenew = 1; /* Process non-halo vertices */
       vertnum < grafptr->vnohnnd; vertnum ++, vertnew ++) {
    int                       degrval;
    int                       edgenum;

    degrval = vendtax[vertnum] - verttax[vertnum];
    petax[vertnew]   = edgenew;
    lentax[vertnew]  = degrval;
    elentax[vertnew] = 0;
    nvtax[vertnew]   = (velotax != NULL) ? velotax[vertnum] : 1;

    for (edgenum = verttax[vertnum]; edgenum < vendtax[vertnum]; edgenum ++, edgenew ++)
      iwtax[edgenew] = edgetax[edgenum] + vertadj;
  }
  for ( ; vertnum < grafptr->s.vertnnd; vertnum ++, vertnew ++) { /* Process halo vertices */
    int                       degrval;
    int                       edgenum;

    degrval = grafptr->s.verttax[vertnum] - grafptr->s.vendtax[vertnum]; /* Negative degree */
    petax[vertnew]   = edgenew;
    lentax[vertnew]  = (degrval != 0) ? degrval : (-1 - grafptr->s.vertnbr);
    elentax[vertnew] = 0;
    nvtax[vertnew]   = (velotax != NULL) ? velotax[vertnum] : 1;

    for (edgenum = grafptr->s.verttax[vertnum];
         edgenum < grafptr->s.vendtax[vertnum]; edgenum ++, edgenew ++)
      iwtax[edgenew] = edgetax[edgenum] + vertadj;
  }

  pfree = edgenew;                             /* Set index to first free area */

  mumps_hamf4_ (&norig, &n, &nbelts, &nbbuck,
                &iwlen, petab, &pfree, lentab, iwtab, nvtab, elentab,
                lasttab, &ncmpa, degtab, wftab, nexttab, wtab, headtab);

  if (ncmpa < 0) {
    errorPrint ("mumps_hamf4: internal error");
    return     (1);
  }

  memFree (iwtab);
  memFree (headtab);
  memFree (wtab);
  memFree (nexttab);
  memFree (wftab);
  memFree (degtab);
  memFree (lasttab);
  memFree (elentab);
  memFree (nvtab);
  memFree (lentab);
  memFree (petab);

  return (0);
}

int
test_hamf_graph2 (
SCOTCH_Graph *              grafptr)
{
  Arch                archdat;
  ArchDom             domntab[3];
  Gnum                vflowgttab[2] = { 0, 0 };
  Bgraph              bipgrafdat;
  BgraphBipartGpParam bipstradat = { 3 };
  Gnum                bipvertnum;
  Gnum *              indvnumtab;
  Gnum                indvnumnbr;
  Hgraph              halgrafdat;
  Hgraph              indgrafdat;
  int                 o;

  SCOTCH_archCmplt ((SCOTCH_Arch *) &archdat, 2);
  archDomFrst (&archdat, &domntab[0]);
  archDomBipart (&archdat, &domntab[0], &domntab[1], &domntab[2]);

  bgraphInit (&bipgrafdat, (Graph *) grafptr, &archdat, &domntab[1], vflowgttab);

  if (bgraphBipartGp (&bipgrafdat, &bipstradat) != 0) {
    errorPrint ("test_hamf_graph2: cannot compute bipartition");
    bgraphExit (&bipgrafdat);
    return (1);
  }

  if ((indvnumtab = memAlloc (bipgrafdat.compsize0 * sizeof (Gnum))) == NULL) {
    errorPrint ("test_hamf_graph2: out of memory (3)");
    return (1);
  }
  for (bipvertnum = bipgrafdat.s.baseval, indvnumnbr = 0; bipvertnum < bipgrafdat.s.vertnnd; bipvertnum ++) {
    if (bipgrafdat.parttax[bipvertnum] == 0)
      indvnumtab[indvnumnbr ++] = bipvertnum;
  }
  if (indvnumnbr != bipgrafdat.compsize0) {
    errorPrint ("test_hamf_graph2: internal error");
    return (1);
  }

  halgrafdat.s = bipgrafdat.s;                    /* Create cloned halo graph from initial graph */
  halgrafdat.s.flagval &= ~GRAPHFREETABS;
  halgrafdat.vnohnbr    = halgrafdat.s.vertnbr;
  halgrafdat.vnohnnd    = halgrafdat.s.vertnnd;
  halgrafdat.vnhdtax    = halgrafdat.s.vendtax;
  halgrafdat.vnlosum    = halgrafdat.s.velosum;
  halgrafdat.enohnbr    = halgrafdat.s.edgenbr;
  halgrafdat.enlosum    = halgrafdat.s.edlosum;
  halgrafdat.levlnum    = 0;

  bgraphExit (&bipgrafdat);

  if (hgraphInduceList (&halgrafdat, indvnumnbr, indvnumtab, halgrafdat.s.vertnbr - indvnumnbr, &indgrafdat) != 0) {
    memFree (indvnumtab);
    return  (1);
  }
  memFree (indvnumtab);

  o = test_hamf_graph1 (&indgrafdat);

  hgraphExit (&indgrafdat);

  return (o);
}

int
test_hamf_graph3 (
SCOTCH_Graph *              grafptr)
{
  Hgraph              halgrafdat;
  int                 o;

  halgrafdat.s = *((Graph *) grafptr);            /* Create cloned halo graph from initial graph */
  halgrafdat.s.flagval &= ~GRAPHFREETABS;
  halgrafdat.vnohnbr    = halgrafdat.s.vertnbr;
  halgrafdat.vnohnnd    = halgrafdat.s.vertnnd;
  halgrafdat.vnhdtax    = halgrafdat.s.vendtax;
  halgrafdat.vnlosum    = halgrafdat.s.velosum;
  halgrafdat.enohnbr    = halgrafdat.s.edgenbr;
  halgrafdat.enlosum    = halgrafdat.s.edlosum;
  halgrafdat.levlnum    = 0;

  o = test_hamf_graph1 (&halgrafdat);

  hgraphExit (&halgrafdat);

  return (o);
}

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
  SCOTCH_Graph        grafdat;
  FILE *              fileptr;
  int                 flagval;

  flagval = 0;
  if ((argc < 2) || (argc > 3)) {
    SCOTCH_errorPrint ("usage: %s [-i] graph_file", argv[0]);
    exit (EXIT_FAILURE);
  }
  if (strcmp (argv[1], "-i") == 0)
    flagval = 1;

  if (SCOTCH_graphInit (&grafdat) != 0) {
    SCOTCH_errorPrint ("main: cannot initialize graph (1)");
    exit (EXIT_FAILURE);
  }

  if ((fileptr = fopen (argv[1 + flagval], "r")) == NULL) {
    SCOTCH_errorPrint ("main: cannot open file");
    exit (EXIT_FAILURE);
  }

  if (SCOTCH_graphLoad (&grafdat, fileptr, -1, 0) != 0) { /* Read source graph */
    SCOTCH_errorPrint ("main: cannot load graph");
    exit (EXIT_FAILURE);
  }

  fclose (fileptr);

  if (flagval != 0)
    test_hamf_graph2 (&grafdat);
  else
    test_hamf_graph3 (&grafdat);

  SCOTCH_graphExit (&grafdat);

  exit (EXIT_SUCCESS);
}

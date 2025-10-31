/* Copyright 2008-2011,2023,2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : kdgraph_map_st.c                        **/
/**                                                        **/
/**   AUTHOR     : Jun-Ho HER (v6.0)                       **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module contains the strategy and   **/
/**                method tables for the parallel          **/
/**                multi-way static mapping routines.      **/
/**                                                        **/
/**   DATES      : # Version 5.1  : from : 16 jun 2008     **/
/**                                 to   : 14 apr 2011     **/
/**                # Version 7.0  : from : 20 jan 2023     **/
/**                                 to   : 07 nov 2024     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "parser.h"
#include "graph.h"
#include "dgraph.h"
#include "dgraph_coarsen.h"
#include "arch.h"
#include "mapping.h"
#include "dmapping.h"
#include "bdgraph.h"
#include "bdgraph_bipart_st.h"
#include "kgraph.h"
#include "kgraph_map_st.h"
#include "kdgraph.h"
#include "kdgraph_map_rb.h"
#include "kdgraph_map_st.h"

/*
**  The static and global variables.
*/

static union {
  KdgraphMapRbParam         param;
  StratNodeMethodData       padding;
} kdgraphmapstdefaultrb = { { &stratdummy, &stratdummy, 0.05 } };

static StratMethodTab       kdgraphmapstmethtab[] = { /* Mapping methods array */
                              { KDGRAPHMAPSTMETHRB, "r",  (StratMethodFunc) kdgraphMapRb, &kdgraphmapstdefaultrb },
                              { -1,                 NULL, (StratMethodFunc) NULL,         NULL } };

static StratParamTab        kdgraphmapstparatab[] = { /* Method parameter list */
                              { KDGRAPHMAPSTMETHRB,  STRATPARAMSTRAT,  "sep",
                                (byte *) &kdgraphmapstdefaultrb.param,
                                (byte *) &kdgraphmapstdefaultrb.param.stratsep,
                                (void *) &bdgraphbipartststratab },
                              { KDGRAPHMAPSTMETHRB,  STRATPARAMSTRAT,  "seq",
                                (byte *) &kdgraphmapstdefaultrb.param,
                                (byte *) &kdgraphmapstdefaultrb.param.stratseq,
                                (void *) &kgraphmapststratab },
                              { KDGRAPHMAPSTMETHRB,  STRATPARAMDOUBLE, "bal",
                                (byte *) &kdgraphmapstdefaultrb.param,
                                (byte *) &kdgraphmapstdefaultrb.param.kbalval,
                                NULL },
                              { KDGRAPHMAPSTMETHNBR, STRATPARAMINT,    NULL,
                                NULL, NULL, NULL } };

static StratParamTab        kdgraphmapstcondtab[] = { /* Graph condition parameter table */
                              { STRATNODENBR,        STRATPARAMINT,    NULL,
                                NULL, NULL, NULL } };

StratTab                    kdgraphmapststratab = { /* Strategy tables for graph mapping methods */
                              kdgraphmapstmethtab,
                              kdgraphmapstparatab,
                              kdgraphmapstcondtab };

/****************************************/
/*                                      */
/* This is the generic mapping routine. */
/*                                      */
/****************************************/

/* This routine computes the given
** mapping according to the given
** strategy.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
kdgraphMapSt (
Kdgraph * restrict const      grafptr,            /*+ Mapping graph    +*/
Kdmapping * restrict const    mappptr,            /*+ Dynamic mapping  +*/
const Strat * restrict const  straptr)            /*+ Mapping strategy +*/
{
  StratTest           testdat;
  int                 o;

#ifdef SCOTCH_DEBUG_KDGRAPH2
  if (sizeof (Gnum) != sizeof (INT)) {
    errorPrint ("kdgraphMapSt: invalid type specification for parser variables");
    return (1);
  }
  if ((sizeof (KdgraphMapRbParam) > sizeof (StratNodeMethodData))) {
    errorPrint ("kdgraphMapSt: invalid type specification");
    return (1);
  }
#endif /* SCOTCH_DEBUG_KDGRAPH2 */
#ifdef SCOTCH_DEBUG_KDGRAPH1
  if ((straptr->tablptr != &kdgraphmapststratab) &&
      (straptr          != &stratdummy)) {
    errorPrint ("kdgraphMapSt: invalid parameter (1)");
    return (1);
  }
#endif /* SCOTCH_DEBUG_KDGRAPH1 */

  o = 0;
  switch (straptr->typeval) {
    case STRATNODECONCAT :
      o = kdgraphMapSt (grafptr, mappptr, straptr->data.concdat.stratab[0]); /* Apply first strategy          */
      if (o == 0)                                 /* If it worked all right                                   */
        o |= kdgraphMapSt (grafptr, mappptr, straptr->data.concdat.stratab[1]); /* Then apply second strategy */
      break;
    case STRATNODECOND :
      o = stratTestEval (straptr->data.conddat.testptr, &testdat, (void *) grafptr); /* Evaluate expression */
      if (o == 0) {                               /* If evaluation was correct */
#ifdef SCOTCH_DEBUG_KDGRAPH2
        if ((testdat.testval != STRATTESTVAL) ||
            (testdat.nodeval != STRATPARAMLOG)) {
          errorPrint ("kdgraphMapSt: invalid test result");
          o = 1;
          break;
        }
#endif /* SCOTCH_DEBUG_KDGRAPH2 */
        if (testdat.data.val.vallog == 1)         /* If expression is true */
          o = kdgraphMapSt (grafptr, mappptr, straptr->data.conddat.stratab[0]); /* Apply first strategy    */
        else {                                    /* Else if expression is false                            */
          if (straptr->data.conddat.stratab[1] != NULL)  /* And if there is an else statement               */
            o = kdgraphMapSt (grafptr, mappptr, straptr->data.conddat.stratab[1]); /* Apply second strategy */
        }
      }
      break;
    case STRATNODEEMPTY :
      break;
    case STRATNODESELECT :
      errorPrint ("kdgraphMapSt: selection operator not implemented for k-way strategies");
      return (1);
#ifdef SCOTCH_DEBUG_KDGRAPH1
    case STRATNODEMETHOD :
#else  /* SCOTCH_DEBUG_KDGRAPH1 */
    default :
#endif /* SCOTCH_DEBUG_KDGRAPH1 */
      return (((KdgraphMapFunc) straptr->tablptr->methtab[straptr->data.methdat.methnum].funcptr)
              (grafptr, mappptr, (const void * const) &straptr->data.methdat.datadat));
#ifdef SCOTCH_DEBUG_KDGRAPH1
    default :
      errorPrint ("kdgraphMapSt: invalid parameter (2)");
      return (1);
#endif /* SCOTCH_DEBUG_KDGRAPH1 */
  }
  return (o);
}

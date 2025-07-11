/* Copyright 2007,2008,2023,2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : hdgraph_order_st.c                      **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This module is the generic call to the  **/
/**                distributed graph ordering module,      **/
/**                using a given strategy.                 **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 15 apr 2006     **/
/**                                 to   : 21 aug 2006     **/
/**                # Version 5.1  : from : 11 nov 2008     **/
/**                                 to   : 11 nov 2008     **/
/**                # Version 7.0  : from : 19 jan 2023     **/
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
#include "order.h"
#include "hgraph.h"
#include "hgraph_order_st.h"
#include "dgraph.h"
#include "dorder.h"
#include "hdgraph.h"
#include "hdgraph_order_nd.h"
#include "hdgraph_order_si.h"
#include "hdgraph_order_sq.h"
#include "hdgraph_order_st.h"
#include "vdgraph.h"
#include "vdgraph_separate_st.h"

/*
**  The static and global variables.
*/

static Hdgraph              hdgraphorderstgraphdummy; /* Dummy graph for offset computations */

static union {                                    /* Default parameters for nested dissection method */
  HdgraphOrderNdParam       param;
  StratNodeMethodData       padding;
} hdgraphorderstdefaultnd = { { &stratdummy, &stratdummy, &stratdummy } };

static union {                                    /* Default parameters for sequential method */
  HdgraphOrderSqParam       param;
  StratNodeMethodData       padding;
} hdgraphorderstdefaultsq = { { &stratdummy } };

static StratMethodTab       hdgraphorderstmethtab[] = { /* Graph ordering methods array */
                              { HDGRAPHORDERSTMETHND, "n",  (StratMethodFunc) hdgraphOrderNd, &hdgraphorderstdefaultnd },
                              { HDGRAPHORDERSTMETHSI, "s",  (StratMethodFunc) hdgraphOrderSi, NULL },
                              { HDGRAPHORDERSTMETHSQ, "q",  (StratMethodFunc) hdgraphOrderSq, &hdgraphorderstdefaultsq },
                              { -1,                   NULL, (StratMethodFunc) NULL,           NULL } };

static StratParamTab        hdgraphorderstparatab[] = { /* The method parameter list */
                              { HDGRAPHORDERSTMETHND,  STRATPARAMSTRAT,  "sep",
                                (byte *) &hdgraphorderstdefaultnd.param,
                                (byte *) &hdgraphorderstdefaultnd.param.sepstrat,
                                (void *) &vdgraphseparateststratab },
                              { HDGRAPHORDERSTMETHND,  STRATPARAMSTRAT,  "ole",
                                (byte *) &hdgraphorderstdefaultnd.param,
                                (byte *) &hdgraphorderstdefaultnd.param.ordstratlea,
                                (void *) &hdgraphorderststratab },
                              { HDGRAPHORDERSTMETHND,  STRATPARAMSTRAT,  "ose",
                                (byte *) &hdgraphorderstdefaultnd.param,
                                (byte *) &hdgraphorderstdefaultnd.param.ordstratsep,
                                (void *) &hdgraphorderststratab },
                              { HDGRAPHORDERSTMETHND,  STRATPARAMSTRAT,  "osq",
                                (byte *) &hdgraphorderstdefaultnd.param,
                                (byte *) &hdgraphorderstdefaultnd.param.ordstratseq,
                                (void *) &hgraphorderststratab },
                              { HDGRAPHORDERSTMETHSQ,  STRATPARAMSTRAT,  "strat",
                                (byte *) &hdgraphorderstdefaultsq.param,
                                (byte *) &hdgraphorderstdefaultsq.param.ordstratseq,
                                (void *) &hgraphorderststratab },
                              { HDGRAPHORDERSTMETHNBR, STRATPARAMINT,    NULL,
                                NULL, NULL, NULL } };

static StratParamTab        hdgraphorderstcondtab[] = { /* Graph condition parameter table */
                              { STRATNODECOND,         STRATPARAMINT,    "edge",
                                (byte *) &hdgraphorderstgraphdummy,
                                (byte *) &hdgraphorderstgraphdummy.s.edgeglbnbr,
                                NULL },
                              { STRATNODECOND,         STRATPARAMINT,    "levl",
                                (byte *) &hdgraphorderstgraphdummy,
                                (byte *) &hdgraphorderstgraphdummy.levlnum,
                                NULL },
                              { STRATNODECOND,         STRATPARAMINT,    "load",
                                (byte *) &hdgraphorderstgraphdummy,
                                (byte *) &hdgraphorderstgraphdummy.s.veloglbsum,
                                NULL },
                              { STRATNODECOND,         STRATPARAMDOUBLE, "mdeg",
                                (byte *) &hdgraphorderstgraphdummy,
                                (byte *) &hdgraphorderstgraphdummy.s.degrglbmax,
                                NULL },
                              { STRATNODECOND,         STRATPARAMINT,    "proc",
                                (byte *) &hdgraphorderstgraphdummy,
                                (byte *) &hdgraphorderstgraphdummy.s.procglbnbr,
                                NULL },
                              { STRATNODECOND,         STRATPARAMINT,    "rank",
                                (byte *) &hdgraphorderstgraphdummy,
                                (byte *) &hdgraphorderstgraphdummy.s.proclocnum,
                                NULL },
                              { STRATNODECOND,         STRATPARAMINT,    "vert",
                                (byte *) &hdgraphorderstgraphdummy,
                                (byte *) &hdgraphorderstgraphdummy.s.vertglbnbr,
                                NULL },
                              { STRATNODENBR,          STRATPARAMINT,    NULL,
                                NULL, NULL, NULL } };

StratTab                    hdgraphorderststratab = { /* Strategy tables for graph ordering methods */
                              hdgraphorderstmethtab,
                              hdgraphorderstparatab,
                              hdgraphorderstcondtab };

/************************************/
/*                                  */
/* This routine is the entry point  */
/* for the graph ordering routines. */
/*                                  */
/************************************/

/* This routine computes an ordering
** with respect to a given strategy.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

int
hdgraphOrderSt (
Hdgraph * restrict const      grafptr,            /*+ Subgraph to order       +*/
DorderCblk * restrict const   cblkptr,            /*+ Current column block    +*/
const Strat * restrict const  straptr)            /*+ Graph ordering strategy +*/
{
  StratTest           testdat;
  int                 o;

  if (grafptr->s.vertglbnbr == 0)                 /* Return immediately if nothing to do */
    return (0);

  o = 0;
  switch (straptr->typeval) {
    case STRATNODECONCAT :
      errorPrint ("hdgraphOrderSt: concatenation operator not available for graph ordering strategies");
      return (1);
    case STRATNODECOND :
      o = stratTestEval (straptr->data.conddat.testptr, &testdat, (void *) grafptr); /* Evaluate expression */
      if (o == 0) {                               /* If evaluation was correct */
#ifdef SCOTCH_DEBUG_HDGRAPH2
        if ((testdat.testval != STRATTESTVAL) &&
            (testdat.nodeval != STRATPARAMLOG)) {
          errorPrint ("hdgraphOrderSt: invalid test result");
          o = 1;
          break;
        }
#endif /* SCOTCH_DEBUG_HDGRAPH2 */
        if (testdat.data.val.vallog == 1)         /* If expression is true */
          o = hdgraphOrderSt (grafptr, cblkptr, straptr->data.conddat.stratab[0]); /* Apply first strategy    */
        else {                                    /* Else if expression is false                              */
          if (straptr->data.conddat.stratab[1] != NULL)  /* And if there is an else statement                 */
            o = hdgraphOrderSt (grafptr, cblkptr, straptr->data.conddat.stratab[1]); /* Apply second strategy */
        }
      }
      break;
    case STRATNODEEMPTY :
      hdgraphOrderSi (grafptr, cblkptr);          /* Always maintain a consistent ordering */
      break;
    case STRATNODESELECT :
      errorPrint ("hdgraphOrderSt: selection operator not available for graph ordering strategies");
      return (1);
#ifdef SCOTCH_DEBUG_HDGRAPH2
    case STRATNODEMETHOD :
#else /* SCOTCH_DEBUG_HDGRAPH2 */
    default :
#endif /* SCOTCH_DEBUG_HDGRAPH2 */
      return (((HdgraphOrderFunc) (straptr->tablptr->methtab[straptr->data.methdat.methnum].funcptr))
              (grafptr, cblkptr, (const void * const) &straptr->data.methdat.datadat));
#ifdef SCOTCH_DEBUG_HDGRAPH2
    default :
      errorPrint ("hdgraphOrderSt: invalid parameter");
      return (1);
#endif /* SCOTCH_DEBUG_HDGRAPH2 */
  }
  return (o);
}

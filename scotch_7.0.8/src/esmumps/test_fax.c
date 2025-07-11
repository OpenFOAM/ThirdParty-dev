/* Copyright 2022,2024 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : test_fax.c                              **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : This is the test module for the         **/
/**                symbolic factorization routine.         **/
/**                                                        **/
/**   DATES      : # Version 7.0  : from : 21 apr 2022     **/
/**                                 to   : 10 sep 2024     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "module.h"
#include "common.h"
#include "scotch.h"
#include "graph.h"
#include "order.h"
#include "dof.h"
#include "symbol.h"
#include "fax.h"

/******************************/
/*                            */
/* This is the main function. */
/*                            */
/******************************/

int
main (
int                         argc,
char *                      argv[])
{
  Graph               grafdat;                    /* Graph to load */
  Order               ordedat;
  SymbolMatrix        symbdat;
  Dof                 deofdat;
  FILE *              stream;
  double              nonzval;
  double              opcoval;

  if (argc != 3) {
    errorPrint ("test_fax: usage: test_fax graph_file ordering_file");
    exit       (EXIT_FAILURE);
  }

  graphInit (&grafdat);
  if ((stream = fopen (argv[1], "r")) == NULL) {
    errorPrint ("test_fax: cannot open graph file");
    graphExit  (&grafdat);
    exit       (EXIT_FAILURE);
  }
  if (graphLoad (&grafdat, stream, -1, 3) != 0) { /* Graph with untouched base value and loads */
    errorPrint ("test_fax: cannot load graph file");
    graphExit  (&grafdat);
    exit       (EXIT_FAILURE);
  }
  fclose (stream);

  orderInit (&ordedat);
  if ((stream = fopen (argv[2], "r")) == NULL) {
    errorPrint ("test_fax: cannot open ordering file");
    orderExit  (&ordedat);
    graphExit  (&grafdat);
    exit       (EXIT_FAILURE);
  }
  if (orderLoad (&ordedat, stream) != 0) {
    errorPrint ("test_fax: cannot load ordering file");
    orderExit  (&ordedat);
    graphExit  (&grafdat);
    exit       (EXIT_FAILURE);
  }
  fclose (stream);

  symbolInit (&symbdat);
  if (symbolFaxGraph (&symbdat, &grafdat, &ordedat) != 0) {
    errorPrint ("test_fax: error in symbolic factorization");
    exit       (EXIT_FAILURE);
  }

  dofInit  (&deofdat);
  dofGraph (&deofdat, &grafdat, 1, ordedat.peritab);

  if (symbolCost (&symbdat, &deofdat, SYMBOLCOSTLDLT, &nonzval, &opcoval) != 0) {
    errorPrint ("test_fax: error in symbolic cost computation");
    exit       (EXIT_FAILURE);
  }

  printf ("NNZ: %le\nOPC: %le\n", nonzval, opcoval);

  dofExit    (&deofdat);
  symbolExit (&symbdat);
  orderExit  (&ordedat);
  graphExit  (&grafdat);

  exit (EXIT_SUCCESS);
}

!* Copyright 2021, IPB, Universite de Bordeaux, INRIA & CNRS
!*
!* This file is part of the Scotch software package for static mapping,
!* graph partitioning and sparse matrix ordering.
!*
!* This software is governed by the CeCILL-C license under French law
!* and abiding by the rules of distribution of free software. You can
!* use, modify and/or redistribute the software under the terms of the
!* CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
!* URL: "http://www.cecill.info".
!*
!* As a counterpart to the access to the source code and rights to copy,
!* modify and redistribute granted by the license, users are provided
!* only with a limited warranty and the software's author, the holder of
!* the economic rights, and the successive licensors have only limited
!* liability.
!*
!* In this respect, the user's attention is drawn to the risks associated
!* with loading, using, modifying and/or developing or reproducing the
!* software by the user in light of its specific status of free software,
!* that may mean that it is complicated to manipulate, and that also
!* therefore means that it is reserved for developers and experienced
!* professionals having in-depth computer knowledge. Users are therefore
!* encouraged to load and test the software's suitability as regards
!* their requirements in conditions enabling the security of their
!* systems and/or data to be ensured and, more generally, to use and
!* operate it in the same conditions as regards security.
!*
!* The fact that you are presently reading this means that you have had
!* knowledge of the CeCILL-C license and that you accept its terms.
!
!***********************************************************
!*                                                        **
!*   NAME       : library_metis_f.h                       **
!*                                                        **
!*   AUTHOR     : Marc FUENTES                            **
!*                                                        **
!*   FUNCTION   : FORTRAN declaration file for the        **
!*                libmetis options vector.                **
!*                                                        **
!*   DATES      : # Version 6.1  : from : 18 may 2021     **
!*                                 to   : 27 may 2021     **
!*                                                        **
!***********************************************************

!* Size of option vector

        INTEGER METIS_NOPTIONS
        PARAMETER (METIS_NOPTIONS = 40)

!* Definition of enum OPTIONS

        INTEGER METIS_OPTION_PTYPE
        INTEGER METIS_OPTION_OBJTYPE
        INTEGER METIS_OPTION_CTYPE
        INTEGER METIS_OPTION_IPTYPE
        INTEGER METIS_OPTION_RTYPE
        INTEGER METIS_OPTION_DBGLVL
        INTEGER METIS_OPTION_NITER
        INTEGER METIS_OPTION_NCUTS
        INTEGER METIS_OPTION_SEED
        INTEGER METIS_OPTION_NO2HOP
        INTEGER METIS_OPTION_MINCONN
        INTEGER METIS_OPTION_CONTIG
        INTEGER METIS_OPTION_COMPRESS
        INTEGER METIS_OPTION_CCORDER
        INTEGER METIS_OPTION_PFACTOR
        INTEGER METIS_OPTION_NSEPS
        INTEGER METIS_OPTION_UFACTOR
        INTEGER METIS_OPTION_NUMBERING
        PARAMETER (METIS_OPTION_PTYPE =      1)
        PARAMETER (METIS_OPTION_OBJTYPE =    2)
        PARAMETER (METIS_OPTION_CTYPE =      3)
        PARAMETER (METIS_OPTION_IPTYPE =     4)
        PARAMETER (METIS_OPTION_RTYPE =      5)
        PARAMETER (METIS_OPTION_DBGLVL =     6)
        PARAMETER (METIS_OPTION_NITER =      7)
        PARAMETER (METIS_OPTION_NCUTS =      8)
        PARAMETER (METIS_OPTION_SEED =       9)
        PARAMETER (METIS_OPTION_NO2HOP =    10)
        PARAMETER (METIS_OPTION_MINCONN =   11)
        PARAMETER (METIS_OPTION_CONTIG =    12)
        PARAMETER (METIS_OPTION_COMPRESS =  13)
        PARAMETER (METIS_OPTION_CCORDER =   14)
        PARAMETER (METIS_OPTION_PFACTOR =   15)
        PARAMETER (METIS_OPTION_NSEPS =     16)
        PARAMETER (METIS_OPTION_UFACTOR =   17)
        PARAMETER (METIS_OPTION_NUMBERING = 18)

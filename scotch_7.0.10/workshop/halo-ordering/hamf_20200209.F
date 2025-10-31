C Version generated on October 29th 2019
C
C This file includes various modifications of an original
C LGPL/ CeCILL-C compatible
C code implementing the AMD (Approximate Minimum Degree) ordering
C   Patrick Amestoy, Timothy A. Davis, and Iain S. Duff,
C    "An approximate minimum degree ordering algorithm,"
C    SIAM J. Matrix Analysis  vol 17, pages=886--905 (1996)
C    MUMPS_ANA_H is based on the original AMD code:
C
C    AMD, Copyright (c), 1996-2016, Timothy A. Davis,
C    Patrick R. Amestoy, and Iain S. Duff.  All Rights Reserved.
C    Used in MUMPS under the BSD 3-clause license.
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
C    CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
C    BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
C    FITNESS FOR A PARTICULAR PURPOSE
C    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS BE LIABLE FOR
C    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
C    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
C    OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
C    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
C    STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
C    IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
C    POSSIBILITY OF SUCH DAMAGE.
C
C   MUMPS_HAMF4 is a major modification of AMD
C   since the metric used to select pivots in not anymore the
C   degree but an approximation of the fill-in.
C   In this approximation
C   all cliques of elements adjacent to the variable are deducted.
C   MUMPS_HAMF4 also enables to take into account a halo in the graph.
C   The graph is composed is partitioned in two types of nodes
C   the so called internal nodes and the so called halo nodes.
C   Halo nodes cannot be selected the both the initial degrees
C   and updated degrees of internal node should be taken
C   into account.
C   Written by Patrick Amestoy between 1999 and 2019.
C   and used by F. Pellegrini in SCOTCH since 2000.
C
C   Unique version to order both graph of variables and
C   graphs with both elements and variables.
C
C   Notations used:
C   Let us refer to as
C     Gv a graph with only variables
C     Ge a graph with both variables and elements
C
C   Let V be the set of nodes in the graph
C       V = Ve + V0 + V1
C           V0 = Set of variable nodes (not in halo)
C           V1 = Set of variable nodes (in halo)
C           Ve = Set of element nodes |Ve|=nbelts
C
C   All 3 sets are disjoint, Ve and V1 can be empty
C   If nbelts>0 then a bipartite graph bewteen
C      (V0 U V1) and Ve is provided on entry.

C   A graph of elements and variables is a bipartite
C   graph between the set of variables (that may
C   include halo variables) and the set of elements.
C   Thus variables are only adjacent to elements and
C   in the element list we only have variables.
C   Elements are "considered" as already eliminated and
C   are provided here only to describe the adjacency between
C   variables. Only variables in V0 need to be eliminated.
C
C Comments relative to 64/32 bits integer representation:
C Restrictive integer 64 bit variant :
C     it is assumed that IW array size can exceed 32-bit integer
C     and thus iwlen is INTEGER(8) and pe in an INTEGER(8) array
C     Graphe size n must be smaller than 2x10^9 but number of
C     edges is a 64-bit integer.
C
C
      SUBROUTINE MUMPS_HAMF4(NORIG, N, NBELTS, NBBUCK,
     &                   IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
     &                   LAST, NCMPA, DEGREE, WF, NEXT, W, HEAD
     &                   )
      IMPLICIT NONE
C
C Parameters
C    Input not modified
      INTEGER, INTENT(IN)    :: NORIG, N, NBELTS, NBBUCK
      INTEGER(8), INTENT(IN) :: IWLEN
C     Input undefined on output
      INTEGER, INTENT(INOUT)  :: LEN(N), IW(IWLEN)
C     NV also meaningful as input to encode compressed graphs
      INTEGER, INTENT(INOUT)  :: NV(N)
C
C     Output only
      INTEGER, INTENT(OUT)   :: NCMPA
      INTEGER, INTENT(OUT)   :: ELEN(N), LAST(N)
C
C     Input/output
      INTEGER(8), INTENT(INOUT) :: PFREE
      INTEGER(8), INTENT(INOUT) :: PE(N)
C
C     Internal Workspace only
C       Min fill approximation one extra array of size NBBUCK+2
C       is also needed
      INTEGER     :: NEXT(N), DEGREE(N), W(N)
      INTEGER     :: HEAD(0:NBBUCK+1), WF(N)
C
C  Comments on the OUTPUT:
C  ----------------------
C  Let V= V0 U V1 the nodes of the initial graph (|V|=n).
C  The assembly tree corresponds to the tree
C    of the supernodes (or supervariables). Each node of the
C    assembly tree is then composed of one principal variable
C    and a list of secondary variables. The list of
C    variable of a node (principal + secondary variables) then
C    describes the structure of the diagonal bloc of the
C    supernode.
C  The elimination tree denotes the tree of all the variables(=node) and
C    is therefore of order n.
C
C  The arrays NV(N) and PE(N) give a description of the
C  assembly tree.
C
C   1/ Description of array nv(N) (on OUTPUT)
C    nv(i)=0 i is a secondary variable
C    nv(i) >0 i is a principal variable, nv(i) holds the
C                  the number of elements in column i of L (true degree of i)
C    with compressed graph (nv(1).ne.-1 on input),
C    nv(i) can be greater than N since degree can be as large as NORIG.
C
C   2/ Description of array PE(N) (on OUTPUT)
C       Note that on
C       pe(i) = -(father of variable/node i) in the elimination tree:
C       If nv (i) .gt. 0, then i represents a node in the assembly tree,
C       and the parent of i is -pe (i), or zero if i is a root.
C       If nv (i) = 0, then (i,-pe (i)) represents an edge in a
C       subtree, the root of which is a node in the assembly tree.
C
C   3/ Example:
C      Let If be a root node father of Is in the assembly tree.
C      If is the principal
C      variable of the node If and let If1, If2, If3 be the secondary variables
C      of node If.
C      Is is the principal
C      variable of the node Is and let Is1, Is2 be the secondary variables
C      of node Is.
C
C      THEN:
C        NV(If1)=NV(If2)=NV(If3) = 0  (secondary variables)
C        NV(Is1)=NV(Is2) = 0  (secondary variables)
C        NV(If) > 0  ( principal variable)
C        NV(Is) > 0  ( principal variable)
C        PE(If)  = 0 (root node)
C        PE(Is)  = -If (If is the father of Is in the assembly tree)
C        PE(If1)=PE(If2)=PE(If3)= -If  ( If is the principal variable)
C        PE(Is1)=PE(Is2)= -Is  ( Is is the principal variable)
C
C
C
C HALOAMD_V1: (September 1997)
C **********
C Initial version designed to experiment the numerical (fill-in) impact
C of taking into account the halo. This code should be able
C to experiment no-halo, partial halo, complete halo.
C DATE: September 17th 1997
C
C HALOAMD is designed to process a gragh composed of two types
C            of nodes, V0 and V1, extracted from a larger gragh.
C            V0^V1 = {},
C
C            We used Min. degree heuristic to order only
C            nodes in V0, but the adjacency to nodes
C            in V1 is taken into account during ordering.
C            Nodes in V1 are odered at last.
C            Adjacency between nodes of V1 need not be provided,
C            however |len(i)| must always corresponds to the number of
C            edges effectively provided in the adjacency list of i.
C          On input :
c          ********
C            Nodes INODE in V1 are flagged with len(INODE) = -degree
C                           if len(i) =0 and i \in V1 then
C                           len(i) must be set on input to -NORIG-1
C          ERROR return (negative values in ncmpa)
C          ************
C            negative value in ncmpa indicates an error detected
C               by HALOAMD.
C
C            The graph provided MUST follow the rule:
C             if (i,j) is an edge in the gragh then
C             j must be in the adjacency list of i AND
C             i must be in the adjacency list of j.
C    REMARKS
C    -------
C
C        1/  Providing edges between nodes of V1 should not
C            affect the final ordering, only the amount of edges
C            of the halo should effectively affect the solution.
C            This code should work in the following cases:
C              1/ halo not provided
C              2/ halo partially provided
C              3/ complete halo
C              4/ complete halo+interconnection between nodes of V1.
C
C              1/ should run and provide identical results (w.r.t to current
C               implementation of AMD in SCOTCH).
C             3/ and 4 should provide identical results.
C
C        2/ All modifications of the AMD initial code are indicated
C           with begin HALO .. end HALO
C
C
C   Given a representation of the nonzero pattern of a symmetric matrix,
C       A, (excluding the diagonal) perform an approximate minimum
C       fill-in heuristic. Aggresive absorption is
C       used to tighten the bound on the degree.  This can result an
C       significant improvement in the quality of the ordering for
C       some matrices.
C-----------------------------------------------------------------------
C INPUT ARGUMENTS (unaltered):
C-----------------------------------------------------------------------
C n:    number of nodes in the complete graph including halo
C       n = size of V0 U V1 U Ve
C       Restriction:  n .ge. 1
C
C norig   if compressed graph (nv(1).ne.-1) then
C            norig is the sum(nv(i)) for i \in [1:n]
C            and could ne larger than n
C         else norig = n
C
C nbbuck should be greater than norig
C        advised value is 2*norig
C
C nbelts number of elements (size of Ve)
C            =0 if Gv (graph of variables)
C            >0 if Ge
C         nbelts > 0 extends/changes the meaning of
C                   len/elen on entry (see below)
C
C
C iwlen:        The length of iw (1..iwlen).  On input, the matrix is
C       stored in iw (1..pfree-1).  However, iw (1..iwlen) should be
C       slightly larger than what is required to hold the matrix, at
C       least iwlen .ge. pfree + n is recommended.  Otherwise,
C       excessive compressions will take place.
C       *** We do not recommend running this algorithm with ***
C       ***      iwlen .lt. pfree + n.                      ***
C       *** Better performance will be obtained if          ***
C       ***      iwlen .ge. pfree + n                       ***
C       *** or better yet                                   ***
C       ***      iwlen .gt. 1.2 * pfree                     ***
C       *** (where pfree is its value on input).            ***
C       The algorithm will not run at all if iwlen .lt. pfree-1.
C
C       Restriction: iwlen .ge. pfree-1
C-----------------------------------------------------------------------
C INPUT/OUPUT ARGUMENTS:
C-----------------------------------------------------------------------
C pe:   On input, pe (i) is the index in iw of the start of row i, or
C       zero if row i has no off-diagonal non-zeros.
C
C       During execution, it is used for both supervariables and
C       elements:
C
C       * Principal supervariable i:  index into iw of the
C               description of supervariable i.  A supervariable
C               represents one or more rows of the matrix
C               with identical nonzero pattern.
C       * Non-principal supervariable i:  if i has been absorbed
C               into another supervariable j, then pe (i) = -j.
C               That is, j has the same pattern as i.
C               Note that j might later be absorbed into another
C               supervariable j2, in which case pe (i) is still -j,
C               and pe (j) = -j2.
C       * Unabsorbed element e:  the index into iw of the description
C               of element e, if e has not yet been absorbed by a
C               subsequent element.  Element e is created when
C               the supervariable of the same name is selected as
C               the pivot.
C       * Absorbed element e:  if element e is absorbed into element
C               e2, then pe (e) = -e2.  This occurs when the pattern of
C               e (that is, Le) is found to be a subset of the pattern
C               of e2 (that is, Le2).  If element e is "null" (it has
C               no nonzeros outside its pivot block), then pe (e) = 0.
C
C       On output, pe holds the assembly tree/forest, which implicitly
C       represents a pivot order with identical fill-in as the actual
C       order (via a depth-first search of the tree).
C
C       On output:
C       If nv (i) .gt. 0, then i represents a node in the assembly tree,
C       and the parent of i is -pe (i), or zero if i is a root.
C       If nv (i) = 0, then (i,-pe (i)) represents an edge in a
C       subtree, the root of which is a node in the assembly tree.
C
C pfree:        On input, the matrix is stored in iw (1..pfree-1) and
C       the rest of the array iw is free.
C       During execution, additional data is placed in iw, and pfree
C       is modified so that components  of iw from pfree are free.
C       On output, pfree is set equal to the size of iw that
C       would have been needed for no compressions to occur.  If
C       ncmpa is zero, then pfree (on output) is less than or equal to
C       iwlen, and the space iw (pfree+1 ... iwlen) was not used.
C       Otherwise, pfree (on output) is greater than iwlen, and all the
C       memory in iw was used.
C
C nv:   On input, encoding of compressed graph:
C        if nv(1) = -1 then graph is not compressed otherwise
C        nv(I) holds the weight of node i.
C       During execution, abs (nv (i)) is equal to the number of rows
C       that are represented by the principal supervariable i.  If i is
C       a nonprincipal variable, then nv (i) = 0.  Initially,
C       nv (i) = 1 for all i.  nv (i) .lt. 0 signifies that i is a
C       principal variable in the pattern Lme of the current pivot
C       element me.  On output, nv (e) holds the true degree of element
C       e at the time it was created (including the diagonal part).
C begin HALO
C       On output, nv(I) can be used to find node in set V1.
C       Not true anymore : ( nv(I) = N+1 characterizes nodes in V1.
C                 instead nodes in V1 are considered as a dense root node )
C end HALO
C-----------------------------------------------------------------------
C INPUT/MODIFIED (undefined on output):
C-----------------------------------------------------------------------
C len:  On input, len (i)
C           positive or null (>=0) : i \in V0 U Ve and
C                   if (nbelts==0) then (graph of variables)
C                     len(i) holds the number of entries in row i of the
C                     matrix, excluding the diagonal.
C                   else (graph of elements+variables)
C                     if i\in V0 then len(i) = nb of elements adjacent to i
C                     if i\in Ve then len(i) = nb of variables adjacent to i
C                   endif
C           negative (<0) : i \in V1, and
C                   if (nbelts==0) then (graph of variables)
C                     -len(i) hold the number of entries in row i of the
C                      matrix, excluding the diagonal.
C                      len(i) = - | Adj(i) | if i \in V1
C                              or -N -1 if  | Adj(i) | = 0 and i \in V1
C                   else  (graph of elements+variables)
C                     -len(i) nb of elements adjacent to i
C                   endif
C       The content of len (1..n) is undefined on output.
C
C elen: defined on input only if nbelts>0
C       if e \in Ve then elen (e) = -N-1
C       if v \in V0 then elen (v) = External degree of v
C                                   that should thus be provided
C                                   on entry to initialize degree
C
C       if v \in V1 then elen (v) = 0
C
C iw:   On input, iw (1..pfree-1) holds the description of each row i
C       in the matrix.  The matrix must be symmetric, and both upper
C       and lower triangular parts must be present.  The diagonal must
C       not be present.  Row i is held as follows:
C
C               len (i):  the length of the row i data structure
C               iw (pe (i) ... pe (i) + len (i) - 1):
C                       the list of column indices for nonzeros
C                       in row i (simple supervariables), excluding
C                       the diagonal.  All supervariables start with
C                       one row/column each (supervariable i is just
C                       row i).
C               if len (i) is zero on input, then pe (i) is ignored
C               on input.
C
C               Note that the rows need not be in any particular order,
C               and there may be empty space between the rows.
C
C       During execution, the supervariable i experiences fill-in.
C       This is represented by placing in i a list of the elements
C       that cause fill-in in supervariable i:
C
C               len (i):  the length of supervariable i
C               iw (pe (i) ... pe (i) + elen (i) - 1):
C                       the list of elements that contain i.  This list
C                       is kept short by removing absorbed elements.
C               iw (pe (i) + elen (i) ... pe (i) + len (i) - 1):
C                       the list of supervariables in i.  This list
C                       is kept short by removing nonprincipal
C                       variables, and any entry j that is also
C                       contained in at least one of the elements
C                       (j in Le) in the list for i (e in row i).
C
C       When supervariable i is selected as pivot, we create an
C       element e of the same name (e=i):
C
C               len (e):  the length of element e
C               iw (pe (e) ... pe (e) + len (e) - 1):
C                       the list of supervariables in element e.
C
C       An element represents the fill-in that occurs when supervariable
C       i is selected as pivot (which represents the selection of row i
C       and all non-principal variables whose principal variable is i).
C       We use the term Le to denote the set of all supervariables
C       in element e.  Absorbed supervariables and elements are pruned
C       from these lists when computationally convenient.
C
C       CAUTION:  THE INPUT MATRIX IS OVERWRITTEN DURING COMPUTATION.
C       The contents of iw are undefined on output.
C
C-----------------------------------------------------------------------
C OUTPUT (need not be set on input):
C-----------------------------------------------------------------------
C elen: See the description of iw above.  At the start of execution,
C       elen (i) is set to zero.  During execution, elen (i) is the
C       number of elements in the list for supervariable i.  When e
C       becomes an element, elen (e) = -nel is set, where nel is the
C       current step of factorization.  elen (i) = 0 is done when i
C       becomes nonprincipal.
C
C       For variables, elen (i) .ge. 0 holds
C       For elements, elen (e) .lt. 0 holds.
C
C last: In a degree list, last (i) is the supervariable preceding i,
C       or zero if i is the head of the list.  In a hash bucket,
C       last (i) is the hash key for i.  last (head (hash)) is also
C       used as the head of a hash bucket if head (hash) contains a
C       degree list (see head, below).
C
C ncmpa:        The number of times iw was compressed.  If this is
C       excessive, then the execution took longer than what could have
C       been.  To reduce ncmpa, try increasing iwlen to be 10% or 20%
C       larger than the value of pfree on input (or at least
C       iwlen .ge. pfree + n).  The fastest performance will be
C       obtained when ncmpa is returned as zero.  If iwlen is set to
C       the value returned by pfree on *output*, then no compressions
C       will occur.
C begin HALO
C        on output ncmpa <0 --> error detected during HALO_AMD:
C           error 1: ncmpa = -NORIG , ordering was stopped.
C end HALO
C
C-----------------------------------------------------------------------
C LOCAL (not input or output - used only during execution):
C-----------------------------------------------------------------------
C degree:       If i is a supervariable, then degree (i) holds the
C       current approximation of the external degree of row i (an upper
C       bound).  The external degree is the number of nonzeros in row i,
C       minus abs (nv (i)) (the diagonal part).  The bound is equal to
C       the external degree if elen (i) is less than or equal to two.
C       We also use the term "external degree" for elements e to refer
C       to |Le \ Lme|.  If e is an element, then degree (e) holds |Le|,
C       which is the degree of the off-diagonal part of the element e
C       (not including the diagonal part).
C begin HALO
C       degree(I) = n+1 indicates that i belongs to V1
C end HALO
C
C head: head is used for degree lists.  head (deg) is the first
C       supervariable in a degree list (all supervariables i in a
C       degree list deg have the same approximate degree, namely,
C       deg = degree (i)).  If the list deg is empty then
C       head (deg) = 0.
C
C       During supervariable detection head (hash) also serves as a
C       pointer to a hash bucket.
C       If head (hash) .gt. 0, there is a degree list of degree hash.
C               The hash bucket head pointer is last (head (hash)).
C       If head (hash) = 0, then the degree list and hash bucket are
C               both empty.
C       If head (hash) .lt. 0, then the degree list is empty, and
C               -head (hash) is the head of the hash bucket.
C       After supervariable detection is complete, all hash buckets
C       are empty, and the (last (head (hash)) = 0) condition is
C       restored for the non-empty degree lists.
C next: next (i) is the supervariable following i in a link list, or
C       zero if i is the last in the list.  Used for two kinds of
C       lists:  degree lists and hash buckets (a supervariable can be
C       in only one kind of list at a time).
C w:    The flag array w determines the status of elements and
C       variables, and the external degree of elements.
C
C       for elements:
C          if w (e) = 0, then the element e is absorbed
C          if w (e) .ge. wflg, then w (e) - wflg is the size of
C               the set |Le \ Lme|, in terms of nonzeros (the
C               sum of abs (nv (i)) for each principal variable i that
C               is both in the pattern of element e and NOT in the
C               pattern of the current pivot element, me).
C          if wflg .gt. w (e) .gt. 0, then e is not absorbed and has
C               not yet been seen in the scan of the element lists in
C               the computation of |Le\Lme| in loop 150 below.
C
C       for variables:
C          during supervariable detection, if w (j) .ne. wflg then j is
C          not in the pattern of variable i
C
C       The w array is initialized by setting w (i) = 1 for all i,
C       and by setting wflg = 2.  It is reinitialized if wflg becomes
C       too large (to ensure that wflg+n does not cause integer
C       overflow).
C
C wf : integer array  used to store the already filled area of
C      the variables adajcent to current pivot.
C      wf is then used to update the score of variable i.
C
C-----------------------------------------------------------------------
C LOCAL INTEGERS:
C-----------------------------------------------------------------------
      INTEGER :: DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     &        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     &        LENJ, LN, ME, MINDEG, NEL,
     &        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X,
     &        NBFLAG, NREAL, LASTD, NELME, WF3, WF4, N2, PAS
       INTEGER KNT1_UPDATED, KNT2_UPDATED
       INTEGER(8) :: MAXMEM, MEM, NEWMEM
       INTEGER    :: MAXINT_N
       INTEGER(8) :: HASH, HMOD
       DOUBLE PRECISION RMF, RMF1
       DOUBLE PRECISION dummy
       INTEGER idummy
       LOGICAL COMPRESS
C deg:        the degree of a variable or element
C degme:      size, |Lme|, of the current element, me (= degree (me))
C dext:       external degree, |Le \ Lme|, of some element e
C dmax:       largest |Le| seen so far
C e:          an element
C elenme:     the length, elen (me), of element list of pivotal var.
C eln:        the length, elen (...), of an element list
C hash:       the computed value of the hash function
C hmod:       the hash function is computed modulo hmod = max (1,n-1)
C i:          a supervariable
C ilast:      the entry in a link list preceding i
C inext:      the entry in a link list following i
C j:          a supervariable
C jlast:      the entry in a link list preceding j
C jnext:      the entry in a link list, or path, following j
C k:          the pivot order of an element or variable
C knt1:       loop counter used during element construction
C knt2:       loop counter used during element construction
C knt3:       loop counter used during compression
C lenj:       len (j)
C ln:         length of a supervariable list
C maxint_n:   large integer to test risk of overflow on wflg
C maxmem:     amount of memory needed for no compressions
C me:         current supervariable being eliminated, and the
C                     current element created by eliminating that
C                     supervariable
C mem:        memory in use assuming no compressions have occurred
C mindeg:     current minimum degree
C nel:        number of pivots selected so far
C newmem:     amount of new memory needed for current pivot element
C nleft:      n - nel, the number of nonpivotal rows/columns remaining
C nvi:        the number of variables in a supervariable i (= nv (i))
C nvj:        the number of variables in a supervariable j (= nv (j))
C nvpiv:      number of pivots in current element
C slenme:     number of variables in variable list of pivotal variable
C we:         w (e)
C wflg:       used for flagging the w array.  See description of iw.
C wnvi:       wflg - nv (i)
C x:          either a supervariable or an element
C wf3:  off diagoanl block area
C wf4:  diagonal block area
C mf : Minimum fill
C begin HALO
C nbflag:     number of flagged entries in the initial gragh.
C nreal :     number of entries on which ordering must be perfomed
C             (nreel = N- nbflag)
C nelme number of pivots selected when reaching the root
C lastd index of the last row in the list of dense rows
C end HALO
C-----------------------------------------------------------------------
C LOCAL POINTERS:
C-----------------------------------------------------------------------
      INTEGER(8) :: P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2,
     &              PN, PSRC
C             Any parameter (pe (...) or pfree) or local variable
C             starting with "p" (for Pointer) is an index into iw,
C             and all indices into iw use variables starting with
C             "p."  The only exception to this rule is the iwlen
C             input argument.
C p:          pointer into lots of things
C p1:         pe (i) for some variable i (start of element list)
C p2:         pe (i) + elen (i) -  1 for some var. i (end of el. list)
C p3:         index of first supervariable in clean list
C pdst:       destination pointer, for compression
C pend:       end of memory to compress
C pj:         pointer into an element or variable
C pme:        pointer into the current element (pme1...pme2)
C pme1:       the current element, me, is stored in iw (pme1...pme2)
C pme2:       the end of the current element
C pn:         pointer into a "clean" variable, also used to compress
C psrc:       source pointer, for compression
C-----------------------------------------------------------------------
C  FUNCTIONS CALLED:
C-----------------------------------------------------------------------
      INTRINSIC max, min, mod, huge
      INTEGER TOTEL
C=======================================================================
C  INITIALIZATIONS
C=======================================================================
C     HEAD (0:NBBUCK+1)
C
C idummy holds the largest integer - 1
C dummy  = dble (idummy)
      idummy = huge(idummy) - 1
      dummy = dble(idummy)
C     variable with degree equal to N2 are in halo
C     bucket NBBUCK+1 used for HALO variables
      N2 = -NBBUCK-1
C Distance betweeen elements of the N, ..., NBBUCK entries of HEAD
C
      PAS = max((NORIG/8), 1)
      WFLG = 2
      MAXINT_N=huge(WFLG)-NORIG
      NCMPA = 0
      NEL = 0
      HMOD = int(max (1, NBBUCK-1),kind=8)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      MINDEG = 0
C
      NBFLAG = 0
      LASTD  = 0
      HEAD(0:NBBUCK+1) = 0
      DO 10 I = 1, N
        LAST(I) = 0
        W(I) = 1
   10 CONTINUE
      IF(NV(1) .LT. 0) THEN
C        uncompress graph
C        NV should be set to 1
         DO I=1, N
          NV(I) = 1
         ENDDO
         COMPRESS = .FALSE.
      ELSE
         COMPRESS = .TRUE.
      ENDIF
      TOTEL = 0
      IF (NBELTS.EQ.0) THEN
C        =========================
C        Graph with only variables
C        =========================
         DO I=1,N
          ELEN(I) = 0
         ENDDO
         DO I=1,N
            IF (LEN(I).LT.0) THEN
               DEGREE (I) = N2
               NBFLAG     = NBFLAG +1
               IF (LEN(I).EQ.-NORIG-1) THEN
C     variable in V1 with empty adj list
                  LEN (I)    = 0
C     Because of compress, we force skipping this
C     entry which is anyway empty
                  PE (I)     = 0_8
               ELSE
                  LEN (I)    = - LEN(I)
               ENDIF
C       end HALO V3
            ELSE
               TOTEL = TOTEL + NV(I)
               IF (.NOT.COMPRESS) THEN
                 DEGREE(I) = LEN(I)
               ELSE
                 DEGREE(I) = 0
                 DO P= PE(I) , PE(I)+int(LEN(I)-1,8)
                  DEGREE(I) = DEGREE(I) + NV(IW(P))
                 ENDDO
               ENDIF
            ENDIF
         ENDDO
      ELSE
C        ===============================================
C        Bipartite graph between variables and elements
C        ===============================================
        DO I=1,N
           IF (LEN(I).LT.0) THEN
C          i \in V1
               DEGREE (I) = N2
               NBFLAG     = NBFLAG +1
               IF (LEN(I).EQ.-NORIG-1) THEN
C     variable in V1 with empty adj list
                  LEN (I)    = 0
C     Because of compress, we force skipping this
C     entry which is anyway empty
                  PE (I)     = 0_8
                  ELEN(I)    = 0
               ELSE
                  LEN (I)    = - LEN(I)
C     only elements adjacent to a variable
                  ELEN (I)   = LEN(I)
               ENDIF
           ELSE
              IF (ELEN(I).LT.0) THEN
C             i \in Ve
                NEL     = NEL + NV(I)
                ELEN(I) = -NEL
                IF (.NOT.COMPRESS) THEN
                 DEGREE(I) = LEN(I)
                ELSE
                 DEGREE(I) = 0
                 DO P= PE(I) , PE(I)+int(LEN(I)-1,8)
                  DEGREE(I) = DEGREE(I) + NV(IW(P))
                 ENDDO
                ENDIF
                DMAX = MAX(DMAX, DEGREE(I))
              ELSE
C             i \in V0
                TOTEL = TOTEL + NV(I)
                DEGREE(I) = ELEN(I)
                ELEN(I)   = LEN(I)
              ENDIF
           ENDIF
        ENDDO

      ENDIF

#if defined(SCOTCH_DEBUG_ORDER2)
      IF (NBELTS.NE.NEL) THEN
        WRITE(6,*) " Error 8Dec2003"
      ENDIF
#endif
C
C
C     number of entries to be ordered.
      NREAL = N - NBFLAG
C     ----------------------------------------------------------------
C     initialize degree lists and eliminate rows with no off-diag. nz.
C     ----------------------------------------------------------------
      DO 20 I = 1, N

C       Skip elements
        IF (ELEN(I).LT.0) CYCLE

        DEG = DEGREE (I)
        IF (DEG.EQ.N2) THEN
C           V1 variables
C           DEG = N2 (flagged variables are stored
C                  in the degree list of NBBUCK + 1
C                  (safe: because max
C                         max value of degree is NBBUCK)
C
             DEG = NBBUCK + 1
             IF (LASTD.EQ.0) THEN
C              degree list is empty
               LASTD     = I
               HEAD(DEG) = I
               NEXT(I)   = 0
               LAST(I)   = 0
             ELSE
               NEXT(LASTD) = I
               LAST(I)     = LASTD
               LASTD       = I
               NEXT(I)     = 0
             ENDIF
C
        ELSE IF (DEG .GT. 0) THEN
CN--------------
CNstrat0
CNC           DEG   = int( anint (
CNC     &             (dble(DEG)*dble(DEG-1)) / dble(2) ) )
CNC           DEG = max (DEG,1)
CNstat2
          WF(I) = DEG
C   version 1 FIXME (see Patch 04/01/04)
C          compute an estimation of fill
C            DEG = ( DEG * (DEG-1) - l*(l-1))/2
C            avec l  largest element adjacent to I
C            DEG = max (DEG, 0)
           IF (DEG.GT.NORIG) THEN
            DEG = min(((DEG-NORIG)/PAS) + NORIG , NBBUCK)
           ENDIF
C           Note that if deg=0 then
C           No fill-in will occur,
C           but one variable is adjacent to I
C          ----------------------------------------------------------
C          place i in the degree list corresponding to its degree
C          ----------------------------------------------------------
           INEXT = HEAD (DEG)
           IF (INEXT .NE. 0) LAST (INEXT) = I
           NEXT (I) = INEXT
           HEAD (DEG) = I
        ELSE
C         ----------------------------------------------------------
C         we have a variable that can be eliminated at once because
C         there is no off-diagonal non-zero in its row.
C         ----------------------------------------------------------
          NEL = NEL + NV(I)
          ELEN (I) = -NEL
          PE (I) = 0_8
          W (I) = 0
        ENDIF
C=======================================================================
C
   20 CONTINUE
C=======================================================================
C  WHILE (selecting pivots) DO
C=======================================================================
C     if elements provided (NBELTS>0) they are eliminated
      NLEFT = TOTEL-NEL
C=======================================================================
C =====================================================================
   30 IF (NEL .LT. TOTEL) THEN
C =====================================================================
C  GET PIVOT OF MINIMUM DEGREE
C=======================================================================
C       -------------------------------------------------------------
C       find next supervariable for elimination
C       -------------------------------------------------------------
        DO 40 DEG = MINDEG, NBBUCK
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
        IF (ME.LE.0) THEN
CNOCOVBEG
        write (*,*) ' INTERNAL ERROR in  MUMPS_HAMF4 '
        write(6,*) ' NEL, DEG=', NEL,DEG
        NCMPA = -NORIG
CNOCOVEND
        ENDIF
       IF (DEG.GT.N) THEN
C        -------------------------------
C        Linear search to find variable
C        with best score in the list
C        -------------------------------
C        While end of list list not reached
C         NEXT(J) = 0
         J = NEXT(ME)
         K = WF(ME)
   55    CONTINUE
         IF (J.GT.0) THEN
          IF (WF(J).LT.K) THEN
           ME = J
           K  = WF(ME)
          ENDIF
          J= NEXT(J)
          GOTO 55
         ENDIF
         ILAST = LAST(ME)
         INEXT = NEXT(ME)
         IF (INEXT .NE. 0) LAST (INEXT) = ILAST
         IF (ILAST .NE. 0) THEN
           NEXT (ILAST) = INEXT
         ELSE
C          me is at the head of the degree list
           HEAD (DEG) = INEXT
         ENDIF
C
        ELSE
C         -------------------------------------------------------------
C         remove chosen variable from link list
C         -------------------------------------------------------------
          INEXT = NEXT (ME)
          IF (INEXT .NE. 0) LAST (INEXT) = 0
          HEAD (DEG) = INEXT
        ENDIF
C       -------------------------------------------------------------
C       me represents the elimination of pivots nel+1 to nel+nv(me).
C       place me itself as the first in this set.  It will be moved
C       to the nel+nv(me) position when the permutation vectors are
C       computed.
C       -------------------------------------------------------------
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
C=======================================================================
C  CONSTRUCT NEW ELEMENT
C=======================================================================
C       -------------------------------------------------------------
C       At this point, me is the pivotal supervariable.  It will be
C       converted into the current element.  Scan list of the
C       pivotal supervariable, me, setting tree pointers and
C       constructing new list of supervariables for the new element,
C       me.  p is a pointer to the current position in the old list.
C       -------------------------------------------------------------
C       flag the variable "me" as being in Lme by negating nv (me)
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
C         ----------------------------------------------------------
C         construct the new element in place
C         ----------------------------------------------------------
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + LEN (ME) - 1
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
C             ----------------------------------------------------
C             i is a principal variable not yet placed in Lme.
C             store i in new list
C             ----------------------------------------------------
              DEGME = DEGME + NVI
C             flag i as being in Lme by negating nv (i)
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
              IF (DEGREE(I).NE.N2) THEN
C             ----------------------------------------------------
C             remove variable i from degree list. (only if i \in V0)
C             ----------------------------------------------------
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
C               i is at the head of the degree list
CNversion0              DEG = min(WF(I),NBBUCK)
CNversion1
                IF (WF(I).GT.NORIG) THEN
                 DEG = min(((WF(I)-NORIG)/PAS) + NORIG, NBBUCK)
                ELSE
                 DEG = WF(I)
                ENDIF
                HEAD (DEG) = INEXT
              ENDIF
              ENDIF
            ENDIF
   60     CONTINUE
C         this element takes no new memory in iw:
          NEWMEM = 0
        ELSE
C         ----------------------------------------------------------
C         construct the new element in empty space, iw (pfree ...)
C         ----------------------------------------------------------
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          KNT1_UPDATED = 0
          DO 120 KNT1 = 1, ELENME + 1
            KNT1_UPDATED = KNT1_UPDATED +1
            IF (KNT1 .GT. ELENME) THEN
C             search the supervariables in me.
              E = ME
              PJ = P
              LN = SLENME
            ELSE
C             search the elements in me.
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
C           -------------------------------------------------------
C           search for different supervariables and add them to the
C           new list, compressing when necessary. this loop is
C           executed once for each element in the list and once for
C           all the supervariables in the list.
C           -------------------------------------------------------
            KNT2_UPDATED = 0
            DO 110 KNT2 = 1, LN
              KNT2_UPDATED = KNT2_UPDATED+1
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
C               -------------------------------------------------
C               compress iw, if necessary
C               -------------------------------------------------
                IF (PFREE .GT. IWLEN) THEN
C                 prepare for compressing iw by adjusting
C                 pointers and lengths so that the lists being
C                 searched in the inner and outer loops contain
C                 only the remaining entries.
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1_UPDATED
C                 Reset KNT1_UPDATED in case of recompress
C                 at same iteration of the loop 120
                  KNT1_UPDATED = 0
C                 Check if anything left in supervariable ME
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0_8
                  PE (E) = PJ
                  LEN (E) = LN - KNT2_UPDATED
C                 Reset KNT2_UPDATED in case of recompress
C                 at same iteration of the loop 110
                  KNT2_UPDATED = 0
C                 Check if anything left in element E
                  IF (LEN (E) .EQ. 0) PE (E) = 0_8
                  NCMPA = NCMPA + 1
C                 store first item in pe
C                 set first entry to -item
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = int(IW (PN),8)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
C                 psrc/pdst point to source/destination
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
C                 while loop:
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
C                   search for next negative entry
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = int(PE (J))
                      PE (J) = PDST
                      PDST = PDST + 1_8
C                     copy from source to destination
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
C                 move the new partially-constructed element
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
C               -------------------------------------------------
C               i is a principal variable not yet placed in Lme
C               store i in new list
C               -------------------------------------------------
                DEGME = DEGME + NVI
C               flag i as being in Lme by negating nv (i)
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
              IF (DEGREE(I).NE.N2) THEN
C               -------------------------------------------------
C               remove variable i from degree link list
C                            (only if i in V0)
C               -------------------------------------------------
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
CNversion0              DEG = min(WF(I),NBBUCK)
CNversion1
                  IF (WF(I).GT.NORIG) THEN
                   DEG = min(((WF(I)-NORIG)/PAS) + NORIG, NBBUCK)
                  ELSE
                   DEG = WF(I)
                  ENDIF
C                 i is at the head of the degree list
                  HEAD (DEG) = INEXT
                ENDIF
              ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
C             set tree pointer and flag to indicate element e is
C             absorbed into new element me (the parent of e is me)
              PE (E) = int(-ME,8)
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
C         this element takes newmem new memory in iw (possibly zero)
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = max (MAXMEM, MEM)
        ENDIF
C       -------------------------------------------------------------
C       me has now been converted into an element in iw (pme1..pme2)
C       -------------------------------------------------------------
C       degme holds the external degree of new element
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = int(PME2 - PME1 + 1_8)
C       -------------------------------------------------------------
C       make sure that wflg is not too large.  With the current
C       value of wflg, wflg+n must not cause integer overflow
C       -------------------------------------------------------------
        IF (WFLG .GT. MAXINT_N) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
C=======================================================================
C  COMPUTE (w (e) - wflg) = |Le\Lme| FOR ALL ELEMENTS
C=======================================================================
C       -------------------------------------------------------------
C       Scan 1:  compute the external degrees of previous elements
C       with respect to the current element.  That is:
C            (w (e) - wflg) = |Le \ Lme|
C       for each element e that appears in any supervariable in Lme.
C       The notation Le refers to the pattern (list of
C       supervariables) of a previous element e, where e is not yet
C       absorbed, stored in iw (pe (e) + 1 ... pe (e) + iw (pe (e))).
C       The notation Lme refers to the pattern of the current element
C       (stored in iw (pme1..pme2)).   If (w (e) - wflg) becomes
C       zero, then the element e will be absorbed in scan 2.
C       -------------------------------------------------------------
        DO 150 PME = PME1, PME2
          I = IW (PME)
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
C           note that nv (i) has been negated to denote i in Lme:
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + int(ELN - 1,8)
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
C               unabsorbed element e has been seen in this loop
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
C               e is an unabsorbed element
C               this is the first we have seen e in all of Scan 1
                WE = DEGREE (E) + WNVI
                WF(E) = 0
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
C=======================================================================
C  DEGREE UPDATE AND ELEMENT ABSORPTION
C=======================================================================
C       -------------------------------------------------------------
C       Scan 2:  for each i in Lme, sum up the degree of Lme (which
C       is degme), plus the sum of the external degrees of each Le
C       for the elements e appearing within i, plus the
C       supervariables in i.  Place i in hash list.
C       -------------------------------------------------------------
        DO 180 PME = PME1, PME2
          I = IW (PME)
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0_8
          DEG  = 0
          WF3  = 0
          WF4  = 0
          NVI  = -NV(I)
C         ----------------------------------------------------------
C         scan the element list associated with supervariable i
C         ----------------------------------------------------------
          DO 160 P = P1, P2
            E = IW (P)
C           dext = | Le \ Lme |
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              IF ( WF(E) .EQ. 0 ) THEN
C              First time we meet e : compute wf(e)
C              which holds the surface associated to element e
C              it will later be deducted from fill-in
C              area of all variables adjacent to e
               WF(E) = DEXT * ( (2 * DEGREE(E))  -  DEXT - 1)
              ENDIF
              WF4 = WF4 + WF(E)
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E, kind=8)
            ELSE IF (DEXT .EQ. 0) THEN
C             aggressive absorption: e is not adjacent to me, but
C             the |Le \ Lme| is 0, so absorb it into me
              PE (E) = int(-ME,8)
              W (E) = 0
            ENDIF
  160     CONTINUE
C         count the number of elements in i (including me):
          ELEN (I) = int(PN - P1 + 1_8)
C         ----------------------------------------------------------
C         scan the supervariables in the list associated with i
C         ----------------------------------------------------------
          P3 = PN
          DO 170 P = P2 + 1_8, P1 + int(LEN (I) - 1,8)
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
C             j is unabsorbed, and not in Lme.
C             add to degree and add to new list
              DEG = DEG + NVJ
              WF3 = WF3 + NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + int(J,kind=8)
            ENDIF
  170     CONTINUE
C
          IF (DEGREE(I).EQ.N2) DEG = N2
C         ----------------------------------------------------------
C         update the degree and check for mass elimination
C         ----------------------------------------------------------
          IF (DEG .EQ. 0) THEN
C           -------------------------------------------------------
C           mass elimination
C           -------------------------------------------------------
C           There is nothing left of this node except for an
C           edge to the current pivot element.  elen (i) is 1,
C           and there are no variables adjacent to node i.
C           Absorb i into the current pivot element, me.
            PE (I) = int(-ME,8)
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
C           -------------------------------------------------------
C           update the upper-bound degree of i
C           -------------------------------------------------------
C           the following degree does not yet include the size
C           of the current element, which is added later:
            IF (DEGREE(I).NE.N2) THEN
C                I does not belong to halo
CNCCC             DEG        = min (DEG, NLEFT)
                 IF ( DEGREE (I).LT.DEG ) THEN
C                  Our appox degree is loose.
C                  we keep old value. Note that in
C                  this case we cannot substract WF(I)
C                  for min-fill score.
                   WF4 = 0
                   WF3 = 0
                 ELSE
                   DEGREE(I)  = DEG
                 ENDIF
            ENDIF
C
C           compute WF(I) taking into account size of block 3.0
            WF(I)      = WF4 + 2*NVI*WF3
C           -------------------------------------------------------
C           add me to the list for i
C           -------------------------------------------------------
C           move first supervariable to end of list
            IW (PN) = IW (P3)
C           move first element to end of element part of list
            IW (P3) = IW (P1)
C           add new element to front of list.
            IW (P1) = ME
C           store the new length of the list in len (i)
            LEN (I) = int(PN - P1 + 1)
            IF (DEG.NE.N2) THEN
C           -------------------------------------------------------
C           place in hash bucket.  Save hash key of i in last (i).
C           -------------------------------------------------------
            HASH = mod (HASH, HMOD) + 1_8
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
C             the degree list is empty, hash head is -j
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
C             degree list is not empty
C             use last (head (hash)) as hash head
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = int(HASH,kind=kind(LAST))
            ENDIF
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
C       -------------------------------------------------------------
C       Clear the counter array, w (...), by incrementing wflg.
C       -------------------------------------------------------------
        DMAX = max (DMAX, DEGME)
        WFLG = WFLG + DMAX
C       make sure that wflg+n does not cause integer overflow
        IF (WFLG .GT. MAXINT_N) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
C       at this point, w (1..n) .lt. wflg holds
C=======================================================================
C  SUPERVARIABLE DETECTION
C=======================================================================
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF ( (NV (I) .LT. 0) .AND. (DEGREE(I).NE.N2) ) THEN
C           i is a principal variable in Lme
C           -------------------------------------------------------
C           examine all hash buckets with 2 or more variables.  We
C           do this by examing all unique hash keys for super-
C           variables in the pattern Lme of the current element, me
C           -------------------------------------------------------
            HASH = int(LAST (I),kind=8)
C           let i = head of hash bucket, and empty the hash bucket
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
C             degree list is empty
              I = -J
              HEAD (HASH) = 0
            ELSE
C             degree list is not empty, restore last () of head
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
C           while loop:
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
C             ----------------------------------------------------
C             this bucket has one or more variables following i.
C             scan all of them to see if i can absorb any entries
C             that follow i in hash bucket.  Scatter i into w.
C             ----------------------------------------------------
              LN = LEN (I)
              ELN = ELEN (I)
C             do not flag the first element in the list (me)
              DO 210 P = PE (I) + 1_8, PE (I) + int(LN - 1,8)
                W (IW (P)) = WFLG
  210         CONTINUE
C             ----------------------------------------------------
C             scan every other entry j following i in bucket
C             ----------------------------------------------------
              JLAST = I
              J = NEXT (I)
C             while loop:
  220         CONTINUE
              IF (J .NE. 0) THEN
C               -------------------------------------------------
C               check if j and i have identical nonzero pattern
C               -------------------------------------------------
C               jump if i and j do not have same size data structure
                IF (LEN (J) .NE. LN) GO TO 240
C               jump if i and j do not have same number adj elts
                IF (ELEN (J) .NE. ELN) GO TO 240
C               do not flag the first element in the list (me)
                DO 230 P = PE (J) + 1_8, PE (J) + int(LN - 1,8)
C                 jump if an entry (iw(p)) is in j but not in i
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
C               -------------------------------------------------
C               found it!  j can be absorbed into i
C               -------------------------------------------------
                PE (J) = int(-I,8)
                WF(I)  = max(WF(I),WF(J))
C               both nv (i) and nv (j) are negated since they
C               are in Lme, and the absolute values of each
C               are the number of variables in i and j:
                NV (I) = NV (I) + NV (J)
                NV (J) = 0
                ELEN (J) = 0
C               delete j from hash bucket
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
C               -------------------------------------------------
  240           CONTINUE
C               j cannot be absorbed into i
C               -------------------------------------------------
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
C             ----------------------------------------------------
C             no more variables can be absorbed into i
C             go to next i in bucket and clear flag array
C             ----------------------------------------------------
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
C=======================================================================
C  RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVAR. FROM ELEMENT
C=======================================================================
        P = PME1
        NLEFT = TOTEL - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
C           i is a principal variable in Lme
C           restore nv (i) to signify that i is principal
            NV (I) = NVI
            IF (DEGREE(I).NE.N2) THEN
C           -------------------------------------------------------
C           compute the external degree (add size of current elem)
C           -------------------------------------------------------
CNamdonly            DEG = min (DEGREE (I) + DEGME - NVI, NLEFT - NVI)
C--------------------------
C--------------------------
            DEG = min (DEGREE (I) + DEGME - NVI, NLEFT - NVI)
            IF (DEGREE (I) + DEGME .GT. NLEFT ) THEN
C
              DEG = DEGREE(I)
              RMF1  = dble(DEG)*dble( (DEG-1) + 2*DEGME )
     &              - dble(WF(I))
              DEGREE(I) = NLEFT - NVI
              DEG       = DEGREE(I)
              RMF = dble(DEG)*dble(DEG-1)
     &         -  dble(DEGME-NVI)*dble(DEGME-NVI-1)
              RMF = min(RMF, RMF1)
            ELSE
              DEG = DEGREE(I)
              DEGREE(I) = DEGREE (I) + DEGME - NVI
C             All previous cliques taken into account (AMF4)
              RMF  = dble(DEG)*dble( (DEG-1) + 2*DEGME )
     &              - dble(WF(I))
            ENDIF
C
            RMF =  RMF / dble(NVI+1)
C
            IF (RMF.LT.dummy) THEN
             WF(I) = int ( anint( RMF ))
            ELSEIF (RMF / dble(N) .LT. dummy) THEN
             WF(I) = int ( anint( RMF/dble(N) ))
            ELSE
             WF(I) = idummy
            ENDIF
            WF(I) = max(1,WF(I))
CN
CN
CN
CN
CN
CN
            DEG = WF(I)
            IF (DEG.GT.NORIG) THEN
              DEG = min(((DEG-NORIG)/PAS) + NORIG, NBBUCK)
            ENDIF
            INEXT = HEAD (DEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (DEG) = I
C           -------------------------------------------------------
C           save the new degree, and find the minimum degree
C           -------------------------------------------------------
            MINDEG = min (MINDEG, DEG)
C begin HALO
              ENDIF
C end HALO
C           -------------------------------------------------------
C           place the supervariable in the element pattern
C           -------------------------------------------------------
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
C=======================================================================
C  FINALIZE THE NEW ELEMENT
C=======================================================================
        NV (ME) = NVPIV + DEGME
C       fill_est = fill_est + nvpiv * (nvpiv + 2 * degme)
C       nv (me) is now the degree of pivot (including diagonal part)
C       save the length of the list for the new element me
        LEN (ME) = int(P - PME1)
        IF (LEN (ME) .EQ. 0) THEN
C         there is nothing left of the current pivot element
          PE (ME) = 0_8
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
C         element was not constructed in place: deallocate part
C         of it (final size is less than or equal to newmem,
C         since newly nonprincipal variables have been removed).
          PFREE = P
          MEM = MEM - NEWMEM + int(LEN (ME),8)
        ENDIF
C=======================================================================
C       END WHILE (selecting pivots)
      GO TO 30
      ENDIF
C=======================================================================
C begin HALO V2
C     IF (NEL.LT.N+NBELTS) THEN
      IF (NBFLAG.GT.0) THEN
C
C     All possible pivots (not flagged have been eliminated).
C     We amalgamate all flagged variables at the root and
C     we finish the elimination tree.
C          1/ Go through all
C          non absorbed elements (root of the subgraph)
C          and absorb in ME
C          2/ perform mass elimination of all dense rows
           DO DEG = MINDEG, NBBUCK+1
             ME = HEAD (DEG)
             IF (ME .GT. 0) GO TO 51
           ENDDO
   51      MINDEG = DEG
           NELME    = -(NEL+1)
           DO X=1,N
            IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
C            X is an unabsorbed element
             PE(X) = int(-ME,8)
C            W(X) = 0 could be suppressed ?? check it
            ELSEIF (DEGREE(X).EQ.N2) THEN
C            X is a dense row, absorb it in ME (mass elimination)
             NEL   = NEL + NV(X)
             PE(X) = int(-ME,8)
             ELEN(X) = 0
C            Correct value of NV is (secondary variable)
             NV(X) = 0
            ENDIF
           ENDDO
C          ME is the root node
           ELEN(ME) = NELME
C          Correct value of NV is (principal variable)
           NV(ME)   = N-NREAL   ! =NBFLAG
           PE(ME)   = 0_8
CN
#if defined(SCOTCH_DEBUG_ORDER2)
C       IF (NEL.NE.N) THEN
        IF (NEL.NE.NORIG) THEN
         NCMPA = -NORIG - 1
         WRITE(6,*) " Warning 30Oct2019"
        ENDIF
#endif
      ENDIF
C end HALO
C
C=======================================================================
C  COMPUTE THE PERMUTATION VECTORS
C=======================================================================
C     ----------------------------------------------------------------
C     The time taken by the following code is O(n).  At this
C     point, elen (e) = -k has been done for all elements e,
C     and elen (i) = 0 has been done for all nonprincipal
C     variables i.  At this point, there are no principal
C     supervariables left, and all elements are absorbed.
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
C     compute the ordering of unordered nonprincipal variables
C     ----------------------------------------------------------------
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
C         ----------------------------------------------------------
C         i is an un-ordered row.  Traverse the tree from i until
C         reaching an element, e.  The element, e, was the
C         principal supervariable of i and all nodes in the path
C         from i to when e was selected as pivot.
C         ----------------------------------------------------------
          J = int(-PE (I))
C         while (j is a variable) do:
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = int(-PE (J))
              GO TO 270
            ENDIF
            E = J
C           ----------------------------------------------------------
C           get the current pivot ordering of e
C           ----------------------------------------------------------
            K = -ELEN (E)
C           ----------------------------------------------------------
C           traverse the path again from i to e, and compress the
C           path (all nodes point to e).  Path compression allows
C           this code to compute in O(n) time.  Order the unordered
C           nodes in the path, and place the element e at the end.
C           ----------------------------------------------------------
            J = I
C           while (j is a variable) do:
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = int(-PE (J))
              PE (J) = int(-E,8)
              IF (ELEN (J) .EQ. 0) THEN
C               j is an unordered row
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
C         leave elen (e) negative, so we know it is an element
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
C     ----------------------------------------------------------------
C     reset the inverse permutation (elen (1..n)) to be positive,
C     and compute the pivot order (last (1..n)).
C     ----------------------------------------------------------------
C=======================================================================
C  RETURN THE MEMORY USAGE IN IW
C=======================================================================
C     If maxmem is less than or equal to iwlen, then no compressions
C     occurred, and iw (maxmem+1 ... iwlen) was unused.  Otherwise
C     compressions did occur, and iwlen would have had to have been
C     greater than or equal to maxmem for no compressions to occur.
C     Return the value of maxmem in the pfree argument.
 500  PFREE = MAXMEM
      RETURN
      END SUBROUTINE MUMPS_HAMF4

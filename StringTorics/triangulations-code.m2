---------------------------------
-- Triangulation code -----------
---------------------------------
-*
findMaxRegularStarTriangulation = method(Options=>{MonomialOrder=>null})
findMaxRegularStarTriangulation Matrix := opts -> (A) -> (
    -- steps: 1. Augment A
    --        2. Compute the toric ideal
    --        3. do gfan of that
    --        4. get triangulations
    --        5. restrict to star ones (i.e. ones for which the n-th one corresp to origin appears in all simplices)
    --        return these
    n := numColumns A;
    d := numRows A;
    A1 := matrix{{matrix{{n:1}}, 1},{A,0}};
    D := symbol D;
    R := QQ[D_0..D_n];
    time IA := toBinomial(toricMarkov A1, R);
    time interiors := flatten select(IA_*/terms//flatten/support, t -> # t == 1);
    time interiors = interiors/index//set; -- select(interiors/index, i -> i < n);
    time wt1 := splice{n:1, 0};
    time wt2 := for i from 0 to n list if member(i,interiors) then 0 else 1;
    --R1 := QQ[gens R, Weights=>{wt1,wt2}];
    time monorder := if opts.MonomialOrder===null then {Weights=>wt1, Weights=>wt2}
                else {Weights=>wt1, Weights=>wt2, opts.MonomialOrder};
    time R1 := ZZ/101[gens R, MonomialOrder=>monorder];
    IA1 := gb sub(IA,R1);
    inComplex := radical monomialIdeal leadTerm IA1;
    fac := ideal facets simplicialComplex inComplex;
    sort for m in fac_* list drop((support m)/index,-1)
    )
*-

-*
Subdivision = new Type of HashTable
Triangulation - new Type of Subdivision
triangulation = method()
triangulation(Matrix, List) := (A, tri) -> (
    T := new Triangulation from { 
        symbol Matrix => A, 
        symbol AugmentedMatrix => augment A, 
        symbol Facets => tri, 
        symbol Dimension => rank A,
        cache => new CacheTable
        };
    -- check some basics:
    -- A is a matrix over ZZ or QQ
    -- tri is a set of subsets, all of size dim+1, of 0..numcols A - 1.
    T
    )

isWellDefined Triangulation := (T) -> (
    -- check that T is indeed a triangulation of the convex hull of the columns of T.Matrix.
    )

-- what is not yet written?
-- isWellDefined, fine to fine star triangulation (does method always work?)

isWellDefined Triangulation
flips Triangulation
isFine
isStar
isRegular
isFRST
regularTriangulationWeights
flip(Triangulation, List)
generateTriangulations(Triangulation, opts)
volume
volumeVector, or gkzVector.
regularFineTriangulation A -- from topcom
allTriangulations A -- from topcom
fineRegularStarTriangulation A
*-


-- TODO: starting with a point configuration, and the origin, 
--   create a regular, fine, star, triangulation of this set.
--   does this always exist?

-- TODO: this routine should be submitted to Polyhedra as a bug fix.
  regularSubdivision (Matrix,Matrix) := (M,w) -> (
      n := numColumns M;
      M = M ** QQ;
      -- Checking for input errors
      if numColumns w != numColumns M or numRows w != 1 then 
          error("The weight must be a one row matrix with number of points many entries");
      P := convexHull(M||w,matrix (toList(numRows M:{0})|{{1}}));
      vertsP := transpose drop(entries (vertices P), -1);
      LPs := transpose entries M;
      LPH := hashTable for i from 0 to #LPs - 1 list LPs#i => i;
      F := select(faces (1,P), f -> #(f#1) == 0);
      sort apply (F, f -> sort apply(f#0, v -> LPH#(vertsP#v)))
      )

-- TODO/BUG: this ASSUMES (A, tri) is a triangulation.
isFine = method()
isFine(Matrix, List) := (A, tri) -> (
    numcols A == tri//flatten//unique//length
    )

-- TODO/BUG: this ASSUMES (A, tri) is a triangulation.
isStar = method()
isStar(Matrix, List) := (A, tri) -> (
    -- assumption?  last column of A is the zero element?  (or the one the star is taken with respect to).
    origin := numcols A - 1;
    all(tri, t -> member(origin, t))
    )

delaunayWeights = method()
delaunayWeights Matrix := (A) -> (
    matrix{for i from 0 to numcols A - 1 list ((transpose A_{i}) * A_{i})_(0,0)}
    )

delaunaySubdivision = method()
delaunaySubdivision Matrix := A -> elapsedTime regularSubdivision(A, elapsedTime delaunayWeights A)

-- isRegularTriangulation is defined in Topcom.

isWellDefinedTriangulation = method() -- WRITE THIS! -- perhaps this is in Topcom.m2...

-- TODO: why codim2 and codim2s?
codim2 = (tri) -> (
    -- find all pairs whose intersection is all but one of the vertices.
    n := #tri_0;
    select(subsets(tri, 2), v -> length toList ((set (v#0)) * (set (v#1))) == n-1)
    )

codim2s = (tri) -> (
    n := #tri_0;
    C2 := unique apply(subsets(tri, 2), v -> sort toList (set v#0 + set v#1));
    select(C2, x -> #x == n+1)
    )

affineCircuits = method()

-- probably remove this one.
affineCircuits(List, List) := (pts, tri) -> (
    -- pts is a list of points (in homogeneous form?)
    -- tri is a triangulation of this configuration
    -- returns: a minimal list of (sorted) affine circuits
    Amat := transpose matrix pts;
    c2 := codim2s tri;
    for c in c2 list (
        z := flatten entries syz Amat_c;
        {c_(positions(z, zi -> zi > 0)), c_(positions(z, zi -> zi < 0))}
        )
    )

-- remove this one
affineCircuits(Matrix, List) := (Amat, tri) -> (
    -- tri: a set of subsets of size d+1 from 0..#columns(Amat)-1, that form
    --  a triangulation of conv(columns of Amat).
    -- Amat: d x n matrix over ZZ or QQ of the points, with the d rows independent.
    --   OR of size(d+1) x n, with the last row all 1's. (or the rows are rank d).
    -- output: 
    c2 := codim2s tri;
    for c in c2 list (
        z := flatten entries syz Amat_c;
        {c_(positions(z, zi -> zi > 0)), c_(positions(z, zi -> zi < 0))}
        )
    )

-- TODO: RENAME to flips, or possibleFlips
affineCircuits(Matrix, List) := (Amat, tri) -> (
    -- tri: a set of subsets of size d from 0..#columns(Amat)-1, that form
    --  a triangulation of conv(columns of Amat).
    -- Amat: d-1 x n matrix over ZZ or QQ of the points, with the d-1 rows independent.
    --   OR of size d x n, with the last row all 1's. (or the rows are rank (d-1)).
    d := #(tri#0);
    assert all(tri, t -> #t == d);
    if numrows Amat == d-1 then
      Amat = Amat || matrix{ toList(numcols Amat : 1) };
    c2 := codim2s tri;
    unique for c in c2 list (
        z := flatten entries syz Amat_c;
        rsort {c_(positions(z, zi -> zi > 0)), c_(positions(z, zi -> zi < 0))}
        )
    )

-- TODO: write an `affineCircuits Amat` function.
-- Allow Homogenize option?  Use Topcom here?

link(List, List) := (tau, triangulation) -> (
    S := select(triangulation, t -> isSubset(tau, t));
    sort for s in S list sort toList (set s - set tau)
    )

flip(List,List) := (tri, affineCircuit) -> (
    -- returns a new triangulation
    -- input: tri: a list of lists of integers (a triangulation).
    --        affineCircuit: a pair of lists of integers.  Each is disjoint, the negative/positive part of the affine circuit.
    -- a. start with an affine circuit
    -- 1. compute 2 choices of triangulation for the circuit
    whole := sort unique flatten affineCircuit;
    T1 := for i in affineCircuit#0 list sort toList (set whole - set {i});
    T2 := for i in affineCircuit#1 list sort toList (set whole - set {i});
    -- 2. which is contained in tri?
    -- XXXXXXXXXXXX
    S1 := unique for tau in T1 list link(tau, tri);
    S2 := unique for tau in T2 list link(tau, tri);
    if #S1 > 1 or #S2 > 1 then (
        --<< "link is not unique, here they are: S1: " << S1 << " S2: " << S2 << endl;
        return null;
        );
    S1 = first S1;
    S2 = first S2;
    if #S1 > 0 and #S2 > 0 then (
        error "my logic is wrong: somehow both links exist...";
        return null;
        );
    if #S1 == 0 then (
        S1 = S2;
        (T1,T2) = (T2,T1);
        );
    -- Now we can change the triangulation:
    T := set tri;
    outgoing := set flatten for s in S1 list for tau in T1 list sort join(s,tau);
    incoming := set flatten for s in S1 list for tau in T2 list sort join(s,tau);
    --<< "outgoing: " << sort toList outgoing << endl;
    --<< "incoming: " << sort toList incoming << endl;    
    sortTriangulation toList(T - outgoing + incoming)
    )

generateTriangulations = method(Options => {Limit=>infinity, Regular=>false})
generateTriangulations(Matrix, List) := opts -> (Amat, TRI) -> (
    TRI = sortTriangulation TRI;
    allT := new MutableHashTable;
    allT#TRI = true;
    TODO := {TRI};
    while #TODO > 0 and #(keys allT) < opts.Limit do (
        nextTRI := TODO#0;
        TODO = drop(TODO,1);
        flips := select(affineCircuits(Amat, nextTRI), z -> #z#0 > 1 and #z#1 > 1);
        fliptris := for f in flips list flip(nextTRI, f);
        newT := select(fliptris, x -> x =!= null);
        for T in newT do (
            if not allT#?T then (
                --<< "new triangulation: " << T << endl;
                if not opts.Regular or isRegularTriangulation(Amat, T) then (
                    allT#T = true;
                    TODO = append(TODO, T);
                    );
                ));
        if debugLevel > 0 then 
            << "todo = " << #TODO << " and #triang = " << #(keys allT) << endl;
        );
    keys allT
    )
generateTriangulations Matrix := opts -> Amat -> (
    generateTriangulations(Amat, regularFineTriangulation Amat, opts))

-- Not working yet.  Make sure this is correct: is it generating the correct kind of triangulations?
-- 15 May 2018.  REMOVE? (MES, 5/1/2020)
generateTriangulations Polyhedron := opts -> P2 -> (
    (LP,tri) := regularStarTriangulation(dim P2-2,P2);
    Amat := transpose matrix LP;
    allT := new MutableHashTable;
    allT#tri = true;
    TODO := {tri};
    while #TODO > 0 and #(keys allT) < opts.Limit do (
        nextTRI := TODO#0;
        TODO = drop(TODO,1);
        flips := select(affineCircuits(Amat, nextTRI), z -> #z#0 > 1 and #z#1 > 1);
        fliptris := for f in flips list flip(nextTRI, f);
        newT := select(fliptris, x -> x =!= null);
        for T in newT do (
            if not allT#?T then (
                --<< "new triangulation: " << T << endl;
                if not opts.Regular or isRegularTriangulation(Amat, T) then (
                    allT#T = true;
                    TODO = append(TODO, T);
                    );
                ));
        << "todo = " << #TODO << " and #triang = " << #(keys allT) << endl;
        );
    keys allT
    )

-- TODO: this function is not working correctly, I believe.  Remove it.
--  replace it with a function `isTriangulation(A, tri)`.
checkFan = (Amat, tri) -> (
    -- returns null if this is a (complete) triangulation...
    -- otherwise gives an error
    F := fan for t in tri list posHull Amat_t;
    if not isComplete F then error "triangulation is a fan, but is not complete";
    )

volumeVector = method()
volumeVector(Matrix, List) := (Amat, tri) -> (
    if #tri == 0 then error "expected at least one simplex";
    nelems := #tri#0;
    d := nelems-1;
    if not all(tri, f -> #f == nelems)
    then error "expected a triangulation";
    if numrows Amat =!= nelems then Amat = Amat || matrix{{(numcols Amat):1}};
    if numrows Amat =!= nelems then error "triangulation not compatible with matrix";
    H := hashTable for t in tri list t => (abs det Amat_t)/d!;
    if debugLevel > 0 then
    << "Volume = " << (sum values H)/d! << endl;
    for i from 0 to numColumns Amat - 1 list (
        T := select(tri, t -> member(i,t));
        sum for t in T list H#t
        )
    )

-- NOTE!!!!! This code is being migrated to the package Triangulations.
-- 
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


-- isRegularTriangulation is defined in Triangulations and Topcom.

isWellDefinedTriangulation = method() -- WRITE THIS! -- perhaps this will be in Topcom or Triangulations.

-- TODO: this function is not working correctly, I believe.  Remove it.
--  replace it with a function `isTriangulation(A, tri)`.
checkFan = (Amat, tri) -> (
    -- returns null if this is a (complete) triangulation...
    -- otherwise gives an error
    F := fan for t in tri list posHull Amat_t;
    if not isComplete F then error "triangulation is a fan, but is not complete";
    )


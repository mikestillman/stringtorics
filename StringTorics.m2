-- TODO (Nov 2019)
-- 1. documentation
-- 2. Extra polyhedral functions (do we still need these?  Probably...)
--    Make sure they work with non reflexive polytopes too!
--    Or say they onyl work for reflexives?
-- 3. Triangulations
--    Computing FRST's
-- 4. Cohomology of line bundles on a toric variety
-- 5. Complete Intersections in torics
-- 6. Intersection rings (especially intersection numbers).
--      In particular, I think there is a bug in computing intersection numbers.
--      Also, handle non-favorable case.
-- 7. Mori cones

newPackage(
        "StringTorics",
        Version => "0.5", 
        Date => "26 May 2021",
        Authors => {
            {Name => "Mike Stillman", 
            Email => "mike@math.cornell.edu", 
            HomePage => "http://pi.math.cornell.edu/~mike"}
        },
        Headline => "toric variety functions for string theory",
        DebuggingMode => true,
        AuxiliaryFiles => true,
        PackageExports => {
            "FourTiTwo", -- where is this used?
            "SimplicialComplexes",
            "NormalToricVarieties",
            "Schubert2",
            "Polyhedra",
            "ReflexivePolytopesDB",
            "CohomCalg",
            "Topcom"
            },
        PackageImports => {
            --"Graphs", 
            "LLLBases"}
        )

export {
    -- Types defined here
    "ReflexivePolytope",
    "CalabiYauInToric",
    "TopologicalDataOfCY3",
    
    -- The following should be placed into ReflexivePolytopesDB:
--    "matrixFromKSEntry", -- replace this with matrixFromKS.
--    "matrixFromKS",

    -- Extra polyhedral facilities, for lattice points and faces of a Polyhedron
    -- how much of this shoiuld be exported??
    "vertexMatrix",
    "vertexList",
    "faceDimensionHash",
    "faceList",
    "dualFace",
    "minimalFace",
    "latticePointList",
    "latticePointHash",
    "interiorLatticePointList",
    "annotatedFaces",

    -- Triangulation code
    "Origin",
    "pointConfiguration",
    "regularStarTriangulation",
    "generateTriangulations",
    "isFine",
    "isStar",
    "volumeVector",
    "delaunaySubdivision",
    "delaunayWeights",
    "findAllFRSTs",
    
    -- This set maybe should be included in NormalToricVarieties?
    "singularCones",
    "singularLocusInToric",
    "normalToricVarietyFromGLSM",
    -- These should stay here
    "reflexiveToSimplicialToricVariety",
    "reflexiveToSimplicialToricVarietyCleanDegrees",
    "allZeros",
    "augmentWithOrigin",

    "readSageTriangulations",
    "sortTriangulation",
    "matchNonZero",
    "applyPermutation",
    "affineCircuits",
    "checkFan",
    "sageTri",
    "Regular",

    -- Reflexive polytope code
    -- Uses ReflexivePolytope, a wrapper over Polyhedra package, and containing the info we want/need.
    "reflexivePolytope",
    "basisIndices",
    "topologicalData",
    
    -- gvInvariants
    "toricMoriCone",
    "moriCone",
    "moriConeByGV",
    "intersectionNumbersOfCY", -- possibly not for export
    "gvInvariants",
    "invariants", -- really in the topology section...
    
    -- CompleteIntersectionInToric's
    "completeIntersection",
    "CompleteIntersectionInToric",
        "Ambient",
        "CI",

    -- Cohomology
    "toricCohomologySetup",
    "cohomologyBasis",
    "CohomologySetup",
    "cohomologyMatrix",
    "cohomologyMatrixRank",
    "genericCohomologyMatrix",
    
    -- new attempt at faster cohomology of line bundles in toric varieties.
    -- want the Laurent monomials.
    "toricOrthants",
    "H1orthants",
    "normalDegrees",
    "normalDegree",
    "setOfColumns",
    "getFraction",
    "findTope",
    "cohomologyFractions",
    
    -- cohomology for complete intersections in torics
    --TODO clean up this interface!
    "cohomologyVector",  -- currently calls CohomCalg, but shouold take an optional argument?
    "cohomologyFromLES",
    "collectLineBundles",
    "cohomologyOmega1",
    "removeUnusedVariables",
    "example0912'3524",
    "exampleP111122'44",
    "cohom",
    "nextBundle",
    "basicCohomologies",
    "hodgeDiamond",
    "hodgeVector",
    
    -- CY hypersurface cohomology info from polyhedral information
    "subcomplex",
    "complexWithInteriorFacesRemoved",
    -- The following are the ones we care about
    "hodgeCY", -- debugging only?
    "hodgeOfCYToricDivisor",
    "hodgeOfCYToricDivisors",
    "h11OfCY",
    "h21OfCY",
    "isFavorable",
    -- new formula
    "hodgeVectorViaTheorem", -- TODO: is likely not correct currently.
    "tentativeHodgeVector", -- deprecated
    
    -- Topology of a CY3-fold.
    -- input: reflexive poytope, and triangulation.
    -- Some functions to have here:
    -- type: CY3ToricHypersurface.
    -- h11 X
    -- h12, h21 X
    -- intersectionNumbers X
    -- topology X (returns an object of TopologicalDataOfCY3)
    
    -- Gopakumar-Vafa invariants
    -- computed using Andres' computeGV code in C++
    
    -- Flop chains, Mori cones

    "FilePrefix",
    "Executable",
    "Mori"    
    }

--- kludge to access parts of the 'Core'
hasAttribute = value Core#"private dictionary"#"hasAttribute";
getAttribute = value Core#"private dictionary"#"getAttribute";
ReverseDictionary = value Core#"private dictionary"#"ReverseDictionary";

------------------------------------
-- New types -----------------------
------------------------------------

ReflexivePolytope = new Type of HashTable -- contains: data about a reflexive polytope.

CalabiYauInToric = new Type of HashTable -- currently a hypersurface, eventually a CI.
  -- contains ReflexivePolytope, and a triangulation.  

-- deprecate this one
ToricHypersurface = new Type of HashTable
  -- contains ReflexivePolytope, and a triangulation.  
  -- TODO: better name? perhaps TriangulatedReflexivePolytope?

TopologicalDataOfCY3 = new Type of HashTable
  -- contains h11, h21, c2, cubic intersection form


load (currentFileDirectory | "StringTorics/MyPolyhedra.m2")
load (currentFileDirectory | "StringTorics/ToricCompleteIntersections.m2")
load (currentFileDirectory | "StringTorics/triangulations-code.m2")

protect nextVar
protect nextFinalVar
---------------------------------------------------
singularLocusInToric = method()
singularLocusInToric(Ideal, Ideal) := (I, B) -> (
    -- I: an ideal in the Cox ring, where B is the irrelevant ideal.
    -- return: a list of (codim Sing locus, groebner basis) for each max cone
    -- TODO: this I think only works in the smooth case so far.  FIX THIS.
    -- find singular locus on each open set.
    c := codim I;
    Is := for b in B_* list (
      sub(I, for x in support b list x => 1)
      );
    for J in Is list (
        singJ := J + minors(c, jacobian J);
        (codim singJ, gens gb singJ)
        )
    )

---------------------------------------------------

allzeros0 = (F, topvar, pt, lo, hi) -> (
  -- F is a polynomial in a poly ring (e.g. QQ[a,b,c])
  -- topvar
  if topvar == -1 then (return if F == 0 then {pt} else {});
  R := ring F;
  flatten for v from lo to hi list (
      G := sub(F, R_topvar => v);
      pt1 := prepend(v, pt);
      allzeros0(G, topvar-1, pt1, lo, hi)
      )
  )
allzeros1 = (F, topvar, pt, lo, hi, f) -> (
  -- F is a polynomial in a poly ring (e.g. QQ[a,b,c])
  -- topvar
  if topvar == -1 then (if F == 0 then f pt; return null);
  R := ring F;
  for v from lo to hi do (
      G := sub(F, R_topvar => v);
      pt1 := prepend(v, pt);
      allzeros1(G, topvar-1, pt1, lo, hi, f)
      )
  )
allZeros = method()
allZeros(RingElement, Sequence) := (F, lohi) -> (
    (lo,hi) := lohi;
    allzeros0(F, numgens ring F - 1, {}, lo, hi)
    )
allZeros(RingElement, ZZ, Sequence) := (F, topvar, lohi) -> (
    (lo,hi) := lohi;
    allzeros0(F, topvar, {}, lo, hi)
    )
allZeros(RingElement, ZZ, Sequence, Function) := (F, topvar, lohi, f) -> (
    (lo,hi) := lohi;
    allzeros1(F, topvar, {}, lo, hi, f)
    )

reflexiveToSimplicialToricVariety = method(Options => {
        CoefficientRing => QQ, 
        Variable =>  getSymbol "x"
        })

reflexiveToSimplicialToricVariety Polyhedron := opts -> (P1) -> (
    -- P1 is a reflexive polytope in the M lattice.
    -- Creates a simplicial toric variety via a triangulation
    -- of the normal fan of P1, using all of lattice points of (polar P1).
    -- Returns: a normal toric variety, whch uses all
    -- lattice points in the triangulation.
    P2 := polar P1;
    (LP,tri) := regularStarTriangulation(dim P2-2,P2);
    normalToricVariety(LP,tri,opts)
    )

-- Subroutines for reflexiveToSimplicialToricVarietyCleanDegrees
-- These might be generally useful too?
  findFirstUnitVectors = method()
  findFirstUnitVectors Matrix := List => (M) -> (
      -- returns the indices of the the columns which are unit vectors,
      -- if there are two or more columns corresponding to the same unit vector, take the first.
      -- the result is in increasing list of integers.
      e := entries transpose M;
      unitposition := (col) -> if all(col, a -> a >= 0) and sum col === 1 then position(col, a -> a === 1) else null;
      H := partition(c -> unitposition e_c, toList(0..numcols M-1));
      ks := select(keys H, k -> k =!= null);
      sort for k in ks list min(H#k)
      )

  -- extend the nice unit vectors to a list p of n columns (n = numrows M)
  -- so that det M_p = 1 or -1.
  -- compute the inverse A of this n x n matrix.
  -- then compute A^-1 * M, and return this list of columns.
  -- TODO: don't use subsets...
  findInvertibleSubmatrix = method(Options => {Limit => 10000})
  findInvertibleSubmatrix(Matrix, List) := List => opts -> (M, p) -> (
      -- M is a matrix over the integers
      -- p is a list of column indices of M (indices run from 0 to numcols M - 1)
      -- Find a list q containing p (if possible), of column indices,
      --   s.t. det(M_q) = 1 or -1.
      -- Return q, or null, if one cannot be found.
      S := set toList(0..numcols M - 1);
      others := sort toList(S - set p);
      needed := numrows M - #p;
      if binomial(#others, needed) > opts.Limit then return null;
      -- really: do some random choices of q containing p.
      trythese := subsets(others, needed);
      tried := 0;
      for q1 in trythese do (
          q := join(q1,p);
          tried = tried+1;
          if abs det(M_q) == 1 then (
              if debugLevel >= 1 then (
                  << "found suitable set of columns after " << tried << " step" 
                  << if tried == 1 then "" else "s" << endl;
                  );
              return sort q;
              );
          );
      << "tried all determinants: none had unit determinant" << endl;
      null
      )

reflexiveToSimplicialToricVarietyCleanDegrees = method(Options => options reflexiveToSimplicialToricVariety)
reflexiveToSimplicialToricVarietyCleanDegrees Polyhedron := Sequence => opts -> (P1) -> (
    -- returns (V, q) where
    -- V is essentially reflexiveToSimplicialToricVariety P1
    --   except that the degrees of the ring have been cleaned up
    -- AND 
    -- a list of column indices (or of rays) in ascending order, 
    -- whose degrees form a basis, in fact the identity matrix.
    -- null is returned if no such q can be found, or the algorithm
    -- doesn't find it.
    -- Note: we don't really need ring V0 here...
    P2 := polar P1;
    (LP,tri) := regularStarTriangulation(dim P2-2,P2);
    D := transpose syz matrix transpose LP;
    print D;
    p := findFirstUnitVectors D;
    q := findInvertibleSubmatrix(D, p);
    if q === null then null
    else (
        A := D_q;
        D' := A^-1 * D;
        (normalToricVariety(LP, tri, WeilToClass => D', opts), q)
        )
    )

augmentWithOrigin = method()
augmentWithOrigin Matrix := (A) -> (
    -- A is a matrix over ZZ
    -- add column of 0's, then add in a first row of 1's.
    n := numColumns A;
    zeros := matrix {numRows A : {0}};
    ones := matrix {{(n+1) : 1}};
    ones || (A | zeros)
    )

augmentBack = (A) -> (
    -- A is a matrix over ZZ
    -- add in a last row of 1's.  Also add in one column of zeros
    n := numColumns A;
    zeros1 := matrix{numRows A : {0}};
    zeros0 := matrix {1+numRows A : {0}};
    ones := matrix {{n+1 : 1}};
    ((A | zeros1) || ones) --| zeros0
    )

readSageTriangulations = method()
readSageTriangulations String := (str) -> (
   --  "     [[array([0, 1, 2, 3]), array([0, 1, 3, 4]), array([1, 2, 3, 4]), array([2, 3, 4, 5]), array([2, 3, 5, 8]), array([2, 5, 6, 7]), array([2, 5, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 7]), array([0, 2, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 4, 6]), array([2, 4, 5, 6])], [array([0, 1, 2, 3]), array([0, 1, 3, 4]), array([1, 2, 3, 5]), array([1, 3, 4, 5]), array([2, 3, 5, 8]), array([2, 5, 6, 7]), array([2, 5, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 7]), array([0, 2, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 5, 6]), array([1, 4, 5, 6])], [array([0, 1, 2, 3]), array([0, 1, 3, 4]), array([1, 2, 3, 5]), array([1, 3, 4, 5]), array([2, 3, 5, 8]), array([2, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 5, 6]), array([1, 4, 5, 6])], [array([0, 1, 2, 3]), array([0, 1, 3, 4]), array([1, 2, 3, 4]), array([2, 3, 4, 5]), array([2, 3, 5, 8]), array([2, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 4, 6]), array([2, 4, 5, 6])], [array([0, 1, 2, 3]), array([0, 1, 3, 4]), array([1, 2, 3, 5]), array([1, 3, 4, 5]), array([2, 3, 5, 6]), array([3, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 6]), array([0, 3, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 5, 6]), array([1, 4, 5, 6])], [array([0, 1, 2, 5]), array([0, 1, 4, 5]), array([0, 2, 3, 5]), array([0, 3, 4, 5]), array([2, 3, 5, 8]), array([2, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 5, 6]), array([1, 4, 5, 6])], [array([0, 1, 2, 3]), array([0, 1, 3, 4]), array([1, 2, 3, 4]), array([2, 3, 4, 5]), array([2, 3, 5, 6]), array([3, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 6]), array([0, 3, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 4, 6]), array([2, 4, 5, 6])], [array([0, 1, 2, 4]), array([0, 2, 3, 5]), array([0, 2, 4, 5]), array([0, 3, 4, 5]), array([2, 3, 5, 8]), array([2, 5, 6, 7]), array([2, 5, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 7]), array([0, 2, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 4, 6]), array([2, 4, 5, 6])], [array([0, 1, 2, 5]), array([0, 1, 4, 5]), array([0, 2, 3, 5]), array([0, 3, 4, 5]), array([2, 3, 5, 6]), array([3, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 6]), array([0, 3, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 5, 6]), array([1, 4, 5, 6])], [array([0, 1, 2, 5]), array([0, 1, 4, 5]), array([0, 2, 3, 5]), array([0, 3, 4, 5]), array([2, 3, 5, 8]), array([2, 5, 6, 7]), array([2, 5, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 7]), array([0, 2, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 5, 6]), array([1, 4, 5, 6])], [array([0, 1, 2, 4]), array([0, 2, 3, 5]), array([0, 2, 4, 5]), array([0, 3, 4, 5]), array([2, 3, 5, 8]), array([2, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 4, 6]), array([2, 4, 5, 6])], [array([0, 1, 2, 4]), array([0, 2, 3, 5]), array([0, 2, 4, 5]), array([0, 3, 4, 5]), array([2, 3, 5, 6]), array([3, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 6]), array([0, 3, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 4, 6]), array([2, 4, 5, 6])]] "
    s1 := replace("array", "", str);
    s2 := replace("\\(", "", s1);    
    s3 := replace("\\)", "", s2);
    s4 := replace("\\[", "{", s3);
    s5 := replace("\\]", "}", s4);
    value s5
    )

sortTriangulation = method()
sortTriangulation List := (T) -> sort for t in T list sort t

matchNonZero = method()
matchNonZero(Matrix, Matrix) := (A,B) -> (
    eA := entries transpose A;
    Ah := hashTable for i from 0 to #eA-1 list (
        if all(eA#i, x -> x == 0) then continue else eA#i => i
        );
    eB := entries transpose B;
    Bh := hashTable for i from 0 to #eB-1 list (
        if all(eB#i, x -> x == 0) then continue else eB#i => i
        );
    if set keys Ah =!= set keys Bh then (
        error "non-zero columns are different";
        );
    AtoB := for e in eA list if Bh#?e then Bh#e else -1;
    BtoA := for e in eB list if Ah#?e then Ah#e else -1;
    (AtoB, BtoA)
    )

applyPermutation = method()
applyPermutation(List, ZZ) := (P, i) -> if i >= 0 and i < #P then P#i else -1
applyPermutation(List, List) := (P, L) -> sort for f in L list applyPermutation(P,f)

---------------------------------------
-- Code for toric varieties packages --
---------------------------------------

  singularCones = method()
  singularCones(ZZ, NormalToricVariety) := (coneDim, X) -> (
      rys := transpose matrix rays X;
      which := orbits(X, coneDim);
      select(which, f -> (
              I := trim minors(#f, rys_f);
              I != 1))
      )
  singularCones(ZZ, NormalToricVariety, Thing) := (coneDim, X, notused) -> (
      -- This version also returns the multiplicity, for non smooth cones.
      rys := transpose matrix rays X;
      which := orbits(X, coneDim);
      for f in which list (
          I := trim minors(#f, rys_f);
          if I == 1 then continue;
          f => I_0
          )
      )

  normalToricVarietyFromGLSM = method(Options=>options normalToricVariety)
  normalToricVarietyFromGLSM(Matrix, List) := opts -> (GLSM, maxCones) -> (
      -- create rays from GLSM degrees (which will then be the output of 
      rays := entries syz GLSM;
      normalToricVariety(rays, maxCones, opts, WeilToClass => GLSM)
      )
  normalToricVarietyFromGLSM(Matrix, MonomialIdeal) := opts -> (GLSM, SRIdeal) -> (
      -- this switch to dual should be much simpler than this...
      n := numColumns GLSM; -- also should be the number of variables of the ring SRIdeal
      allvars := set toList(0..n-1);
      maxcones := (dual SRIdeal)_*/(f -> sort toList(allvars - (support f)/index//set));
      normalToricVarietyFromGLSM(GLSM, maxcones, opts)
      )

  -- This function isn't used anywhere??
  pointConfiguration = method(Options => {Homogenize=>false, Origin => false});
  pointConfiguration(ZZ,Polyhedron) := opts -> (maxdim, P2) -> (
      LP := select(latticePointList P2, lp -> dim(P2, minimalFace(P2, lp)) <= maxdim);
      A := transpose matrix LP;
      if opts.Origin then 
        A = A | matrix for i from 1 to numRows A list {0};
      if opts.Homogenize then
        A = A || matrix {for i from 1 to numColumns A list 0};
      A
      )
  
regularStarTriangulation = method()
regularStarTriangulation Polyhedron := (P2) -> (
    LP := drop(latticePointList P2, -1);
    A := transpose matrix LP;
    tri := regularFineTriangulation A;
    -- Now, this is not a star triangulation...  But we think we can make it so in the
    -- following way.
    facetlist := (faceList(dim P2-1, P2))/(f -> latticePointList(P2,f));
    newtri := sort flatten for f in tri list (
        fac := select(facetlist, g -> #((set f) * (set g)) == dim P2);
        for g in fac list sort toList ((set f) * (set g))
        );
    (LP,newtri)
    )

regularStarTriangulation(ZZ,Polyhedron) := (maxdim, P2) -> (
    LP := select(latticePointList P2, lp -> dim(P2, minimalFace(P2, lp)) <= maxdim);
    A := transpose matrix LP;
    A = A | matrix for i from 1 to numRows A list {0};
    tri := regularFineTriangulation A;
    -- Now, this is not a star triangulation...  But we think we can make it so in the
    -- following way.
    H := latticePointHash P2;
    fullset := set for lp in LP list H#lp;
    facetlist := faceList(dim P2-1, P2);
    facetlist = facetlist/(f -> (
            sort toList(set latticePointList(P2,f) * fullset)
            ));
    newtri := sort flatten for f in tri list (
        fac := select(facetlist, g -> #((set f) * (set g)) == dim P2);
        for g in fac list sort toList ((set f) * (set g))
        );
    (LP,newtri)
    )

  sageTri = "     [[array([0, 1, 2, 3]), array([0, 1, 3, 4]), array([1, 2, 3, 4]), array([2, 3, 4, 5]), 
array([2, 3, 5, 8]), array([2, 5, 6, 7]), array([2, 5, 7, 8]), 
array([0, 2, 3, 8]), array([0, 2, 6, 7]), array([0, 2, 7, 8]), 
array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), 
array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), 
array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), 
array([1, 2, 4, 6]), array([2, 4, 5, 6])], [array([0, 1, 2, 3]), 
array([0, 1, 3, 4]), array([1, 2, 3, 5]), array([1, 3, 4, 5]), 
array([2, 3, 5, 8]), array([2, 5, 6, 7]), array([2, 5, 7, 8]), 
  array([0, 2, 3, 8]), array([0, 2, 6, 7]), array([0, 2, 7, 8]),   
  array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), 
  array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), 
  array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), 
  array([1, 2, 5, 6]), array([1, 4, 5, 6])], [array([0, 1, 2, 3]), 
  array([0, 1, 3, 4]), array([1, 2, 3, 5]), array([1, 3, 4, 5]), 
  array([2, 3, 5, 8]), array([2, 5, 6, 8]), array([5, 6, 7, 8]), 
  array([0, 2, 3, 8]), array([0, 2, 6, 8]), array([0, 6, 7, 8]), 
  array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), 
  array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), 
  array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), 
  array([1, 2, 5, 6]), array([1, 4, 5, 6])], [array([0, 1, 2, 3]), 
  array([0, 1, 3, 4]), array([1, 2, 3, 4]), array([2, 3, 4, 5]), 
  array([2, 3, 5, 8]), array([2, 5, 6, 8]), array([5, 6, 7, 8]), 
  array([0, 2, 3, 8]), array([0, 2, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), 
  array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), 
  array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 4, 6]), array([2, 4, 5, 6])], [array([0, 1, 2, 3]), 
  array([0, 1, 3, 4]), array([1, 2, 3, 5]), array([1, 3, 4, 5]), array([2, 3, 5, 6]), array([3, 5, 6, 8]), 
  array([5, 6, 7, 8]), array([0, 2, 3, 6]), array([0, 3, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), 
  array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), 
  array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 5, 6]), array([1, 4, 5, 6])], 
[array([0, 1, 2, 5]), array([0, 1, 4, 5]), array([0, 2, 3, 5]), array([0, 3, 4, 5]), array([2, 3, 5, 8]), 
    array([2, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 8]), array([0, 6, 7, 8]), 
    array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), 
    array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 5, 6]), 
    array([1, 4, 5, 6])], [array([0, 1, 2, 3]), array([0, 1, 3, 4]), array([1, 2, 3, 4]), array([2, 3, 4, 5]), 
    array([2, 3, 5, 6]), array([3, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 6]), array([0, 3, 6, 8]), 
    array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), 
    array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), 
    array([1, 2, 4, 6]), array([2, 4, 5, 6])], [array([0, 1, 2, 4]), array([0, 2, 3, 5]), array([0, 2, 4, 5]), 
    array([0, 3, 4, 5]), array([2, 3, 5, 8]), array([2, 5, 6, 7]), array([2, 5, 7, 8]), array([0, 2, 3, 8]), 
    array([0, 2, 6, 7]), array([0, 2, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), 
    array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), 
    array([0, 1, 2, 6]), array([1, 2, 4, 6]), array([2, 4, 5, 6])], [array([0, 1, 2, 5]), array([0, 1, 4, 5]), 
    array([0, 2, 3, 5]), array([0, 3, 4, 5]), array([2, 3, 5, 6]), array([3, 5, 6, 8]), array([5, 6, 7, 8]), 
    array([0, 2, 3, 6]), array([0, 3, 6, 8]), array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), 
    array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), 
    array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 5, 6]), array([1, 4, 5, 6])], [array([0, 1, 2, 5]), 
    array([0, 1, 4, 5]), array([0, 2, 3, 5]), array([0, 3, 4, 5]), array([2, 3, 5, 8]), array([2, 5, 6, 7]), 
    array([2, 5, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 7]), array([0, 2, 7, 8]), array([0, 1, 4, 7]), 
    array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), array([4, 5, 6, 7]), 
    array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 5, 6]), array([1, 4, 5, 6])], 
[array([0, 1, 2, 4]), array([0, 2, 3, 5]), array([0, 2, 4, 5]), array([0, 3, 4, 5]), array([2, 3, 5, 8]), 
    array([2, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 8]), array([0, 2, 6, 8]), array([0, 6, 7, 8]), 
    array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), array([1, 4, 6, 7]), 
    array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), array([1, 2, 4, 6]), 
    array([2, 4, 5, 6])], [array([0, 1, 2, 4]), array([0, 2, 3, 5]), array([0, 2, 4, 5]), array([0, 3, 4, 5]), 
    array([2, 3, 5, 6]), array([3, 5, 6, 8]), array([5, 6, 7, 8]), array([0, 2, 3, 6]), array([0, 3, 6, 8]), 
    array([0, 6, 7, 8]), array([0, 1, 4, 7]), array([0, 1, 6, 7]), array([0, 3, 4, 8]), array([0, 4, 7, 8]), 
    array([1, 4, 6, 7]), array([4, 5, 6, 7]), array([3, 4, 5, 8]), array([4, 5, 7, 8]), array([0, 1, 2, 6]), 
    array([1, 2, 4, 6]), array([2, 4, 5, 6])]] "

--------------------------------------------------------
-- Code for complete intersections in toric varieties --
-- more functionality is in CohomCalg.m2              --
--------------------------------------------------------
CompleteIntersectionInToric = new Type of HashTable

completeIntersection = method()
completeIntersection(NormalToricVariety, List) := (Y,CIeqns) -> (
    if not all(CIeqns, d -> instance(d, ToricDivisor))
    then error "expected a list of toric divisors";
    if not all(CIeqns, d -> variety d === Y)
    then error "expected a list of toric divisors on the given toric variety";
    new CompleteIntersectionInToric from {
        symbol Ambient => Y,
        symbol CI => CIeqns,
        symbol cache => new CacheTable
        }
    )
dim CompleteIntersectionInToric := (X) -> dim X.Ambient - #X.CI
ambient CompleteIntersectionInToric := (X) -> X.Ambient

abstractVariety(CompleteIntersectionInToric, AbstractVariety) := opts -> (X,B) -> (
    if not X.cache#?(abstractVariety, B) then X.cache#(abstractVariety, B) = (
        aY := abstractVariety(ambient X, B);
        -- Question: how best to define F??
        bundles := X.CI/(d -> OO d);
        F := bundles#0;
        for i from 1 to #bundles-1 do F = F ++ bundles#i;
        aF := abstractSheaf(ambient X, B, F);
        sectionZeroLocus aF
        );
    X.cache#(abstractVariety, B)
    )

load (currentFileDirectory | "StringTorics/LineBundleCohomology.m2")

-----------------------
-- Examples -----------
-----------------------
example0912'3524 = () -> (value /// () -> (
    R := QQ[v1, v2, v3, v4, v5, v6, v1s, v7, v8, v9, v10];
    SR := {{v3,v9},{v5,v9},{v7,v10},{v1,v2,v3},
          {v4,v1s,v8},{v4,v7,v8},{v4,v8,v9},
          {v5,v6,v1s},{v5,v6,v10},{v1,v2,v6,v1s}};
    SR = monomialIdeal(SR/product);
    GLSM := transpose matrix {{3, 3, 3, 3, 0}, {2, 2, 2, 2, 0},
      {1, 0, 0, 0, 0}, {0, 0, 1, 0, 0}, {0, 0, 0, 1, 0},
      {0, 1, 0, 0, 0}, {0, 1, 1, 0, 0}, {0, 0, 1, 0, 1},
      {0, 0, 1, 0, 0}, {0,-1,-1, 1,-1}, {0, 0, 0, 0, 1}};
    normalToricVarietyFromGLSM(GLSM, SR)
    )
  ///
  )

exampleP11222'8 = () ->  (value /// () -> (
        -- P^4_{11222}[8]
        -- Actually: we desingularize the ZZ/2 - singularity of this
        rays := {{-1,-2,-2,-2}, 
                 {1,0,0,0},
                 {0,1,0,0},
                 {0,0,1,0},
                 {0,0,0,1},
                 {0,-1,-1,-1}};
        maxcones := {
            {0, 2, 3, 4}, 
            {0, 2, 3, 5}, 
            {0, 2, 4, 5}, 
            {0, 3, 4, 5}, 
            {1, 2, 3, 4}, 
            {1, 2, 3, 5}, 
            {1, 2, 4, 5}, 
            {1, 3, 4, 5}};
        GLSM := transpose matrix {
            {1,0}, {1,0}, {2,1}, {2,1}, {2,1}, {0,1}};
        normalToricVariety(rays, maxcones, WeilToClass => GLSM)
        )
    ///
    )

exampleP111122'44 = () -> (value /// () -> (
        rays := {{-1,-1,-1,-2,-2}, 
                 {1,0,0,0,0},
                 {0,1,0,0,0},
                 {0,0,1,0,0},
                 {0,0,0,1,0},
                 {0,0,0,0,1}};
        maxcones := {
            {0, 1, 2, 3, 4},
            {0, 1, 2, 3, 5},
            {0, 1, 2, 4, 5},
            {0, 1, 3, 4, 5},
            {0, 2, 3, 4, 5},
            {1, 2, 3, 4, 5}};
        GLSM := transpose matrix {
            {1}, {1}, {1}, {1}, {2}, {2}};
        normalToricVariety(rays, maxcones, WeilToClass => GLSM)
        )
    ///
    )

 findAllFRSTs = method()
 findAllFRSTs NormalToricVariety := (V) -> findAllFRSTs transpose matrix rays V
 findAllFRSTs Matrix := List => (A) -> (
     A1 := A | map(target A, (ring source A)^1, 0);
     Ts := allTriangulations(A1, Fine => true, RegularOnly => true);
     if #Ts === 0 or #Ts#0 == 0 then (
         count := 0;
         while count < 10 and (#Ts === 0 or #Ts#0 == 0) do (
             Ts = allTriangulations(A1, Fine => true, RegularOnly => true);
             count = count + 1;
             );
         --if #Ts == 0 then error "no triangulation could be found";
         << "WARNING: TOPCOM failed to find triangulations, then found them after " << 
           count << " attempt(s)" << endl;
         );
     --<< "Ts = (before selection): " << netList Ts << endl;
     Ts1 := select(Ts, isStar_A1);
     --<< "Ts1 = (after): " << netList Ts1 << endl;
     assert all(Ts1, tri -> all(tri, s -> s#-1 == numcols A));
     Ts1/(t -> (entries transpose A, t/(s -> drop(s, -1))))
     )
 findAllFRSTs Polyhedron := List => (P) -> (
     L := latticePointList P;
     assert all(L#-1, a -> a == 0);
     L = drop(L, -1);
     A := transpose matrix L;
     findAllFRSTs A
     )

-------------------------------------------------
-- Intersection numbers for Calabi-Yau 3-folds --
-------------------------------------------------
-- currently is functional only for hypersurfaces
-- in a toric 4-fold.  Also, the polytope must be favorable.
-- TODO: remove favorable hypothesis
-- TODO: allow 4D CY's?
-- TODO: allow CI CY's in a Fano toric.

load (currentFileDirectory | "StringTorics/IntersectionNumbers.m2")

------ end of Intersection Numbers code ---------


---------------------------------------------
-- ReflexivePolytope ------------------------
---------------------------------------------
reflexivePolytope = method()

-- The following is meant to be an internal method.
reflexivePolytope(List, List, List) := ReflexivePolytope => (latticePoints, GLSM, basisIndices) -> (
    -- What should be checked here to validate the input data?
    new ReflexivePolytope from {
        symbol cache => new CacheTable,
        "rays" => latticePoints,
        "glsm" => GLSM,
        "basis indices" => basisIndices
        }
    )

-- TODO: this function needs TLC.  It is useful, but can't always find the 
-- basis indices, or compute D_q^-1...
reflexivePolytope Polyhedron := ReflexivePolytope => (P2) -> (
    LP := select(latticePointList P2, lp -> dim(P2, minimalFace(P2, lp)) <= 2);
    mLP := transpose matrix LP;
    D := transpose syz mLP;
    p := findFirstUnitVectors D; -- TODO: p,q computation can be slow!
    q := findInvertibleSubmatrix(D, p);
    if q === null then error ("oops, can't find a good GLSM matrix"); -- hasn't happened yet. HAS NOW!!
    GLSM := (D_q)^-1 * D;
    -- the rays of each triangulation should match LP.
    result := reflexivePolytope(LP, entries transpose GLSM, q);
    result.cache#"N polytope" = P2;
    result.cache#"M polytope" = polar P2;
    result
    )

-- TODO: the following function often fails, due to TODO on above function.
reflexivePolytope Matrix := ReflexivePolytope => (A) -> (
    reflexivePolytope polar convexHull A
    )

degrees ReflexivePolytope := List => P -> P#"glsm"
rays ReflexivePolytope := List => P -> P#"rays" -- maybe call this something else?
dim ReflexivePolytope := ZZ => P -> dim polytope P
basisIndices = method()
basisIndices ReflexivePolytope := List => P -> P#"basis indices" -- returns the indices into `rays P` 

polar ReflexivePolytope := ReflexivePolytope => P -> reflexivePolytope polytope(P, "M")

polytope ReflexivePolytope := Polyhedron => P -> polytope(P, "N")
polytope(ReflexivePolytope, String) := Polyhedron => (P, which) -> (
    -- which is either "M" or "N"
    if which === "M" then P.cache#"M polytope"
    else if which === "N" then P.cache#"N polytope"
    else error "expected \"M\" or \"N\""
    )

annotatedFaces ReflexivePolytope := (cacheValue symbol annotatedFaces)(P -> (
    annotatedFaces polytope(P, "N")
    ))

annotatedFaces(ZZ, ReflexivePolytope) := (i,P) -> (
    -- i is the dimension of the face on polytope on the N side...
    A := annotatedFaces P;
    for f in A list if f#0 == i then drop(f, 1) else continue
    )

-- internal function to set h11, h21 in P.cache.
setH11H21 = P -> (
    -- this version is only for CY 3-fold hypersurfaces...
    -- P:ReflexivePolytope
    A := annotatedFaces P; -- polytope on N side.
    A0 := annotatedFaces(0, P);
    A1 := annotatedFaces(1, P);
    A2 := annotatedFaces(2, P);
    A3 := annotatedFaces(3, P);
    npM := A/(x -> x#4)//sum + 1;
    npN := A/(x -> x#3)//sum; -- origin is included in the dim 4 face.
    -- points in facets (on M side) -- this is part of h21
    -- points in facets (on N side) -- this is part of h11
    facetInteriorsM := A0/(v -> v#3)//sum;
    facetInteriorsN := A3/(v -> v#2)//sum;
    -- points interior to 2-faces (times their genus) (on M-side)
    -- points interior to 2-faces (times their genus) (on N-side)
    twoFacesM := A1/(v -> v#2 * v#3)//sum;
    twoFacesN := A2/(v -> v#2 * v#3)//sum;
    -- now set the h11, h21.
    h11 := npN - 5 - facetInteriorsN + twoFacesN;
    h21 := npM - 5 - facetInteriorsM + twoFacesM;
    P.cache#"h11" = h11;
    P.cache#"h21" = h21;
    P.cache#"favorable" = (twoFacesN == 0);
    (h11, h21)
    )

h11OfCY ReflexivePolytope := ZZ => P -> (
    if not P.cache#?"h11" then elapsedTime setH11H21 P;
    P.cache#"h11"
    )

h21OfCY ReflexivePolytope := ZZ => P -> (
    if not P.cache#?"h21" then elapsedTime setH11H21 P;
    P.cache#"h21"
    )

isFavorable ReflexivePolytope := Boolean => P -> (
    if not P.cache#?"favorable" then elapsedTime setH11H21 P;
    P.cache#"favorable"
    )

-- TODO: test that the lattice points, vertices match up...
-- i.e. annotatedFaces align with other aspects(?) of these polyhedra.
TEST ///
  restart
  needsPackage "StringTorics"
  
  topes = kreuzerSkarke(3, Limit => 50);    
  A = matrix topes_30
  convexHull A
  elapsedTime h11OfCY oo
  polar convexHull A
  elapsedTime P = reflexivePolytope A
  elapsedTime assert(h11OfCY P == 3)
  elapsedTime assert(h21OfCY P == 69)
  elapsedTime assert(isFavorable P)
  netList elapsedTime annotatedFaces P
  annotatedFaces(4, P)
  annotatedFaces(3, P)
  annotatedFaces(2, P)
  annotatedFaces(1, P)
  annotatedFaces(0, P)
  
  
  annotatedFaces(0, polytope P)
  elapsedTime latticePointList polytope P
  assert(dim P == dim polytope P)
  rays P
  latticePoints polytope P
  elapsedTime transpose matrix latticePointList polytope P -- these should match annotatedFaces.  TEST THIS.
  transpose matrix rays P -- this matches the non-zero lattice points of P.  Or maybe those not interior to 3-faces.
  annotatedFaces polytope P
  elapsedTime assert(h11OfCY P == 3)
  elapsedTime assert(h21OfCY P == 69)
  elapsedTime assert isFavorable P
  -- Q = polar P  -- doesn't work yet for this example.
///
----------------------------------------------------------------


----------------------------------------------------------------
-- FRST Triangulations (Fine, regular, star triangulations) ----
----------------------------------------------------------------
-- 
  -- Here, we only consider a triangulation of a reflexive polytope
  
CalabiYauInToric.synonym = "normal toric variety"
CalabiYauInToric.GlobalAssignHook = globalAssignFunction
CalabiYauInToric.GlobalReleaseHook = globalReleaseFunction
expression CalabiYauInToric := X -> if hasAttribute (X, ReverseDictionary) 
    then expression getAttribute (X, ReverseDictionary) else 
    (describe X)#0
describe CalabiYauInToric := X -> Describe (expression CalabiYauInToric) (
      expression "a" , expression 3)
--    expression rays X, expression max X)
  
makeCYInToric = (reflexivePolytope, FRSTtriangulation) -> (
    -- TODO: consistency check for data.
    new CalabiYauInToric from {
        symbol cache => new CacheTable,
        "polytope" => reflexivePolytope,
        "triangulation" => FRSTtriangulation
        }
    )

net CalabiYauInToric := X -> "a Calabi-Yau hypersurface in a simplicial toric variety of dimension " | dim polytope X

-- net CalabiYauInToric := T -> net T#"triangulation"
-- CalabiYauInToric#{Standard,AfterPrint} = X -> (
--      << endl;				  -- double space
--      << concatenate(interpreterDepth:"o") << lineNumber << " : "
--      << "a Calabi-Yau hypersurface in a simplicial toric variety of dimension " << 4 << endl;
--      )

findAllFRSTs ReflexivePolytope := List => P -> (
    T := findAllFRSTs transpose matrix rays P;
    for t in T do if last t === {} then error "what?!"; -- this should not happen?
    for t in T list makeCYInToric(P, last t) -- t is a pair: list of vertices, list of list of indices
    )

findOneCY = method()
findOneCY ReflexivePolytope := CalabiYauInToric => P -> (
    P2 := polytope P;
    (LP,tri) := regularStarTriangulation(dim P2-2,P2);
    if rays P =!= LP then error "I have a lattice point mismatch";
    makeCYInToric(P, tri)
    )
    
reflexiveToSimplicialToricVariety Polyhedron := opts -> (P1) -> (
    -- P1 is a reflexive polytope in the M lattice.
    -- Creates a simplicial toric variety via a triangulation
    -- of the normal fan of P1, using all of lattice points of (polar P1).
    -- Returns: a normal toric variety, whch uses all
    -- lattice points in the triangulation.
    P2 := polar P1;
    (LP,tri) := regularStarTriangulation(dim P2-2,P2);
    normalToricVariety(LP,tri,opts)
    )

normalToricVariety CalabiYauInToric := opts -> X -> (
    if not X.cache.?NormalToricVariety then X.cache.NormalToricVariety = (
        P := X#"polytope";
        T := X#"triangulation";
        GLSM := transpose matrix P#"glsm";
        normalToricVariety(rays P, T, opts, WeilToClass => matrix GLSM)
        );
    X.cache.NormalToricVariety
    -- TODO: this fails if the class group is torsion! (Fails: later it gives an inscrutable error...)
    )

ambient CalabiYauInToric := X -> normalToricVariety X
polytope CalabiYauInToric := X -> X#"polytope"
basisIndices CalabiYauInToric := X -> basisIndices polytope X
abstractVariety CalabiYauInToric := opts -> X -> (
    -- Store this with X.
    V := ambient X;
    aX := completeIntersection(V, {-toricDivisor V});
    abstractVariety(aX, base())
    )
abstractVariety(CalabiYauInToric, AbstractVariety) := opts -> (X, pt) -> (
    -- Store this with X?
    -- Check: pt is of dimension zero?
    V := ambient X;
    aX := completeIntersection(V, {-toricDivisor V});
    abstractVariety(aX, pt)
    )

topologicalData = method()
topologicalData(CalabiYauInToric, Ring) := TopologicalDataOfCY3 => (X, RZ) -> (
    V := ambient X;
    P := polytope X;
    data := elapsedTime topologyOfCY3(V, P#"basis indices");
    -- this data above computes intersection numbers for all toric divisors. 
    -- So we consider only the ones whose indices are contained in basis indices:
    new TopologicalDataOfCY3 from {
        "h11" => elapsedTime h11OfCY P,
        "h21" => elapsedTime h21OfCY P,
        "c2" => sub(data_3, vars RZ),
        "cubic intersection form" => sub(data_2, vars RZ)
        }
    )

c2 = method()
cubicForm = method()

hh(Sequence, TopologicalDataOfCY3) := (pq, T) -> (
    (p,q) := pq;
    if p > q then (p, q) = (q, p);
    if p == 0 then (
        if q == 3 or q == 0 then 1 else 0
        )
    else if p == 1 then (
        if q == 1 then T#"h11"
        else if q == 2 then T#"h21"
        else 0
        )
    else if p == 2 then (
        if q == 2 then T#"h11" else 0
        )
    else if p == 3 then (
        if q == 3 then 1
        else 0
        )
    )

--h11 CalabiYauInToric := X -> h11OfCY polytope X
--h21 CalabiYauInToric := X -> h21OfCY polytope X
c2 TopologicalDataOfCY3 := T -> T#"c2"
cubicForm TopologicalDataOfCY3 := T -> T#"cubic intersection form"

TEST ///
-- XXX
  restart
  needsPackage "StringTorics"
  
  topes = kreuzerSkarke(3, Limit => 50);    
  A = matrix topes_30
  P = reflexivePolytope A
  findAllFRSTs P
  X = first oo
  
  h11OfCY polytope P
  elapsedTime polar polytope P  
  h21OfCY polytope P  
  RZ = ZZ[x,y,z]
  elapsedTime topologicalData(X, RZ) -- cache this result?
  
  dim X
  ambient X -- give the normal toric variety.  Works now.
  abstractVariety X -- give the abstract variety
  abstractVariety(X, base(a,b,c)) -- give the abstract variety
  inheritedMoriCone X
  gvInvariants X

  aX = abstractVariety X -- TODO: should stash the value...
  IX = intersectionRing aX
  
  debug StringTorics
  intersectionNumbers X  

  
  -- TODO: line bundles on X, and their cohomology.
  
  topes = kreuzerSkarke(6, Limit => 50);    
  A = matrix topes_30
  P = reflexivePolytope A
  assert(# findAllFRSTs P == 8)

  A = matrix topes_40
  P = reflexivePolytope A
  elapsedTime findAllFRSTs P -- takes longer than I would like...
  assert(#oo == 36) -- not sure if this is correct, but checking for possible change.

  topes = kreuzerSkarke(8, Limit => 50);    
  A = matrix topes_30
  P = reflexivePolytope A
  --elapsedTime findAllFRSTs P; -- takes some time...  How many are there?  Improve this time...

///

----------------------------------------------------------------  
-- gvInvariants ------------------------------------------------
-- Uses computeGV.cpp from CYtools -----------------------------
----------------------------------------------------------------
gvInvariants = method(Options => {
    Mori => null, -- null means: compute rays of the Mori cone of V (in ZZ^(h11))
    Heft => null, -- null means: compute it
    DegreeLimit => infinity,
    Precision => 150,
    FilePrefix => "foo",
--    Executable => "~/src/git-from-others/cytools-private/external/gv/computeGV-good/computeGV",
    Executable => "~/src/git-from-others/cytools-private/external/gv/computeGV",
    KeepFiles => true
    })

intersectionNumbersOfCY = method()
intersectionNumbersOfCY(NormalToricVariety, List) := (V, basisIndices) -> (
    X := completeIntersection(V, {-toricDivisor V});
    Xa := abstractVariety(X, base());
    IX := intersectionRing Xa;
    intersectionNumbers(IX, basisIndices)
    )

-- The function to write the data needed by the computeGV program
gvInput = (moriGenerators, heftval, GLSM, intersectionnums, degreelimit, prec) -> (
    -- moriGenerators: list of lists
    -- heftval: list of ints
    -- GLSM: list of list of ints
    -- intersectionnums: list of triples of ints
    -- degreelimit: infinity or positive integer
    -- prec: positive integer
    str1 := toString moriGenerators;
    str3 := toString heftval;
    str4 := toString GLSM;
    str5 := toString intersectionnums;
    str6 := toString ({
            if degreelimit === infinity then -1 else degreelimit, 
            prec,
            0,
            300000
            });
    concatenate between("\n", {str1, toString {}, str3, str4, toString {}, str5, str6})
    )

moriCone = method()
moriCone NormalToricVariety := List => (V) -> (
    IV := intersectionRing (abstractVariety V);
    Cs := matrix for x in orbits(V, 1) list (
        c := product(x, i -> IV_i);
        for d in gens IV list integral(c*d)
        );
    M := posHull transpose lift(Cs, QQ);
    GLSM := matrix degrees ring V;
    entries transpose((rays M) // GLSM)
    )

toricMoriCone = method()
toricMoriCone CalabiYauInToric := Cone => X -> (
    posHull transpose matrix moriCone X
    )

gvInvariants(NormalToricVariety, List) := HashTable => opts -> (V, basisIndices) -> (
    -- Compute intersection numbers for X in V (using this basis)
    -- Compute mori cone (if needed) (?? requires basis too...)
    -- Compute a vector which dots positively with all these generators.
    -- Then write the file
    -- Execute the command
    -- Read the results, and return them
    intersectionnums := for t in intersectionNumbersOfCY(V, basisIndices) list append(t#0, t#1);
    -- X := completeIntersection(V, {-toricDivisor V});
    -- Xa := abstractVariety(X, base());
    -- IX := intersectionRing Xa;
    -- intersectionnums := for x in pairs intersectionNumbers(IX, basisIndices) list append(x#0, x#1);
    -- H := hashTable for i from 0 to #basisIndices-1 list basisIndices#i => i;
    -- intersectionnums := for x in pairs CY3NonzeroMultiplicities V list (
    --     if isSubset(x#0, basisIndices) then
    --         append(sort for a in x#0 list H#a, x#1)
    --     else 
    --         continue
    --     );
    mori := if opts.Mori =!= null then opts.Mori else moriCone V;
    heft := if opts.Heft =!= null then opts.Heft else (
      sum entries transpose rays dualCone posHull transpose matrix mori
      );
    -- OK, now we have computed everything we need.  Write it to a file
    infile := opts.FilePrefix | "-input";
    outfile := opts.FilePrefix | "-output";
    infile << gvInput(mori, heft, transpose degrees ring V, intersectionnums,
        opts.DegreeLimit, opts.Precision) << close;
    inputLine := opts.Executable | " <" | infile | " >" | outfile;
    print inputLine;
    run inputLine;
    -- Get the output, package as a hash table
    (lines get outfile)/value//hashTable
    )

gvInvariants CalabiYauInToric := HashTable => opts -> X -> (
    gvInvariants(ambient X, basisIndices polytope X, opts)
    )

TEST ///
-- XXX
  restart
  needsPackage "StringTorics"
  
  topes = kreuzerSkarke(3, Limit => 50);    
  A = matrix topes_30
  A = matrix topes_31
  P = reflexivePolytope A
  findAllFRSTs P
  X = first oo
  polytope X

  V = ambient X
  intersectionNumbersOfCY(V, basisIndices polytope X)
  elapsedTime T = topologicalData(X, ZZ[a,b,c])
  hh^(1,1) T
  hh^(1,2) T
  hh^(1,1) X
  hh^(1,2) X
  hodgeDiamond X

  gv = gvInvariants(X, DegreeLimit => 30);
  C = posHull transpose matrix ((keys gv)/toList)
  rays C
  moriX = entries transpose rays C -- not exactly moriX..
  partition(f -> gv#(toSequence f), moriX)

  moriConeByGV(X, DegreeLimit => 30)
  moriConeByGV(X, DegreeLimit => 20)
  moriConeByGV(X, DegreeLimit => 10)
  moriConeByGV(X, DegreeLimit => 5)
///  
----------------------------------------------------------------

factors = method()
factors RingElement := (F) -> (
     facs := factor F;
     facs//toList/toList/reverse
     )

invariants = method()
invariants List := (f) -> (
    RQ := QQ[gens ring first f];
    facs := select((factors f_1 )/toList/last, g -> support g != {});
    l1 := sub(f_0, RQ);
    f1 := sub(f_1, RQ);
    d := dim saturate ideal jacobian f1;
    nc := # decompose ideal(l1, f1);
    singZ := flatten entries gens gb saturate(ideal(f_1) + ideal jacobian f_1);
    badp := select(singZ, a -> support leadTerm a === {});
    badp = if badp === {} then 0 else first badp;
    {badp, (trim content f_0)_0, (trim content f_1)_0, #facs, d, nc, f_2, f_3}
    )

moriConeByGV = method(Options => options gvInvariants)
moriConeByGV CalabiYauInToric := HashTable => opts -> X -> (
    gv := gvInvariants(X, opts);
    C := posHull transpose matrix ((keys gv)/toList);
    moriX := entries transpose rays C; -- not exactly moriX..
    partition(f -> gv#(toSequence f), moriX)
    )

----------------------------------------------------------------

beginDocumentation()

-- Problems
-- . how to determine favorable?
-- . which toric variety do we want?
-- . want a smooth one
-- . want to compute h^11(X) using cohomcalg, but
--   examples can be too big
-- . triangulations can be too big
-- . what else can be too big?


load (currentFileDirectory | "StringTorics/doc.m2")
load (currentFileDirectory | "StringTorics/test.m2")



end--

restart
needsPackage "StringTorics"

restart
uninstallPackage "StringTorics"
restart
installPackage "StringTorics"

restart
needsPackage "StringTorics"
check oo

-- Generation of some examples
L = kreuzerSkarke(20, Limit => 100, Access=>"wget")
L = kreuzerSkarke(15, Limit => 100)
L = kreuzerSkarke(13, Limit => 100)
L = kreuzerSkarke(10, Limit => 100)
A = matrix L_50
P = convexHull A
isReflexive P
X = reflexiveToSimplicialToricVariety P
  max X
  # rays X
  sr = dual monomialIdeal X
  for i from 0 to #L-1 list (
      elapsedTime X = reflexiveToSimplicialToricVariety convexHull matrix L_i;
      if # rays X > 64 then continue;
      elapsedTime S = try ring X else null;
      if S === null then continue;
      sr = elapsedTime dual monomialIdeal X;
      if numgens sr > 64 then continue;
      i
      )
  Xs = for i from 0 to #L-1 list (
      elapsedTime X = reflexiveToSimplicialToricVariety convexHull matrix L_i
      )

  tally apply(Xs, X -> # rays X)
  select(Xs, X -> elapsedTime try (ring X; true) else false) -- h11=13, h11=15: all ok here...
  tally apply(Xs, X -> numgens elapsedTime dual monomialIdeal X)
  tally apply(Xs, X -> # rays X)
  positions(Xs, X -> numgens elapsedTime dual monomialIdeal X < 64)
  max Xs_28
  debug CohomCalg
  print toCohomCalg Xs_28  -- h11 = 15
  
  -- h11 = 15:
  positions(Xs, X -> numgens elapsedTime dual monomialIdeal X == 39) -- h11=15 indices: 51, 66.
  max Xs_51      
  rays Xs_51
  print toCohomCalg Xs_51

  -- h11 = 15:
  positions(Xs, X -> numgens elapsedTime dual monomialIdeal X == 24) -- h11=10 indices: 6.
  max Xs_6
  rays Xs_6
  print toCohomCalg Xs_6
    

str = getKreuzerSkarke(4,100,Limit=>10)
time polytopes = parseKS str;
#polytopes
time for p in polytopes list (p_0, matrixFromString p_1);


str = getKreuzerSkarke(50,100,Limit=>4000);
str = getKreuzerSkarke(50,100,Limit=>0)
str = getKreuzerSkarke(4,4,Limit=>10)
str = getKreuzerSkarke(4,8,Limit=>10)
for i from 1 to 20 list parseKS getKreuzerSkarke(4,i,Limit=>10)
for i from 21 to 30 list parseKS getKreuzerSkarke(4,i,Limit=>10)
for i from 31 to 40 list parseKS getKreuzerSkarke(4,i,Limit=>10)
join oo
flatten oo
netList oo

-- Experiments for Batyrev formula
restart
needsPackage "StringTorics"
getKreuzerSkarke(3,57)
polytopes = parseKS oo;
netList polytopes
polystr = polytopes_2
A = matrixFromString polystr_1
P = convexHull A
P2 = polar P

-- Part 1:
# latticePoints(P2) == 8

-- Part2: facets of P2
faces1 = faces(1, P2)

F = faces1_0
vertices F
latticePoints F

for f in faces1 list ((# latticePoints f) - numColumns vertices F)
F = faces1_3
hyperplanes F
(first (hyperplanes F)) * vertices F
H = for f in faces(1,F) list hyperplanes f
latticePoints F
for p in oo list H_0_0 * p

-- Need a function: given a lattice point, and a polytope, is the point in the polytope?
  -- Is it on the boundary?
L = latticePoints F  
interiorLatticePoints F
for p in L list contains(F,p)
for f in faces(1,P2) list interiorLatticePoints f
for f in faces(2,P2) list interiorLatticePoints f
F22 = faces(2,P2)
first oo
polar oo
vertices oo

F2 = faces(2,P)
#F22
#F2
faces(4,P)
faces(1,P2)
# faces(2,P2)
# faces(3,P)
faces(3,P2)
faces(2,P)

-- Ben's code for morigenerators, translated into M2.
moriGenerators = method()
moriGenerators(List, List) := (rays, cones) -> (
    -- rays: a list of points in ZZ^n
    -- cones: a list of lists of size n, each element is a subset of 0..n-1,
    --   representing the simplicial cones of a fan with rays 'rays'
    n := #(rays#0);
    cones := cones/set;
    -- assert: all rays have length n, and all cones have n as well
    walls := select(subsets(n,2), x -> # toList (cones#(x_0) * cones#(x_1)) == n-1);
    commons := for x in walls list toList (cones#(x_0) + cones#(x_1));
    walls
    )
rays25 = {{2, -1, -1, -1}, {-1, 0, -1, -1}, {-1, 1, -1, -1}, {-1, 0, 
    4, -1}, {-1, 1, 0, -1}, {-1, 0, -1, 0}, {-1, 1, -1, 4}, {-1, 0, 4,
     0}, {-1, 1, 0, 4}, {-1, 1, 0, 3}, {-1, 1, -1, 3}, {-1, 1, 0, 
    2}, {-1, 1, -1, 2}, {-1, 1, 0, 1}, {-1, 1, -1, 1}, {-1, 1, 0, 
    0}, {-1, 1, -1, 0}, {0, 0, 1, 1}, {0, 0, 0, 1}, {0, 0, -1, 1}, {0,
     0, 1, 0}, {0, 0, -1, 0}, {0, 0, 1, -1}, {0, 0, 0, -1}, {0, 
    0, -1, -1}, {-1, 0, 3, 0}, {-1, 0, 2, 0}, {-1, 0, 1, 0}, {-1, 0, 
    0, 0}, {-1, 0, 3, -1}, {-1, 0, 2, -1}, {-1, 0, 1, -1}, {-1, 0, 
    0, -1}};
cones25 = {{0, 7, 9, 11}, {0, 6, 9, 11}, {6, 7, 9, 11}, {0, 4, 10, 
    12}, {0, 5, 10, 12}, {4, 5, 10, 12}, {0, 4, 7, 15}, {0, 4, 6, 
    15}, {4, 6, 7, 15}, {0, 2, 4, 16}, {0, 2, 5, 16}, {2, 4, 5, 
    16}, {6, 7, 8, 17}, {7, 8, 9, 17}, {6, 7, 8, 9}, {0, 6, 8, 
    17}, {0, 8, 9, 17}, {0, 6, 8, 9}, {0, 6, 7, 17}, {0, 7, 9, 
    17}, {0, 6, 7, 18}, {3, 4, 5, 6}, {5, 6, 10, 19}, {4, 5, 6, 
    10}, {5, 6, 18, 19}, {0, 4, 6, 10}, {0, 6, 10, 19}, {0, 6, 18, 
    19}, {0, 5, 10, 19}, {0, 5, 18, 19}, {7, 11, 13, 20}, {6, 7, 11, 
    13}, {7, 13, 15, 20}, {6, 7, 13, 15}, {0, 11, 13, 20}, {0, 6, 11, 
    13}, {0, 13, 15, 20}, {0, 6, 13, 15}, {0, 7, 11, 20}, {0, 7, 15, 
    20}, {5, 12, 14, 21}, {4, 5, 12, 14}, {5, 14, 16, 21}, {4, 5, 14, 
    16}, {0, 4, 12, 14}, {0, 12, 14, 21}, {0, 4, 14, 16}, {0, 14, 16, 
    21}, {0, 5, 12, 21}, {0, 5, 16, 21}, {2, 3, 4, 22}, {2, 3, 4, 
    5}, {3, 4, 7, 22}, {3, 4, 6, 7}, {0, 2, 4, 22}, {0, 4, 7, 22}, {0,
     2, 3, 22}, {0, 3, 7, 22}, {0, 2, 3, 23}, {1, 2, 5, 24}, {1, 2, 
    23, 24}, {0, 2, 5, 24}, {0, 2, 23, 24}, {0, 1, 5, 24}, {0, 1, 23, 
    24}, {6, 7, 18, 25}, {3, 6, 7, 25}, {0, 7, 18, 25}, {0, 3, 7, 
    25}, {6, 18, 25, 26}, {3, 6, 25, 26}, {0, 18, 25, 26}, {0, 3, 25, 
    26}, {6, 18, 26, 27}, {3, 6, 26, 27}, {0, 18, 26, 27}, {0, 3, 26, 
    27}, {6, 18, 27, 28}, {3, 6, 27, 28}, {0, 18, 27, 28}, {0, 3, 27, 
    28}, {5, 6, 18, 28}, {3, 5, 6, 28}, {0, 5, 18, 28}, {0, 3, 5, 
    28}, {2, 3, 23, 29}, {2, 3, 5, 29}, {0, 3, 5, 29}, {0, 3, 23, 
    29}, {2, 23, 29, 30}, {2, 5, 29, 30}, {0, 5, 29, 30}, {0, 23, 29, 
    30}, {2, 23, 30, 31}, {2, 5, 30, 31}, {0, 5, 30, 31}, {0, 23, 30, 
    31}, {2, 23, 31, 32}, {2, 5, 31, 32}, {0, 5, 31, 32}, {0, 23, 31, 
    32}, {1, 2, 23, 32}, {1, 2, 5, 32}, {0, 1, 5, 32}, {0, 1, 23, 
    32}};

-- Investigate cohomology cones
-- fano(4,10):
B = ideal V
S = (coefficientRing ring B)[gens ring B, DegreeRank=>numgens ring B]
B = sub(B,S)
Ext^2(comodule B, S)
  -- get x_0*x_1*x_6, x_0*x_5*x_6
  -- remember: allow these to be as negative as desired, rest of the variables must be positive.
  -- Still: get a cone in RR^n (n is 7 here).
  degs = degrees ring V
  for i from 0 to numgens ring V list if (x_0*x_1*x_6)%x_i == 0 then - degs_i else degs_i
  

---------------------------------------------------
-- remove code below this line (18 May 2019) !!! --
---------------------------------------------------

-- This changes the function from ReflexivePolytopesDB
parseKS String := (str) -> (
    -- result is a List of pairs of strings.
    locA := regex("<b>Result:</b>\n", str);
    locB := regex("#NF", str);
    if locB === null then locB = regex("Exceeded", str);
    --if locA === null or locB === null then error "data not in correct Kreuzer-Skarke format";
    firstloc := if locA === null then 0 else locA#0#0 + locA#0#1;
    lastloc := if locB === null then #str else locB#0#0;
    --firstloc := locA#0#0 + locA#0#1;
    --lastloc := locB#0#0;
    cys := substring(firstloc, lastloc-firstloc, str);
    cys = lines cys;
    cys = select(cys, s -> #s > 0);
    starts := positions(cys, s -> s#0 != " ");
    starts = append(starts, #cys);
    for i from 0 to #starts-2 list (
        cys_(starts#i) | " id:" | i | "\n" |  demark("\n", cys_{starts#i+1 .. starts#(i+1)-1})
        )
    )

-- expect that str is something like:
str = ///4 5  M:53 5 N:9 5 H:3,43 [-80] id:0
   1   0   2   4 -10
   0   1   3   5  -9
   0   0   4   0  -4
   0   0   0   8  -8
   ///

matrixFromKS = method()
matrixFromKS String := (str) -> (
    matrixFromString concatenate between("\n", drop(lines str, 1))
    )

matrixFromKSEntry = method()
matrixFromKSEntry String := (str) -> (
    matrixFromString concatenate between("\n", drop(lines str, 1))
    )

///
  assert(matrixFromKS str == matrix {{1, 0, 2, 4, -10}, {0, 1, 3, 5, -9}, {0, 0, 4, 0, -4}, {0, 0, 0, 8, -8}})
///
  
  

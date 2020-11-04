-- TODO (Nov 2019)
-- 1. documentation
-- 2. Extra polyhedral functions (do we still need these?  Probably...)
--    Make sure they work with non reflexive polytopes too!
--    Or say they onyl work for reflexives?
-- 3. Triangulations
--    Computing FRST's
-- 4. Cohomology of line bundles on a toric variety
-- 5. Complete Intersections in torics
-- 6. Intersection rings
-- 7. Mori cones

newPackage(
        "StringTorics",
        Version => "0.1", 
        Date => "",
        Authors => {{Name => "Mike Stillman", 
                  Email => "", 
                  HomePage => ""}},
        Headline => "toric variety functions for string theory",
        DebuggingMode => true,
        AuxiliaryFiles => true,
        PackageExports => {
            "FourTiTwo",
            "SimplicialComplexes",
            "NormalToricVarieties",
            "Schubert2",
            "Polyhedra",
            "ReflexivePolytopesDB",
            "CohomCalg",
            "Topcom"
--            "AbstractToricVarieties"
            },
        PackageImports => {"Graphs"}
        )

export {
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
    
    -- This set maybe should be included in NormalToricVarieties?
    "singularCones",
    "singularLocusInToric",
    "normalToricVarietyFromGLSM",
    -- These should stay here
    "reflexiveToSimplicialToricVariety",
    "allZeros",
    "augment",

    "readSageTriangulations",
    "sortTriangulation",
    "matchNonZero",
    "applyPermutation",
    "findMaxRegularStarTriangulation",
    "affineCircuits",
    "checkFan",
    "sageTri",
    "Regular",

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
    "tentativeHodgeVector" -- deprecated
    }



load (currentFileDirectory | "StringTorics/MyPolyhedra.m2")
load (currentFileDirectory | "StringTorics/ToricCompleteIntersections.m2")
load (currentFileDirectory | "StringTorics/triangulations-code.m2")
--<< "-------------------------------------------------------------------------" << endl;
--<< "-- WARNING: StringTorics is still experimental, and not complete.  The --" << endl;
--<< "--   interface might change, and it is not well documented yet.        --" << endl;
--<< "-------------------------------------------------------------------------" << endl;


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

augment = (A) -> (
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

----------------------------------------------------------------
-- Cohomology of line bundles: obtaining a basis of fractions --
----------------------------------------------------------------
toricCohomologySetup = method()
toricCohomologySetup NormalToricVariety := (X) -> (
    if not X.cache#?CohomologySetup then (
        -- part 1: compute free res of B^* (over ZZ^n)
        R := ring X;
        R1 := (coefficientRing R)(monoid[gens R, DegreeRank=>numgens R]);
        B1 := monomialIdeal sub(ideal X, vars R1);
        B1' := dual B1;
        C := resolution comodule B1';
        -- part 2: take the (fine) degrees in there, and place them depending on HH^i, in a hash table
        -- part 3: create multi-graded rings for each of these
        -- finally, stash this info into the cache table of X.
        -- outside of this: make a function that takes a degree (or an OO(D)?), and i, and computes the fractions for HH^i(V, OO(V)).
        rawDegs := sort flatten for j from 0 to length C list (degrees C_j)/(x -> (sum x - j,x));
        rawDegsWithRings := for x in rawDegs list (
            -- x is (i, I), where I has length n is an array of 0 and 1's.
            Ra := if all(x#1, i -> i == 0) then R else (
                ds := degrees R;
                degs := for i from 0 to #x#1-1 list if x#1#i == 0 then ds_i else -ds_i;
                (coefficientRing R)(monoid[generators R, Degrees=>degs])
                );
            (x#0, x#1, Ra)
            );
        H := partition(x -> x#0, rawDegsWithRings);
        X.cache#CohomologySetup = applyPairs(H, x -> (x#0, (x#1/(y -> drop(y,1)))));
        );
    X.cache#CohomologySetup
    )
cohomologyBasis = method()
cohomologyBasis(ZZ, NormalToricVariety, List) := (i,X,L) -> (
    H := toricCohomologySetup X;
    L = -L;
    KR := frac ring X;
    if not H#?i then return {};
    degs := H#i;
    flatten for x in degs list (
        -- x is (I, R), I is a vector of 0,1's of length n = numgens R, and R is a multigraded ring
        (I,Ra) := x;
        deg := degree(Ra_I);
        ---M := transpose matrix degrees Ra | transpose matrix{-L+deg};
        ---result := (normaliz(M, 5))#"gen";
        ---elems0 := if result === null then {} else entries result;
        ---elems1 := select(elems0, v -> v#-1 == 1);
        ---elems2 := elems1/(v -> drop(v,-1));
        ---exps := elems2;
        --<< "L=" << L << " deg=" << deg << " -L-deg=" << -L-deg << endl;
        elems := flatten entries basis(-L-deg, Ra);
        --print elems;
        exps := elems/exponents/first;
        -- now need to place these into frac R.
        for e in exps list (
            e1 := for i from 0 to #I-1 list if I#i == 0 then e#i else - e#i - 1;
            KR_e1
            )
        )
    )
genericCohomologyMatrix = method()

-- Compute the map HH^i(V, OO_V(D+K)) --> HH^i(V, OO_V(D))
-- as a mutable matrix, over a finite field.
-- given V, i, and D.
-- the map is given by multiplication by a random form F of HH^0(OO_V(K)).
-- Caveat: it is extemely unlikely, but could happen that the generic rank
--   does not occur: the actual generic rank could potentially be larger than 
--   the rank of this matrix.
-- 
genericCohomologyMatrix(ZZ, NormalToricVariety, List) := (i,V,D) -> (
    K := degree toricDivisor V;
    basis1 := cohomologyBasis(1, V, D+K);
    basis2 := cohomologyBasis(1, V, D);
    basis0 := first entries basis(-K, ring V);
    hash1 := hashTable for i from 0 to #basis1-1 list basis1#i => i;
    hash2 := hashTable for i from 0 to #basis2-1 list basis2#i => i;
    M := mutableMatrix(ZZ/32003, #basis2, #basis1);
    for p in basis0 do (
        c := random ring M;
        for k in keys hash1 do (
            m := p * k;
            if hash2#?m then M_(hash2#m, hash1#k) = c;
            )
        );
    M
    )

cohomologyMatrix = method()
cohomologyMatrix(ZZ, NormalToricVariety, List, RingElement) := (i,V,D,F) -> (
    -- Create the map HH^i(V, D - degree F) --> HH^i(V, D), induced by mult by F (in Cox ring).
    kk := coefficientRing ring F;
    deg := degree F; -- F should be an element of the Cox ring.
    basis1 := cohomologyBasis(i, V, D-deg);
    basis2 := cohomologyBasis(i, V, D);
    basis0 := terms F;
    hash1 := hashTable for i from 0 to #basis1-1 list basis1#i => i;
    hash2 := hashTable for i from 0 to #basis2-1 list basis2#i => i;
    M := mutableMatrix(kk, #basis2, #basis1);
    for tm in basis0 do (
        c := leadCoefficient tm;
        p := leadMonomial tm;
        for k in keys hash1 do (
            m := p * k;
            if hash2#?m then M_(hash2#m, hash1#k) = c;
            )
        );
    (M, basis2, basis1)
    )

cohomologyMatrixRank = method()
cohomologyMatrixRank(ZZ, NormalToricVariety, List, RingElement) := (i,V,deg,F) -> (
    z := cohomologyMatrix(i,V,deg,F);
    nrows := #z#1;
    ncols := #z#2;
    rk := rank z#0;
    (nrows, ncols, rk)
    )

-- The above code doesn't work on larger examples.  Let's see if we can do better.
toricOrthants = method()
toricOrthants NormalToricVariety := (X) -> (
        -- part 1: compute free res of B^* (over ZZ^n)
        R := ring X;
        R1 := (coefficientRing R)(monoid[gens R, DegreeRank=>numgens R]);
        B1 := monomialIdeal sub(ideal X, vars R1);
        B1' := dual B1;
        C := resolution comodule B1';
        -- part 2: take the (fine) degrees in there, and place them depending on HH^i, in a hash table
        -- part 3: create multi-graded rings for each of these
        -- finally, stash this info into the cache table of X.
        -- outside of this: make a function that takes a degree (or an OO(D)?), and i, and computes the fractions for HH^i(V, OO(V)).
        rawDegs := sort flatten for j from 0 to length C list (degrees C_j)/(x -> (sum x - j,x));
        rawDegs/(x -> (x#0, positions(x#1, i -> i != 0)))
    )

H1orthants = method()
H1orthants NormalToricVariety := (X) -> (
        R := ring X;
        R1 := (coefficientRing R)(monoid[gens R, DegreeRank=>numgens R]);
        B1 := monomialIdeal sub(ideal X, vars R1);
        B1' := dual B1;
        R2 := (coefficientRing R)(monoid[gens R]);
        B1'' := sub(B1', R2);
        I := ideal select(flatten entries gens B1'', m -> degree m === {2});
        C := resolution(I, DegreeLimit=>1);
        C1 := for i from 1 to length C list sub(C.dd_i, R1);
        C2 := new MutableList;
        C2#0 = map(target C1_0,, C1_0);
        for i from 1 to length C1-1 do 
          C2#i = map(source C2#(i-1),, C1#i);
        degs := flatten for phi in toList C2 list degrees source phi;
        return for d in degs list positions(d, a -> a != 0);
        -- part 2: take the (fine) degrees in there, and place them depending on HH^i, in a hash table
        -- part 3: create multi-graded rings for each of these
        -- finally, stash this info into the cache table of X.
        -- outside of this: make a function that takes a degree (or an OO(D)?), and i, and computes the fractions for HH^i(V, OO(V)).
        rawDegs := sort flatten for j from 0 to length C list (degrees C_j)/(x -> (sum x - j,x));
        rawDegs/(x -> (x#0, positions(x#1, i -> i != 0)))
    )

-- goal: given V, find a permutation P, a d x d invertible (over ZZ) matrix Q, s.t. QAP = [-I | C]
-- where A == transpose matrix degrees ring V (assuming we can compute that ring!  If not, we need to fix that later).
normalDegrees = method()
normalDegrees NormalToricVariety := (V) -> (
    if not V.cache.?normalDegrees then V.cache.normalDegrees = (
        S := ring V;
        A := transpose matrix degrees S;
        d := numRows A;
        n := numColumns A; -- also numgens S.
        good := select(1, max V, a -> (d := det A_a; d == 1 or d == -1));
        if good === {} then error "cannot find maximal cone with volume 1";
        cols := good#0; -- ordered list of d column indices of the (d x n) matrix A.
        others := sort toList (set (0..n - 1) - set cols);
        perm := join(cols, others); -- or is the inverse permutation!?
        permInv := hashTable for i from 0 to #perm-1 list perm#i => i;
        permInv = for i from 0 to #perm-1 list permInv#i;
        Q := A_cols;
        Anew := Q^-1 * A_perm;
        C := Anew_{d..n-1};
        -- newdeg: take a degree in ZZ^d (old basis, using A), and change it to 
        --  a degree in ZZ^d, wrt the new basis.
        newdeg := (deg1) -> flatten entries (Q^-1 * transpose matrix {deg1});
        -- newneg: take a subset of {0..n-1} (i.e. columns of A), and return the sorted
        --  list of the corresponding columns of Anew.
        newneg := (neg1) -> sort for x in neg1 list permInv#x;
        -- mon is a row vector  0..n-1 w.r.t the new set of columns in Anew.
        -- this function translates it back to a row vector w.r.t A,
        -- and then provides a monomial in S (all with respect to the original
        -- set of variables).
        tofrac := (mon) -> (
            v1 := flatten entries (mon^permInv); 
            product for i from 0 to #v1-1 list (
                if v1_i > 0 then 
                    S_i^(v1_i) 
                else if v1_i < 0 then 
                    (1/S_i^(-v1_i)) 
                else 1_S
            ));
        -- A is a map ZZ^n --> ZZ^d
        -- Anew is as well.
        -- These induce an isomorphism ZZ^d (old) --> ZZ^d (new).
        {Anew, perm, Q, newdeg, newneg, tofrac, permInv}
        );
    V.cache.normalDegrees#0
    )

normalDegree = method()
normalDegree(NormalToricVariety, List) := (V, deg1) -> (
    normalDegrees V; -- we don't use the value, just the cached values.
    V.cache.normalDegrees#3 deg1
    )

setOfColumns = method()
setOfColumns(NormalToricVariety, List) := (V, cols) -> (
    normalDegrees V; -- we don't use the value, just the cached values.
    V.cache.normalDegrees#4 cols
    )

getFraction = method()
getFraction(NormalToricVariety, List) := (V, mon) -> (
    normalDegrees V; -- we don't use the value, just the cached values.
    V.cache.normalDegrees#5 mon
    )

-- This one is now correct, I think.
findTope = method()
findTope(Matrix, List, List) := (normalA, normalNeg, normalDeg) -> (
    n := numColumns normalA;
    d := numRows normalA;
    C := normalA_{d..n-1};
    I := id_(ZZ^(n-d));
    alpha := if #normalNeg == 0 then normalDeg else
      normalDeg + flatten entries sum(normalNeg, p -> normalA_p);
    -- create hyperplanes matrix M, RHS b, polytope will be Mx <= b
    hypers := mutableMatrix(C || -I);
    RHS := mutableMatrix transpose matrix{join(alpha, toList(n-d: 0))};
    --return (hypers, RHS);
    for p in normalNeg do (
        rowMult(hypers, p, -1);
        rowMult(RHS, p, -1);
        );
    (polyhedronFromHData(matrix hypers, matrix RHS), hypers, RHS)
    )

cohomologyFractions = method()
cohomologyFractions(NormalToricVariety, List, List) := (V, negativeSet, deg) -> (
    -- returns a list of fractions in Ring S.
    Anew := normalDegrees V;
    d := numRows Anew;
    n := numColumns Anew;
    normalNeg := setOfColumns_V negativeSet;
    beta := normalDegree(V, deg);
    normalNegDegree := (sum for i in normalNeg list flatten entries Anew_i);
    gamma := beta + normalNegDegree;
    P := first findTope(Anew, normalNeg, beta);
    if not isCompact P then error "Your negative set doesn't correspond to a cohomology cone";
    if isEmpty P then return {};
    LP := latticePoints P;
    C := Anew_{d..n-1};
    for lp in LP list (
        mon := (transpose matrix{gamma} - C * lp) || lp;
        (V.cache.normalDegrees#5 mon)/(product(negativeSet, i -> (ring V)_i))
        )
    )
--------------------------------------------------------------
-- Cohomology for complete intersections in toric varieties --
--------------------------------------------------------------

cohomologyVector = method()
cohomologyVector(NormalToricVariety, ToricDivisor) := (Y,D) -> cohomologyVector(Y, degree D)
cohomologyVector(NormalToricVariety, List) := (Y, D) -> (cohomCalg(Y, {D}); first Y.cache.CohomCalg#D)

-------------------------------------------
-- Cohomology from short exact sequences --
-------------------------------------------
cohomologyFromLES = method()
cohomologyFromLES Matrix := (M) -> (
    les := flatten entries M;
    eqns := trim splitLES les;
    M % eqns
    )

splitLES = L -> (
    alternatingSum := (P) -> (sign := 1; sum for i from 0 to #P-1 list (result := sign * P#i; sign = -sign; result));
    pos := positions(L, x -> x == 0);
    if #pos == 0 then return ideal alternatingSum L;
    pos = {-1}|pos|{#L};
    ideal for i from 1 to #pos-1 list (
        alternatingSum L_{pos#(i-1)+1..pos#i-1}
        )
    )
removeUnusedVariables = method()
removeUnusedVariables List := (L) -> (
    -- L is a list of matrices
    R := ring L#0;
    keepthese := L/support/set//sum//toList//sort;
    A := (coefficientRing R)[keepthese];
    phi := map(A,R);
    for m in L list phi m
    )

collectLineBundles = method()
collectLineBundles(NormalToricVariety, List, List) := (Y, CI, Ds) -> (
    -- Y is a projective normal toric variety
    -- CI is a list of multidegrees, i.e. elements in Cl(Y).
    -- Ds is another such list
    -- Result: a list of (r, {...}) of cohomologies that need to be computed, where r is the number of
    --  CI divisors being used (0 <= r <= #CI).
    N := #CI;
    degRank := # rays Y - dim Y;
    koszuls := for v in subsets CI list (
        if #v > 0 then degree(-sum v) else toList(degRank:0)
        );
    rsort unique flatten for D in Ds list (
        d := if instance(D, ToricDivisor) then degree D else D;
        for h in koszuls list (h+d)
        )
    )

--------------------------------------------
-- Complete Intersection in Toric Variety --
--------------------------------------------
initializeCohomologies = method(Options => {Symbol => getSymbol "a"})
initializeCohomologies(CompleteIntersectionInToric, ZZ, ZZ) := opts -> (X, nlinebundles, nbundles) -> (
    -- remember that for each line bundle, we need (#CI-1) * (dim Y + 1) + (dim X + 1) variables
    -- and for each new bundle, we need (dim X + 1) variables
    nvars := nlinebundles * ((#X.CI - 1) * (dim X.Ambient + 1) + dim X + 1) + nbundles * (dim X + 1);
    A := QQ(monoid [VariableBaseName=>opts#Symbol, Variables => 3*nvars+1200]);
    X.cache.ring = A;
    X.cache.cohom = new MutableHashTable;
    X.cache.nextVar = 0; -- before a call to one of the routines here, nextVar == nextFinalVar, but this
      -- is used to 'reserve' variables while determining LES relations.
      -- After a set of cohomology vectors is completed, then they are 'squashed down': the variables appearing
      -- (that are not finalized) are placed starting at nextFinalVar.
    X.cache.nextFinalVar = 0; -- variables before this are finalized, and appear in some cohomology
    )
basicCohomologies = method(Options=>{Limit=>null})
basicCohomologies(CompleteIntersectionInToric) := opts -> (X) -> (
    -- these include all the D_i's, and 00_X too
    -- all that is needed on Y to get Omega_X^1...
    N := #X.CI;
    Y := ambient X;
    dimY := dim Y;
    dimX := dim X;
    if not X.cache.?ring then 
      initializeCohomologies(X, #(rays Y) + 1, 2);
    ndegrees := #(rays Y) - dimY;
    Ds := prepend(Y_0-Y_0, for i from 0 to #(rays Y)-1 list (-Y_i));
    Ds = join(Ds, for h in X.CI list -h);
    Ds = Ds/degree;
    linebundles := collectLineBundles(Y, X.CI, Ds);
    cohomCalg(Y, linebundles); -- actually computes the cohomologies needed over Y, cohoms is a MutableHashTable.
    for d in Ds do cohomologyVector(X, d);
    X.cache.cohom
    )
nextLineBundle = method()
nextLineBundle CompleteIntersectionInToric := (X) -> (
    -- returns a list of (#CI-1) cohom vectors of length: dim Y + 1
    --  and the last one has the same length, but is zero above dim X.
    Y := X.Ambient;
    CI := X.CI;
    nvecs := #CI-1;
    A := X.cache.ring;
    firstvecs := for i from 0 to #CI-2 list (
        result := flatten entries genericMatrix(A, A_(X.cache.nextVar), 1, dim Y + 1);
        X.cache.nextVar = X.cache.nextVar + dim Y + 1;
        result
        );
    lastvec := flatten entries genericMatrix(A, A_(X.cache.nextVar), 1, dim X + 1);
    lastvec = join(lastvec, toList(#CI : 0));
    X.cache.nextVar = X.cache.nextVar + dim X + 1;
    append(firstvecs, lastvec)
    )
nextBundle = method()
nextBundle CompleteIntersectionInToric := (X) -> (
    -- returns a single list, of length: dim X + 1, representing the
    -- cohomoogies of a sheaf or vector bundle on X.
    A := X.cache.ring;
    result := flatten entries genericMatrix(A, A_(X.cache.nextVar), 1, dim X + 1);
    X.cache.nextVar = X.cache.nextVar + dim X + 1;
    result
    )
finalizeVariables = (X, L) -> (
    -- X is a CompleteIntersectionInToric
    -- L: List of Matrix
    --  each one should be a matrix over X.cache.ring.
    -- Returned value: L': a list of matrices over the same ring, 
    --   where the variables have been moved
    --   up to not waste ring indeterminants.
    -- This updates X.cache.nextFinalVar, X.cache.nextVar as well.
    -- YYY
    R := ring L#0;
    firstvar := X.cache.nextFinalVar;
    -- We make a list of the variables that actually occur here
    keepthese := L/support/set//sum;
    keepthese = keepthese - set for i from 0 to firstvar-1 list X.cache.ring_i;
    keepthese = sort toList keepthese;
    -- now we create a ring map that maps these existing variables to their compacted counterpart
    phi := map(R,R,for v from 0 to #keepthese-1 list (keepthese#v => R_(firstvar + v)));
    X.cache.nextVar = X.cache.nextFinalVar = firstvar + #keepthese;
    for m in L list phi m
    )
-- Interface for cohomology in toric complete intersections

cohomologyVector(CompleteIntersectionInToric, List) := (X, degreeD) -> (
    if not X.cache.?cohom then basicCohomologies X;
    if not X.cache.cohom#?degreeD then X.cache.cohom#degreeD = (
        CI := X.CI;
        N := #CI;
        Y := ambient X;
        dimY := dim Y;
        linebundles := collectLineBundles(Y, CI, {degreeD});
        cohoms := cohomCalg(Y, linebundles); -- actually computes the cohomologies needed over Y, cohoms is a MutableHashTable.
        -- Next step is to create the short exact sequences we need.
        H0 := partition(s -> #s, subsets CI);
        H := applyPairs(H0, (r,v) -> (r, if r > 0 then apply(v, v0 -> degreeD - degree sum v0) else apply(v, v0->degreeD)));
        -- H is a hash table: keys are 0..#CI, values: lists of degrees of line bundles at that step in Koszul complex.
        -- Compute all of these cohomologies:
        H1 := applyPairs(H, (r,v) -> (r, sum for v1 in v list cohomologyVector(Y,v1)));
        -- Now create a ring with #CI*(dim Y + 1) number of variables
        -- L contains the ansatz cohomology vectors for the syzygy vector bundles
        L := reverse nextLineBundle X;
        ses := prepend(matrix transpose{H1#N, H1#(N-1), L#(N-1)},
            reverse for i from 0 to N-2 list matrix transpose{L#(i+1), H1#i, L#i});
        J := trim (sum for s in ses list splitLES flatten entries s);
        ses1 := for s in ses list (s%J);
        if debugLevel > 0 then << "exact sequences: " << ses1 << endl;
        M := (last ses1)_{2};
        M = first finalizeVariables(X,{M});
        result := take(flatten entries M, dim X + 1);
        if all(result, x -> liftable(x, ZZ))
        then result = result/(x -> lift(x,ZZ));
        result
        );
    X.cache.cohom#degreeD
    )
cohomologyVector(CompleteIntersectionInToric, ToricDivisor) := (X, D) -> cohomologyVector(X, degree D)
cohomologyVector CompleteIntersectionInToric := (X) -> cohomologyVector(X, 0 * (ambient X)_0)

cohomologyOmega1 = method()
cohomologyOmega1(CompleteIntersectionInToric) := (X) -> (
    -- Need to make 2 ses's, and two cohom vectors for dim X
    -- 
    if not X.cache.?cohom then basicCohomologies X;
    if not X.cache.cohom#?"omega1" then X.cache.cohom#"omega1" = (
        dimX := dim X;
        Y := ambient X;
        CI := X.CI;
        ndegrees := #(rays Y) - dim Y;
        cohoms := basicCohomologies X;
        cohomOOX := cohomologyVector(X, degree (0*Y_0));
        spot1 := sum for i from 0 to #(rays Y)-1 list cohomologyVector(X, - degree Y_i);
        spot2 := sum for h in CI list cohomologyVector(X, - degree h);
        -- note: we computed all the cohomologies above, now we consider vec1, vec2
        -- these are created after the above 4 lines, so that when we call finalizeVariables 
        --in cohomologyVector, they don't conflict with our choice of vec1 and and vec2 variables.
        vec1 := nextBundle X;
        vec2 := nextBundle X;
        ses1 := {
            vec1,
            spot1,
            ndegrees * cohomOOX};
        ses2 := {
            spot2,
            vec1,
            vec2
            };
        sess := {transpose matrix ses1, transpose matrix ses2};
        J1 := ideal(vec2_0 - cohomOOX_1, vec2_(dimX) - cohomOOX_(dimX-1));
        J := trim (J1 + (sum for s in sess list splitLES flatten entries s));
        sess = for s in sess list (s%J);
        if debugLevel > 0 then << "exact sequences: " << sess << endl;
        M := (last sess)_{2};
        M = first finalizeVariables(X,{M});
        vec := flatten entries M;
        if all(vec, x -> liftable(x,ZZ))
        then vec = vec/(x -> lift(x,ZZ));
        vec
        );
    X.cache.cohom#"omega1"
    )

hodgeDiamond = method()
hodgeDiamond CompleteIntersectionInToric := (X) -> (
    -- Assumptions: Y is smooth?
    --  Certainly want: X is smooth.
    -- if dimX <= 3, then we only need Omega1_X
    -- if dimX == 4, then we can either use Omega2_X, or the topological Euler characteristic
    dimX := dim X;
    if dimX >= 4 then <<  "warning: not yet implemented for dimension >= 4, -1's mean not computed" << endl;
    Y := ambient X;
    vec1 := cohomologyOmega1 X;
    vec0 := cohomologyVector(X, degree(0*Y_0));
    matrix for p from 0 to dimX list for q from 0 to dimX list (
        if p == 0 then vec0_q 
        else if q == 0 then vec0_p
        else if p == 1 then vec1_q
        else if q == 1 then vec1_p
        else if p == dimX then vec0_(dimX-q)
        else if q == dimX then vec0_(dimX-p)
        else if p == dimX-1 then vec1_(dimX-q)
        else if q == dimX-1 then vec1_(dimX-p)
        else -1
        )
    )

-- Ds are toric divisors in a toric variety
-- Computes the cohomology vector of the intersection
-- of these divisors.
hodgeVector = method()
hodgeVector List := (Ds) -> (
    if #Ds == 0 then error "expected at least one divisor";
    V := variety Ds#0;
    if not all(Ds, D -> instance(D,ToricDivisor) and variety D === V)
      then error "expected divisors all on the same toric variety";
    D := completeIntersection(V, Ds);
    cohomologyVector(D, degree (V_0-V_0))
    )

-- The following uses cohomology matrix and bases to compute
-- the cohomology dimensions.  This gets exact values, based on specific
-- polynomials in the Cox ring, however it needs specific polynomials.

makeCohomMatrix = method()
makeCohomMatrix(ZZ, NormalToricVariety, RingElement, RingElement) := (i,V,F,G) -> (
   (m1, base1, base2) := cohomologyMatrix(i,V,-degree F, G);
   (m2, base3, base4) := cohomologyMatrix(i,V,-degree G, F);
   if base2 != base4 then error "hmm, expected bases to be identical for columns";
   (matrix m1 || matrix m2, join(base1, base3), base2)
   )
cohomVectorOfCodim2 = method()
cohomVectorOfCodim2(NormalToricVariety, RingElement, RingElement) := (V,F,G) -> (
    maps := for i from 0 to dim V list makeCohomMatrix(i,V,F,G);
    rks := for i from 0 to dim V list (#maps_i_1, rank maps_i_0, #maps_i_2);
    if rks_0 != (0,0,0) then << "warning: expected H^0 part to be all zero" << endl;
    if rks_1_2 != rks_1_1 then << "warning: expected H^1 map to be an inclusion" << endl;
    if rks_(dim V)_1 != rks_(dim V)_0 then << "warning: expected H^" << dim V << " map to be a surjection" << endl;
    for i from 0 to dim V - 2 list (
        rks_(i+1)_0 + rks_(i+2)_2 - rks_(i+1)_1 - rks_(i+2)_1 + if i == 0 then 1 else 0
        )
    )
hodgeVector(List, List) := (Ds, Fs) -> (
    -- Ds: list of (effective) divisor classes on a toric variety
    -- Fs: list of polynomials, one for each class.
    -- currently: these are limited to #Ds == #Fs == 2
    -- invariant: degree Ds_i == degree Fs_i, all i.
    if #Ds =!= 2 or #Fs =!= 2 then error "expected codimension 2 complete intersection in a toric variety";
    if degree Ds_0 =!= degree Fs_0 or degree Ds_1 =!= degree Fs_1 then
      error "expected polynomials to be in the given divisor classes";
    V := variety Ds_0;
    cohomVectorOfCodim2(V, Fs_0, Fs_1)
    )

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
installPackage "StringTorics" -- ERROR at end.  Not sure what it is...

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
  
  
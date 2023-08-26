-- TODO: compute toricMoriCone without intersection ring?
-- 

----------------------------------------------------------------  
-- gvInvariants ------------------------------------------------
-- Uses computeGV.cpp from CYtools -----------------------------
----------------------------------------------------------------
-- moriCone = method()
-- -- Not functional...
-- moriCone NormalToricVariety := List => (V) -> (
--     IV := intersectionRing (abstractVariety V);
--     Cs := matrix for x in orbits(V, 1) list (
--         c := product(x, i -> IV_i);
--         for d in gens IV list integral(c*d)
--         );
--     M := posHull transpose lift(Cs, QQ);
--     GLSM := matrix degrees ring V;
--     entries transpose((rays M) // GLSM)
--     )


toricMoriCone(NormalToricVariety, List) := Cone => (V, basisIndices) -> (
    IV := intersectionRing (abstractVariety V);
    Cs := matrix for x in orbits(V, 1) list (
        c := product(x, i -> IV_i);
        for j in basisIndices list integral(c * IV_j)
        );
    posHull transpose lift(Cs, QQ) -- TODO: lift to ZZ?
    )

toricMoriCone CalabiYauInToric := Cone => X -> (
    -- TODO: handle toric mori cones of non-favorables
    if isFavorable X then toricMoriCone(ambient X, basisIndices X)
    )

hilbertBasisGenerators = method()
hilbertBasisGenerators Cone := List => C -> (
    for x in hilbertBasis C list flatten entries x
    )

-- This function returns a very large heft vector.  Not so good!
heft CalabiYauInToric := List => X -> (
    C := toricMoriCone X;
    sum entries transpose rays dualCone C
    )

gvInvariants = method(Options => {
    Mori => null, -- null means: compute rays of the Mori cone of V (in ZZ^(h11))
    Heft => null, -- null means: compute it
    DegreeLimit => infinity,
    Precision => 150,
    FilePrefix => "foo",
--    Executable => "~/src/git-from-others/cytools-private/external/gv/computeGV-good/computeGV",
    Executable => "~/src/git-from-others/cytools-private/external/gv/computeGV-good/computeGV",
    KeepFiles => true
    })

-- The function to write the data needed by the computeGV program
gvInput = (moriGenerators, heftval, GLSM, intersectionnums, degreelimit, prec) -> (
    -- moriGenerators: list of lists. Hilbert basis of the cone of 
    --   irreducible curves induced from the toric variety.
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

-- TODO: make the gvInvariants code not go through NormalToricVarieties.
--  and use only info obtained from data we have in CalabiYauInToric.
--  requires: intersectionNumbersOfCY
--            toricMoriCone

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
    mori := if opts.Mori =!= null then 
                opts.Mori 
            else 
                hilbertBasisGenerators toricMoriCone(V, basisIndices);
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
    if not isFavorable X then return null;
    intersectionnums := for t in intersectionNumbers X list append(t#0, t#1);
    mori := if opts.Mori =!= null then 
                opts.Mori 
            else 
                hilbertBasisGenerators toricMoriCone(ambient X, basisIndices X);
    heft := if opts.Heft =!= null then opts.Heft else (
      sum entries transpose rays dualCone posHull transpose matrix mori
      );
    -- OK, now we have computed everything we need.  Write it to a file
    infile := opts.FilePrefix | "-input";
    outfile := opts.FilePrefix | "-output";
    infile << gvInput(mori, heft, transpose degrees X, intersectionnums,
        opts.DegreeLimit, opts.Precision) << close;
    inputLine := opts.Executable | " <" | infile | " >" | outfile;
    print inputLine;
    run inputLine; -- TODO: run this as a program and if it crashes, return something reasonable.
    -- Get the output, package as a hash table
    contents := get outfile;
    if #contents == 0 then return null;
    (lines contents)/value//hashTable
    )

gvRay = method(Options => options gvInvariants)
gvRay(CalabiYauInToric, List) := HashTable => opts -> (X, curveClass) -> (
    -- This doesn't seem to be correct
    if not isFavorable X then return null;
    return gvInvariants(X,Mori => {curveClass})
    )

gvCone = method(Options => options gvInvariants)
gvCone CalabiYauInToric := Cone => opts -> X -> (
    if not isFavorable X then return null;
    gv := gvInvariants(X, opts);
    if gv === null then return null;
    posHull transpose matrix ((keys gv)/toList)
    )

gvInvariantsAndCone = method(Options => options gvInvariants)
gvInvariantsAndCone(CalabiYauInToric, ZZ) := Sequence => opts -> (X, D) -> (
    -- D is the degree bound to start with.  We could start with 5, or DegreeLimit/2 or DegreeLimit/4, or ...
    if not isFavorable X then return null;
    degvec := heft X;
    gv := gvInvariants(X, opts);
    if gv === null then return null;
    keysgv := keys gv;
    H := hashTable for k in keysgv list k => dotProduct(k, degvec);
    firstSet := select(keys H, k -> H#k <= D);
    if debugLevel > 0 then << "The number of curves in the first set: " << #firstSet << endl;
    C := posHull transpose matrix (firstSet);
    Cdual := dualCone C;
    HC := transpose rays Cdual;
    curves := for k in keys H list if H#k > D then transpose matrix {k} else continue;
    set2 := select(curves, c  -> any(flatten entries (HC * c), a -> a < 0));
    if debugLevel > 0 then << "The number of curves not in the first cone: " << #set2 << endl;
    C2 := if #set2 == 0 then C else posHull (rays C | matrix{set2});
    if debugLevel > 0 and #set2 == 0 then (
        << "CY " << label X << " C = " << rays C  << endl
        )
    else
        << "*differs* CY " << label X << " C1 = " << rays C << " and C2 = " << rays C2 << endl;
    (gv, C2)
    )

partitionGVConeByGV = method(Options => options gvInvariants)
partitionGVConeByGV CalabiYauInToric := HashTable => opts -> X -> (
    -- return null if we cannot computr GV invariants (i.e. if non-favorable).
    if not isFavorable X then return null;
    gv := gvInvariants(X, opts); -- TODO: stash this?
    if gv === null then return null;
    C := posHull transpose matrix ((keys gv)/toList);
    gvX := entries transpose rays C;
    partition(f -> if gv#?(toSequence f) then gv#(toSequence f) else 0, gvX)
    )

partitionGVConeByGV(CYToolsCY3, ZZ) := HashTable => opts -> (X, D) -> (
    -- return null if we cannot computr GV invariants (i.e. if non-favorable).
    (gv, C) := gvInvariantsAndCone(X, D, opts);
    gvX := entries transpose rays C;
    partition(f -> if gv#?f then gv#f else 0, gvX)
    )

-- TODO: move to Topology.m2? file?
findLinearMaps = method()
findLinearMaps(HashTable, HashTable) := List => (gv1, gv2) -> (
    -- gv1, gv2: result of partitionGVConeByGV
    if sort keys gv1 =!= sort keys gv2 then return {};
    for k in keys gv1 do if #gv1#k =!= #gv2#k then return {};
    n := # (first values gv1)_0; -- we should check if all the values are lists of integers of this size.
    t := symbol t;
    T := QQ[t_(1,1)..t_(n,n)];
    M := genericMatrix(T, n, n);
    -- now we make the ideals for each key, and each permutation.
    ids := for k in keys gv1 list (
        perms := permutations(#gv1#k);
        mat1 := transpose matrix gv1#k;
        mat2 := transpose matrix gv2#k;
        for p in perms list (
            I := trim ideal (M * mat1 - mat2_p); 
            if I == 1 then continue else I
            )
        );
    topval := ids/(x -> #x - 1);
    zeroval := ids/(x -> 0);
    fullIdeals := for a in zeroval .. topval list (
        J := trim sum for i from 0 to #ids-1 list ids#i#(a#i);
        if J == 1 then continue else J
        );
    Ms := for i in fullIdeals list M % i;
    --newMs := select(Ms, m -> (d := det m; d == 1 or d == -1));
    --if any(newMs, m -> support m =!= {}) then << "some M is not reduced to a constant" << endl;
    Ms
    )

----------------------------
-- new code ----------------
----------------------------
-- GVInvariants class
--  has X in cache, or way to rerun gv's at higher invarisnts, rays, etc.
--  

-- Design: What should a GVInvariants class look like?
--  1. Has hash table as it does now.
--  2. Knows its degree limit, and grading vector.
--  3. Can compute "infinity cone": actually, should be done for 2 or 3 different degrees,
--       then compare them?
--  4. Compute ray of GV values out some distance.
--    This should use special features of the code?  Does it work on non-extremal rays?
-- The following is perhaps not part of this class...
--  5. Determine what kind of extremal ray a ray is:
--    1. nilpotent (type I)
--    2. nilpotent (type II0, type IIg)
--    3. potent ray.
--    4. is a ray in the closure of the infinity cone?  Or can we not consider this possibility?
-- For non-general CY3's it is possible for a curve to be effective, but have gv ray all 0's.

-- GVInvariants = new Type of HashTable
-- gvInvariantsObject = method(Options => {
--         Mori => null, 
--         Heft => null,
--         DegreeLimit => 5,
--         Precision => 150,
--         FilePrefix => "foo",
--         Executable => "~/src/git-from-others/cytools-private/external/gv/computeGV",
--         KeepFiles => true
--     })
--gvInvariantsObject()
 -- need: intersection numbers
 --       mori cone hilbert basis gens
 --       degrees

moriConeGVs = method()
moriConeGVs(CalabiYauInToric, ZZ) := (X, deglimit) ->(
    -- deglimit that the GV invariants were computed to.
    -- loop thru the toric mori cone cap generators, and for each,
    -- look at the gv ray. -- then return a hash table whose keys are among
    --   {POTENT, {gv vals on ray}
    --   {FLOP, {gv vals on ray}, 
    --   {TYPEII0, {gv vals on ray},
    --   {TYPEIIg, {gv vals on ray}}
    --   {ZERO} -- this means that 
    --   {...}, C^perp has D^3 = 0.  Not sure what the gv invariants are in this case...
    -- and whose values are the list of curve classes with that type.
    )

gvRay(HashTable, List, ZZ, List) := opts -> (GVHash, C, deglimit, degvector) -> (
    contentC := gcd C;
    if contentC =!= 1 then C = C // contentC;
    d := dotProduct(degvector, C);
    rayC := for i from 1 to floor(deglimit/d) list (
        Cseq := toSequence(i*C);
        if GVHash#?Cseq then GVHash#Cseq else 0
        );
    rayC
    )

count = 0;

classifyExtremalCurve = method()
classifyExtremalCurve(HashTable, List, ZZ, List) := (GVHash, C, deglimit, degvector) -> (
    rayC := gvRay(GVHash, C, deglimit, degvector);
    if #rayC <= 2 then return {"OTHER", rayC};
    if all(2..#rayC-1, i -> rayC#i == 0) then (
        -- only first two, possibly, are non-zero.
        if rayC#0 == 0 and rayC#1 == 0 then (count=count+1; return {"ZERO", count});
        if rayC#0 == -2 or rayC#1 == -2 then return {"TYPEII0", {rayC#0, rayC#1}};
        if rayC#0 >= 0 and rayC#1 >= 0 then return {"FLOP", {rayC#0, rayC#1}};
        if rayC#0 < 0 or rayC#1 < 0 then return {"TYPEIIg", {rayC#0, rayC#1}};
        )
    else return {"POTENT", {rayC#0, rayC#1, rayC#2, "..."}}
    )



///
  restart
  debug needsPackage "StringTorics" -- the debug is because some functions are not yet exported.
  DB3 = "../Databases/cys-ntfe-h11-3.dbm"
  RZ = ZZ[a,b,c]
  RQ = QQ (monoid RZ);
  (Qs, Xs) = readCYDatabase(DB3, Ring => RZ);

  moris = for k in sort keys Xs list (
      X = Xs#k;
      if not isFavorable X then continue; -- only handles favorable polytopes and CY3's.
      classifyExtremalCurves(X, Verbose => 2)
      );

  moris/keys/sort//unique
      degvec = heft X;
      deglimit = max for c in toricMoriConeCap X list 3 * dotProduct(degvec, c);
      << "degree limit for " << k << " is " << deglimit << endl;
      gvX = gvInvariants(X, DegreeLimit => deglimit);
      moriClass = sort for c in toricMoriConeCap X list
          classifyExtremalCurve(gvX, c, deglimit, degvec);
      print netList (ans := prepend(k, moriClass));
      moriClass
      )     

  restart
  debug needsPackage "StringTorics" -- the debug is because some functions are not yet exported.
  DB4 = "../Databases/cys-ntfe-h11-4.dbm"
  RZ = ZZ[a,b,c,d]
  RQ = QQ (monoid RZ);
  (Qs, Xs) = readCYDatabase(DB4, Ring => RZ);
  
  X = Xs#(137,0)
  for k in sort keys Xs list (
  deglimit = 14
  degvec = heft X
  gvX = gvInvariants(X, DegreeLimit => deglimit);

  netList for c in toricMoriConeCap X list
    sort classifyExtremalCurve(gvX, c, deglimit, degvec)

  for k in sort keys Xs list (
      X = Xs#k;
      if not isFavorable X then continue; -- only handles favorable polytopes and CY3's.
      deglimit = 14;
      degvec = heft X;
      gvX = gvInvariants(X, DegreeLimit => deglimit);
      moriClass = sort for c in toricMoriConeCap X list
          classifyExtremalCurve(gvX, c, deglimit, degvec);
      print moriClass;
      moriClass
      )     
  moricap = toricMoriConeCap X
  for c in toricMoriConeCap X list
    c => gvRay(gvX, c, 22, heft X)
  netList oo
  

///

gvTopMoriConeCapDegree = method()
gvTopMoriConeCapDegree CalabiYauInToric := X -> (
    if not isFavorable X then error "expected a favorable polytope";
    degvec := heft X;
    max for c in toricMoriConeCap X list dotProduct(degvec, c)
    )

classifyExtremalCurves = method(Options => {
        Verbose => 0,
        DegreeLimit => null,
        MoriHilbertGens => null
        })
classifyExtremalCurves(HashTable, List, ZZ, List) := (GVHash, Cs, deglimit, degvector) -> (
    partition(c -> classifyExtremalCurve(GVHash, c, deglimit, degvector), Cs)
    )

classifyExtremalCurves CalabiYauInToric := opts -> X -> (
    if not isFavorable X then error "expected favorable CY3-fold";
    mori := if opts.MoriHilbertGens === null then toricMoriConeCap X else opts.MoriHilbertGens;
    deglimit := 3 * gvTopMoriConeCapDegree X;
    if opts.Verbose > 1 then << "*** mori cone cap degree limit is " << deglimit << " ***" << endl;
    degvec := if opts.DegreeLimit === null then heft X else opts.DegreeLimit;
    gvX := gvInvariants(X, DegreeLimit => deglimit);
    partition(c -> classifyExtremalCurve(gvX, c, deglimit, degvec), mori)
    )

///
  restart
  debug needsPackage "StringTorics"
  DB4 = "../Databases/cys-ntfe-h11-4.dbm"
  RZ = ZZ[a,b,c,d]
  RQ = QQ (monoid RZ);
  (Qs, Xs) = readCYDatabase(DB4, Ring => RZ);

  moris = for k in sort keys Xs list elapsedTime (
      X = Xs#k;
      if not isFavorable X then continue; -- only handles favorable polytopes and CY3's.
      mori1 = classifyExtremalCurves(X, Verbose => 2);
      print netList {mori1};
      mori1
      );


  X = Xs#(34,0)
  assert isFavorable X
  classifyExtremalCurves(Xs#(34,0), Verbose => 2)
  classifyExtremalCurves(Xs#(35,0), Verbose => 2)
  mori = toricMoriConeCap X;
  deglimit = 3 * gvTopMoriConeCapDegree X;
  deglimit
  degvec = heft X;
  gvX = gvInvariants(X, DegreeLimit => deglimit)
  partition(c -> classifyExtremalCurve(gvX, c, deglimit, degvec), mori)
///
-- GVInvariantsTable = new Type of HashTable

-- makeGVInvariantsTable = (intersectionNums, moriGens, GLSM, heftvec, deglimit) -> (
--     new GVInvariantsTable from {
--         IntersectionNumbers => intersectionNums,
--         MoriHilbertGens => moriGens, -- a list of curve classes (each a list of n integers).
--         Degrees => GLSM, -- format: a list if all the toric degrees.
--         Heft => heftvec, -- format: a list of integers, of same length as each GLSM degree
--         cache => new CacheTable from {
--             "GV" =>  new MutableHashTable -- keys: (degreelimit, precision), value: a hash table c => gv.
--             }
--         }
--     )

-- gvInvariantsTable = method()

-- gvInvariantsTable CalabiYauInToric := GVInvariantsTable => X -> (
--     makeGVInvariantsTable (
--     )

moriConeCapGVInvariants = method()
moriConeCapGVInvariants CalabiYauInToric := X -> (
    )
    )

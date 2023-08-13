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
    Executable => "~/src/git-from-others/cytools-private/external/gv/computeGV",
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
    run inputLine;
    -- Get the output, package as a hash table
    (lines get outfile)/value//hashTable
    )

gvCone = method(Options => options gvInvariants)
gvCone CalabiYauInToric := Cone => opts -> X -> (
    gv := gvInvariants(X, opts);
    posHull transpose matrix ((keys gv)/toList)
    )

gvInvariantsAndCone = method(Options => options gvInvariants)
gvInvariantsAndCone(CalabiYauInToric, ZZ) := Sequence => opts -> (X, D) -> (
    -- D is the degree bound to start with.  We could start with 5, or DegreeLimit/2 or DegreeLimit/4, or ...
    if not isFavorable X then return null;
    degvec := heft X;
    gv := gvInvariants(X, opts);
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

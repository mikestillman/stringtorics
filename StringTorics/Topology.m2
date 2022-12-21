-- Code to compute the topology of a CY 3-fold.
-- Also contains code to help determine whether two CY3's are homeomorphic
-- (this appears to be equivalent to being diffeomorphic).

topologicalData = method()
topologicalData CYData := TopologicalDataOfCY3 => X -> (
    -- TODO: this does not consider torision in H_2(X, ZZ) or H_3(X, ZZ)
    elapsedTime new TopologicalDataOfCY3 from {
        "h11" => hh^(1,1) cyPolytopeData X,
        "h21" => hh^(2,1) cyPolytopeData X,
        "c2" => c2 X,
        "intersection numbers" => intersectionNumbers X
        }
    )

-- this is the older, alternate version of this function.
-- this is to be removed.
topologicalData(CYData, Ring) := TopologicalDataOfCY3 => (X, RZ) -> (
    V := ambient X;
    Q := X#"polytope data";
    P := polytope Q;
    data := elapsedTime topologyOfCY3(V, basisIndices X);
    -- this data above computes intersection numbers for all toric divisors. 
    -- So we consider only the ones whose indices are contained in basis indices:
    new TopologicalDataOfCY3 from {
        "h11" => elapsedTime hh^(1,1) Q,
        "h21" => elapsedTime hh^(2,1) Q,
        "c2" => sub(data_3, vars RZ),
        "cubic intersection form" => sub(data_2, vars RZ)
        }
    )

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

pointCount = method()
pointCount(RingElement, ZZ) := ZZ => (F, p) -> (
    -- F is a polynomial in 3 variables (FIXME: any number of variables)
    -- p is a prime number
    R := (ZZ/p) (monoid ring F);
    Fp := sub(F, R);
    # for x in (0,0,0)..(p-1,p-1,p-1) list if sub(Fp, matrix{{x}}) == 0 then x else continue
    )

invariants CYData := List => X -> (
    L := c2Form X;
    F := cubicForm X;
    h11 := hh^(1,1) X;
    h12 := hh^(1,2) X;
    RZ := ring L;
    if ring F =!= RZ then 
        error "expected same rings for c2 and cubic form";
    if coefficientRing RZ =!= ZZ then 
        error "expected c2 and cubic form ring to be a polynomial ring over ZZ";
    RQ := QQ[gens RZ];
    facs := select((factors F)/toList/last, g -> support g != {});
    LQ := sub(L, RQ);
    FQ := sub(F, RQ);
    d := dim saturate ideal jacobian FQ;
    nc := # decompose ideal(LQ, FQ);
    singZ := flatten entries gens gb saturate(ideal(F) + ideal jacobian F);
    badp := select(singZ, a -> support leadTerm a === {});
    badp = if badp === {} then 0 else sub(first badp, ZZ);
    ptcounts := for p in {2, 3, 5, 7, 11} list pointCount(F, p);
    {h11, h12, ptcounts, badp, (trim content L)_0, (trim content F)_0, #facs, d, nc}
    )

invariants1 = method()
invariants1 CYData := List => X -> (
    L := c2Form X;
    F := cubicForm X;
    h11 := hh^(1,1) X;
    h12 := hh^(1,2) X;
    contentL := (trim content L)_0;
    contentF := (trim content F)_0;
    ptcounts := for p in {2, 3, 5, 7, 11} list pointCount(F, p);
    {h11, h12} |  ptcounts | {contentL, contentF}
    )

invariants2 = method()
invariants2 CYData := List => X -> (
    L := c2Form X;
    F := cubicForm X;
    h11 := hh^(1,1) X;
    h12 := hh^(1,2) X;
    RZ := ring L;
    if ring F =!= RZ then 
        error "expected same rings for c2 and cubic form";
    if coefficientRing RZ =!= ZZ then 
        error "expected c2 and cubic form ring to be a polynomial ring over ZZ";
    RQ := QQ[gens RZ];
    facs := select((factors F)/toList/last, g -> support g != {});
    LQ := sub(L, RQ);
    FQ := sub(F, RQ);
    d := dim saturate ideal jacobian FQ;
    nc := # decompose ideal(LQ, FQ);
    singZ := flatten entries gens gb saturate(ideal(F) + ideal jacobian F);
    badp := select(singZ, a -> support leadTerm a === {});
    badp = if badp === {} then 0 else sub(first badp, ZZ);
    {h11, h12, badp, (trim content L)_0, (trim content F)_0, #facs, d, nc}
    )

mapIsIsomorphism = method()
mapIsIsomorphism(Matrix, CYData, CYData) := Boolean => (M, X1, X2) -> (
    -- M is a matrix over the base field, a possible map giving
    -- an isomorphism of topologies.
    -- T1, T2 are two topologies.
    F1 := cubicForm X1;
    F2 := cubicForm X2;
    RZ := ring F1;
    if RZ =!= ring F2 then error "expected the same picard ring";
    phi := map(RZ, RZ, M);
    (phi c2Form X1 == c2Form X2) and (phi F1 == F2)
    )

partitionByTopology = method()
partitionByTopology List := LGVs -> (
    -- LGVs is a list of X => gvPartition.
    -- output: a hashtable, keys are labels, values are lists of {label, matrix}
    labels := for x in LGVs list label first x;
    hashXs := hashTable for y in LGVs list (label first y) => y;
    distinctTops := new MutableHashTable; -- label => list of labels.
    for i in labels do (
        Xi := first hashXs#i;
        GVi := last hashXs#i;
        << "trying " << i << endl;
        prev := keys distinctTops;
        isFound := false;
        for j in prev do (
            Xj := first hashXs#j;
            GVj := last hashXs#j;
            -- compare CY's i, j
            Ms := findLinearMaps(GVj, GVi);
            if Ms === {} then continue;
            Ms = Ms/(m -> lift(m, ZZ));
            Ms = select(Ms, m -> (d := det m; d === 1 or d === -1));
            isIsos := Ms/(m -> mapIsIsomorphism(m, Xj, Xi));
            if any(isIsos, x -> true) then (
                mi := position(isIsos, x -> true);
                << "found isomorphism" << endl;
                distinctTops#j = append(distinctTops#j, {i, Ms#mi});
                isFound = true;
                break;
                ));
        if not isFound then (
            distinctTops#i = {};
            << "found new top: " << topologicalData Xi << endl;
            );
        );
    new HashTable from distinctTops
    )

partitionByTopology(List, HashTable, ZZ) := (Ls, Xs, degreelimit) -> (
    -- LGVs is a list of X => gvPartition.
    -- output: a hashtable, keys are labels, values are lists of {label, matrix}
    labels := Ls;
    if #Ls === 1 then return hashTable {Ls#0 => {}};
    hashXs := Xs;
    GVs := hashTable for lab in labels list lab => partitionGVConeByGV(hashXs#lab, DegreeLimit => degreelimit);
    distinctTops := new MutableHashTable; -- invariants => list of labels.
    for i in labels do (
        Xi := hashXs#i;
        GVi := GVs#i;
        << "trying " << i << endl;
        prev := keys distinctTops;
        isFound := false;
        for j in prev do (
            Xj := hashXs#j;
            GVj := GVs#j;
            -- compare CY's i, j
            Ms := findLinearMaps(GVj, GVi);
            if Ms === {} then continue;
            Ms = Ms/(m -> lift(m, ZZ));
            Ms = select(Ms, m -> (d := det m; d === 1 or d === -1));
            isIsos := Ms/(m -> mapIsIsomorphism(m, Xj, Xi));
            if any(isIsos, x -> true) then (
                mi := position(isIsos, x -> true);
                << "found isomorphism" << endl;
                distinctTops#j = append(distinctTops#j, {i, Ms#mi});
                isFound = true;
                break;
                ));
        if not isFound then (
            distinctTops#i = {};
            << "found new top: " << i << endl;
            );
        );
    new HashTable from distinctTops
    )

-- TODO: remove the following code (any reason to keep it?)
topologyOfCY3 = method(Options => {
        Variable => "x",
        Ring => null
        })

topologyOfCY3(NormalToricVariety, List) := opts -> (V, basisIndices) -> (
    -- input: 
    --   V: a simplicial resolution of a Fano toric 4-fold
    --      X is a (general) anti-canonical section of V.
    --   basisIndices: list of integer indicesas to which V_i will be in the 
    --      basis of Pic X that you choose.
    -- output: a hash table containing:
    --  a. triple intersection numbers (a hash table, H#{a,b,c}, with 0 <= a <= b <= c < h11(X))
    --  b. the h11 numbers: c2(X) . D_i, 0 <= i < h11
    --  c. the integers h11, h12
    --  d. the cubic form C(x,y,z) in a polynomial ring ZZ[x_0, ..., x_(h11-1)]
    --  e. a linear form L(x,y,z) in the same ring, representing c2(X).D_i
    --
    P := convexHull transpose matrix rays V;
    h11 := h21OfCY P; -- we want h11 of `polar P`.
    h21 := h11OfCY P;
    if #basisIndices != h11 then error("expected "|h11|" indices");
    H := CY3NonzeroMultiplicities V;
    -- basisInv := new MutableHashTable;
    -- for i from 0 to #basisIndices-1 do basisInv#(basisIndices#i) = i;
    -- H3 := hashTable for x in keys H list (
    --     if isSubset(x, basisIndices) then (
    --         x' := apply(x, i -> basisInv#i);
    --         x' => H#x 
    --         ) else continue
    --     );
    -- Now let's get the cubic form and the linear form directly from the intersection theory.
    -- For larger h11, this method will need to change.
    x := getSymbol opts.Variable;
    pt := base(x_0..x_(h11-1));
    A := intersectionRing pt; -- over QQ
    R := if opts#Ring =!= null then opts#Ring else ZZ (monoid A);
    if numgens R =!= h11 then error("expected a ring with "|toString h11|" variables");

    X := completeIntersection(V, {-toricDivisor V});
    Xa := abstractVariety(X, pt);
    IX := intersectionRing Xa;
    h := sum(h11, i -> A_i * IX_(basisIndices#i));
    C := sub(integral(h^3), vars R);
    L := integral((chern_2 tangentBundle Xa) * h);
    L = sub(L, vars R);
    (h11, h21, C, L)
    )

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

c2 TopologicalDataOfCY3 := T -> T#"c2"
cubicForm TopologicalDataOfCY3 := T -> T#"cubic intersection form"
----- end of removing code TODO -----------------------------

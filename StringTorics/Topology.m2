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

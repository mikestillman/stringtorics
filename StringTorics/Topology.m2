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

allPoints = (p, n) -> (
    -- all points in kk = ZZ//p in kk^n
    pts := for a from 0 to p-1 list {a};
    if n == 1 then return pts;
    if n <= 0 then error "internal logic error";
    b := allPoints(p, n-1);
    flatten for a from 0 to p-1 list (b/(b1 -> prepend(a, b1)))
    )

pointCount = method()
pointCount(RingElement, ZZ) := ZZ => (F, p) -> (
    -- F is a polynomial in 3 variables (FIXME: any number of variables)
    -- p is a prime number
    kk := ZZ/p;
    R := kk (monoid ring F);
    Fp := sub(F, R);
    allpts := allPoints(p, numgens ring F);
    allmaps := allpts/(pt -> map(kk, R, pt));
    ans1 := # for a in allpts list (phi := map(kk, R, a); if phi Fp == 0 then a else continue);
    --ans2 := # for x in (0,0,0)..(p-1,p-1,p-1) list if sub(Fp, matrix{{x}}) == 0 then x else continue;
    --if ans1 != ans2 then << "My previous code was incorrect" << endl;
    ans1
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

invariants3 = method()
invariants3 CYData := List => X -> (
    -- these are some invariants only involving the cubic form, not the c2 form...
    -- (but currently also involving h11 and h12.
    F := cubicForm X;
    h11 := hh^(1,1) X;
    h12 := hh^(1,2) X;
    RZ := ring F;
    if coefficientRing RZ =!= ZZ then 
        error "expected c2 and cubic form ring to be a polynomial ring over ZZ";
    RQ := QQ[gens RZ];
    facs := select((factors F)/toList/last, g -> support g != {});
    FQ := sub(F, RQ);
    d := dim saturate ideal jacobian FQ;
    singFZ := ideal gens gb saturate(ideal(F) + ideal jacobian F);
    sing := (F) -> ideal F + ideal jacobian F;
    linearcontent := (I) -> (
        if I == 0 then return 0;
        lins := select(I_*, f -> f != 0 and first degree f <= 1);
        if #lins == 0 then return 0;
        gcd for ell in lins list (trim content ell)_0
        );
    lincontent := linearcontent saturate sing F;
    singZ := flatten entries gens gb saturate(ideal(F) + ideal jacobian F);
    badp := select(singZ, a -> support leadTerm a === {});
    badp = if badp === {} then 0 else sub(first badp, ZZ);
    {h11, h12, badp, (trim content F)_0, #facs, d, lincontent}
    )

invariants4 = method()
invariants4 CYData := List => X -> (
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
    sing := (F) -> ideal F + ideal jacobian F;
    linearcontent := (I) -> (
        if I == 0 then return 0;
        lins := select(I_*, f -> f != 0 and first degree f <= 1);
        if #lins == 0 then return 0;
        gcd for ell in lins list (trim content ell)_0
        );
    lincontent := linearcontent saturate sing F;
    nc := # decompose ideal(LQ, FQ);
    singZ := flatten entries gens gb saturate(ideal(F) + ideal jacobian F);
    badp := select(singZ, a -> support leadTerm a === {});
    badp = if badp === {} then 0 else sub(first badp, ZZ);
    {h11, h12, badp, (trim content L)_0, (trim content F)_0, #facs, d, nc, lincontent}
    )

-- This one contains the best info we have to date (which isn't quite good enough).
polynomialContent = method()

-- Assumption: F is a polynomial over ZZ (not over a field).
-- Outout: the integer content of F.
polynomialContent RingElement := F -> (trim content F)_0

integerPart = method()
integerPart Ideal := (I) -> (
    -- expected: I is an ideal in a polynomial ring over ZZ.
    gs := select(flatten entries gens gb I, f -> support f === {});
    if #gs == 0 then 0 
    else if #gs == 1 then lift(gs_0, ZZ) 
    else error "internal error: somehow have two generators in ZZ in this GB"
    )

hessian = method()
hessian RingElement := F -> diff(vars ring F, diff(transpose vars ring F, F))

--factorShape = method()
-- This one is WRONG: lift(xxx, ZZ) could be positive or negative.  Those cannot be different.
-- factorShape RingElement := F -> (
--     facs := factors F;
--     sort for x in facs list if support x#1 == {} then 
--             {0, lift(x#1, ZZ)}
--         else
--             {first degree x#1, x#0}
--     )

factorShape = method()
factorShape RingElement := F -> (
    facs := factors F;
    sort for x in facs list if support x#1 == {} then 
            {0, abs lift(x#1, ZZ)}
        else
            {first degree x#1, x#0}
    )

-- factorShape RingElement := F -> (
--     facs := factors F;
--     con := trim content F;
--     con = con_0;
--     posfactors := sort for x in facs list if support x#1 == {} then 
--                      continue
--                    else
--                      {first degree x#1, x#0};
--     prepend({0, con},  posfactors)
--     )


invariantsAll = method()
invariantsAll(RingElement, RingElement, ZZ, ZZ) := (L, F, h11, h12) -> (
    RZ := ring L;
    if ring F =!= RZ then 
        error "expected same rings for c2 and cubic form";
    if coefficientRing RZ =!= ZZ then 
        error "expected c2 and cubic form ring to be a polynomial ring over ZZ";
    RQ := QQ[gens RZ];
    toQQ := F -> sub(F, vars RQ);
    sing := (cod, I) -> trim(I + minors(cod, jacobian I));
    linearcontent := (I) -> (
        if I == 0 then return 0;
        lins := select(I_*, f -> f != 0 and first degree f <= 1);
        if #lins == 0 then return 0;
        gcd for ell in lins list (trim content ell)_0
        );
    FQ := toQQ F;
    LQ := toQQ L;
    inv0 := polynomialContent L;
    inv1 := polynomialContent F;
    -- dimension and degree of each component of the singular loci over QQ.
    inv2 := sort for c in decompose sing_1 ideal FQ list {codim c, degree c};
    inv3 := sort for c in decompose sing_2 ideal(LQ, FQ) list {codim c, degree c};
    inv4 := betti res saturate sing_1 ideal FQ;
    -- inverse system of FQ
    inv5 := betti res inverseSystem FQ; -- not clear this one is worthwhile
    -- integer parts of singular loci.
    conductF := integerPart saturate sing_1 ideal F;
    conductLF := integerPart saturate sing_2 ideal(L,F);
    inv6 := conductF;
    inv7 := conductLF;
    inv8 := factorShape det hessian F;
    inv9 := linearcontent saturate sing_1 F;
    hashTable {"h11" => h11,
     "h12" => h12,
     "c(L)" => inv0, 
     "c(F)" => inv1, 
     "comps sing FQ" => inv2, 
     "comps sing LFQ" => inv3, 
     "bettti sing LFQ" => inv4,
     "betti inv F" => inv5,
     "conduct(F)" => inv6,
     "conduct(L,F)}" => inv7,
     "hessian shape" => inv8,
     "lincontent sing F" => inv9
     }
    )

invariantsAll CYData := X -> invariantsAll(c2Form X, cubicForm X, hh^(1,1) X, hh^(1,2) X)

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
            --Ms = Ms/(m -> lift(m, ZZ));
            Ms = for m in Ms list try lift(m, ZZ) else continue;
            Ms = select(Ms, m -> (d := det m; d === 1 or d === -1));
            isIsos := Ms/(m -> mapIsIsomorphism(m, Xj, Xi));
            if any(isIsos, x -> x == true) then (
                mi := position(isIsos, x -> x == true);
                << "found isomorphism from " << j << " to " << i << ": " << Ms#mi << endl;
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

-- Code to help determine equivalences between cubic forms.
-- 1. DONE Create ring T, TR, and (general) matrix A, and perhaps a phi: TR --> TR
--   corresponding to A.
-- 2. Given a list of linear forms that must map to other linear forms, find the constraint ideal
--   in T.
-- 3. Given points that map to points (in the contravariant map between affine spaces), 
--   find the constraint ideal.
-- 4. Given a list of linear forms that must match, up to sign, another set of linear forms.
--   
genericLinearMap = method(Options => {Variable => null})
genericLinearMap Ring := opts -> R -> (
    -- R should be a polynomial ring in n variables.
    n := numgens R;
    K := coefficientRing R;
    t := if opts.Variable === null then getSymbol "t" else opts.Variable;
    T := K[t_(1,1)..t_(n,n)];
    TR := T [gens R, Join => false];
    A := map(T^n,,transpose genericMatrix(T, T_0, n, n));
    phi := map(TR, TR, transpose A);
    (A, phi)
    )

TEST ///
  R = ZZ/101[a..d]
  (A, phi) = genericLinearMap R
  TR = target phi
  assert(source phi === TR)
  assert(ring A === coefficientRing TR)
  for i from 0 to 3 do 
    assert(phi TR_i == (A^{i} * (transpose vars TR))_(0,0))

  R = ZZ[a..d]
  (A, phi) = genericLinearMap R
  TR = target phi
  assert(source phi === TR)
  assert(ring A === coefficientRing TR)
  for i from 0 to 3 do 
    assert(phi TR_i == (A^{i} * (transpose vars TR))_(0,0))

  R = QQ[a..e]
  (A, phi) = genericLinearMap R
  TR = target phi
  assert(source phi === TR)
  assert(ring A === coefficientRing TR)
  for i from 0 to numgens R - 1 do 
    assert(phi TR_i == (A^{i} * (transpose vars TR))_(0,0))
///

linearEquationConstraints = method()
linearEquationConstraints(Matrix, RingMap, List, List) := Sequence => (A, phi, Ls, pts) -> (
    -- each entry of Ls is a list/sequence of length 2: {F, G}
    -- where F, G are polynomials in a ring R, (A, phi) are obtained from
    -- genericLinearMap.  We return the ideal of constraints in T
    -- for which phi(F) = G, for all pairs {F,G} in Ls.
    -- We also return A0, phi0 corresponding to these constraints.
    T := ring A;
    A0 := A;
    TR := target phi;
    I := sum for L in Ls list ideal sub(last coefficients(phi L_0 - L_1), T);
    if I != 0 then A0 = sub(A, T) % (trim I);
    J := sum for pq in pts list minors(2, 
        (A0 * (transpose matrix {pq#1})) | transpose matrix {pq#0});
    if J != 0 then A0 = A0 % (trim J);
    (A0, map(TR, TR, transpose A0))
    )

TEST ///
  R = QQ[a..d]
  (A, phi) = genericLinearMap R
  TR = target phi
  assert(source phi === TR)
  assert(ring A === coefficientRing TR)
  for i from 0 to 3 do 
    assert(phi TR_i == (A^{i} * (transpose vars TR))_(0,0))
    
  linearEquationConstraints(A, phi, {}, {
          {{1,0,0,0}, {1,1,0,0}},
          {{0,1,0,0}, {1,1,3,7}},
          {{0,0,1,0}, {5,6,-2,8}},
          {{0,0,0,1}, {0,1,0,0}}}
      )

  F1 = -2*a^3-6*a^2*b+6*b^2*c-12*a^2*d+12*a*b*d+12*b^2*d+36*b*c*d+30*a*d^2+60*b*d^2+54*c*d^2+76*d^3
  F7 = -2*a^3+6*a^2*b-6*a*b^2+2*b^3-6*a*c^2-4*c^3+6*a^2*d-6*a*d^2+2*d^3

  (A0, phi0) = linearEquationConstraints(A, phi, {
          {b+2*d, a-d},
          {b+3*d, a+c},
          {a+b+2*d, a-b},
          {F1, F7}
          }, {
          }
      )
  phi0 F1 == F7

  -- this one isn't correct yet.
  (A0, phi0) = linearEquationConstraints(A, phi, {
          {b+2*d, a-d},
          {b+3*d, a+c},
          {a+b+2*d, a-b}
          }, {
          }
      )
          {{0,0,1,0}, {1,1,-1,1}}

  trim ideal last coefficients(phi0 F1 - F7)


///

findMaps = (top1, top2, A, phi, RQ) -> (
    TR := target phi;
    n := numgens TR;
    T := coefficientRing TR;
    toTR := f -> sub(f, TR);
    toQQ := f -> sub(f, RQ);
    evalphi := (F,G) -> trim sub(ideal last coefficients((phi toTR F) - toTR G), T);
    (L1, F1, h11, h12) := toSequence top1;
    (L2, F2, l11, l12) := toSequence top2;
    RZ := ring L1;
    if RZ =!= ring F1 or RZ =!= ring L2 or RZ =!= ring F2 then error "expected polynomials over the same ring";
    if h11 != l11 or h12 != l12 then return null;
    I := (evalphi(L1, L2) + evalphi(F1, F2));
    if I == 1 then return null;
    -- first see if there is a unique solution.
    -- if codim I === n*n and degree I === 1 then (
    --     A0 := A % I;
    --     if support A0 === {} then (
    --         A0 = lift(A0, QQ);
    --         phi0 := map(RQ, RQ, transpose A0);
    --         if phi0 toQQ L1 != toQQ L2 or phi0 toQQ F1 != toQQ F2 then error "map is not correct!";
    --         return (A0, phi0)
    --         );
    --     );
    -- now let's look through all of the components for a smooth point.
    compsI := decompose I;
    As := for c in compsI list A % c;
    As = for a in As list try lift(a, ZZ) else continue; -- grab the ones that lift.
    As = select(As, a -> (d := det a; d == 1 or d == -1));
    if #As > 0 then (
        A0 := As#0;
        phi0 := map(RZ, RZ, transpose A0);
        if phi0 L1 != L2 or phi0 F1 != F2 then error "map is not correct!";
        (A0, phi0)
        )
    else (
        if any(compsI, c -> codim c < n*n or degree c =!= 1) then (
            << "warning: there might be a map in this case!" << endl;
            << netList compsI << endl;
            << "----------------------------------" << endl;
            compsI 
            )
        else null
        )
    )

partitionH113sByTopology = method()
partitionH113sByTopology(List, HashTable, Ring) := HashTable => (Ls, Ts, RQ) -> (
    -- Ls is a list of labels to separate.
    -- Ts is a hash table of label => {c2, cubicform, h11, h12}
    -- This function first separates these by the invariants: invariantsAll.
    -- The for each pair in each set, it attempts to find a map between them.
    -- output: a hashtable, keys are labels, values are lists of {label, matrix}
    (A, phi) := genericLinearMap RQ;
    if #Ls === 1 then return hashTable {Ls#0 => {}};
    H := partition(lab -> invariantsAll toSequence Ts#lab, Ls);
    distinctTops := new MutableHashTable; -- label => list of {label, matrix}, those with the same topology
    for i in Ls do (
        Ti := Ts#i;
        -- now we attempt to match this with each key of distinctTops
        << "trying " << i << endl;
        prev := keys distinctTops;
        isFound := false;
        for j in prev do (
            Tj := Ts#j;
            ans := findMaps(Tj, Ti, A, phi, RQ);
            --if j == (115,0) and i == (120,0) then error "debug me";
            if ans === null then (
                -- Ti is distinct from Tj
                )
            else if class first ans === Matrix then (
                -- we have a match!
                (A0, phi0) := ans;
                if all(flatten entries A0, a -> liftable(a, ZZ))
                then (
                    isFound = true;
                    distinctTops#j = append(distinctTops#j, {i, lift(A0, ZZ)});
                    break;
                    )
                )
            else (
                << (i,j) << " might be the same, might not" << endl;
                )
            );
        if not isFound then (
            distinctTops#i = {};
            << "found new top: " << i << endl;
            );
        );
    new HashTable from distinctTops
    )


-------------------------------------------------------------------------
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

newPackage(
    "IntegerEquivalences",
    Version => "0.1",
    Date => "13 Jan 2024",
    Headline => "finding invertible integral matrices preserving points, linear forms and ideals",
    Authors => {{ Name => "", Email => "", HomePage => ""}},
    AuxiliaryFiles => false,
    DebuggingMode => true
    )

export {
    "MatchingData",
    "equivalenceIdeal",
    "extendToMatrix", -- extendToMatrix(List of integers) ==> Matrix (over ZZ).
    "factorsByType",
    "idealsByBetti",
    "genericLinearMap", -- genericLinearMap(R).  Constructs two new rings, T, U, a matrix A over T nxn, n = numgens R, and phi = map(U, U, transpose A).
    "invertibleMatrixOverZZ",
    "matches",
    "selectLinear",
    "matchingData",
    "allSigns",
    "signedPermutations",
    "cartesian",
    "singularPoints",
    "singularPointMatches",
    "RowVector",
    "ColumnVector",
    "Unknown",
    "CONSISTENT",
    "INCONSISTENT",
    "INDETERMINATE",
    "SignedPermutations",
    "Permutations"
    }

importFrom_"LLLBases"{"gcdLLL"};

extendToMatrix = method()
extendToMatrix List := Matrix => L -> (
    if not all(L, a -> instance(a, ZZ))
    then error "expected a list of integers";
    (g, A) := gcdLLL L; -- coming from LLLBases.
    transpose A
    -- TODO: should this insure that the matrix has determinant 1 (not -1)?
    -- I don't really need that...
    )

genericLinearMap = method(Options => {Variable => null})
genericLinearMap Ring := Sequence => opts -> R -> (
    -- R should be a polynomial ring in n >= 1 variables.
    n := numgens R;
    if n == 0 then (
        kk := coefficientRing R; -- TODO: if none, this should give a better error message
        A := map(kk^0, kk^0, {});
        return (A, id_R);
        );
    K := coefficientRing R;
    t := if opts.Variable === null then getSymbol "t" else opts.Variable;
    T := K[t_(1,1)..t_(n,n)];
    U := T [gens R, Join => false];
    A = map(T^n,,transpose genericMatrix(T, T_0, n, n));
    phi := map(U, U, transpose A);
    (A, phi)
    )

-- MatchingData: is a list of elements each f the form
--  L => {M0, M1, ..., Ms}
-- where L is a list of:
--   RingElement: a polynomial in the original ring RZ or RQ.
--   Ideal: an ideal in the original ring.
--   Matrix: either a row vector or column vector, over ZZ or QQ (or RZ or RQ)
-- and each list Mi has the same length as L, and the same types of its elements.
MatchingData = new Type of List

matchingData = method()
matchingData List := LMs -> (
    if not all(LMs, x -> instance(x, Option) or (instance(x, List) and #x === 3))
    then error "expected each element to be Item => Item, or a list {type, ItemList, ItemList}";
    ans := new MatchingData from LMs;
    if not isWellDefined ans then error "expected matching data to match.  Set `debugLevel=1` to investigate";
    ans
    )

-- helper function for validMatchingItem.  
--   Input: either source or one of the targets of the matching data.
--   Output: a list of types.
itemType = L -> (
    -- L is a list or a single element of the following form
    -- returns a list of RingElement, Ideal, RowVector, ColumnVector.
    if not instance(L, List) then L = {L};
    for elem in L list (
        if instance(elem, RingElement) then RingElement
        else if instance(elem, Ideal) then Ideal
        else (
            if instance(elem, Matrix) then (
                if numrows elem === 1 then RowVector
                else if numcols elem === 1 then ColumnVector
                else Unknown
            ) else 
                Unknown
        ))
    )

-- helper function for (isWellDefined, MatchingData)
--   Input: one element (an Option, source and possible targets) of the matching data
--          which index this data sits at (used for error messages if debugLevel > 0)
--   Output: Boolean, whether this item is valid.
validMatchingItem = (LM,i) -> (
    -- TODO: if typ is null, then L should be a RingElement, Ideal or row or column Matrix, and M should be the same type
    if instance(LM, Option) then (
        if not (instance(LM#0, RingElement) or instance(LM#0, Ideal) or instance(LM#0, Matrix))
          or class LM#0 =!= class LM#1
          then (
              if debugLevel > 0 then 
                << "excepted each item of Option to be the same type: an Ideal, RingElement or row or column matrix" << endl;
                return false;
              );
          -- todo: make sure matrices are row or column matrices.
          return true;
        );
    (typ, L, M) := (LM#0, LM#1, LM#2);
    if typ =!= SignedPermutations and typ =!= Permutations then (
        if debugLevel > 0 then 
            << "expected first entry of list to be SignedPermutations or Permutations, instead, received "
                << toString typ << endl;
        return false;                
        );

    if not instance(L, List) or not instance(M, List) then (
        if debugLevel > 0 then 
          << "expected elements in source " << i << " to be lists" << endl;
          return false;
        );
    -- M should be a list, same length and types as L.
    Ltype := itemType L;
    Mtype := itemType M;
    if any(Ltype, x -> x === Unknown) then (
        if debugLevel > 0 then << "one element in source is unknown" << endl;
        return false;
        );
    if Mtype =!= Ltype then (
        if debugLevel > 0 then << "in source " << i << ", target doesn't match source = " << Ltype << " obtaining instead " << Mtype << endl;
        return false;
        );
    true
    )

isWellDefined MatchingData := Boolean => LMs -> (
    for i from 0 to #LMs-1 do (
        LM := LMs#i;
        if not instance(LM, Option) and not (instance(LM, List) and #LM === 3) then (
            if debugLevel > 0 then << "elements of list must be of the form {type, List, List}" << endl;
            return false;
            );
        if not validMatchingItem(LM,i) then return false;
        );
    true
    )

-- Used in creating MatchingData: use permutations or signedPermutations.
allSigns = method()
allSigns List := L -> (
    if #L <= 0 then return {{}};
    if #L == 1 then return {{L#0}, {-L#0}};
    flatten for q in allSigns(drop(L, 1)) list {prepend(L#0, q), prepend(-L#0, q)}
    )

signedPermutations = method()
signedPermutations List := List =>  L -> (
    flatten for p in permutations L list allSigns p
    )

cartesian = method()
cartesian List := (Ls) -> (
    -- cartesian product of Ls: one element from each
    -- so the result is a list of lists.
    if #Ls == 1 then return for p in Ls#0 list {p};
    Ls1 := cartesian drop(Ls,1);
    flatten for p in Ls#0 list for q in Ls1 list prepend(p, q)
    )

-- XXX this is being tested now.
matches = method()
matches MatchingData := List => (MD) -> (
    -- for each element of the MatchingData, we make the list of all possible targets
    src := flatten for elem in MD list
        if instance(elem, Option) then elem#0 else elem#1;
    -- now create targets.  This means making all permutations, signed or not, and taking cartesian products.
    targs := cartesian for elem in MD list (
        if instance(elem, Option) then {elem#1}
        else if elem#0 === SignedPermutations then signedPermutations elem#2
        else if elem#0 === Permutations then permutations elem#2
        );
    targs = targs/flatten;
    (src, targs)
    )

selectLinear = method()
selectLinear MatchingData := MatchingData => (MD) -> (
    -- ASSUMPTION: the base ring for Ideals and RingElement's is standard graded polynomial ring.
    -- we keep only the ring elements that are linear.
    -- for each ideal, we take only the linear elements.
    -- We keep all row and column matrices.
    matchingData for elem in MD list (
        if instance(elem, Option) then (
            if instance(elem#0, RingElement) then (
                if degree elem#0 === {1} then elem else continue
            ) else if instance(elem#0, Matrix) then elem
        ) else (
            if instance(elem#1#0, RingElement) then (
                if degree elem#1#0 === {1} then elem else continue
                )
            else if instance(elem#1#0, Ideal) then (
                elem1 := for x in elem#1 list
                    select(x_*, f -> degree f === {1});
                if #elem1#0 === 0 then continue;
                elem2 := for x in elem#2 list
                    select(x_*, f -> degree f === {1});
                {elem#0, elem1/ideal, elem2/ideal}
                )
            else elem
        ))
    )


invertibleMatrixOverZZ = method()
invertibleMatrixOverZZ(Matrix, Ideal) := Sequence => (A, J) -> (
    -- returns (determinacy, A0), or (INCONSISTENT, null) or ...?
    if J == 1 then 
        (INCONSISTENT, null)
    else (
        A0 := A % J;
        detA0 := (det A0) % J;
        suppA0 := support A0;
        if suppA0 === {} then (
            -- In this case we either have an integer matrix, or a rational matrix.
            if liftable(detA0, ZZ) and all(flatten entries A0, f -> liftable(f, ZZ)) then
                return (CONSISTENT, sub(A0, ZZ));
            return (INCONSISTENT, sub(A0, QQ));
            );
        jc := decompose J;
        if isPrime J then return (INDETERMINATE, jc);
        possibles := for j in jc list (
            ans := invertibleMatrixOverZZ(A0, j);
            if ans#0 == CONSISTENT then return ans else ans
            );
        if all(possibles, a -> a#0 === INCONSISTENT) then (
            ans := select(1, possibles, a -> instance(a#1, Matrix));
            if #ans > 0 then return ans#0 else return (INCONSISTENT, null);
            )
        else
            return (INDETERMINATE, jc);
        )
    )

equivalenceIdeal = method()
equivalenceIdeal(List, List, Ring, Sequence) := Ideal => (List1, List2, RQ, Aphi) -> (
    if itemType List1 =!= itemType List2 then 
        error("expected two lists to have the same list of types, they are: " 
            | toString itemType List1 | " and " | toString itemType List2);
    (A,phi) := Aphi;
    U := target phi;
    if U =!= source phi then error "expected ring map with same source and target";
    if #List1 == 0 then return ideal(0_U);
    B := coefficientRing U;
    toU := map(U, RQ, vars U);
    toB := map(B, U);
    List1 = List1/(I -> if ring I =!= RQ then toU (sub(I, RQ)) else toU I);
    List2 = List2/(I -> if ring I =!= RQ then toU (sub(I, RQ)) else toU I);
    -- List1 and List2 are lists with the same length, consisting of RingElement's, Ideal's, Matrices.
    -- List1 and List2 should each have RingElement's and Ideal's in the same spot.
    ids := for i from 0 to #List1-1 list (
        if instance(List1#i, RingElement) then (
            ideal toB (last coefficients(phi List1#i - List2#i))
            )
        else if instance(List1#i, Ideal) then (
            ideal toB (last coefficients((gens phi List1#i) % List2#i))
            )
        else if instance(List1#i, Matrix) then (
            rowvec := (numrows List1#i === 1);
            -- if rowvec is false, then this must be a column vector.
            if rowvec then
                ideal toB last coefficients sub(List2#i * (transpose A) - List1#i, U)
            else
                ideal toB last coefficients sub((transpose A) * List1#i - List2#i, U)
            )
        );
    sum ids
    )

tryEquivalences = method()
tryEquivalences(MatchingData, Ring, Sequence) := (MD, RQ, Aphi) -> (
    (A,phi) := Aphi;
    badJs := {};
    inconsistentMatrix := null;
    (src, tar) := matches MD;
    for i from 0 to #tar-1 do (
        << "doing " << i << endl;
        J := trim equivalenceIdeal(src, tar#i, RQ, Aphi);
        ans := invertibleMatrixOverZZ(A, J);
        if ans#0 == CONSISTENT then return ans;
        if ans#0 == INCONSISTENT and instance(ans#1, Matrix) then inconsistentMatrix = ans#1;
        if ans#0 == INDETERMINATE then (
            badJs = join(badJs, ans#1);
            );
        );
    if #badJs > 0 then return (INDETERMINATE, badJs);
    (INCONSISTENT, inconsistentMatrix)
    )

-------------------------------
-- Finding matching data of (L1,F1), (L2,F2)
-- using singular loci and hessian factorizations
-------------------------------
factors = method()
factors RingElement := (F) -> (
     facs := factor F;
     facs//toList/toList/reverse
     )

factorsByType = method()
factorsByType RingElement := HashTable => F -> (
    facs := factors F;
    faclist := for fx in facs list (fx#0, sum first exponents fx#1, fx#1);
    H := partition(x -> {x#0, x#1}, faclist);
    hashTable for k in keys H list k => for x in H#k list x_2
    )   

idealsByBetti = method()
idealsByBetti(List, List) := List => (J1s, J2s) -> (
    H1 := partition(J -> betti gens J, J1s);
    H2 := partition(J -> betti gens J, J2s);
    if sort keys H1 =!= sort keys H2 then return {};
    for k in sort keys H1 list (
        H1#k => permutations H2#k
        )
    )

-- TODO: place in M2 Core?
hessian = method()
hessian RingElement := F -> diff(vars ring F, diff(transpose vars ring F, F))

hessianMatches = method()
hessianMatches(Sequence, Sequence) := MatchingData => (LF1, LF2) -> (
    (L1, F1) := LF1;
    (L2, F2) := LF2;
    fac1 := factorsByType(det hessian F1);
    fac2 := factorsByType(det hessian F2);
    keys1 := sort select(keys fac1, k -> k =!= {1,0}); -- remove constant
    keys2 := sort select(keys fac2, k -> k =!= {1,0}); -- remove constant
    if keys1 =!= keys2 then return null; -- no matches.
    matchingData for k in keys1 list {SignedPermutations, fac1#k, fac2#k}
    )

-- Note: the returned value over a finite field is fine, and over QQ, all denominators and lcms have been cleared
singularPoints = method()
singularPoints RingElement := List => F -> (
    if not isHomogeneous F then error "expected homogeneous polynomial";
    R := ring F;
    kk := coefficientRing R;
    n := numgens R;
    singlocus := trim saturate(ideal F + ideal jacobian F);
    if singlocus == 1 then return {};
    comps := (decompose singlocus);
    comps0 := select(comps, c -> codim c == n-1 and degree c === 1); -- zero-dimensional rational points
    comps1 := select(comps, c -> not(codim c == n-1 and degree c === 1)); -- the rest
    comps1 = comps1/trim;
    if #comps1 != 0 then 
        << "other singular components: " << netList comps1 << endl;
    for c in comps0 list (
        pt := (vars R) % c;
        vs := support pt;
        if #vs != 1 then error "internal error: number of variables is not 1";
        pt = sub(pt, vs#0 => 1_kk);
        pt = flatten entries pt;
        if kk === QQ then (
            g := gcd pt;
            pt = 1/g * pt
            );
        transpose matrix{pt}
        )
    )
singularPointMatches = method()
singularPointMatches(RingElement, RingElement) := MatchingData => (F1, F2) -> (
    pts1 := singularPoints F1;
    pts2 := singularPoints F2;
    if #pts1 =!= #pts2 then return null; -- no matches.
    matchingData{{SignedPermutations, pts1, pts2}}
    )

singularMatches = method()
singularMatches(RingElement, RingElement) := MatchingData => (F1, F2) -> (
    sing1 := trim saturate(ideal F1 + ideal jacobian F1);
    sing2 := trim saturate(ideal F2 + ideal jacobian F2);
    if sing1 == 1 then return null;
    comps1 := (decompose sing1)/trim;
    comps2 := (decompose sing2)/trim;
    idealsByBetti(comps1, comps2)
    )

TEST ///
-- These 3 forms were generated from h11=5 database, with:
-- {(2249, 0), (2255, 0), (2270, 0)}
-- L1 = c2Form Xs#(2249,0)
-- L2 = c2Form Xs#(2255,0)
-- L3 = c2Form Xs#(2270,0)
-- F1 = cubicForm Xs#(2249,0)
-- F2 = cubicForm Xs#(2255,0)
-- F3 = cubicForm Xs#(2270,0)

-*
  restart
  needsPackage "IntegerEquivalences"
*-
-- This is supposed to test formation of matchingData
  debug needsPackage "IntegerEquivalences"
  RZ = ZZ[a,b,c,d,e]
  RQ = QQ (monoid RZ)
  (A,phi) = genericLinearMap RQ
  use RQ
  L1 = 12*a+10*b+34*c-4*d+4*e
  L2 = -4*a+12*b+18*c+28*d+4*e
  L3 = -4*a+12*b+26*c+10*d+20*e
  F1 = b^3-6*a^2*c-3*b^2*c+3*b*c^2+c^3+6*a^2*d+12*a*c*d+6*c^2*d-12*a*d^2-
      12*c*d^2+8*d^3-6*a^2*e+12*a*c*e+12*a*e^2-12*c*e^2-8*e^3
  F2 = 8*a^3-12*a^2*b+6*a*b^2-12*a^2*c+12*a*b*c-6*b^2*c+6*a*c^2-3*c^3-
      12*a^2*d+12*a*b*d-6*b^2*d+12*a*c*d-6*c^2*d+6*a*d^2-6*c*d^2-
      2*d^3+12*c*d*e+6*d^2*e-12*c*e^2-8*e^3
  F3 = 8*a^3-12*a^2*b+6*a*b^2-12*a^2*c+12*a*b*c-6*b^2*c+6*a*c^2-c^3+
      3*c^2*d-3*c*d^2+d^3-6*c*e^2-4*e^3
  FT1 = factorsByType det hessian F1
  FT2 = factorsByType det hessian F2
  FT3 = factorsByType det hessian F3

  srcs = flatten {FT1#{1,1}, FT1#{2,1}, L1}
  targs = flatten \ (cartesian {signedPermutations FT2#{1,1}, signedPermutations FT2#{2,1}, {L2}})
  for t in targs list (
    A % trim equivalenceIdeal(srcs, t, RQ, (A,phi))
    )

  srcs = flatten {FT1#{1,1}, FT1#{2,1}, L1}
  targs = flatten \ (cartesian {signedPermutations FT3#{1,1}, signedPermutations FT3#{2,1}, {L3}})
  for t in targs list (
    A0 = A % trim equivalenceIdeal(srcs, t, RQ, (A,phi));
    A0 = try lift(A0, QQ) else continue;
    A0 = try lift(A0, ZZ) else continue;
    f := map(RQ, RQ, (transpose A0)**QQ);
    {A0, f(F1) - F3, f(L1) - L3}
    )

  srcs = flatten {FT2#{1,1}, FT2#{2,1}, L2}
  targs = flatten \ (cartesian {signedPermutations FT3#{1,1}, signedPermutations FT3#{2,1}, {L3}})
  for t in targs list (
    A0 = A % trim equivalenceIdeal(srcs, t, RQ, (A,phi));
    A0 = try lift(A0, QQ) else continue;
    A0 = try lift(A0, ZZ) else continue;
    f := map(RQ, RQ, (transpose A0)**QQ);
    {A0, f(F2) - F3, f(L2) - L3}
    )

  singularPointMatches(F1, F2)
  singularMatches(F1, F2)
  equivalenceIdeal(first first oo, first last first oo, RQ, (A,phi))
  equivalenceIdeal({transpose matrix{{2,0,0,1,1}}}, {transpose matrix{{0,0,-2,2,1}}}, RQ, (A,phi))
  equivalenceIdeal({transpose matrix{{2,0,0,1,1}}}, {transpose matrix{{0,0,2,-2,-1}}}, RQ, (A,phi))
///

-- NOT AT ALL FUNCTIONAL YET (14 Jan 2024)
singularAndHessianMatches = method()
singularAndHessianMatches(Sequence, Sequence, Sequence) := MatchingData => (LF1, LF2, Aphi) -> (
    (L1, F1) := LF1;
    (L2, F2) := LF2;
    (A,phi) := Aphi;
    RQ := QQ (monoid ring L1); -- TODO: allow ZZ/p as well...
    FQ1 := sub(F1, RQ);
    FQ2 := sub(F2, RQ);
    fac1 := factorsByType(det hessian F1);
    fac2 := factorsByType(det hessian F2);
    keyset1 := sort for k in sort keys fac1 list if k === {1,0} then continue else k;
    keyset2 := sort for k in sort keys fac2 list if k === {1,0} then continue else k;
    if keyset1 =!= keyset2 then error "expected same factors and their degrees";
    list1 := flatten for k in keyset1 list fac1#k;
    list2 := for k in keyset2 list fac2#k;
    -- list2perms := cartesian for fs in list2 list signedPermutations fs;
    -- hessPart := (list1, list2perms);
    sing1 := trim saturate(ideal FQ1 + ideal jacobian FQ1);
    sing2 := trim saturate(ideal FQ2 + ideal jacobian FQ2);
    if sing1 == 1 then return {};
    comps1 := (decompose sing1)/trim;
    comps2 := (decompose sing2)/trim;
    list2perms := null;
    (list1, list2perms) = idealsByBetti(comps1, comps2);
    singPart := (join({L1, F1, sing1}, list1), cartesian{{{L2, F2, sing2}}, list2perms});
    --secondPart := flatten for a in singPart#1 list for b in hessPart#1 list (a | b);
    --((singPart#0 | hessPart#0), secondPart)
    )
-- XXX


TEST ///
-*
  restart
  needsPackage "IntegerEquivalences"
*-
  A = extendToMatrix{10,15,6}
  assert(A * transpose matrix{{10,15,6}} == transpose matrix{{0,0,1}})
  assert(abs det A == 1)

  A = extendToMatrix{10,15,1}
  assert(A * transpose matrix{{10,15,1}} == transpose matrix{{0,0,1}})
  assert(abs det A == 1)

  A = extendToMatrix{10,15,0}
  assert(A * transpose matrix{{10,15,0}} == transpose matrix{{0,0,5}}) -- 5 is gcd!
  assert(abs det A == 1)
///

TEST ///
-*
  restart
  needsPackage "IntegerEquivalences"
*-
  -- trivial case
  R = ZZ[]
  (A, phi) = genericLinearMap R
  assert(numRows A == 0 and numcols A == 0)
  assert(ring A === coefficientRing target phi)
  assert(source phi === target phi)

  R = ZZ[a..d]
  (A, phi) = genericLinearMap R
  assert(numRows A == 4 and numcols A == 4)
  assert(ring A === coefficientRing target phi)
  assert(source phi === target phi)

  (A, phi) = genericLinearMap(R, Variable => symbol s)
  assert(numRows A == 4 and numcols A == 4)
  assert(ring A === coefficientRing target phi)
  assert(source phi === target phi)
  use ring A
  assert(A_(0,0) == s_(1,1))  
  
  RQ = QQ (monoid R)
  (A, phi) = genericLinearMap RQ
  assert(numRows A == 4 and numcols A == 4)
  assert(ring A === coefficientRing target phi)
  assert(source phi === target phi)
  
  R = ZZ[a]  
  (A, phi) = genericLinearMap R
  assert(numRows A == 1 and numcols A == 1)
  assert(ring A === coefficientRing target phi)
  assert(source phi === target phi)

  R = ZZ/101[a..d]
  (A, phi) = genericLinearMap R
  U = target phi
  assert(source phi === U)
  assert(ring A === coefficientRing U)
  for i from 0 to 3 do 
    assert(phi U_i == (A^{i} * (transpose vars U))_(0,0))

  R = ZZ[a..d]
  (A, phi) = genericLinearMap R
  U = target phi
  assert(source phi === U)
  assert(ring A === coefficientRing U)
  for i from 0 to 3 do 
    assert(phi U_i == (A^{i} * (transpose vars U))_(0,0))

  R = QQ[a..e]
  (A, phi) = genericLinearMap R
  U = target phi
  assert(source phi === U)
  assert(ring A === coefficientRing U)
  for i from 0 to numgens R - 1 do 
    assert(phi U_i == (A^{i} * (transpose vars U))_(0,0))
///

TEST ///
-*
  restart
  needsPackage "IntegerEquivalences"
*-
  RZ = ZZ[a,b,c]
  F1 = 3*a^2*b-3*a*b^2+b^3+6*a^2*c-6*a*c^2+2*c^3
  L1 = 24*a+10*b+8*c
  F2 = F1
  L2 = L1
  M1 = a
  M2 = b
  M3 = c
  N1 = a+b
  N2 = a+c
  N3 = 2*b+c

  md = matchingData {
    F1 => F2, 
    L1 => L2,
    {Permutations, {ideal M1, ideal M2}, {ideal N1, ideal N2}},
    {SignedPermutations, {M1,M2,M3}, {N1,N2,N3}},
    matrix{{1,2,3}} => matrix{{1,-2,1}}
    }
  assert isWellDefined md
  matches md
  netList first matches md
  netList last matches md

  md = matchingData {
    F1 => F2, 
    L1 => L2,
    {Permutations, {ideal(M1, M2^2), ideal(M2,M3^2)}, {ideal(N1,N2^2), ideal(N2,N3^2)}},
    {SignedPermutations, {M1,M2,M3}, {N1,N2,N3}},
    matrix{{1,2,3}} => matrix{{1,-2,1}}
    }
  assert isWellDefined md
  matches md
  netList first matches md
  netList last matches md

  md1 = selectLinear md
  
  assert try (matchingData {
    {F1, F2},
    {Permutations, {L1}, {L2}},
    {SignedPermutations, {M1,M2,M3}, {N1,N2,N3}},
    {matrix{{1,2,3}}, matrix{{1,-2,1}}}
    }; false
  ) else true

  assert try (
      matchingData {
          F1 => F2,
          L1 => L2,
          {M1,M2,M3} => {{N1,N2,N3,N1}},
          matrix{{1,2,3}} => matrix{{1,-2,1}}
          };
      false
  ) else true;
///

TEST ///
-*
  restart
  needsPackage "IntegerEquivalences"
*-
  assert(# signedPermutations{1,2,3} === 48)
  assert(# permutations{1,2,3} === 6)
  assert(cartesian{{1,2}} == {{1}, {2}})
  assert(cartesian{{1,2}, {3,4}} === {{1, 3}, {1, 4}, {2, 3}, {2, 4}})
  assert(cartesian{{1,2},{3,4},{5,6}} === {{1, 3, 5}, {1, 3, 6}, {1, 4, 5}, {1, 4, 6}, {2, 3, 5}, {2, 3, 6}, {2, 4, 5}, {2, 4, 6}})
  assert(cartesian{{1,2},{3},{5,6}} === {{1, 3, 5}, {1, 3, 6}, {2, 3, 5}, {2, 3, 6}})
  assert(cartesian{{1,2},{},{5,6}} === {})
///

TEST ///
-*
  restart
  needsPackage "IntegerEquivalences"
*-
  RZ = ZZ[a,b,c]
  RQ = QQ (monoid RZ)
  (A,phi) = genericLinearMap RQ

  -- Here is how we construct these polynomials (using h11=3 database).  
  --(L1, F1) = (c2Form Xs#(26,0), cubicForm Xs#(26,0))
  --(L2, F2) = (c2Form Xs#(37,0), cubicForm Xs#(37,0))
  
  (L1, F1) = (10*a+28*b+26*c,a^3-3*a^2*b-3*a*b^2-2*b^3-3*a^2*c+6*a*b*c+6*b^2*c+3*a*c^2+6*b*c^2-c^3)
  (L2, F2) = (16*a+10*b+26*c,-2*a^3-3*a^2*b-3*a*b^2+b^3+6*a*b*c-3*b^2*c+6*a*c^2+3*b*c^2-c^3)
  
  MD = matchingData {
      L1 => L2,
      F1 => F2
      }
  (src, tars) = matches MD
  J = equivalenceIdeal(src, tars#0, RQ, (A,phi))
  A % J
  assert(first invertibleMatrixOverZZ(A, J) == INCONSISTENT)
///

TEST ///
-*
  restart
  needsPackage "IntegerEquivalences"
*-
  RZ = ZZ[a,b,c]
  RQ = QQ (monoid RZ)
  (A,phi) = genericLinearMap RQ

  -- Here is how we construct these polynomials (using h11=3 database).  
  --(L1, F1) = (c2Form Xs#(26,0), cubicForm Xs#(26,0))
  --(L2, F2) = (c2Form Xs#(37,0), cubicForm Xs#(37,0))

  (L1, F1) = (8*a-4*b+36*c,2*a^3-3*a^2*b-3*a*b^2+8*b^3-6*a^2*c+6*a*b*c-6*b^2*c+6*a*c^2)
  (L2, F2) = (-4*a+8*b+36*c,8*a^3-3*a^2*b-3*a*b^2+2*b^3-6*a^2*c+6*a*b*c-6*b^2*c+6*b*c^2)
  --factorsByType det hessian F1
  
  MD = matchingData {
      L1 => L2,
      F1 => F2
      }
  (src, tars) = matches MD
  J = equivalenceIdeal(src, tars#0, RQ, (A,phi))
  A % J
  assert(first invertibleMatrixOverZZ(A, J) == CONSISTENT)
///

TEST ///
-*
  restart
  needsPackage "IntegerEquivalences"
*-
///

end--

-* Documentation section *-
beginDocumentation()

doc ///
Key
  IntegerEquivalences
Headline
Description
  Text
  Tree
  Example
  CannedExample
Acknowledgement
Contributors
References
Caveat
SeeAlso
Subnodes
///

doc ///
Key
Headline
Usage
Inputs
Outputs
Consequences
  Item
Description
  Text
  Example
  CannedExample
  Code
  Pre
ExampleFiles
Contributors
References
Caveat
SeeAlso
///

-* Test section *-
TEST /// -* [insert short title for this test] *-
-- test code and assertions here
-- may have as many TEST sections as needed
///

end--

-* Development section *-
restart
needsPackage "IntegerEquivalences"
check "IntegerEquivalences"

uninstallPackage "IntegerEquivalences"
restart
installPackage "IntegerEquivalences"
viewHelp "IntegerEquivalences"

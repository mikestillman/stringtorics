-- In this Macaulay2 file, we (try to) determine the inequivalent CY3's of h11=3
-- First, we find all possible topologies for toric hypersurfaces.
-- Then we will also consider hopefully all non-toric flops as well.

---------------------------------------------------------
-- Step 1.  Generate the database file for h11=3 CY3's --
---------------------------------------------------------
-- This we do only once.

  restart
  needsPackage "StringTorics"
  topes = kreuzerSkarke(3, Limit => 1000);
  assert(#topes == 244)
  elapsedTime createPolytopeDatabase("foo-cys-ntfe-h11-3.dbm", topes) -- 83 seconds

  -- Now let's add in all the CY's total, including all triangulations
  -- which are distinct when restricted to 2-faces (NTFE => true says don't use all triangulations).
  Qs = readCYPolytopes "foo-cys-ntfe-h11-3.dbm";
  elapsedTime for Q in values Qs do addToCYDatabase("foo-cys-ntfe-h11-3.dbm", Q, NTFE => true); -- 11 seconds

  -- this creates a database whose keys are integers 0, 1, ..., 243
  -- (one for each polytope in the KS database, in the same order, and
  -- (polytope#, triangulation#).  See step 2 for grabbing data from
  -- there.

----------------------------------------------------------------------------
-- Step 2.  Regenerate from data base all polytopes (Qs), all CY3's (Xs) ---
----------------------------------------------------------------------------
  -- First we grab all of the polytopes.  This returns a hash table: keys are integers 0..243, 
  -- and the corresponding value is the corresponding "CYPolytopeData"

  restart
  debug needsPackage "StringTorics"

  Qs = readCYPolytopes "foo-cys-ntfe-h11-3.dbm";
  sort keys Qs -- all 0, 1, ..., 243.

  RZ = ZZ[a,b,c]
  Xs = readCYs("foo-cys-ntfe-h11-3.dbm", Qs, Ring => RZ);
  sort keys Xs
  
  # keys Xs

  -- some simple checking.
  for k in sort keys Qs do assert instance(Qs#k, CYPolytopeData)
  for k in sort keys Xs do assert instance(Xs#k, CYData)

  -- example of use of Q:
  Q = Qs#100
    degrees Q
    isFavorable Q
    hh^(1,1) Q
    hh^(1,2) Q
    label Q
    netList annotatedFaces Q
    peek Q.cache -- this data is stored in Q, so it doesn't need to be recomputed.
  X = Xs#(100,0)
    intersectionNumbers X
    c2 X
    cubicForm X
    c2Form X
    label X
    -- the toric varity for X can be constructed:
    max X
    rays X
    -- playing well with NormalToricVarieties, Schubert2.
    X' = normalToricVariety X
    Xa = abstractVariety X
    peek X.cache

----------------------------------------------------------------------------
-- Step 3.  Enumerate all possible topologies ------------------------------
----------------------------------------------------------------------------
  #(keys Qs) == 244
  #(keys Xs) == 306 -- this is the number of possibly different CY3's which are toric, h11=3.
  (keys Xs)/(lab -> topologicalData(Xs#lab))//unique;
  #oo == 292 -- so only 14 "exact" duplicates
  partition(lab -> topologicalData(Xs#lab), keys Xs);
  
  -- currently, we can only compute GV invariants for those whose normal toric variety has non-torsion class group,
  -- and also must be favorable.
  -- I hope to remove these restrictions soon!
  
  nontorsionfrees = for lab in sort keys Xs list (
      X := Xs#lab;
      V := normalToricVariety(rays X, max X); 
      if classGroup V != ZZ^3 then lab else continue
      )
  nontorsionfrees == {(0, 0), (9, 0), (10, 0), (55, 0), (62, 0), (232, 0)}
  torsionfrees = for lab in sort keys Xs list (
      X := Xs#lab;
      V := normalToricVariety(rays X, max X); 
      if classGroup V == ZZ^3 then lab else continue
      );
  #torsionfrees == 300

  invariantsAll Xs#(100,0)
  -- so FOR NOW, we ignore the 6 X's which are in the nontorsionfree list.
  
    elapsedTime H = partition(lab -> invariantsAll Xs#lab, torsionfrees); -- 19 seconds
    #(keys H) == 170 -- used to be == 195 -- this is the minimum number of different topologies.
    (keys H)/(k -> #H#k)//tally
    elapsedTime INV = for k in keys H list k => partitionByTopology(H#k, Xs, 15); -- 34 seconds
    INV/(x -> #x#1)//tally
    INVGV = INV;
    INVGV = (INVGV/last)
    -- so the maximum number of distinct topologies is 159 + 18 + 6 = 183
    -- used to be (before fixing at least 3 bugs!): so the maximum number of distinct topologies is 186 + 16 + 3 == 205

  REPS = for k in INV list (
      H := last k; -- a hash table
      K := keys H;
      if #K  == 1 then continue
      else sort K)
  netList REPS


  ALLREPS = for k in INV list (
      H := last k; -- a hash table
      K := keys H;
      sort K)
  netList sort (ALLREPS)
  DIFFTOPS = sort flatten ALLREPS -- these are all, I believe, distinct topologies (mod worries about torsion)

----------------------------------------------------------------------------
-- Step 3A  Enumerate all possible topologies, not using GV's --------------
----------------------------------------------------------------------------
-- The idea here: first separate all topologies via invariantsAll.
-- Then for each of those, use the ansatz phi(L1)=L2, phi(F1)=F2 to
-- find an isomorphism.  The key thing used here: for all or almost
-- all of the cases here, this ansatz produces a unique map.  In this
-- case, if it is integral, we have an isomorphism, if not, we don't.
  nontorsionfrees = for lab in sort keys Xs list (
      X := Xs#lab;
      V := normalToricVariety(rays X, max X); 
      if classGroup V != ZZ^3 then lab else continue
      )
  nontorsionfrees == {(0, 0), (9, 0), (10, 0), (55, 0), (62, 0), (232, 0)}
  torsionfrees = for lab in sort keys Xs list (
      X := Xs#lab;
      V := normalToricVariety(rays X, max X); 
      if classGroup V == ZZ^3 then lab else continue
      );
  #torsionfrees == 300
  favorables = for lab in sort keys Xs list if first lab === 232 then continue else lab
  
Ts = hashTable for lab in favorables list (
    X := Xs#lab;
    lab => {c2Form X, cubicForm X, hh^(1,1) X, hh^(1,2) X}
    )
#Ts == 305

elapsedTime H = partition(lab -> invariantsAll toSequence (Ts#lab), favorables); -- 19 sec
#keys H == 170
-- now for each set, we separate via finding isomorphisms...
RQ = QQ[a,b,c]
elapsedTime INV = for k in keys H list k => partitionH113sByTopology(H#k, Ts, RQ); -- 40 sec now 170 seconds...  Fix that!
    INV/(x -> #x#1)//tally
    INVMES = INV/last; -- 186 different topologies (not counting the non-favorable (232,0))
  REPS = for k in INV list (
      H := last k; -- a hash table
      K := keys H;
      if #K  == 1 then continue
      else sort K)
  netList REPS


  ALLREPS = for k in INV list (
      H := last k; -- a hash table
      K := keys H;
      sort K)
  netList sort (ALLREPS)

  DIFFTOPS = sort flatten ALLREPS -- these are all, I believe, distinct topologies (mod worries about torsion)

  -- Now we do the above, with the DIFFTOPS and all the ones from Richard's h11=3 database.
  load "~/src/stringtorics/m2-examples/richard-example/richard-db.m2"
  RICH = flatten flatten for h21 in findH21s DIR list 
    for pol in findPolys(DIR, h21) list
      for cy in findCYs(DIR, h21, pol) list 
        (h21, pol, cy) => toList join(getTopology(h21, pol, cy, RZ), {3, h21});

  allTs = hashTable join(
      for lab in DIFFTOPS list lab => Ts#lab,
      for x in RICH list x
      );
  #keys allTs == 969
  assert(#keys allTs == #DIFFTOPS + #RICH)

  elapsedTime H = partition(lab -> invariantsAll toSequence (allTs#lab), keys allTs); 
  #keys H == 279
  -- now for each set, we separate via finding isomorphisms...
  RQ = QQ[a,b,c]
  elapsedTime INV = for k in keys H list k => partitionH113sByTopology(H#k, allTs, RQ); --

  ALLREPS = for k in INV list (
      H := last k; -- a hash table
      K := keys H;
      sort K)
  netList sort (ALLREPS)

  DIFFTOPS = sort flatten ALLREPS -- these are all, I believe, distinct topologies (mod worries about torsion)
  #DIFFTOP == 301
    
-- testing:
-- these all have the same invariants.  How many are isomorphic?
Ls = {(138, 0), (141, 0), (142, 0), (143, 0), (143, 1), (143, 2), (149, 0), (149, 1), (149, 2), (149, 3), (149, 4), (149, 5), (152, 0), (152, 1)}
Ls = {(115,0), (120,0)}
P1 = partitionByTopology(Ls, Xs, 15)
P2 = partitionH113sByTopology(Ls, Ts, RQ)

  X1 = Xs#(115,0)
  X2 = Xs#(120,0)
  L1 = Ts#(115,0)#0
  F1 = Ts#(115,0)#1
  L2 = Ts#(120,0)#0
  F2 = Ts#(120,0)#1

  X1 = Xs#(138,0)
  X2 = Xs#(143,0)
  L1 = Ts#(138,0)#0
  F1 = Ts#(138,0)#1
  L2 = Ts#(143,0)#0
  F2 = Ts#(143,0)#1
  cubicForm X1 == F1
  c2Form X1 == L1

  A0 = P1#(138,0)#0#1
  A1 = P2#(138,0)#0#1
  phi1 = map(RZ, RZ, A0)
  phi1 L2 - L1
  phi1 = map(RZ, RZ, transpose A1)
  phi1 L1 - L2
  phi1 F1 - F2

  gv1 = partitionGVConeByGV X1
  gv2 = partitionGVConeByGV X2
  Ms = findLinearMaps(gv1, gv2)
  Ms = for m in Ms list try lift(m, ZZ) else false
  Ms = select(Ms, m -> (d := det m; d === 1 or d === -1))
  isIsos = Ms/(m -> mapIsIsomorphism(m, X1, X2))
  Ms_1  
  phi1 = map(RZ, RZ, Ms_1)
  phi1 F1 - F2
  phi1 L1 - L2
  det Ms_1
  Ms_1 == A0

  mapIsIsomorphism(sub(A0s_0, QQ), X2, X1)
  mapIsIsomorphism(sub(A0s_1, QQ), X2, X1)

  Ms = findLinearMaps(gv2, gv1)
  Ms = Ms/(m -> lift(m, ZZ))
  Ms = select(Ms, m -> (d := det m; d === 1 or d === -1));
  isIsos := Ms/(m -> mapIsIsomorphism(m, X2, X1))

            Ms = for m in Ms list try lift(m, ZZ) else false
            isIsos = Ms/(m -> mapIsIsomorphism(m, X1, X2))
  Ms_1  
  
Ts#(142,0)
Ts#(152,0)
(A, phi) = genericLinearMap RQ
findMaps(Ts#(142,0), Ts#(152,0), A, phi, RQ)
kH = keys H;
H#(kH#0)
H#(kH#1)
H#(kH#2)
(A, phi) = genericLinearMap RQ
partitionH113sByTopology(H#(kH#2), Ts, RQ)
findMaps(Ts#(227,0), Ts#(229,0), A, phi, RQ)

----------------------------------------------------------------------------
-- Step 4.  Analyze these ones that are possibly still equivalent-----------
----------------------------------------------------------------------------


  X1 = Xs#(141,0)
  X2 = Xs#(142,0)
  X3 = Xs#(143,0)
  F1 = cubicForm X1
  F2 = cubicForm X2
  F3 = cubicForm X3
  L1 = c2Form X1
  L2 = c2Form X2
  L3 = c2Form X3

  factor det hessian F1
  factor det hessian F2
  factor det hessian F3
  
  R = QQ[a,b,c]
  (A, phi) = genericLinearMap R
  TR = target phi
  T = coefficientRing TR
  phi
  A
  assert(source phi === TR)
  assert(ring A === coefficientRing TR)

   use TR
    (A0, phi0) = linearEquationConstraints(A, phi, {
            {sub(L1,TR), sub(L2,TR)},
            {a+c, a+b+3*c},
            {a, b+3*c}}, {})

    (A0, phi0) = linearEquationConstraints(A, phi, {
            {sub(L1,TR), sub(L2,TR)},
            {a+c, -(a+b+3*c)},
            {a, b+3*c}}, {})

    (A0, phi0) = linearEquationConstraints(A, phi, {
            {sub(L1,TR), sub(L2,TR)},
            {a+c, (a+b+3*c)},
            {a, -(b+3*c)}}, {})

    (A0, phi0) = linearEquationConstraints(A, phi, {
            {sub(L1,TR), sub(L2,TR)},
            {a+c, -(a+b+3*c)},
            {a, -(b+3*c)}}, {})

  phi0 sub(F1, TR) == sub(F2, TR)
  phi0 sub(L1, TR) == sub(L2, TR)
A0


   use TR
    (A0, phi0) = linearEquationConstraints(A, phi, {
            {sub(L1,TR), sub(L3,TR)},
            {a+c, (b+c)},
            {a, (a+2*b+c)}}, {})

    (A0, phi0) = linearEquationConstraints(A, phi, {
            {sub(L1,TR), sub(L3,TR)},
            {a+c, (b+c)},
            {a, -(a+2*b+c)}}, {})

  phi0 sub(F1, TR) == sub(F3, TR)
  phi0 sub(L1, TR) == sub(L3, TR)

    (A0, phi0) = linearEquationConstraints(A, phi, {
            {a+c, (b+c)},
            {a, -(a+2*b+c)}}, {})

  I = trim ideal last coefficients (phi0 sub(F1, TR) - sub(F3, TR))
  A0 % sub(I, T)


  ----------------------------------
  -- checking ones that invariants can't tell apart, but GV thinks are different.
  RQ = QQ[a,b,c]
  toQQ = f -> sub(f, RQ)
  (A, phi) = genericLinearMap R
  TR = target phi
  T = coefficientRing TR

  X1 = Xs#(34,0)
  X2 = Xs#(37,0) -- these are different, equiv over QQ.

  X1 = Xs#(183,0)
  X2 = Xs#(195,0) -- these are different, equiv over QQ.

  -- these two have a 1D family of matrices over QQ that map F1-->F2, L1-->L2.
  -- two (over QQ) have det 1 or -1.  No integral matrices do it.
  X1 = Xs#(63,1)
  X2 = Xs#(72,0) -- these are different, equiv over QQ

  X1 = Xs#(138,0)
  X2 = Xs#(152,0) -- these are different, equiv over QQ

  X1 = Xs#(227,0)
  X2 = Xs#(229,0) -- these are different, equiv over QQ

  X1 = Xs#(63,0)
  X2 = Xs#(72,0) -- these are different, not equiv over QQ (well, that maps cubic and c2)

  -- these two have a 1D family of matrices over QQ that map F1-->F2, L1-->L2.
  -- two (over QQ) have det 1 or -1.  No integral matrices do it.
  X1 = Xs#(87,0)
  X2 = Xs#(88,1) -- these are different, equiv over QQ

  X1 = Xs#(141,0)
  X2 = Xs#(142,0) -- different, equiv over QQ.

  X1 = Xs#(141,0)
  X2 = Xs#(143,0) -- different, equiv over QQ.

  X1 = Xs#(142,0)
  X2 = Xs#(143,0) -- different, equiv over QQ.
  
  -- having tried all these pairs, we now have 205 different topologies at h11=3.
  -- but we must consider the 6 torsion examples too.

  F1 = cubicForm X1
  F2 = cubicForm X2
  L1 = c2Form X1
  L2 = c2Form X2
  factor det hessian F1
  factor det hessian F2

  I1 = sub(trim ideal last coefficients(phi sub(L1, TR) - sub(L2, TR)), coefficientRing TR)
  I2 = sub(trim ideal last coefficients(phi sub(F1, TR) - sub(F2, TR)), coefficientRing TR)
  I = trim sub(I1 + I2, T)
  gbI = ideal gens gb I;
  A0 = A % gbI
  phi0 = map(RQ, RQ, transpose lift(A0, QQ))
  phi0 toQQ F1 - toQQ F2
  phi0 toQQ L1 - toQQ L2

  A0 = A0 % (det A0 - 1) -- over QQ
  A0 = A0 % (det A0 + 1) -- over QQ
  
-----------------------------------
-- Step 5. Consider all the 205, together with the reminaing 6.
-- Let's analyze how many different cubics there are...

  cubics = hashTable(DIFFTOPS/(lab -> lab => cubicForm Xs#lab))

  -- let's separate out the ones that factor.
  partition(lab -> drop(factorShape cubics#lab, 1), keys cubics)
  partition(lab -> factorShape cubics#lab, keys cubics)

  GROUPS = partition(lab -> (
          F := cubics#lab;
          singF := saturate ideal jacobian toQQ F;
          compsF := decompose singF;
          {
              drop(factorShape F, 1),
              dim singF,
              degree singF,
              compsF/((i -> degree i)),
              drop(factorShape det hessian F, 1)}
          ),
          keys cubics
          )


  netList for k in sort keys GROUPS list append(k, #GROUPS#k)

  GROUPSZZ = partition(lab -> (
          F := cubics#lab;
          singF := saturate (ideal F + ideal jacobian F);
          cond := select(flatten entries gens gb singF, f -> support f === {});
          cond = trim sub(ideal cond, ZZ);
          cond = if numgens cond === 0 then 0 else cond_0;
          compsF := decompose toQQ singF;
          {
              cond,
              factorShape F,
              dim toQQ singF,
              degree toQQ singF,
              compsF/((i -> degree i)),
              factorShape det hessian F}
          ),
          keys cubics
          )
  netList for k in sort keys GROUPSZZ list append(k, #GROUPSZZ#k)

  select(values GROUPSZZ, v -> #v == 9)
  
  SMOOTHS = select(keys cubics, lab -> (
          F := cubics#lab;
          singF := saturate ideal jacobian toQQ F;
          singF == 1))

  partition(lab -> (
          F := cubics#lab;
          singF := saturate ideal jacobian toQQ F;
          compsF := decompose singF;
          {drop(factorShape F, 1),
              dim singF,
              degree singF,
              compsF/((i -> degree i))
              }
          ),
          keys cubics
          )

-- Compare j-invariants of all smooth examples
  needsPackage "EllipticCurves"
  toWeierstrass(RingElement, List) := (F, pt) -> (
      R := ring F;
      kk := coefficientRing R;
      if numgens R =!= 3 or #pt != 3 then error "expected ring in 3 variables and point in 3-space";
      mons3 := basis(3, R);
      coeffsF := flatten entries last coefficients(F, Monomials => mons3);
      toWeierstrass(coeffsF, pt, kk)
      )
  
  findPoint = (F, bound) -> (
      for p in {-bound,-bound,-bound}..{bound,bound,bound} do (
          if p === {0,0,0} then continue;
          if sub(F, matrix{p}) == 0 then return p
          );
      null
      )
  NOPOINTYET = SMOOTHS
  NOPOINTYET1 = {}

  FOUND2 = for lab in NOPOINTYET list (
      pt := findPoint(cubics#lab, 7);
      if pt =!= null then (lab => pt) 
      else (
          NOPOINTYET1 = append(NOPOINTYET1, lab);
          continue)
      )
  NOPOINTYET = NOPOINTYET1
  #FOUND2 == 99
  #NOPOINTYET == 18

  FOUND = FOUND2 

  JINV = hashTable for x in FOUND list x#0 => jInvariant toWeierstrass(toQQ cubics#(x#0), x#1)  
  partition(lab -> JINV#lab, keys JINV)
  netList for k in sort keys oo list {k, sort oo#k}

  -- Let's create from FOUND a file so that Richard can determine j-invariants of these curves
  -- This is the line I used to create (via cut and paste) smooth-cubics-for-richard.txt.
  for x in FOUND/first//sort do print (toString cubics#x | " # "|toString x)

  for x in sort NOPOINTYET do print (toString cubics#x | " # "|toString x)
  pts = hashTable FOUND  
  for x in FOUND/first//sort list {x, jInvariant toWeierstrass(toQQ cubics#x, pts#x)}
  
  jInvariant toWeierstrass(toQQ o311, {0, -2, 0})
  jInvariant toWeierstrass(F1, {2, 3, -4})
  jInvariant toWeierstrass(F2, {1, 1, -1})

  findPoint(cubics#(5,0), 3)
  toWeierstrass(toQQ cubics#(5,0), {1,0,0})
needsPackage "RationalPoints2"
rationalPoints(ideal toQQ F, Projective => true, Bound => 200, KeepAll => true)
Rp = ZZ/1009[a,b,c]
rationalPoints(ideal sub(F, Rp), Projective => true, KeepAll => true)

--------------------------------------
-- Step: take the different cubics, which have been partitioned into categories 
--   which (I hope) will keep ones that are eqivalent over ZZ or QQ in the same bin.
--   Figure out the ZZ-equivalence classes and the QQ-equivalence classes of each set.
-----
  -- here are the groups:
  {{{1, 1}, {1, 1}, {1, 1}}, 1, 3, {1, 1, 1}, {{1, 1}, {1, 1}, {1, 1}}}
  {{{1, 1}, {2, 1}}, 1, 2, {1, 1}, {{1, 1}, {2, 1}}}
  {{{1, 1}, {2, 1}}, 1, 2, {2}, {{1, 1}, {2, 1}}}
  {{{1, 1}, {2, 1}}, 1, 3, {1}, {{1, 3}}}
  {{{3, 1}}, -1, 0, {}, {{1, 1}, {1, 1}, {1, 1}}}
  {{{3, 1}}, -1, 0, {}, {{1, 1}, {2, 1}}}
  {{{3, 1}}, -1, 0, {}, {{3, 1}}}
  {{{3, 1}}, 1, 1, {1}, {{3, 1}}}
  {{{3, 1}}, 1, 2, {1}, {{1, 1}, {1, 2}}}

  -- let's do the smooth ones whose Hessian splits into linear factors
  SET = sort GROUPS #  {{{3, 1}}, -1, 0, {}, {{1, 1}, {1, 1}, {1, 1}}}
  #SET == 50
  SETcubics = SET/(lab -> lab => cubics#lab)
  #(SETcubics/last//unique) -- they are all different cubics (possibly (linearly) equiv over ZZ or QQ though)
  SETcubics/(f -> factor det hessian f#1)//netList

  -- Assumes A, phi
  findMapIdeal = (F1, F2, A, phi) -> (
        H1 := factors det hessian F1;
        H2 := factors det hessian F2;
        pt1 := H1#0#1;
        ps := select(H1, x -> support x#1 =!= {} and first degree x#1 === 1);
        qs := select(H2, x -> support x#1 =!= {} and first degree x#1 === 1);
        ps = ps/last;
        qs = qs/last;
        if #ps =!= 3 or #qs =!= 3 then << "number of factors of Hessian is not correct?";
        -- each element of qs will give us two maps ps#0 --> qs#i or - qs#i.
        (ps, qs)
        )

  -- SETcubics#0, SETcubics#1: not equivalent over QQ, but are so over QQ(cube root of 2).
  --              SETcubics#2: same
  (ps, qs) = findMapIdeal(F1 = SETcubics#0#1, F2 = SETcubics#4#1, A, phi)

  evalphi = (F,G) -> trim sub(ideal last coefficients(phi sub(F, TR) - sub(G, TR)), coefficientRing TR);
  evalIdeal = (I1, I2) -> trim sub(ideal last coefficients((phi sub(I1, TR)) % sub(I2, TR)), coefficientRing TR)

  (i,j) = (0,4) -- no match
  (i,j) = (0,5) -- over CC
  (i,j) = (0,6) -- over ZZ
  (i,j) = (0,7) -- over ZZ
  (i,j) = (0,8) -- no match
  (i,j) = (0,9) -- no match
  (i,j) = (0,10) -- over CC
  (i,j) = (0,11) -- no match
  (i,j) = (0,12) -- over CC
  (i,j) = (0,13) -- no match
  (i,j) = (0,14) -- over CC
  (i,j) = (0,15) -- over CC
  (i,j) = (0,16) -- no match
  (i,j) = (0,17) -- no match
  (i,j) = (0,18) -- no match
  (i,j) = (0,19) -- over CC
  (i,j) = (0,20) -- no match
  (i,j) = (0,21) -- over CC
  (i,j) = (0,22) -- no match
  (ps, qs) = findMapIdeal(F1 = SETcubics#i#1, F2 = SETcubics#j#1, A, phi)
  I0 = evalphi(F1,F2)
  possibles = flatten for q in {qs#0, qs#1, qs#2, -qs#0, -qs#1, -qs#2} list
    decompose ideal gens gb (I0 + evalphi(ps#0, q));
  netList possibles


  I0 = trim sub(ideal last coefficients(phi sub(F1, TR) - sub(F2, TR)), coefficientRing TR);
  I1 = trim sub(ideal last coefficients(phi sub(ps#0, TR) - sub(-qs#0,  TR)), coefficientRing TR);
  I = trim sub(I0 + I1, T)
  gbI = ideal gens gb I
  compsgbI = decompose gbI
  netList compsgbI
  A0 = A % gbI
  phi0 = map(RQ, RQ, transpose lift(A0, QQ))
  phi0 = map(TR, TR, transpose A0)
  A0 = o588
  phi0 sub(F1, TR) - sub(F2, TR)
  phi0 toQQ F1 - toQQ F2
  phi0 toQQ L1 - toQQ L2

      )

-- The following have the same invariants on GROUPSZZ.  Which are equivalent over ZZ?
  labs = {(108, 0), (126, 0), (178, 0), (179, 0), (202, 0), (206, 0), (217, 0), (220, 0), (228, 0)}
  cubes = labs/(lab -> cubics#lab)

  factor det hessian cubes#0
  factor det hessian cubes#1
  decompose ideal toQQ jacobian cubes#0   
  decompose ideal toQQ jacobian cubes#1

  F1 = cubes#0
  F2 = cubes#1
  evalphi(F1,F2)

  Tp = (ZZ/32003) (monoid T)
  I = sub(evalphi(F1,F2), Tp)
  gbTrace=1
  gens gb I;
  compsI = decompose I;
  netList (compsI/trim)
  
  
  findPoints = (F, bound) -> (
      for p in {-bound,-bound,-bound}..{bound,bound,bound} list (
          if p === {0,0,0} then continue;
          if sub(F, matrix{p}) == 0 then p else continue
          )
      )

  cubes/(f -> findPoints(f, 3))

  
  I0 = evalphi(cubes#0, cubes#1)
  I1 = evalIdeal(gens ideal(a,b), ideal(a, b+c))
  decompose ideal gens gb (I0+I1)

----------------------------------------------------------
-- Obtaining Richard's examples, of flopped h11=3 CY3's --
----------------------------------------------------------

load "~/src/stringtorics/m2-examples/richard-example/richard-db.m2"
getTopology(103, 159, 1, RZ)

-- example: h12=111
difftops111 = select(DIFFTOPS, (lab -> hh^(1,2) Xs#lab == 111))
difftops111 == {(182, 0), (183, 0), (184, 0), (185, 0), (186, 0), (187, 0), (188, 0), (189, 0), (194, 2), (195, 0), (198, 0), (199, 0), (200, 0)}

-- get Richard's
RICH = flatten flatten for h21 in findH21s DIR list 
  for pol in findPolys(DIR, h21) list
    for cy in findCYs(DIR, h21, pol) list 
      (h21, pol, cy) => toList join(getTopology(h21, pol, cy, RZ), {3, h21});

RICH111 = select(RICH, x -> x#0#0 == 111)


invariantsAll toSequence RICH111#0#1
topsToric = difftops111/(lab -> lab => invariantsAll Xs#lab)
#topsToric == 13
RICH111 = select(RICH111, x -> class x#1#0 =!= String)
topsRich = RICH111/(x -> (x#0 => invariantsAll toSequence x#1))
#unique topsToric == 9

tops111 = join(topsToric, topsRich)
hashTop111 = hashTable tops111
H = partition(x -> x#1, tops111);
H1 = for k in keys H list (H#k/first)
for lab in H1#1 list hashTop111#lab

-- NOT DONE ...
-- TODO:
--  1. grab all of the examples of Richard's
--  2. separate by invariantsAll
--  3. separate topologies, as done in 3A above.
--    There are perhaps 3 times as many than in 3A, so it should take longer...
-- get Richard's
load "~/src/stringtorics/m2-examples/richard-example/richard-db.m2"
RICH = flatten flatten for h21 in findH21s DIR list 
  for pol in findPolys(DIR, h21) list
    for cy in findCYs(DIR, h21, pol) list 
      (h21, pol, cy) => toList join(getTopology(h21, pol, cy, RZ), {3, h21});



------------------------------------------------------
-- Code to take 2 (L,F)'s and see if there are any isomorphisms over ZZ bwtween them.
findMaps

load "../m2-examples/find-maps.m2"
  {c2Form Xs#(185,0), cubicForm Xs#(185,0)}
  hRICH111 = hashTable RICH111
  drop(hRICH111#(111, 186, 1), -2)
  (L1,F1) = (c2Form Xs#(185,0), cubicForm Xs#(185,0))
  (L2,F2) = toSequence drop(hRICH111#(111, 186, 3), -2)
  ans = findMaps((L1,F1), (L2,F2))
  (A0, phi0) = ans
  phi0 toQQ L1 == toQQ L2  
  phi0 toQQ F1 == toQQ F2  

hRICH111

-- consider h12=111
labs111 = select(keys Xs, lab -> hh^(1,2) Xs#lab == 111)
h111 = hashTable(labs111/(lab -> (X := Xs#lab; lab => {c2Form X, cubicForm X})))
r111 = hashTable((keys hRICH111)/(x -> x => {hRICH111#x#0, hRICH111#x#1}))
#keys h111
#keys r111

findMaps(h111#(first keys h111), r111#((keys r111)#1))
findMaps(h111#((keys h111)#1), r111#((keys r111)#2))
findMaps(h111#((keys h111)#1), r111#((keys r111)#3))

hr111 = hashTable join(pairs h111, pairs r111)
H = partition(lab -> invariantsAll(hr111#lab#0, hr111#lab#1, 3, 111), keys hr111);
#(keys H) == 11
netList for k in keys H list {k, #H#k}

invset1 = first select(keys H, k -> #H#k == 15)
SET1 = rsort H#invset1

for x from 0 to #SET1-1 list for y from x+1 to #SET1-1 list
(
  ans = (x,y) => findMaps(hr111#(SET1#x), hr111#(SET1#y));
  print ans;
  ans)
  

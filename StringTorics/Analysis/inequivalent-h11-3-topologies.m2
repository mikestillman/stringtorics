-- We compute all of the different topologies for h11=3 toric CY3 hypersurfaces
-- Caveats:
--  1 example is non favorable: #232
--  N examples have Cl V torsion.  Currently our methods do not handle these 
--    examples.

restart
debug needsPackage "StringTorics" -- the debug is because some functions are not yet exported.
DB3 = "../Databases/cys-ntfe-h11-3.dbm"
DB3 = "/Users/mike/src/stringtorics/StringTorics/Databases/cys-ntfe-h11-3.dbm"

R = ZZ[a,b,c]
RZ = R
RQ = QQ (monoid R);
elapsedTime (Qs, Xs) = readCYDatabase(DB3, Ring => R);

-- We collect the nontorsion, torsion, favorable, nonfavorable's.
   torsions = for k in keys Qs list (
      istor := prune coker matrix rays Qs#k;
      if not isFreeModule istor then k else continue
      )
   nonfavorables = for k in keys Qs list (
       if not isFavorable Qs#k then k else continue
      )

torsionCYs = sort select(keys Xs, lab -> member(first lab, torsions))
nonfavorableCYs = sort select(keys Xs, lab -> member(first lab, nonfavorables))
favorableXs = sort select(keys Xs, lab -> not member(first lab, nonfavorables))
assert(torsionCYs == {(0, 0), (9, 0), (10, 0), (55, 0), (62, 0)})
assert(nonfavorableCYs == {(232,0)})

-- how many of these have nonfavorable dual polytopes
-- these are the ones that are not nec general in moduli
-- select(sort keys Qs, lab -> (ans := not isFavorable polar Qs#lab; print ans; ans))

---------------------------------------------------------------
-- Next step: How many of these 306 are distinct topologies? --
---------------------------------------------------------------
  allXs = sort keys Xs

  allXs = favorableXs -- Numbers refer to this situation.
  allT = topologySet(allXs, Xs);
  info allT -- 274 possibly different topologies

  allT = separateIfDifferent(allT, invariantsH11H12)
  info allT

  allT = separateIfDifferent(allT, hubschInvariants)
  info oo

  PC = pointCounter(RZ, Projective => true);
  allT = separateIfDifferent(allT, pointCounts_PC)
  info allT
  --partProj = partition(lab -> pointCounts_PC Xs#lab, sort keys Xs);

  allT = separateIfDifferent(allT, hessianInvariants)
  info oo
  --partition(lab -> hessianInvariants Xs#lab, sort keys Xs);

  -- This gives the same benefit as hessianInvariants, if we first do: pointcounts, hubsch, h11h12.
  allT = separateIfDifferent(allT, X -> polynomialContent det hessian cubicForm X)
  info oo

  allT = separateIfDifferent(allT, cubicConductorInvariants)
  info oo -- on its own, this breaks into 67 different cases. (but 145 are in one class!)
    -- still at 174 if we do the ones above this too.

  allT = separateIfDifferent(allT, cubicLinearConductorInvariants) -- this is a good one!
  info oo -- on its own, breaks into 168 classes, BUT no new ones after using above

  allT = separateIfDifferent(allT, singularContents) -- this is a good one!
  info oo -- on its own, breaks into 168 classes, BUT no new ones after using above

  newsets = for reps in representatives allT list (
      partition(lab -> classifyExtremalCurves Xs#lab, reps)
      )

  netList extremalCurveInvariant(Xs#(2,0))
  netList extremalCurveInvariant(Xs#(2,1))
  
  newsets = for reps in representatives allT list (
      partition(lab -> netList extremalCurveInvariant Xs#lab, reps)
      )
  
  -- -- Don't use affine points.  They are likely equivalent to projective points, and *much* slower.
  -- PC = pointCounter(RZ, Projective => false); -- very pricy...
  -- allTA = separateIfDifferent(allT, pointCounts_PC)
  -- info allTA
  -- partAffine = partition(lab -> pointCounts_PC Xs#lab, sort keys Xs);
  -- (values partProj)/sort//sort
  -- (values partAffine)/sort//sort
  -- oo === ooo
  
  allT1 = combineIfSame(allT, X -> (c2Form X, cubicForm X))
  equivalences allT1
  info allT1 

  elapsedTime allT2 = separateIfDifferent(allT, invariantsAll) -- 18 seconds
  info allT2 
  netList representatives allT2
  equivalences allT2

  elapsedTime allT2 = separateIfDifferent(allT1, invariantsAll) -- 18 seconds
  info allT2 
  netList representatives allT2
  equivalences allT2

  -- We have two ways to proceed here.

  -- VERSION #1: use GV invariants to find equivalences

  elapsedTime allT3a = combineByGV(allT2, DegreeLimit => 5); -- 14 sec
    info allT3a
    netList representatives allT3a
  elapsedTime allT3b = combineByGV(allT3a, DegreeLimit => 10); -- 4.4 sec
    info allT3b
    netList representatives allT3b
    equivalences allT3b
  allT3 = allT3b

  for x in representatives allT3 list (
      if #x != 2 then continue;
      ans := x => findIsomorphism((toSequence x)/(lab -> Xs#lab));
      print ans;
      ans
      )

  for x in representatives allT3 list (
      if #x == 2 then continue;
      tries := subsets(x, 2);
      for t in tries list (
        ans := t => findIsomorphism((toSequence t)/(lab -> Xs#lab));
        print ans;
        )
      )


  -- find all 186 representatives.  If possible (it is), choose favorable, and
  -- also if possible choose Q s.t. polar Q is favorable (i.e. no missing complex structure.
  
  for s in allT3#"Sets" list for L in s list (
      newa := for a in L list if instance(a, List) then first a else a;
      qs := newa/first//unique//sort;
      j := select(qs, lab -> isFavorable polar Qs#lab);
      if #j == 0 then newa#0 else (
        j1 := position(newa, lab -> first lab == j#0);
        newa#j1
        )
      )
  

    Ts = for lab in allXs list (
        X := Xs#lab;
        lab => {c2Form X, cubicForm X, hh^(1,1) X, hh^(1,2) X}
        );
    assert(#Ts == 306)
    hashTs = hashTable Ts;
    keyTs = Ts/first//sort -- 306 of these.

    reps = representatives allT3    
    partitionH113sByTopology(reps#11, hashTs, RQ)

    (A, phi) = genericLinearMap RQ
    findMaps(hashTs#(232,0), hashTs#(233,0), A, phi, RQ)
    equivs = equivalences(allT3, IgnoreSingles=>false)
    result = flatten for L in representatives(allT3, IgnoreSingles=>false) list (
        newset0 := partitionH113sByTopology(L, hashTs, RQ);
        for x in keys newset0 list (
            -- we include x, and its value, and for each element in there, attach its rest also.
            val := newset0#x;
            {{x}|val|equivs#x|for v in val list equivs#(first v)}
            )
        )
    new TopologySet from {
        "Sets" => result,
        "CYHash" => allT3#"CYHash"
        }

  -- here we check that the equivalent ones are really equivalent.
  equivs = sort flatten for x in allT3a#"Sets" list (
      y1 := for y in x list if #y > 1 then y else continue;
      if #y1 == 0 then continue else y1)
  equivs_0

  
  isEquivalent(Xs#(1,0), Xs#(2,0), equivs_0_1_1)
  -- These seem to all be equivalent.
  for e in equivs do (
      X1 = Xs#(e#0);
      for i from 1 to #e-1 do (
          if instance(e#i, Sequence) then (
              -- check equality of c2, cubic
              )
          else if not isEquivalent(X1, Xs#(e#i#0), e#i#1)
            then << e#0 << " and " << e#i#0 << " should be equivalent, but seem not to be" << endl;
          )
      )

  elapsedTime allT3b = combineByGV(allT3a, DegreeLimit => 10);

  equivs = sort flatten for x in allT3b#"Sets" list (
      y1 := for y in x list if #y > 1 then y else continue;
      if #y1 == 0 then continue else y1)
  equivs_0

  
  isEquivalent(Xs#(1,0), Xs#(2,0), equivs_0_1_1)
  -- These seem to all be equivalent.
  for e in equivs do (
      X1 = Xs#(e#0);
      for i from 1 to #e-1 do (
          if instance(e#i, Sequence) then (
              -- check equality of c2, cubic
              )
          else if not isEquivalent(X1, Xs#(e#i#0), e#i#1)
            then << e#0 << " and " << e#i#0 << " should be equivalent, but seem not to be" << endl;
          )
      )

  elapsedTime allT3c = combineByGV(allT3b, DegreeLimit => 15);

  equivs = sort flatten for x in allT3c#"Sets" list (
      y1 := for y in x list if #y > 1 then y else continue;
      if #y1 == 0 then continue else y1)
  equivs_0

  elapsedTime allT3 = combineByGV(allT3b, DegreeLimit => 20);


  elapsedTime allT3 = combineByGV allT2; -- HERE XXX
    elapsedTime allT3 = separateByGV allT2 -- 45 sec
    info allT3
    #allT3a#"Sets"
    #allT3b#"Sets"
    allT3a#"Sets"/(x -> x/length//sum)//sum
    allT3b#"Sets"/(x -> x/length//sum)//sum
    -- Tally{1 => 159}
    --       2 => 9
    --       3 => 2
    -- UPSHOT: There are at least 170 topologies, but at most 170+9+2*2 = 183 different topologies.

    -- these are the possible pairs that could be equivalent.
    -- if the GV code is giving the correct answer (not clear, as there is a degree bound given)  
    -- then these should all not have any solution over the integeres.
    onestocheck = flatten for x in allT3#"Sets" list if #x == 1 then continue else (
        subsets(x/first, 2)
        )
    elapsedTime for x in onestocheck list (
        (J,A) := getEquivalenceIdeal(x#0, x#1, Xs);
        x => for j in decompose J list A % j
        )
    -- these are all the unit ideal, so none of these have any equivalence
    -- over QQ or RR, let alone ZZ.

    sort (flatten allT3#"Sets")/(x -> sort prepend(first x, (drop(x, 1))/first))
  -- Upshot: There are 187 different topologies for CY3 which are have h11=3 and are toric hypersurfaces (a ll Batyrev construction)

  -- VERSION #2: Use ansatz to find equivalences, DOES NOT use GV invariants
  -- This is somewhat older code, but hopefully it gives the same answer!
  -- Actually: this version can handle the torsion examples, just not the non-favorable one.
    T1s = for lab in allXs list (
        X := Xs#lab;
        lab => {c2Form X, cubicForm X, hh^(1,1) X, hh^(1,2) X}
        );
    assert(#T1s == 306)
  
    -- now, how many of these are the same?
    285 == # unique values hashTable T1s -- ?? 286 + 5 -- the 5 are the torsion but favorables.

    Ts = T1s; -- these are all the label => topology pairs we have (291 here).
    hashTs = hashTable Ts;
    keyTs = Ts/first//sort -- 306 of these.

    elapsedTime H = partition(lab -> invariantsAll toSequence (hashTs#lab), keyTs); -- 
    #keys H == 173
    (keys H)/(k -> #H#k)//tally
    
    elapsedTime INV = for k in keys H list k => elapsedTime partitionH113sByTopology(H#k, hashTs, RQ); -- 196 seconds
    tally for x in INV list #(keys x#1) -- 162 have one group, 9 have 2, 2 have 3.
    162 + 18 + 6
    INV/last
    netList oo -- torsions? (0,0),
    torsionCYs -- (9,0), (10,0) are equivalent; (55,0), (62,0) are equivalent.
    INV/last/(x -> for k in keys x list {k}|x#k)//flatten
    -- UPSHOT:
    -- found 186 different topologies
    -- 3 of these correspond to torsion examples (the 5 appear in 3 sets).
    -- 183 correspond to torsion-free examples.  This matches the code above.
    -- UNKNOWN YET: the nonfavorable example: is that by itself, or equivalent to another?
    -- here are the h12=165 examples.  The last 5 are all the same topology.
    --   {(232, 0), (233, 0), (234, 0), (235, 0), (236, 0), (237, 0)}
    --   Richhard's database indicates that there are no flops of these examples.
    -- the first is possibly a new topology.  But is it?
    
-- QUESTIONS/TODO:
--  1. which invariants do we really need here?
--  2. can we use GV code for flopped examples?
--  3. how can we use the first code above, using GV invariants, with the flopped examples.

  -- VERSION #3.  Do the same, but on these examples, conbined with Richard's database (includes flopped CY3s)
  -- How to handle these?
  -- How to compute the GV cone directly?
    load "../../m2-examples/richard-example/richard-db.m2"
    DIR = "~/Dropbox/Collaboration/Physics-Liam/Inequivalent CYs/h11_3NontoricNewFlopCode/"
    T2s = flatten flatten for h21 in findH21s DIR list 
      for pol in findPolys(DIR, h21) list
        for cy in findCYs(DIR, h21, pol) list 
          (h21, pol, cy) => toList join(getTopology(DIR, h21, pol, cy, R), {3, h21});
    assert(#T2s == 487)
    T2s/first/toList/isToric_DIR//tally -- 220 are not toric, 267 are toric. QUESTION: why 267??
    Ts = join(T1s, T2s); -- ours and Richard's.  Ours are here to calibrate the toric examples.
    hashTs = hashTable Ts;
    keyTs = Ts/first//sort
    #keyTs == 792

    elapsedTime H = partition(lab -> invariantsAll toSequence (hashTs#lab), keyTs); -- 
    #keys H == 277
    (keys H)/(k -> #H#k)//tally
    
    elapsedTime INV = for k in keys H list k => elapsedTime partitionH113sByTopology(H#k, hashTs, RQ); -- 
    tally for x in INV list #(keys x#1) -- 257 have one group, 18 have 2, 2 have 3.
    topList = INV/last/(x -> for k in keys x list {k}|(x#k)/first)//flatten
    select(topList, x -> all(x, x1 -> #x1 == 3))


-- Upshot: I am pretty sure it is 186 toric topologies.
--  the one non-favorable (232,0) is equivalent to (233,0).
--  the torsion ones form 3 separate topologies.
--  non-toric phases: 113
-- 186 + 113 = 299.

(J, A) = getEquivalenceIdeal((232,0),(233,0),Xs)
-- This shows they have the same intersection and c2 forms.
(L1,F1) = (c2Form Xs#(232,0), cubicForm Xs#(232,0))
(L2,F2) = (c2Form Xs#(233,0), cubicForm Xs#(233,0))
phi = map(R, R, transpose matrix{{1,0,3},{0,1,1},{0,0,1}})
phi L1 == L2
phi F1 == F2

-- now let's try the equivalenceIdeal.  It should find it too!
(A,phi) := genericLinearMap RQ;
T = target phi
(L1,F1) = (sub(L1,T), sub(F1,T))
(L2,F2) = (sub(L2,T), sub(F2,T))
phi L1 - L2

-------------------------------------------
-- Automorphisms --------------------------
-------------------------------------------
automorphisms Qs#0
automorphisms Qs#1
netList annotatedFaces Qs#0
topenum = 0
select(sort keys Xs, lab -> first lab == topenum)

Xs = findAllCYs Qs#2
tri1 = restrictTriangulation(2, Xs#0)
tri2 = restrictTriangulation(2, Xs#1)
tri3 = restrictTriangulation(2, Xs#2)
tri1 === tri2
tri1 === tri3
tri2 === tri3
netList restrictTriangulation(Xs#0)
findAllFRSTs(Qs#2)

--------------------------------------
-- (166,0), (171,0)
A = matrix{{0,-3,1},{0,1,0},{1,-1,2}}
X1= Xs#(166, 0)
X2 = Xs#(171, 0)
phi = map(RZ, RZ, A)
phi vars RZ
phi cubicForm X1 == cubicForm X2
phi c2Form X1 == c2Form X2
classifyExtremalCurves X1
classifyExtremalCurves X2
A * transpose matrix{{1,0,0}} 

-- curves are transformed via A^T
matrix{{1,0,0}} * (transpose A)
matrix{{1,1,0}} * (transpose A)
matrix{{-2,0,1}} * (transpose A)

-- divisors are transformed via A^-1
matrix{{0,1,0}} * A^-1
matrix{{0,0,1}} * A^-1
matrix{{1,-1,2}} * A^-1

-- The ring map is: a,b,c -> matrix{{a,b,c}} * A.
-- or, a\\b\\c --> (transpose A) * a\\b\\c
A * transpose matrix{{a,b,c}}
-- matrix{{a,b,c}} --> matrix{{a,b,c}} * A
phi.matrix

-- if L = a*x + b*y +c*z, what is phi(L) in terms of the
-- row or column vector {x,y,z}.

-----------------------------------------------
-- Given allT2, i.e. first separate by invariantsAll, now use IntegerEquivalences to
-- see if they are the same.
reps = {{(1, 0), (2, 0), (2, 1)},
    {(4, 0), (6, 0)},
    {(9, 0), (10, 0)},
    {(11, 0), (13, 0)},
    {(12, 0), (15, 0)}, {(20, 0), (21, 0)},
    {(23, 0), (24, 0), (24, 1), (24, 2), (31, 0)},
    {(26, 0), (34, 0), (37, 0)}, {(27, 0), (32, 0)},
    {(28, 0), (40, 0)}, {(29, 0), (35, 0)}, {(30, 0), (36, 0)},
    {(47, 0), (48, 0)}, {(53, 0), (53, 1)}, {(54, 0), (59, 0)},
    {(55, 0), (62, 0)}, {(56, 0), (74, 0)}, {(57, 0), (67, 0), (73, 0)},
    {(60, 0), (64, 0)}, {(61, 0), (77, 0)}, {(63, 0), (72, 1)}, {(63, 1), (72, 0)},
    {(81, 0), (85, 0)}, {(87, 0), (88, 1)}, {(90, 0), (92, 0)}, {(94, 0), (104, 1)},
    {(94, 1), (104, 0)}, {(95, 0), (103, 0)}, {(96, 0), (98, 0)}, {(97, 0), (102, 0)},
    {(115, 0), (120, 0)}, {(122, 0), (123, 0), (125, 0)}, {(130, 0), (132, 0), (133, 0)},
    {(134, 0), (136, 0), (137, 0)}, {(138, 0), (141, 0), (142, 0), (143, 0),
        (143, 1), (143, 2), (149, 0), (152, 0)},
    {(139, 0), (144, 0), (145, 0), (148, 0), (150, 0), (154, 0)},
    {(140, 0), (151, 0), (153, 0)}, {(146, 0), (157, 0)}, {(159, 0), (162, 0), (163, 0)},
    {(166, 0), (171, 0), (177, 0)}, {(167, 0), (172, 0)}, {(168, 0), (174, 0), (176, 0)},
    {(170, 0), (173, 0), (175, 0)}, {(178, 0), (179, 0), (180, 0)}, {(182, 0), (188, 0), (196, 0)},
    {(183, 0), (194, 1), (195, 0), (197, 0)}, {(184, 0), (192, 0)}, {(185, 0), (200, 0)},
    {(186, 0), (190, 0), (193, 0)}, {(187, 0), (200, 1)}, {(189, 0), (191, 0), (199, 0)},
    {(201, 0), (205, 0), (208, 1)}, {(201, 1), (208, 0)}, {(202, 0), (206, 0)},
    {(203, 0), (204, 0)}, {(207, 0), (210, 0)}, {(211, 0), (214, 0)},
    {(215, 0), (221, 0)}, {(216, 0), (218, 0), (222, 0)}, {(217, 0), (220, 0)},
    {(219, 0), (223, 0)}, {(227, 0), (229, 0)}, {(230, 0), (231, 0)},
    {(233, 0), (234, 0), (235, 0), (236, 0), (237, 0)},
    {(239, 0), (240, 0)}, {(241, 0), (242, 0)}}
needsPackage "IntegerEquivalences"
findEquivalence (Sequence, Sequence) := (lab1, lab2) -> (
    X1 := Xs#lab1;
    X2 := Xs#lab2;
    LF1 := {c2Form X1, cubicForm X1};
    LF2 := {c2Form X2, cubicForm X2};
    findEquivalence(LF1, LF2)
    )
findEquivalence((4,0), (6,0))
findEquivalence((9,0), (10,0))
findEquivalence((11,0), (13,0))
rep2s = flatten(reps/(x -> subsets(x,2)))
for x in rep2s list (
    result := findEquivalence(x#0, x#1);
    print result;
    x => result
    )

findEquivalences List := (LFs) -> (
    -- return a pair: hash table: keys are (some of) the indices of LFs
    -- value is a list of {index, matrix}.
    -- list of inconclusives.
    )

for x in reps list (
    x2 := subsets(x, 2);
    for y in x2 list findEquivalence(y#0, y#1)
    )


netList oo

for x in reps list (
    for y in x list classifyExtremalCurves Xs#y
    )

reps#1
X1 = Xs#(4,0)
use ring c2Form X1
classifyExtremalCurves X1
LF1 = {c2Form X1, cubicForm X1}
LF2 = {c2Form X1 + 2*4*a, cubicForm X1 - 4*a^3}
LF3 = {c2Form X1 + 2*(-a+b-c), cubicForm X1 - (-a+b-c)^3}
LF4 = {c2Form X1 + 2*(-b+2*c), cubicForm X1 - (-b+2*c)^3}
factorsByType det hessian LF1_1
factorsByType det hessian LF2_1
factorsByType det hessian LF3_1 -- different...
factorsByType det hessian LF4_1 -- different
findEquivalence(LF1, LF2) -- same, square of the matrix is 1.
A1 = oo_1
findEquivalence(LF1, LF3) -- diff
findEquivalence(LF1, LF4) -- diff
findEquivalence(LF3, LF4) -- diff.
classifyExtremalCurves X1

transpose matrix degrees X1
toricMoriConeCap X1
dualCone posHull transpose matrix oo
rays oo
-- need to determine nef cone, Mori cone of flopped CY3 LF2.  This is equivalent to X1.
-- but we need to get the maps correct...
effV = rays posHull transpose matrix degrees X1
moriX = transpose matrix toricMoriConeCap X1 -- tentative?
nefX = rays dualCone posHull moriX
for a in (-3,-3,-3)..(3,3,3) list (b := hh^0 OO_X1 a; if b > 0 then a => b else continue)
oo/first/toList//matrix//transpose
posHull oo
rays oo

(transpose A1) * nefX

hh^0 OO_X1(2,2,1)

effV
A1 * effV

hh^0 OO_X1(1,0,0)
hh^0 OO_X1(-1,1,1)

hh^0 OO_X1(0,1,0)
hh^0 OO_X1(-1,1,1)

hh^0 OO_X1(-2,1,1)
hh^0 OO_X1(2,-1,-1)

hh^0 OO_X1(-2,1,1)
hh^0 OO_X1(-1,3,2) -- not the same


(transpose A1) * nefX


-- {(11, 0), (13, 0)}
reps#3
X1 = Xs#(reps#3#0)
for x in subsets(reps#3, 2) list findEquivalence(x#0, x#1)
use ring c2Form X1
classifyExtremalCurves X1
LF1 = {c2Form X1, cubicForm X1}
classifyExtremalCurves X1
LF2 = {c2Form X1 + 2*6*b, cubicForm X1 - 6*b^3}
LF3 = {c2Form X1 + 2*6*a, cubicForm X1 - 6*a^3}
factorsByType det hessian LF1_1
factorsByType det hessian LF2_1
factorsByType det hessian LF3_1

findEquivalence(LF1, LF2) -- same, square of the matrix is 1.
A1 = oo_1
findEquivalence(LF1, LF3) -- same
A2 = oo_1
findEquivalence toSequence reps#3
A0 = oo_1

findEquivalence(LF2, LF3) -- flip a, b.  this is A2*A1
A1^4 == 1
A2^2 == 1
(A1*A2)^2 == 1
(A1^2*A2)
A0^4

-- it looks like all flops go to "Weyl" flops (their 

X = Xs#(12,0)
X = Xs#(15,0)
X = Xs#(20,0)
X = Xs#(21,0)
X = Xs#(23,0)
X = Xs#(24,0)
X = Xs#(24,1)
X = Xs#(24,2)
X = Xs#(31,0)

X = Xs#(26,0)
X = Xs#(34,0)
X = Xs#(37,0)

X = Xs#(27,0)
X = Xs#(32,0)

X = Xs#(28,0) -- this one has an interesting symmetric flop.
X = Xs#(40,0)

findEquivalence((28,0), (40,0))

-----------------------------------------------------
-- Example: CY3 with an elliptic ruled surface
--  whose flop goes to same topology CY3, but whose effective surfaces change.
-----------------------------------------------------
-- Trying to find my bug or logic error in this example.
-- Problem: I get a matrix A corresponding to what I think is a flop.
--   This matrix implements the isomorphism on divisor classes/curves
--   of X and its flopped X', which is isomorphic to X.
-- BUT: this means that A^-1 times any effective class should be effective.
--   actually, I think it means that all of q, A^-1*q, A^-2*q, ....
--   should all have the same o-th cohomology...
--   This seems to not be the case.
restart
debug needsPackage "StringTorics" -- the debug is because some functions are not yet exported.
DB3 = "../Databases/cys-ntfe-h11-3.dbm"
DB3 = "./StringTorics/Databases/cys-ntfe-h11-3.dbm"

R = ZZ[a,b,c]
RZ = R
RQ = QQ (monoid R);
elapsedTime (Qs, Xs) = readCYDatabase(DB3, Ring => R);

X = Xs#(28,0) -- this one has an interesting symmetric flop.
classifyExtremalCurves X
flopset = flatten for p in pairs classifyExtremalCurves X list (
    if p#0#0 =!= "FLOP" then continue;
    nv := p#0#1#0; -- only one that is non-zero??!
    for x in p#1 list {nv, sum for i from 0 to #x-1 list x#i * R_i}
    )
for p in flopset list (
    LF2 := {c2Form X + 2 * p#0 * p#1, cubicForm X - p#0 * p#1 ^3};
    findEquivalence({c2Form X, cubicForm X}, LF2) 
    )
A = last last oo  -- from (28,0)
(A - 1)^3 == 0
det A
A^2
A^3

-- this is the same A as above (for X = Xs#(28,0)
A == matrix {{-1, 0, -4}, {-1, 1, 2}, {1, 0, 3}}
label X == (28,0)

V = normalToricVariety(X, CoefficientRing => ZZ/32003)
X' = completeIntersection(V, {-toricDivisor V}, Basis => {0,1,4}, Variables => {symbol a, symbol b, symbol c})
hh^0(OO_X'(3,5,-1))
hh^0(OO_X(3,5,-1))
hh^0(OO_X'(5,14,-2))
hh^0(OO_X(5,14,-2))

Q = transpose matrix degrees X -- these are the inherited effective divisors
hh^* OO_X(-2,-1,1)
A^-1 * Q
for q in entries transpose Q list hh^0 OO_X toSequence q
for q in entries transpose (A^-1*Q) list hh^0 OO_X toSequence q
for q in entries transpose (A^-2*Q) list hh^0 OO_X toSequence q
-- (-2,-1,1) is effective on V, therefore X.
-- On X', isom to X''=X, one thinks that (-2, -5, 1) should be effective on X'' too.
-- The isomorphism sends the effective divisor (-2,-1,1) to (-2,-5, 1), but this is not
-- effective on X.
-- Wnat is the surface (-2,-1,1) on X?
hh^* OO_X(-2, -1, 1) -- {1, 1, 0, 0}
  -- is this an elliptic ruled surface? This is D_6.
  -- How do I tell?
  
basis({-2,-1,1}, ring V)

-- Let's analyze this surface E = D_6 \cap X.

rays V
-- can we nail down what the Mori cone is?
toricMoriConeCap X
Q = transpose matrix degrees X
rays posHull transpose matrix degrees X -- effV...
nefXsmall = rays dualCone posHull transpose matrix toricMoriConeCap X; -- potential nef X...
 -- pos hull of the rays {{0, 1, 0}, {1, 1, 0}, {0, 1, 1}}

for q in entries transpose Q list hh^0 OO_X toSequence q
for q in entries transpose (A^-1*Q) list hh^0 OO_X toSequence q
A^-1 * Q
A^-2 * Q
for q in entries transpose (A^-2*Q) list hh^0 OO_X toSequence q
for q in entries transpose (A*Q) list hh^0 OO_X toSequence q
for q in entries transpose (A^2*Q) list hh^0 OO_X toSequence q

A^-1 * nefX
 
for a in (-3,-3,-3)..(3,3,3) list (b := hh^0 OO_X a; if b > 0 then a => b else continue)
for a in (-3,-3,-3)..(3,3,3) list (b := hh^* OO_X a; if b === {1,0,0,0} then a else continue)
for a in (-3,-3,-3)..(3,3,3) list (b := hh^0 OO_X a; if b == 1 then toList a else continue)
oo/first/toList//matrix//transpose
posHull oo
rays oo
-- nef X seems to be gen by rays {{1, 0, 0}, {0, 1, 0}, {-2, -1, 1}} -- NO!!
hh^0 OO_X(1,0,0)
hh^0 OO_X(-1,-1,1)
hh^0 OO_X(-3,2,2)

(transpose A) * transpose matrix{{1,0,0}}
hh^0 OO_X(-1,0,4)

A * transpose matrix{{1,0,0}}
A^2 * transpose matrix{{1,0,0}}

transpose matrix{{-1,1,-1}}
(transpose A)^-1 * oo

(transpose A) * transpose matrix{{1,0,0}}
hh^0(OO_X(1,0,0))
hh^0(OO_X(-1,-1,1))

hh^0(OO_X(0,0,1))
hh^0(OO_X(-4,2,3))
hh^0(OO_X(1,0,3))
hh^0(OO_X(4,6,-1))
hh^0(OO_X(3,5,-1))

hh^0(OO_X(5,14,-2))
hh^0(OO_X(8,20,-3))

-- Question: is the curve {1,0,0} an extremal ray in the Mori cone of X?
--
gvX = gvInvariants(X, DegreeLimit => 24)
(keys gvX)/toList//matrix//transpose
moriX = rays posHull oo -- columns are curve classes
nefX = rays dualCone posHull moriX
for k in keys gvX list (
    if minors(2, matrix{toList k, {0,0,1}}) == 0 then - toList k else toList k)
transpose matrix oo
rays posHull oo
for k in keys gvX list (
    matrix{toList k, {0,0,1}}
    )
keys gvX



-- which curves are not potent?
--
hf = heft X
for k in keys gvX list (
    if dotProduct(hf, toList k) <= 8 then toList k else continue
    )
nilpcurves = for k in oo list (
    if gvX#?(toSequence (3 * k)) then continue else k
    )
for k in nilpcurves list gvX#(toSequence k)
(transpose A) * transpose matrix nilpcurves
(transpose A^-1) * transpose matrix nilpcurves
pt = base(a,b,c)
Xa = abstractVariety(X, pt)
IX = intersectionRing Xa
use IX
E = -2*t_0 - t_1 + t_4
(chern_2 tangentBundle Xa) * E == 0 -- I think that any ruled surface on X must have this be true?
c2Form X -- this is >= 0 on the effective cone?  Or just on the nef cone?
-- what is the class of a fiber on E?

flatten for i from 0 to numgens IX - 2 list for j from i+1 to numgens IX-1 list (
    D := IX_i * IX_j;
    {i,j} => (integral(D * IX_0), integral(D * IX_1), integral(D * IX_4))
    )
V = normalToricVariety X
Va = abstractVariety(V, pt)
IV = intersectionRing Va
for x in subsets(splice{0..numgens IV-1}, 3) list (
    D := IV_(x#0) * IV_(x#1) * IV_(x#2);
    x => (integral(D * IV_0), integral(D * IV_1), integral(D * IV_4))
    )

flatten for i from 0 to numgens IV - 3 list for j from i+1 to numgens IX-1 list (
    D := IX_i * IX_j;
    {i,j} => (integral(D * IX_0), integral(D * IX_1), integral(D * IX_4))
    )



c1 = transpose matrix{{0,0,1}}

(transpose A) * c1
(transpose A^2) * c1
gvX#(1,0,3)
gvX#(2,0,5)
gvX#(3,0,7)
select(keys gvX, k -> gvX#k == 48)
(transpose A^3) * c1
gvX#(-1,0,-1)

(keys gvX)/toList//matrix//transpose
moriX = rays posHull oo -- columns are curve classes
nefX = rays dualCone posHull moriX
classifyExtremalCurves X
gvRay
gvX#(1,0,1)
gvX#(2,0,2)
gvX#(3,0,3)
gvX#(4,0,4)

-- So: the Mori cone of X seems likely to be generated by {1,1,0}, {0,0,1} (gv=48), {-1,1,-1} (gv=1).
-- The nef cone then is generated by rays {0,1,0}, {1,1,0}, {-1,0,1}.

toricMoriConeCap X
rays posHull transpose matrix oo
hh^0 OO_X(1,2,1) -- 19 sections...
hh^0 OO_X(-1,1,2) -- 19 sections...
hh^0 OO_X(0,1,1) -- 9.
I = intersectionRing abstractVariety(X, base(a,b,c))
integral((a*t_0 + b*t_1 + c*t_4) * (t_0 + 2*t_1 + t_4)^2)


------------------------------------------------------------------
-- Another try May 27, 2024 to find non-equivalent CY3's at h11=3 in KS.
restart
debug needsPackage "StringTorics" -- the debug is because some functions are not yet exported.
DB3 = "/Users/mike/src/stringtorics/StringTorics/Databases/cys-ntfe-h11-3.dbm"

RZ = ZZ[a,b,c]
elapsedTime (Qs, Xs) = readCYDatabase(DB3, Ring => RZ);
allXs = sort keys Xs
diffsets = elapsedTime partition(lab -> invariants Xs#lab, allXs);
tally apply(keys diffsets, k -> #diffsets#k)
-- Tally{1 => 100}
--       2 => 43
--       3 => 20
--       4 => 1
--       5 => 1
--       6 => 2
--       8 => 1
set2 = for k in sort keys diffsets list if #diffsets#k == 2 then diffsets#k else continue
for x in set2 list (
    X1 := Xs#(x#0);
    X2 := Xs#(x#1);
    ans := findEquivalence({c2Form X1, cubicForm X1}, {c2Form X2, cubicForm X2});
    << x << "  " << ans << endl;
    x => ans
    )
netList oo

set3 = for k in sort keys diffsets list if #diffsets#k == 3 then diffsets#k else continue
for y in set3 list (
    << "--- " << y << " -----" << endl;
    for x in subsets(y, 2) list (
      X1 := Xs#(x#0);
      X2 := Xs#(x#1);
      ans := findEquivalence({c2Form X1, cubicForm X1}, {c2Form X2, cubicForm X2});
      << x << "  " << ans << endl;
      x => ans
      )
  )
netList oo

set4 = for k in sort keys diffsets list if #diffsets#k == 4 then diffsets#k else continue
for y in set4 list (
    << "--- " << y << " -----" << endl;
    for x in subsets(y, 2) list (
      X1 := Xs#(x#0);
      X2 := Xs#(x#1);
      ans := findEquivalence({c2Form X1, cubicForm X1}, {c2Form X2, cubicForm X2});
      << x << "  " << ans << endl;
      x => ans
      )
  )
netList oo

set5 = for k in sort keys diffsets list if #diffsets#k == 5 then diffsets#k else continue
for y in set5 list (
    << "--- " << y << " -----" << endl;
    for x in subsets(y, 2) list (
      X1 := Xs#(x#0);
      X2 := Xs#(x#1);
      ans := findEquivalence({c2Form X1, cubicForm X1}, {c2Form X2, cubicForm X2});
      << x << "  " << ans << endl;
      x => ans
      )
  )
netList oo
netList transpose ooo -- all the same

set6 = for k in sort keys diffsets list if #diffsets#k == 6 then diffsets#k else continue
for y in set6 list (
    << "--- " << y << " -----" << endl;
    for x in subsets(y, 2) list (
      X1 := Xs#(x#0);
      X2 := Xs#(x#1);
      ans := findEquivalence({c2Form X1, cubicForm X1}, {c2Form X2, cubicForm X2});
      << x << "  " << ans << endl;
      x => ans
      )
  )
netList oo
netList transpose ooo -- each set consists of all equivalenct tops.  So get 2 new tops here.

set8 = for k in sort keys diffsets list if #diffsets#k == 8 then diffsets#k else continue
for y in set8 list (
    << "--- " << y << " -----" << endl;
    for x in subsets(y, 2) list (
      X1 := Xs#(x#0);
      X2 := Xs#(x#1);
      ans := findEquivalence({c2Form X1, cubicForm X1}, {c2Form X2, cubicForm X2});
      << x << "  " << ans << endl;
      x => ans
      )
  )
netList transpose ooo -- need to see how many different ones there are.


  allXs = sort keys Xs
  allT = topologySet(allXs, Xs);
  info allT -- 275 possibly different topologies

  allT = separateIfDifferent(allT, invariantsH11H12)
  info allT

  allT = separateIfDifferent(allT, hubschInvariants)
  info oo

  PC = pointCounter(RZ, Projective => true);
  allT = separateIfDifferent(allT, pointCounts_PC)
  info allT
  --partProj = partition(lab -> pointCounts_PC Xs#lab, sort keys Xs);

  allT = separateIfDifferent(allT, hessianInvariants)
  info oo
  --partition(lab -> hessianInvariants Xs#lab, sort keys Xs);

  elapsedTime allT' = separateIfDifferent(allT, x -> elapsedTime invariants x)

  representatives allT
  for x in representatives allT list (
      X1 := Xs#(x#0);
      for y in drop(x, 1) list (
          X2 := Xs#y;
          ans := findEquivalence({c2Form X1, cubicForm X1}, {c2Form X2, cubicForm X2});
          {x#0, y} => ans
      ))

-- We compute all of the different topologies for h11=4 toric CY3 hypersurfaces
-- Caveats:
--  non-favorables:  These are handled now, except not for GV invariants.
--   I believe the code works for these, just never finds isomorphisms using GV invariants.
--  N examples have Cl V torsion.  Currently our methods do not handle these 
--    examples.  Actually, they do handle these.

restart
debug needsPackage "StringTorics" -- the debug is because some functions are not yet exported.
DB4 = "../Databases/cys-ntfe-h11-4.dbm"

R = ZZ[a,b,c,d]
RQ = QQ (monoid R);
(Qs, Xs) = readCYDatabase(DB4, Ring => R);

-- We collect the nontorsion, torsion, favorable, nonfavorable's.
   torsions = for k in keys Qs list (
      istor := prune coker matrix rays Qs#k;
      if not isFreeModule istor then k else continue
      )
   nonfavorables = for k in keys Qs list (
       if not isFavorable Qs#k then k else continue
      )
   #torsions == 6
   #nonfavorables == 12

   #(keys Xs) == 2014 -- this is the maximum number of inequivalent topologies

   nonfavorableQs = {796, 800, 803, 1059, 1060, 1064, 1065, 1134, 1135, 1151, 1153, 1155}
   nonfavorableXs = sort select(keys Xs, k -> member(first k, nonfavorableQs))
   nonfavorableXs = {
       (796, 0), -- potentially same as (798,1)? (which is equiv to 5 others) SAME YES the same!
       (796, 1), -- potentially same as (795,1)? (which is equiv to 5 others) -- seemingly not same
       (800, 0), -- potentially same as (844,0)? (which is equiv to 2 others) SAME 
       (803, 0), -- potentially same as (807,0)? (which is equiv to 5 others) SAME YES the same!
       (1059, 0), (1060, 0), (1064, 0), (1065, 0), -- these 4 have same c2, cubic (exactly same) 
         -- possibly same as (1071,0)? (which is itself equiv to 10 others) SAME

       (1134, 0), -- potentially same as (1136,0)? SAME
       (1135, 0), -- potentially same as (1139,0)? SAME
       (1151, 0), (1153, 0), -- potentially same each other, and as (1156,0)? (which is equiv to a bunch of others).  hessian doesn't factor...
       (1155, 0), -- UNIQUE this one is unique in its invariants bucket. bu hessian does factor.
       (1155, 1)} -- potentially same as (1168,0). hessian doesn't factor

  torsions = {0, 3, 4, 5, 12, 15}
  torsionCYs = sort select(keys Xs, k -> member(first k, torsions))
  assert(torsionCYs == {(0, 0), (3, 0), (4, 0), (5, 0), (12, 0), (15, 0)})

-- There are 6 polytopes whose toric variety has torsion class group.
-- Each of these has 1 NTFE triangulation.
-- These 6 form three topological classes.
-- (0,0) -- is alone, not equivalent to any other in Xs
-- (5,0), (4,0), (3,0) -- have equivalent c2, cubic forms,  not equivalent to any other in Xs
-- (15,0), (12,0) -- have equivalent c2, cubic forms, not equivalent to any other in Xs
-- These are distinct from all others at the invariant level.

factorsByType = method()
factorsByType RingElement := HashTable => F -> (
    facs := factors F;
    faclist := for fx in facs list (fx#0, sum first exponents fx#1, fx#1);
    H := partition(x -> {x#0, x#1}, faclist);
    hashTable for k in keys H list k => for x in H#k list x_2
    )   

factorsByDegree = method()
factorsByDegree RingElement := HashTable => F -> (
    facs := factors F;
    faclist := for fx in facs list (sum first exponents fx#1, fx#1);
    H := partition(x -> x#0, faclist);
    hashTable for k in keys H list k => for x in H#k list x_1
    )   
///
  X1 = Xs#(796,0)
  X2 = Xs#(798,1)
  F1 = cubicForm X1
  F2 = cubicForm X2
  L1 = c2Form X1
  L2 = c2Form X2

  facs1 = factorsByType det hessian F1
  facs2 = factorsByType det hessian F2
  facs1 = factorsByDegree det hessian F1
  facs2 = factorsByDegree det hessian F2

  RZ = R
  (A,phi) = genericLinearMap RQ
  T = source phi
  Ps = permutations facs2#{1,1}
  trythem = () -> (
    A0ZZ := null;
    for p in Ps do (
      (A0, phi0, J0) := linearEquationConstraints(A, phi, 
          prepend({sub(L1, T), sub(L2, T)},
            for j from 0 to #facs1#{1,1} - 1 list {sub(facs1#{1,1}#j, T), sub(p#j, T)}),
          {});
      -- if A0 is over ZZ, invertible, then we return it.
      mapshouldwork := try (
            A0ZZ = lift(A0, ZZ);
            print (A0, A0ZZ);
            abs det A0ZZ === 1
          ) else false;
      if mapshouldwork then (
          if mapIsIsomorphism(transpose A0ZZ, X1, X2) then 
              return A0ZZ;
          );
      )
  )

  trythem()      

  findEquivalenceViaHessians = method()
  findEquivalenceViaHessians(CalabiYauInToric, CalabiYauInToric) := (X1, X2) -> (
      F1 := cubicForm X1;
      F2 := cubicForm X2;
      L1 := c2Form X1;
      L2 := c2Form X2;
      facs1 := factorsByType det hessian F1;
      facs2 := factorsByType det hessian F2;
      RZ := ring F1;
      RQ := QQ (monoid RZ);
      (A,phi) := genericLinearMap RQ;
      T := source phi;
      Ps := permutations facs2#{1,1};
      possibleSigns := flatten for i1 in {1,-1} list flatten for i2 in {-1,1} list flatten for i3 in {-1,1} list for i4 in {-1,1} list {i1,i2,i3,i4};
      trythem = () -> (
          A0ZZ := null;
          for sgn in possibleSigns do
          for p in Ps do (
              (A0, phi0, J0) := linearEquationConstraints(A, phi, 
                  prepend({sub(L1, T), sub(L2, T)},
                      for j from 0 to #facs1#{1,1} - 1 list {sub(facs1#{1,1}#j, T), sgn#j * sub(p#j, T)}),
                  {});
              print A0;
              -- if A0 is over ZZ, invertible, then we return it.
              mapshouldwork := try (
                  A0ZZ = lift(A0, ZZ);
                  abs det A0ZZ === 1
                  ) else false;
              if mapshouldwork then (
                  if mapIsIsomorphism(transpose A0ZZ, X1, X2) then 
                  return A0ZZ;
                  );
              )
          );
      trythem()
      )

  findEquivalenceViaHessians(CalabiYauInToric, CalabiYauInToric) := (X1, X2) -> (
      F1 := cubicForm X1;
      F2 := cubicForm X2;
      L1 := c2Form X1;
      L2 := c2Form X2;
      facs1 := factorsByDegree det hessian F1;
      facs2 := factorsByDegree det hessian F2;
      RZ := ring F1;
      RQ := QQ (monoid RZ);
      (A,phi) := genericLinearMap RQ;
      T := source phi;
      Ps := permutations facs2#1;
      possibleSigns := flatten for i1 in {1,-1} list flatten for i2 in {-1,1} list flatten for i3 in {-1,1} list for i4 in {-1,1} list {i1,i2,i3,i4};
      trythem = () -> (
          A0ZZ := null;
          for sgn in possibleSigns do
          for p in Ps do (
              (A0, phi0, J0) = linearEquationConstraints(A, phi, 
                  prepend({sub(F1, T), sub(F2, T)},
                  prepend({sub(L1, T), sub(L2, T)},
                      for j from 0 to #facs1#1 - 1 list {sub(facs1#1#j, T), sgn#j * sub(p#j, T)})),
                  {});
              << A0 << endl << endl;
              -- if A0 is over ZZ, invertible, then we return it.
              mapshouldwork := try (
                  A0ZZ = lift(A0, ZZ);
                  abs det A0ZZ === 1
                  ) else false;
              if mapshouldwork then (
                  if mapIsIsomorphism(transpose A0ZZ, X1, X2) then 
                  return A0ZZ;
                  );
              )
          );
      trythem()
      )

findEquivalenceViaHessians(Xs#(796,0), Xs#(798,1))
findEquivalenceViaHessians(Xs#(796,1), Xs#(795,1)) -- doesn't find one? I think these cannot be equivalent.
ret = findEquivalenceViaHessians(Xs#(800,0), Xs#(844,0)) -- SAME
(803, 0)
factor det hessian cubicForm Xs#(803,0)
factor det hessian cubicForm Xs#(807,0)
ret = findEquivalenceViaHessians(Xs#(803,0), Xs#(807,0)) -- SAME
ret = findEquivalenceViaHessians(Xs#(1071,0), Xs#(1059,0)) -- SAME
       (1134, 0), -- potentially same as (1136,0)?
       
ret = findEquivalenceViaHessians(Xs#(1134,0), Xs#(1136,0)) -- SAME
ret = findEquivalenceViaHessians(Xs#(1135,0), Xs#(1139,0)) -- SAME
///
---------------------------------------------------------------
-- Next step: How many of these 2014 are distinct topologies? --
---------------------------------------------------------------
  allXs = sort keys Xs
  --allXs = torsionfrees -- these are the ones we consider
  allT = topologySet(allXs, Xs);
  info allT -- 2014 possibly different topologies
  
  allT1 = combineIfSame(allT, X -> (c2Form X, cubicForm X))

  identicals = sort first for x in allT1#"Sets" list (
      for x1 in x list if #x1 > 1 then x1 else continue
      )
  -- Question: are there any torsions or nonfavorables in here?
  -- Torsions: none on this list.
  -- Nonfavorables: (1059, 0), (1060, 0), (1064, 0), (1065, 0) are all the same (on the nose).
  
  info allT1 
  allT1#"Sets"#0/length//tally 
  netList select(allT1#"Sets"#0, x -> #x > 1) -- there are 1945 seemingly different (69 are same as some other one).
  
  elapsedTime allT2 = separateIfDifferent(allT1, invariantsAll) -- 260 seconds

  for x in allT2#"Sets" list (
      if any(flatten x, y -> member(first y, nonfavorables)) then x else continue
      )
  info allT2
  #allT2#"Sets" == 1130
  allT2#"Sets"/length//tally

  elapsedTime allT3a = combineByGV(allT2, DegreeLimit => 10); -- 220 sec
  info allT3a
    -- Total number of objects considered:         2014
    -- Number of known different topologies:       1130
    -- Maximum possible # of different topologies: 1251
    -- Largest number in one set:                  51
  elapsedTime allT3b = combineByGV(allT3a, DegreeLimit => 15); -- 
  elapsedTime allT3c = combineByGV(allT3b, DegreeLimit => 20); -- 

  elapsedTime allT3 = combineByGV(allT3c, DegreeLimit => 25); -- 
  info allT3
  -- Let it run at this point. XXX overnight run.
  -- Next steps:
  --  Are the torsions all on their own?
  --  How many nonfavorables are still around? (most, I would guess...)
  --  For those with linear hessians, we can match linear forms to linear forms (or negatives):
  --    for 4 linear factors, we have 4! * 2^3 different possible ansatz'.
  --    Use these to find (or show they do not exist) equivalences.
  --  This can probably work unless the hessian is irreducible.  What do we do then?

  elapsedTime allT3 = separateByGV allT2 -- 837 sec
  elapsedTime allT3a = separateByGV allT3 -- using degree limit 25.
  info allT3
  #allT3#"Sets" == 1130
  allT3#"Sets"/length//tally
  -- Tally{1 => 1060}
  --       2 => 58
  --       3 => 10
  --       4 => 2
  -- UPSHOT: There are at least 1130 topologies, but at most 
  -- 
  1060 + 2*58 + 3*10 + 2*4  == 1214 -- This is the current upper bound on number of toric phases.
  1060 + 58 + 10 + 2 == 1130
  
  info allT3a
  #allT3a#"Sets" == 1130
  allT3a#"Sets"/length//tally
  
  for a in allT3#"Sets" list (
      a/first
      )

  netList for x in allT3#"Sets" list (
      if any(flatten x, y -> member(first y, nonfavorables)) then x else continue
      )

  netList for x in allT3#"Sets" list (
      if any(flatten x, y -> member(first y, torsions)) then x else continue
      )

-- preliminary numbers:
--  #polytopes
--  #CY's (NTFE)
--  # nonfavorables not equivalent to favorable
--  same for torsions.
--  

for x in allT3#"Sets" list (
    for x1 in x list (
        print x1;
        if any(x1, x2 -> member(first x2, nonfavorables))
        then x1 else continue
        )
    )

netList for x in allT3a#"Sets" list (
    x' := for x1 in x list (
        if any(x1, x2 -> member(first x2, nonfavorables))
        then x1 else continue
        );
    if #x' == 0 then continue else x'
    )

for x in allT3a#"Sets" list (
    x' := for x1 in x list (
        if any(x1, x2 -> member(first x2, torsions))
        then x1 else continue
        );
    if #x' == 0 then continue else x'
    )
netList oo

    
    -- these are the possible pairs that could be equivalent.
    -- if the GV code is giving the correct answer (not clear, as there is a degree bound given)  
    -- then these should all not have any solution over the integeres.
    onestocheck = flatten for x in allT3#"Sets" list if #x == 1 then continue else (
        subsets(x/first, 2)
        )
    #onestocheck == 100 -- this is the number of things still needed to check.
      -- if both are favorable, it is unlikely to give the same topology.
      
    for x in onestocheck list (isFavorable cyPolytope Xs#(x#0), isFavorable cyPolytope Xs#(x#1))
    

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

kk = ZZ/32003
-- (1135,0) is nonfavorable, (1139,0) is favorable, and they are the same topology
(J, A) = getEquivalenceIdeal((1135,0),(1139,0),Xs)

(lab1, l1b2) = ((1135, 0), (1139, 0)) -- same
(lab1, lab2) = ((800,0), (846,0)) -- same
(lab1, lab2) = ((796, 1), (823, 3)) -- same
(lab1, lab2) = ((337, 0), (339, 2))  -- not done yet: probably different.

(lab1, lab2) = toSequence onestocheck#0 -- exists over ZZ[1/2]
(lab1, lab2) = toSequence onestocheck#1 -- EQUAL
(lab1, lab2) = toSequence onestocheck#2 -- probably different
(lab1, lab2) = toSequence onestocheck#3 -- 


o271_0
o271_1
o271_3

(J, A) = getEquivalenceIdeal(lab1, lab2, Xs);
-- This shows they have the same intersection and c2 forms.
(L1,F1) = (c2Form Xs#lab1, cubicForm Xs#lab1)
(L2,F2) = (c2Form Xs#lab2, cubicForm Xs#lab2)

Tp = kk (monoid ring J)
Jp = sub(J, Tp);
gbJp = groebnerBasis(Jp, Strategy => "F4");
--decompose ideal gbJp
gbJ = ideal sub(gbJp, ring J); -- maybe not quite gb of J
see gbJ
use ring A
A0 = A % (gbJ + ideal(t_(4,4)))
A0 = A % (gbJ)
phi = (map(R, R, transpose lift(A0, ZZ))) 
phi F1 - F2 == 0
phi L1 - L2 == 0
det A0

gv1 = partitionGVConeByGV(Xs#lab1, DegreeLimit => 20)
gv2 = partitionGVConeByGV(Xs#lab2, DegreeLimit => 20)
options separateByGV
methods separateByGV
code 0
findLinearMaps(gv1, gv2)
tryAtLargerDegreeLimit = for x in onestocheck list if isFavorable Qs#(first x#0) and isFavorable Qs#(first x#1) then x else continue
for x in onestocheck list if isFavorable Qs#(first x#0) and isFavorable Qs#(first x#1) then continue else x

for a in tryAtLargerDegreeLimit list (
    (lab1, lab2) = toSequence a;
    gv1 = partitionGVConeByGV(Xs#lab1, DegreeLimit => 20);
    gv2 = partitionGVConeByGV(Xs#lab2, DegreeLimit => 20);
    ans := findLinearMaps(gv1, gv2);
    << "#### " << a << " " << ans << endl;
    a => ans
    )

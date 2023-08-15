-- We compute all of the different topologies for h11=3 toric CY3 hypersurfaces
-- Caveats:
--  1 example is non favorable: #232
--  N examples have Cl V torsion.  Currently our methods do not handle these 
--    examples.

restart
debug needsPackage "StringTorics" -- the debug is because some functions are not yet exported.
  -- list of these functions:
  
--DB3 = "../Databases/test2-cys-ntfe-h11-3.dbm"
--DB3 = "../StringTorics/Databases/test2-cys-ntfe-h11-3.dbm"
DB3 = "../Databases/cys-ntfe-h11-3.dbm"

R = ZZ[a,b,c]
RZ = R
RQ = QQ (monoid R);
(Qs, Xs) = readCYDatabase(DB3, Ring => R);

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
assert(torsionCYs == {(0, 0), (9, 0), (10, 0), (55, 0), (62, 0)})
assert(nonfavorableCYs == {(232,0)})

-- how many of these have nonfavorable dual polytopes
-- these are the ones that are not nec general in moduli
select(sort keys Qs, lab -> (ans := not isFavorable polar Qs#lab; print ans; ans))

---------------------------------------------------------------
-- Next step: How many of these 306 are distinct topologies? --
---------------------------------------------------------------
  allXs = sort keys Xs
  allT = topologySet(allXs, Xs);
  info allT -- 306 possibly different topologies

  -- allT = separateIfDifferent(allT, invariantsH11H12)
  -- info allT

  -- allT = separateIfDifferent(allT, hubschInvariants)
  -- info oo

  -- PC = pointCounter RZ;
  -- allT = separateIfDifferent(allT, pointCounts_PC)
  -- info allT

  allT1 = combineIfSame(allT, X -> (c2Form X, cubicForm X))
  equivalences allT1
  info allT1 

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

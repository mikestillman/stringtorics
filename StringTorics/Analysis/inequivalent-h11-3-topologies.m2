-- We compute all of the different topologies for h11=3 toric CY3 hypersurfaces
-- Caveats:
--  1 example is non favorable: #232
--  N examples have Cl V torsion.  Currently our methods do not handle these 
--    examples.

restart
debug needsPackage "StringTorics" -- the debug is because some functions are not yet exported.
  -- list of these functions:
  
DB3 = "../Databases/cys-ntfe-h11-3.dbm"

R = ZZ[a,b,c]
RQ = QQ (monoid R);
(Qs, Xs) = readCYDatabase(DB3, Ring => R);

-- We collect the nontorsion, torsion, favorable, nonfavorable's.
   torsions = for k in keys Qs list (
      istor := prune coker matrix rays Qs#k != ZZ^3;
      if istor then k else continue
      )
torsionCYs = sort select(keys Xs, lab -> member(first lab, torsions))
assert(torsionCYs == {(0, 0), (9, 0), (10, 0), (55, 0), (62, 0), (232, 0)})

-- UPSHOT: There are 6 topologies that we are not yet prepared to analyze.
-- Each of these 6 polytopes has exactly one NTFE triangulation (non-2-face-equivalent)
-- Note that (232,0) is not favorable, but is also torsion.

torsionfrees = sort toList (set keys Xs - set torsionCYs)
favorables = sort toList (set keys Xs - set {(232,0)})
assert(#torsionfrees == 300)

---------------------------------------------------------------
-- Next step: How many of these 300 are distinct topologies? --
---------------------------------------------------------------
  allXs = torsionfrees -- these are the ones we consider
  allT = topologySet(allXs, Xs);
  info allT -- 300 possibly different topologies
  
  allT1 = combineIfSame(allT, X -> (c2Form X, cubicForm X))
  info allT1 
  -- TODO: this should be a method
  allT1#"Sets"#0/length//tally -- 14 of these are identical, leaving 286 possibly different tops
  select(allT1#"Sets"#0, x -> #x > 1) -- these have the duplicates
  
  elapsedTime allT2 = separateIfDifferent(allT1, invariantsAll) -- 

  info allT2
  -- This divides the 300 (really, 286) topologies into 170 different groups, each group 
  -- consisting of CY3's with the same invariants (in invariantsAll).
  #allT2#"Sets" == 170
  allT2#"Sets"/length//tally -- 108 of these have a unique CY3 in them (so these are not equivalent to anything else)
  -- Tally{1 => 108}
  --       2 => 35
  --       3 => 19
  --       4 => 4
  --       5 => 1
  --       7 => 1
  --       10 => 1
  --       13 => 1
  
  -- The largest set has 13 potentially the same topology
  
  -- We have two ways to proceed here.

  -- VERSION #1: use GV invariants to find equivalences
  -- TODO: I think this removes the duplicates found equivalent in allT1.  Fix that.
  --   But the result still 
    elapsedTime allT3 = separateByGV allT2 -- 45 sec
    info allT3
    allT3#"Sets"/length//tally
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
    for x in onestocheck list getEquivalenceIdeal(x#0, x#1, Xs)
    -- these are all the unit ideal, so none of these have any equivalence
    -- over QQ or RR, let alone ZZ.

  -- Upshot: There are 183 different topologies for CY3 which are have h11=3 and are toric hypersurfaces (a ll Batyrev construction)

  -- VERSION #2: Use ansatz to find equivalences, DOES NOT use GV invariants
  -- This is somewhat older code, but hopefully it gives the same answer!
  -- Actually: this version can handle the torsion examples, just not the non-favorable one.
    T1s = for lab in favorables list (
        X := Xs#lab;
        lab => {c2Form X, cubicForm X, hh^(1,1) X, hh^(1,2) X}
        );
    assert(#T1s == 305)
  
    -- now, how many of these are the same?
    291 == # unique values hashTable T1s -- 286 + 5 -- the 5 are the torsion but favorables.

    Ts = T1s; -- these are all the label => topology pairs we have (291 here).
    hashTs = hashTable Ts;
    keyTs = Ts/first//sort -- 305 of these.

    elapsedTime H = partition(lab -> invariantsAll toSequence (hashTs#lab), keyTs); -- 
    #keys H == 173
    (keys H)/(k -> #H#k)//tally
    
    elapsedTime INV = for k in keys H list k => elapsedTime partitionH113sByTopology(H#k, hashTs, RQ); -- 
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


-- Upshot:
--  toric topologies from torsionfrees: 183
--  toric topologies from torsions: 3
--  toric from unfavorable: 1 (or 0 new possibly!)
--  non-toric phases: 113
-- 186 + 113 = 299.  Or maybe 300...

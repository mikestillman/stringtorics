-- We compute all of the different topologies for h11=4 toric CY3 hypersurfaces
-- Caveats:
--  non-favorables:  These are handled now, except not for GV invariants.
--   I believe the code works for these, just never finds isomorphisms using GV invariants.
--  N examples have Cl V torsion.  Currently our methods do not handle these 
--    examples.  Actually, they do handle these.

restart
debug needsPackage "StringTorics" -- the debug is because some functions are not yet exported.
DBNAME = "../Databases/cys-ntfe-h11-5.dbm"

RZ = ZZ[a,b,c,d,e]
RQ = QQ (monoid RZ);
(Qs, Xs) = readCYDatabase(DBNAME, Ring => RZ);

-- Considering invariants (not coming from GV invariants):
  allXs = sort keys Xs;
  allT = topologySet(allXs, Xs);
  info allT -- 1 bucket. (11847 possibly different topologies)

  elapsedTime allT = separateIfDifferent(allT, invariantsH11H12) -- 3 sec
  info allT -- 113 different buckets

  elapsedTime allT = separateIfDifferent(allT, hubschInvariants); -- 127 sec
  info oo -- after H11H12: 411 different buckets

  PC = pointCounter(RZ, "Primes" => {2,3,5,7,11,13}, Projective => true);
  --PC = pointCounter(RZ, "Primes" => {2,3,5,7,11,13,(2,2),(3,2),(2,3)}, Projective => true);
  elapsedTime allT1 = separateIfDifferent(allT, pointCounts_PC);
  info allT1

  -- At this point we have XXXX different buckets.

  -- By itself: gets it to 1111.
  -- After invariantsH11H12, hubschInvariants, this gives: 1114 different.

  allT = separateIfDifferent(allT, hessianInvariants)
  info oo -- goes from 1114 to 1123, but is much faster too than point counts.
  -- by itself: 100 different classes.
  -- after h11h12, hubsch: 660 different classes
  --partition(lab -> hessianInvariants Xs#lab, sort keys Xs);

  -- This gives the same benefit as hessianInvariants, if we first do: pointcounts, hubsch, h11h12.
  allT = separateIfDifferent(allT, X -> polynomialContent det hessian cubicForm X)
  info oo

  allT = separateIfDifferent(allT, cubicConductorInvariants)
  info oo -- to 1128.
  -- by itself: 645.

  allT = separateIfDifferent(allT, cubicLinearConductorInvariants) -- this is a good one!
  info oo -- to 1135.  
  -- by itself: 1049 classes


------------------------------
-- Try separating given h12 --
------------------------------
X35s = select(sort keys Xs, lab -> hh^(1,2) Xs#lab == 35)  
  allT = topologySet(X35s, Xs);
  info allT -- 1 bucket. 31 different possible.

  elapsedTime allT = separateIfDifferent(allT, hubschInvariants); -- 
  info oo -- 3 different buckets

  --PC = pointCounter(RZ, "Primes" => {2,3,5,7,11,13}, Projective => true);
  PC = pointCounter(RZ, "Primes" => {2,3,5,7,11,13,(2,2),(3,2),(2,3),(2,4),17}, Projective => true);
  elapsedTime allT1 = separateIfDifferent(allT, pointCounts_PC);
  info allT1 -- 17 different buckets
  netList representatives allT1
  allT = allT1

  elapsedTime allT = separateIfDifferent(allT, singularContents)
  info allT -- still 17 -- no change.

  elapsedTime allT = separateIfDifferent(allT, hessianInvariants)
  info allT -- still 17 -- no change.

  for S in allT#"Sets" list if #S == 1 then continue else (
      partition(lab -> hessianInvariants Xs#lab, S/first)
      )
netList oo

  netList representatives allT  
  toricMoriConeCap Xs#(9,0)
  toricMoriConeCap Xs#(19,2)
  heft  Xs#(9,0)
  heft  Xs#(19,2)
-- We collect the nontorsion, torsion, favorable, nonfavorable's.
   torsions = for k in keys Qs list (
      istor := prune coker matrix rays Qs#k;
      if not isFreeModule istor then k else continue
      )
   nonfavorables = for k in keys Qs list (
       if not isFavorable Qs#k then k else continue
      )
   favorablesXs = sort for k in keys Xs list (
       if not isFavorable Qs#(first k) then continue else k
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


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
elapsedTime (Qs, Xs) = readCYDatabase(DBNAME, Ring => RZ);

-- Considering invariants (not coming from GV invariants):
  allXs = sort keys Xs;
  
  nonfavorables = for k in keys Qs list (
      if not isFavorable Qs#k then k else continue
      )
  
  nonfavorableXs = sort select(keys Xs, lab -> member(first lab, nonfavorables))
  #nonfavorables == 93
  #nonfavorableXs == 134

  favorableXs = sort for k in keys Xs list if isFavorable Qs#(first k) then k else continue;
  #favorableXs == 11713

  allT = topologySet(allXs, Xs);
  allT = topologySet(favorableXs, Xs);
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


----------------------------------
-- Try separating given h12 = 35--
----------------------------------
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

  flatten flatten select(allT#"Sets", s -> #s == 8)
  set8 = {(8, 0), (11, 0), (12, 0), (12, 1), (12, 2), (21, 0), (21, 1), (21, 2)}
  set8/(lab -> isFavorable Xs#lab)
  set8/(lab -> isFavorable polar Qs#(first lab))
  gvs = hashTable(set8/(lab -> lab => elapsedTime gvInvariants(Xs#lab, DegreeLimit => 16)));
  for lab in set8 list (
      X := Xs#lab;
      degvec := heft X;
      print degvec;
      for c in toricMoriConeCap X list classifyExtremalCurve(gvs#lab, c, 16, degvec)
      )

  transpose matrix toricMoriConeCap Xs#(8,0)
  transpose matrix toricMoriConeCap Xs#(11,0)
  f = (lab) -> (X := Xs#lab; degvec := heft X; for c in toricMoriConeCap X list classifyExtremalCurve(

  C120 = transpose matrix toricMoriConeCap Xs#(12,0)
  C121 = transpose matrix toricMoriConeCap Xs#(12,1)
  A = C120^-1 * C121_{1,0,2,3,4}
  A = C120^-1 * C121_{1,0,3,2,4}
  A = C120^-1 * C121_{0,1,3,2,4}
  phi = map(RZ, RZ, A)
  phi c2Form Xs#(12,0)
  c2Form Xs#(12,1)
  partitionByTopology(set8, Xs, 16)
  
  set8/(lab -> partitionGVConeByGV(Xs#lab, DegreeLimit => 5))
  
-- Here is the group of 66 that might still be equivalent after point counts ---
-- And the second largest group: of 30.
L = value get "allT-h11-5-after-pointCounts";
netList select(L, x -> #x > 9)
for x in L list {
    M := flatten x;
    M = select(M, lab -> isFavorable Xs#lab);
    allLi = topologySet(M, Xs);
    combineByGV(allLi, DegreeLimit => 10);
    print info allLi
    }

set66 = {(1817, 0), (1818, 0), (1821, 0), (1821, 1), (1821, 2), (1821, 3), (1821, 4), (1821, 5), (1821, 6), (1821, 7), (1821, 8), (1821, 9), (1829, 0), (1829, 1), (1829, 2), (1829, 3), (1829, 4), (1838, 0), (1838, 1), (1838, 2), (1838, 3), (1838, 4), (1838, 5), (1838, 6), (1838, 7), (1838, 8), (1838, 9), (1838, 10), (1838, 11), (1838, 12), (1838, 13), (1838, 14), (1839, 0), (1839, 1), (1839, 2), (1839, 3), (1839, 4), (1839, 5), (1839, 6), (1839, 7), (1839, 8), (1839, 9), (1843, 0), (1843, 1), (1843, 2), (1843, 3), (1843, 4), (1843, 5), (1843, 6), (1843, 7), (1843, 8), (1843, 9), (1849, 0), (1849, 1), (1849, 2), (1858, 0), (1858, 1), (1858, 2), (1858, 3), (1923, 0), (1923, 1), (1923, 2), (1923, 3), (1923, 4), (1923, 5), (2009, 0)}  
set30 = {(4425, 0), (4426, 0), (4428, 0), (4429, 0), (4438, 0), (4439, 0), (4440, 0), (4442, 0), (4444, 0), (4445, 0), (4446, 0), (4450, 0), (4454, 0), (4460, 0), (4461, 0), (4464, 0), (4465, 0), (4471, 0), (4476, 0), (4477, 0), (4478, 0), (4479, 0), (4481, 0), (4485, 0), (4488, 0), (4498, 0), (4499, 0), (4503, 0), (4504, 0), (4506, 0)}
partitionByTopology(set66, Xs, 5) -- these 66 are all the same.
partitionByTopology(set30, Xs, 5) -- 19 different classes left

-- 8 of the 30 are not favorable, 22 are.  The 22 are 
all30 = topologySet(set30, Xs);
  elapsedTime all30a = combineByGV(all30, DegreeLimit => 10); -- sec
  elapsedTime all30b = combineByGV(all30a, DegreeLimit => 15); -- sec
  --elapsedTime all30c = combineByGV(all30b, DegreeLimit => 18); -- didn't seem to work?!

  -- 9 different ones left.
  -- can we separate them by various invariants?
  separateIfDifferent(all30c, invariantsAll) -- didn't separate!
  -- let's do more point counts.
  
  elapsedTime PC = pointCounter(RZ, "Primes" => {(2,2),(3,2),(2,3), 17, 19, 23}, Projective => true);
  elapsedTime all30d = separateIfDifferent(all30b, pointCounts_PC);
  info all30d -- still 9 different classes!
  
  the9 = flatten representatives all30d
  toricMoriConeCap(Xs#(the9_0))

  set22 = select(set30, lab -> isFavorable Xs#lab)
  all22 = topologySet(set22, Xs)
  info all22
  elapsedTime all22a = combineByGV(all22, DegreeLimit => 5); -- sec
  info all22a -- <= 11 different.
  elapsedTime all22b = combineByGV(all22a, DegreeLimit => 10); -- sec
  info all22b  -- <= 2 different.
  elapsedTime all22c = combineByGV(all22b, DegreeLimit => 15); -- sec
  info all22c -- 1 DIFFERENT TOP HERE.

-- read in the different sets after point counts, from heaviside.
allS = value get "allT-h11-5-after-pointCounts";
allS = select(allS, x -> #x > 1)
tally (allS/(x -> #x)) -- 1662 pairs of 2 each.
-- Let's attack those.
allS2 = select(allS, x -> #x == 2)
gv90 = gvInvariants(Xs#(9,0), DegreeLimit => 16)
gv192 = gvInvariants(Xs#(19,2), DegreeLimit => 16)

degvec = heft Xs#(9,0)
for c in toricMoriConeCap Xs#(9,0) list dotProduct(c, degvec)
mori = posHull transpose matrix toricMoriConeCap Xs#(9,0)
matrix{hilbertBasis dualCone mori}
(transpose rays mori) * matrix{hilbertBasis dualCone mori}
classifyExtremalCurves(gv90, toricMoriConeCap Xs#(9,0), 16, heft Xs#(9,0))
classifyExtremalCurves(gv192, toricMoriConeCap Xs#(19,2), 16, heft Xs#(19,2))
toricMoriConeCap Xs#(19,2)

--------------- BELOW THIS IS PROBABLY FROM h11=4... ------------------------------
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


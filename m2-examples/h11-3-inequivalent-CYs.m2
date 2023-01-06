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

  -- so FOR NOW, we ignore the 6 X's which are in the nontorsionfree list.
  
    elapsedTime H = partition(lab -> invariantsAll Xs#lab, torsionfrees); -- 19 seconds
    #(keys H) == 195 -- this is the minimum number of different topologies.
    (keys H)/(k -> #H#k)//tally
    elapsedTime INV = for k in keys H list k => partitionByTopology(H#k, Xs, 15); -- 34 seconds
    INV/(x -> #x#1)//tally
    -- so the maximum number of distinct topologies is 186 + 16 + 3 == 205

  REPS = for k in INV list (
      H := last k; -- a hash table
      K := keys H;
      if #K  == 1 then continue
      else sort K)
  netList REPS
----------------------------------------------------------------------------
-- Step 4.  Analyze these ones that are possibly still equivalent-----------
----------------------------------------------------------------------------


  X1 = Xs#(227,0)
  X2 = Xs#(229,0)
  F1 = cubicForm X1
  F2 = cubicForm X2
  factor det hessian F1
  factor det hessian F2
  
  R = QQ[a,b,c]
  (A, phi) = genericLinearMap R
  TR = target phi
  T = coefficientRing TR
  phi
  A
  assert(source phi === TR)
  assert(ring A === coefficientRing TR)

    linearEquationConstraints(A, phi, {
            {a,b}
            }, {})


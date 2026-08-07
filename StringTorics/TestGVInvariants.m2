-*
  restart
  needsPackage "StringTorics"
*-
TEST /// -- WORKING ON THIS ONE, AND OTHER GV invariants code
  debug needsPackage "StringTorics" -- displayRays
  Q = reflexivePolytope(
      {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
       {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
      ID => 7)
  R = ZZ[a,b,c];
  X = makeCY(Q, PicardRing => R, ID => 0);
  elapsedTime GVT = gvInvariants(X, DegreeLimit => 12);
  GVT
  rays GVT
  displayRays GVT
  displayRays(GVT, "OneOnly" => true)

  C = gvCone GVT -- gives cone, and the weight vector non-negative on this cone.
  rays C

  curveGVs = extremalCurves GVT
  curveGVs2 = extremalCurves(X, entries transpose rays C)
  assert(curveGVs === curveGVs2)

  netList nilpotentCurves gvRayTable GVT

  partitionGVConeByGV(X, DegreeLimit => 12)
  partitionGVConeByGV(X, DegreeLimit => 7)
///

-*
  restart
  needsPackage "StringTorics"
*-
TEST /// -- tests flop code that is in GVInvariants (i.e. doesn't refer to CY3 directly).
  Q = reflexivePolytope(
      {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
       {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
      ID => 7)
  R = ZZ[a,b,c];
  X = makeCY(Q, PicardRing => R, ID => 0);
  elapsedTime GVT = gvInvariants(X, DegreeLimit => 12);

  -- gvRayTable looks good.
  GVR = gvRayTable GVT
  assert instance(GVR, GVRayTable)
  hashTable gvRayTable GVT
  assert(gvRayTable(GVT, set{}) === gvRayTable GVT)
  gvRayTable(GVT, set{{-1,1,0}})

  -- gvCone using gvRayTable, looks good
  C0 = gvCone GVT
  C1 = gvCone(GVT, set {})
  C2 = gvCone(GVT, set{{-1,1,0}})
  C1' = gvCone(GVR, degreeLimit GVT)
  assert(C0 == C1)
  assert(C1' == C1)
  rays C1
  rays C2

  -- extremalCurves, looks good.
  extremalCurves(GVT)
  extremalCurves(GVT, set{})
  extremalCurves(GVT, set{{-1,1,0}})
  extremalCurves(GVR, degreeLimit GVT)

  C1 = moriCone GVT -- doesn't exist
  rays (C1' = moriCone(GVT, set{})) -- OK
  C2' = moriCone(GVT, set{{-1,1,0}})
  assert(C2 == C2')
  assert(C1 == C1')

  assert(
      set entries transpose rays moriCone GVT ==
      set {{-1, 1, 0}, {1, 2, -1}, {0, -2, 1}}
      )
///

-*
  restart
  needsPackage "StringTorics"
*-
TEST ///
  Q = reflexivePolytope(
      {{-1,-1,-1,-1},{-1,-1,-1,0},{-1,-1,0,2},
       {-1,0,-1,-1},{0,-1,-1,-1},{1,-1,0,-1},{1,2,2,2}},
      ID => 7)
  R = ZZ[a,b,c];
  X = makeCY(Q, PicardRing => R, ID => 0);
  elapsedTime GVT = gvInvariants(X, DegreeLimit => 12);

  X1 = makeCY3(X, DegreeLimit => 16)
  assert instance(gvTable X1, GVTable)
  extremalCurves X1
  rays gvCone X1
  rays moriCone X1
  -- Now let's do a flop.  There are 2 that can be done from X: {-1,1,0}, and {1,2,-1}.

  X2 = performFlop(X1, {-1,1,0})
  moriCone X2
  rays oo
  extremalCurves X2

  X3 = performFlop(X1, {1,2,-1})
  extremalCurves X3

  extremalCurves X2
  X4 = performFlop(X2, {-1,-1,1})
  extremalCurves X4
  
  extremalCurves X2
  X5 = performFlop(X2, {0,3,-1})
  extremalCurves X5

  extremalCurves X3
  X6 = performFlop(X3, {-1,-2,1})
  extremalCurves X6

  extremalCurves X4

  X6 = performFlop(X3, {-1,-2,1})
  extremalCurves X6

  
  makeCY3(X1, DegreeLimit => 20) -- can't change the degree yet...
  extremalCurves X1
  flops X1

  
  
  X1 = makeCY3(X, GVT)
  X1 = makeCY3(X, GVTable => GVT) -- X a CalabiYauInToric
  X1 = makeCY3(X, DegreeLimit => 16) -- X a CalabiYauInToric OR a CY3 (either creates a GVTable, or replaces it).
  rays moriCone X1 -- works for no flops, at least.
  
  flops X1 -- list of (extremal) curve classes that are flops curves.
  extremalCurves gvTable X1 -- {curve class, deg, gvs, classificaiton}
--  nilpotentCurves gvTable X1
    
  performFlop X1

  nilpotentCurves X1 -- needed?

  -- apply a invertible integral map to c2 form, cubic form, and the curve classes in GVT.
  changeBasis(X1, A) -- A is an n x n invertible integral matrix, n = h^(1,1) X.

  curveGVs = extremalCurves GVT
  curveGVs2 = extremalCurves(X, entries transpose rays C)
  assert(curveGVs === curveGVs2)

  netList nilpotentCurves GVT

///

/// -- test: which X of h11=3 have seemingly infinitelg generated Mori cone?  effective cone?
  restart
  debug needsPackage "StringTorics"
  databaseLOC -- ~/utah/CYToolsM2/Databases, in my init.m2 file
  DBNAME = databaseLOC | "/cy3-h11-3.dbm"
  R = ZZ[a,b,c]
  RQ = QQ (monoid R);
  (Qs, Xs) = readCYDatabase(DBNAME, Ring => R);
  allXs = sort keys Xs

  lab = allXs#19 -- allXs#18: likely infinite flop chains: 17, 18, (not 15, 16).
  X = Xs#lab 
  heft X
  GVT = gvInvariants(X, DegreeLimit => 40);
  displayRays(GVT, "OneOnly" => true)
  -- rays gvCone GVT
  extremalCurves GVT
  netList nilpotentCurves GVT
///



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
  elapsedTime GVT = gvInvariantsNew(X, DegreeLimit => 12);
  GVT
  rays GVT
  displayRays GVT
  displayRays(GVT, "OneOnly" => true)

  C = gvCone GVT -- gives cone, and the weight vector non-negative on this cone.
  rays C

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
  GVT = gvInvariantsNew(X, DegreeLimit => 40);
  displayRays(GVT, "OneOnly" => true)
  -- rays gvCone GVT
  extremalCurves GVT
  netList nilpotentCurves GVT
///



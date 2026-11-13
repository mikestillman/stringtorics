-*
  restart
  needsPackage "StringTorics"
*-
/// -- MES: In process (chosen by opus 5.5) -- test of: simplexes and polar dual simplexes
  debug needsPackage "StringTorics"

  simplexVerts = d -> append(entries id_(ZZ^d), toList(d:-1));
  dualSimplexVerts = d -> entries transpose lift(vertices polar convexHull transpose matrix simplexVerts d, ZZ)

  -- Creation of CYPolytopeNew's for simplices and dual simplices.
  set1 = hashTable for d from 1 to 9 list elapsedTime d => simplexVerts d
  set2 = hashTable for d from 1 to 9 list elapsedTime d => dualSimplexVerts d

  simp1 = hashTable for d from 1 to 9 list elapsedTime d => cyPolytopeNew set1#d
  simp2 = hashTable for d from 1 to 9 list elapsedTime d => cyPolytopeNew set2#d

  for d from 1 to 9 list elapsedTime vertices simp1#d
  for d from 1 to 9 list elapsedTime vertices simp2#d

  for d from 1 to 9 list elapsedTime halfspaces simp1#d
  for d from 1 to 9 list elapsedTime halfspaces simp2#d

  for d from 1 to 9 list elapsedTime dim simp1#d
  for d from 1 to 9 list elapsedTime dim simp2#d
  
  
  -- the old `latticePoints` code is way slow: `# latticePoints set2#8` takes 225 sec.
  -- before changing old lattice point code.
  for d from 1 to 9 list elapsedTime numcols latticePoints simp1#d
  for d from 1 to 9 list elapsedTime numcols latticePoints simp2#d

  -- This one is pretty slow...
  for d from 8 to 8 list elapsedTime latticePointListNew simp1#d
  for d from 8 to 8 list elapsedTime latticePointListNew simp2#d;

  polyt1 = hashTable for d from 1 to 9 list d => elapsedTime polytope simp1#d
  polyt2 = hashTable for d from 1 to 9 list d => elapsedTime polytope simp2#d

  -- list of all faces isn't fast, but at d=9, still under a second (.83 sec simp1#9, .89 sec sim2#9))
  polyt1 = hashTable for d from 1 to 9 list d => elapsedTime faces polytope simp1#d;
  polyt2 = hashTable for d from 1 to 9 list d => elapsedTime faces polytope simp2#d;

  -- the main time in heere appears to be in `faces Polyhedron`.  Can we speed that up?
  -- TODO: I want this associated to the CYPolytopeNew, I think...  for example, for:
  -- BUG: these are not necessarily for the order of vertices in the CYPolytopeNew
  hashTable for d from 1 to 9 list d => elapsedTime faceDimensionHash polytope simp1#d;
  hashTable for d from 1 to 9 list d => elapsedTime faceDimensionHash polytope simp2#d;

  -- the following is just an indication for amount of time it takes to do `vertices` for Polyhedron.
  polytp1 = hashTable for d from 1 to 9 list d => elapsedTime polytope(simp1#d, "M")
  polytp2 = hashTable for d from 1 to 9 list d => elapsedTime polytope(simp2#d, "M")
  for d from 1 to 9 list elapsedTime vertices polytp1#d
  for d from 1 to 9 list elapsedTime vertices polytp2#d

  faceDimensions simp2#3
  latticePoints simp2#3
  -- The key thing to speed up, or determine exactly what we need: where the lattice points sit.
  -- This seems to take real time.

  -- automorphisms.  This is slow in higher dimension (simp1#6: 2.5 sec)
  -- TODO: speed this up?  Do we need it for larger dimensions?
  for d from 1 to 9 list print (elapsedTime # automorphisms simp1#d);
  for d from 1 to 6 list print (elapsedTime # automorphisms simp2#d); -- simp2#5: 2.7 sec

  -- elapsedTime boundedLatticePoints(A,b, Bound => max flatten entries vP, Limit => 10000000, NodeLimit => 10000000);

///


-*
  restart
  needsPackage "StringTorics"
*-
/// -- test of larger examples in KreuzerSkarke, dim=3
  debug needsPackage "StringTorics"

  eg1 = KSEntry "4 9  M:68 9 N:341 9 H:250,46 [408] id:96
    1    0    0    0  -21  -39  -43  -49  -63
    0    1    0    0  -13  -26  -28  -32  -42
    0    0    1    0   -6  -10  -12  -14  -18
    0    0    0    1   -1   -2   -2   -2   -2"
  eg2 = KSEntry "4 10  M:68 10 N:342 10 H:250,46 [408] id:97
    1    0    0    0  -19   -1   -3   -4   -4  -25
    0    1    0    0  -13    0   -1   -1   -2  -16
    0    0    1    0   -5    0   -1   -2   -1   -8
    0    0    0    1   -1    1    2    2    2    2"
  eg3 = KSEntry "4 11  M:68 11 N:346 11 H:250,46 [408] id:98
    1    0    0    0  -19   -1  -21   -3   -5   -8  -17
    0    1    0    0  -13    0  -14   -1   -2   -4  -10
    0    0    1    0   -5    0   -6   -1   -2   -3   -6
    0    0    0    1   -1    1    0    2    2    4    4"
  eg4 = KSEntry "4 11  M:68 11 N:346 11 H:250,46 [408] id:99
    1    0    0    0  -19   -1  -21   -3   -5   -6  -16
    0    1    0    0  -13    0  -14   -1   -2   -3   -9
    0    0    1    0   -5    0   -6   -1   -2   -2   -6
    0    0    0    1   -1    1    0    2    2    3    5"

  eg300a = KSEntry "4 9  M:26 9 N:415 10 H:300,14 [572] id:0
    1    0    0   -7    0  -23   -1   -3   -7
    0    1    0   -4    0  -15   -1   -2   -5
    0    0    1   -2    0   -6    0   -1   -2
    0    0    0    0    1   -1    1    1    1"
  eg300b = KSEntry "4 9  M:27 9 N:411 9 H:300,14 [572] id:1
    1    0    0   -3    0  -23   -1   -3   -7
    0    1    0   -1    0  -15   -1   -2   -5
    0    0    1   -1    0   -6    0   -1   -2
    0    0    0    0    1   -1    1    1    1"

  eg491 = KSEntry "4 5  M:26 5 N:680 5 H:491,11 [960] id:0
    1    0    0    0  -42
    0    1    0    0  -28
    0    0    1    0  -12
    0    0    0    1   -1"

  set1 = {eg1, eg2, eg3, eg4, eg300a, eg300b, eg491}

  polys = for ks in set1 list elapsedTime cyPolytopeNew ks;

  for p in polys list elapsedTime vertices p
  for p in polys list elapsedTime polar p;
  for p in polys list elapsedTime dim p;
  for p in polys list elapsedTime numcols latticePoints p
  for p in polys list elapsedTime faces polytope p; -- < .02 each, not so bad.
  for p in polys list elapsedTime faces polar polytope p; -- < .09 sec each
  for p in polys list elapsedTime annotatedFaces p; -- expensive: over .33 seconds each
  for p in polys list elapsedTime hh^(2,1) p -- not fast, almost all of cost is in annotatedFaces.
  for p in polys list elapsedTime # automorphisms p -- none have automorphisms, I guess.  When alot, it takes a lot longer.
  for p in polys list elapsedTime isFavorable p -- 
  for p in polys list elapsedTime isFavorable polar p -- 
  
///

-*
  restart
  needsPackage "StringTorics"
*-
"TEST" -- add back in when CYPolytopeNew is incorporated.
/// -- test of weight systems of larger dimension
  debug needsPackage "StringTorics"
  eg1 = {8,1,1,1,1,2,2}
  eg2 = {40,5,6,7,7,7,8}
  set5 = {eg1, eg2}

  polys5 = for ws in set5 list elapsedTime cyPolytopeNew palpNVertices ws;

  for p in polys5 list elapsedTime dim p
  for p in polys5 list elapsedTime vertices p
  for p in polys5 list elapsedTime polar p;

  for p in polys5 list elapsedTime numcols latticePoints p
  for p in polys5 list elapsedTime faces polytope p; -- < .09 each, not so bad yet
  for p in polys5 list elapsedTime faces polar polytope p; -- < .09 sec each
  for p in polys5 list elapsedTime annotatedFaces p; -- expensive: over .33 seconds each

  for p in polys5 list elapsedTime faceDimensionHash polytope p; -- < .09 each, not so bad yet
  for p in polys5 list elapsedTime faceDimensions p; 


  for p in polys5 list elapsedTime hh^(2,1) p -- not fast, almost all of cost is in annotatedFaces.
  for p in polys5 list elapsedTime # automorphisms p -- none have automorphisms, I guess.  When alot, it takes a lot longer.
  for p in polys5 list elapsedTime isFavorable p -- 
  for p in polys5 list elapsedTime isFavorable polar p -- 

  ws = weightSystem{50, 1, 1, 1, 2, 4, 16, 25}
  palpNefPartitions(ws, 3)
  elapsedTime wss = weightSystems(6, Degrees => 50, Reflexive => true); -- 3.3 sec, 874 of these.
  for w in wss list (
      elapsedTime n := palpNefPartitions(w, 3);
      if #n#"Partitions" > 0 then print n#"Partitions";
      w => n)
  Q = cyPolytopeNew palpNVertices ws

  elapsedTime dim Q
  elapsedTime vertices Q
  Q' = elapsedTime polar Q

  elapsedTime latticePointList polytope Q
  elapsedTime latticePointList polytope Q'; -- 62495 lattice points
  elapsedTime latticePoints Q';

  vQ = matrix vertices Q
  vQ' = transpose matrix vertices Q'  
  vQ * vQ'  
  elapsedTime annotatedFaces Q -- quite long, due to latticePoints...
  


  -- bench("ws 5D {1,1,1,1,2,2}", () -> reflexivePolytope palpNVertices {8,1,1,1,1,2,2});
  -- bench("ws 5D {40,5,6,7,7,7,8}", () -> reflexivePolytope palpNVertices {40,5,6,7,7,7,8});
///



-*
  restart
  needsPackage "StringTorics"
*-
/// -- test of: reflexivePolytope, polar.
  -- this example is from this:
  tope = KSEntry "4 6  M:57 6 N:9 8 H:4,52 [-96] id:40
    1   -1    1    1    1   -3
    0    0    2    2    2   -6
    0    0    0    4    0   -4
    0    0    0    0    4   -4
  "

  -- method #1 from a list
  verts = {{-1, 0, 0, 0}, {-1, 0, 0, 1}, {-1, 0, 1, 0},
      {-1, 2, -1, -1}, {1, -1, 0, 0}, {1, -1, 0, 1},
      {1, -1, 1, 0}, {1, 1, -1, -1}}
  Q = reflexivePolytope(verts, ID => 40) -- from h11=4, id=40 Kreuzer Skarke
  assert isReflexive polytope Q
  assert isWellDefined Q

  -- method #2 from a matrix
  Q2 = reflexivePolytope(transpose matrix verts, ID => 40)
  assert isWellDefined Q2
  assert(Q2 === Q)
  Q2.cache === Q.cache -- this should be false?  Here is the actual check:
  assert not (hashTable pairs Q2.cache === hashTable pairs Q.cache)
  assert(vertices Q === vertices Q2)

  -- method #3: using polar
  Q' = polar Q -- It should be that Q stashes Q' (and vice versa?)
  Q3 = polar Q'
  isWellDefined Q3
  assert(Q3 === Q)
  hashTable pairs Q.cache === hashTable pairs Q3.cache
  -- TODO: Q' only stashes the polytope object, not the ReflexivePolytope object.
  --   fix this.

  -- method #4: using a Polyhedron object
  P4 = convexHull transpose matrix verts
  Q4 = reflexivePolytope(P4, ID => 40)
  isWellDefined Q4
  assert(Q4 === Q)

  -- method #5: Using Kreuzer-Skarke database
  Q5 = reflexivePolytope tope -- takes ID from tope.
  isWellDefined Q5
  assert(Q5 === Q)
///

-*
  restart
  needsPackage "StringTorics"
*-
/// -- test of: dump, reflexivePolytope String (inverses of each other
  -- this example is from this:
  verts = {{-1, 0, 0, 0}, {-1, 0, 0, 1}, {-1, 0, 1, 0},
      {-1, 2, -1, -1}, {1, -1, 0, 0}, {1, -1, 0, 1},
      {1, -1, 1, 0}, {1, 1, -1, -1}}
  Q = reflexivePolytope(verts, ID => 40) -- from h11=4, id=40 Kreuzer Skarke
  assert isWellDefined Q

  computeBasics Q
  str = dump Q
  Q2 = reflexivePolytope str
  assert(Q === Q2)
  str2 = dump Q2
  assert(str === str2)
///  

-*
  restart
  needsPackage "StringTorics"
*-
/// -- test of: findAllFRVTs
  topes = kreuzerSkarke 5;
  tope = topes_101
  Q = reflexivePolytope tope
  tri = findOneFRST Q
  computeBasics Q
  elapsedTime # findAllFRSTs Q
  elapsedTime # findAllFRVTs Q
  tris = findAllFRVTs Q
  assert all(for t in tris list (
      V := normalToricVariety(rays Q, t);
      result := isWellDefined V;
      << result << endl;
      result
      ))
///

-*
  restart
  needsPackage "StringTorics"
*-
///
  --topes = kreuzerSkarke(5, Limit => 100)
  --tope = topes_55
  tope = KSEntry "4 11  M:48 11 N:10 8 H:5,41 [-72] id:55
   1   0   0   0  -3   3   1  -3  -4  -1  -4
   0   1   0   1   2  -4   0   2   3   0   3
   0   0   1   0  -1   2   0  -2  -2  -2  -2
   0   0   0   3   3  -6   2   2   4   2   5
  "
  Q = reflexivePolytope tope

  tri = findOneFRST Q

  --findOneFRVT Q -- WRONG!!

  elapsedTime ptris = findAllFRSTs Q; -- relatively expensive
  assert(#ptris == 18)
  for tri in ptris list (
      V := normalToricVariety(rays Q, tri);
      isWellDefined V and isProjective V
      )
  assert all(oo, x  -> x)
  
  vtris = sort findAllFRVTs Q;
  for tri in vtris list (
      V := normalToricVariety(rays Q, tri);
      isWellDefined V and isProjective V
      )
  tally oo -- 18 false, 121 true.
  assert all(oo, x  -> x)
  assert(#vtris == 103)
  assert(#ptris == 18)

  ---------------------------------------------------------
  -- Testing with vector triangulations
  -- Now let's try to maneuver in the space of triangulations
  -- starting with tri above.
  cQ = chirotope(transpose matrix rays Q, Homogenize => false)
  "topcomfoo.in" << toString cQ << endl << close
  ctris = for L in lines get ("!chiro2finetriangs <topcomfoo.in") list (
      matches := regex("\\{\\{[0-9,\\{\\}]*", L);
      if matches === null or #matches != 1 then error "my logic is wrong";
      value substring(matches_0, L)
      );
  ctris = sort for x in ctris list (x/sort//sort);
  #ctris == 139
  toList(set vtris + set ptris - set ctris) -- every one in vtris and ptris is in ctris.
  extras = toList (set ctris - (set vtris + set ptris))
  for tri in extras list (
      V := normalToricVariety(rays Q, tri);
      isWellDefined V and isProjective V
      )

  get "!chiro2circuits <topcomfoo.in"
  circuits cQ
  circs = flipCandidates(transpose matrix rays Q, tri)

  bistellarFlip(tri, circs_3)
  T = triangulation(transpose matrix rays Q, tri)
  flipCandidates(transpose matrix rays Q, tri)
  debugLevel = 1
  newtris = generateTriangulations(T, Limit => 200, Homogenize => false)
  newtris = generateTriangulations(transpose matrix rays Q, tri, Limit => 200, Homogenize => false)
  T = triangulation(transpose matrix rays Q, vtris_0)
  newtris = generateTriangulations(T, Homogenize => false)
  set(newtris/max) - set ptris
  set(newtris/max) - (set vtris + set ptris)
  set(newtris/max) - set ctris
  -- why do the others not show up?
  o28 - vtris
  
  for c in circs list (
      newtri := bistellarFlip(tri, c);
      if newtri === null then continue;
      {c, newtri}
      )
  #ptris
  vtris
  ctris
  -- end of testing with vector triangulations, etc
  ----------------------------------------------------------

  
  vtris = findAllSimplicialFans(transpose matrix rays Q);
  #vtris == 121
  ptris = sort for t in vtris list if isTriangulationOfPolytope(Q,t) then max t else continue;
  assert(tris1 === ptris)

  alltris = findAllSimplicialFans(transpose matrix rays Q, Fine => false);
  #alltris == 369
  allptris = sort for t in alltris list if isTriangulationOfPolytope(Q,t) then max t else continue; -- 30
  
  -- findAllFRSTs works as follows currently:
  -- 1. find all fine regular triangulations of the point set, including origin
  -- 2. take those that are star.  In this case, these are the ones that
  --    all subsets contain the last index.
  -- 3. remove the origin from each simplex.

  -- Another way: compute all vector triangulations
  tris = topcomAllTriangulations(transpose matrix rays Q, Homogenize => false, Fine => true); -- 121 of these
  #tris == 121
  
  topcomAllTriangulations(A1, Fine => true)

  elapsedTime findAllFRSTs Q;
  options  allTriangulations
  options topcomAllTriangulations
  A = transpose matrix rays Q
  topcomAllTriangulations A 
  topcomAllTriangulations(A, Fine => true)
  topcomAllTriangulations(A, Homogenize => false, Fine => true)

  dump Q
  
///


-*
  restart
  needsPackage "StringTorics"
*-
TEST ///
  --topes = kreuzerSkarke(10, Limit => 1000);
  -- this is topes_57 (which has interior facet lattice points).
  tope = KSEntry "4 11  M:29 11 N:17 10 H:10,22 [-24] id:57
     1   0   0   0   0   3   3  -1  -2  -3  -4
     0   1   0   0   1  -2  -3   2   2   0   1
     0   0   1   0   0  -1  -1  -1  -1   0   0
     0   0   0   1   1  -1  -1  -1  -1   1   1
  "  

  Q = reflexivePolytope tope

  -- we need to exert some care that
  -- vertices, latticePoints, rays, annotatedFaces all refer to the same indices!
  elapsedTime assert(
      latticePoints Q
      ===
      {{-1, -1, -1, 0}, {-1, -1, -1, 1}, {-1, -1, 0, 0}, {-1, 0, -1, -1},
          {0, -1, -1, 0}, {0, 0, -1, -1}, {0, 0, -1, 2}, {0, 0, 2, -1},
          {0, 1, -1, -1}, {1, 1, -1, 2}, {0, 0, 0, -1}, {0, 0, 0, 1},
          {0, 0, 1, -1}, {0, 0, 1, 0}, {0, 0, -1, 0}, {0, 0, -1, 1},
          {0, 0, 0, 0}}
      )
  elapsedTime assert(
      rays Q
      ===
      (latticePoints Q)_{0..13}
      )
  elapsedTime assert(
      vertices Q
      ===
      (latticePoints Q)_{0..9}
      )
  assert(faceDimensions Q
      ===
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 3, 3, 4}
      )
  assert(faceDimension(Q, 10) == 1)

  -- the following are cached and computed together
  -- degrees, basisIndices, toricBasisIndices
  elapsedTime degrees Q -- currently a bit expensive
  elapsedTime transpose matrix degrees Q
  assert(
      basisIndices Q
      ===
      {0, 1, 2, 3, 4, 5, 6, 7, 10, 11}
      )

  -- the following are cached and computed together
  elapsedTime assert(hh^(1,1) Q == 10)
  assert(hh^(1,2) Q == 22)
  assert isFavorable Q

  -- the following are cached and computed together
  -- automorphisms, autPermutations
  -- TODO: add these back in
  elapsedTime automorphisms Q -- in this case, just the identity
  debug StringTorics
  elapsedTime automorphismsAsPermutations Q

  -- point and vector configurations
  -- TODO.
  -- elapsedTime findAllFRSTs Q; -- too many for this h11=9 example??
  -- options  allTriangulations
  -- options topcomAllTriangulations
  -- A = transpose matrix rays Q
  -- topcomAllTriangulations A -- 105
  -- topcomAllTriangulations(A, Fine => true) -- 75
  -- topcomAllTriangulations(A, Homogenize => false, Fine => true) -- 121
  -- dump Q
///

end--
-- We want ReflexivePolytope's to be easy to construct in dimensions up to 6, 7, 8...
Q = reflexive KSEntry "4 8  M:10 8 N:55 15 H:50,5 [90]
    1    0    0   -1    1   -1    0    2
    0    1    0    0   -1    0    0   -1
    0    0    1   -1    1    1    0    0
    0    0    0    0    0    0    1   -1"
degrees Q

Q = reflexive KSEntry"4 6  M:11 6 N:119 7 H:100,4 [192]
    1    0    0   -3    0   -4
    0    1    0   -2    0    0
    0    0    1    0    0   -2
    0    0    0    0    1   -1"
elapsedTime degrees Q -- 133 sec.  90 x 94...  IMPROVE!

-- TODO: reflexive should take a WeightSystem...

methods weightSystems
ws5 = weightSystems(5, Degrees => 40, Reflexive => true);
Q = reflexive palpNVertices {40, 5, 6, 7, 7, 7, 8}
dim Q
vertices Q
latticePoints Q
faceDimensions Q
degrees Q
basisIndices Q

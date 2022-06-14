TEST ///
-*
  restart
  needsPackage "StringTorics"
*-
  str4 =         "         1   0   0   0   1   1   0  -1  -1  -2  -4
                  0   1   1   0  -2   2   3  -1  -4   1  -1
                  0   0   2   0  -2   4   4  -1  -4  -2  -4
                  0   0   0   1   0  -2  -2   2   2   0   2  
                  "
  A = matrixFromString str4
  P = convexHull A
  V = reflexiveToSimplicialToricVariety P
  X = completeIntersection(V, {-toricDivisor V})
  pt = base(a,b,c)
  Xa = abstractVariety(X, pt)
  I = intersectionRing Xa
  a*t_1
  hodgeDiamond X
  assert(h11OfCY P == 6)
  assert(h21OfCY P == 50)
  a*t_1 -- should still work...
///


TEST ///
-*
  restart
  needsPackage "StringTorics"
*-
  -- augment
    mat = "   1   0   0   1  -3   3   3  -3   5
            0   1   0   0   2  -4  -6   4  -8  
            0   0   1   0   2  -2  -2   0  -4 
            0   0   0   2   0  -4  -6   6  -6  "
  A = matrixFromString mat;
  assert(
      augment A 
      == 
      matrix {
          {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
          {1, 0, 0, 1, -3, 3, 3, -3, 5, 0}, 
          {0, 1, 0, 0, 2, -4, -6, 4, -8, 0}, 
          {0, 0, 1, 0, 2, -2, -2, 0, -4, 0},
          {0, 0, 0, 2, 0, -4, -6, 6, -6, 0}
          }
      )
///


TEST ///
  -- readSageTriangulation
  Ts = readSageTriangulations sageTri
  
  -- let's also switch to a different choice of rays
  -- Bstr is the lattice points (a.k.a. rays) used by Cody's code in Sage.
    Bstr = "[ 1 -1 -1  1 -1 -1 -1  1  1  0]
[ 0  1  0 -1  1  0  0 -1 -1  0]
[-1  0  0  0  1  1  0  0  0  0]
[ 2  0  0  1 -1 -1 -1 -1  0  0]"
  -- rays coming from M2
  Amat = matrix {{-1, -1, -1, -1, -1, 1, 1, 1, 1}, {0, 0, 0, 1, 1, -1, -1, -1, 0}, {0, 0, 1, 0, 1, 0, 0, 0, -1}, {-1, 0, -1, 0, -1, -1, 0, 1, 2}}
  Bmat = matrixFromString Bstr
  
  (fromM2, toM2) = matchNonZero(Amat, Bmat)
  applyPermutation(toM2, Ts)
  
  assert(
      applyPermutation(toM2, {0,1,2,3}) 
      == 
      {1,3,7,8}
      )
///

-- commenting this out.  It should really be a test in ReflexivePolytopesDB
///
  L1 = kreuzerSkarke(5,57, Access => "wget");
  assert(#L1 == 197)

  L2 = kreuzerSkarke(5,57, Access => "curl");
  assert(#L2 == 197)
  
  assert(L1 === L2)
///

TEST ///
  -- XXX  
  -- Test functionality of triangulations, part 1. Basic tests
  -- 
-*  
restart
needsPackage "StringTorics"
*-
  -- WARNING: currently, we use the polar dual of a convex polytope to determine
  -- minimal faces, etc.  However, for this to work, the convex polytope
  -- MUST contain the origin in the interior.
  -- Here are the functions that won't work if this is not the case:
  --polar P -- doesn't work as expected, since the origin is not an interior point of the polytope.
  --  minimalFace(P, pts)
  --  dualFace(P,f)
  --  latticePointList(P,f)
  --  genus(P,f)
  --  annotatedFaces P
  --  annotatedFaces(i,P)
  -- TODO: it would be possible to get all of these to work, except genus, and dualFace,
  --  by implementing them a bit differently.  Is it worth it? Most polytopes of interest will contain
  -- the origin in the interior?

  A = transpose matrix{{0,0,0},{0,0,1},{0,1,0},{0,1,1},{1,0,0},{1,0,1},{1,1,0},{1,1,1}}
  A1 = matrix{{8:1}} || A
  P = convexHull A
  assert(dim P == 3)

  vertices P
  latticePoints P
  faces P
  faces(1,P)

  LP = latticePointList P
  assert(set LP === set {{1, 1, 1}, {0, 1, 1}, {1, 0, 1}, {1, 1, 0}, {0, 0, 1}, {0, 1, 0}, {1, 0, 0}, {0, 0, 0}})
  assert(set vertexList P === set {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}})
  assert(entries transpose vertexMatrix P == vertexList P)

  assert(# faceList P == 27)
  faceDimensionHash P
  for f in faceList P do assert (2^(dim(P,f)) === #f)
  assert(faceList(0,P) == for i from 0 to 7 list {i})

  Amat = transpose matrix LP
  tri = regularFineTriangulation Amat
  naiveIsTriangulation(Amat, tri)
  topcomIsTriangulation(Amat, tri)
  
  -- check what happens if Amat is homoogenized:
  AmatH = Amat || matrix{{8:1}}
  naiveIsTriangulation(AmatH, tri) -- this should be false...?
  assert not topcomIsTriangulation(AmatH, tri) -- good! it complains that the index sets are not full dimensional (I think that is good?)
  
  affineCircuits(Amat, tri)
  for x in affineCircuits(Amat, tri) list flip(tri, x)
  flips(Amat, tri)
  
  generateTriangulations(Amat, tri)
  allTriangulations(Amat, RegularOnly => false, ConnectedToRegular => false, Fine => false)
  
  generateTriangulations(Amat, tri, Regular => true)
  allTriangulations Amat
  
  -- let's check that this is a triangulation.
  -- part of what we are checking: calls relative to homogenization are correct, and types make sense.
  -- part 1: for each oriented circuit
  circs = orientedCircuits transpose matrix LP
  assert(circs == {
          {{0, 4}, {1, 2}}, {{0, 4, 5}, {1, 6}}, {{0, 4, 6}, {2, 5}}, {{0, 5}, {1, 3}}, 
          {{0, 5, 6}, {3, 4}}, {{0, 6}, {2, 3}}, {{0, 7}, {1, 2, 3}}, {{0, 7}, {1, 6}}, 
          {{0, 7}, {2, 5}}, {{0, 7}, {3, 4}}, {{0, 7}, {4, 5, 6}}, {{1, 2, 7}, {3, 4}}, 
          {{1, 3, 7}, {2, 5}}, {{1, 6}, {2, 3, 7}}, {{1, 6}, {2, 5}}, {{1, 6}, {3, 4}}, 
          {{1, 7}, {4, 5}}, {{2, 5}, {3, 4}}, {{2, 7}, {4, 6}}, {{3, 7}, {5, 6}}
          }
      )
  -- it is possible that another triangulation would be output.
  assert(tri == {{0, 1, 2, 3}, {1, 2, 3, 4}, {1, 3, 4, 5}, {2, 3, 4, 6}, {3, 4, 5, 6}, {4, 5, 6, 7}})
  for c in circs list (
      n1 := # select(tri, t -> isSubset(c#0, t));
      n2 := # select(tri, t -> isSubset(c#1, t));
      ok := n1 == 0 or n2 == 0 or n1 == #tri or n2 == #tri;
      (n1,n2,ok)
      )

  walls = tri/(x -> subsets(x, #x-1))//flatten
  nfacets = tally walls
  facs = (faces(1,P))/first
  walls = partition(k -> nfacets#k, keys nfacets)
  facs = for f in facs list latticePointList(P, f)
  for w in walls#1 list (
      # select(facs, f -> isSubset(w, f))
      )
  for w in walls#2 list (
      # select(facs, f -> isSubset(w, f))
      )
  
  walls#2 -- 6 walls here.  Compute the vector for each.
  matrix {for w in walls#2 list (
      --w = {2,3,4} -- a wall
      circ := select(tri, t -> isSubset(w, t));
      others := circ/(c -> toList(set c - set w));
      elems := (flatten others) | w;
      print elems;
      id_(ZZ^8)_elems * syz AmatH_elems
      )}
      --id_(ZZ^8)_{1,6,2,3,4} * syz AmatH_{1,6,2,3,4} 
  
///

TEST /// 
  -- XXX  
  -- simple test of the polyhedral functions here: on the cube with the origin as its only 
  -- interior point
  needsPackage "StringTorics"
  P = hypercube 3
  assert isCompact P -- These functions fail if P is not a polytope.
  assert(
    vertices P  == matrix(QQ, {
        {-1, 1, -1, 1, -1, 1, -1, 1}, 
        {-1, -1, 1, 1, -1, -1, 1, 1}, 
        {-1, -1, -1, -1, 1, 1, 1, 1}
        })
    )
  vertices2 = vertexList P
  vertices3 = entries transpose vertexMatrix P
  assert(vertices2 == vertices3)
  
  LP = latticePointList P
  LP2 = transpose entries lift(matrix {latticePoints P}, ZZ)
  assert(set LP === set LP2)
  assert(27 == # LP)
  assert(27 == # LP2)

  assert(# faceList P == 27)
  faceDimensionHash P
  for f in faceList P do assert (2^(dim(P,f)) === #f)
  assert(faceList(0,P) == for i from 0 to 7 list {i})

  hashTable for pt in LP list pt => minimalFace(P, pt)      
  for pt in LP do (
      f := minimalFace(P, pt);
      assert(2^(# select(pt, v -> v == 0)) == #f)
      )
  P2 = polar P
  for f in faceList P list (f, dualFace(P, f))
  faces P
  faceList P
  assert(set((flatten values faces P)/first) === set faceList P)
  for f in faceList P list (
      latticePointList(P, f)
      )
  
  -- annotatedFaces
  netList annotatedFaces P
  for p in annotatedFaces P do (
      f := p#1;
      assert(p#0 == dim(P,f));
      lps := latticePointList(P, f);
      assert(p#2 == lps);
      assert(p#3 == # interiorLatticePointList(P,f));
      assert(p#4 == genus(P,f))
      )  
  
  regularFineTriangulation vertices P
  regularStarTriangulation P -- it is a shame that this returns something different from regularFineTriangulation.
///

TEST ///
  -- Test of triangulation code
-*
  restart
  needsPackage "StringTorics"
*-  
  -- XXX
  topes = kreuzerSkarke 3;
  A = matrix topes_30
  A
  P = reflexivePolytope A
  Amat = transpose matrix latticePointList polytope P
  regularSubdivision(Amat, matrix{{0,0,1,3,6,9,20,30}}) -- seems incorrect.
  TRI = regularFineTriangulation Amat -- is this including the origin automatically?
  wts = regularTriangulationWeights(Amat, TRI)
  -- check that this is a triangulation!
  TRI2 = regularSubdivision(Amat, matrix{{2, 4, 2, 0, 0, 0, 0, 0}}) -- good!
  assert(sortTriangulation TRI === sortTriangulation TRI2) -- works!

  -- let's check 'affineCircuits'
  C = affineCircuits(Amat, TRI)  
  4! * volumeVector(Amat, TRI)
  flip(TRI, C_0) === null
  flip(TRI, C_1) === null
  TRI2 = flip(TRI, C_2)
  wts2 = regularTriangulationWeights(Amat, TRI2)
  wts2 == {-2, 6, 4, 0, 0, 0, 0, 0}    
  TRI2' = regularSubdivision(Amat, matrix{{-2, 6, 4, 0, 0, 0, 0, 0}}) -- good!
  assert(TRI2' == TRI2)
  flip(TRI, C_3) === null
  4! * volumeVector(Amat, flip(TRI, C_4))
  flip(TRI, C_5) === null
  flip(TRI, C_6)
  4! * volumeVector(Amat, flip(TRI, C_6))
  
  findAllFRSTs P
///

TEST ///
  -- Test functionality of triangulations, part 2. reflexive polytope in 4D
  -- We try one with h11=3
  needsPackage "StringTorics"
  tope = "4 13  M:64 13 N:8 7 H:3,51 [-96]
   1   0   0   2  -2   0  -1  -1   0   2   2  -3  -3
   0   1   1   3  -5   2   0   0   2   4   4  -6  -6
   0   0   4   0  -4   5   4  -1   0   1   0  -4  -5
   0   0   0   4  -4   1  -1  -1   1   5   5  -5  -5"
  A = matrix first kreuzerSkarke tope
  P = convexHull A
  P2 = polar P
  V = reflexiveToSimplicialToricVariety P
  assert isFavorable P
  assert isCompact P
  assert isCompact P2
  -- regularFineTriangulation (from topcom)
  -- regularFineStarTriangulation
  -- isRegularTriangulation (from topcom)
  nonzeroLP = transpose matrix drop(latticePointList P2,-1)
  LP = transpose matrix latticePointList P2
  chirotope nonzeroLP

  regularFineTriangulation nonzeroLP
  regularFineTriangulation LP
  regularStarTriangulation P2
///

TEST ///
  -- this is an example of a polytope with 200 lattice points.
  tope = "4 12  M:72 12 N:276 12 H:200,50 [300]
          1    0    0    0  -14   -1  -17   -4   -5   -9  -15  -29
          0    1    0    0   -9   -1  -11   -2   -4   -6  -10  -20
          0    0    1    0   -3    0   -4   -2   -2   -4   -6  -10
          0    0    0    1   -1    1   -1    1    2    2    2    2"
  A = matrix first kreuzerSkarke tope
  P = convexHull A
  P2 = polar P
  elapsedTime faces P2;
  elapsedTime halfspaces P2
  LP = latticePoints P2
  # faces(1,P2)
  elapsedTime regularFineTriangulation matrix transpose latticePointList P2;
  V = reflexiveToSimplicialToricVariety P
  rays V
  max V
///

TEST ///
  -- testing triangulations of fans
  -- Our plan: start with example #26 from Kreuzer-Skarke with h11=5, h12=57
  --  this one has a number of triangulations.
  -- How do we create Amat?  This is the way:
  --      4 11  M:58 11 N:10 8 H:5,51 [-92]
  -- XXXXXXXXXX This test is failing May 2022.
-*
  restart
*-
  needsPackage "StringTorics"
  mat = "  1   1   1  -1   0   1   1  -1  -3  -1  -3
         0   2   0   0   0   0   2  -2  -2  -2  -4
         0   0   2  -2   0  -2   2   2  -2   4   2
         0   0   0   0   1  -2   0   2   0   2   2"
  A = matrixFromString mat
  P1 = convexHull A -- (extremal) vertices in RR^4
  P2 = polar P1
  LP = latticePoints P2 -- these will be the origin, the extremal vertices of P2 and possibly some more.
  Amat = matrix {select(LP, x -> x != 0)}
  elapsedTime   allTRIS = generateTriangulations Amat; -- removed in commit 33e77a592d2890c7ebf134e13b95d5915a624039

  -- now let's change these to sage indexing
    Bstr = "[ 1 -1 -1  1 -1 -1 -1  1  1  0]
            [ 0  1  0 -1  1  0  0 -1 -1  0]
            [-1  0  0  0  1  1  0  0  0  0]
            [ 2  0  0  1 -1 -1 -1 -1  0  0]"
  -- rays coming from M2
  Bmat = matrixFromString Bstr
  
  (fromM2, toM2) = matchNonZero(Amat, Bmat)

  Ts = readSageTriangulations sageTri
  --elapsedTime for T in Ts do time checkFan(Bmat, T) -- this takes a while (24 seconds), too long for testing
  elapsedTime checkFan(Bmat, Ts_5)
  applyPermutation(fromM2, allTRIS)

  -- the following are all in this list
  time TRI = regularFineStarTriangulation Amat -- this appears to not be returning regular triangulations?

  -- make a toric variety from one of the triangulations:
  elapsedTime assert({(true, true, true)} === 
      unique for T in allTRIS list (
      X = normalToricVariety(entries transpose Amat, T);
      time (isSimplicial X, isComplete X, isProjective X)
      )
  )
///

TEST ///
  -- analyze polytopes with h11=30, h12=50.
  -- grabbed 3000 examples, but there are more!
-*
  restart
  needsPackage "StringTorics"
*-
  -- one of the "favorable" examples
  --  4 7  M:51 7 N:45 7 H:30,50 [-40]
  mat = "    1    0    0    0    0   -4  -10
    0    1    0    0    2    2  -10
    0    0    1    0   -3   -7    5
    0    0    0    1    2    2   -4"

  A = matrixFromString mat
  P1 = convexHull A -- (extremal) vertices in RR^4
  P2 = polar P1
  LP = latticePoints P2 -- these will be the origin, the extremal vertices of P2 and possibly some more.
  Amat = matrix {select(LP, x -> x != 0)}

  elapsedTime regularStarTriangulation P2;
  elapsedTime   TRI = regularFineStarTriangulation Amat;
  assert(sort unique flatten TRI == toList(0..numcols Amat - 1))
///

TEST ///
      -- 4 10  M:105 10 N:71 10 H:50,80 [-60]
     mat = "     1    0    0    0   -1   -1   -1   -5  -15  -15
          0    1    0    0    0   -2   -2   -6   -8  -16
          0    0    1    0    1    1    0   -2   -6   -6
          0    0    0    1   -1    1    2    4    0    8"

  A = matrixFromString mat
  P1 = convexHull A -- (extremal) vertices in RR^4
  P2 = polar P1
  LP = latticePoints P2 -- these will be the origin, the extremal vertices of P2 and possibly some more.
  Amat = matrix {select(LP, x -> x != 0)}

  elapsedTime   TRI = regularFineStarTriangulation Amat;
  assert(sort unique flatten TRI == toList(0..numcols Amat - 1))
///

TEST ///
  -- creating simplicial toric varieties from data base
      -- 4 10  M:105 10 N:71 10 H:50,80 [-60]
     mat = "     1    0    0    0   -1   -1   -1   -5  -15  -15
          0    1    0    0    0   -2   -2   -6   -8  -16
          0    0    1    0    1    1    0   -2   -6   -6
          0    0    0    1   -1    1    2    4    0    8"

  A = matrixFromString mat
  P1 = convexHull A
  P2 = polar P1
  (LP,tri) = regularStarTriangulation P2
  V = normalToricVariety(LP,tri)  
  assert isSimplicial V
  assert not isSmooth V
  --assert elapsedTime isWellDefined V -- this currently takes some time

  (LP,tri) = regularStarTriangulation(2,P2)
  V = normalToricVariety(LP,tri)  
  assert isSimplicial V
  assert not isSmooth V
  --assert isWellDefined V -- this currently takes some time
  
  V1 = reflexiveToSimplicialToricVariety P1
  assert isSimplicial V
  assert not isSmooth V
///

TEST ///
  -- toric complete intersection cohomology code
  -- YYY
-*
  restart
  needsPackage "StringTorics"
*-
  mat = "  1   1   1  -1   0   1   1  -1  -3  -1  -3
         0   2   0   0   0   0   2  -2  -2  -2  -4
         0   0   2  -2   0  -2   2   2  -2   4   2
         0   0   0   0   1  -2   0   2   0   2   2"
  A = matrixFromString mat
  P = convexHull A
  V = reflexiveToSimplicialToricVariety P
  KV = toricDivisor V
  X = completeIntersection(V, {-KV})
  hodgeDiamond X
  assert(h11OfCY P == 5)
  assert(h21OfCY P == 51)
  assert(cohomologyVector X == {1,0,0,1})
  cohomologyVector(X, degree(V_0+V_1+V_2))
  allsums = drop(subsets for i from 0 to #rays V - 1 list V_i, 1);
  allsums = allsums/sum;
  elapsedTime for i from 0 to 100 list cohomologyVector(X, degree allsums_i)
  for i from 0 to 10 list cohomologyVector(X, 2 * degree allsums_i)
///

-- example: regular star triangulations
///
-*
  restart
*-
  needsPackage "StringTorics"
  mat = "  1   1   1  -1   0   1   1  -1  -3  -1  -3
         0   2   0   0   0   0   2  -2  -2  -2  -4
         0   0   2  -2   0  -2   2   2  -2   4   2
         0   0   0   0   1  -2   0   2   0   2   2"
  A = matrixFromString mat
  P = convexHull A
  V = reflexiveToSimplicialToricVariety P
  pts = transpose matrix rays V
  pts = pts || matrix{{numColumns pts : 1}}
  T = max V -- triangulation
  -- is this list correct, or do we need to "homogenize 'pts'?
  annotatedFaces polar P
  ac = select(affineCircuits(pts,T), x -> #x#0 > 1 and #x#1 > 1)
  ac = unique(ac/sort//sort)
  netList oo
  volumeVector(pts, T)
  ac#0
  T1 = flip(T,ac#0)
  volumeVector(pts, T1)
  checkFan(pts, T1)
  checkFan(pts, T)
  isRegularTriangulation(pts,T)
  isRegularTriangulation(pts,T1)
  triS = new MutableHashTable from {T=>true}
  ac = select(affineCircuits(pts,T), x -> #x#0 > 1 and #x#1 > 1)
  newT = for a in ac list (t := flip(T,a); if t === null then continue else t)
  Ts = join({T},newT)
  Ts/(t -> volumeVector(pts,t))/sum
  Ts/(t -> elapsedTime checkFan(pts,t))
  Ts/(t -> sort unique flatten t)
  Ts/(t -> isRegularTriangulation(pts,t))
  unique oo
///

TEST ///
-*
  restart
*-
  -- from Kreuzer-Skarke database
  -- 4 9  M:32 9 N:11 8 H:6,30 [-48]
  polystr = "   1   0   1   1  -1   1   0  -1  -2
    0   1   0   0   0  -2  -2   2   2
    0   0   2   0  -2  -2  -2   2   2
    0   0   0   2  -2   2   1   0  -1"

  A = matrixFromString polystr
  P1 = convexHull A
  P2 = polar P1

  elapsedTime LP1 = latticePointList P1
  V1 = vertexList P1
  assert(take(LP1,#V1) == V1)

  elapsedTime LP2 = latticePointList P2
  V2 = vertexList P2
  assert(take(LP2,#V2) == V2)

  elapsedTime assert(faceList(0,P1) == for i from 0 to 8 list {i})
  assert(# faceList(1,P1) == 22)
  assert(# faceList(2,P1) == 21)
  assert(# faceList(3,P1) == 8)

  elapsedTime assert(faceList(0,P2) == for i from 0 to 7 list {i})
  assert(# faceList(1,P2) == 21)
  assert(# faceList(2,P2) == 22)
  assert(# faceList(3,P2) == 9)

-*
  faceList(4,P1) -- what should this do?
  assert(faceList(-1,P1) == {{}}) -- not correct yet
  faceList(4,P2) -- what should this do?
  assert(faceList(-1,P2) == {{}}) -- not correct yet
*-

  -- now test dual faces...    
  dualfaces = for f in faceList(1,P2) list dualFace(P2,f)
  origfaces = for g in dualfaces list dualFace(P1,g)
  assert(origfaces == faceList(1,P2))
  for g in dualfaces do assert(2 == dim(P1,g))
  
  -- now test lattice point containment
  -- for each face, want the lattice points on that face
  faceList(1,P2)  
  for f in faceList(1,P2) list f => latticePointList(P2,f)
  for f in faceList(2,P2) list f => latticePointList(P2,f)
  -- now test interior lattice points
  -- need to run through all lattice points, and find max face containing it
  -- check this against Polyhedra code
  H = latticePointHash P2

  (vertexMatrix P2)_{0}
  Q = convexHull oo
  vertexList Q
  vertexMatrix Q
  
  latticePoints Q
  -- TODO: fix the following bug.  This results when using these functions in cases when 
  -- the polytope isn't e.g. reflexive, or really, doesn't have the origin in the interior...
  -- latticePointList Q -- fails, since polar Q isn't really what should be used here...

  for f in (faceList P2)_{0} do (
      Q := convexHull (vertexMatrix P2)_f;
      lp := (latticePoints Q)/(m -> flatten entries m);
      lpi := latticePointList(P2,f);
      lpi2 := lp/(p -> H#p);
      assert(lpi == sort lpi2)
      )

  H = latticePointHash P1
  for f in faceList P1 do (
      Q := convexHull (vertexMatrix P1)_f;
      lp := (latticePoints Q)/(m -> flatten entries m);
      lpi := latticePointList(P1,f);
      lpi2 := lp/(p -> H#p);
      assert(lpi == sort lpi2)
      )
  
  -- test interior lattice point code
  lpi = (latticePointHash P1)#{-2, 2, 2, -1}
  assert(minimalFace(P1, {-2,2,2,-1}) == {0})
  
  assert(interiorLatticePointList(P1, {0,2,3,5,7}) == {24, 26, 28})

  hashTable for f in faceList(2,P2) list (
      f => {dim(P2,f), 
          latticePointList(P2,f), 
          # interiorLatticePointList(P2,f), 
          # interiorLatticePointList(P1, dualFace(P2,f))}
      )
  allinfo = sort for f in faceList P2 list (
      {dim(P2,f), 
          f, 
          latticePointList(P2,f), 
          # interiorLatticePointList(P2,f), 
          # interiorLatticePointList(P1, dualFace(P2,f))
          }
      )
  allinfo'ans = {
      {0, {0}, {0}, 1, 0}, 
      {0, {1}, {1}, 1, 0}, 
      {0, {2}, {2}, 1, 0}, 
      {0, {3}, {3}, 1, 0}, 
      {0, {4}, {4}, 1, 0}, 
      {0, {5}, {5}, 1, 0}, 
      {0, {6}, {6}, 1, 0}, 
      {0, {7}, {7}, 1, 0}, 
      {1, {0, 1}, {0, 1}, 0, 1}, 
      {1, {0, 2}, {0, 2}, 0, 0}, 
      {1, {0, 3}, {0, 3}, 0, 0}, 
      {1, {0, 5}, {0, 5}, 0, 0}, 
      {1, {0, 6}, {0, 6}, 0, 0}, 
      {1, {1, 2}, {1, 2}, 0, 0}, 
      {1, {1, 4}, {1, 4, 8}, 1, 0}, 
      {1, {1, 5}, {1, 5}, 0, 1}, 
      {1, {1, 6}, {1, 6}, 0, 0}, 
      {1, {1, 7}, {1, 7}, 0, 1}, 
      {1, {2, 3}, {2, 3}, 0, 0}, 
      {1, {2, 4}, {2, 4}, 0, 0}, 
      {1, {2, 5}, {2, 5}, 0, 1}, 
      {1, {3, 4}, {3, 4}, 0, 0}, 
      {1, {3, 6}, {3, 6}, 0, 0}, 
      {1, {4, 5}, {4, 5}, 0, 0}, 
      {1, {4, 6}, {4, 6}, 0, 0}, 
      {1, {4, 7}, {4, 7}, 0, 1}, 
      {1, {5, 6}, {5, 6, 9}, 1, 3}, 
      {1, {5, 7}, {5, 7}, 0, 1}, 
      {1, {6, 7}, {6, 7}, 0, 1}, 
      {2, {0, 1, 2}, {0, 1, 2}, 0, 0}, 
      {2, {0, 1, 3, 4}, {0, 1, 3, 4, 8}, 0, 1}, 
      {2, {0, 1, 5}, {0, 1, 5}, 0, 0}, 
      {2, {0, 1, 6}, {0, 1, 6}, 0, 1}, 
      {2, {0, 2, 3}, {0, 2, 3}, 0, 1}, 
      {2, {0, 2, 5}, {0, 2, 5}, 0, 0}, 
      {2, {0, 3, 6}, {0, 3, 6}, 0, 1}, 
      {2, {0, 5, 6}, {0, 5, 6, 9}, 0, 1}, 
      {2, {1, 2, 4}, {1, 2, 4, 8}, 0, 1}, 
      {2, {1, 2, 5}, {1, 2, 5}, 0, 0}, 
      {2, {1, 4, 7}, {1, 4, 7, 8}, 0, 1}, 
      {2, {1, 5, 6}, {1, 5, 6, 9}, 0, 0}, 
      {2, {1, 5, 7}, {1, 5, 7}, 0, 0}, 
      {2, {1, 6, 7}, {1, 6, 7}, 0, 0}, 
      {2, {2, 3, 4}, {2, 3, 4}, 0, 1}, 
      {2, {2, 3, 5, 6}, {2, 3, 5, 6, 9}, 0, 1}, 
      {2, {2, 4, 5}, {2, 4, 5}, 0, 1}, 
      {2, {3, 4, 6}, {3, 4, 6}, 0, 1}, 
      {2, {4, 5, 6}, {4, 5, 6, 9}, 0, 0}, 
      {2, {4, 5, 7}, {4, 5, 7}, 0, 0}, 
      {2, {4, 6, 7}, {4, 6, 7}, 0, 0}, 
      {2, {5, 6, 7}, {5, 6, 7, 9}, 0, 1}, 
      {3, {0, 1, 2, 3, 4}, {0, 1, 2, 3, 4, 8}, 0, 1}, 
      {3, {0, 1, 2, 5}, {0, 1, 2, 5}, 0, 1}, 
      {3, {0, 1, 3, 4, 6, 7}, {0, 1, 3, 4, 6, 7, 8}, 0, 1}, 
      {3, {0, 1, 5, 6}, {0, 1, 5, 6, 9}, 0, 1}, 
      {3, {0, 2, 3, 5, 6}, {0, 2, 3, 5, 6, 9}, 0, 1}, 
      {3, {1, 2, 4, 5, 7}, {1, 2, 4, 5, 7, 8}, 0, 1}, 
      {3, {1, 5, 6, 7}, {1, 5, 6, 7, 9}, 0, 1}, 
      {3, {2, 3, 4, 5, 6}, {2, 3, 4, 5, 6, 9}, 0, 1}, 
      {3, {4, 5, 6, 7}, {4, 5, 6, 7, 9}, 0, 1},
      {4, {0, 1, 2, 3, 4, 5, 6, 7}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 1, 0}
      }
  assert(allinfo == allinfo'ans)
  assert(allinfo == annotatedFaces P2)  

  hodgeOfCYToricDivisors P1
  hodgeOfCYToricDivisors P2
  elapsedTime hashTable for p in latticePointList P2 list p => hodgeCY(P1, p)
  elapsedTime hashTable for p in latticePointList P1 list p => hodgeCY(P2, p)
  
  netList annotatedFaces(0,P1)
  netList annotatedFaces(3,P2)
  netList annotatedFaces(1,P1)
  netList annotatedFaces(2,P1)
  netList annotatedFaces(3,P1)
///

TEST ///
  -- Test of h11 and h21 formulae for an example from Kreuzer-Skarke database.
  -- as well as minimalFace, dim, genus, isFavorable.
-*
  restart
*-
  -- from Kreuzer-Skarke database
  -- 4 9  M:32 9 N:11 8 H:6,30 [-48]
  polystr = "   1   0   1   1  -1   1   0  -1  -2
    0   1   0   0   0  -2  -2   2   2
    0   0   2   0  -2  -2  -2   2   2
    0   0   0   2  -2   2   1   0  -1"

  A = matrixFromString polystr
  P1 = convexHull A
  P2 = polar P1

  elapsedTime assert(h11OfCY P1 == 6)
  elapsedTime assert(h11OfCY P2 == 30)
  elapsedTime assert(h21OfCY P1 == 30)
  elapsedTime assert(h21OfCY P2 == 6)

  assert isFavorable P1  
  assert not isFavorable P2 
  
  for f in faceList P2 list {dim(P2,f), genus(P2,f)}
  minfaces = for lp in latticePointList P2 list dim(P2,minimalFace(P2,lp))
  assert(minfaces == {0,0,0,0,0,0,0,0,1,1,4})
///

TEST ///
  -- We work on one example in 4 dimensions, where we know the answers (or have computed them elsewhere).
  -- Second polytope (index 1) on h11=3 Kreuzer-Skarke list of 4d reflexive polytopes for h11=3.
  -- 4 5  M:48 5 N:8 5 H:3,45 [-84]
  str = "  1   0   2   4  -8
         0   1   5   3  -9
         0   0   6   0  -6
         0   0   0   6  -6
  "
  M = matrixFromString str
  assert(M == matrix {{1, 0, 2, 4, -8}, {0, 1, 5, 3, -9}, {0, 0, 6, 0, -6}, {0, 0, 0, 6, -6}})

  P = convexHull M  
  assert(# latticePointList P == 48)
  assert(# latticePointList polar P == 8)
  assert(h11OfCY P == 3)
  assert(h21OfCY P == 45)
  assert(h11OfCY polar P == 45)
  assert(h21OfCY polar P == 3)
  
  -- Now compute all of the cohomologies of the (irreducible) toric divisors 
  LP = latticePointList polar P
  assert(#LP == 8)
  LP = drop(LP, -1)
  cohoms = for v in LP list hodgeOfCYToricDivisor(P, v)
  cohomH = hodgeOfCYToricDivisors P
  cohoms1 = for v from 0 to #LP-1 list cohomH#v -- last LP is the origin, which isn't one of these divisors
  assert(cohoms == cohoms1)
  assert isFavorable P

  -- Now compute all of the cohomologies of the (irreducible) toric divisors for the polar dual
  LP = latticePointList P
  LP = drop(LP, -1)
  cohoms = for v in LP list hodgeOfCYToricDivisor(polar P, v)
  cohomH = hodgeOfCYToricDivisors polar P
  cohoms1 = for v from 0 to #LP-1 list cohomH#v
  assert(cohoms == cohoms1)
  assert not isFavorable polar P
///

TEST ///
  -- id=0 h11=11
  -- 4 13  M:23 13 N:16 13 H:11,18 [-14]
  str = "    1    0    0    0   -1    1    0    0    0    1   -1    1   -2
    0    1    0    0    1   -1    0    1    1   -1    1   -2    0
    0    0    1    0    1   -1    0    1    0    0   -1   -2    2
    0    0    0    1   -1    1   -1   -1   -1    1    1    0   -1
    "

  P = convexHull matrixFromString str -- last first eg11
  P2 = polar P
  assert(h11OfCY P == 11)
  assert(h21OfCY P == 18)
  assert(h11OfCY polar P == 18)
  assert(h21OfCY polar P == 11)

  -- Now compute all of the cohomologies of the (irreducible) toric divisors 
  LP = latticePointList polar P  
  LP = drop(LP, -1)
  assert(#LP == 15)
  
  --elapsedTime cohoms = for v in LP list hodgeOfCYToricDivisor(P, v)
  cohomH = hodgeOfCYToricDivisors P
  cohoms1 = for v from 0 to #LP-2 list cohomH#v

  -- the following line is too slow, and 
  cohoms = for v in drop(LP,-1) list hodgeOfCYToricDivisor(P, v)
  assert(cohoms == cohoms1)
  assert isFavorable P

  -- Now compute all of the cohomologies of the (irreducible) toric divisors for the polar dual
  LP = latticePointList P
  LP = drop(LP, -1)
  assert(#LP == 22)
  elapsedTime cohoms = for v in LP list hodgeOfCYToricDivisor(polar P, v)
  cohomH = hodgeOfCYToricDivisors polar P
  cohoms1 = for v from 0 to #LP-1 list cohomH#v
  assert(cohoms == cohoms1)
  -- the following line is too slow, and 
--  cohoms = for v in LP list hodgeOfCYToricDivisor(P, v) -- actually, gives an error here...
--  assert(cohoms == cohoms1)
  assert isFavorable polar P
///

TEST ///
-*
restart
*-
str = "    1    0    0    0  -11   -3   -1   -3
    0    1    0    0   -6   -2   -2   -6
    0    0    1    0   -2   -2   -2   -2
    0    0    0    1   -2    2    4    6 "
M = matrixFromString str
P = convexHull M
vertexMatrix P
lp = latticePointList polar P
sort for v in vertexList P list (latticePointHash P)#v

assert try (hodgeOfCYToricDivisor(P,{0,0,0,1}); false) else true -- {0,0,0,1} is not a lattice point
for p in lp list p => hodgeOfCYToricDivisor(P,p)
hodgeOfCYToricDivisors P

assert(h11OfCY(P) == 90)
assert(h21OfCY(P) == 50)
assert(h11OfCY(polar P) == 50)
assert(h21OfCY(polar P) == 90)
///

----------------------------
-- tests from MyPolyhedra --
----------------------------

TEST ///  
  A = transpose matrix {{-1,-1,2},{-1,0,1},{-1,1,1},{0,-1,2},{0,1,1},{1,-1,3},{1,0,-1},{1,1,-2}}
  tri = regularFineTriangulation A
  volumeVector(augment A,tri)
  P = convexHull A

  A = transpose matrix {{-1, 0, -1, -1}, {-1, 0, 0, -1}, {-1, 1, 2, -1}, {-1, 1, 2, 0}, {1, -1, -1, -1}, {1, -1, -1, 1}, {1, 0, -1, 2}, {1, 0, 1, 2}}
  C = transpose matrix latticePointList polar convexHull A
  tri = regularFineTriangulation C

  elapsedTime delaunaySubdivision C -- takes 20 seconds?!
  isRegularTriangulation(C, tri)
  P2 = polar convexHull A
  regularStarTriangulation P2
  volume P2
///


TEST ///
  A = transpose matrix {{1, 0, 0, 0}, {1, 2, 0, 0}, {1, 0, 2, 0}, {0, 0, 0, 1}, {0, 4, 0, 1}, {0, 0, 4, 1}, {-2, -4, -6, -3}, {-2, -6, -6, -3}, {-2, -6, -4, -3}}
  P = convexHull A
  P2 = polar P

  -- triangulations
  A1 = transpose matrix latticePointList P
  tri = regularFineTriangulation A1
  isRegularTriangulation(A1,tri)

  debugLevel = 3
  tri = regularStarTriangulation P
  normalToricVariety(drop(latticePointList P,-1), tri_1)
  -- elapsedTime isWellDefined oo -- ouch, pretty long (4/27/20) (96 sec)

  -- vertices
  vertexList P
  vertexMatrix P
  vertexList P2
  vertexMatrix P2

  -- latticePoints
  latticePointList P
  latticePointList P2
  
  -- interior lattice points
  debug StringTorics
  hash1 = hashTable for f in faceList P list f => interiorLatticePointList(P, f)
  hash2 = hashTable select(pairs hash1, (k,v) -> #v > 0)
  hash3 = P.cache.TCIInteriorLatticeHash
  assert(hash2 === hash3)

  -- faces
  faceDimensionHash P
  faceList P
  faceList(0, P)
  faceList(1, P)
  faceList(2, P)
  faceList(3, P)
  faceList(4, P)

  FL1 = for f in faceList P list dualFace(P,f)
  FL2 = for f in FL1 list dualFace(polar P,f)  
  FL3 = for f in FL2 list dualFace(P,f)  
  assert(FL1 == FL3)
  assert(faceList P == FL2)

  time for f in faceList P list dim(P,f)
  time for f in faceList P list latticePointList(P,f)
  time for lp in latticePointList P list minimalFace(P,lp)
  time for f in faceList P list genus(P,f)

  netList annotatedFaces P
  for i from 0 to 4 list netList annotatedFaces(i,P)
///


TEST /// -- test of functions here on the square and the cube
-*
  restart
  needsPackage "StringTorics"
*-
  square = transpose matrix{{1,1},{-1,1},{-1,-1},{1,-1}}
  
  regularFineTriangulation square
  assert(# allTriangulations square == 2)
    
  -- Now consider all of the lattice points of the square
  sq9 = transpose matrix latticePointList convexHull square
  assert(sq9 == matrix {{-1, -1, 1, 1, -1, 0, 0, 1, 0}, {-1, 1, -1, 1, 0, -1, 1, 0, 0}})

  -- test function from Topcom.
  t1 = regularFineTriangulation sq9
  regularTriangulationWeights(sq9, t1)
  fineStarTriangulation(sq9, t1) -- has central element removed.  Don't do that?
  delaunaySubdivision sq9 -- not a triangulation (4 squares).
  orientedCircuits sq9 -- many of these are not useful when considering only fine triangulations.

  -- Polyhedra command, bug fix currently in StringTorics.
  t2 = regularSubdivision(sq9, matrix{{4,4,4,4,1,1,1,1,0}})
  assert(
    sortTriangulation t2 
    == 
    {{0, 4, 5}, {1, 4, 6}, {2, 5, 7}, {3, 6, 7}, {4, 5, 8}, {4, 6, 8}, {5, 7, 8}, {6, 7, 8}}
    )

  t3 = regularSubdivision(sq9, matrix{{4,4,4,4,1,1,1,1,-4}})
  t4 = regularSubdivision(sq9, matrix{{1,1,1,1,1,1,1,1,-4}})
  
  affineCircuits(sq9, t3)
  elapsedTime affineCircuits(sq9, t1) -- TODO: change the name of this function? Why? affineOrientedCircuits?
  
  elapsedTime tris = allTriangulations sq9; -- computes all triangulations, Fine => false, regular.
  assert(387 == #tris) -- check this number. 
  elapsedTime tris1 = allTriangulations(sq9, Fine => true);
  assert(64 == #tris1)
  elapsedTime tris2 = generateTriangulations(sq9, t1);
  elapsedTime tris2a = generateTriangulations(sq9, t1, Regular => true); -- slow...
  #tris2 == #tris2a

  -- via topcom:
  elapsedTime assert(387 == # select(tris, t -> isRegularTriangulation(sq9, t))) -- slow
  -- via Polyhedra/
  elapsedTime assert(387 == # select(tris, t -> null =!= regularTriangulationWeights(sq9, t))) -- slow, same as Polyhedra, I think (same code is being used)
  -- 64 fine triangulations
  assert(64 == # select(tris, t -> isFine(sq9, t)))
  assert(16 == # select(tris, t -> isStar(sq9, t)))
  assert(1 == # select(tris, t -> isStar(sq9, t) and isFine(sq9, t)))

  -- generate only fine 
  elapsedTime finetris = allTriangulations(sq9, Fine => true);
  -- elapsedTime finetris1 = generateTriangulations(sq9, t1, Fine => true); -- TODO
  elapsedTime finetris1 = generateTriangulations(sq9, t1);
  assert(sort finetris == sort finetris1)

  finetris/volumeVector_(augment sq9) -- TODO: augment?  what to do with that?  
  --finetris/volumeVector_sq9 -- TODO: allow this?  In any case: this fails now...
  -- TODO: augment, in TopCom, change to a method.
  
  -- TODO: audit/allow "A" matrices to be homogeneous or not?
  --   one way: default is not, and Homogenize => true will add the extra row.
  --   use this for most routines taking point configurations or vector configurations.

  -- note: given two regular, fine, star, triangulations, they can be connected via 
  -- bistellar flips that keep the triangulation fine and star (and regular, I believe).
///

TEST ///
  -- triangulation code, test 1.  
  -- some triangulation code is in Topcom, Polyhedra
  -- (comes from KS database, h11=3, #232.
-*
  restart
*-  
  needsPackage "StringTorics"
  A0 = matrix {{1, 1, 1, 1, -11}, {0, 2, 2, 2, -10}, {0, 0, 4, 4, -8}, {0, 0, 0, 12, -12}}

  P2 = polar convexHull A0
  LP = latticePointList P2
  A  = transpose matrix LP
  
  -- We want triangulations of the point configuration given by the columns of A.
  T = regularFineTriangulation A
  assert isFine(A, T)
  assert isStar(A, T)
  assert isRegularTriangulation(A, T)
      
  -- second: is the union the same as the polytope? How do we tell?
  -- consider all facets of each cell.  Each must occur once or twice.  Internal ones
  -- must occur twice, ones on the boundary must occur once.
  facetsA = T/(t -> subsets(t, #t-1))//flatten//tally
  outerfacets = select(keys facetsA, k -> facetsA#k == 1)
  facetsP2 = (annotatedFaces(3, P2))/(x -> x#1) -- sets of LP's on each facet.
  all(outerfacets, f -> any(facetsP2, g -> isSubset(f, g))) -- this checks 
  -- now make sure no interior simplices overlap in their interior.
  volume P2
  sum for t in T list volume convexHull A_t

  elapsedTime tris = allTriangulations A;
  FRST = first select(tris, t -> isFine(A, t) and isStar(A,t) and isRegularTriangulation(A,t))
  volumeVector(augment A, FRST) -- WRONG?? Somehow an extra element is being considered...
  elapsedTime tris = generateTriangulations(A, T);

  wts = regularTriangulationWeights(A, T)
  assert(regularSubdivision(A, matrix{wts}) == T//sort)
  assert isRegularTriangulation(A, T) -- from Topcom.
  -- also check Topcom's code...
  flips(A, T) -- not unpacked... TODO: need to understand and document what this gives
  --A1 = A || splice matrix{{numColumns A: 1}}
  A1 = augment A
  --volumeVector(A, T) -- need to homogenize A? at bottom or at top?
  volumeVector(A1, T) -- need to homogenize A? at bottom or at top?
  -- check that T is a triangulation?
  circs = affineCircuits(A, T) -- NEEDS CLEANING UP (? why)
  flip(T, {{5}, {0, 1, 2}})  
  for c in circs list flip(T, c)

  for v in subsets(T, 2) list {v#0, v#1, sort unique flatten v}
  select(oo, v -> #v#2 == 6)

  -- now generate new triangulations from these.
  -- check: they are triangulations, are they regular?
  -- go from fine regular --> fine star regular
  -- can we get all affine circuits of A?  Yes: from topcom...  also, can just do it?
  
  -- flips, generating more triangulations.  Are they regular?
  -- volume of a polytope.
///

TEST /// -- medium size (h^11 = 15) example
-*
  restart
  needsPackage "StringTorics"
*-
  topes = kreuzerSkarke(15, Limit=>10, Access=>"wget")
  A1 = matrix topes_8
  P = convexHull A1
  P2 = polar P
  A = transpose matrix latticePointList P2  

  elapsedTime tri = regularStarTriangulation P2;
  A2 = transpose matrix first tri
  assert(A2 == submatrix(A, 0..numcols A-2))
  tri = tri_1/(x -> append(x, numcols A - 1))
  isFine(A, tri)
  isStar(A, tri)
  wts = regularTriangulationWeights(A, tri)
  elapsedTime regularSubdivision(A, matrix{wts}) -- this is slower than we would like
  assert(sortTriangulation oo == tri)

  circs = affineCircuits(A, tri)
  circs0 = select(circs, x -> not member(numcols A - 1, flatten x))  
  for c in circs0 list flip(tri, c)
  flip(tri, circs0_1)
  
  elapsedTime tris = generateTriangulations(A, tri, Limit => 50);
  tris/isFine_A//tally
  tris/isStar_A//tally
  elapsedTime(tris/regularTriangulationWeights_A);
  
  elapsedTime tri = regularFineTriangulation A; -- topcom, fast.
  --  elapsedTime tris = allTriangulations(A, Fine => true); -- pretty long, how many are there?

  ans = {{0, 2, 3, 5}, {0, 2, 3, 9}, {0, 2, 5, 15}, {0, 2, 6, 9}, {0, 2, 6, 13}, 
      {0, 2, 7, 13}, {0, 2, 7, 15}, {0, 3, 5, 11}, {0, 3, 9, 11}, {0, 5, 6, 11}, 
      {0, 5, 6, 13}, {0, 5, 13, 15}, {0, 6, 9, 11}, {0, 7, 13, 15}, {1, 2, 3, 9}, 
      {1, 2, 3, 10}, {1, 2, 4, 10}, {1, 2, 4, 12}, {1, 2, 6, 9}, {1, 2, 6, 12}, 
      {1, 3, 5, 10}, {1, 3, 5, 11}, {1, 3, 9, 11}, {1, 4, 5, 10}, {1, 4, 5, 12}, 
      {1, 5, 6, 11}, {1, 5, 6, 12}, {1, 6, 9, 11}, {2, 3, 5, 10}, {2, 4, 5, 10}, 
      {2, 4, 5, 16}, {2, 4, 12, 16}, {2, 5, 7, 15}, {2, 5, 7, 16}, {2, 6, 12, 13}, 
      {2, 7, 12, 13}, {2, 7, 12, 16}, {4, 5, 12, 14}, {4, 5, 14, 16}, {4, 8, 12, 14}, 
      {4, 8, 12, 16}, {4, 8, 14, 16}, {5, 6, 12, 14}, {5, 6, 13, 14}, {5, 7, 13, 14}, 
      {5, 7, 13, 15}, {5, 7, 14, 16}, {6, 8, 12, 13}, {6, 8, 12, 14}, {6, 8, 13, 14}, 
      {7, 8, 12, 13}, {7, 8, 12, 16}, {7, 8, 13, 14}, {7, 8, 14, 16}}
  assert(ans == regularFineStarTriangulation A)
///

TEST ///
-- This test is failing: May 2022.  It isn't a complete test anyway...
-- Test of intersection number computations.
-- This requires that V be favorable?
-*
  restart
  needsPackage "StringTorics"
*-
  debug StringTorics

  topes = kreuzerSkarke(3, Limit => 50);    
  A = matrix topes_30
  P = convexHull A
  (V, basisElems) = reflexiveToSimplicialToricVarietyCleanDegrees(P, CoefficientRing => ZZ/32003)
  basisElems -- for the moment, we ignore this, and write down all of the elements...
  GLSM = transpose matrix degrees ring V
  X = completeIntersection(V, {-toricDivisor V})
  Xa = abstractVariety(X, base())
  IX = intersectionRing Xa
  elemsToConsider = toList(0..numcols GLSM-1);
  triples = (subsets(elemsToConsider, 3))/sort//sort;
  Htriples = hashTable for a in triples list (
      (i,j,k) := toSequence a;
      val := integral(IX_i * IX_j * IX_k);
      if val == 0 then continue else {i,j,k} => val
    )


  -- Now we try the code above
  (singles, doubles, triples) = toSequence possibleNonZeros V

  Htriples = hashTable for a in join(singles, doubles, triples) list (
      (i,j,k) := toSequence a;
      val := integral(IX_i * IX_j * IX_k);
      if val == 0 then continue else {i,j,k} => val
    )

  assert(triples === (keys Htriples)/sort//sort) -- failing.
  
  elapsedTime tripleProductsCY V
  elapsedTime CY3NonzeroMultiplicities V -- much slower for small h11...
  
  ans = toSequence topologyOfCY3(V, basisElems)
  (h11, h21, H3, C, L) = ans
  hashTable for x in keys H3 list (
      if isSubset(x, basisElems) then x => H3#x else continue
      )


  (h11, h21, C, L) = toSequence topologyOfCY3(V, basisElems)  
  (h11', h21', C', L') = toSequence topologyOfCY3(V, {0,1,2}, Ring => ring C)  
  A = ring C;
  M = GLSM_{0,1,2}  
  phi1 = map(A, A, flatten entries((M) * transpose vars A))
  L' == phi1 L
  C' == phi1 C
  netList {L, L', phi1 L, phi1 L'}
  
///

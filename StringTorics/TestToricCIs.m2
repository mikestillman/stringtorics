-*
  restart
  needsPackage "StringTorics"
*-
TEST ///
  debug needsPackage "StringTorics" -- remove 'debug'
  -- Let's consider smooth toric surfaces.
  for i from 0 to 4 list (
      V = smoothFanoToricVariety(2, 1);
      picardGroup V
      )
  V = smoothFanoToricVariety(2, 2)
  transpose matrix degrees ring V
  X = completeIntersection(V, {-toricDivisor V}, Basis => {0,3}, Variables => {symbol a, symbol b})
  linearForm X
  intersectionForm X
  linearForm X
  hh^* OO_X(0,0) -- X is an elliptic curve
  hh^*(OO_X(-1,1))

  L = OO_X(-1,1)
  assert(hh^0(L) == 1)
  assert(hh^1(L) == 0)
  assert(hh^* L == {1,0})
  intersectionForm X
  chern_1 tangentBundle abstractVariety X
  c2Form X
///

-*
  restart
  needsPackage "StringTorics"
*-
TEST ///
  debug needsPackage "StringTorics" -- remove 'debug', used for linearForm.

  -*    -- code to generate this example
    topes = kreuzerSkarke(3, Limit => 50);
    A = matrix topes_30
    P = convexHull A
    (V, basisElems) = reflexiveToSimplicialToricVarietyCleanDegrees(P, CoefficientRing => ZZ/32003)
  *-

  verts = {{-1, -1, 0, 0}, {-1, -1, 0, 1}, {-1, -1, 2, 0}, {-1, 0, 0, 0}, {1, -1, -1, 1}, {1, 2, -1, -1}, {-1, -1, 1, 0}}
  maxcones = {{0, 1, 3, 4}, {0, 1, 3, 6}, {0, 1, 4, 6}, {0, 3, 4, 5}, {0, 3, 5, 6}, {0, 4, 5, 6}, {1, 2, 3, 5}, {1, 2, 3, 6}, {1, 2, 4, 5}, {1, 2, 4, 6}, {1, 3, 4, 5}, {2, 3, 5, 6}, {2, 4, 5, 6}}
  V = normalToricVariety(verts, maxcones, CoefficientRing => ZZ/32003)
  glsm = transpose matrix degrees ring V
  basiselems = {0, 5, 6}
  X = completeIntersection(V, {-toricDivisor V}, Basis => {0, 5, 6}, Variables => {symbol a, symbol b, symbol c})
  assert(dim X == 3)
  intersectionRing X
  h = linearForm X
  integral(h^3)
  assert(intersectionForm X == a^3-6*a^2*b+6*a*b^2+12*b^3-9*a^2*c+24*a*b*c+21*a*c^2-24*b*c^2-16*c^3)
  assert(c2Form X == 10*a + 60*b + 8*c)

  X = completeIntersection(V, {-toricDivisor V}, Basis => {0, 5, 6}, Variables => symbol a)

  D = completeIntersection(V, {-toricDivisor V, V_0}, Basis => {0, 5, 6}, Variables => {symbol a, symbol b, symbol c})
  dim D == 2
  saturate(ideal equations D, ideal V)
  intersectionRing D
  intersectionForm D
  c2Form D
  chern_2 tangentBundle abstractVariety D

  L = OO_X(1,1,2)
  assert(hh^* L == {4, 15, 0, 0})
  assert(variety L === X)
  assert(degree L == {1,1,2})

  S = ring V
  M = S^1/(ideal equations X) ** S^{{1,1,2}}
  rank HH^0(V, sheaf M) == 4
  --elapsedTime rank HH^1(V, sheaf M) -- takes a long time...
  --elapsedTime rank HH^2(V, sheaf M) == 0 -- 61 sec
  elapsedTime rank HH^3(V, sheaf M) == 0 -- quick
  -- can we get the euler characteristic of this?
  euler OO_V(1,1,2) -- -11
  euler OO_V(-3,-3,0) -- 0
  -- euler L (on X) is -11.  So we can deduce rank HH^1(V, sheaf M) == 15.
  
  -- TODO: ADD BACK IN once DanilovKhovanskii more functional
  -- needsPackage "DanilovKhovanskii"
  -- computeHodgeDeligne(-toricDivisor V)
  -- oo#1
  -- matrix for i from 0 to dim X list for j from 0 to dim X list (-1)^(i+j) * oo#(i,j)

  -- want to be able to turn this into a CalabiYauInToric...
  -- then we can check computations against each other too.

///

-*
  restart
  needsPackage "StringTorics"
*-
TEST ///
  debug needsPackage "StringTorics" -- remove 'debug'
  V = kleinschmidt(3, {2,1}, CoefficientRing => ZZ/101)
  rays V
  max V
  picardGroup V
  isSmooth V
  transpose matrix degrees ring V
  X = completeIntersection(V, {-toricDivisor V}, Basis => {3,0}, Variables => {a,b})
  dim X == 2
  hh^*(OO_X(0,0)) == {1, 0, 1}
  F = first equations X
  saturate(ideal F + ideal jacobian F, ideal V) -- X is smooth
  -- hh^*(OO_V(-1,3)) -- ouch!  needs to work...
///

-*
  restart
  needsPackage "StringTorics"
*-
"TEST"
///
  debug needsPackage "StringTorics" -- remove 'debug'?

  -- idea: create a toric simplicial variety of dim 6, corresponding to
  --ws = {12,1,1,2,2,2,2,2}
  ws = {20, 1, 1, 1, 2, 3, 6, 6}  
  -- now triangulate it using all of the lattice points not in facets
  palpNVertices ws
  palpNPoints ws
  Q = reflexive palpNVertices ws
  Ts = allTriangulations palpNPoints ws
  tri = first select(Ts, t -> isStar t and isFine t and isRegularTriangulation t)
  maxT = (max tri)/(d -> drop(d, -1)) -- this is the triangulation we want of the 6D point set, pts 0..7
  pts = drop(entries transpose palpNPoints ws, -1)
  V = normalToricVariety(pts, maxT)
  isWellDefined V
  isSimplicial V
  isSmooth V -- false...
  -- Now, given V, let's check nef partitions...
  --  {{0, 5}, {1, 2}, {3, 4, 6, 7}}
  Vs = for i from 0 to # rays V - 1 list V_i
  Vs/isNef
  -- this is *not* a nef partition!!
  isNef(V_0 + V_5)
  isNef(V_1 + V_2)
  isNef(V_3 + V_4 + V_5 + V_6 + V_7)
  isNef (-toricDivisor V)
  X = completeIntersection(V, {V_0 + V_5, V_1 + V_2, V_3 + V_4 + V_5 + V_6 + V_7}, Basis => {0,1,2}, Variables => {a,b,c})
  hodgeDiamond X -- this appears wrong!

  -- this one is!
  isNef(V_0 + V_4)
  isNef(V_3 + V_6)
  isNef(V_1 + V_2 + V_5 + V_7 + V_8)
  isNef (-toricDivisor V)
  X = completeIntersection(V, {V_0 + V_4, V_3 + V_6, V_1 + V_2 + V_5 + V_7 + V_8}, Basis => {0,1,2}, Variables => {a,b,c})
  hodgeDiamond X -- this appears correct! (And no variables...)

  D0 = completeIntersection(V, {V_0 + V_4, V_3 + V_6, V_1 + V_2 + V_5 + V_7 + V_8, V_0});
  hodgeDiamond D0
  D1 = completeIntersection(V, {V_0 + V_4, V_3 + V_6, V_1 + V_2 + V_5 + V_7 + V_8, V_1});
  hodgeDiamond D1
  D2 = completeIntersection(V, {V_0 + V_4, V_3 + V_6, V_1 + V_2 + V_5 + V_7 + V_8, V_2});
  hodgeDiamond D2
  D3 = completeIntersection(V, {V_0 + V_4, V_3 + V_6, V_1 + V_2 + V_5 + V_7 + V_8, V_3});
  hodgeDiamond D3
  D4 = completeIntersection(V, {V_0 + V_4, V_3 + V_6, V_1 + V_2 + V_5 + V_7 + V_8, V_4});
  hodgeDiamond D4
  D5 = completeIntersection(V, {V_0 + V_4, V_3 + V_6, V_1 + V_2 + V_5 + V_7 + V_8, V_5});
  hodgeDiamond D5
  D6 = completeIntersection(V, {V_0 + V_4, V_3 + V_6, V_1 + V_2 + V_5 + V_7 + V_8, V_6});
  hodgeDiamond D6
  D7 = completeIntersection(V, {V_0 + V_4, V_3 + V_6, V_1 + V_2 + V_5 + V_7 + V_8, V_7});
  hodgeDiamond D7  -- D7 has 2 components, probably both P^2's!
  D8 = completeIntersection(V, {V_0 + V_4, V_3 + V_6, V_1 + V_2 + V_5 + V_7 + V_8, V_8});
  hodgeDiamond D8  -- intersection here appears to be 0??

  intersectionRing X
  
  
  intersectionRing X
  faceDimensions Q
  annotatedFaces Q
  P = polytope Q
  makeCY
  methods triangulations
  verts = {{0, 0, 0, 0, 0, 1}, {-1, -1, -2, -3, -6, -6}, {0, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0}, {0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 1, 0}, {1, 0, 0, 0, 0, 0}, {0, 0, 0, -1, -2, -2}}
  vertsM = {{-1, -1, -1, -1, -1, -1}, {11, -1, -1, -1, -1, -1}, {-1, 5, -1, -1, -1, -1}, {-1, -1, 5, -1, -1, -1}, {-1, -1, -1, 5, -1, -1}, {-1, -1, -1, -1, 5, -1}, {-1, -1, -1, -1, -1, 5}}
  P = convexHull transpose matrix vertsM
  isReflexive P
  V0 = normalToricVariety P
  dim P == 6
  isSimplicial V0
  rays V0
  max V0
  latticePoints polar P
  vertices polar P
  Q = reflexive P
  makeCY Q
///  

debug needsPackage "ReflexivePolytopesDB"
kreuzerSkarkeURL ZZ := URL => opts -> (h11) -> (
    getopt := (ksStr, m2Key, opts) -> (
        if opts#m2Key === null then ""
        else (ksStr | "=" | toString opts#m2Key)
        );
    valV := getopt("V", Vertices, opts);
    valM := getopt("M", LatticePoints, opts);
    valN := getopt("N", DualLatticePoints, opts);
    valF := getopt("F", Facets, opts);
    valH21 := getopt("h12", H12, opts);
    valL := getopt("L", Limit, opts);
    str := "https://quark.itp.tuwien.ac.at/cgi-bin/cy/cydata.cgi?h11="|toString h11;
    ksopts := {valV, valM, valN, valF, valH21, valL};
    ksopts = select(ksopts, x -> x =!= "");
    URL concatenate between("&", prepend(str,ksopts))
    )

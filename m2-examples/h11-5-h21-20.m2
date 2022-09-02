-- Consider all examples with h11=5, h21=20.
-- There is only one polytope with this data.
-- 142 triangulations.
-- 2 possibly different cubic forms...
--   both of these seem to be smooth cubics (over QQ), and might be equivalent.

restart
needsPackage "StringTorics"

-- First, lets find the polytope (CYPolytopeData) we will use.
  F = openDatabase "polytopes-h11-5.dbm"
  elapsedTime Qs = for i from 0 to 4989 list cyPolytopeData F#(toString i); -- 3.25 seconds to read this...
  close F
  assert(1 == # select(Qs, Q -> h21OfCY Q == 20))
  Q = first select(Qs, Q -> h21OfCY Q == 20)
  dump Q
  -- CYPolytopeData
  --   rays:{{-1, -1, 0, 0}, {-1, -1, 0, 1}, {-1, 0, -1, 0}, {-1, 1, -1, 0}, {0, -1, 0, -1}, {0, -1, 1, -1}, {0, 0, -1, -1}, {1, 0, -1, -1}, {1, 1, 1, 1}}
  --   face dimensions:{0, 0, 0, 0, 0, 0, 0, 0, 0}
  --   id:0
  --   favorable:true
  --   h11:5
  --   h21:20
  --   basis indices:{0, 1, 2, 4, 5}
  --   glsm:{{1, 0, 0, 0, 0}, {0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}, {0, -1, -1, 1, 2}, {0, 0, 0, 1, 0}, {0, 0, 0, 0, 1}, {1, 5, 2, -2, -5}, {0, -2, -1, 1, 3}, {1, 2, 1, 0, -1}}
  --   annotated faces:{{0, {0}, {0}, 1, 0}, {0, {1}, {1}, 1, 0}, {0, {2}, {2}, 1, 0}, {0, {3}, {3}, 1, 0}, {0, {4}, {4}, 1, 0}, {0, {5}, {5}, 1, 0}, {0, {6}, {6}, 1, 0}, {0, {7}, {7}, 1, 0}, {0, {8}, {8}, 1, 0}, {1, {0, 1}, {0, 1}, 0, 0}, {1, {0, 2}, {0, 2}, 0, 0}, {1, {0, 3}, {0, 3}, 0, 0}, {1, {0, 4}, {0, 4}, 0, 0}, {1, {0, 5}, {0, 5}, 0, 0}, {1, {1, 2}, {1, 2}, 0, 0}, {1, {1, 3}, {1, 3}, 0, 0}, {1, {1, 4}, {1, 4}, 0, 0}, {1, {1, 5}, {1, 5}, 0, 0}, {1, {1, 7}, {1, 7}, 0, 1}, {1, {1, 8}, {1, 8}, 0, 2}, {1, {2, 3}, {2, 3}, 0, 0}, {1, {2, 6}, {2, 6}, 0, 0}, {1, {2, 7}, {2, 7}, 0, 0}, {1, {3, 5}, {3, 5}, 0, 1}, {1, {3, 6}, {3, 6}, 0, 0}, {1, {3, 7}, {3, 7}, 0, 0}, {1, {3, 8}, {3, 8}, 0, 2}, {1, {4, 5}, {4, 5}, 0, 0}, {1, {4, 6}, {4, 6}, 0, 0}, {1, {4, 7}, {4, 7}, 0, 0}, {1, {5, 6}, {5, 6}, 0, 0}, {1, {5, 7}, {5, 7}, 0, 0}, {1, {5, 8}, {5, 8}, 0, 2}, {1, {6, 7}, {6, 7}, 0, 0}, {1, {7, 8}, {7, 8}, 0, 2}, {2, {0, 1, 2}, {0, 1, 2}, 0, 0}, {2, {0, 1, 3}, {0, 1, 3}, 0, 0}, {2, {0, 1, 4}, {0, 1, 4}, 0, 0}, {2, {0, 1, 5}, {0, 1, 5}, 0, 1}, {2, {0, 2, 3}, {0, 2, 3}, 0, 0}, {2, {0, 2, 4, 6}, {0, 2, 4, 6}, 0, 0}, {2, {0, 3, 5}, {0, 3, 5}, 0, 0}, {2, {0, 4, 5}, {0, 4, 5}, 0, 0}, {2, {1, 2, 3}, {1, 2, 3}, 0, 1}, {2, {1, 2, 7}, {1, 2, 7}, 0, 0}, {2, {1, 3, 8}, {1, 3, 8}, 0, 0}, {2, {1, 4, 5}, {1, 4, 5}, 0, 0}, {2, {1, 4, 7}, {1, 4, 7}, 0, 0}, {2, {1, 5, 8}, {1, 5, 8}, 0, 0}, {2, {1, 7, 8}, {1, 7, 8}, 0, 0}, {2, {2, 3, 6}, {2, 3, 6}, 0, 0}, {2, {2, 3, 7}, {2, 3, 7}, 0, 0}, {2, {2, 6, 7}, {2, 6, 7}, 0, 0}, {2, {3, 5, 6}, {3, 5, 6}, 0, 0}, {2, {3, 5, 8}, {3, 5, 8}, 0, 0}, {2, {3, 6, 7}, {3, 6, 7}, 0, 1}, {2, {3, 7, 8}, {3, 7, 8}, 0, 0}, {2, {4, 5, 6}, {4, 5, 6}, 0, 0}, {2, {4, 5, 7}, {4, 5, 7}, 0, 1}, {2, {4, 6, 7}, {4, 6, 7}, 0, 0}, {2, {5, 6, 7}, {5, 6, 7}, 0, 0}, {2, {5, 7, 8}, {5, 7, 8}, 0, 0}, {3, {0, 1, 2, 3}, {0, 1, 2, 3}, 0, 1}, {3, {0, 1, 2, 4, 6, 7}, {0, 1, 2, 4, 6, 7}, 0, 1}, {3, {0, 1, 3, 5, 8}, {0, 1, 3, 5, 8}, 0, 1}, {3, {0, 1, 4, 5}, {0, 1, 4, 5}, 0, 1}, {3, {0, 2, 3, 4, 5, 6}, {0, 2, 3, 4, 5, 6}, 0, 1}, {3, {1, 2, 3, 7, 8}, {1, 2, 3, 7, 8}, 0, 1}, {3, {1, 4, 5, 7, 8}, {1, 4, 5, 7, 8}, 0, 1}, {3, {2, 3, 6, 7}, {2, 3, 6, 7}, 0, 1}, {3, {3, 5, 6, 7, 8}, {3, 5, 6, 7, 8}, 0, 1}, {3, {4, 5, 6, 7}, {4, 5, 6, 7}, 0, 1}, {4, {0, 1, 2, 3, 4, 5, 6, 7, 8}, {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, 1, 0}}
  elapsedTime Xs = findAllCYs Q;
  twoD = partition(restrictTriangulation, Xs);
  # keys oo -- 2
  netList annotatedFaces Q

  X1 = first twoD#((keys twoD)#0)
  X2 = first twoD#((keys twoD)#1)
  gv1 = partitionGVConeByGV(X1, DegreeLimit => 40)
  gv2 = partitionGVConeByGV(X2, DegreeLimit => 40)
  Ms = findLinearMaps(gv1, gv2)

  -- Now we need to see if these two topologies are the same.
  Ms = Ms/(m -> lift(m, ZZ))
  Ms/det
  phis = for m in Ms list map(RZ, RZ, m)

  RZ = ZZ[a,b,c,d,e]
  T1 = topologicalData(X1, RZ)
  T2 = topologicalData(X2, RZ)
  (L1, F1) = (c2 T1, cubicForm T1)
  (L2, F2) = (c2 T2, cubicForm T2)
  phis_0 L1 - L2
  phis_0 F1 - F2

  phis_1 L1 - L2 -- 0
  phis_1 F1 - F2 -- 0: so both topologies are the same!

  phis_2 L1 - L2 -- 0
  phis_2 F1 - F2 -- 0: so both topologies are the same!

  phis_3 L1 - L2 -- 0
  assert(phis_3 F1 - F2 != 0) -- this doesn't give an isomorphism...

  -- it seems from the stuff below, that maybe I thought they were different before?

  elapsedTime Ts = findAllFRSTs Q;
  Xs = for i from 0 to #Ts - 1 list cyData(Q, Ts#i, ID => i);
  assert(#Xs == 142) -- There are 142 triangulations 
  -- These 142 triangulations have at most 2 different topologies.
  RZ = ZZ[a,b,c,d,e]
  H = partition(X -> topologicalData(X, RZ), Xs);
  # keys H
  H#(first keys H)/(X -> X.cache#"id") -- 0, 1, 4, 5, 6, ...
  H#(last keys H)/(X -> X.cache#"id") -- 2, 3, 13, ...
-- Here we pick out two of these (Xs#0, Xs#2), one of each possibly different "topology".  

restart
debug needsPackage "StringTorics"
Qrays = {{-1, -1, 0, 0}, {-1, -1, 0, 1}, {-1, 0, -1, 0}, {-1, 1, -1, 0}, 
    {0, -1, 0, -1}, {0, -1, 1, -1}, {0, 0, -1, -1}, {1, 0, -1, -1}, {1, 1, 1, 1}}
basisInd = {0, 1, 2, 4, 5}
GLSM = {{1, 0, 0, 0, 0}, {0, 1, 0, 0, 0}, {0, 0, 1, 0, 0}, {0, -1, -1, 1, 2}, 
    {0, 0, 0, 1, 0}, {0, 0, 0, 0, 1}, {1, 5, 2, -2, -5}, {0, -2, -1, 1, 3}, {1, 2, 1, 0, -1}}
T1 = {{0, 1, 2, 3}, {0, 1, 2, 4}, {0, 1, 3, 5}, {0, 1, 4, 5}, {0, 2, 3, 4}, 
    {0, 3, 4, 5}, {1, 2, 3, 7}, {1, 2, 4, 6}, {1, 2, 6, 7}, {1, 3, 5, 8}, 
    {1, 3, 7, 8}, {1, 4, 5, 7}, {1, 4, 6, 7}, {1, 5, 7, 8}, {2, 3, 4, 6}, 
    {2, 3, 6, 7}, {3, 4, 5, 6}, {3, 5, 6, 7}, {3, 5, 7, 8}, {4, 5, 6, 7}}
T2 = {{0, 1, 2, 3}, {0, 1, 2, 6}, {0, 1, 3, 5}, {0, 1, 4, 5}, {0, 1, 4, 6}, 
    {0, 2, 3, 6}, {0, 3, 5, 6}, {0, 4, 5, 6}, {1, 2, 3, 7}, {1, 2, 6, 7}, 
    {1, 3, 5, 8}, {1, 3, 7, 8}, {1, 4, 5, 7}, {1, 4, 6, 7}, {1, 5, 7, 8}, 
    {2, 3, 6, 7}, {3, 5, 6, 7}, {3, 5, 7, 8}, {4, 5, 6, 7}}
Q = cyPolytopeData convexHull transpose matrix Qrays
assert(degrees Q === GLSM)
assert(basisIndices Q == basisInd)
X1 = cyData(Q, T1)
X2 = cyData(Q, T2)

mori1 = toricMoriCone X1
  assert(transpose entries rays mori1 == {
          {1, 0, -1, -1, 0}, {-1, 1, -1, 1, 0}, {0, -1, 2, 1, 0}, {
              0, -1, 0, 5, -3}, {0, 0, 0, 0, -1}, {-1, 0, 2, -1, 1}})
mori2 = toricMoriCone X2
  assert(transpose entries rays mori2 == {
          {0, 1, -2, 0, 0}, {2, -1, 0, -1, 0}, {0, 0, 0, 0, -1}, 
          {-1, 0, 0, 3, -1}, {-2, 0, 3, 0, 1}, {-1, 1, 0, -1, 1}})
RZ = ZZ[a,b,c,d,e]
topologicalData(X1, RZ)
topologicalData(X2, RZ)

-- now find heft vectors (weight vectors) for the curve classes...
-- Here is the lazy way
weights1 = sum entries transpose rays dualCone posHull transpose matrix rays mori1
weights2 = sum entries transpose rays dualCone posHull transpose matrix rays mori2

-- By hand, find a somewhat better wts1.
CM1 = posHull transpose matrix rays mori1
CM1 = mori1
DM1 = dualCone CM1
(transpose rays DM1) * (rays CM1)
DM = entries transpose rays DM1
rays DM1
wts1 = DM_1 + DM_2 + DM_6 + DM_3
assert(wts1 == {20, 30, 14, 5, -2})
-- or -- wts1 = DM_1 + DM_2 + DM_6 + DM_5
matrix{ wts1} * rays CM1
--wts1 = {22, 32, 15, 6, -1}


-- By hand, find a somewhat better wts2 for mori2
CM2 = mori2
DM2 = dualCone CM2
(transpose rays DM2) * (rays CM2)
DM = entries transpose rays DM2
wts2 = DM_0 + DM_1 + DM_2  + DM_5
wts2 = DM_0 + DM_1 + DM_4  + DM_5
matrix{ wts2} * rays CM2
assert(wts2 == {16, 25, 12, 6, -1})
hmori1 = hilbertBasisGenerators mori1
hmori2 = hilbertBasisGenerators mori2

gv1 = gvInvariants(X1, Heft => wts1, DegreeLimit => 30);
gv2 = gvInvariants(X2, Heft => wts2, DegreeLimit => 30);


gv1 = gvInvariants(ambient X1, basisInd, Mori => hmori1, Heft => wts1, DegreeLimit => 30);
gv2 = gvInvariants(ambient X2, basisInd, Mori => hmori2, Heft => wts2, DegreeLimit => 30);

gv1 = gvInvariants(ambient X1, basisInd, Heft => wts1, DegreeLimit => 30); -- not yet working...  
gv2 = gvInvariants(ambient X2, basisInd, Heft => wts2, DegreeLimit => 30);

gvCone1 = (keys gv1)/toList//matrix//transpose//posHull//rays
gvCone2 = (keys gv2)/toList//matrix//transpose//posHull//rays
mv1 = partition(x -> gv1#(toSequence x), entries transpose gvCone1)
mv2 = partition(x -> gv2#(toSequence x), entries transpose gvCone2)

T = QQ[t_(1,1) .. t_(5,5)]
M = genericMatrix(T, 5, 5)
mv1#-2/(x -> transpose matrix {x})
mv2#-2/(x -> transpose matrix {x})

vecs1 = mv1#1/(x -> transpose matrix {x})
vecs2 = mv2#1/(x -> transpose matrix {x})

vecs1b = mv1#-2/(x -> transpose matrix {x})
vecs2b = mv2#-2/(x -> transpose matrix {x})

findmat = (p,q) -> (
    trim(sum for i from 0 to #p-1 list ideal( M * (vecs1#i) - vecs2#(p#i))
    +
    sum for i from 0 to #q-1 list ideal( M * (vecs1b#i) - vecs2b#(q#i))
    ))
findmat({0,1,2,3,4}, {0,1})
findmat({0,1,2,3,4}, {1,0})
for p in subsets(7, 5) list findmat(p, {0,1}) -- all 1's
for p in subsets(7, 5) list findmat(p, {1,0}) -- all 1's.  Therefore, it appears that these two
  -- topologies are different...

dualCone posHull transpose matrix mori2
rays oo
(transpose o25) * o22
transpose (o25_{0} + o25_{1} + o25_{2} + o25_{5}) * o22
GV1 = gvInvariants(ambient Xs#0, basisIndices Xs#0, Mori => mori1, DegreeLimit => 40)
GV2 = gvInvariants(ambient Xs#2, basisIndices Xs#2, Mori => mori2, DegreeLimit => 40, Heft => {11, 17, 8, 4, -1})
rays posHull transpose matrix (keys GV1/toList)
rays posHull transpose matrix (keys GV2/toList)
F1 = cubicForm first keys H
F2 = cubicForm last keys H
factor(F1-F2)


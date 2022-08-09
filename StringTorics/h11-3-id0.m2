-- This example fails several things, due to its class group being torsion.
restart
loadPackage "StringTorics"

-- First polytope in h11 = 3 list.
-- 4 5  M:53 5 N:9 5 H:3,43 [-80]
polystr = /// 
         1   0   2   4 -10
         0   1   3   5  -9
         0   0   4   0  -4
         0   0   0   8  -8
///

A = matrixFromString polystr
P = convexHull A
P2 = polar P
dim P
dim P2

-- Analyze the polytope
isSimplicial P

-- Triangulate it
(LP,tri) = regularStarTriangulation(3,P2)
(LP,tri) = regularStarTriangulation(2,P2)
regularStarTriangulation P2
-- Basics of the normal toric variety
V = normalToricVariety(LP, tri)
debugLevel = 2
assert isWellDefined V
assert(not isSmooth V) -- no
assert isSimplicial V -- true
assert isProjective V -- true
assert(rays V == LP)
picardGroup V == ZZ^3
classGroup V 

-- Analyze Hodge numbers of the corresponding CY3-fold X
assert(3 == h11OfCY P)
assert(43 == h21OfCY P)
assert isFavorable P
assert(hodgeOfCYToricDivisors P
    ===
    hashTable {0 => {1,0,0}, 1 => {1,0,0}, 2 => {1,0,0}, 3 => {1,0,0},
        4 => {1,0,9}, 5 => {1,2,0}, 6 => {1,2,0}, 7 => null}
    )
netList annotatedFaces(0,P2)
netList annotatedFaces(1,P2)
netList annotatedFaces(2,P2)
netList annotatedFaces(3,P2)

-- the following doesn't work, since we can't do gradings by torsion groups yet
for i from 0 to 6 list for j from 0 to 3 list HH^j(V, OO V_i)
-- so let's try cohomcalg:
KV = toricDivisor V
X = completeIntersection(V, {-KV})
hodgeDiamond X -- fails
for i from 0 to 6 list cohomologyVector(X, degree(KV-KV)) -- not working yet
cohomCalgVector (KV-KV) -- not working yet.
for i from 0 to 6 list hodgeVector({-KV, V_i})

-- Hodge-Deligne numbers for X

-- 

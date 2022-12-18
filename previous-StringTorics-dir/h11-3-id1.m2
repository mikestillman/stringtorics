restart
needsPackage "StringTorics"
-- id=1 polytope in h11 = 3 list.
-- 4 5  M:48 5 N:8 5 H:3,45 [-84]
polystr = ///  1   0   2   4  -8
         0   1   5   3  -9
         0   0   6   0  -6
         0   0   0   6  -6
///

A = matrixFromString polystr
P = convexHull A
P2 = polar P
dim P
dim P2

-- Analyze the polytope
isSimplicial P2

V = reflexiveToSimplicialToricVariety P
debugLevel = 2
assert isWellDefined V
debugLevel = 0
assert(not isSmooth V) -- no
assert isSimplicial V -- true
assert isProjective V -- true
pic V == ZZ^3
assert(cl V == pic V) -- so this example has no torsion

-- Analyze Hodge numbers of the corresponding CY3-fold X
assert(3 == h11OfCY P)
assert(45 == h21OfCY P)
assert isFavorable P
assert(hodgeOfCYToricDivisors P
    ===
    hashTable {
        0 => {1,0,1}, 
        1 => {1,0,1}, 
        2 => {1,0,0}, 
        3 => {1,0,0},
        4 => {1,0,4},
        5 => {1,4,0}, 
        6 => {1,4,0}
        }
    )
netList annotatedFaces(0,P2)
netList annotatedFaces(1,P2)
netList annotatedFaces(2,P2)
netList annotatedFaces(3,P2)

-- the following doesn't work, since we can't do gradings by torsion groups yet
-- but in any case, it isn't quite what we want right now anyway
matrix for i from 0 to 6 list for j from 0 to 3 list rank HH^j(V, OO V_i)

-- so let's try cohomcalg:
KV = toricDivisor V
X = completeIntersection(V, {-KV})
assert(hodgeDiamond X == matrix{
    {1,0,0,1},
    {0,3,45,0},
    {0,45,3,0},
    {1,0,0,1}}
    )
for i from 0 to 6 list cohomologyVector(X, degree(V_i)) -- ok
for i from 0 to 6 do assert(hodgeVector({-KV, V_i}) == (hodgeOfCYToricDivisors P)#i)

-- Intersection theory on X
    pt = base(symbol a, symbol b, symbol c, symbol d)
    Va = abstractVariety(V, pt)
    Xa = abstractVariety(X, pt)
    IX = intersectionRing Xa
    use IX
    F1 = chi OO(a*t_2 + b*t_3 + c*t_6)
    time zeroset1 = allZeros(F1-1, 2, (-20, 20))
    hashTable for d in zeroset1 list d => cohomologyVector(X, d) -- fails, BUG, due to #variables needed....
    cohomologyVector(X, zeroset1_0)    
    cohomologyVector(X, zeroset1_1)
    cohomologyVector(X, zeroset1_2)
    cohomologyVector(X, zeroset1_3)    
    cohomologyVector(X, zeroset1_4)
    cohomologyVector(X, zeroset1_5)
    cohomologyVector(X, zeroset1_6) -- rigid
    cohomologyVector(X, zeroset1_7) -- rigid
    cohomologyVector(X, zeroset1_8) -- rigid {0,1,1}
    cohomologyVector(X, zeroset1_9) -- rigid {0,1,1}
    cohomologyVector(X, zeroset1_10) -- rigid {-2,3,2}
    cohomologyVector(X, zeroset1_11) -- rigid {-3,4,2}
    cohomologyVector(X, zeroset1_12) -- 
    cohomologyVector(X, zeroset1_13) -- 
    cohomologyVector(X, zeroset1_14) --
    cohomologyVector(X, zeroset1_15) 
    cohomologyVector(X, zeroset1_16) -- FAILS, too many variables in use...

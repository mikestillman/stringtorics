restart
--load "h11-template.m2"
needsPackage "StringTorics"
--polystrs = parseKS get "h11-6-KS.txt";
polystrs = parseKS getKreuzerSkarke 6;
A = matrixFromString last polystrs_14
P = convexHull A
P2 = polar P

V = reflexiveToSimplicialToricVariety P
Va = abstractVariety V
A = intersectionRing Va
KV = toricDivisor V

netList annotatedFaces P2

isWellDefined V
Ds = for i from 0 to #rays V - 1 list V_i
matrix for D in Ds list cohomologyVector(V, D)
matrix for D in Ds list cohomologyVector(V, D+KV)
select(subsets(Ds,2), e -> (c := cohomologyVector(V, sum e + KV); sum c > 0))
select(subsets(Ds,4), e -> (c := cohomologyVector(V, sum e + KV); sum c > 0))
for e in oo list e => cohomologyVector(V, sum e + KV)
hashTable for e in subsets(for i from 0 to #rays V - 1 list V_i, 2) list e#0 + e#1 => (cohomologyVector(V, e#0 + e#1), cohomologyVector(V, e#0 + e#1 + KV))
hashTable for e in subsets(for i from 0 to #rays V - 1 list V_i, 3) list sum e => cohomologyVector(V, sum e)
hashTable for e in subsets(for i from 0 to #rays V - 1 list V_i, 8) list sum e => cohomologyVector(V, -sum e)
subcomplex(V,{4,5,6,7})
D = completeIntersection(V, {-KV, V_1})
    cohomologyVector(D, degree(V_0-V_0))
cohomologyVector(X, degree(V_7+V_8))
cohomologyVector(V, degree(V_7+V_8))
cohomologyVector(V, degree(V_7+V_8+KV))
--for i from 0 to 2 list HH^i(V, MS)
cohomologyVector(V, -KV)
isAmple(-KV)
hashTable for e in subsets(for i from 0 to #rays V - 1 list V_i, 2) list sum e => isAmple(-4*KV + sum e)
for i from 0 to numgens A - 1 list integral( (2*( (-A_8-A_9) + 2 * sum gens A)^3) * A_i)
isAmple((2*-V_8-V_9) - 2*KV)
isNef(2*(-V_8-V_9) - 2*KV)

-- nef cone "by hand":
-- use V_4, V_5, V_6, V_7, V_8, V_9
orbits(V,1)
matrix {for f in orbits(V,1) list f/(i -> A_i)//product}
entries((transpose matrix{{A_4, A_5, A_6, A_7, A_8, A_9}}) * oo)
oo/(L -> L/integral)
matrix oo
M = lift(oo,QQ)
C = dualCone posHull M
MC = rays C
isNef(V_6+V_7)
assert isNef(V_6+V_7+V_9)
assert isNef(V_5+V_6+V_7+V_9)
assert isNef(V_5+V_6+V_7+V_9+V_4+V_8)
(transpose MC) * M
D = sum {3*V_4, 6*V_5, 13*V_6, 13*V_7, 4*V_8, 12*V_9}
D = sum {3*V_4, 6*V_5, 13*V_6, 13*V_7, 4*V_8, 12*V_9} - sum{V_4,V_5,V_6,V_7,V_8,V_9}
cohomologyVector(V,D)
select(1..40, i ->  isCartier (i*D))
isCartier(24*D)
isAmple(24*D) -- yes
isAmple D-- incorrect....
isVeryAmple(24*D)
{3*V_4, 6*V_5, 13*V_6, 13*V_7, 4*V_8, 12*V_9}/(d -> cohomologyVector(V,d))
restart
--load "h11-template.m2"
needsPackage "StringTorics"
--polystrs = parseKS get "h11-6-KS.txt";
polystrs = parseKS getKreuzerSkarke 6;
A = matrixFromString last polystrs_14
P = convexHull A
P2 = polar P

V = reflexiveToSimplicialToricVariety P
Va = abstractVariety V
A = intersectionRing Va
KV = toricDivisor V

netList annotatedFaces P2

isWellDefined V
matrix for i from 0 to #rays V - 1 list cohomologyVector(V, V_i)
matrix for i from 0 to #rays V - 1 list cohomologyVector(V, V_i+KV)
hashTable for e in subsets(for i from 0 to #rays V - 1 list V_i, 2) list e#0 + e#1 => (cohomologyVector(V, e#0 + e#1), cohomologyVector(V, e#0 + e#1 + KV))
hashTable for e in subsets(for i from 0 to #rays V - 1 list V_i, 3) list sum e => cohomologyVector(V, sum e)
hashTable for e in subsets(for i from 0 to #rays V - 1 list V_i, 8) list sum e => cohomologyVector(V, -sum e)
D = completeIntersection(V, {-KV, V_1})
    cohomologyVector(D, degree(V_0-V_0))
cohomologyVector(X, degree(V_7+V_8))
cohomologyVector(V, degree(V_7+V_8))
cohomologyVector(V, degree(V_7+V_8+KV))
--for i from 0 to 2 list HH^i(V, MS)
cohomologyVector(V, -KV)
isAmple(-KV)
hashTable for e in subsets(for i from 0 to #rays V - 1 list V_i, 2) list sum e => isAmple(-4*KV + sum e)
for i from 0 to numgens A - 1 list integral( (2*( (-A_8-A_9) + 2 * sum gens A)^3) * A_i)
isAmple((2*-V_8-V_9) - 2*KV)
isNef(2*(-V_8-V_9) - 2*KV)

-- nef cone "by hand":
-- use V_4, V_5, V_6, V_7, V_8, V_9
orbits(V,1)
matrix {for f in orbits(V,1) list f/(i -> A_i)//product}
entries((transpose matrix{{A_4, A_5, A_6, A_7, A_8, A_9}}) * oo)
oo/(L -> L/integral)
matrix oo
M = lift(oo,QQ)
C = dualCone posHull M
rays C
isNef(V_6+V_7+V_9)

tally values oo
S = ring V
basis(degree(x_8), Ext^1(ideal V, S))

SM = newRing(ring V, DegreeRank=>numgens ring V)
BM = sub(ideal V, SM)
basis({0,0,0,0,0,0,0,0,1,0},Ext^1(BM^[4],SM))
basis({0,0,0,0,0,0,0,0,0,1},Ext^1(BM^[4],SM))
F = random(-degree KV, S)
M = comodule ideal(F, x_0)
MS = sheaf(V, M)
HH^0(V, MS)
HH^1(V, MS)
--elapsedTime HH^2(V, MS)

M = comodule saturate(ideal(F, x_1), ideal V)
numColumns basis({0,0,0,0,0},Ext^1(ideal V, M))
numColumns basis({0,0,0,0,0},Ext^1(module (ideal V)^[2], M))
numColumns basis({0,0,0,0,0},Ext^1(module (ideal V)^[3], M))
numColumns basis({0,0,0,0,0},Ext^2(comodule (ideal V)^[4], M))
numColumns basis({0,0,0,0,0},time Ext^2(comodule (ideal V)^[5], M))
X = completeIntersection(V, {-KV})
hashTable for i from 0 to 9 list i => cohomologyVector(X, degree V_i)
hashTable for i from 0 to 9 list i => {cohomologyVector(X, - degree V_i), cohomologyVector(X, - 2 * degree V_i), cohomologyVector(X, - 3 * degree V_i)}
isAmple D
tally values oo
S = ring V
basis(degree(x_8), Ext^1(ideal V, S))

SM = newRing(ring V, DegreeRank=>numgens ring V)
BM = sub(ideal V, SM)
basis({0,0,0,0,0,0,0,0,1,0},Ext^1(BM^[4],SM))
basis({0,0,0,0,0,0,0,0,0,1},Ext^1(BM^[4],SM))
F = random(-degree KV, S)
M = comodule ideal(F, x_0)
MS = sheaf(V, M)
HH^0(V, MS)
HH^1(V, MS)
--elapsedTime HH^2(V, MS)

M = comodule saturate(ideal(F, x_1), ideal V)
numColumns basis({0,0,0,0,0},Ext^1(ideal V, M))
numColumns basis({0,0,0,0,0},Ext^1(module (ideal V)^[2], M))
numColumns basis({0,0,0,0,0},Ext^1(module (ideal V)^[3], M))
numColumns basis({0,0,0,0,0},Ext^2(comodule (ideal V)^[4], M))
numColumns basis({0,0,0,0,0},time Ext^2(comodule (ideal V)^[5], M))
X = completeIntersection(V, {-KV})
hashTable for i from 0 to 9 list i => cohomologyVector(X, degree V_i)
hashTable for i from 0 to 9 list i => {cohomologyVector(X, - degree V_i), cohomologyVector(X, - 2 * degree V_i), cohomologyVector(X, - 3 * degree V_i)}

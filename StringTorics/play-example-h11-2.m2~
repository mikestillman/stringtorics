needsPackage "StringTorics"
Ps = kreuzerSkarke(2, Access => "wget")

A = matrix Ps_2
P1 = convexHull A
P2 = polar P1
V = reflexiveToSimplicialToricVariety P1
rays V -- 6 rays
ring V
transpose matrix degrees ring V

h11OfCY P1
h21OfCY P1

h11OfCY P2
h21OfCY P2

-- too hard for the moment.
V2 = reflexiveToSimplicialToricVariety P2
--ring V2;

-- Cohomology on V itself.
HH^1(V, OO_V(1,0))
cohoms = deg -> for i from 0 to 3 list rank HH^i(V, OO_V deg)
cohoms(1,0)
cohoms(2,1)
L = flatten for a from -10 to 10 list for b from -10 to 10 list (a,b) => cohoms(a,b);
select(L, x -> x#1#0 == 1)
netList oo

KV = toricDivisor V
X = completeIntersection(V, {-KV})
cohomologyVector(X, {1,0})
for d in degrees ring V list d => cohomologyVector(X, d)
L = flatten for a from -10 to 10 list for b from -10 to 10 list (a,b) => cohomologyVector(X, {a,b});

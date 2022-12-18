path=append(path,"/mnt/d/git_repo/stringtorics")
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
dual monomialIdeal (ideal V)
ideal V
HH^1(V,OO_V(-4,-1))
HH^1(V,OO_V(0,1))
posHull transpose matrix {{-1,0},{2,1},{-3,-1},{0,1}} 
rays oo
matrix for i from -10 to 0 list for j from -4 to 4 list rank HH^1(V,OO_V(i,j))
matrix for i from -10 to 0 list for j from -4 to 4 list rank HH^2(V,OO_V(i,j))
matrix for i from -10 to 0 list for j from -4 to 4 list rank HH^3(V,OO_V(i,j))
matrix for i from -10 to 0 list for j from -4 to 4 list rank HH^4(V,OO_V(i,j))

-- too hard for the moment.
V2 = reflexiveToSimplicialToricVariety P2

--ring V2;

-- Cohomology on V itself.
HH^1(V, OO_V(1,0))
cohoms = deg -> for i from 0 to 3 list rank HH^i(V, OO_V deg)
cohoms(1,0)
cohoms(2,1)
L= flatten for a from -10 to 10 list for b from -10 to 10 list (a,b) => sum cohoms(a,b);
L= matrix for a from -10 to 10 list for b from -10 to 10 list sum cohoms(a,b);

matrix (L/last)
select(L, x -> x#1#0 == 1)
netList oo

KV = toricDivisor V
X = completeIntersection(V, {-KV})
cohomologyVector(X, {1,0})
for d in degrees ring V list d => cohomologyVector(X, d)
L = flatten for a from -10 to 10 list for b from -10 to 10 list (a,b) => cohomologyVector(X, {a,b});

-- Example usage of StringTorics, and other facilities in Macaulay2
-- for toric varieties, and Calabi-Yau hypersurfaces on them.

-- Packages installed include:
--  NormalToricVarieties
--  Polyhedra
--  Triangulations
--  Schubert2
--  ReflexivePolytopesDB
--  CohomCalg
--  Topcom

restart
needsPackage "StringTorics"

------------------------------
-- A. Reflexive polytopes ----
-- The main type is `ReflexivePolytope`
-- Also: Polyhedron
------------------------------

-- We use the incredibly useful Kreuzer-Skarke database!

topes = kreuzerSkarke(4, Limit => 10000);
#topes
tope = topes_50
Q = reflexive tope
transpose matrix vertices Q -- "polar dual" to convexHull of tope
isFavorable Q
rays Q -- Q doesn't have "rays"!  These are rays of the corresponding fan in N lattice
transpose matrix oo
degrees Q
transpose matrix oo
hh^(1,1) Q
hh^(1,2) Q

P = polytope(Q, "N") 
vertices P -- permuted!
isCompact P
isReflexive P
latticePoints P
latticePoints polar P

------------------------------
-- B. Triangulations ---------
-- main type: Triangulation
------------------------------
A = transpose matrix rays Q
options allTriangulations
methods allTriangulations
Ts = allTriangulations(A | transpose matrix{{0,0,0,0}}, Fine => true) -- triangulations of the point set.
#Ts
Ts = select(Ts, isStar)
#Ts -- 2 of these
t = triangulation(A, apply(max Ts_0, x -> drop(x, -1)))

Vs = allTriangulations(A, Homogenize => false)
Vs = select(Vs, isFine)
#Vs -- 28

member(t, Vs)

------------------------------
-- C. Making CalabiYau 3-fold hypersurfaces
-- main type: CalabiYauInToric
------------------------------

options findAllCYs -- finds all "Batyrev" Calabi-Yaus with this polytope.
Xs = findAllCYs(Q, PicardRing => ZZ[u_0..u_3])
X = first oo
dim X
hh^(1,1) X
hh^(1,2) X
peek X
max X
rays X

-- Same data as a NormalToricVariety!
-- (except must come from reflexive polytope)
V = normalToricVariety X
max V
rays V
isSmooth V
isSimplicial V
isProjective V
picardGroup V

c2Form X
cubicForm X
RQ = QQ[x,y,z,w]
F = sub(cubicForm X, vars RQ)
jacobian F
saturate ideal jacobian F -- F=0 is smooth!

-- Question: F=0 must be the blowup of 6 points.  Which points can correspond to the cubic form
-- of a h11=4 Calabi-Yau 3-fold

-- Let's grab some h11=3 examples
topes = kreuzerSkarke 3;
#topes
Q = reflexive topes_60
hh^(1,1) Q, hh^(1,2) Q

V = smoothFanoToricVariety(4, 30)
picardGroup V
Q = reflexive rays V
set rays Q == set rays V

A = transpose matrix rays Q
needsPackage "PALPInterface"
normalForm A

Xs = findAllCYs Q
X = first Xs
F = cubicForm X
factor F
RQ = QQ[gens ring F]
FQ = sub(F, RQ)
saturate ideal jacobian FQ -- singular at a point

------------------------
-- D. CY3, Gopakumar-Vafa/Gromov-Witten numbers.
-- type: CY3 (topological + GV information for a CY3).
------------------------

GVT = gvInvariants(X, DegreeLimit => 20)
GVR = gvRayTable GVT

X1 = makeCY3(X, GVTable => GVT, DegreeLimit => 20)
extremalCurves X1

X2 = performFlop(X1, {0, 0, 1})
extremalCurves X2

findEquivalence({c2Form X1, cubicForm X1}, {c2Form X2, cubicForm X2})

X3 = performFlop(X2, {0,1,3})
extremalCurves X3

--------------------------
-- E. Line bundle cohomology
-- type: CompleteIntersectionInToric
--------------------------
Xci = completeIntersection(V, {-toricDivisor V})
hodgeDiamond Xci
hh^* OO_Xci(-1,1,0)
hh^0 OO_Xci(-1,1,0)
effV = posHull transpose matrix degrees ring V
rays oo -- effective cone of V

nonzeros = for x in (-1,-1,-1)..(1,1,1) list (
    r := hh^* OO_Xci x;
    if r#0 > 0 then x => r else continue
    )
netList nonzeros
effX = posHull transpose matrix(nonzeros/first/toList)
contains(effV, effX) -- no autotochthonous divisors found yet
effV == effX

nonzeros = for x in (-2,-2,-2)..(2,2,2) list (
    r := hh^* OO_Xci x;
    if r#0 > 0 then x => r else continue
    )
effFound = transpose matrix(nonzeros/first/toList)
(halfspaces effV) * effFound -- all are in this cone.

------------------------------
-- F. AbstractVariety --------
-- type: AbstractVariety
------------------------------
X = makeCY Q
transpose matrix degrees Q
hh^(1,1) X
hh^(1,2) X
Xa = abstractVariety(X, base(a_0..a_2))
IX = intersectionRing Xa
tangentBundle Xa
p = chern_2 tangentBundle Xa
basisIndices X
integral(p * t_0)
integral(p * t_1)
integral(p * t_2)
c2 X
chi OO(t_0) -- sum (-1)^i hh^i(OO_X(X_0))
chi OO(t_1) -- sum (-1)^i hh^i(OO_X(X_0))
chi OO(t_1-t_0) -- sum (-1)^i hh^i(OO_X(X_0)) 


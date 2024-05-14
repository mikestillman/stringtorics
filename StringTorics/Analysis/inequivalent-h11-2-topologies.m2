-- I get:
-- 36 different polytopes
-- 39 different NTFE CY3's
-- 29 distinct topologies

restart
debug needsPackage "StringTorics" -- the debug is because some functions are not yet exported.
DBNAME = "../Databases/cys-ntfe-h11-2.dbm"

R = ZZ[a,b]
RQ = QQ (monoid R);
(Qs, Xs) = readCYDatabase(DBNAME, Ring => R);
sort keys Xs
tally for lab in keys Qs list #automorphisms Qs#lab

-- 39 different examples.

-- We collect the nontorsion, torsion, favorable, nonfavorable's.
   torsions = for k in keys Qs list (
      istor := prune coker matrix rays Qs#k;
      if not isFreeModule istor then k else continue
      )
   nonfavorables = for k in keys Qs list (
       if not isFavorable Qs#k then k else continue
      )
   favorablesXs = sort for k in keys Xs list (
       if not isFavorable Qs#(first k) then continue else k
      )
   assert(#torsions == 2)
   assert(#nonfavorables == 0)

---------------------------------------------------------------
-- Next step: How many of these 39 are distinct topologies? --
---------------------------------------------------------------
  allXs = sort keys Xs
  allT = topologySet(allXs, Xs);
  info allT -- 39 possibly different topologies, so far
  
  allT1 = combineIfSame(allT, X -> (c2Form X, cubicForm X))
  info allT1 -- down to 35 different

  equivalences allT1

  -- Question: are there any torsions or nonfavorables in here?
  -- Torsions: none on this list.
  -- Nonfavorables: (1059, 0), (1060, 0), (1064, 0), (1065, 0) are all the same (on the nose).
  
  elapsedTime allT2 = separateIfDifferent(allT1, invariantsAll) -- 1 second
  info allT2
  equivalences allT2
  representatives allT2
  
  -- Now we take the buckets consisting of more than one set, and try to make equivalences.
  elapsedTime allT3a = combineByGV(allT2, DegreeLimit => 20); -- 
  info allT3a  
  representatives allT3a
  
  -- only two pairs left to check:
  -- {(9, 0), (10, 0)}
  -- {(28, 0), (31, 0)}
  -- below shows that the data is equivalent over QQ for each pair, but not over ZZ.
  
  (A, phi) = genericLinearMap RQ  
  T = target phi

  Ls = {(9, 0), (10, 0)}
  (X1, X2) = toSequence (Ls/(lab -> Xs#lab))
  (L1, F1) = (c2Form X1, cubicForm X1)
  (L2, F2) = (c2Form X2, cubicForm X2)
  -- check X1,X2: SAME
  JL = ideal last coefficients(phi sub(L1, T) - sub(L2, T))
  JF = ideal last coefficients(phi sub(F1, T) - sub(F2, T))
  J = trim(JL + JF)
  trim J
  decompose oo -- not equivalent, no integral points.
  
  Ls = {(28, 0), (31, 0)}
  (X1, X2) = toSequence (Ls/(lab -> Xs#lab))
  (L1, F1) = (c2Form X1, cubicForm X1)
  (L2, F2) = (c2Form X2, cubicForm X2)
  factor det hessian F1
  factor det hessian F2
  -- check X1,X2: SAME
  Ja = ideal last coefficients(phi sub(b, T) % ideal sub(a, T))
  JL = ideal last coefficients(phi sub(L1, T) - sub(L2, T))
  JF = ideal last coefficients(phi sub(F1, T) - sub(F2, T))
  J = trim(Ja + JL + JF)
  decompose J -- no integral points.
  -- somewhat easier:
  JL = ideal last coefficients(phi sub(L1, T) - sub(L2, T))
  JF = ideal last coefficients(phi sub(F1, T) - sub(F2, T))
  J = trim(JL + JF)
  decompose J

--------------------------------------------------------
-- Now let's add in the ones that come from non-Batyrev fans
-- These arise from searching for all fans, ignoring rays interior to facets (ALLOWED??)
--------------------------------------------------------
restart
debug needsPackage "StringTorics" -- the debug is because some functions are not yet exported.
DBNAME = "../Databases/cys-ntfe-h11-2.dbm"

R = ZZ[a,b]
RQ = QQ (monoid R);
(Qs, Xs) = readCYDatabase(DBNAME, Ring => R);
sort keys Xs

prevTops = for lab in sort keys Xs list {hh^(1,1) Xs#lab, hh^(1,2) Xs#lab, c2Form Xs#lab, cubicForm Xs#lab}
labs128 = for lab in sort keys Xs list if hh^(1,2) Xs#lab == 128 then lab else continue
tops128 = for lab in labs128 list {hh^(1,1) Xs#lab, hh^(1,2) Xs#lab, c2Form Xs#lab, cubicForm Xs#lab}
--tops128 = select(prevTops, x -> x#1 == 128)
netList oo
labs128 == {(27, 0), (28, 0), (29, 0), (30, 0), (31, 0)}

newTVs = {
    {10, {{-1, -1, -1, 0}, {-1, -1, -1, 1}, {-1, -1, 0, 0}, {-1, 0, -1, 1}, {0, -1, -1, 0}, {2, 2, 2, -1}},
        {{0, 2, 3, 4}, {0, 2, 3, 5}, {0, 2, 4, 5}, {0, 3, 4, 5}, {1, 2, 3, 4}, {1, 2, 3, 5}, {1, 2, 4, 5}, {1, 3, 4, 5}}},
    {11, {{-1, 0, -1, 0}, {-1, 0, -1, 1}, {-1, 0, 0, 0}, {-1, 1, -1, 2}, {1, -1, -1, -1}, {1, 0, 2, -1}},
        {{0, 2, 3, 4}, {0, 2, 3, 5}, {0, 2, 4, 5}, {0, 3, 4, 5}, {1, 2, 3, 4}, {1, 2, 3, 5}, {1, 2, 4, 5}, {1, 3, 4, 5}}},
    {13, {{-1, -1, 0, -1}, {-1, -1, 0, 0}, {-1, -1, 1, -1}, {-1, 2, -1, 3}, {0, -1, 0, -1}, {2, 0, 0, -1}},
        {{0, 2, 3, 4}, {0, 2, 3, 5}, {0, 2, 4, 5}, {0, 3, 4, 5}, {1, 2, 3, 4}, {1, 2, 3, 5}, {1, 2, 4, 5}, {1, 3, 4, 5}}},
    {30, {{-1, 0, -1, 0}, {-1, 0, -1, 1}, {-1, 0, 0, 0}, {-1, 1, -1, 1}, {-1, 2, 4, -1}, {1, -1, -1, 0}},
        {{0, 2, 3, 4}, {0, 2, 3, 5}, {0, 2, 4, 5}, {0, 3, 4, 5}, {1, 2, 3, 4}, {1, 2, 3, 5}, {1, 2, 4, 5}, {1, 3, 4, 5}}}}

newVs = for x in newTVs list (
    V0 := normalToricVariety(toSequence drop(x, 1), CoefficientRing => ZZ/32003);
    Q := Qs#(x#0);
    degrees Q;
    V0.cache#"basis indices" = basisIndices Q; -- TODO: only works for favorables
    V0
    )

for V in newVs list isWellDefined V
for V in newVs list isSmooth V -- 3 true, #30 is false
newXs = for V in newVs list (
    X0 := completeIntersection(V, {-toricDivisor V});
    X0.cache#"basis indices" = V.cache#"basis indices";
    X0)

pt = base(symbol a, symbol b)
Vas = for V in newVs list abstractVariety(V, pt)
Xas = for X in newXs list abstractVariety(X, pt)
Is = for X in Xas list intersectionRing X
for V in newVs list classGroup V -- none have torsion.
-- we need bases for each of these
for V in newVs list transpose matrix degrees ring V
for X in newXs list X.cache#"basis indices"
for X in newXs list hodgeDiamond X -- #10,11,13 seem to have h11=2, h12=86.  #30 seems to have h12 = 128.

for i in toList(0..#newXs-1) list (
    X := newXs#i;
    IX := Is#i;
    A := coefficientRing IX;
    basIndices := X.cache#"basis indices";
    c2 := chern_2 tangentBundle Xas#i;
    h := sum for j from 0 to #basIndices - 1 list A_j * IX_(basIndices#j);
    print h;
    L := sub(integral(h*c2), R);
    F := sub(integral(h^3), R);
    (L, F)
    )
new128 = last oo
new128 = {2, 128, new128#0, new128#1}
tops128 = append(tops128, new128)
netList tops128
for i from 0 to #tops128 - 1 list (
    x := tops128#i;
    invariantsAllX(x#2, x#3, x#0, x#1) => i
    )
partition(x -> x#0, oo)
for k in keys oo list k => (oo#k/(x -> x#1))
tops128/(x -> invariantsAllX(x#2, x#3, x#0, x#1))

-- the two that may be the same: (29,0), new128
netList tops128_{2,5}
+-+---+---------+-----------+
| |   |         |  2        |
|2|128|36a + 24b|6a b       |
+-+---+---------+-----------+
| |   |         |    2     3|
|2|128|24a + 12b|6a*b  - 6b |
+-+---+---------+-----------+

-- check equivalence of these.
  RQ = QQ (monoid R)
  (A, phi) = genericLinearMap RQ  
  T = target phi
  (L1, F1) = (sub(tops128_2_2, T), sub(tops128_2_3, T))
  (L2, F2) = (sub(tops128_5_2, T), sub(tops128_5_3, T))
  I1 = ideal last coefficients ((phi L1) - L2)
  I2 = ideal last coefficients ((phi F1) - F2)
  I = I1 + I2
  trim I
  A0 = A % sub(ideal gens gb I, ring A)
  phi0 = map(R, R, lift(A0, ZZ))

  (L1, F1) = (tops128_2_2, tops128_2_3)
  (L2, F2) = (tops128_5_2, tops128_5_3)
  
  phi0 L1 - L2
  phi0 F1 - F2

  
--  let's look at morio cones of these two.
-- First, what curves on V lie on X: F=0?

---------------------
-- Analyze Mori cone/nef cone and corresponding maps to P^N's
---------------------
-- First we do this for the "new" #30. XXX
V = newVs_3
X = newXs_3
FX = first equations X
B = ideal V
S = ring V
for c in orbits(V, 1) list c => saturate(ideal FX + ideal for c1 in c list S_c1, B)
-- one toric curve appears to lie on X: x2 = x3 = x4 = 0.

Xa = abstractVariety(X, pt)
IX = intersectionRing Xa
X.cache#"basis indices" -- {0,2}
transpose  matrix degrees ring V
effV = posHull transpose  matrix degrees ring V
rays effV
use coefficientRing IX
e = a*IX_0 + b*IX_2
integral(e * IX_4^2)
integral(e^3)
F = sub(oo, R)
factor F

-- Let's compute the image of FX=0 -- under D4 (of degree {-1,0}
-- phi : X --> PP^5
hh^0(OO_X(-1,0)) -- 6 sections
-- is this map a morphism?
trim(ideal basis(degree V_4, ring V) + ideal FX)
saturate(oo, B) == 1 -- yes, the map is base point free.
assert isNef V_4 -- this is nef on V, therefore on X too.
-- There is some chance this is ample on X, is that possible?
integral(IX_4^2 * (-a * IX_4 + b * IX_5))
integral(IX_4^2 * (-8 * IX_4 + 3 * IX_5)) -- degree {8,3}
-- is this surface effective?
HH^0(V, OO_V(8,3)) -- nothing...
hh^0(OO_X(8,3)) -- nothing...
HH^0(V, OO_V(-8,-3)) -- nothing...
hh^0(OO_X(-8,-3)) -- nothing...
-- It appears it is not effective.  Interesting, as last night I thought it was...
-- How about a curve whose image is a point?
-- For C to be a curve class which gets contracted under phi, need V_4 * C == 0.
-- i.e.
C = a*t_4*t_5 + b*t_5^2
integral(IX_4 * IX_2 * IX_3) == 0 -- so this curve (x2=x3=FX=0) should be contracted to a point?
integral(IX_5 * IX_2 * IX_3) == -2
IX_2*IX_3
use S
trim ideal(FX, x_2, x_3) -- looks like a double line?
IX_4 * IX_2
IX_4 * IX_3
-- actually, it looks like C = (x_2 = x_3 = x_5 = 0) lies on $X$ (but it looks double on $X$?)
-- C is likely a P^1.

-- if this is true, nef cone for X is generated by degrees (-1,0),(-3,-1).
-- and the Mori cone is generated by (0,-2), (1,0)
-- for (1,0) or (-1,0), what is the curve?
-- curve class
use coefficientRing IX
for x in subsets(splice{0..5}, 2) list (
    x => {- integral(IX_(x#0) * IX_(x#1) * IX_4), integral(IX_(x#0) * IX_(x#1) * IX_5)}
    )
C = a * t_4*t_5 + b * t_5^2
integral(IX_4 * C) - 1
integral(IX_5 * C)

-- What about the local equations and the morphism itself on these
-- local equations
D = V_4
monomials D
monomialsToLatticePoints(D, monomials D)
QQCartierCoefficients D
debug NormalToricVarieties
cartierCoefficients(3*V_4)

RAYS = rays V
mysigma = last max V
hilbB = hilbertBasis dualCone posHull transpose matrix RAYS_mysigma
transpose matrix o80
D = 3*V_4
msigma = flatten entries last cartierCoefficients(D)
for p in monomialsToLatticePoints(D, monomials D) list p + msigma

t = symbol t
s = symbol s
R = ZZ/32003[t_0..t_3, s_0..s_3, x_0..x_3, MonomialOrder => {8,4}]
R1 = R/(ideal for i from 0 to 3 list s_i * t_i - 1)
I1 = ideal for j from 0 to #hilbB - 1 list (
    e1 := flatten entries hilbB_j;
    x_j - product for i from 0 to #e1-1 list if e1#i < 0 then s_i else if e1#i > 0 then t_i else 1_R1
    )

see I1
DT = -toricDivisor V
msigma = flatten entries last cartierCoefficients DT
for p in monomialsToLatticePoints(DT, monomials DT) list (
    e := p + msigma;
    m := product for i from 0 to #e-1 list if e#i < 0 then s_i^(-e#i) else if e#i > 0 then t_i^(e#i) else 1_R1;
    {e, m % I1})
    

for 

for x in max V list dualCone posHull transpose matrix RAYS_x
for x in oo list hilbertBasis x
last oo
PD = polytope(-toricDivisor V)

# max V
vertices PD
-- IX_4

degree V_0
isNef V_0

-- This is the potential vertex for dual to cone sigma.
-- If QQ-basepointfree, then these are all vertices, but may have some duplication.
-- If no duplication, then D is ample.
msigma = (V, sigmaIndices, D) -> (
    RAYS := transpose matrix rays V;
    RAYS = RAYS ** QQ;
    RHS := transpose matrix {-(entries vector D)_sigmaIndices};
    RHS = RHS ** QQ;
    MYRAYS := RAYS_sigmaIndices;
    msig := RHS // MYRAYS;
    flatten entries msig
    )

isAmple (V_4 + V_5)
isNef V_4
for x in max V list msigma(V, x, V_4)
for x in max V list msigma(V, x, V_0)
polytope V_0

msigma(V, {1,3,4,5}, V_4)
vertices polytope V_4
for x in max V list msigma(V, x, V_4)
msigma(V, {1,3,4,5}, V_4)
msigma(V, {1,3,4,5}, V_0)
polytope V_0
basis(degree V_0, S)
V_0
u0 = transpose matrix{(rays V)#0}

sigmaIndices = last max V
sigma = posHull transpose matrix RAYS_sigmaIndices
sigmaDual = dualCone sigma
rays sigmaDual
hilbertBasis sigmaDual
-- need the corresponding vertex in M.
transpose rays sigmaDual * transpose matrix RAYS
vertices PD
-- For #30, which is a flop from this, most likely.

P5xV = (coefficientRing S)[x_0..x_5, y_0..y_5, Degrees => {1,1,1,2,3,8,1,1,1,1,1,1}, MonomialOrder => Eliminate 6]
B = sub(basis(degree V_4, ring V), P5xV)
I = ideal(sub(FX, P5xV)) + ideal(B - matrix{{y_0..y_5}})
gens gb I;
I = ideal(B - matrix{{y_0..y_5}})
eliminate(I, splice{x_0..x_5})


------------------------
-- Analyze (29,0) in a similar manner
------------------------

X' = Xs#(29,0)
V' = normalToricVariety(X', CoefficientRing => ZZ/32003)
S' = ring V'
F' = random(degree(-toricDivisor V'), S')
B' = ideal V'
for c in orbits(V', 1) list c => saturate(ideal F' + ideal for c1 in c list S'_c1, B')
  -- three toric curves lie on this X.

classifyExtremalCurves X'
classifyExtremalCurves Xs#(30,0)
transpose matrix degrees ring V'
cubicForm X
cubicForm X'
cubicForm Xs#(29,0)
c2Form X'
c2Form Xs#(29,0)

netList for lab in labs128 list lab => classifyExtremalCurves Xs#lab

-- upshot for h12=128:
-- (27,0), (28,0) are Wall equivalent
-- (29,0), (new) are Wall equivalent (assuming no torsion in (new).
-- (30,0) is alone
-- (31,0) is not Wall equivalent to any of these.

-- Let's look at the effective cone of each or (29,0) and (new).
-- (29,0):
--  Mori cone appears to be {0,1}, {1,0}, so nef cone is positive quadrant.
-- (new)?

describe X
cohomologyVector OO_X(1,0)
cohomologyVector OO_X(0,1)
for a  from -1 to 1 list for b from -1 to 1 list cohomologyVector OO_X(a,b)
hh^0 OO_X(1,0) == 0
hh^0 OO_X(1,1) == 12
hh^0 OO_X(1,2)
hh^0 OO_X(2,1)

for a  from -1 to 1 list for b from -1 to 1 list cohomologyVector OO_X'(a,b)
netList oo

hh^0 OO_X'(1,0) == 3
hh^0 OO_X'(0,1) == 2
hh^0 OO_X'(0,2) == 3
hh^0 OO_X'(0,3) == 4
hh^0 OO_X(1,2)
hh^0 OO_X(2,1)


-- Let's compute the nef cone and types of walls of this cone, for X = Xs#(30,0)
X = Xs#(30,0)
V = normalToricVariety(X, CoefficientRing => ZZ/32003)
cubicForm X
Xa = abstractVariety(X, pt)
basisIndices X -- {0,2}
IX = intersectionRing Xa
D0 = IX_0
D2 = IX_2
integral(D0^3)
integral(D0^2*D2)
integral(D0*D2^2)
integral(D2^3)
-- possible nef cone walls: {2,1}, {1,1}
classifyExtremalCurves X
hh^0(OO_X(2,1)) -- 6
hh^0(OO_X(1,1)) -- 3
hh^0(OO_X(1,2))
isNef (2*V_0 + V_2)
isNef (2*(V_0 + V_2))
not isNef (V_0 + 2 * V_2) -- not nef on V, but on X?
HH^0(V, OO_V(2,1))
HH^0(V, OO_V(1,1))
HH^0(V, OO_V(2,2))
ideal basis({2,2}, ring V) + ideal random(degree(-toricDivisor V), ring V)
decompose oo
decompose ideal V
use coefficientRing IX
isNef(V_0 + V_2)
5 * (a*IX_0 + b*IX_1) * IX_0 + 4 * (a*IX_0 + b*IX_1) * IX_2
(IX_0 + IX_2) *  (IX_2) -- so V_2 intersect X maps to a point here, so is contracted.
hh^0(OO_X(2,2))

----------------------------------
-- Xs#(29,0)
X = Xs#(29,0)
V = normalToricVariety(X, CoefficientRing => ZZ/32003)
Xa = abstractVariety(X, pt)
IX = intersectionRing Xa
basisIndices X -- {0,1}
transpose  matrix degrees ring V
h = IX_0
use coefficientRing IX
e = a*IX_0 + b*IX_1
integral(h^2*e) -- b = 0.
for i from 0 to numgens IX - 1 list h * IX_i
for i from 0 to numgens IX - 1 list IX_1 * IX_i
-- D = D_1
h = IX_1
h^2*e -- 0...
HH^0(V, OO_V(0,3))
IX_1^2

-- XXXXX grab from here (30 new)

-- Let's analyze phi_(D0) : X --> PP^1
-- D0 is on (D3 = 0)
hh^0 OO_X(-3, -1)
hh^0 OO_X(-6, -2)
isNef(V_0) -- true
Iphi = trim(ideal basis({-3,-1}, S) + ideal FX)
trim(ideal basis(2 * degree V_0, ring V) + ideal FX)

-- Let's analyze phi_(D4) : X --> PP^5
-- Is this base point free?
isNef V_4
degree V_4
hh^0 OO_X(-1,0) -- 6
Iphi = trim(ideal basis({-1,0}, S) + ideal FX)
saturate(Iphi, ideal V) -- so map is defined everywhere.

-- my guess is this contracts no surface, but does contract a curve.
-- I would like to understand the image of X under this map too.

-- What curves are blown down by this map?
for x in subsets(splice{0..5}, 2) list (
    x => integral(IX_(x#0) * IX_(x#1) * IX_4)
    )
-- the curve x2=x3=FX=0 seems to get blown down.
-- looking at FX, and the basis, we see that this curve is blown down.

P5xV = (coefficientRing S)[x_0..x_5, y_0..y_5, Degrees => {1,1,1,2,3,8, 3,3,3,3,3,3}, MonomialOrder => {6,6}]
B = matrix entries sub(basis(degree V_4, ring V), P5xV)
isHomogeneous B
FX' = sub(FX, P5xV)
isHomogeneous FX'
I = ideal(FX') + ideal(B - matrix{{y_0..y_5}})
isHomogeneous I
gens gb I;
gbI = ideal oo
leadTerm(1, gens gbI)
I = ideal(B - matrix{{y_0..y_5}})
eliminate(I, splice{x_0..x_5})

-- take2: do the projective map
P5xV = (coefficientRing S)[x_0..x_5, y_0..y_5, Degrees => {1,1,1,2,3,8, 1,1,1,1,1,1}, MonomialOrder => {6,6}]
B = matrix entries sub(basis(degree V_4, ring V), P5xV)
isHomogeneous B
IV = sub(ideal V, P5xV)
FX' = sub(FX, P5xV)
isHomogeneous FX'
I = ideal(FX') + minors(2, B || matrix{{y_0..y_5}})
saturate(I, IV)
saturate(oo, ideal(y_0..y_5))
Jgr = trim oo
-- We should check: any surface blown down?
-- it appears that if a surface is blown down, it must be numerically on ray (1,-1).
e = a*IX_0 + b*IX_2
e = -a*IX_4 + b*IX_5
e * IX_4^2
IX_4^2 * (8*IX_4 - 3*IX_5) == 0
hh^0(OO_X(-8,3)) -- definitely something here (dim = 4129)...
-- this seems to contract the surface (-8, 3) to a curve?
-- i.e. this is not a flip !?

degree V_4
use ring FX
decompose trim ideal(x_2, x_3, FX)

-- Let's compute the image and also the locus where the map contracts
-- surfaces or curves, or is not birational
PD = polytope(-toricDivisor V)
vertices PD
latticePoints PD -- 172!!
size FX == 172
matrix rays V
LLL syz transpose oo


phi = map(S/FX, P5, basis(degree V_4, ring V))
ID4 = ker phi;
isNef V_0 -- true
for i from 0 to 5 list isNef V_i
hh^0 OO_X(-1,0)
HH^0(V, OO_V(-3,-1)) -- 2
HH^0(V, OO_V(3,1)) == 0
HH^0(V, OO_V(2,1)) -- 3
HH^0(V, OO_V(-2,-1)) == 0
isNef (-3*V_0 - V_2) -- false
isNef (2*V_0 + V_2) -- true
h = IX_0
use coefficientRing IX
e = a*IX_0 + b*IX_1
integral(h^2*e) -- b = 0.

Xs

labs128 = select(sort keys Xs, lab -> hh^(1,2) Xs#lab == 128)
for lab in labs128 list {lab, classifyExtremalCurves Xs#lab}

labs86 = select(sort keys Xs, lab -> hh^(1,2) Xs#lab == 86)
for lab in labs86 list {lab, classifyExtremalCurves Xs#lab}
netList oo

for lab in labs86 list {lab, invariantsAll Xs#lab}
partition(x -> x#1, oo)
(keys oo)/(k -> k => oo#k/first)

-- h12=86
-- 10, 11, 13, flop to a exotic fan CY3
-- 12 flop that doesn't appear as exotic fan CY3
-- nontoric(11) wall equiv to 8 (not Mori equiv).  Original Gross pair.
-- nontoric(13) wall equiv to 12 (not Mori equiv). (So Gross pair).
-- 6+7 wall+GV equivalent
-- all other torics are distinct

-- let's get all c2, cubics for torics
use R
torics = for lab in labs86 list lab => {c2Form Xs#lab, cubicForm Xs#lab, 2, 86}
nontorics = {(10,1) => {38*a + 24*b, -1*a*a*a + 12*a*a*b + 0*a*b*b + 0*b*b*b, 2, 86},
 (11,1) => {24*a + 68*b, 0*a*a*a + 0*a*a*b + 12*a*b*b + 14*b*b*b, 2, 86},
 (12,1) => {50*a + 56*b, 5*a*a*a + 12*a*a*b + 0*a*b*b + -16*b*b*b, 2, 86},
 (13,1) => {24*a + 74*b, 0*a*a*a + 0*a*a*b + 12*a*b*b + 17*b*b*b, 2, 86}
 }
all86 = join(torics, nontorics)
hashTable oo
partition(x -> invariantsAll toSequence x#1, all86)
all86 = hashTable all86

  (A, phi) = genericLinearMap RQ  
  T = target phi

  Ls = {(9, 0), (10, 0)}
  (X1, X2) = toSequence (Ls/(lab -> Xs#lab))
  (L1, F1) = (c2Form X1, cubicForm X1)
  (L2, F2) = (c2Form X2, cubicForm X2)
  -- check X1,X2: SAME
  JL = ideal last coefficients(phi sub(L1, T) - sub(L2, T))
  JF = ideal last coefficients(phi sub(F1, T) - sub(F2, T))
  J = trim(JL + JF)
  trim J
  decompose oo -- not equivalent, no integral points.

  (A, phi) = genericLinearMap RQ  
  T = target phi

  -- 10', 12
  Ls = {(10,1), (12,0)}
  (L1, F1) = toSequence take(all86#(10,1), 2)
  (L2, F2) = toSequence take(all86#(12,0), 2)
  factor det hessian F1
  factor det hessian F2
  JL = ideal last coefficients(phi sub(L1, T) - sub(L2, T))
  JF = ideal last coefficients(phi sub(F1, T) - sub(F2, T))
  decompose trim(JL + JF)
  A0 = sub(A % sub(first oo, ring A), QQ)
  A0 == matrix {{1, 0}, {-1/2, 1}}
  raysOf12 = {{-1, 0, 0, -1}, {-1, 0, 0, 0}, {-1, 0, 1, -1}, {-1, 1, 0, -1}, {0, 0, 0, -1}, {3, -1, -1, 3}}
  basisOf12 = {0, 1}

  Ls = {(10,1), (13,1)}
  (L1, F1) = toSequence take(all86#(10,1), 2)
  (L2, F2) = toSequence take(all86#(13,1), 2)
  factor det hessian F1
  factor det hessian F2
  JL = ideal last coefficients(phi sub(L1, T) - sub(L2, T))
  JF = ideal last coefficients(phi sub(F1, T) - sub(F2, T))
  decompose trim(JL + JF)

  decompose J -- no integral points.
  -- somewhat easier:
  JL = ideal last coefficients(phi sub(L1, T) - sub(L2, T))
  JF = ideal last coefficients(phi sub(F1, T) - sub(F2, T))
  J = trim(JL + JF)
  decompose J

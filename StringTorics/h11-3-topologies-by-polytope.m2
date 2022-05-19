restart
load "WriteToricData.m2"

kk = ZZ/32003
topes = kreuzerSkarke(3, Limit => 10000); -- 244
topes = kreuzerSkarke(10, Limit => 10000, Access => "wget"); -- 
#topes
assert(#topes == 244)
  -*
    elapsedTime Vs = topes / (P -> elapsedTime reflexiveToSimplicialToricVariety(convexHull matrix P, CoefficientRing => kk));
    torsionFrees = positions(Vs, V -> classGroup V == ZZ^3);
    nonTorsionFrees = positions(Vs, V -> classGroup V != ZZ^3);
  *-
  nonTorsionFrees = {0, 9, 10, 55, 62, 232}
  torsionFrees = sort toList(set(0..#topes-1) - set nonTorsionFrees);
  assert(#torsionFrees == 238)

topes = topes_{100, 200, 300}
    torsionFrees = positions(Vs, V -> classGroup V == ZZ^10);
    torsionFrees = {0,1,2}
    Vs/classGroup
debugLevel = 1
elapsedTime Ts = hashTable for i in torsionFrees list i => (print i; elapsedTime (
    P = reflexivePolytopeData matrix topes_i;
    findAllFRSTs P
    ));

///
A = matrix {{-1, -1, -1, -1, 1, 2, -1}, {-1, -1, -1, 3, -1, -1, 1}, {1, 1, 2, -2, 0, 0, 0}, {0, 1, 1, -1, 0, 0, 0}}
A1 = A | map(target A, (ring source A)^1, 0);
Ts = allTriangulations(A1, Fine => true, RegularOnly => true)
///

  -- Naomi gets: 518 total number of triangulations for these 238 polytopes.
  sum for k in sort keys Ts list #Ts#k
RZ = ZZ[x,y,z]
elapsedTime tops = hashTable for i in keys Ts list i => (
    ts := Ts#i;
    for X in ts list elapsedTime (
        t := topologicalDataOfCY3(X, RZ);
        {t#"c2", t#"cubic intersection form", t#"h11", t#"h21"}
        )
    );

sort keys tops

-- create a list of all of the topologies of all the triangulations, adding in which polytope it comes from
alltops = flatten for i in sort keys tops list (
    count := -1;
    for t in tops#i list (count = count+1; t | {i, count})
    )

byInvariants = partition(invariants, alltops)
netList for k in sort keys byInvariants list {k, netList byInvariants#k}

-- Now lets take each of the values of byInvariants, and see which topologies are equivalent.
#sort keys byInvariants -- 16 sets
ks = sort keys byInvariants

RQ = QQ[x,y,z]
(A, xyz) = makeGLRing RQ
inc = map(ring A, RZ, xyz)

positions(ks, k -> #byInvariants#k > 1) == {0, 1, 2, 4, 5, 6, 8, 9, 11, 13}




for k in sort keys byInvariants list (
    elapsedTime separateTopologies byInvariants#k
    )
sum for a in o19 list #keys a
netList oo





elapsedTime Ts = hashTable for i in torsionFrees list i => (print i; elapsedTime (
    h21OfCY reflexivePolytopeData matrix topes_i
    ));

tally for i in keys Ts list #Ts#i

-- Tally{1 => 88}
--       2 => 93
--       3 => 22
--       4 => 9
--       5 => 15
--       6 => 11

-- Question: For those that have more than one triangulation (150 of these polytopes),
-- which have more than one topology?

-- Compute the {c2, D^3} forms for each polytope.
-- Compute the topologies of all triangulations in Ts
RZ = ZZ[x,y,z]
elapsedTime tops = hashTable for i in keys Ts list i => (
    ts := Ts#i;
    for X in ts list elapsedTime (
        t := topologicalDataOfCY3(X, RZ);
        {t#"c2", t#"cubic intersection form", t#"h11", t#"h21"}
        )
    ); -- (w/o h11, h21): 443 seconds...  (w h11, h21): 590 sec.

hashTable((keys tops)/(i -> i => netList tops#i))

alltops = flatten values tops;
#alltops == 517
# unique alltops == 286
partition(test2, alltops)
(pairs oo)/((k,v) -> k => #v)

uniqueTops = hashTable for k in keys tops list k => unique tops#k
hashTable select(pairs uniqueTops, (k,v) -> #v > 1)
uniqueMultiples = hashTable select(pairs uniqueTops, (k,v) -> #v > 1)
hashTable((pairs uniqueMultiples)/((k,v) -> k => netList v))
hashTable((pairs uniqueTops)/((k,v) -> 
        
test1 = (f) -> {(trim content f_0)_0, (trim content f_1)_0}

RQ = QQ[gens RZ]
test2 = (f) -> (
    facs := select((factors f_1 )/toList/last, g -> support g != {});
    l1 := sub(f_0, RQ);
    f1 := sub(f_1, RQ);
    d := dim saturate ideal jacobian f1;
    nc := # decompose ideal(l1, f1);
    singZ := flatten entries gens gb saturate(ideal(f_1) + ideal jacobian f_1);
    badp := select(singZ, a -> support leadTerm a === {});
    badp = if badp === {} then 0 else first badp;
    {badp, (trim content f_0)_0, (trim content f_1)_0, #facs, d, nc, f_2, f_3}
    )

tops1 = hashTable for k in keys uniqueMultiples list k => (
    partition(test1, uniqueMultiples#k)
    )
tops2 = hashTable for k in keys uniqueMultiples list k => (
    partition(test2, uniqueMultiples#k)
    )


tops2A = hashTable for k in keys tops list k => (
    partition(test2, tops#k)
    );

(pairs tops2A)/((k,v) -> keys v)//flatten//tally
tops2 = hashTable for k in keys uniqueMultiples list k => (
    p := partition(test1, uniqueMultiples#k);
    if any(values p, v -> #v != 1) then p else continue
    )
#keys tops2 -- 37 left here...
tops2A
flatten values tops2A
keys tops2A

tops
topsB = (pairs tops)/((k,v) -> for v1 in v list {v1, k})//flatten
topsB#0#0
topsC = partition(t -> {t#0#0, t#0#1}, topsB)
# keys topsC == 286
(topsC#(first keys topsC))/(v -> {last first v, last v})
topsD = (pairs topsC)/((k,v) -> (k => (v/(v1 -> {last first v1, last v1}))))//hashTable

(values topsD)/unique

(pairs topsC)((k,v) -> print v)
(k => {last first v, last v}))


H = partition(test2, flatten values tops);
H = (pairs H)/((k,v) -> k => unique v)//hashTable

uniquesInH = (select(pairs H, (k,v) -> #v == 1))//hashTable
  #keys uniquesInH == 102
-- The ones to determine if they give different topologies  
H1 = (select(pairs H, (k,v) -> #v > 1))//hashTable
netList for k in sort keys H1 list k => #H1#k  


-- some examples
-- example 1
-- upshot: there are 3 (L,F)'s here, from 2 different polytopes.
  -- they ALL have hessian factoring as l1*l2^2.
  -- using these linear forms, we can find psi that maps between 1st, 2nd,  and 
  -- another from the 1st to 3rd.
  -- Therefore: they all have the same topology!
  fs = H1#{0, 2, 1, 1, 1, 1, 3, 69}
  
  (L1, L2, L3) = toSequence(fs/first)
  (F1, F2, F3) = toSequence(fs/(f -> f#1))
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F1, RQ)))
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F2, RQ)))
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F3, RQ)))
  -- all of these Hessians factor as l1 * l2^2 * 6912.
  
  -- can we find linear transformations which maps these to each other?
  select(keys tops, k -> any(tops#k, v -> v == fs#0))
  select(keys tops, k -> any(tops#k, v -> v == fs#1))
  select(keys tops, k -> any(tops#k, v -> v == fs#2))

  -- find phi s.t. phi(26x+8y+42z) = 8x + 16y + 42z
  --               phi(x+z) = z or -z
  --               phi(y-z) = x-y-z or negative of this.
      
                   
  A = QQ[a_1..a_9, x, y, z]
  M = genericMatrix(A, 3, 3)
  lins = flatten entries (matrix{{x,y,z}} * M)
  phi = map(A, A, {a_1..a_9} | lins)
  J = ideal(phi(26*x+8*y+42*z) - (8*x+16*y+42*z),
        phi(x+z) - z,
        phi(y-z) - (x-y-z))
  trim ideal last coefficients(gens J, Variables => {x,y,z})
  soln = matrix{{x,y,z}} * (M % oo)
  -- {-y, x, y+z})
  phi.matrix % ooo
  
  use RZ
  psi = map(RZ, RZ, sub(soln, RZ))
  psi L1 - L2 == 0
  psi F1 - F2 == 0

  use A
  J = ideal(phi(26*x+8*y+42*z) - (16*x+8*y+42*z),
        phi(x+z) - z,
        phi(y-z) + (x-y+z))
  trim ideal last coefficients(gens J, Variables => {x,y,z})
  soln = matrix{{x,y,z}} * (M % oo)
  -- {-x, y, x+z})
  phi.matrix % ooo
  
  use RZ
  psi = map(RZ, RZ, sub(soln, RZ))
  psi L1 - L3 == 0
  psi F1 - F3 == 0


-- example 2:
  fs = H1#{0, 2, 1, 1, 1, 1, 3, 75}
  (L1, L2, L3) = toSequence(fs/first)
  (F1, F2, F3) = toSequence(fs/(f -> f#1))
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F1, RQ)))
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F2, RQ)))
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F3, RQ)))

  -- TODO: looks to be the same story as example 1.

-- example 3: 13 examples have this data
  fs = H1#{0, 2, 1, 1, 1, 1, 3, 99}
  (L1, L2, L3) = take(toSequence(fs/first), 3)
  (F1, F2, F3) = take(toSequence(fs/(f -> f#1)), 3)
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F1, RQ)))
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F2, RQ)))
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F3, RQ)))
  netList for v in fs list   factor det diff(transpose vars RQ, diff(vars RQ, sub(v_1, RQ)))

-- example 4:
  fs = H1#{0, 2, 1, 1, 1, 1, 3, 111}
  (L1, L2, L3) = take(toSequence(fs/first), 3)
  (F1, F2, F3) = take(toSequence(fs/(f -> f#1)), 3)
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F1, RQ)))
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F2, RQ)))
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F3, RQ)))
  netList for v in fs list   factor det diff(transpose vars RQ, diff(vars RQ, sub(v_1, RQ)))

-- example 5 (Hessian NOT l1*l2^2)
  fs = H1#{0, 2, 1, 2, 1, 2, 3, 78}
  (L1, L2, L3) = take(toSequence(fs/first), 3)
  (F1, F2, F3) = take(toSequence(fs/(f -> f#1)), 3)
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F1, RQ)))
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F2, RQ)))
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F3, RQ)))

-- example 6: Hessians's don't factor here! But: there is a singular point on the g=1 curve
  fs = H1#{0, 4, 1, 1, 1, 1, 3, 103}
  (L1, L2, L3) = take(toSequence(fs/first), 3)
  (F1, F2, F3) = take(toSequence(fs/(f -> f#1)), 3)
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F1, RQ)))
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F2, RQ)))
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F3, RQ)))

-- example 7: Hessians's don't factor here! But: there is a singular point on the g=1 curve
  fs = H1#{0, 4, 1, 1, 1, 1, 3, 111}
  (L1, L2) = take(toSequence(fs/first), 2)
  (F1, F2) = take(toSequence(fs/(f -> f#1)), 2)
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F1, RQ)))
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F2, RQ)))
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F3, RQ)))
  L1
  L2
  decompose sub(ideal(L1, F1), RQ)
  gens gb ideal(L1, F1)
  gens gb ideal(L2, F2)
  decompose sub(ideal F1 + ideal jacobian F1, RQ) -- (x,y)
  decompose sub(ideal F2 + ideal jacobian F2, RQ) -- (y,z)

-- example 8. 2 examples. there are 2 points to match (this then matches the line oo). NOT DONE
  k = (sort keys H1)#7
  k == {0, 4, 1, 2, 1, 3, 3, 75}
  fs = H1#k
  (L1, L2) = take(toSequence(fs/first), 2)
  (F1, F2) = take(toSequence(fs/(f -> f#1)), 2)
  decompose sub(ideal(L1, F1), RQ)
  decompose sub(ideal(L2, F2), RQ)

-- example 9.  7 different (L,F)'s.
  k = (sort keys H1)#8
  k == {0, 4, 2, 1, 1, 1, 3, 99}
  fs = H1#k
  netList oo
  (L1, L2, L3, L4, L5, L6, L7) = toSequence(fs/first)
  (F1, F2, F3, F4, F5, F6, F7) = toSequence(fs/(f -> f#1))
  G1 = det diff(transpose vars RQ, diff(vars RQ, sub(F1, RQ)))
  G1 = G1//2160
  G2 = det diff(transpose vars RQ, diff(vars RQ, sub(F2, RQ)))
  G2 = G2/2160
  G3 = det diff(transpose vars RQ, diff(vars RQ, sub(F3, RQ)))
  G3 = G3 // 2160
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F4, RQ)))
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F5, RQ)))
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F6, RQ)))
  factor det diff(transpose vars RQ, diff(vars RQ, sub(F7, RQ)))

  -- determined this from the singular locus of the Hessians of F1, F2.
  -- they map (x,z) to (x,z).
  phi = map(RZ, RZ, {x, y-z, z})
  phi L1 - L2
  phi F1 == F2
  F1-F2
  factor(F2-F3)
  factor(F3-F4)
  factor(F2-F4)
  factor(F5-F4)

-- {0, 8, 2, 1, 1, 1, 3, 69}
  fs = H1#((sort keys H1)#20)
  netList fs
  k = {0, 8, 2, 1, 1, 1, 3, 69}
  positions(sort keys H1, k1 -> k1 == k)

  (L1, F1) = (fs#0#0, fs#0#1)
  (L2, F2) = (fs#1#0, fs#1#1)
  (L3, F3) = (fs#2#0, fs#2#1)
  F1
  decompose ideal jacobian sub(F1, RQ)
  decompose ideal jacobian sub(F2, RQ)
  decompose ideal jacobian sub(F3, RQ)
  
/opt/homebrew/bin/points2finetriangs --regular --heights -v < /var/folders/wb/8v4mm0j52pq9pf5gkr8f23z40000gr/T/M2-48780-0/398.in


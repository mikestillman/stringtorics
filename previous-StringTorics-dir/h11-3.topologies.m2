restart
load "WriteToricData.m2"

kk = ZZ/32003
topes = kreuzerSkarke(3, Limit => 10000); -- 244
assert(#topes == 244)
  -*
    elapsedTime Vs = topes / (P -> elapsedTime reflexiveToSimplicialToricVariety(convexHull matrix P, CoefficientRing => kk));
    torsionFrees = positions(Vs, V -> classGroup V == ZZ^3);
    nonTorsionFrees = positions(Vs, V -> classGroup V != ZZ^3);
  *-
  nonTorsionFrees = {0, 9, 10, 55, 62, 232}
  torsionFrees = sort toList(set(0..#topes-1) - set nonTorsionFrees);
  assert(#torsionFrees == 238)

Ts = hashTable for i in torsionFrees list i => (print i; elapsedTime (
    P = reflexivePolytope matrix topes_i;
    findAllFRSTs P
    ))

for i in keys Ts list #Ts#i

{45, 51, 57, 59, 63, 65, 66, 67, 69, 71, 72, 73, 75, 76, 
    77, 78, 79, 81, 83, 84, 85, 87, 89, 91, 93, 95, 99, 
    103, 105, 107, 111, 115, 119, 123, 127, 131, 141, 
    165, 195, 231, 243}

-- Organize all (torsion free) polyhedra via h21 of the CY3)
elapsedTime h21s = partition(i -> (
      P = reflexivePolytope matrix topes_i;
      h21OfCY P
      ),
    torsionFrees) -- 121 seconds.
  assert(sort keys h21s == {45, 51, 57, 59, 63, 65, 66, 67, 69, 71, 72, 73, 75, 76, 
    77, 78, 79, 81, 83, 84, 85, 87, 89, 91, 93, 95, 99, 
    103, 105, 107, 111, 115, 119, 123, 127, 131, 141, 
    165, 195, 231, 243}
    )

-- Compute the topologies of all triangulations in Ts
RZ = ZZ[x,y,z]
elapsedTime tops = hashTable for i in keys Ts list i => (
    ts := Ts#i;
    for X in ts list elapsedTime (
        t := topologicalDataOfCY3(X, RZ);
        {t#"c2", t#"cubic intersection form"}
        )
    ); -- 443 seconds...

topsByH21 = hashTable for h21 in keys h21s list h21 => (
    for i in h21s#h21 list tops#i
    );

hashTable for h21 in keys h21s list h21 => (
    netList for i in h21s#h21 list tops#i
    )

uniqueFormsByH21 = hashTable for h21 in keys h21s list h21 => (
    unique flatten for i in h21s#h21 list tops#i
    )

-- consider all contents for c2, D^3 and #components of C=0
hashTable for h21 in keys h21s list h21 => (
    tally for f in uniqueFormsByH21#h21 list (
        facs := factors f_1;
        facs = for x in facs list (g := last x; if support g == {} then continue else g);
        {(trim content f_0)_0, (trim content f_1)_0, #facs}
        )
    )

-- topologies for all h11=3 triangulations.
-- tops: key is h21
--       value is a hashtable with key: number of the polytope, value: list of (c2, D^3) forms, one for each triangulation.
top3 = hashTable for h21 in keys h21s list h21 => (
    hashTable for i in h21s#h21 list i => (
        -- now we list all of the topological data for this polytope
        P = reflexivePolytope matrix topes_i;
        Ts = findAllFRSTs P;
        << "i = " << i << " triangulations: " << netList Ts << endl;
        topD = Ts/(X -> topologicalDataOfCY3(X, RZ));
        forms = topD/(t -> {t#"c2", t#"cubic intersection form"});
        forms
        );
    )

findTopologies = (h21) -> (
  unique flatten for i in h21s#h21 list (
    P = reflexivePolytope matrix topes_i;
    Ts = findAllFRSTs P;
    << "i = " << i << " triangulations: " << netList Ts << endl;
    topD = Ts/(X -> topologicalDataOfCY3(X, RZ));
    forms = topD/(t -> {t#"c2", t#"cubic intersection form"});
    forms
    ))

findTopologies1 = (h21) -> (
  for i in h21s#h21 list (
    P = reflexivePolytope matrix topes_i;
    Ts = findAllFRSTs P;
    << "i = " << i << " triangulations: " << netList Ts << endl;
    topD = Ts/(X -> topologicalDataOfCY3(X, RZ));
    forms = topD/(t -> {t#"c2", t#"cubic intersection form"});
    forms
    ))


----------------------------------
-- Analyze h11=3, h21=81
-- Number of polytopes: 11.
-- Total # of triangulations: 24
-- Removing uniqueness separately for each polytope: 14 different.
  assert(h21s#81 == {94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104})
  netList topsByH21#81
  netList(topsByH21#81/(fs -> unique fs))
  flatten topsByH21#81
  tops0 = unique flatten topsByH21#81
  for t in tops0 list ((trim content t_0)_0, (trim content t_1)_0)

  byContent = partition(t -> ((trim content t_0)_0, (trim content t_1)_0),
      tops0)
  (keys byContent)/(k -> {k, #byContent#k})

  -- so at least 6 different topologies.
  netList for k in keys byContent list netList for f in byContent#k list {f_0, factors f_1}
  
    netList tops#95
    netList tops#103
    netList (topsByH21#81/unique)

  -- the following is used to find linear transformation mapping c2,D^3 for #94 tri1, #104 tri1.
  -- we can make this into a routine as well...
    A = QQ[a_1..a_9, x, y, z]
    M = genericMatrix(A, 3, 3)
    lins = flatten entries (matrix{{x,y,z}} * M)
    phi = map(A, A, {a_1..a_9} | lins)
    
    L1 = y+z
    L2 = x-2*z
    L3 = x+y+z
    c2 = x+4*y+6*z

    M1 = z
    M2 = x-y-2*z
    M3 = x+y+z
    d2 = x+3*y+6*z
    
    phi matrix{{L1,L2,L3,c2}} - matrix{{M1, M2, M3, d2}}
    trim ideal last coefficients(oo, Variables => {x,y,z})
    M % oo

  -- (6,3)
    unique topsByH21#81#2
    (L1, F1) = toSequence topsByH21#81#0#0
    (L2, F2) = toSequence topsByH21#81#-1#0
    decompose ideal jacobian sub(F1, RQ) -- smooth cubic
    decompose ideal jacobian sub(F2, RQ) -- smooth cubic
    gens gb saturate ideal jacobian F1 -- both give (9)
    gens gb saturate ideal jacobian F2 -- both give (9)
    det diff(transpose vars RQ, diff(vars RQ, sub(F1, RQ)))
    det diff(transpose vars RQ, diff(vars RQ, sub(F2, RQ)))

    A = QQ[a_1..a_9, x, y, z]
    M = genericMatrix(A, 3, 3)
    lins = flatten entries (matrix{{x,y,z}} * M)
    phi = map(A, A, {a_1..a_9} | lins)
    
    L1 = z
    L2 = y+2*z
    L3 = x-z
    c2 = x+7*y+12*z

    M1 = y+z
    M2 = y+2*z
    M3 = x-z
    d2 = x+6*y+12*z
    
    phi matrix{{L1,L2,L3,c2}} - matrix{{M1, M2, M3, d2}}
    vals = trim ideal last coefficients(oo, Variables => {x,y,z})
    M % vals

    phi.matrix % vals
    use RZ
    phi0 = map(RZ, RZ, {x+y, -y, y+z})
    phi0 L1 == L2 -- original L1, L2, F1, F2.
    phi0 F1 == F2

  -- (4,2)
    netList topsByH21#81
    netList(topsByH21#81/(fs -> unique fs))
    flatten topsByH21#81
    (L1, F1) = toSequence (topsByH21#81/(fs -> unique fs))#2#0
    (L2, F2) = toSequence (topsByH21#81/(fs -> unique fs))#4#1
    decompose ideal jacobian sub(F1, RQ)
    decompose ideal jacobian sub(F2, RQ)
    gens gb saturate ideal jacobian F1 
    gens gb saturate ideal jacobian F2 
    det diff(transpose vars RQ, diff(vars RQ, sub(F1, RQ)))
    det diff(transpose vars RQ, diff(vars RQ, sub(F2, RQ)))
    decompose sub(ideal(F1, L1), RQ)
      unique topsByH21#81#2
      tops
    L1 = L1//4
    L2 = L2//4
    L1 = sub(L1, RQ)
    L2 = sub(L2, RQ)
    L1 = sub(L1, A)
    L2 = sub(L2, A)
    trim ideal last coefficients(phi L1 - L2, Variables => {x,y,z})
    -- NOT DONE WITH THIS CASE YET!
    
 -- (2,1)
    (L1, F1) = toSequence byContent#(2,1)#0
    (L2, F2) = toSequence byContent#(2,1)#1
    (L3, F3) = toSequence byContent#(2,1)#2
    decompose ideal jacobian sub(F1, RQ)
    decompose ideal jacobian sub(F2, RQ)
    decompose ideal jacobian sub(F3, RQ)
    det diff(transpose vars RQ, diff(vars RQ, sub(F1, RQ)))
    det diff(transpose vars RQ, diff(vars RQ, sub(F2, RQ)))
    det diff(transpose vars RQ, diff(vars RQ, sub(F3, RQ)))
    gens gb saturate(ideal F1 + ideal jacobian F1)
    gens gb saturate(ideal F2 + ideal jacobian F2)
    gens gb saturate(ideal F3 + ideal jacobian F3)
-- End of analyze h21=81 ---------
----------------------------------

findTopologies 45
netList findTopologies 51
netList findTopologies 57
netList findTopologies 57
netList findTopologies 59
netList findTopologies 63
netList findTopologies 65 -- only 1.
netList findTopologies 66
netList findTopologies 67
tops69 = findTopologies 69
tops169 = findTopologies1 69

tops71 = findTopologies 71 -- 3 different topologies
tops72 = findTopologies 72 -- 
tops73 = findTopologies 73 -- 
#tops69 == 31
netList tops73

for f in tops71 list {(trim content f_0)_0, (trim content f_1)_0}|(factors f_1)
for f in tops73 list {(trim content f_0)_0, (trim content f_1)_0}|{gens gb saturate(ideal jacobian f_1)}|(factors f_1) 
for f in tops72 list decompose sub(ideal jacobian f_1, RQ)
netList oo
for f in tops71 list factors f_1
tops69'21 = for f in tops69 list (
    if {(trim content f_0)_0, (trim content f_1)_0} == {2,1}
    then f else continue)
for f in tops69'21 list (
    selectInSubring(1, gens gb saturate ideal jacobian f_1)
    )
for f in tops72 list (
    gens gb saturate ideal jacobian f_1
    )
tops72/(f -> factor f_1)//netList
RQ = QQ[gens RZ]
for f in tops69'21 list (
    gens gb saturate sub(ideal jacobian f_1, RQ)
    )


netList oo
netList oo
alltopes = hashTable for h21 from 1 to 300 list (
    elapsedTime vals := kreuzerSkarke(3, h21);
    if #vals == 0 then continue;
    h21 => vals
    );

ks = select(keys alltopes, k -> #alltopes#k > 0)
alltopes = hashTable for k in ks list k => alltopes#k;
for k in keys alltopes list #alltopes#k
topes = kreuzerSkarke(3);

for k in sort keys alltopes list (
    P = reflexivePolytope matrix topes_i;
    Ts = findAllFRSTs P;
    << "i = " << i << " triangulations: " << netList Ts << endl;
    topD = Ts/(X -> topologicalDataOfCY3(X, RZ));
    forms = topD/(t -> {t#"c2", t#"cubic intersection form"});
    forms = unique forms;
    if #forms != 2 then continue else i => forms
    )

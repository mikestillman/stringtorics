-- Goal: given h11, or a set of polytopes:
--  construct all triangulations.
--  index each as {polytope#, triangulation#}
--  construct the list of such polytopes.

debug needsPackage "StringTorics"

end--

restart
uninstallAllPackages()

///
-- Create the h11=3 database.
restart
needs "./finding-all-topologies.m2"
topes = kreuzerSkarke(3, Limit => 10000);
assert(#topes == 244)
createCYDatabase("polytopes-h11-3.dbm", topes)
///

///
  -- Create the h11=4 database.
  restart
  needsPackage "StringTorics"
  topes = kreuzerSkarke(4, Limit => 10000);
  assert(#topes == 1197)
  elapsedTime createCYDatabase("polytopes-h11-4.dbm", topes) -- 730 sec
///

///
  restart
  needsPackage "StringTorics"
  topes = kreuzerSkarke(5, Limit => 10000);
  assert(#topes == 4990)
  elapsedTime createCYDatabase("polytopes-h11-5.dbm", topes)
///



///
-- check that it is working: find all the (3,h21) pairs in the table.
restart
needsPackage "StringTorics"
F = openDatabase "polytopes-h11-3.dbm"
F#"1"
Qs = for i from 0 to 243 list cyPolytope F#(toString i);
tally for Q in Qs list hh^(2,1) Q -- 42 different h12's...
///


-- check that it is working: find all the (4,h21) pairs in the table.
restart
needsPackage "StringTorics"
F = openDatabase "polytopes-h11-4.dbm"
Qs = for i from 0 to 1196 list cyPolytope F#(toString i);
tally for Q in Qs list hh^(2,1) Q -- 87 different h12's...
///




topes = kreuzerSkarke(5, Limit => 20);
elapsedTime V = cyPolytope(topes_4, ID => 4)
str = dump V
V = cyPolytope(str, ID => 4)
assert(h11OfCY V == 5)
assert(h21OfCY V == 29)
assert isFavorable V
V1 = cyPolytope(dump V, ID => 4)
assert(V === V1) -- note that the cache's differ.
dump V

makeCY V

Ts = findAllFRSTs V
X = cyData(V, Ts_0, ID => 0)
dump X
cyData(dump X, i -> V)


X#"polytope data"
dump X#"polytope data"
createCYDatabase("test-polytopes-h11-5.dbm", topes)
F = openDatabase "test-polytopes-h11-5.dbm" -- or also open it for writing?
F#"10"
peek cyPolytope F#"10"
close F

F = openDatabase "polytopes-h11-5.dbm" -- or also open it for writing?
F#"10"
close F

F = openDatabase "polytopes-h11-5.dbm" -- or also open it for writing?
F#"3000"
Q = cyPolytope F#"3000"
findAllFRSTs Q
close F

F = openDatabase "polytopes-h11-5.dbm" -- or also open it for writing?
elapsedTime Qs = for i from 0 to 4989 list cyPolytope F#(toString i); -- 8 sec to read them all... now 3.25 secon...
tally for Q in Qs list (h11OfCY Q, h21OfCY Q)
close F
Ts = Qs/(Q -> elapsedTime findAllFRSTs Q);

-- Analyze one set of triangulations
F = openDatabase "polytopes-h11-5.dbm"
elapsedTime Qs = for i from 0 to 4989 list cyPolytope F#(toString i); -- 8 sec to read them all... now 3.25 secon...
close F
Q = first select(Qs, Q -> h21OfCY Q == 20)
Ts = findAllFRSTs Q;
Xs = for i from 0 to #Ts - 1 list cyData(Q, Ts#i, ID => i);
assert(#Xs == 142)

X = Xs_0
RZ = ZZ[a,b,c,d,e]
isFavorable Q
topologicalData (X, RZ)
cubicForm oo
elapsedTime topXs = unique for X in Xs list topologicalData(X, RZ);
H = partition(X -> topologicalData(X, RZ), Xs);
F1 = cubicForm topXs#0
F2 = cubicForm topXs#1
ideal gens gb saturate(ideal jacobian F1, ideal(a,b,c,d,e))
ideal gens gb saturate(ideal jacobian F2, ideal(a,b,c,d,e))
RQ = QQ(monoid RZ)
saturate ideal jacobian sub(F1, RQ) -- smooth
saturate ideal jacobian sub(F2, RQ) -- smooth
minimalBetti inverseSystem sub(F1, RQ)
minimalBetti inverseSystem sub(F2, RQ)
-- use the following one to get everything working
topes = kreuzerSkarke(5, Limit => 30);

--elapsedTime Ps = topes/matrix/reflexivePolytope;
elapsedTime PCYs = topes/matrix/convexHull/polar/cyPolytopeData;
for i from 0 to #PCYs - 1 do (PCYs#i)#"id" = i
elapsedTime for x in PCYs list annotatedFaces x;

for x in PCYs list isFavorable x

PCYs/degrees


elapsedTime for x in oo list value x;
elapsedTime for x in oo list polytope x;

dbmData = elapsedTime for x in PCYs list writeCYPolytope x;

-- create the polytope data base for a specific set of polytopes (e.g. h11=5)
F = openDatabaseOut "temp-h11-5-cy-polytope-data"
F#"0" = dbmData#0
for i from 0 to #dbmData - 1 do F#(toString i) = writeCYPolytope PCYs#i
close F

F = openDatabase "temp-h11-5-cy-polytope-data"
scanKeys(F, print)
keys F
F#"5"
for i in 0..29 list F#(toString i);
readCYPolytope F#"10"
viewHelp Database
close F

F = openDatabaseOut "temp-h11-5-cy-polytope-data"
keys F
for i from 0 to #dbmData - 1 do F#(toString i) = writeCYPolytope PCYs#i
F#"10"
writeCYPolytope PCYs_0
readCYPolytope oo
cyMakePolytope oo
peek oo
cySetGLSM PCYs_0
cySetH11H21 PCYs_0
PCYs_0
writeCYPolytope PCYs_0
cyMakePolytope oo
oo === o8
CY2 = value writeCYPolytope PCYs_0
assert(polytope CY2 == polytope PCYs_0)
value oo
elapsedTime PCYs = topes/matrix/convexHull/polar/latticePointList; -- much slower than latticePoints...
Ps_0
favorables = positions(Ps, isFavorable)
unfavorables = sort toList(set splice {0.. #topes - 1} - favorables)
unfavorables == {1,2}

debug StringTorics
writePolytope Ps#0
-- For now, we want to construct all of the CalabiYau's coming from all
-- triangulations of these.  Each is indexed by (tope#, triangulation#).
-- tope# is in range 0..#topes-1 (unfavorables are ignored.

Xs = elapsedTime flatten for i in favorables list (
    P = Ps#i;
    elapsedTime xs = findAllFRSTs P;
    for j from 0 to #xs-1 list (i,j) => xs#j
    )
Xs = hashTable Xs
Xs#(29,4)

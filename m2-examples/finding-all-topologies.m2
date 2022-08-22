-- Goal: given h11, or a set of polytopes:
--  construct all triangulations.
--  index each as {polytope#, triangulation#}
--  construct the list of such polytopes.

debug needsPackage "StringTorics"


-- Function: write out all polytopes.
--           read in all polytopes.
--           write out all triangulations
--           read in all triangulations.

createPolytopeDatabase = method()
createPolytopeDatabase(String, List) := (dbfilename, topes) -> (
    -- open data base file
    F := openDatabaseOut dbfilename;
    -- F["info"] = "4990 reflexive polytopes of h11=5"
    -- F["topes"] = toString topes;
    -- loop through topes, create CYPolytopeData, populate it, write it to data base.
    elapsedTime for i from 0 to #topes - 1 do elapsedTime (
        << "computing for polytope " << i << endl;
        V := cyPolytopeData(topes#i, ID => i); -- NOT correct!!! gives the dual...
        -- now fill it with data we want
        basisIndices V; -- compute them
        isFavorable V; -- compute h11, h21, favorability.
        annotatedFaces V; -- compute annotated faces
        -- now write it
        F#(toString i) = dump V;
        );
    close F;
    )
end--

restart
uninstallAllPackages()

restart
needs "./finding-all-topologies.m2"
topes = kreuzerSkarke(5, Limit => 10000);
assert(#topes == 4990)
createPolytopeDatabase("polytopes-h11-5.dbm", topes)

topes = kreuzerSkarke(5, Limit => 20);
elapsedTime V = cyPolytopeData(topes_4, ID => 4)
str = dump V
V = cyPolytopeData(str, ID => 4)
assert(h11OfCY V == 5)
assert(h21OfCY V == 29)
assert isFavorable V
V1 = cyPolytopeData(dump V, ID => 4)
assert(V === V1) -- note that the cache's differ.
dump V

makeCY V

Ts = findAllFRSTs V
X = cyData(V, Ts_0, ID => 0)
dump X
cyData(dump X, i -> V)


X#"polytope data"
dump X#"polytope data"
createPolytopeDatabase("test-polytopes-h11-5.dbm", topes)
F = openDatabase "test-polytopes-h11-5.dbm" -- or also open it for writing?
F#"10"
peek cyPolytopeData F#"10"
close F

F = openDatabase "polytopes-h11-5.dbm" -- or also open it for writing?
F#"10"
close F

F = openDatabase "polytopes-h11-5.dbm" -- or also open it for writing?
F#"3000"
Q = readCYPolytopeData F#"3000"
findAllFRSTs Q
close F

F = openDatabase "polytopes-h11-5.dbm" -- or also open it for writing?
elapsedTime Qs = for i from 0 to 4989 list cyPolytopeData F#(toString i); -- 8 sec to read them all... now 3.25 secon...
tally for Q in Qs list (h11OfCY Q, h21OfCY Q)
close F
Ts = Qs/(Q -> elapsedTime findAllFRSTs Q);

-- Analyze one set of triangulations
F = openDatabase "polytopes-h11-5.dbm"
elapsedTime Qs = for i from 0 to 4989 list cyPolytopeData F#(toString i); -- 8 sec to read them all... now 3.25 secon...
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

dbmData = elapsedTime for x in PCYs list writeCYPolytopeData x;

-- create the polytope data base for a specific set of polytopes (e.g. h11=5)
F = openDatabaseOut "temp-h11-5-cy-polytope-data"
F#"0" = dbmData#0
for i from 0 to #dbmData - 1 do F#(toString i) = writeCYPolytopeData PCYs#i
close F

F = openDatabase "temp-h11-5-cy-polytope-data"
scanKeys(F, print)
keys F
F#"5"
for i in 0..29 list F#(toString i);
readCYPolytopeData F#"10"
viewHelp Database
close F

F = openDatabaseOut "temp-h11-5-cy-polytope-data"
keys F
for i from 0 to #dbmData - 1 do F#(toString i) = writeCYPolytopeData PCYs#i
F#"10"
writeCYPolytopeData PCYs_0
readCYPolytopeData oo
cyMakePolytopeData oo
peek oo
cySetGLSM PCYs_0
cySetH11H21 PCYs_0
PCYs_0
writeCYPolytopeData PCYs_0
cyMakePolytopeData oo
oo === o8
CY2 = value writeCYPolytopeData PCYs_0
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

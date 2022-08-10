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
        A := matrix topes#i;
        P := convexHull A;
        P2 := polar P;
        X := cyPolytopeData P2;
        -- now fill it with data we want
        cySetGLSM X;
        cySetH11H21 X;
        annotatedFaces X;
        F#(toString i) = writeCYPolytopeData X;
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

F = openDatabase "polytopes-h11-5.dbm" -- or also open it for writing?
F#"10"
close F

F = openDatabase "polytopes-h11-5.dbm" -- or also open it for writing?
F#"3000"
Q = readCYPolytopeData F#"3000"
findAllFRSTs Q
close F

F = openDatabase "polytopes-h11-5.dbm" -- or also open it for writing?
elapsedTime Qs = for i from 0 to 4989 list readCYPolytopeData F#(toString i); -- 8 sec to read them all...
close F
Ts = Qs/(Q -> elapsedTime findAllFRSTs Q);


F#"3000"

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

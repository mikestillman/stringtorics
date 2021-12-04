-- Currently this is an "add-on" to string torics, not part of the package yet.
debug needsPackage "StringTorics" 
  -- debug: findFirstUnitVectors, findInvertibleSubmatrix.
  -- note, these functions are dirt slow, I think.

PolytopeData = new Type of HashTable
TopologicalDataOfCY3 = new Type of HashTable

polytopeData = method()
polytopeData Matrix := PolytopeData => (A) -> (
    P2 := polar convexHull A;
    LP = select(latticePointList P2, lp -> dim(P2, minimalFace(P2, lp)) <= 2);
    mLP = transpose matrix LP;
    D := transpose syz mLP;
    p := findFirstUnitVectors D;
    q := findInvertibleSubmatrix(D, p);
    if q === null then error ("oops, can't find a good GLSM matrix");
    GLSM := (D_q)^-1 * D;
    T := findAllFRSTs mLP;
    -- the rays of each triangulation should match LP.
    new PolytopeData from {
        "h11" => h21OfCY P2,
        "h21" => h11OfCY P2,
        "rays" => LP,
        "glsm" => entries transpose GLSM,
        "basis indices" => q,
        "triangulations" => T/last
        }
    )

topologicalDataOfCY3 = method()
topologicalDataOfCY3(PolytopeData, List, Ring) := TopologicalDataOfCY3 => (P, triang, RZ) -> (
    V := normalToricVariety(P#"rays", triang);
    data := topologyOfCY3(V, P#"basis indices");
    new TopologicalDataOfCY3 from {
        "h11" => P#"h11",
        "h21" => P#"h21",
        "c2" => sub(data_3, vars RZ),
        "cubic intersection form" => sub(data_2, vars RZ)
        }
    )

writeDataFormat = method()
writeDataFormat String := (filename) -> (
    F := openOutAppend filename;
    F << "-- Each triangulation entry: " << endl;
    F << "  (h11, h12)" << endl;
    F << "  rays of the fan" << endl;
    F << "  GLSM charges" << endl;
    F << "  basis indices for divisors" << endl;
    F << "  max simplices in the triangulation" << endl;
    F << "  linear form (dot with c2(X))" << endl;
    F << "  cubic intersection form, using basis of divisors above" << endl;
    close F;
    )

writeCY3Data = method()
writeCY3Data(String, ZZ, PolytopeData, List) := (filename, idx, P, Ts) -> (
    if #Ts != # P#"triangulations" then 
        error "expected topological data for each triangulation";
    if not all(Ts, t -> instance(t, TopologicalDataOfCY3)) then 
        error "expected a list of top data for each triangulation";
    F := openOutAppend filename;
    F << "-- polytope " << idx << endl;
    for i from 0 to #Ts - 1 list (
        F << "-- triangulation " << i << endl;
        topdata := Ts#i;
        F << "  " << (P#"h11", P#"h21") << " -- (h11, h21)" << endl;
        F << "  " << P#"rays" << endl;
        F << "  " << P#"glsm" << endl;
        F << "  " << P#"basis indices" << endl;
        F << "  " << P#"triangulations"#i << endl;
        F << "  " << toString(topdata#"c2") << endl;
        F << "  " << toString(topdata#"cubic intersection form") << endl;
        );
    close F
    )

computeAndWriteTopologicalData = method()
computeAndWriteTopologicalData(String, Matrix, ZZ) := (filename, A, idx) -> (
    polydata := polytopeData A;
    nvars := polydata#"h11";
    RZ := ZZ[vars(0..nvars-1)];
    elapsedTime topdata = for t in polydata#"triangulations" list
      elapsedTime topologicalDataOfCY3(polydata, t, RZ);
    writeCY3Data(filename, idx, polydata, topdata)
    )

end--

restart
load "WriteToricData.m2"
kk = ZZ/32003
RZ = ZZ[x,y,z,w]
topes = kreuzerSkarke(4, Access => "wget", Limit => 10000); -- 1197
assert(#topes == 1197)

-*
elapsedTime Vs = topes / (P -> elapsedTime reflexiveToSimplicialToricVariety(convexHull matrix P, CoefficientRing => kk));
torsionFrees = positions(Vs, V -> classGroup V == ZZ^4);
nonTorsionFrees = positions(Vs, V -> classGroup V != ZZ^4);
*-
nonTorsionFrees = {0, 3, 4, 5, 12, 15, 796, 800, 803, 1059, 1060, 1064, 1065, 1134, 1135, 1151, 1153, 1155}
torsionFrees = sort toList(set(0..#topes-1) - set nonTorsionFrees);
assert(#torsionFrees == 1179) -- not 1197!! -- so 18 are torsion...

for x in torsionFrees do (
    << "doing " << x << endl;
    A = matrix topes_x;
    computeAndWriteTopologicalData("foo4", A, x)
    )

A = matrix topes_100
writeDataFormat("foo4")
computeAndWriteTopologicalData("foo4", A, 100)

elapsedTime polydata = polytopeData A
elapsedTime topdata = for t in polydata#"triangulations" list
  elapsedTime topologicalDataOfCY3(polydata, t, RZ)
netList topdata
writeCY3Data("foo", 13, polydata, topdata)
elapsedTime toricdata = computeToricData A
topdatas = for t in last toricdata list elapsedTime computeTopologicalData(RZ, toricdata_2, t, toricdata_4)
use RZ


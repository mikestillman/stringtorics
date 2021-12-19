-- Currently this is an "add-on" to string torics, not part of the package yet.
debug needsPackage "StringTorics" 
  -- debug: findFirstUnitVectors, findInvertibleSubmatrix.
  -- note, these functions are dirt slow, I think.

PolytopeData = new Type of HashTable -- contains: data about a reflexive polytope.
  -- contains reflexive polytope data for P2 (i.e. in N lattice).
  -- in particular: rays, glsm (which it computes), basis indices.
  -- in cache: h11, h21, triangulations?
  --           annotated faces, perhaps?

ReflexivePolytopeData = new Type of HashTable -- contains: data about a reflexive polytope.
  -- TODO: better name?  perhaps ReflexivePolytope

ToricHypersurface = new Type of HashTable
  -- contains PolytopeData, and a triangulation.  TODO: better name? perhaps TriangulatedReflexivePolytope?

TopologicalDataOfCY3 = new Type of HashTable
  -- contains h11, h21, c2, cubic intersection form

dim ReflexivePolytopeData := ZZ => P -> # P#"rays"#0 -- dimension of the polytope
rays ReflexivePolytopeData := List => P -> P#"rays"
degrees ReflexivePolytopeData := List => P -> P#"glsm"
h11OfCY ReflexivePolytopeData := ZZ => P -> (
    if not P.cache#?"h11" then P.cache#"h11" = h21OfCY convexHull transpose matrix rays P;
    P.cache#"h11"
    )
h21OfCY ReflexivePolytopeData := ZZ => P -> (
    if not P.cache#?"h11" then P.cache#"h21" = h11OfCY convexHull transpose matrix rays P;
    P.cache#"h21"
    )

reflexivePolytopeData = method()
reflexivePolytopeData Matrix := ReflexivePolytopeData => (A) -> (
    P2 := polar convexHull A;
    LP = select(latticePointList P2, lp -> dim(P2, minimalFace(P2, lp)) <= 2);
    mLP = transpose matrix LP;
    D := transpose syz mLP;
    p := findFirstUnitVectors D;
    q := findInvertibleSubmatrix(D, p);
    if q === null then error ("oops, can't find a good GLSM matrix"); -- hasn't happened yet
    GLSM := (D_q)^-1 * D;
    -- the rays of each triangulation should match LP.
    new ReflexivePolytopeData from {
        symbol cache => new CacheTable,
        "rays" => LP,
        "glsm" => entries transpose GLSM,
        "basis indices" => q
        }
    )

-- TODO: allow options, e.g. limit the number.
-- TODO: Another routine: find a random one?
-- TODO: Another routine: start with one, then do bistellar flip.
toricHypersurface = method()
toricHypersurface(ReflexivePolytopeData, List) := ToricHypersurface => (P, triang) -> (
    -- TODO: do some basic checking:
    --   triang should be a list of lists of indices from the rays of P,
    --   each of length dim P.
    --   BUT: TODO: really want to be able to make sure it is a triangulation.
    new ToricHypersurface from {
        symbol cache => new CacheTable,
        "polytope" => P,
        "triangulation" => triang
        }
    )

-- Input: ReflexivePolytopeData
-- Output: List of ToricHypersurface's.
findAllFRSTs ReflexivePolytopeData := List => P -> (
    T := findAllFRSTs transpose matrix rays P;
    for t in T list toricHypersurface(P, last t) -- t is a pair: list of vertices, list of list of indices
    )

topologicalDataOfCY3 = method()
topologicalDataOfCY3(ToricHypersurface, List, Ring) := TopologicalDataOfCY3 => (P, triang, RZ) -> (
    V := normalToricVariety(P#"rays", triang);
    data := topologyOfCY3(V, P#"basis indices");
    -- this data above computes intersection numbers for all toric divisors. 
    -- So we consider only the ones whose indices are contained in basis indices:
    new TopologicalDataOfCY3 from {
        "h11" => P#"h11",
        "h21" => P#"h21",
        "c2" => sub(data_3, vars RZ),
        "cubic intersection form" => sub(data_2, vars RZ)
        }
    )

normalToricVariety ToricHypersurface := opts -> (X) -> (
    P := X#"polytope";
    T := X#"triangulation";
    GLSM := P#"glsm";
    normalToricVariety(rays P, T, WeilToClass => transpose matrix GLSM) -- TODO: this fails if the class group is torsion! (Fails: later it gives an inscrutable error...)
    )

gvInputFromToric = method(Options => {
    Heft => null,
    DegreeLimit => infinity,
    Precision => 150
    }
)

gvInputFromToric(NormalToricVariety, List, List) := opts -> (V, basisIndices, moriGenerators) -> (
    << "in gvInput" << endl;
    H := hashTable for i from 0 to #basisIndices-1 list basisIndices#i => i;
    intersectionnums := for x in pairs CY3NonzeroMultiplicities V list (
        if isSubset(x#0, basisIndices) then
            append(sort for a in x#0 list H#a, x#1)
        else 
            continue
        );
    str1 := toString moriGenerators;
    str2 := "{}";
    str3 := toString if opts.Heft === null then heft ring V else opts.Heft; -- this might not be correct
    str4 := toString transpose degrees ring V;
    str5 := toString intersectionnums;
    str6 := toString ({
            if opts.DegreeLimit === infinity then -1 else opts.DegreeLimit, 
            opts.Precision,
            0,
            300000
            });
    concatenate between("\n", {str1, str2, str3, str4, toString {}, str5, str6})
    )

------------------------------
-- Original code below here --
------------------------------
polytopeData = method()
polytopeData Matrix := PolytopeData => (A) -> (
    P2 := polar convexHull A;
    LP = select(latticePointList P2, lp -> dim(P2, minimalFace(P2, lp)) <= 2);
    mLP = transpose matrix LP;
    D := transpose syz mLP;
    p := findFirstUnitVectors D;
    q := findInvertibleSubmatrix(D, p);
    if q === null then error ("oops, can't find a good GLSM matrix"); -- hasn't happened yet
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

TEST ///
  topes = kreuzerSkarke 3;
  P = reflexivePolytopeData matrix topes_12
  -- check consistency:
  transpose matrix rays P
  GLSM = transpose matrix degrees P
  assert(dim P == 4)
  assert(GLSM_(P#"basis indices") == 1)

  -- Now create triangulations  
  Xs = findAllFRSTs P
  
  -- Now compute topology
  V = normalToricVariety Xs_0
  elapsedTime CY3NonzeroMultiplicities V
  for x in pairs oo list append(x#0, x#1)
  elapsedTime tripleProductsCY V -- returns monomials => value.  Bit annoying...
  P#"basis indices"
  topologyOfCY3(V, P#"basis indices")

  "foo1" << gvInputFromToric(V, P#"basis indices", {{1,0,0},{0,1,0},{0,0,1}}) << endl << close;
  -- Now create computeGV input, then get GV invariants.
  
///

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

-- By hand:
P = reflexivePolytopeData matrix topes_10
Ts = findAllFRSTs P
netList Ts

rays P


-- Below this still works (hopefully).



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


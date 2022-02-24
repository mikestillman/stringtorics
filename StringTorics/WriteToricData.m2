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
    if not P.cache#?"h21" then P.cache#"h21" = h11OfCY convexHull transpose matrix rays P;
    P.cache#"h21"
    )

reflexivePolytopeData = method()
reflexivePolytopeData(List, List, List) := ReflexivePolytopeData => (latticePoints, GLSM, basisIndices) -> (
    -- What should be checked here to validate the input data?
    new ReflexivePolytopeData from {
        symbol cache => new CacheTable,
        "rays" => latticePoints,
        "glsm" => GLSM,
        "basis indices" => basisIndices
        }
    )

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
    reflexivePolytopeData(LP, entries transpose GLSM, q)
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

normalToricVariety ToricHypersurface := opts -> (X) -> (
    if not X.cache.?NormalToricVariety then X.cache.NormalToricVariety = (
        P := X#"polytope";
        T := X#"triangulation";
        GLSM := transpose matrix P#"glsm";
        normalToricVariety(rays P, T, opts, WeilToClass => matrix GLSM)
        );
    X.cache.NormalToricVariety
    -- TODO: this fails if the class group is torsion! (Fails: later it gives an inscrutable error...)
    )

topologicalDataOfCY3 = method()
topologicalDataOfCY3(ToricHypersurface, Ring) := TopologicalDataOfCY3 => (X, RZ) -> (
    V := normalToricVariety X;
    P := X#"polytope";
    data := topologyOfCY3(V, P#"basis indices");
    -- this data above computes intersection numbers for all toric divisors. 
    -- So we consider only the ones whose indices are contained in basis indices:
    new TopologicalDataOfCY3 from {
        "h11" => h11OfCY P,
        "h21" => h21OfCY P,
        "c2" => sub(data_3, vars RZ),
        "cubic intersection form" => sub(data_2, vars RZ)
        }
    )

moriCone = method()
moriCone NormalToricVariety := Cone => (V) -> (
    IV := intersectionRing (abstractVariety V);
    Cs := matrix for x in orbits(V, 1) list (
        c := product(x, i -> IV_i);
        for d in gens IV list integral(c*d)
        );
    posHull transpose lift(Cs, QQ)
    )

gvInputFromToric = method(Options => {
    Heft => null,
    DegreeLimit => infinity,
    Precision => 150
    }
)

gvInputFromToric(NormalToricVariety, List, List) := opts -> (V, basisIndices, moriGenerators) -> (
    << "in gvInput" << endl;
    X := completeIntersection(V, {-toricDivisor V});
    Xa := abstractVariety(X, base());
    IX := intersectionRing Xa;
    intersectionnums := for x in pairs intersectionNumbers(IX, basisIndices) list append(x#0, x#1);
    -- H := hashTable for i from 0 to #basisIndices-1 list basisIndices#i => i;
    -- intersectionnums := for x in pairs CY3NonzeroMultiplicities V list (
    --     if isSubset(x#0, basisIndices) then
    --         append(sort for a in x#0 list H#a, x#1)
    --     else 
    --         continue
    --     );
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

gvInput = (moriGenerators, heftval, GLSM, intersectionnums, degreelimit, prec) -> (
    -- moriGenerators: list of lists
    -- heftval: list of ints
    -- GLSM: list of list of ints
    -- intersectionnums: list of triples of ints
    -- degreelimit: infinity or positive integer
    -- prec: positive integer
    str1 := toString moriGenerators;
    str3 := toString heftval;
    str4 := toString GLSM;
    str5 := toString intersectionnums;
    str6 := toString ({
            if degreelimit === infinity then -1 else degreelimit, 
            prec,
            0,
            300000
            });
    concatenate between("\n", {str1, toString {}, str3, str4, toString {}, str5, str6})
    )

gvInvariants = method(Options => {
    Mori => null, -- null means: compute rays of the Mori cone of V in 
    Heft => null, -- null means: compute it
    DegreeLimit => infinity,
    Precision => 150,
    FilePrefix => "foo",
    Executable => "~/src/git-from-others/cytools-private/external/gv/computeGV",
    KeepFiles => true
    })

intersectionNumbersOfCY = method()
intersectionNumbersOfCY(NormalToricVariety, List) := (V, basisIndices) -> (
    X := completeIntersection(V, {-toricDivisor V});
    Xa := abstractVariety(X, base());
    IX := intersectionRing Xa;
    intersectionNumbers(IX, basisIndices)
    )

gvInvariants(NormalToricVariety, List) := HashTable => opts -> (V, basisIndices) -> (
    -- Compute intersection numbers for X in V (using this basis)
    -- Compute mori cone (if needed) (?? requires basis too...)
    -- Compute a vector which dots positively with all these generators.
    -- Then write the file
    -- Execute the command
    -- Read the results, and return them
    intersectionnums := for t in intersectionNumbersOfCY(V, basisIndices) list append(t#0, t#1);
    -- X := completeIntersection(V, {-toricDivisor V});
    -- Xa := abstractVariety(X, base());
    -- IX := intersectionRing Xa;
    -- intersectionnums := for x in pairs intersectionNumbers(IX, basisIndices) list append(x#0, x#1);
    -- H := hashTable for i from 0 to #basisIndices-1 list basisIndices#i => i;
    -- intersectionnums := for x in pairs CY3NonzeroMultiplicities V list (
    --     if isSubset(x#0, basisIndices) then
    --         append(sort for a in x#0 list H#a, x#1)
    --     else 
    --         continue
    --     );
    mori := if opts.Mori =!= null then opts.Mori else (
        M := moriCone V;
        GLSM := matrix degrees ring V;
        entries transpose((rays M) // GLSM)
        );
    heft := if opts.Heft =!= null then opts.Heft else (
      sum entries transpose rays dualCone posHull transpose matrix mori
      );
    -- OK, now we have computed everything we need.  Write it to a file
    infile := opts.FilePrefix | "-input";
    outfile := opts.FilePrefix | "-output";
    infile << gvInput(mori, heft, transpose degrees ring V, intersectionnums,
        opts.DegreeLimit, opts.Precision) << close;
    inputLine := opts.Executable | " <" | infile | " >" | outfile;
    print inputLine;
    run inputLine;
    -- Get the output, package as a hash table
    (lines get outfile)/value//hashTable
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

-*
-- Older version
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
*-

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
--
  restart
  needs "WriteToricData.m2"
--
  topes = kreuzerSkarke 3;
  P = reflexivePolytopeData matrix topes_12
  V = reflexiveToSimplicialToricVariety convexHull matrix topes_12
  -- check consistency:
  transpose matrix rays P
  GLSM = transpose matrix degrees P
  assert(dim P == 4)
  assert(GLSM_(P#"basis indices") == 1)

  -- Now create triangulations  
  Xs = findAllFRSTs P

  V = normalToricVariety Xs_0
  gvInvariants(V, {0, 2, 3}, DegreeLimit => 10, Precision => 150)
  gvInvariants(V, {2, 3, 6}, DegreeLimit => 10, Precision => 150)
    
  -- Now compute topology
  elapsedTime CY3NonzeroMultiplicities V
  for x in pairs oo list append(x#0, x#1)
  elapsedTime tripleProductsCY V -- returns monomials => value.  Bit annoying...
  P#"basis indices"
  topologyOfCY3(V, P#"basis indices")

///

end--

TEST ///
  -- example: mirror of (hypersurface in) P(1,1,6,9)

  restart
  needs "WriteToricData.m2"
  A = matrix"1,0,0,0,-1;0,1,0,0,-1;0,0,1,0,-6;0,0,0,1,-9"

  P2 = convexHull A
  vertices  P2
  latticePoints P2
  isReflexive P2

  V = reflexiveToSimplicialToricVariety(polar P2, CoefficientRing => ZZ/101)

  Ts = findAllFRSTs transpose matrix rays V -- only one.
  t = Ts#0#1
  
  GLSM = transpose matrix degrees ring V -- matches Andres'

  -- XXX
  gvInvariants(V, {0,5}, DegreeLimit => 10, Precision => 150)

  -- Used to check the above:
  -- now we need to compute the intersection numbers and Mori cone.
  pt = base(a,b)
  X = completeIntersection(V, {-toricDivisor V})
  Va = abstractVariety(V, pt)
  Xa = abstractVariety(X, pt)
  basisIndices = {0,5}
  use Xa
  -- Intersection numbers of X:
  integral(t_0^3) == 0
  integral(t_0^2*t_5) == 1
  integral(t_0*t_5^2) == -3
  integral(t_5^3) == 9
  
  --
  rays moriCone V
  mori = entries transpose ((rays moriCone V) // (transpose GLSM))
  -- gives [[1,-3], [0,1]]



  -- TODO: write triangulate.
  --   which vertices to take
  --   how many
  --   keeps origin in it too (or should?)

  --reflexiveToSimplicialToricVariety(P2, Lattice => "N") -- this might be nice...

  -- t = triangulate(P2, ....) -- which lattice points to use?
///


TEST ///
  -- example from Naomi (poly_111): Start this over...
  restart
  needs "WriteToricData.m2"

  rys = transpose matrix{{-1,-1,0,0},
         {-1,0,-1,1},
         {0,-1,-1,0},
         {-1,-1,-1,0},
         {-1,-1,-1,1},
         {1,1,1,-1},
         {1,1,1,0}}

  A = transpose LLL syz GLSM
  P2 = convexHull A
  isReflexive P2
  P = polar P2
  assert(h11OfCY P === 2)
  assert(h21OfCY P === 272)
  
  P = reflexivePolytopeData matrix topes_10
  Ts = findAllFRSTs P
  netList Ts

  -- want a normal toric variety with the rays in the same location as GLSM
  -- with this GLSM matrix.
  
  -- First, I need the triangulations...?
///


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

V = normalToricVariety(P#"rays", Ts_0#"triangulation", WeilToClass => transpose matrix P#"glsm", CoefficientRing => kk)
topologicalDataOfCY3(Ts_0, RZ)
V = normalToricVariety(Ts_0, CoefficientRing => ZZ/32003)
Va = abstractVariety(V, base())
IV = intersectionRing Va
for i in {0,2,3,4} list integral(IV_i * IV_0 * IV_1 * IV_2)
transpose matrix for x in orbits(V, 1) list (
  for i in {0,2,3,4} list integral(IV_i * IV_(x#0) * IV_(x#1) * IV_(x#2))
)
lift(oo, QQ)
posHull oo
rays oo
X = completeIntersection(V, {-toricDivisor V})
pt = base(symbol x, symbol y, symbol z, symbol w)
Xa = abstractVariety(X, pt)
IX = intersectionRing Xa

cohomologyVector(X, {1,0,0,0})
cohomologyVector(X, {0,1,0,0})
c2X = chern_2 tangentBundle Xa
integral(c2X * t_0)
integral(c2X * t_2) -- this should be, by Friedman, -6...
integral(c2X * t_3)    
integral(c2X * t_4)

cohomologyVector(X, {0,0,0,1})
cohomologyVector(completeIntersection(V, {-toricDivisor V, V_0})) -- (1,0,1)
cohomologyVector(completeIntersection(V, {-toricDivisor V, V_1})) -- (1,0,1)
cohomologyVector(completeIntersection(V, {-toricDivisor V, V_2})) -- (1,0,0)
cohomologyVector(completeIntersection(V, {-toricDivisor V, V_3})) -- (1,0,1)
cohomologyVector(completeIntersection(V, {-toricDivisor V, V_4})) -- (1,0,0)
cohomologyVector(completeIntersection(V, {-toricDivisor V, V_5})) -- (1,0,0)
cohomologyVector(completeIntersection(V, {-toricDivisor V, V_6})) -- (1,0,2)
cohomologyVector(completeIntersection(V, {-toricDivisor V, V_7})) -- (1,0,0)

h = x*t_0 + y*t_2 + z*t_3 + w*t_4

integral(c2X * t_2) == 10
integral(c2X * t_4) == -4
integral(c2X * t_5) == 16
integral(c2X * t_7) == -4

-- Try to understand these surfaces?
F = random(-degree toricDivisor V, ring V)
use ring F
F4 = sub(F, x_4 => 0)
factor F4
-- Below this still works (hopefully). (it doesn't!)
exponents (F4 // x_7^2)
oo_{0,1,2,3,5,6}
convexHull transpose oo
vertices oo
latticePoints ooo
-- How can we investigate the surface which is F = x_4 = 0.  (and so x_7 = 1)?
-- Method #1: look locally on the toric variety, at the codim 2 locus.
-- in this example, we consider the generator of ideal V: x_0 x_2 x_3 x_4
-- this turns out to be a smooth cone.
G = sub(F, {x_1 => 1, x_5 => 1, x_6 => 1, x_7 => 1})
G0 = sub(G, {x_4 => 0}) -- surface in CC^3.
singG = trim(ideal G0 + ideal jacobian G0)
codim singG
support G0
decompose singG
-- G: F = x_4 = 0, on the patch (x_0 x_2 x_3 x_4), appears to be singular
-- at x_2 = x_3 = 0.
trim(ideal G0 + ideal(x_2, x_3))
singularCones(4, V)
singularCones(3, V)
singularCones(2, V)
singularCones(1, V)
singularCones(0, V)


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

-------------------------------------------------
-- todo: given a toric variety:
--   1. give (consistent with GLSM charges): curve classes generating all (effective) curves.
--   2. write out the input to computeGV
--   3. read in the output from computeGV.

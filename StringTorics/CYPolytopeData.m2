---------------------------------------
-- CYPolytopeData ---------------------
---------------------------------------
-- This type can be written to disk, and tries to retain computations computed already.
-- It does not retain Polyhedron objects, but hopefully it recreates these quickly.
-- 

CYPolytopeFields = {
    -- first entry: true means it must exist and be in the main hash table
    --   false: it might exist, and is in the cache table.
    "rays" => {value, toString, List},
    "face dimensions" => {value, toString, List}
    }

-- These are the cache fields that we write to a string via 'dump'
CYPolytopeCache = {
    -- first entry: true means it must exist and be in the main hash table
    --   false: it might exist, and is in the cache table.
    "id" => {value, toString, ZZ},
    "favorable" => {value, toString, Boolean},
    "h11" => {value, toString, ZZ},
    "h21" => {value, toString, ZZ},
    "basis indices" => {value, toString, List},
    "glsm" => {value, toString, List},
    "annotated faces" => {value, toString, List}
    }

cyPolytopeData = method(Options => {ID => null})

cyPolytopeData Polyhedron := opts -> P2 -> (    
    LP := latticePointList P2;
    LPdim := for lp in LP list dim(P2, minimalFace(P2, lp));
    -- now remove the ones that are in facets (or the origin):
    LP = for i from 0 to #LP-1 list if LPdim#i <= 2 then LP#i else continue;
    LPdim = for i from 0 to #LP-1 list if LPdim#i <= 2 then LPdim#i else continue;
    cyData := new CYPolytopeData from {
        symbol cache => new CacheTable,
        "rays" => LP,
        "face dimensions" => LPdim
        };
    if opts.ID =!= null then cyData.cache#"id" = opts.ID;
    cyData
    )
-- vertices: Matrix whose columns are the vertices of the reflexive polytope.
cyPolytopeData Matrix := CYPolytopeData => opts -> vertices -> (
    P2 := convexHull vertices;
    cyPolytopeData(P2, opts)
    )
-- vertices: A list of the integer coordinates (also a list) of the vertices of the polytope
cyPolytopeData List := CYPolytopeData => opts -> vertices -> (
    cyPolytopeData(transpose matrix vertices, opts)
    )
cyPolytopeData KSEntry := CYPolytopeData => opts -> tope -> (
    -- KSEntry is a Kreuzer-Skarke polytope entry, returned from
    --   ReflexivePolytopesDB functions.
    P1 := convexHull matrix tope;
    P2 := polar P1;
    cyPolytopeData(P2, opts)
    )

cyPolytopeData String := CYPolytopeData => opts -> str -> (
    L := lines str;
    if L#0 != "CYPolytopeData" then error "string is not in proper format";
    fields := hashTable for i from 1 to #L-1 list getKeyPair L#i;
    -- First get the main elements (these are required!):
    required := for field in CYPolytopeFields list (
        k := field#0;
        readFcn := field#1#0;
        if fields#?k then k => readFcn fields#k else error("expected key "|k)
        );
    cyData := new CYPolytopeData from prepend(symbol cache => new CacheTable, required);
    -- now read in the cache values (including "id" value, if any)
    for field in CYPolytopeCache do (
        k := field#0;
        readFcn := field#1#0;
        if fields#?k then cyData.cache#k = readFcn fields#k;
        );
    if opts.ID =!= null then cyData.cache#"id" = opts.ID; -- just for compatibility with other constructors...
    cyData
    )

-- todo: translation function: {1, 2, 3, 6} ==> "1 2 3 6" (and viceversa)
-- todo: translation function: {{1,3},{4,7},{6,7,8}} ==> "1 3;4 7;6 7 8;" or "1 3;4 7;6 7 8" (white space is not relevant after or before a ;)
-- Format
-- CYPolytopeData
--   rays: 1 0 0; 1 0 -1; 1 1 1
--   face dimensions: 0 0 0
--   id: 12
--   favorable: true
--   h11: 5
--   h21: 20
--   basis indices: 0 1 2 3
--   glsm: 1 1 1; 1 2 3

-- Then need to be able to set fields
--
-- Need a isWellFormed function.  Checks that the correct fields are
-- present, and the lengths of the various integer vectors and lists
-- are compatible.

dump CYPolytopeData := String => {} >> opts -> (Q) -> (
    s1 := "CYPolytopeData\n";
    strs := for field in CYPolytopeFields list (
        k := field#0;
        writerFunction := field#1#1;
        if not Q#?k then error("expected key: "|k#0);
        "  " | k | ":" | writerFunction(Q#k) | "\n"
        );
    strs2 := for field in CYPolytopeCache list (
        k := field#0;
        writerFunction := field#1#1;
        if not Q.cache#?k then continue;
        "  " | k | ":" | writerFunction(Q.cache#k) | "\n"
        );
    strs = join({s1}, strs, strs2);
    concatenate strs
    )

getKeyPair = method()
getKeyPair String := Sequence => str -> (
    str1 := replace("^ *", "", str);
    result := separate(" *: *", str1); -- separate at colon, ignoring white space around colon.
    if #result != 2 then error("expected a key and a value for "|str);
    toSequence result
    )

cySetGLSM = method()
cySetGLSM CYPolytopeData := (cyData) -> elapsedTime (
    if cyData.cache#?"glsm" then return;
    mLP := transpose matrix cyData#"rays";
    D := transpose syz mLP;
    p := findFirstUnitVectors D; -- TODO: p,q computation can be slow!
    q := findInvertibleSubmatrix(D, p);
    if q === null then error ("oops, can't find a good GLSM matrix"); -- hasn't happened yet. HAS NOW!!
    GLSM := (D_q)^-1 * D;
    cyData.cache#"glsm" = entries transpose GLSM;
    cyData.cache#"basis indices" = q
    )

cySetH11H21 = cyData -> (
    -- this version is only for CY 3-fold hypersurfaces...
    -- P:ReflexivePolytope
    -- P := polytope cyData;
    elapsedTime A := annotatedFaces cyData; -- polytope on N side.
    A0 := for x in A list if x#0 == 0 then drop(x,1) else continue; -- annotatedFaces(0, P);
    A1 := for x in A list if x#0 == 1 then drop(x,1) else continue; -- annotatedFaces(1, P);
    A2 := for x in A list if x#0 == 2 then drop(x,1) else continue; -- annotatedFaces(2, P);
    A3 := for x in A list if x#0 == 3 then drop(x,1) else continue; -- annotatedFaces(3, P);
    npM := A/(x -> x#4)//sum + 1;
    npN := A/(x -> x#3)//sum; -- origin is included in the dim 4 face.
    -- points in facets (on M side) -- this is part of h21
    -- points in facets (on N side) -- this is part of h11
    facetInteriorsM := A0/(v -> v#3)//sum;
    facetInteriorsN := A3/(v -> v#2)//sum;
    -- points interior to 2-faces (times their genus) (on M-side)
    -- points interior to 2-faces (times their genus) (on N-side)
    twoFacesM := A1/(v -> v#2 * v#3)//sum;
    twoFacesN := A2/(v -> v#2 * v#3)//sum;
    -- now set the h11, h21.
    h11 := npN - 5 - facetInteriorsN + twoFacesN;
    h21 := npM - 5 - facetInteriorsM + twoFacesM;
    cyData.cache#"h11" = h11;
    cyData.cache#"h21" = h21;
    cyData.cache#"favorable" = (twoFacesN == 0);
    (h11, h21)
    )

rays CYPolytopeData := List => cyData -> cyData#"rays"
dim CYPolytopeData := List => cyData -> dim polytope(cyData, "N")
degrees CYPolytopeData := List => cyData -> (
    if not cyData.cache#?"glsm" then cySetGLSM cyData;
    cyData.cache#"glsm"
    )
basisIndices = method()
basisIndices CYPolytopeData := List => cyData -> (
    if not cyData.cache#?"basis indices" then cySetGLSM cyData;
    cyData.cache#"basis indices"
    )
-- h11OfCY CYPolytopeData := ZZ => cyData -> (
--     if not cyData.cache#?"h11" then cySetH11H21 cyData;
--     cyData.cache#"h11"
--     )
-- h21OfCY CYPolytopeData := ZZ => cyData -> (
--     if not cyData.cache#?"h21" then cySetH11H21 cyData;
--     cyData.cache#"h21"
--     )
isFavorable CYPolytopeData := Boolean => cyData -> (
    if not cyData.cache#?"favorable" then cySetH11H21 cyData;
    cyData.cache#"favorable"
    )
annotatedFaces CYPolytopeData := cyData -> (
    if not cyData.cache#?"annotated faces" then
      cyData.cache#"annotated faces" = annotatedFaces polytope(cyData, "N");
    cyData.cache#"annotated faces"
    )
polytope(CYPolytopeData, String) := Polyhedron => (cyData, which) -> (
    if which === "N" then (
        if not cyData.cache#?"N polytope" then (
            LP := cyData#"rays";
            LPdim := cyData#"face dimensions";
            verts := for i from 0 to #LP - 1 list if LPdim#0 == 0 then LP#i else continue;
            cyData.cache#"N polytope" = convexHull transpose matrix verts
            );
        cyData.cache#"N polytope"
        )
    else if which === "M" then (
        if not cyData.cache#?"M polytope" then (
            cyData.cache#"M polytope" = polar polytope(cyData, "N");
            );
        cyData.cache#"M polytope"
        )
    else
      error "expected second argument to be either \"M\" or \"N\""
    )
polytope CYPolytopeData := Polyhedron => cyData -> polytope(cyData, "N")

polar CYPolytopeData := cyData -> cyPolytopeData polytope(cyData, "M")

findAllFRSTs CYPolytopeData := List => cyData -> (findAllFRSTs(transpose matrix rays cyData))/last

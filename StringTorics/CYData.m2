--------------------------------------------------------------
-- CalabiYauInToric (soon to change back to CalabiYauInToric? ----------
--------------------------------------------------------------

CalabiYauInToric.synonym = "Calabi-Yau in a normal toric variety"
CalabiYauInToric.GlobalAssignHook = globalAssignFunction
CalabiYauInToric.GlobalReleaseHook = globalReleaseFunction
expression CalabiYauInToric := X -> if hasAttribute (X, ReverseDictionary) 
    then expression getAttribute (X, ReverseDictionary) else 
    (describe X)#0
describe CalabiYauInToric := X -> Describe (expression CalabiYauInToric) (
      expression "a" , expression 3)
--    expression rays X, expression max X)

CYDataFields = {
    -- first entry: true means it must exist and be in the main hash table
    --   false: it might exist, and is in the cache table.
    "polytope data" => {value, Q -> toString Q.cache#"id", CYPolytopeData},
    "triangulation" => {value, toString, List}
    }

-- These are the cache fields that we write to a string via 'dump'
CYDataCache = {
    -- first entry: true means it must exist and be in the main hash table
    --   false: it might exist, and is in the cache table.
    "id" => {value, toString, ZZ},
    "c2" => {value, toString, List},
    "intersection numbers" => {value, toString, List},
    "toric intersection numbers" => {value, toString, List}
    }

cyData = method(Options => {ID => null, Ring => null})
cyData(CYPolytopeData, List) := opts -> (Q, triang) -> (
    X := new CalabiYauInToric from {
        symbol cache => new CacheTable,
        "polytope data" => Q,
        "triangulation" => triang
        };
    if opts.ID =!= null then X.cache#"id" = opts.ID;
    if opts#Ring =!= null then X.cache#"pic ring" = opts#Ring; -- Note: we should check that it is over ZZ, has h^(1,1) variables
    X
    )
cyData(String, Function) := CalabiYauInToric => opts -> (str, F) -> (
    -- F is a function which takes an id of a CYPolytopeData and returns the object.
    L := lines str;
    if L#0 != "CYData" then error "string is not in proper format";
    fields := hashTable for i from 1 to #L-1 list getKeyPair L#i;
    -- First get the main elements (these are required!):
    polytopeid := value fields#"polytope data";
    required := for field in CYDataFields list (
        k := field#0;
        if k === "polytope data" then (
            "polytope data" => F polytopeid
            )
        else (
            readFcn := field#1#0;
            if fields#?k then k => readFcn fields#k else error("expected key "|k)
        ));
    X := new CalabiYauInToric from prepend(symbol cache => new CacheTable, required);
    -- now read in the cache values (including "id" value, if any)
    for field in CYDataCache do (
        k := field#0;
        readFcn := field#1#0;
        if fields#?k then X.cache#k = readFcn fields#k;
        );
    if opts.ID =!= null then X.cache#"id" = opts.ID; -- just for compatibility with other constructors...
    if opts#Ring =!= null then X.cache#"pic ring" = opts#Ring; -- Note: we should check that it is over ZZ, has h^(1,1) variables
    X
    )

dump CalabiYauInToric := String => {} >> opts -> X -> (
    s1 := "CYData\n";
    strs := for field in CYDataFields list (
        k := field#0;
        writerFunction := field#1#1;
        if not X#?k then error("expected key: "|k#0);
        "  " | k | ":" | writerFunction(X#k) | "\n"
        );
    strs2 := for field in CYDataCache list (
        k := field#0;
        writerFunction := field#1#1;
        if not X.cache#?k then continue;
        "  " | k | ":" | writerFunction(X.cache#k) | "\n"
        );
    strs = join({s1}, strs, strs2);
    concatenate strs
    )

makeCY = method(Options => {ID => null, Ring => null})
makeCY CYPolytopeData := CalabiYauInToric => opts -> Q -> (
    P2 := polytope Q;
    (LP,tri) := regularStarTriangulation(dim P2-2,P2);
    if rays Q =!= LP then error "I have a lattice point mismatch";
    cyData(Q, tri, opts)
    )    

normalToricVariety CalabiYauInToric := opts -> X -> (
    if not X.cache.?NormalToricVariety then X.cache.NormalToricVariety = (
        Q := X#"polytope data";
        T := X#"triangulation";
        GLSM := transpose matrix degrees Q;
        normalToricVariety(rays Q, T, opts, WeilToClass => matrix GLSM)
        );
    X.cache.NormalToricVariety
    -- TODO: this fails if the class group is torsion! (Fails: later it gives an inscrutable error...)
    )

rays CalabiYauInToric := X -> rays cyPolytopeData X
max CalabiYauInToric := X -> X#"triangulation"

-- TODO: triangulation is used with 2 different pieces of data:
--  with, without cone point!  Change this to use only one point.
-- Also: there are 4 matrices one can imagine: A, A0 (A with origin), Ah, A0h...
-- We need to be consistent about these!
-- TODO: do we really need this?
triangulation CalabiYauInToric := Triangulation => opts -> X -> (
    if not opts.Homogenize then error "Homogenize flag is not used in this method";
    if not X.cache#?"triangulation" then (
        rys := X#"polytope data"#"rays";
        d := #rys#0;
        B := (transpose matrix rys) | matrix{d:{0}};
        X.cache#"triangulation" = triangulation(B, for t in X#"triangulation" list append(t, #rys));
        );
    X.cache#"triangulation"
    )

cyPolytopeData CalabiYauInToric := opts -> X -> X#"polytope data"
dim CalabiYauInToric := X -> dim ambient X - 1
polytope CalabiYauInToric := X -> polytope cyPolytopeData X
polytope(CalabiYauInToric, String) := (X, which) -> polytope(cyPolytopeData X, which)
basisIndices CalabiYauInToric := List => X -> basisIndices cyPolytopeData X
degrees CalabiYauInToric := List => X -> degrees cyPolytopeData X

ambient CalabiYauInToric := X -> normalToricVariety X

label = method()
label CYPolytopeData := Q -> if Q.cache#?"id" then Q.cache#"id" else ""
label CalabiYauInToric := X -> (label cyPolytopeData X, if X.cache#?"id" then X.cache#"id" else "")

hh(Sequence, CalabiYauInToric) := (pq, X) -> hh^pq cyPolytopeData X

abstractVariety CalabiYauInToric := opts -> X -> (
    -- Store this with X.
    V := ambient X;
    aX := completeIntersection(V, {-toricDivisor V});
    abstractVariety(aX, base())
    )
abstractVariety(CalabiYauInToric, AbstractVariety) := opts -> (X, pt) -> (
    -- Store this with X?
    -- Check: pt is of dimension zero?
    V := ambient X;
    aX := completeIntersection(V, {-toricDivisor V});
    abstractVariety(aX, pt)
    )

-- restrictTriangulation: returns a List of
--   {2-face indices, 
--    all indices of points in in this 2-face, 
--    the triangles in this 2-face, 
--    genus of this face}
-- TODO: need also a function which returns just: triangles, genus information.
restrictTriangulation = method()
restrictTriangulation CalabiYauInToric := List => (X) -> (
    -- given X, we use its annotated faces and its triangulation, to write down the triangulations of the 2-faces
    -- of the corresponding reflexive polytope in the N lattice side.
    Q := cyPolytopeData X;
    F := annotatedFaces Q;
    twofaces := for x in F list if x#0 =!= 2 then continue else {x#1, x#2, x#4};
    T := max X; -- triangulation
    for t2 in twofaces list (
        a := set t2#1; -- these are the indices we want.
        atri := sort unique for t in T list (
            b := sort toList(a * set t);
            if #b == 3 then b else continue
            );
        {t2#0, t2#1, atri, t2#2}
        )
    )

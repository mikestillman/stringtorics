newPackage(
    "PALPInterface",
    Version => "0.1",
    Date => "",
    Headline => "An interface to parts of the PALP polyhedra software",
    Authors => {
        {Name => "Mike Stillman", 
            Email => "mike@math.cornell.edu", 
            HomePage => "http://pi.math.cornell.edu/~mike"}
        },
    AuxiliaryFiles => false,
    DebuggingMode => true,
    PackageImports => {"ReflexivePolytopesDB", "Polyhedra", "NormalToricVarieties"}
    )

export {
    -- weight systems
    "getWSFromDim",
    "getVerticesFromWS",
    "weightSystemPolytope",
    
    -- nef partition code
    "palpNefPartitions",
    "parsePALPNefPartitions",
      --"getNefPartitionInfo",

    -- normal form
    "normalForm",

    -- some basic palp functions
    "palpMVertices",
    "palpMPoints",
    "palpNVertices",
    "palpNPoints",
    --"palpVertices",

    "stringToList",
    --"getPartitionInfo", -- nef partition info from nef.x output
    --"readPartitions",
    "formDivisor",
    "getPolyhedralInfo",
    "toPalp",
    "fromPalpMatrix",
    "runPoly", -- TODO: change name, or remove this function, once we understand output of PALP better.
    "runNEF",
    "Normal"
    }

exportMutable {
    "palpVerbosity"
    }

palpVerbosity = 0

POLYX = findProgram("poly.x", "poly.x -h")
NEFX = findProgram("nef.x", "nef.x -h")
CWSX = findProgram("cws.x", "cws.x -h")
polyexe =  "!"|POLYX#"path"| POLYX#"name" -- Don't use this

NEF11D = findProgram("nef-11d.x", "nef-11d.x -h")
--programPaths#"poly.x" = executableDir; -- Don't use this.

nextInt = 0 -- used as a hack for now to create file input names...

runPALP = method()
runPALP(Program, String, String) := (P, args, hereis) -> (
    thisInt := nextInt;
    nextInt = nextInt + 1;
    filename := "run-palp"|thisInt;
    filename << hereis << close;
    result := runProgram(P, args | " " | filename);
    if palpVerbosity > 1 then print result#"output";
    result#"output"
    )

-- translate Matrix or a (combined) weight system list into a String suitable for input to PALP
-- warning: often, the convex hull of the columns must be reflexive.
-- we need to document when this is necessary
-- Additionally, the matrix must satisfy that the polytope is full dimensional
toPalp = method()
toPalp Matrix := String => M -> (
    if ring M =!= ZZ then error "expected integer matrix";
    header := toString(numRows M | " " | numcols M);
    e := entries M;
    es := for e1 in e list for a in e1 list (toString a | " ");
    s := concatenate between("\n", es);
    header | "\n" | s | "\n"
    )
toPalp List := String => wlist -> (
    w1 := for i in wlist list(toString(i)|" ");
    w := concatenate(drop(w1,-1),toString(wlist_(-1)));
    w = w | "\n"
    )

-- From output of palp functions, we might need to grab only some of the lines
-- to form a matrix.
fromPalpMatrix = method()
fromPalpMatrix String := Matrix => str -> (
    matrix KSEntry str -- this will give a somewhat inscrutable error if the format is not correct
    )

normalForm = method()
normalForm Matrix :=
normalForm List := Matrix => M ->
    fromPalpMatrix runPALP(POLYX, "-N", toPalp M)  -- -N: normal form

getVerticesFromWS = method()
getVerticesFromWS List := Matrix => wlist -> (
    -- wlist is of the form {sum, a1, a2, a3, ..., an}
    -- e.g. getVerticesFromWS {10,1,2,3,4}
    w1 := for i in wlist list(toString(i)|" ");
    w := concatenate(drop(w1,-1),toString(wlist_(-1)));
    w = w | "\n";
    fromPalpMatrix runPALP(POLYX, "-v", w)
    )

palpMVertices = method()
palpNVertices = method()
palpMPoints = method()
palpNPoints = method()

palpMVertices Matrix := 
palpMVertices List := Matrix => A -> 
  fromPalpMatrix runPALP(POLYX, "-v", toPalp A)  -- -v: M-lattice vertices

palpMPoints Matrix := 
palpMPoints List := Matrix => A -> 
  fromPalpMatrix runPALP(POLYX, "-p", toPalp A)  -- -p: M-lattice points
  
palpNVertices Matrix := 
palpNVertices List := Matrix => A -> 
  fromPalpMatrix runPALP(POLYX, "-e", toPalp A)  -- -v: N-lattice vertices

palpNPoints Matrix := 
palpNPoints List := Matrix => A -> (
  output := runPALP(POLYX, "-d", toPalp A);  -- -d: N-lattice points
  if output === "" then error "the convex hull of the given matrix is likely not reflexive (or there was some other error";
  fromPalpMatrix output
  )
  

---------------------------
-- Weight systems code ----
---------------------------

-- this function doesn't seem to be correct at all?
weightSystemPolytope = method()
weightSystemPolytope List := (wts) -> (
    -- Note: wts#0 is the sum of the rest.
    -- create it from the definition, using Polyhedra.
    -- used to test the output of PALP...
    --
    wts0 := drop(wts, 1);
    H := matrix{wts0};
    d := #wts0-1;
    halfspaces := -id_(ZZ^(d+1));
    ones := matrix{(d+1):{1}};
    Nabla := polyhedronFromHData(halfspaces, ones, H, matrix{{0}});
    P := polar Nabla;
    vP := vertices P;
    if vP^{0} != 0 then error "my logic is wrong";
    P0 := convexHull ((vP)^{1..numrows vP - 1});
    convexHull matrix{latticePoints P0}
    )

getWSFromDim = method(Options => {Degrees => null})
getWSFromDim ZZ := List => opts -> d -> (
    drange := opts.Degrees;
    if drange =!= null then (
        -- check consistency of opts.Degrees
        if instance(drange, ZZ) then drange = {drange, drange};
        if not instance(drange, BasicList) or #drange =!= 2 then error "expected Degrees => {lo, hi} or Degrees => deg";
        if not all(drange, a -> instance(a, ZZ)) then error "expected list of integers";
        );
    str0 := toString(d);
    str1 := if drange === null then "!cws.x -w"|str0
    else "!cws.x -w"|str0|" "|toString(drange_0)|" "|toString(drange_1);

    PALPOutput := get str1;
    L := lines PALPOutput;
    L1 := if opts.Degrees === null then L else drop(L,-1);
    M := for ell in L1 list(
	L0 := separate(" +",ell);
	L0Mod := take(L0, {0,d+1});
	for x in L0Mod list value x
	)
    )

------------------------------
-- NEF partition code --------
------------------------------

-- read the output of the nef partition code
palpNefPartitions = method()
palpNefPartitions(Matrix, ZZ) :=
palpNefPartitions(List, ZZ) := (M, cod) -> (
    runPALP(NEF11D, "-Lp -c"|cod, toPalp M)
    )

-- deprecate this.
getNefPartitionInfo = method()
getNefPartitionInfo(Matrix, ZZ) :=
getNefPartitionInfo(List, ZZ) := String => (M, cod) -> (
    runPALP(NEF11D, "-Lp -c"|cod, toPalp M)
    )

-- This needs some work, as it is not grabbing the entire parts in cod=3 (or higher).
parsePALPNefPartitions = method()
parsePALPNefPartitions String := output1 -> (
    --reg1 := "V[0-9]*:([0-9 ]+ {1})";
    reg1 := "V[0-9]*:([0-9 ]+ {1})";
    for thisline in (lines output1) list (
        regres := regex(reg1, thisline);
        if not (regres === null) then (
            << thisline << " " << thisline_(regres_1) <<"\n";
            -- for eachparti in thisline_(regres_1) list (regex("^[0-9]$", eachparti))
            stringToList(thisline_(regres_1))
        )
        else continue
        )
    )

-- Remove this after parsePALPNefPartitions is working correctly.
readPartitions = (output1) -> (
    reg1 := "P:[0-9 ]+V:([0-9 ]+ {1})";
    for thisline in (lines output1) list(
        -- reg1 = "M:(.*)N(.*)codim(.*)part*";
        
        -- (for x in L0 list if x!="" then value x else continue)
        regres := regex(reg1, thisline);
        -- << "start: " << thisline;
        -- << "reg: ";
        -- << thisline;
        -- << "\n";
        -- -- if not (regres === null) then << regres;
        -- << regex(reg1, thisline);
        -- << "\n";
        -- << "\n";
        -- << output1_(219, 4);
        if not (regres === null) then (
            << thisline << " " << thisline_(regres_1) <<"\n";
            -- for eachparti in thisline_(regres_1) list (regex("^[0-9]$", eachparti))
            stringToList(thisline_(regres_1))
        )
        else continue
    )
)


-- probably remove?
runPoly = method()
runPoly(Matrix, String) := String => (M, opts) -> (
    "foo" << toPalp M << close;
    cmd := " -" | opts | " foo";
    result := runProgram(POLYX, cmd);
    result#"output"
    )

-- probably remove?
runNEF = method()
runNEF(Matrix, String) := String => (M, opts) -> (
    "foo" << toPalp M << close;
    cmd := " -" | opts | " foo";
    result := runProgram(NEFX, cmd);
    result#"output"
    )

-- two versions here?
-- keep these?
getPolyhedralInfo = method(Options => {Normal => false})
getPolyhedralInfo Matrix := opts -> M -> (
    -- idea: create the matrix of M.
    -- create the call to PALP
    -- call PALP (poly.x here)
    -- parse the output and return the answer as a hash table of desired info.
    w := toPalp M;
    -- this needs to be changed depending on desired info
    str1 := "!poly.x -v  << FOO\n";
    --str1 := " -v << FOO\n";
    str2 := "\nFOO\n";
    str3 := str1|w|str2;
    PALPOutput := get str3;
    -- the following obtains the matrix remaining after removing the given lines.
    keeplines := select(lines PALPOutput, s -> (
            not match("^Degrees", s) and not match("^Type", s) and not match("or", s)));
    str := concatenate between("\n", keeplines);
    matrix KSEntry str
    )
getPolyhedralInfo Matrix := opts -> M -> (
    -- idea: create the matrix of M.
    -- create the call to PALP
    -- call PALP (poly.x here)
    -- parse the output and return the answer as a hash table of desired info.
    w := toPalp M;
    -- this needs to be changed depending on desired info
    str1 := "!poly.x -g  << FOO\n";
    --str1 := " -v << FOO\n";
    str2 := "\nFOO\n";
    str3 := str1|w|str2;
    PALPOutput := get str3;
    return PALPOutput;
    -- the following obtains the matrix remaining after removing the given lines.
    keeplines := select(lines PALPOutput, s -> (
            not match("^Degrees", s) and not match("^Type", s) and not match("or", s)));
    str := concatenate between("\n", keeplines);
    matrix KSEntry str
    )

stringToList = method()
stringToList String := List => stemp -> (
    L0 := separate(" +", stemp);
    (for x in L0 list if x!="" then value x else continue)
)


-- TODO: What is this function doing?
formDivisor = (V, partlist) -> (
    raylist := rays V;
    -- << # raylist;
    -- -- << raylist_0;
    -- << raylist;
    -- << "\n";
    -- << partlist;
    -- << "\n";
    fulllist := toList (0..(# raylist - 1));
    partlist = reverse partlist;
    complist := fulllist;

    for thispart in partlist do complist = drop(complist, {thispart, thispart});
    -- << complist;
    D1list := for thispart in partlist list (V_thispart);
    -- << "D1list: " << for thispart in partlist list ("V_"|toString(thispart));
    D1 := sum(D1list);
    -- << D1;
    
    D2list := for thispart in complist list (V_thispart);
    -- << "D2list: " << for thispart in complist list ("V_"|toString(thispart));
    D2 := sum(D2list);
    -- << D2;

    -- << isNef(D1);
    -- << isNef(D2);
    if not isNef(D1) then (<< "D1 false");
    if not isNef(D2) then (<< "D2 false");
    (D1, D2)
    -- completeIntersection(V, {D1, D2})

    -- test5 = (for thispart in partlist list drop(fulllist, {thispart, thispart}));
    -- test5
    -- complist = 
    -- for thispart in partlist do << thispart;


)

-- Not using palp. Maybe: naiveNefPartitions...?
-- anyway: this is codim=2 only...
findNefPartitions = method()
findNefPartitions(ZZ, NormalToricVariety) := (cod, V) -> (
    subsetRays := subsets(splice{0..#rays V-1});
    subsetRays = drop(subsetRays, 1); -- remove the empty set
    subsetRays = drop(subsetRays, -1); -- remove the full set
    partition(v -> isNef sum(v, i -> V_i), subsetRays, {true, false})
    )

beginDocumentation()

doc ///
Key
  PALPInterface
Headline
  interface to the polyhedral program PALP
Description
  Text
    The program PALP was designed by ... to help compute all reflexive polytopes in dimensions 3 and 4.
    However, it has some generally useful functionality beyond that.  This package interfaces to that
    functionality.
///

doc ///
  Key
    (getVerticesFromWS, List)
  Headline
    vertices of polytope of a weight system
  Usage
    getVerticesFromWS q
  Inputs
    q:List
      of integers: the first is the sum of the rest, and all but the first are positive integers in
      ascending order
  Outputs
    :Matrix
      over the integers, the columns are the vertices of $\Delta(q)$.
  Description
    Text
        This function returns the vertices of the polytope $\Delta(q)$
        (which appears only defined up to an integer change of basis).
    Example
        q = {3,4,5,14,21}
          -- q = {40,41,486,1134,1701}
        ws = prepend(sum q, q)
        verts = getVerticesFromWS ws
        needsPackage "Polyhedra"
        P = convexHull verts
        assert isReflexive P
        assert(dim P == 4)
        numcols  vertices P == 18
        length latticePoints P == 54
        nfverts = normalForm verts
    Example
        A = transpose LLL syz matrix{q}
        B = convexHull latticePoints polar convexHull A
        C = polar B
        assert not isReflexive convexHull A 
        assert isReflexive B
        assert isReflexive C
        nflll = normalForm lift(vertices B, ZZ)
        assert(nfverts == nflll)
    Text
        Let's compare two weight systems to see if they are the same.
    Example
        ws1 = {3402, 40, 41, 486, 1134, 1701}
        ws2 = {3486, 41, 42, 498, 1162, 1743}
        M1 = normalForm getVerticesFromWS ws1
        M2 = normalForm getVerticesFromWS ws2
        M1 == normalForm ws1
        M2 == normalForm ws2
        M1 == M2
  SeeAlso
    getWSFromDim
    normalForm
///

-*
  restart
  needsPackage "PALPInterface"
*-
TEST ///
  -- test of toPalp, fromPalpMatrix
  M = matrix{{1,1,1,1},{0,1,2,3}}
  str1 = toPalp M
  assert(M == fromPalpMatrix str1)

  str2 = toPalp transpose M
  assert(transpose M == fromPalpMatrix str2)

  -- from smoothFanoToricVariety(3, 5)
  M = transpose matrix {{1, 0, 0}, {-1, 0, 1}, {0, 1, 0}, {0, -1, 1}, {0, 0, 1}, {0, 0, -1}}
  assert(M == fromPalpMatrix toPalp M)
///

-*
  restart
  needsPackage "PALPInterface"
*-
TEST ///
  -- test of palpNPoints, palpNVertices, palpMPoints, palpNpoints
  ws = {10,1,2,3,4}
  M = getVerticesFromWS ws

  -- Do these work for non-reflexive polytopes?
  mv = palpMVertices M -- default input is M lattice.
  mp = palpMPoints M
  nv = palpNVertices M
  np = palpNPoints M

  assert(M == mv)
  assert(mv == palpMVertices ws)
  mp == palpMPoints ws -- not the same order as mp! (At least generally).
  assert(set entries mp === set entries palpMPoints ws)
  assert(nv == palpNVertices ws)
  np == palpNPoints ws -- true here, doesn't need to be.
  assert(set entries np === set entries palpNPoints ws)
///

-*
  restart
  needsPackage "PALPInterface"
*-
TEST ///
  --needsPackage "QuillenSuslin"
  -- A = ZZ[x]
  -- A1 = completeMatrix matrix(A, {{6,10,3*26, 3*39}})
  -- assert(A1 == matrix(A, {{6, 10, 78, 117}, {0, 3, 0, 35}, {0, 0, 1, 0}, {1, 0, 0, 0}}))

  -- another way to generate invertible integer matrices
  -- needsPackage "IntegerEquivalences"
  -- A = extendToMatrix {6,10,3*26, 3*39}

  A = matrix {{0, 0, -3, 2}, {-5, 3, 0, 0}, {-3, -6, 1, 0}, {-3, -2, -1, 1}}  
  M = fromPalpMatrix "4 12  M:24 12 N:15 12 H:11,19 [-16] id:13
   1   0   0   1   1  -1   0  -2   4  -2   0   2
   0   1   0   0  -1   0   0   3  -4   1  -1  -3
   0   0   1  -1   0   0   0   1  -4   3  -1  -3
   0   0   0   0   0   0   1  -1   1  -1   1   1
   "
   assert(numrows M == 4 and numcols M == 12)
   assert(normalForm M == M)
   assert(normalForm (A*M) == M)
///

TEST ///
-- This test is not quote correct yet.  weightSystemPolytope returns a polytope with full dim?
  -- this polytope is in a hyperplane in one higher dimension...
  vertices convexHull matrix{latticePoints polar weightSystemPolytope{10, 1,2,3,4}}
  normalForm lift(oo, ZZ)

    normalForm getVerticesFromWS{10, 1,2,3,4}
  assert(oo == ooo)


  ws = {30, 4, 4, 6, 7, 9}
  vertices convexHull matrix{latticePoints polar weightSystemPolytope ws}
  normalForm lift(oo, ZZ)
  normalForm getVerticesFromWS ws
  assert(oo == ooo)

  P = polar weightSystemPolytope {3,1,1,1}

  vertices P
  isReflexive P
  assert(#latticePoints P == 10)
  --assert(vertices P == matrix(QQ, {{2, -1, -1}, {-1, 2, -1}, {-1, -1, 2}})) -- TODO: need better test here.
  vertices polar P

  normalForm getVerticesFromWS({3,1,1,1})
  
  -- TODO: what are these really supposed to be?
  getVerticesFromWS({10,1,2,3,4})
  Q = polar weightSystemPolytope {10,1,2,3,4}
  vertices Q
  latticePoints Q
///

-*
  restart
  needsPackage "PALPInterface"
*-
TEST ///
  -- which operations are fine with non-reflexive polytopes?
  M = transpose matrix{{1,-1,-1}, {1,-1,1}, {1,1,-1}, {1,1,2}, {-1,0,0}}
  P = convexHull M
  isReflexive P
  vertices polar P
  interiorLatticePoints P -- this is an IP polytope...
  -- Now, let's see what palp does with this

  nM = normalForm M
  Mpermuted = M_{4,1,3,0,2} 
  basechange = nM_{0,1,2} * (Mpermuted_{0,1,2} ** QQ)^-1
  basechange = lift(basechange, ZZ)
  assert(nM == basechange * M_{4,1,3,0,2})

  -- anyway, normalForm seems to be fine here.
  assert(palpMVertices M == M)
  assert(numcols palpMPoints M == # latticePoints convexHull M)
  assert(set entries transpose palpMPoints M === set entries transpose matrix {latticePoints convexHull M})

  -- how about N? Yes, the equtions are different though!  They are equations over ZZ, but have an extra column
  -- POSSIBLE TODO: handle the output in a different way?
  vertices polar convexHull M
  palpNVertices M

  -- lattice points inside the dual doesn't work though...
  ans = trap palpNPoints M
  assert(ans#0 === null and instance(ans#1, Error))
///


-*
-- XXX
  restart
  needsPackage "PALPInterface"
*-
TEST /// -- nef partition code
  -- first, let's do codim=2
  ws = {10, 1, 1, 2, 2, 2, 2}
  M = getVerticesFromWS ws
  -- want nef partitions from this
  debug PALPInterface -- getNefPartitionInfo not there yet.
  palpNefPartitions(ws, 2)
  parsePALPNefPartitions oo


  palpNefPartitions(ws, 2)
  parsePALPNefPartitions oo  
  palpNefPartitions(ws, 3)
  parsePALPNefPartitions oo
  
  getNefPartitionInfo(M, 2)
  elapsedTime getNefPartitionInfo(ws, 2) -- not so cheap (0.5 sec)

  elapsedTime getNefPartitionInfo(ws, 3) -- not so cheap (0.5 sec)
  readPartitions oo

  ws = {5, 1, 1, 1, 1, 1, 0, 0,  10, 2, 2, 2, 2, 0, 1, 1}
  elapsedTime getNefPartitionInfo(ws, 2)
  readPartitions oo

  -- codim=3 example
  options getWSFromDim
  wss = getWSFromDim(6, Degrees => {14,14})
  getNefPartitionInfo({14, 1, 2, 2, 2, 2, 2, 3}, 3) -- can't do that with NEFX

  runPALP(NEFX, "-Lp -c2", toPalp ws)
    elapsedTime runPALP(NEFX, "-c2", toPalp ws)
    elapsedTime runPALP(NEFX, "-p -c2", toPalp ws)
    elapsedTime runPALP(NEFX, "-h", toPalp ws)
    elapsedTime runPALP(NEFX, "-D -P -c2", toPalp ws)
    elapsedTime runPALP(NEFX, "-y -c2", toPalp ws) -- has nontrivial nef partition.
    elapsedTime runPALP(NEFX, "-H -c2", toPalp ws) -- display out the full Hodge diamond
///
----------------------------------------------------------------------
----------------------------------------------------------------------
--- Below this line has not been re-vetted ---------------------------
--- Todo: get everything below this line back into the package
----------------------------------------------------------------------
----------------------------------------------------------------------


-*
  restart
  needsPackage "PALPInterface"
*-
TEST ///
  -- 3d reflexive example
  -- from smoothFanoToricVariety(3, 5)
  M = transpose matrix {{1, 0, 0}, {-1, 0, 1}, {0, 1, 0}, {0, -1, 1}, {0, 0, 1}, {0, 0, -1}}
  assert(M == fromPalpMatrix toPalp M)
  normalForm M
  runPoly(M, "t") -- t means give info about the normal form process.

  needsPackage "Polyhedra"
  P = convexHull M
  P' = polar P
  M1 = lift(vertices P', ZZ)

  M = transpose matrix {{1, 0, 0}, {-1, 0, 1}, {0, 1, 0}, {0, -1, 1}, {0, 0, 1}, {0, -1, -2}}
  assert(M == fromPalpMatrix toPalp M)
  normalForm M
  runPoly(M, "t")

  
  needsPackage "IntegerEquivalences"
  A = extendToMatrix{3,4,5}
  A = matrix {{1, -2, 1}, {2, 1, -2}, {-1, 1, 0}}
  assert(det A == 1)
  assert(normalForm (A * M1) == normalForm M1)

  runPoly(M1, "p")
  runPoly(M, "p")

  runPoly(M1, "v")
  runPoly(M, "v")

  runPoly(M, "e")
  runPoly(M1, "e")  

  runPoly(M, "m")
  runPoly(M1, "m")  

  runPoly(M, "g") -- what does output mean? -- M:7 6 N:29 8 Pic:17 Cor:0
  runPoly(M1, "g") -- M:29 8 N:7 6 Pic:3 Cor:0

  runPoly(M, "a")
  runPoly(M1, "a")

  runPoly(M, "i")
  runPoly(M1, "i")

  runPoly(M, "I")
  runPoly(M1, "I")
  runPoly(matrix{{1,1,1,1},{0,1,2,3}}, "I") == "No IP\n"
  runPoly(transpose matrix {{1, 0, 0}, {-1, 0, 1}, {0, 1, 0}, {0, -1, 1}, {0, 0, 1}, {0, 0, -1}}, "I")

  runPoly(M, "S")
  runPoly(M1, "S")

  runPoly(M, "vT")
  runPoly(M, "pT")
  runPoly(M1, "vT")
  runPoly(M1, "pT")

  runPoly(M, "V")
  runPoly(M1, "V")

  runPoly(M, "B")
  runPoly(M1, "B")

  runPoly(M, "F")
  runPoly(M1, "F")

  runPoly(M, "A")
  runPoly(M1, "A")

  runPoly(M, "G")
  runPoly(2*M1, "G")
  ///

-*
  restart
  needsPackage "PALPInterface"
  -- TEST of generating polytope from a weight system
*-
TEST ///
  needsPackage "Polyhedra"
-- this polytope is in a hyperplane in one higher dimension...
  P2 = weightSystemPolytope {2,1,1}
  -- getVerticesFromWS {2,1,1}   -- hmm, this fails... too small?

  -- what is P??
  Q = convexHull getVerticesFromWS {3,1,1,1}
  vertices Q

  -- We would like to make sure we can easily go from P to Q  

  getWSFromDim 3
  
  ws = {3,1,1,1}
  M = getVerticesFromWS ws
  wt = transpose matrix{drop(ws, 1)}
  M2 = transpose LLL syz transpose wt  

  normalForm M
  normalForm M2
  -- Q should be projection from P, and it is.

  ws = {66, 5, 6, 22, 33}  -- last one on dim=3 list.
  M = getVerticesFromWS ws
  Q = convexHull M
  vertices Q
  latticePoints Q
  P = weightSystemPolytope ws
  vertices P

///
-- str = get "!poly.x -v -r <<FOO
-- 10 1 2 3 4
-- FOO
-- "

-- L = lines str
-- L = drop(drop(L, 3), -2)
-- netList L
-- L0 = separate(" +", L_0)
-- for x in drop(L0, 1) list value x
-- M = matrix for ell in L list (
--     L0 := separate(" +", ell);
--     for x in drop(L0, 1) list value x
--     )
-- needsPackage "Polyhedra"
-- needsPackage "StringTorics"
-- P = convexHull M
-- vertices P
-- isReflexive P

-- from Matrix to palp form of a matrix.
-- from palp form of a matrix to Matrix.
-- from ws or cws to matrix whose columns are the vertices of the corresponding reflexive polytope
-- from reflexive polytope vertices to ws or cws... (?)
-- optional: HodgeNumbers: h11, h12, ..., h1(dim X - 1).
-- nef partitions
-- 



///
-*
  restart
  needsPackage "PALPInterface"
*-
  wss = getWSFromDim 3
  assert(#wss == 95)
  ws = wss_90
  assert(ws == {44, 4, 5, 13, 22})  

  M = getVerticesFromWS ws
  wt = transpose matrix{drop(ws, 1)}
  M * wt
  M2 = transpose LLL syz transpose wt  

  normalForm M2
  normalForm M
  
  wss = getWSFromDim(5, Degrees => (20,20))
  assert(#wss == 61)
  V = getVerticesFromWS wss_58
  needsPackage "Polyhedra"
  P = convexHull V
  isReflexive P
  latticePoints P
  isSimplicial polar P
///

-- template for doc nodes for methods/functions
///
Key
Headline
Usage
Inputs
Outputs
Description
  Text
  Example
SeeAlso
///

-*
  restart
  needsPackage "PALPInterface"
*-
TEST ///
  M = matrix({{ -1,  -1,  -1,  -1,   1},
       { -1,  -1,   0,   3,  -1},
       { -1,  -1,   4,   0,  -1},
       { -1,   2,  -1,   2,  -1},
       { -1,   6,   0,  -1,  -1},
       {  1,   2,  -1,   2,  -1},
       {  1,   6,   0,  -1,  -1},
       { -1,  -1,  -1,  -1,  -1},
       { -1,  -1,  -1,   3,  -1},
       { -1,  -1,   5,  -1,  -1},
       { -1,   0,  -1,   3,  -1},
       { -1,   2,   3,  -1,  -1},
       { -1,   3,   1,   0,  -1},
       { -1,   7,  -1,  -1,  -1},
       {  1,   2,   3,  -1,  -1},
       {  5,   0,  -1,   3,  -1},
       {  7,  -1,   5,  -1,  -1},
       {  7,   7,  -1,  -1,  -1},
       { 23,  -1,  -1,   3,  -1},
       {151,  -1,  -1,  -1,  -1}})
  needsPackage "Polyhedra"
  P = convexHull transpose M
  assert(dim P == 5)
  isReflexive P
  runPoly(M, "gD")

  runNEF(M, "c2") -- nothing??

  wss = getWSFromDim(5, Degrees => (30,30))
  #wss == 354
  ws = wss#10
  for ws in take(wss,10) list (
    M = getVerticesFromWS(ws);
    if M === null then continue;
    if not isReflexive convexHull M then << "note: " << ws << " does not give reflexive" << endl;
    nefs = runNEF(M, "c2");
    print nefs;
    ws => nefs)

  wss = getWSFromDim(5, Degrees => (10,10))
  #wss === 5
  for ws in wss list (
    M = getVerticesFromWS(ws);
    if M === null then continue;
    if not isReflexive convexHull M then << "note: " << ws << " does not give reflexive" << endl;
    nefs = runNEF(M, "c2");
    print nefs;
    ws => nefs)

///

-*
  restart
  needsPackage "PALPInterface"
*-
TEST ///
  ws = {10, 1, 1, 1, 1, 3, 3}  
  needsPackage "StringTorics"
  M = getVerticesFromWS ws
  P = convexHull M
  dim P
  V = reflexiveToSimplicialToricVariety P
  isSimplicial V
  isSmooth V
  transpose matrix rays V
  max V
  assert isWellDefined V
  normalForm transpose matrix rays V

  P2 = polar P
  vertices P2
  normalForm lift(vertices P2, ZZ)
  (LP, tri) = regularStarTriangulation(3, P2)
  V = normalToricVariety(LP, tri)
  isSimplicial V
  isSmooth V -- yes!

  m = lift(matrix vertices P2, ZZ)

  vertsP2 = entries transpose vertices P2
  raysV = rays V
  
  runNEF(m, "N -c2")
  runNEF(matrix rays V, "N -c2")
  -- this gives us the following nef partition of V
  D1 = V_2 + V_3
  D2 = -toricDivisor V - D1
  isNef D1
  isNef D2

  - degree toricDivisor V
  degree D1
  degree D2
  transpose matrix LP

  SR = dual monomialIdeal V
  S1 = ZZ/32003[(gens ring V)]
  minimalBetti sub(SR, S1)
  picardGroup V
  classGroup V

  X = completeIntersection(V, {D1, D2})
  hodgeDiamond X -- h11=3, although induced is only ZZ^2.
  
  -- if  ws = {30, 1, 2, 4, 4, 7, 12}
  -- these don't work so well...
  --elapsedTime HH^1(V, OO_V(0, 0, 1, 0, 1, 0, 1, 0, 0, 0))
  --elapsedTime cohomCalg(V, V_2 + V_4 + V_8)

  elapsedTime cohomCalg(V, V_2 + V_4 + V_5)
///

-*
  restart
  needsPackage "PALPInterface"
*-
TEST ///
  needsPackage "StringTorics"

  ws = {10, 1, 1, 1, 1, 3, 3}

  M = getVerticesFromWS ws
  P = convexHull M
  dim P
  V = reflexiveToSimplicialToricVariety P
  isSimplicial V
  isSmooth V
  transpose matrix rays V
  max V
  assert isWellDefined V
  normalForm transpose matrix rays V

  P2 = polar P
  vertices P2
  normalForm lift(vertices P2, ZZ)
  (LP, tri) = regularStarTriangulation(3, P2)
  V = normalToricVariety(LP, tri)
  isSimplicial V
  isSmooth V -- yes!

  m = lift(matrix vertices P2, ZZ)

  vertsP2 = entries transpose vertices P2
  raysV = rays V
  
  runNEF(m, "N -c2")
  runNEF(matrix rays V, "N -c2")
  -- this gives us the following nef partition of V
  D1 = V_2 + V_3 + V_6
  D2 = -toricDivisor V - D1
  isNef D1
  isNef D2
  X = completeIntersection(V, {D1, D2})
  hodgeDiamond X -- h11=2, same as V.

  posHull transpose matrix degrees ring V
  rays oo
///


///
  wss = getWSFromDim(5, Degrees => (10,10));
  ws = wss_4
  assert(ws == {10, 1, 1, 2, 2, 2, 2})
  
  A1 = getVerticesFromWS ws
  A2 = lift(vertices polar convexHull A1, ZZ)
  normalForm A1
  normalForm A2
  
  myrays = entries transpose A1
  annotatedFaces convexHull A2
  convexHull A2
  isReflexive oo

  needsPackage "StringTorics"
  Q = cyPolytope(A2, ID => 3)
  findAllCYs Q; -- FAILS: dim is too high

  P = convexHull A2
  (verts, tri) = regularStarTriangulation(dim P - 2, P)
  rays Q == verts
  V = normalToricVariety(verts, tri, CoefficientRing => ZZ/101)
  dim V
  assert isSimplicial V
  assert isSmooth V
  assert isProjective V

  dual monomialIdeal V

  H = findNefPartitions(2, V)
  oo#true
  netList oo
  for x in H#true list (
      if #x >= 4 then continue;
      other := sort toList ((set splice{0..#rays V - 1}) - set x);
      if isNef(sum(other, i -> V_i)) then {x, other} else continue
      )
  
  -- I want to take Delta(q), take a triuangulation of it.
  -- Use that.  But I want the one with reasonably small h11...

  matrix for x from -10 to 10 list for y from -4 to 4 list rank HH^0(V, OO_V (x,y))
  matrix for x from -10 to 10 list for y from -4 to 4 list rank HH^1(V, OO_V (x,y))
  matrix for x from -10 to 10 list for y from -4 to 4 list rank HH^2(V, OO_V (x,y))
  matrix for x from -15 to 5 list for y from -10 to 4 list rank HH^3(V, OO_V (x,y))
  matrix for x from -15 to 5 list for y from -10 to 4 list rank HH^4(V, OO_V (x,y))
  matrix for x from -15 to 10 list for y from -10 to 0 list rank HH^5(V, OO_V (x,y))
  ///

end--

-* Development section *-
restart
needsPackage "PALPInterface"
check "PALPInterface"

uninstallPackage "PALPInterface"
restart
installPackage "PALPInterface"
viewHelp "PALPInterface"

restart
debug needsPackage "PALPInterface"
M = matrix{{1,1,1,1},{0,1,2,3}}
toPalp transpose M
getPolyhedralInfo M
getPolyhedralInfo transpose M



"foo" << toPalp M << close
runPoly(M, "v")
normalForm M
viewHelp runProgram

needsPackage "StringTorics"
kss = kreuzerSkarke(7, Limit => 10)

L = vertices polar convexHull matrix kss_7
L = sub(L, ZZ)
L1 = normalForm L

first lines runPoly(L1, "g") -- informational line
runPoly(L1, "p")
runPoly(L1, "v")
runPoly(L1, "e")
runPoly(L1, "m")
runPoly(L1, "d")
runPoly(L1, "a")
runPoly(L1, "D")
runPoly(L1, "i")
runPoly(L1, "I")
runPoly(L1, "S")
runPoly(L1, "Tv")
runPoly(L1, "Tp")
runPoly(L1, "N")
runPoly(L1, "t")

-- Using nef.x
3 1 1 1 0 0 0 0 0  2 0 0 0 1 1 0 0 0  3 0 0 0 0 0 1 1 1
nef.x -h
nef.x -p -c2 -V << FOO
3 1 1 1 0 0 0 0 0  2 0 0 0 1 1 0 0 0  3 0 0 0 0 0 1 1 1
FOO

-- Using cws.x
cws.x -h

10 1 1 1 1 1 5 M:1128 6 N:8 6 H:1,0,976 [5910]

# points are the columns
poly.x -v << FOO
5 1 1 1 1 1
FOO

# points are the rows
poly.x -e << FOO
5 1 1 1 1 1
FOO

# points are the rows
poly.x -v << FOO
45 5 6 7 8 9 10
FOO

run "poly.x -e << FOO
3 1 1 1 0 0 0 3 0 0 0 1 1 1
FOO
"

run "poly.x -e << FOO
36 1 4 4 6 9 12
FOO
"

str = get "!poly.x -v -r << FOO
36 1 4 4 6 9 12
FOO
"

str = get "!poly.x -v -r << FOO
10 1 2 3 4
FOO
"

L = lines str
L = drop(drop(L, 3), -2)
netList L
L0 = separate(" +", L_0)
for x in drop(L0, 1) list value x
M = matrix for ell in L list (
    L0 := separate(" +", ell);
    for x in drop(L0, 1) list value x
    )
needsPackage "Polyhedra"
needsPackage "StringTorics"
P = convexHull M
vertices P
isReflexive P

get "!poly.x -e << FOO
42 2 3 5 5 6 21
FOO
"

get "!poly.x -e << FOO
143233 43 1651 3328 20226 47194 70791
FOO
"

get "!poly.x -v << FOO
143233 43 1651 3328 20226 47194 70791
FOO
"

lines get///!curl "http://rgc.itp.tuwien.ac.at/fourfolds/db/5d_reflexive,h11=100.txt"///;
#oo == 204046
o31_100000





---------- nef-partitions ---------
restart
needsPackage "PALPInterface"
needsPackage "StringTorics"
needsPackage "NormalToricVarieties"
viewHelp NormalToricVarieties

V = smoothFanoToricVariety(3, 4)
rays V
max V
V_0
V_1
V_4
toricDivisor V
S = ring V
describe S
2*V_3 + 10*V_2
degree V_3
degree V_1
degree(V_1 + V_3)
degree (-toricDivisor V)
basis({3,2}, S)
F = random({3,2}, S)
size F
isNef V_0
isNef V_1
isNef V_4

V = smoothFanoToricVariety(5, 6)
rays V
transpose matrix degrees ring V
dim V

nef.x -N -c2 -p << FOO
8 5
-1 0 0 0 0
0 -1 0 0 0
0 0 -1 0 0
0 0 0 -1 0
0 0 0 0 1
0 0 0 0 -1
0 0 0 1 -1
1 1 1 0 -3
FOO

D1 = V_4 + V_6 + V_7
D2 = V_0 + V_1 + V_2 + V_3 + V_5
isNef(V_4 + V_6 + V_7)
degree(V_4 + V_6 + V_7)
isNef(V_0 + V_1 + V_2 + V_3 + V_5)
degree(V_0 + V_1 + V_2 + V_3 + V_5)
for i from 0 to 7 list isNef V_i
for i from 0 to 7 list if i == 5 then continue else isNef(V_i + V_3 +  V_5)

X = completeIntersection(V, {D1, D2})
hodgeDiamond X
pt = base(a,b)
Xa = abstractVariety(X, pt)
IX = intersectionRing Xa
describe IX
basis(1, IX)
h = a * t_6 + b * t_7
integral(h^3)
integral((chern_2 tangentBundle Xa) * h)





---------- nef-partitions ---------
restart
needsPackage "PALPInterface"
needsPackage "StringTorics"
needsPackage "NormalToricVarieties"

testinfo = getPartitionInfo(5, 6, 2) -- remove this line: this grabs from DB in toric varieties or polyhedra...
Y = smoothFanoToricVariety(5, 6)
testPartitionList = readPartitions(testinfo)
-- (testD1, testD2) = formDivisor(smoothFanoToricVariety(5, 6), {4, 6, 7})
divisorsList = for thispartcomb in testPartitionList list (
    (testD1, testD2) = formDivisor(Y, thispartcomb)
)
oo/(x -> (x/isNef))
-- X = completeIntersection(testV, {D1, D2});
X = completeIntersection(Y, toList divisorsList_0);
hodgeDiamond X

--- testing PALP calls with hereis docs in them.
restart
needsPackage "PALPInterface"

M = getVerticesFromWS {10,1,2,3,4}
normalForm M
M2 = getVerticesFromWS {3,1,1,1,0,0,0, 3,0,0,0,1,1,1}
normalForm M2

P = convexHull(M2)
assert isReflexive P
assert(#latticePoints P == 100)



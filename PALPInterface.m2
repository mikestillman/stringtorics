newPackage(
    "PALPInterface",
    Version => "0.2",
    Date => "20 September 2026",
    Headline => "an interface to the PALP package for lattice polytopes",
    Authors => {
        {Name => "Mike Stillman",
            Email => "mike@math.cornell.edu",
            HomePage => "https://mikestillman.github.io"}
        },
    Keywords => {"Interfaces", "Toric Geometry"},
    HomePage => "https://github.com/mikestillman/stringtorics",
    AuxiliaryFiles => false,
    DebuggingMode => true,
    PackageImports => {"Polyhedra", "NormalToricVarieties", "LLLBases"}
    )

export {
    -- locating and running the PALP executables
    "palpProgram",
    "palpRun",
    "toPalp",
    "fromPalpMatrix",

    -- poly.x
    "palpNormalForm",
    "palpMVertices",
    "palpMPoints",
    "palpNVertices",
    "palpNPoints",
    "palpInfo",
    "palpIsReflexive",

    -- weight systems
    "weightSystems",
    "weightSystemVertices",
    "weightSystemPolytope",

    -- nef partitions
    "palpNefPartitions",
    "parsePALPNefPartitions",
    "naiveNefPartitions",
    "nefPartitionDivisors"
    }

exportMutable {
    "palpVerbosity",
    "palpDmax"
    }

palpVerbosity = 0

-------------------------------------------------
-- Locating the PALP executables -----------------
-------------------------------------------------

-- PALP is compiled with a fixed upper bound POLY_Dmax on the dimensions it can
-- handle, and `make all-dims` builds one executable per bound.  A program named
-- e.g. "poly" is therefore available as poly-4d.x, poly-5d.x, poly-6d.x,
-- poly-11d.x, and poly.x (which is the same build as poly-6d.x).
--
-- The requirement is  POLY_Dmax >= dim N + codim - 1,  so codim 1 computations
-- (poly.x, cws.x) are fine with the default build up to dimension 6, while nef
-- partitions run out much sooner: dim 5 with codim 3, or dim 6 with codim 2,
-- already need the 11d build.  We therefore prefer the largest build.
palpDmax = 11

palpDmaxSuffixes = new HashTable from {4 => "-4d", 5 => "-5d", 6 => "", 11 => "-11d"}

palpProgramCache = new MutableHashTable

palpProgram = method(TypicalValue => Program)
palpProgram String := name -> (
    if not palpDmaxSuffixes#?palpDmax then error(
        "palpDmax should be one of " | demark(", ", toString \ sort keys palpDmaxSuffixes));
    if palpProgramCache#?(name, palpDmax) then return palpProgramCache#(name, palpDmax);
    exe := name | palpDmaxSuffixes#palpDmax | ".x";
    P := findProgram(exe, exe | " -h", RaiseError => false);
    if P === null and palpDmaxSuffixes#palpDmax =!= "" then (
        -- the requested variant is not installed; the default build is always 6d
        exe = name | ".x";
        P = findProgram(exe, exe | " -h", RaiseError => false);
        );
    if P === null then error(
        "cannot find the PALP program " | name | ".x" | newline |
        "PALP is at https://hep.itp.tuwien.ac.at/~kreuzer/CY/CYpalp.html, or use" | newline |
        "  brew install macaulay2/tap/palp");
    palpProgramCache#(name, palpDmax) = P
    )

-------------------------------------------------
-- Running a PALP executable ---------------------
-------------------------------------------------

-- true if the string contains no non-whitespace character
isBlankString = s -> not match("[^[:space:]]", s)

-- The compile-time bounds PALP complains about, as they appear in its messages
-- ("Please increase POLY_Dmax to at least 8 ...", "Need SYM_Nmax > 46083 !!").
-- None of these strings occurs in ordinary PALP output.
palpLimits = {"POLY_Dmax", "POINT_Nmax", "VERT_Nmax", "AMBI_Dmax", "SYM_Nmax", "WDIM"}

-- PALP takes its input from a file named on the command line.  We write that
-- file under M2's session temporary directory, which M2 removes on exit, rather
-- than into the current directory.  Set palpVerbosity > 0 to see where it is,
-- or use KeepFiles => true to hang on to it.
runPALP = method(Options => {KeepFiles => false})
runPALP(Program, String, String) := String => opts -> (P, args, input) -> (
    f := temporaryFileName() | ".palp";
    f << input << close;
    if palpVerbosity > 0 then printerr(
        "palp: " | P#"name" | " " | args | " " | f | newline | input);
    result := runProgram(P, args | " " | f, RaiseError => false);
    if opts.KeepFiles
    then printerr("palp: input left in " | f)
    else removeFile f;
    out := result#"output";
    if result#"return value" =!= 0 then error(
        "PALP command failed: " | result#"command" | newline | result#"error");
    if not isBlankString result#"error" then error(
        "PALP wrote to stderr: " | result#"command" | newline | result#"error");
    -- PALP is compiled with fixed bounds on the size of everything it handles,
    -- and reports running out of one of them on stdout while still exiting with
    -- status 0, so this can only be caught by inspecting the output.  POLY_Dmax
    -- bounds the dimension and is the one palpDmax selects; the others are fixed
    -- when PALP is built.  SYM_Nmax, for instance, stops poly.x -N at dimension
    -- eight on a simplex.
    lim := select(palpLimits, v -> match(v, out));
    if #lim > 0 then error(
        P#"name" | " ran out of " | demark(", ", lim) | " on this input:" | newline |
        out |
        (if member("POLY_Dmax", lim)
         then "Set palpDmax to a larger value, or install that variant of PALP."
         else "This is a compile-time bound; PALP must be rebuilt with a larger value."));
    if palpVerbosity > 1 then printerr out;
    out
    )

palpRun = method(Options => {KeepFiles => false})
palpRun(String, String, Matrix) :=
palpRun(String, String, List) := String => opts -> (name, args, input) ->
    runPALP(palpProgram name, args, palpInput input, KeepFiles => opts.KeepFiles)

-------------------------------------------------
-- Translating to and from PALP's format ---------
-------------------------------------------------

-- A matrix becomes "#rows #cols" followed by the rows.  A list of integers is
-- passed through as a (combined) weight system, which PALP distinguishes from a
-- matrix by the shape of its first line.
toPalp = method()
toPalp Matrix := String => M -> (
    if ring M =!= ZZ then error "expected a matrix over ZZ";
    if numRows M === 0 or numColumns M === 0 then error "expected a non-empty matrix";
    concatenate(
        toString numRows M, " ", toString numColumns M, newline,
        between(newline, for row in entries M list demark(" ", toString \ row)),
        newline)
    )
toPalp List := String => wlist -> (
    if #wlist === 0 then error "expected a non-empty weight system";
    if not all(wlist, a -> instance(a, ZZ)) then error "expected a list of integers";
    demark(" ", toString \ wlist) | newline
    )

-- PALP writes a block of numbers preceded by a header line
--     <#lines> <#columns>  <comment>
-- Returns (#lines, #columns, rows), or null if there is no block at all.
palpReadBlock = str -> (
    L := select(lines str, ell -> not isBlankString ell);
    if #L === 0 then return null;
    hdr := select(separate(" +", L#0), x -> x =!= "");
    if #hdr < 2 or not all(take(hdr, 2), x -> match("^[0-9]+$", x))
    then error("expected a PALP matrix header, got: " | L#0);
    nlines := value hdr#0;
    ncols := value hdr#1;
    if #L < nlines + 1 then error(
        "PALP announced " | toString nlines | " rows but printed " | toString(#L - 1));
    rows := for ell in take(L, {1, nlines}) list
        for x in select(separate(" +", ell), y -> y =!= "") list (
            if not match("^-?[0-9]+$", x) then error(
                "expected an integer in PALP output, got: " | x);
            value x);
    (nlines, ncols, rows)
    )

-- The matrix exactly as PALP printed it.  No reorientation: this is the raw
-- reader, and is the inverse of toPalp.
fromPalpMatrix = method()
fromPalpMatrix String := Matrix => str -> (
    b := palpReadBlock str;
    if b === null then error "expected non-empty PALP output";
    (nlines, ncols, rows) := b;
    if not all(rows, r -> #r === ncols) then error(
        "PALP announced " | toString ncols | " columns, but its rows do not agree:"
        | newline | str);
    matrix rows
    )

-- PALP chooses the orientation of its output block per example rather than per
-- option -- poly.x -p prints 23 x 3 for one weight system and 2 x 10 for another
-- -- so we orient every point set the same way, one point per COLUMN.  The
-- dimension is the smaller of the two, since a full dimensional polytope in
-- dimension d has more than d points.
--
-- Returns null when PALP printed nothing (for instance poly.x -d applied to a
-- polytope that is not reflexive), and also when the rows carry an extra
-- inhomogeneous entry, which is how poly.x -e reports the inequalities of a
-- polytope that is not reflexive: those are not dual vertices.
palpPointMatrix = str -> (
    b := palpReadBlock str;
    if b === null then return null;
    (nlines, ncols, rows) := b;
    if not all(rows, r -> #r === ncols) then return null;
    M := matrix rows;
    if nlines <= ncols then M else transpose M
    )

-- PALP reads a matrix as a set of points and takes the smaller of its two
-- dimensions to be the dimension of the lattice, assuming throughout that the
-- points affinely span it.  It does not check.  Given a polytope of lower
-- dimension it will quietly return a wrong answer -- poly.x -v reports three of
-- the four vertices of a square lying in a plane in ZZ^3 -- or die on a signal:
-- poly.x -g on that same square exits with a bus error.  Given fewer points
-- than the dimension it prints "Identical points in Vec_Greater_Than !!".  None
-- of that is worth passing on to the caller, so check first.
palpFullDimensional = M -> (
    n := numColumns M;
    if n <= numRows M then false
    else rank((M_(toList(1 ..< n))) - M_{0} * matrix {toList((n-1) : 1)}) === numRows M
    )

checkPalpInput = method()
checkPalpInput Matrix := M -> (
    if numColumns M <= numRows M then error(
        "PALP needs more points than the dimension, but was given a "
        | toString numRows M | " by " | toString numColumns M | " matrix");
    if not palpFullDimensional M then error(
        "PALP needs the points to span the lattice affinely, but the columns of "
        | "this " | toString numRows M | " by " | toString numColumns M
        | " matrix span only dimension "
        | toString rank((M_(toList(1 ..< numColumns M))) - M_{0} * matrix {toList((numColumns M - 1) : 1)})
        | newline
        | "Change coordinates so that the polytope is full dimensional.");
    )
checkPalpInput List := q -> (
    -- a single weight system {d, q_0, q_1} describes a polytope of dimension one,
    -- on which poly.x fails an assertion in ReadCwsPp and aborts
    if #q === 3 then error(
        "PALP aborts on a weight system of dimension one; give the polytope as a "
        | "matrix instead, for instance matrix{{1,-1}}");
    if #q < 3 then error "expected a weight system {d, q_0, ..., q_n}";
    )

-- everything that hands a polytope to PALP goes through this
palpInput = A -> (checkPalpInput A; toPalp A)

palpPoly = (args, A) -> palpPointMatrix runPALP(palpProgram "poly", args, palpInput A)

-------------------------------------------------
-- The basic poly.x computations -----------------
-------------------------------------------------

-- Each of these accepts either a matrix whose columns span the polytope or a
-- (combined) weight system, and returns a matrix whose columns are the points
-- asked for, or null if PALP produced no answer.

palpNormalForm = method()
palpNormalForm Matrix :=
palpNormalForm List := Matrix => A -> palpPoly("-N", A)

palpMVertices = method()
palpMVertices Matrix :=
palpMVertices List := Matrix => A -> palpPoly("-v", A)

palpMPoints = method()
palpMPoints Matrix :=
palpMPoints List := Matrix => A -> palpPoly("-p", A)

palpNVertices = method()
palpNVertices Matrix :=
palpNVertices List := Matrix => A -> palpPoly("-e", A)

palpNPoints = method()
palpNPoints Matrix :=
palpNPoints List := Matrix => A -> palpPoly("-d", A)

-- poly.x -g prints one line of "Key:value ..." fields, for example
--     M:7 6 N:27 8 Pic:17 Cor:0
--     10 1 1 2 2 2 2 M:378 6 N:8 6 H:2,0,350 [2160]
-- where M and N give (#points, #vertices), H the Hodge numbers, and the value in
-- brackets the Euler characteristic.
palpParseInfo = ell -> (
    H := new MutableHashTable;
    key := null;
    for t in select(separate(" +", ell), x -> x =!= "") do (
        if match("^\\[-?[0-9]+\\]$", t) then (
            H#"Euler" = value substring(t, 1, #t - 2);
            key = null;
            )
        else if match(":", t) then (
            i := first first regex(":", t);
            key = substring(t, 0, i);
            rest := substring(t, i + 1);
            H#key = if rest === "" then {} else for x in separate(",", rest) list value x;
            )
        else if key =!= null and match("^-?[0-9]+$", t)
        then H#key = append(H#key, value t)
        else key = null;
        );
    new HashTable from H
    )

palpInfo = method()
palpInfo Matrix :=
palpInfo List := HashTable => A -> (
    out := runPALP(palpProgram "poly", "-g", palpInput A);
    if isBlankString out then new HashTable from {} else palpParseInfo first lines out
    )

-- Two PALP options answer this.  The summary line of poly.x -g carries an "N:"
-- field for a reflexive polytope and an "F:" field for one that is not, but -g
-- counts the lattice points of the dual and so becomes expensive as the
-- dimension grows: on the simplex in dimension eight it takes three seconds
-- against twenty milliseconds for -e.  So we use -e, which returns the vertices
-- of the dual when there are any and inequalities, which palpPointMatrix
-- reports as null, when there are not.  The two agree wherever both are
-- affordable.
palpIsReflexive = method()
palpIsReflexive Matrix :=
palpIsReflexive List := Boolean => A -> palpNVertices A =!= null

-------------------------------------------------
-- Weight systems --------------------------------
-------------------------------------------------

-- A weight system is written {d, q_0, ..., q_n} with d the degree.  The ones of
-- interest here are those with d = sum of the q_i, which are exactly the ones
-- whose degree d hypersurface in the weighted projective space P(q) has trivial
-- canonical class.  cws.x enumerates them.
weightSystems = method(Options => {Degrees => null})
weightSystems ZZ := List => opts -> d -> (
    if d < 1 then error "expected a positive dimension";
    args := "-w" | toString d;
    drange := opts.Degrees;
    if drange =!= null then (
        if instance(drange, ZZ) then drange = {drange, drange};
        if not instance(drange, BasicList) or #drange =!= 2
        then error "expected Degrees => {lo, hi} or Degrees => deg";
        if not all(drange, a -> instance(a, ZZ)) then error "expected Degrees to be integers";
        args = args | " " | toString drange#0 | " " | toString drange#1;
        );
    -- cws.x -w takes no input file, so we call it directly rather than via runPALP
    result := runProgram(palpProgram "cws", args, RaiseError => false);
    if result#"return value" =!= 0 then error(
        "cws.x failed: " | result#"command" | newline | result#"error");
    if match("POLY_Dmax", result#"output") then error(
        "cws.x was built with too small a POLY_Dmax for dimension " | toString d |
        ":" | newline | result#"output");
    -- Each line is "<d> <q_0> ... <q_n>" followed by flags such as "rt".  When a
    -- degree range is given, cws.x adds a trailing summary line beginning with #.
    for ell in lines result#"output" list (
        if isBlankString ell or match("^#", ell) then continue;
        toks := take(select(separate(" +", ell), x -> x =!= ""), d + 2);
        if #toks =!= d + 2 or not all(toks, x -> match("^[0-9]+$", x))
        then error("unexpected cws.x output line: " | ell);
        for x in toks list value x
        )
    )

-- The vertices of Delta(q).  This is palpMVertices on a weight system; it is
-- named separately because that is how one looks for it, and because a combined
-- weight system is also accepted here.
weightSystemVertices = method()
weightSystemVertices List := Matrix => q -> palpMVertices q

-- Delta(q) = {x : q.x = 0, x_i >= -1}, computed from the definition rather than
-- by calling PALP, and returned in a basis of the lattice ker(q) \cap ZZ^n so
-- that it is full dimensional.
--
-- Note that Delta(q) itself is only a rational polytope: the vertex obtained by
-- setting every coordinate but the j-th to -1 is (d - q_j)/q_j, which is an
-- integer only when q_j divides d.  The lattice polytope that PALP reports, and
-- that this function returns, is the convex hull of its lattice points.
weightSystemPolytope = method()
weightSystemPolytope List := Polyhedron => wts -> (
    if #wts < 3 then error "expected a weight system {d, q_0, ..., q_n}";
    q := drop(wts, 1);
    if not all(q, a -> instance(a, ZZ) and a > 0)
    then error "expected the weights to be positive integers";
    if wts#0 =!= sum q then error(
        "expected the degree " | toString wts#0 | " to be the sum of the weights " |
        toString sum q | newline |
        "(a combined weight system is not handled here; use weightSystemVertices)");
    n := #q;
    B := transpose LLL syz matrix {q};    -- rows: a ZZ-basis of ker(q) \cap ZZ^n
    Delta := polyhedronFromHData(- transpose B, matrix toList(n : {1}));
    convexHull matrix {latticePoints Delta}
    )

-------------------------------------------------
-- Nef partitions --------------------------------
-------------------------------------------------

-- nef.x reports one nef partition per line, beginning with H:, for example
--     H:2 86 [-168] P:0 V:1 5  6   (2 8) (1 4)     0sec  0cpu
--     H:20 [24] P:0 V0:1 5  6  V1:3 4   (2 4 4) (1 2 2)    1sec  0cpu
-- The V fields list the parts of the partition, as indices into the list of
-- points that -Lp prints.  In codimension c there are c-1 of them; the points
-- not mentioned form the remaining part.
palpParsePartitionLine = ell -> (
    groups := {};
    collecting := false;
    for t in select(separate(" +", ell), x -> x =!= "") do (
        if match("^V[0-9]*:", t) then (
            rest := substring(t, 1 + first first regex(":", t));
            groups = append(groups, if rest === "" then {} else {value rest});
            collecting = true;
            )
        else if collecting and match("^[0-9]+$", t)
        then groups = append(drop(groups, -1), append(last groups, value t))
        else collecting = false;
        );
    groups
    )

parsePALPNefPartitions = method()
parsePALPNefPartitions String := List => out -> (
    for ell in lines out list (
        if not match("^H:", ell) then continue;
        palpParsePartitionLine ell
        )
    )

-- the block of N-lattice points that nef.x -Lp prints, one point per column
palpNefPoints = out -> (
    L := lines out;
    i := position(L, ell -> match("Points of Poly in N-Lattice", ell));
    if i === null then null
    else palpPointMatrix concatenate between(newline, drop(L, i))
    )

-- Returns a hash table with
--   "Points"     the N-lattice points, one per column, or null if PALP said nothing
--   "Partitions" one entry per nef partition, each a list of the parts that PALP
--                printed, as column indices into "Points"
-- An empty list of partitions is an answer, not a failure: it says this polytope
-- has no nef partition of that codimension.
palpNefPartitions = method()
palpNefPartitions(Matrix, ZZ) :=
palpNefPartitions(List, ZZ) := HashTable => (A, c) -> (
    if c < 1 then error "expected a positive codimension";
    out := runPALP(palpProgram "nef", "-Lp -c" | toString c, palpInput A);
    new HashTable from {
        "Points" => palpNefPoints out,
        "Partitions" => parsePALPNefPartitions out
        }
    )

-------------------------------------------------
-- Nef partitions without PALP --------------------
-------------------------------------------------

-- A nef partition of codimension c is a partition of the rays of V into c parts,
-- each of whose sum of divisors is nef.  This searches for them directly, which
-- is a useful check on what PALP reports, but it looks at every set partition of
-- the rays into c parts and so is only usable for small numbers of rays.
naiveNefPartitions = method()
naiveNefPartitions(NormalToricVariety, ZZ) := List => (V, cod) -> (
    if cod < 1 then error "expected a positive number of parts";
    n := #rays V;
    if cod > n then return {};
    isNefPart := memoize(b -> isNef sum(b, i -> V_i));
    result := new MutableList;
    -- restricted growth strings visit each set partition of {0, ..., n-1} once
    recurse := null;
    recurse = (i, blocks) -> (
        if #blocks + (n - i) < cod then return;     -- too few elements left
        if i === n then (
            if #blocks === cod and all(blocks, isNefPart)
            then result#(#result) = blocks;
            return;
            );
        for b from 0 to min(#blocks, cod - 1) do recurse(i + 1,
            if b < #blocks then replace(b, append(blocks#b, i), blocks)
            else append(blocks, {i}));
        );
    recurse(0, {});
    toList result
    )

-- The divisors of a nef partition.  The parts may be given as a list of lists of
-- ray indices, or, for codimension two, as the single list of indices making up
-- one part.  Any rays left over form one more part, so this can be applied
-- directly to what palpNefPartitions reports, which names only c-1 of the c
-- parts.
nefPartitionDivisors = method()
nefPartitionDivisors(NormalToricVariety, List) := List => (V, part) -> (
    if #part === 0 then error "expected a nonempty partition";
    n := #rays V;
    parts := if all(part, x -> instance(x, ZZ)) then {toList part}
        else for p in part list (
            if not instance(p, VisibleList)
            then error "expected a list of ray indices, or a list of lists of them";
            toList p);
    if any(parts, p -> #p === 0) then error "expected the parts to be nonempty";
    all0n := flatten parts;
    if not all(all0n, i -> instance(i, ZZ) and 0 <= i and i < n)
    then error("expected ray indices in the range 0 .. " | toString(n - 1));
    if #unique all0n =!= #all0n then error "expected the parts to be disjoint";
    rest := sort toList(set toList(0 .. n - 1) - set all0n);
    if #rest > 0 then parts = append(parts, rest);
    for p in parts list sum(p, i -> V_i)
    )

beginDocumentation()

doc ///
Key
  PALPInterface
Headline
  an interface to the PALP package for lattice polytopes
Description
  Text
    PALP, a Package for Analyzing Lattice Polytopes, was written by Maximilian
    Kreuzer and Harald Skarke in order to classify reflexive polytopes, and in
    particular to produce the list of the 473800776 reflexive polytopes of
    dimension four.  Along the way it acquired a good deal of functionality that
    is useful on its own: normal forms of lattice polytopes, the lattice points
    of a polytope and of its dual, the enumeration of weight systems, and nef
    partitions.  This package makes that functionality available inside
    Macaulay2.
  Text
    The program, its source and its documentation are at
    @HREF "https://hep.itp.tuwien.ac.at/~kreuzer/CY/CYpalp.html"@; the reference
    documentation is reached from there, under "Former PALP wiki" and "Last
    version of the PALP online documentation".  PALP is described in
  Pre
    M. Kreuzer and H. Skarke, PALP: A Package for Analyzing Lattice Polytopes
    with Applications to Toric Geometry, Comput. Phys. Commun. 157 (2004)
    87-106, arXiv:math/0204356.
  Text
    PALP is a separate program, and must be installed before any of these
    functions will run.  On a machine already set up for Macaulay2, use
  Pre
    brew install macaulay2/tap/palp
  Text
    The package itself loads whether or not PALP is present; the executables are
    looked for only when a function needs one.  If they are installed somewhere
    that is not on your path, add that directory to @TO "programPaths"@.
  Text
    A first computation: the polytope of a weight system, its vertices, and its
    normal form.
  Example
    ws = {10, 1, 2, 3, 4}
    V = weightSystemVertices ws
    palpNormalForm V
      P = convexHull V
    isReflexive P
  Text
    The columns of the matrix are the vertices.  PALP itself prints a point set
    sometimes as rows and sometimes as columns, depending on the example, but
    every function here returns one point per column.
  Text
    PALP is compiled with a fixed bound on the dimensions it can handle, and
    ships as several executables with different bounds.  See @TO "palpDmax"@,
    which matters mainly for @TO palpNefPartitions@.
  Text
    PALP assumes its input is a full dimensional polytope and does not check, so
    the functions here check before calling it.  A set of points spanning
    something of lower dimension makes PALP give a wrong answer under some
    options and crash under others, and a weight system of dimension one makes
    it abort; rather than pass any of that on, these are refused with an error.
    The file @TT "PALP-upstream-notes.md"@ in this package's repository has the
    details.
  Example
    M = matrix{{1,-1,-1,1},{1,1,-1,-1},{0,0,0,0}};   -- a square in a plane in ZZ^3
    assert(try (palpMVertices M; false) else true)
SeeAlso
  palpNormalForm
  weightSystems
  palpNefPartitions
  "palpDmax"
  "ReflexivePolytopesDB::ReflexivePolytopesDB"
///

doc ///
Key
  "palpDmax"
Headline
  which build of the PALP executables to use
Description
  Text
    PALP is compiled with a fixed upper bound @TT "POLY_Dmax"@ on the dimensions
    it can handle, and its @TT "make all-dims"@ target builds one executable per
    bound: @TT "poly-4d.x"@, @TT "poly-5d.x"@, @TT "poly-6d.x"@,
    @TT "poly-11d.x"@, and @TT "poly.x"@, which is the same build as the 6d one.
    This variable says which of them to prefer, and may be set to 4, 5, 6 or 11.
  Example
    palpDmax
  Text
    The requirement PALP imposes is
  Pre
    POLY_Dmax >= dim N + codim - 1
  Text
    so the default build is enough for @TO palpMVertices@ and friends up to
    dimension six, while nef partitions run out much sooner: a five dimensional
    polytope in codimension three, or a six dimensional one in codimension two,
    already needs the 11d build.  That is the main reason to have it.
  Example
    ws = {10, 1, 1, 2, 2, 2, 2};
    #(palpNefPartitions(ws, 3))#"Partitions"
  Text
    Asking for a build that is not installed falls back on the default one.  If
    the result is a PALP that cannot handle the input, that is reported as an
    error rather than being silently mistaken for an empty answer, since PALP
    says so on its standard output and still exits successfully.
Caveat
  Changing this value selects a different executable, so any program already
  located is looked up again.
SeeAlso
  palpProgram
  palpNefPartitions
  "palpVerbosity"
///

doc ///
Key
  "palpVerbosity"
Headline
  how much to report about the PALP commands being run
Description
  Text
    At 0, the default, nothing is printed.  At 1, each PALP command is printed
    together with the temporary file holding its input.  At 2, PALP's output is
    printed as well.
  Example
    palpVerbosity = 1;
    palpNormalForm {3, 1, 1, 1}
    palpVerbosity = 0;
SeeAlso
  palpRun
  "palpDmax"
///

doc ///
Key
  palpProgram
  (palpProgram, String)
Headline
  locate one of the PALP executables
Usage
  P = palpProgram name
Inputs
  name:String
    the base name of a PALP program, such as "poly", "cws" or "nef", with no
    dimension suffix and no ".x"
Outputs
  P:Program
    as returned by @TO findProgram@
Description
  Text
    The build selected by @TO "palpDmax"@ is preferred, falling back on the
    default one if it is not installed.  The result is remembered, so the search
    happens once per program and per value of @TO "palpDmax"@.
  Example
    P = palpProgram "poly"
    P#"name"
  Text
    An error is raised, naming where to get PALP, if the program cannot be found
    at all.  Nothing is looked up when this package is loaded, so the package is
    usable without PALP right up to the point where a computation needs it.
SeeAlso
  palpRun
  "palpDmax"
  "programPaths"
///

doc ///
Key
  palpRun
  (palpRun, String, String, Matrix)
  (palpRun, String, String, List)
  [palpRun, KeepFiles]
Headline
  run a PALP program and return its output
Usage
  s = palpRun(name, args, A)
Inputs
  name:String
    the base name of a PALP program, as for @TO palpProgram@
  args:String
    the command line options, for instance "-v" or "-Lp -c2"
  A:Matrix
    over the integers, whose columns span the polytope, or a @TO List@ holding a
    (combined) weight system
  KeepFiles => Boolean
    whether to keep the temporary file holding PALP's input, and print its name
Outputs
  s:String
    everything PALP wrote to its standard output
Description
  Text
    This is the general way in to PALP, for the options that have no wrapper of
    their own.
  Example
    M = transpose matrix {{1,0,0},{-1,0,1},{0,1,0},{0,-1,1},{0,0,1},{0,0,-1}};
    print palpRun("poly", "-g", M)
    print palpRun("poly", "-v", {3, 1, 1, 1})
  Text
    The input is written to a file under Macaulay2's temporary directory, which
    is removed when Macaulay2 exits, and nothing is written to the current
    directory.  Set @TO "palpVerbosity"@ to 1 to see the command and the name of
    that file.
  Text
    PALP reports some failures on its standard output while still exiting
    successfully; those are turned into errors here rather than passed along as
    if they were results.
SeeAlso
  palpProgram
  "palpVerbosity"
  fromPalpMatrix
///

doc ///
Key
  toPalp
  (toPalp, Matrix)
  (toPalp, List)
Headline
  the PALP input format for a matrix or a weight system
Usage
  s = toPalp A
Inputs
  A:Matrix
    over the integers, or a @TO List@ of integers holding a (combined) weight
    system
Outputs
  s:String
    in the format PALP reads
Description
  Text
    A matrix becomes its number of rows and columns followed by its entries.
    PALP reads that as a set of points, taking the smaller of the two numbers to
    be the dimension.
  Example
    M = matrix{{1,1,1,1},{0,1,2,3}}
    print toPalp M
  Text
    A list of integers is passed through as a single line, which is how PALP
    recognises a weight system.  The first entry is the degree and the rest are
    the weights; several such blocks in a row make a combined weight system.
  Example
    print toPalp {10, 1, 2, 3, 4}
    print toPalp {3,1,1,1,0,0,0, 3,0,0,0,1,1,1}
SeeAlso
  fromPalpMatrix
  palpRun
///

doc ///
Key
  fromPalpMatrix
  (fromPalpMatrix, String)
Headline
  read a matrix printed by PALP
Usage
  M = fromPalpMatrix s
Inputs
  s:String
    output of a PALP program beginning with a header line
Outputs
  M:Matrix
    over the integers, exactly as PALP printed it
Description
  Text
    This is the inverse of @TO toPalp@ on matrices, and does not change the
    orientation of what it reads.
  Example
    M = matrix{{1,1,1,1},{0,1,2,3}}
    assert(M == fromPalpMatrix toPalp M)
  Example
    print palpRun("poly", "-v", {3, 1, 1, 1})
    fromPalpMatrix palpRun("poly", "-v", {3, 1, 1, 1})
  Text
    Since PALP does not always print a point set the same way round, the
    functions that return point sets, such as @TO palpMVertices@, reorient what
    they read so that each point is a column.  This function does not: it is the
    raw reader.
Caveat
  An error is raised if PALP printed nothing, or if the rows do not have the
  length that PALP's header announced.  The latter happens for the inequalities
  that @TT "poly.x -e"@ prints for a polytope that is not reflexive, which carry
  an extra inhomogeneous entry.
SeeAlso
  toPalp
  palpMVertices
///

doc ///
Key
  palpNormalForm
  (palpNormalForm, Matrix)
  (palpNormalForm, List)
Headline
  the PALP normal form of a lattice polytope
Usage
  N = palpNormalForm A
Inputs
  A:Matrix
    over the integers, whose columns span the polytope, or a @TO List@ holding a
    (combined) weight system
Outputs
  N:Matrix
    whose columns are the vertices in normal form
Description
  Text
    Two lattice polytopes differ by a change of basis of the lattice exactly when
    they have the same normal form, so this is what to compare in order to decide
    whether two descriptions give the same polytope.
  Example
    M = palpMVertices {10, 1, 2, 3, 4}
    palpNormalForm M
    A = matrix {{1, -2, 1}, {2, 1, -2}, {-1, 1, 0}}
    det A
    assert(palpNormalForm (A * M) == palpNormalForm M)
  Text
    Two weight systems that give the same polytope:
  Example
    palpNormalForm {3402, 40, 41, 486, 1134, 1701}
    palpNormalForm {3486, 41, 42, 498, 1162, 1743}
    assert(oo == ooo)
Caveat
  The result is @TO null@ if PALP printed nothing, which happens for input it
  declines to treat as a polytope, such as a single point.
SeeAlso
  palpMVertices
  weightSystemPolytope
///

doc ///
Key
  palpMVertices
  (palpMVertices, Matrix)
  (palpMVertices, List)
  palpMPoints
  (palpMPoints, Matrix)
  (palpMPoints, List)
  palpNVertices
  (palpNVertices, Matrix)
  (palpNVertices, List)
  palpNPoints
  (palpNPoints, Matrix)
  (palpNPoints, List)
Headline
  vertices and lattice points of a polytope and of its dual
Usage
  palpMVertices A
  palpMPoints A
  palpNVertices A
  palpNPoints A
Inputs
  A:Matrix
    over the integers, whose columns span the polytope, or a @TO List@ holding a
    (combined) weight system
Outputs
  :Matrix
    over the integers, one point per column
Description
  Text
    For a polytope $P$ in the lattice $M$, with dual $P^*$ in $N$, these give
    the vertices of $P$, the lattice points of $P$, the vertices of $P^*$, and
    the lattice points of $P^*$ respectively.
  Example
    ws = {10, 1, 2, 3, 4}
    palpMVertices ws
    palpNVertices ws
    numColumns palpMPoints ws
    numColumns palpNPoints ws
  Text
    These agree with @TO "Polyhedra::Polyhedra"@.
  Example
      P = convexHull palpMVertices ws;
    assert(numColumns palpMPoints ws == #latticePoints P)
    assert(set entries transpose palpNVertices ws
        === set entries transpose lift(vertices polar P, ZZ))
  Text
    Each point is a column, whichever way round PALP chose to print it.  PALP
    varies that from example to example rather than from option to option: for
    the weight system above it prints the lattice points as 23 rows of 3, and
    for the next one as 2 rows of 10.
  Example
    numRows palpMPoints {3, 1, 1, 1}
    numColumns palpMPoints {3, 1, 1, 1}
Caveat
  The result is @TO null@ when PALP has no answer to give.  This is the usual
  outcome of asking about the dual of a polytope that is not reflexive: there
  @TT "poly.x -d"@ prints nothing at all, and @TT "poly.x -e"@ prints
  inequalities rather than vertices.

  PALP requires the polytope to be full dimensional, and so requires more points
  than the dimension.  It assumes both and checks neither: given a polytope of
  lower dimension, @TT "poly.x -v"@ quietly returns too few vertices while
  @TT "poly.x -g"@ exits on a bus error, and given too few points it prints
  @TT "Identical points in Vec_Greater_Than !!"@.  These functions test the
  input first and raise an error rather than pass any of that on.

  @TT "poly.x -p"@, and so @TO palpInfo@ with it, counts the lattice points of
  the dual, which grows quickly: on the simplex it takes about twenty
  milliseconds in dimension six, three seconds in dimension eight, and more than
  a minute in dimension nine, while @TT "-v"@ and @TT "-e"@ stay fast.  PALP is
  chiefly used in low dimensions and is not much exercised above seven.
Description
  Text
    For example:
  Example
    M = transpose matrix{{1,-1,-1}, {1,-1,1}, {1,1,-1}, {1,1,2}, {-1,0,0}};
    isReflexive convexHull M
    palpMVertices M
    palpNPoints M
SeeAlso
  palpNormalForm
  palpInfo
  weightSystemVertices
///

doc ///
Key
  palpInfo
  (palpInfo, Matrix)
  (palpInfo, List)
Headline
  the summary line PALP prints for a polytope
Usage
  H = palpInfo A
Inputs
  A:Matrix
    over the integers, whose columns span the polytope, or a @TO List@ holding a
    (combined) weight system
Outputs
  H:HashTable
Description
  Text
    This is @TT "poly.x -g"@, whose one line of output collects the number of
    lattice points and vertices of the polytope and of its dual, and, when PALP
    knows them, the Hodge numbers and the Euler characteristic.
  Example
    palpInfo {10, 1, 1, 2, 2, 2, 2}
  Text
    Each field printed by PALP as @TT "Key:"@ becomes the list of numbers
    following it, and the value in square brackets, the Euler characteristic,
    becomes an integer under "Euler".  The fields are
  Text
    @UL {
      {TT "M", " -- the number of lattice points and the number of vertices of the polytope"},
      {TT "N", " -- the same two numbers for its dual, present only when the polytope is reflexive"},
      {TT "F", " -- the numbers of faces, printed in place of ", TT "N", " when the polytope is not reflexive"},
      {TT "H", " -- Hodge numbers"},
      {TT "Euler", " -- the Euler characteristic, which PALP prints in square brackets"},
      {TT "Pic", " -- for a three dimensional polytope, the Picard number of the generic K3 hypersurface"},
      {TT "Cor", " -- the correction term in the formula for ", TT "Pic", "; see below"}
      }@
  Text
    Which of them appear depends on the dimension and on whether the polytope is
    reflexive, because PALP prints whatever it computes for that case.  For the
    reflexive simplices in dimensions 1 through 8 one gets
  Pre
    dim 1, 2      M N
    dim 3         M N Pic Cor
    dim 4, 5      M N H Euler
    dim 6, 7, 8   M N H
  Text
    and a polytope that is not reflexive gives @TT "M"@ and @TT "F"@ in any
    dimension.  So do not assume a field is there; test with @TT "H#?\"H\""@
    before reading it.  PALP's own summary of what @TT "-g"@ prints is
  Pre
    P reflexive:     numbers of (dual) points/vertices, Hodge numbers
    P not reflexive: numbers of points, vertices, equations
  Text
    so in the second case @TT "F"@ is the number of equations, standing where
    @TT "N"@ would be.
  Text
    In dimension three the anticanonical hypersurface is a K3 surface, and PALP
    reports the Picard number of the generic member together with the term that
    corrects the naive count.  Writing $\ell$ for the number of lattice points
    of a face and $\ell^*$ for the number in its relative interior, and
    $\theta$ for the face of $\Delta$ dual to a face $\theta^*$ of
    $\Delta^*$, Batyrev's formula reads
  Pre
    Pic = l(Delta*) - 4 - sum over facets theta* of l*(theta*)
                        + sum over edges  theta* of l*(theta*) l*(theta)
  Text
    and @TT "Cor"@ is that second sum.  It is nonzero exactly when some edge of
    $\Delta^*$ and the edge of $\Delta$ dual to it both have interior lattice
    points, which is the part of the answer the point counts @TT "M"@ and
    @TT "N"@ cannot show.
  Example
    H = palpInfo transpose matrix {{1,0,0},{-1,0,1},{0,1,0},{0,-1,1},{0,0,1},{0,0,-1}}
    H#"M"
    H#"Pic"
  Text
    A polytope that is canonical but not reflexive, so that PALP reports
    @TT "F"@ and no @TT "N"@:
  Example
    palpInfo matrix{{1,0,0,-1},{0,1,0,-1},{0,0,1,-2}}
SeeAlso
  palpIsReflexive
  palpMVertices
  palpRun
///

doc ///
Key
  palpIsReflexive
  (palpIsReflexive, Matrix)
  (palpIsReflexive, List)
Headline
  whether PALP considers a polytope reflexive
Usage
  palpIsReflexive A
Inputs
  A:Matrix
    over the integers, whose columns span the polytope, or a @TO List@ holding a
    (combined) weight system
Outputs
  :Boolean
Description
  Text
    This asks PALP for the vertices of the dual, with @TT "poly.x -e"@: a
    reflexive polytope has them, and a polytope that is not reflexive has PALP
    return inequalities instead, which @TO palpNVertices@ reports as
    @TO null@.
  Text
    The summary line of @TT "poly.x -g"@ answers the same question, carrying an
    @TT "N:"@ field for a reflexive polytope and an @TT "F:"@ field for one that
    is not, and @TO palpInfo@ will show you that.  It is the slower route,
    because it counts the lattice points of the dual: in dimension eight it
    takes some seconds against some milliseconds for @TT "-e"@.
  Example
    palpIsReflexive {10, 1, 2, 3, 4}
    palpIsReflexive matrix{{1,0,0,-1},{0,1,0,-1},{0,0,1,-1}}
  Text
    The second polytope below is canonical -- its only interior lattice point is
    the origin -- but it is not reflexive.  It is the fan polytope of the
    weighted projective space $P(1,1,1,2)$, and a simplex of this shape is
    reflexive exactly when every weight divides the sum of the weights.
  Example
    M = matrix{{1,0,0,-1},{0,1,0,-1},{0,0,1,-2}}
    palpIsReflexive M
    isReflexive convexHull M
    interiorLatticePoints convexHull M
Caveat
  This agrees with @TO "Polyhedra::isReflexive"@, but goes through PALP rather
  than computing the dual in Macaulay2.
SeeAlso
  palpInfo
  palpNPoints
///

doc ///
Key
  weightSystems
  (weightSystems, ZZ)
  [weightSystems, Degrees]
Headline
  the weight systems of a given dimension
Usage
  L = weightSystems d
Inputs
  d:ZZ
    the dimension of the resulting polytopes
  Degrees => ZZ
    one degree, or a pair giving a range of degrees; the default @TO null@ asks
    for all of them
Outputs
  L:List
    of lists of integers, each of the form {degree, weights}
Description
  Text
    This is @TT "cws.x -w"@, which enumerates the weight systems whose degree is
    the sum of their weights.  Those are exactly the ones whose hypersurface in
    the weighted projective space has trivial canonical class.
  Example
    L = weightSystems 3;
    #L
    L#0
    L#90
    assert all(L, ws -> ws#0 == sum drop(ws, 1))
  Text
    Restricting the degree:
  Example
    weightSystems(5, Degrees => 10)
  Example
    #weightSystems(5, Degrees => {20, 20})
SeeAlso
  weightSystemVertices
  weightSystemPolytope
///

doc ///
Key
  weightSystemVertices
  (weightSystemVertices, List)
Headline
  the vertices of the polytope of a weight system
Usage
  V = weightSystemVertices q
Inputs
  q:List
    of integers {d, q_0, ..., q_n}, the degree followed by the weights, or
    several such blocks in a row making a combined weight system
Outputs
  V:Matrix
    over the integers, the vertices of the polytope as its columns
Description
  Text
    This is @TO palpMVertices@ applied to a weight system, under the name one
    would look for it by.
  Example
    ws = {3402, 40, 41, 486, 1134, 1701}
    V = weightSystemVertices ws
      P = convexHull V
    assert isReflexive P
    assert(dim P == 4)
  Text
    Combined weight systems are accepted as well, which is why this goes through
    PALP rather than through @TO weightSystemPolytope@.
  Example
    weightSystemVertices {3,1,1,1,0,0,0, 3,0,0,0,1,1,1}
Caveat
  The result is @TO null@ when PALP declines the weight system, for instance
  because the corresponding polytope has no interior lattice point.

  A weight system with just two weights describes a polytope of dimension one,
  on which PALP fails an assertion in @TT "ReadCwsPp"@ and aborts, so that case
  is refused here.  Give the polytope as a matrix instead: the one dimensional
  reflexive polytope is @TT "matrix{{1,-1}}"@.
SeeAlso
  weightSystems
  weightSystemPolytope
  palpMVertices
///

doc ///
Key
  weightSystemPolytope
  (weightSystemPolytope, List)
Headline
  the polytope of a weight system, computed from the definition
Usage
  P = weightSystemPolytope q
Inputs
  q:List
    of integers {d, q_0, ..., q_n}, where d is the sum of the weights
Outputs
  P:Polyhedron
Description
  Text
    For a weight system $q$ of degree $d = \sum q_i$, this is
    $\Delta(q) = \{x : \sum q_i x_i = 0, \ x_i \ge -1\}$, returned in a basis of
    the lattice $\ker(q) \cap \mathbb{Z}^n$ so that it is full dimensional.
  Example
    P = weightSystemPolytope {10, 1, 2, 3, 4}
    dim P
      isReflexive P
    vertices P
  Text
    Nothing here calls PALP, so this is a way of checking PALP's answer against
    the definition.  The two agree:
  Example
    ws = {30, 4, 4, 6, 7, 9};
    assert(palpNormalForm lift(vertices weightSystemPolytope ws, ZZ)
        == palpNormalForm weightSystemVertices ws)
  Text
    Note that $\Delta(q)$ is in general only a rational polytope: setting every
    coordinate but the $j$-th to $-1$ gives the vertex $(d - q_j)/q_j$, which is
    an integer only when $q_j$ divides $d$.  What is returned, and what PALP
    reports, is the convex hull of its lattice points.
  Example
    ws = {10, 1, 2, 3, 4};
    numColumns vertices weightSystemPolytope ws
Caveat
  Only a single weight system whose degree is the sum of its weights is
  accepted; for a combined weight system use @TO weightSystemVertices@.
SeeAlso
  weightSystemVertices
  weightSystems
///

doc ///
Key
  palpNefPartitions
  (palpNefPartitions, Matrix, ZZ)
  (palpNefPartitions, List, ZZ)
Headline
  the nef partitions of a reflexive polytope
Usage
  H = palpNefPartitions(A, c)
Inputs
  A:Matrix
    over the integers, whose columns span the polytope, or a @TO List@ holding a
    (combined) weight system
  c:ZZ
    the codimension, that is, the number of parts
Outputs
  H:HashTable
    with keys "Points" and "Partitions"
Description
  Text
    This is @TT "nef.x -Lp"@.  The value under "Points" is the matrix of
    N-lattice points that PALP printed, one per column, and the value under
    "Partitions" is one entry per nef partition, each a list of the parts that
    PALP named, given as column indices into that matrix.
  Example
    ws = {10, 1, 1, 2, 2, 2, 2};
    H = palpNefPartitions(ws, 2)
    H#"Partitions"
    H#"Points"
  Text
    In codimension $c$, PALP names $c-1$ of the $c$ parts and leaves the
    remaining points to make up the last one.  So in codimension two a single
    part is printed, and in codimension three, two parts:
  Example
    (palpNefPartitions(ws, 3))#"Partitions"
  Text
    That codimension three computation needs a PALP built for dimension at least
    @TT "dim N + codim - 1"@, which here is 7; see @TO "palpDmax"@.
  Text
    A polytope with no nef partition of the given codimension gives an empty
    list, which is an answer rather than a failure.
  Example
    (palpNefPartitions({14, 1, 2, 2, 2, 2, 2, 3}, 3))#"Partitions"
SeeAlso
  parsePALPNefPartitions
  naiveNefPartitions
  nefPartitionDivisors
  "palpDmax"
///

doc ///
Key
  parsePALPNefPartitions
  (parsePALPNefPartitions, String)
Headline
  read the nef partitions out of nef.x output
Usage
  L = parsePALPNefPartitions s
Inputs
  s:String
    the output of a run of @TT "nef.x"@
Outputs
  L:List
    one entry per nef partition, each a list of the parts PALP named
Description
  Text
    @TO palpNefPartitions@ calls this; it is exported for use on output obtained
    some other way, for instance through @TO palpRun@ with different options.
  Example
    s = palpRun("nef", "-Lp -c2", {10, 1, 1, 2, 2, 2, 2});
    parsePALPNefPartitions s
  Text
    Each nef partition is one line of PALP's output beginning with @TT "H:"@,
    and the parts are its @TT "V:"@ or @TT "V0:"@, @TT "V1:"@, ... fields.  Any
    other line is skipped.
  Example
    parsePALPNefPartitions "H:20 [24] P:0 V0:1 5  6  V1:3 4   (2 4 4) (1 2 2)\n"
SeeAlso
  palpNefPartitions
///

doc ///
Key
  naiveNefPartitions
  (naiveNefPartitions, NormalToricVariety, ZZ)
Headline
  search directly for nef partitions of a toric variety
Usage
  L = naiveNefPartitions(V, c)
Inputs
  V:NormalToricVariety
  c:ZZ
    the number of parts
Outputs
  L:List
    of partitions of the rays of V into c parts, each part having nef sum
Description
  Text
    A nef partition of codimension $c$ is a partition of the rays into $c$ parts
    the sum of whose divisors is nef in each part.  This looks for them directly,
    with no help from PALP, which makes it a useful check on
    @TO palpNefPartitions@.
  Example
    needsPackage "NormalToricVarieties"
    V = smoothFanoToricVariety(3, 5);
    #rays V
    L = naiveNefPartitions(V, 2)
    assert all(L, p -> all(p, b -> isNef sum(b, i -> V_i)))
  Text
    In codimension one the only candidate is the anticanonical divisor.
  Example
    naiveNefPartitions(V, 1)
Caveat
  Every set partition of the rays into c parts is examined, so this is only
  practical when there are few rays.
SeeAlso
  palpNefPartitions
  nefPartitionDivisors
///

doc ///
Key
  nefPartitionDivisors
  (nefPartitionDivisors, NormalToricVariety, List)
Headline
  the divisors belonging to a nef partition
Usage
  Ds = nefPartitionDivisors(V, part)
Inputs
  V:NormalToricVariety
  part:List
    of ray indices, or a list of lists of them
Outputs
  Ds:List
    of @TO ToricDivisor@s, one for each part
Description
  Text
    Any rays not named make up one further part, so a codimension two partition
    can be given by naming one of its two parts, and the output of
    @TO palpNefPartitions@, which names all but the last part, can be used
    directly.
  Example
    needsPackage "NormalToricVarieties"
    V = smoothFanoToricVariety(3, 5);
    Ds = nefPartitionDivisors(V, {3, 5})
    assert all(Ds, isNef)
    assert(sum Ds == - toricDivisor V)
  Text
    In higher codimension the parts are given as a list of lists, and again the
    rays left over make up the last part.
  Example
    Es = nefPartitionDivisors(V, {{0, 2}, {1, 3, 4}})
    assert all(Es, isNef)
    assert(sum Es == - toricDivisor V)
  Text
    The divisors always add up to the anticanonical divisor, since the parts
    partition the rays.
SeeAlso
  naiveNefPartitions
  palpNefPartitions
///

TEST ///
  -- toPalp on matrices, and the round trip through fromPalpMatrix
  M = matrix{{1,1,1,1},{0,1,2,3}}
  assert(toPalp M == "2 4\n1 1 1 1\n0 1 2 3\n")
  assert(M == fromPalpMatrix toPalp M)
  assert(transpose M == fromPalpMatrix toPalp transpose M)

  -- from smoothFanoToricVariety(3, 5)
  N = transpose matrix {{1, 0, 0}, {-1, 0, 1}, {0, 1, 0}, {0, -1, 1}, {0, 0, 1}, {0, 0, -1}}
  assert(N == fromPalpMatrix toPalp N)

  -- weight systems are passed through as a single line
  assert(toPalp {10,1,2,3,4} == "10 1 2 3 4\n")

  assert(try (toPalp matrix(QQ, {{1,2}}); false) else true)      -- not over ZZ
  assert(try (toPalp {1, 1/2}; false) else true)                 -- not integers
  assert(try (fromPalpMatrix "   \n  \n"; false) else true)      -- blank output
///

TEST ///
  -- palpProgram picks a variant and caches it
  assert(instance(palpProgram "poly", Program))
  assert(palpProgram "poly" === palpProgram "poly")          -- cached
  assert(palpDmax == 11)
  assert((palpProgram "poly")#"name" == "poly-11d.x")

  palpDmax = 6
  assert((palpProgram "poly")#"name" == "poly.x")
  palpDmax = 11

  assert(try (palpDmax = 7; palpProgram "poly"; false) else true)
  palpDmax = 11
///

TEST ///
  -- runPALP leaves nothing behind in the current directory
  debug PALPInterface
  before = set readDirectory ".";
  out = palpRun("poly", "-v", {10,1,2,3,4})
  assert(match("Vertices of P", out))
  assert(set readDirectory "." === before)
///

TEST ///
  -- PALP announces a too-small POLY_Dmax on stdout with exit status 0;
  -- runPALP has to turn that into an actual error.
  ws = {10, 1, 1, 2, 2, 2, 2}    -- dim N = 5, so codim 3 needs POLY_Dmax >= 7
  palpDmax = 11
  assert(match("#part", palpRun("nef", "-Lp -c3", ws)))
  palpDmax = 6
  assert(try (palpRun("nef", "-Lp -c3", ws); false) else true)
  palpDmax = 11
///

TEST ///
  -- PALP orients its output block per example, not per option: poly.x -p prints
  -- 23 x 3 for {10,1,2,3,4} but 2 x 10 for {3,1,1,1}.  Every wrapper must hand
  -- back one point per column regardless.
  for ws in {{10,1,2,3,4}, {30,4,4,6,7,9}, {3,1,1,1}, {5,1,1,1,1,1}} do (
      mv := palpMVertices ws;
      d := numRows mv;
      assert(d === #ws - 2);                       -- dim of Delta(q)
      assert(numColumns mv > d);
      for M in {palpMPoints ws, palpNVertices ws, palpNPoints ws, palpNormalForm ws} do (
          assert(M =!= null);
          assert(numRows M === d);
          assert(numColumns M >= d);
          );
      -- the M-lattice points really are the lattice points of the polytope
      P := convexHull mv;
      assert(set entries transpose palpMPoints ws
          === set entries transpose matrix {latticePoints P});
      -- and the N-lattice vertices are the vertices of the dual
      assert(set entries transpose palpNVertices ws
          === set entries transpose lift(vertices polar P, ZZ));
      );

  -- the specific shapes that used to come back transposed
  assert((numRows palpMPoints {10,1,2,3,4}, numColumns palpMPoints {10,1,2,3,4}) === (3, 23))
  assert((numRows palpMPoints {3,1,1,1}, numColumns palpMPoints {3,1,1,1}) === (2, 10))
  assert((numRows palpNPoints {30,4,4,6,7,9}, numColumns palpNPoints {30,4,4,6,7,9}) === (4, 21))
///

TEST ///
  -- palpNormalForm is an invariant of the lattice polytope
  A = matrix {{0, 0, -3, 2}, {-5, 3, 0, 0}, {-3, -6, 1, 0}, {-3, -2, -1, 1}}
  assert(det A == 1 or det A == -1)
  M = fromPalpMatrix "4 12
   1   0   0   1   1  -1   0  -2   4  -2   0   2
   0   1   0   0  -1   0   0   3  -4   1  -1  -3
   0   0   1  -1   0   0   0   1  -4   3  -1  -3
   0   0   0   0   0   0   1  -1   1  -1   1   1
   "
  assert(numRows M == 4 and numColumns M == 12)
  assert(palpNormalForm M == M)
  assert(palpNormalForm (A*M) == M)

  -- a matrix and the weight system it comes from give the same normal form
  ws = {10,1,2,3,4}
  assert(palpNormalForm ws == palpNormalForm palpMVertices ws)

  -- two weight systems describing the same polytope
  assert(palpNormalForm {3402, 40, 41, 486, 1134, 1701}
      == palpNormalForm {3486, 41, 42, 498, 1162, 1743})
///

TEST ///
  -- a polytope that is IP but not reflexive: PALP has M-lattice answers for it,
  -- but its dual has no lattice description, and poly.x reports that by printing
  -- nothing (-d) or inequalities with an extra inhomogeneous entry (-e).
  M = transpose matrix{{1,-1,-1}, {1,-1,1}, {1,1,-1}, {1,1,2}, {-1,0,0}}
  assert not isReflexive convexHull M
  assert(palpMVertices M == M)
  assert(numColumns palpMPoints M == #latticePoints convexHull M)
  assert(palpNVertices M === null)
  assert(palpNPoints M === null)

  -- but fromPalpMatrix, the raw reader, still refuses to guess
  assert(try (fromPalpMatrix palpRun("poly", "-e", M); false) else true)
///

TEST ///
  -- palpInfo collects the summary line of poly.x -g
  H = palpInfo {10,1,1,2,2,2,2}
  assert(H#"M" === {378, 6})
  assert(H#"N" === {8, 6})
  assert(H#"H" === {2, 0, 350})
  assert(H#"Euler" === 2160)

  V = transpose matrix {{1,0,0},{-1,0,1},{0,1,0},{0,-1,1},{0,0,1},{0,0,-1}}
  H2 = palpInfo V
  assert(H2#"M" === {7, 6})
  assert(H2#"Cor" === {0})
  assert(not H2#?"Euler")
///

TEST ///
  -- cws.x enumerates the weight systems of a given dimension
  wss = weightSystems 3
  assert(#wss == 95)
  assert(wss#0 == {4, 1, 1, 1, 1})
  assert(wss#90 == {44, 4, 5, 13, 22})
  assert all(wss, ws -> ws#0 == sum drop(ws, 1))     -- degree is the sum of the weights
  assert all(wss, ws -> #ws == 5)

  -- a degree range; cws.x then appends a summary line, which must not be parsed
  assert(#weightSystems(5, Degrees => (20,20)) == 61)
  assert(#weightSystems(5, Degrees => (10,10)) == 5)
  assert(weightSystems(5, Degrees => 10) == weightSystems(5, Degrees => {10,10}))
  assert(member({10, 1, 1, 2, 2, 2, 2}, weightSystems(5, Degrees => 10)))

  assert(try (weightSystems(3, Degrees => {1,2,3}); false) else true)
///

TEST ///
  -- weightSystemPolytope is computed from the definition, with no call to PALP;
  -- it must agree with what poly.x -v reports.
  for ws in {{10,1,2,3,4}, {30,4,4,6,7,9}, {3,1,1,1}, {44,4,5,13,22},
             {3402,40,41,486,1134,1701}, {5,1,1,1,1,1}} do (
      P := weightSystemPolytope ws;
      assert(dim P == #ws - 2);
      assert isReflexive P;
      assert(palpNormalForm lift(vertices P, ZZ) == palpNormalForm weightSystemVertices ws);
      );

  -- weightSystemVertices is poly.x on the weight system
  assert(weightSystemVertices {10,1,2,3,4} == palpMVertices {10,1,2,3,4})

  -- combined weight systems go through PALP, but not through the definition
  cws = {3,1,1,1,0,0,0, 3,0,0,0,1,1,1}
  assert(numRows weightSystemVertices cws == 4)
  assert(try (weightSystemPolytope cws; false) else true)
  assert(try (weightSystemPolytope {10,1,2,3}; false) else true)   -- degree is not the sum
///

TEST ///
  -- ws = {10,1,1,2,2,2,2} has dim N = 5, so codim 3 needs POLY_Dmax >= 7 and
  -- only runs under the 11d build.  This is the main reason for wanting it.
  ws = {10, 1, 1, 2, 2, 2, 2}
  assert(palpDmax == 11)

  H2 = palpNefPartitions(ws, 2)
  assert(numRows H2#"Points" == 5 and numColumns H2#"Points" == 8)
  assert(#H2#"Partitions" == 3)
  assert(H2#"Partitions" == {{{1, 5, 6}}, {{2, 3, 4}}, {{3, 4}}})
  -- in codimension 2 PALP prints a single part; the rest is its complement
  assert all(H2#"Partitions", p -> #p == 1)

  -- codimension 3 prints two parts per partition.  The old parser only ever
  -- picked up the first of them.
  H3 = palpNefPartitions(ws, 3)
  assert(#H3#"Partitions" == 1)
  assert(H3#"Partitions" == {{{1, 5, 6}, {3, 4}}})
  assert all(H3#"Partitions", p -> #p == 2)

  -- every index is a column of the reported point matrix
  assert all(flatten flatten H3#"Partitions", i -> i >= 0 and i < numColumns H3#"Points")

  -- a matrix input gives the same answer as the weight system
  assert((palpNefPartitions(weightSystemVertices ws, 2))#"Partitions" == H2#"Partitions")
///

TEST ///
  -- no nef partitions is an answer, not an error
  H = palpNefPartitions({14, 1, 2, 2, 2, 2, 2, 3}, 3)
  assert(H#"Partitions" === {})
  assert(H#"Points" === null)

  -- the parser on its own
  assert(parsePALPNefPartitions "" === {})
  assert(parsePALPNefPartitions "H:2 86 [-168] P:0 V:1 5  6   (2 8) (1 4)     0sec  0cpu\n"
      === {{{1, 5, 6}}})
  assert(parsePALPNefPartitions "H:20 [24] P:0 V0:1 5  6  V1:3 4   (2 4 4) (1 2 2)  1sec\n"
      === {{{1, 5, 6}, {3, 4}}})
  -- lines that are not partitions are skipped
  assert(parsePALPNefPartitions "5 8  Points of Poly in N-Lattice:\nnp=3 d:0 p:1\n" === {})
///

TEST ///
  -- nef partitions found directly, with no help from PALP
  needsPackage "NormalToricVarieties"
  V = smoothFanoToricVariety(3, 5)
  n = #rays V
  assert(n == 6)

  -- codimension 1: the anticanonical divisor, which is nef since V is Fano
  assert(naiveNefPartitions(V, 1) == {{toList(0..n-1)}})

  P2 = naiveNefPartitions(V, 2)
  assert(#P2 > 0)
  for p in P2 do (
      assert(#p == 2);
      assert(sort flatten p == toList(0..n-1));          -- a partition of the rays
      assert all(p, b -> isNef sum(b, i -> V_i));        -- each part is nef
      );
  assert(#unique P2 == #P2)

  -- more parts than rays is empty, not an error
  assert(naiveNefPartitions(V, n + 1) === {})
  assert(try (naiveNefPartitions(V, 0); false) else true)
///

TEST ///
  -- nefPartitionDivisors completes a partition and adds up to the anticanonical
  needsPackage "NormalToricVarieties"
  V = smoothFanoToricVariety(3, 5)
  n = #rays V

  for p in naiveNefPartitions(V, 2) do (
      Ds := nefPartitionDivisors(V, p);
      assert(#Ds == 2);
      assert all(Ds, isNef);
      assert(sum Ds == - toricDivisor V);
      );

  -- codimension two may be given as the one part, with the rest implied
  D = nefPartitionDivisors(V, {2, 3})
  assert(#D == 2)
  assert(D#0 == V_2 + V_3)
  assert(sum D == - toricDivisor V)
  assert(nefPartitionDivisors(V, {{2, 3}}) == D)

  -- and PALP's output shape, which names all but the last part, works too
  E = nefPartitionDivisors(V, {{0, 1}, {2, 3}})
  assert(#E == 3)
  assert(sum E == - toricDivisor V)

  assert(try (nefPartitionDivisors(V, {0, 0}); false) else true)      -- not disjoint
  assert(try (nefPartitionDivisors(V, {n}); false) else true)         -- out of range
  assert(try (nefPartitionDivisors(V, {}); false) else true)          -- empty
///

TEST ///
-- MES test
  -- here are small examples of canonical polytopes which are not reflexive, in dimensions 3,4,5.
  M3 = matrix{{1,0,0,-1},{0,1,0,-1},{0,0,1,-2}}
  P3 = convexHull M3
  assert not isReflexive P3
  assert(1 === # interiorLatticePoints P3)
  palpInfo M3
  
  M4 = fromPalpMatrix"4 5
    1 0 0 0 -1
    0 1 0 0 -1
    0 0 1 0 -1
    0 0 0 1 -3"
  P4 = convexHull M4
  assert not isReflexive P4
  assert(1 === # interiorLatticePoints P4)
  
  M5 = fromPalpMatrix"5 6
    1 0 0 0 0 -1
    0 1 0 0 0 -1
    0 0 1 0 0 -1
    0 0 0 1 0 -1
    0 0 0 0 1 -2"
  P5 = convexHull M5
  assert not isReflexive P5
  assert(1 === # interiorLatticePoints P5)


 P =  weightSystemPolytope{7,1,1,1,2,2}
 isReflexive P
///

TEST ///
  -- MES test.  d=1 on all functions (not too many such polytopes in d=1!).
  -- Up to a change of basis [-1,1] is the only canonical polytope in dimension
  -- one: any longer interval containing the origin has a second interior point.
  V = matrix{{1,-1}}
  assert(1 === #interiorLatticePoints convexHull V)
  assert(2 === #interiorLatticePoints convexHull matrix{{-1,2}})

  assert(palpMVertices V == matrix{{1,-1}})
  assert(palpNVertices V == matrix{{1,-1}})
  assert(palpNormalForm V == matrix{{1,-1}})
  assert(set entries transpose palpMPoints V === set {{1},{-1},{0}})
  assert(set entries transpose palpNPoints V === set {{1},{-1},{0}})
  assert palpIsReflexive V
  assert(palpInfo V === new HashTable from {"M" => {3,2}, "N" => {3,2}})

  -- everything keeps one point per column here too
  for f in {palpMVertices, palpMPoints, palpNVertices, palpNPoints, palpNormalForm}
  do assert(numRows f V === 1)

  -- there is exactly one weight system in dimension one
  assert(weightSystems 1 == {{2,1,1}})
  assert(lift(vertices weightSystemPolytope {2,1,1}, ZZ) == matrix{{-1,1}})

  -- PALP itself aborts on a dimension one weight system, failing an assertion
  -- in ReadCwsPp, so that case is refused before PALP is called at all
  ans = trap weightSystemVertices {2,1,1}
  assert(ans#0 === null and instance(ans#1, Error))
  assert match("dimension one", toString ans#1)
///

TEST ///
  -- palpIsReflexive asks poly.x -e for the vertices of the dual: a reflexive
  -- polytope has them, one that is not gets inequalities instead, which come
  -- back as null.  poly.x -g would answer too, via its N: and F: fields, but it
  -- counts the dual's lattice points and is far slower in higher dimensions.
  assert palpIsReflexive {10, 1, 2, 3, 4}
  assert palpIsReflexive matrix{{1,0,0,-1},{0,1,0,-1},{0,0,1,-1}}

  -- canonical but not reflexive, in dimensions 3, 4 and 5
  M3 = matrix{{1,0,0,-1},{0,1,0,-1},{0,0,1,-2}}
  M4 = id_(ZZ^4) | transpose matrix {{-1,-1,-1,-3}}
  M5 = id_(ZZ^5) | transpose matrix {{-1,-1,-1,-1,-2}}
  for M in {M3, M4, M5} do (
      assert not palpIsReflexive M;
      assert(palpIsReflexive M === isReflexive convexHull M);
      assert(1 === #interiorLatticePoints convexHull M);
      assert((palpInfo M)#?"F" and not (palpInfo M)#?"N");
      );

  -- the fan polytope of P(1,a_1,..,a_d) is reflexive exactly when every weight
  -- divides the degree 1 + sum a_i
  for a1 from 1 to 4 do for a2 from a1 to 4 do for a3 from a2 to 4 do (
      a := {a1, a2, a3};
      M := id_(ZZ^3) | transpose matrix {-a};
      if dim convexHull M =!= 3 then continue;
      d := 1 + sum a;
      assert(palpIsReflexive M === all(prepend(1, a), q -> d % q === 0));
      );
///

TEST ///
  -- PALP has several compile-time bounds besides POLY_Dmax, and reports running
  -- out of any of them on stdout while exiting with status 0.  poly.x -N hits
  -- SYM_Nmax on the simplex in dimension eight.  This is well above the
  -- dimensions PALP is ordinarily used in; the point of the test is that the
  -- message becomes an error rather than a confusing parse failure.
  V8 = id_(ZZ^8) | transpose matrix {toList(8 : -1)}
  assert(numRows palpMVertices V8 == 8)            -- the cheap calls are fine
  assert(numRows palpNVertices V8 == 8)
  ans = trap palpNormalForm V8
  assert(ans#0 === null and instance(ans#1, Error))
  assert match("SYM_Nmax", toString ans#1)

  -- and POLY_Dmax is still reported, with advice to change palpDmax
  palpDmax = 6
  ans2 = trap palpNefPartitions({10,1,1,2,2,2,2}, 3)
  palpDmax = 11
  assert(ans2#0 === null and instance(ans2#1, Error))
  assert match("POLY_Dmax", toString ans2#1)
  assert match("palpDmax", toString ans2#1)
///

TEST ///
  -- Every wrapper keeps one point per column, in each dimension.  PALP is
  -- mainly used in low dimensions and is not much exercised above seven, so
  -- this stops at six; poly.x -p also slows sharply beyond that, taking three
  -- seconds in dimension eight against twenty milliseconds in dimension six.
  for d from 1 to 6 do (
      V := id_(ZZ^d) | transpose matrix {toList(d : -1)};
      assert(numRows palpMVertices V === d);
      assert(numRows palpMPoints V === d);
      assert(numRows palpNVertices V === d);
      assert palpIsReflexive V;
      assert(numRows palpNormalForm V === d);
      );
///

TEST ///
  -- Input PALP cannot handle is refused here rather than passed on.  See
  -- PALP-upstream-notes.md for what PALP does with each of these.
  debug PALPInterface

  -- not full dimensional: a square in the plane z = 0 of ZZ^3.  poly.x -v
  -- returns three of its four vertices, and poly.x -g exits on a bus error.
  sq = matrix{{1,-1,-1,1},{1,1,-1,-1},{0,0,0,0}}
  assert(4 == numColumns vertices convexHull sq)
  assert not palpFullDimensional sq
  for f in {palpMVertices, palpMPoints, palpNVertices, palpNPoints, palpInfo} do (
      r := trap f sq;
      assert(r#0 === null and instance(r#1, Error));
      assert match("span the lattice affinely", toString r#1);
      );
  assert(try (palpRun("poly", "-v", sq); false) else true)

  -- fewer points than the dimension cannot span it
  for M in {matrix{{1,-1},{0,0},{0,0}}, matrix{{0},{0},{0}}} do (
      r := trap palpMVertices M;
      assert(r#0 === null and instance(r#1, Error));
      assert match("more points than the dimension", toString r#1);
      );

  -- a weight system of dimension one makes PALP abort
  r = trap weightSystemVertices {2,1,1}
  assert(r#0 === null and instance(r#1, Error))
  assert match("dimension one", toString r#1)

  -- and the ordinary cases still go through
  assert palpFullDimensional matrix{{1,0,-1},{0,1,-1}}
  assert(numColumns palpMVertices {10,1,2,3,4} == 7)
///

end--

-* Development section *-
restart
needsPackage "PALPInterface"
check "PALPInterface"

restart
uninstallPackage "PALPInterface"
restart
installPackage "PALPInterface"
check "PALPInterface"

-- Questions to check:
--  weight systems: do the integers need to be in monotone increasing order?
--    {q, n_0, ..., n_r}, sum(n_i) = q, and n_0 <= n_1 <= ... <= n_r.
--    ANSWER: the n_i do NOT need to be numbered.
--    TODO: check that the integers are >= 1, sum to 0th element?  Otherwise we get a strange error...
--    We could also insist that they are in order, but not sure we care?
ws1 = weightSystemVertices {10,1,2,3,4}
ws2 = weightSystemVertices {10,1,2,4,3}
assert(palpNormalForm ws1 == palpNormalForm ws2)
ws3 = weightSystemVertices {10,4,3,2,1}
assert(palpNormalForm ws1 == palpNormalForm ws3)
ws4 = weightSystemVertices {10,2,3,6,0,1}
weightSystemVertices {5,2,2,1}
weightSystemVertices {3,1,1,1}

-- TODO:
-- doc for weightSystems
--   The description is weak.  Define a weight system.  What is the degree vs weight of a weight system?
--     which partitions are removed?
--   The user doesn't know what cws.x -w is.  It is nice to keep it here, but should not be the lead sentence, and should
--     mention that this is a palp command.
--   For dim >= 4, I think we need to limit the degree range automatically.  I don't think it finishes otherwise (I cancelled it
--     after 100 seconds elapsed.
--     and it is annoyingly large in any case.
help weightSystems
elapsedTime weightSystems 3;
# elapsedTime weightSystems(4, Degrees => 30) == 137
# elapsedTime weightSystems(4, Degrees => 100) == 781


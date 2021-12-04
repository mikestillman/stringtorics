-- Code to compute the GV invariants of a set of curve classes on a CY 3-fold.
needsPackage "StringTorics"

GVinvariants = method(Options => {
        Precision => 100, -- number of bits of accuracy
        Bound => infinity, -- max Heft bound
        Heft => List, -- list of positive weights
        })

-*
findFirstUnitVectors = method()
findFirstUnitVectors Matrix := List => (M) -> (
    -- returns the indices of the the columns which are unit vectors,
    -- if there are two or more columns corresponding to the same unit vector, take the first.
    -- the result is in increasing list of integers.
    e := entries transpose M;
    unitposition := (col) -> if all(col, a -> a >= 0) and sum col === 1 then position(col, a -> a === 1) else null;
    H := partition(c -> unitposition e_c, toList(0..numcols M-1));
    ks := select(keys H, k -> k =!= null);
    sort for k in ks list min(H#k)
    )

-- extend the nice unit vectors to a lisp p of n columns (n = numrows M)
-- so that det M_p = 1 or -1.
-- compute the inverse A of this n x n matrix.
-- then compute A^-1 * M, and return this list of columns.

findInvertibleSubmatrix = method(Options => {Limit => 10000})
findInvertibleSubmatrix(Matrix, List) := List => opts -> (M, p) -> (
    -- M is a matrix over the integers
    -- p is a list of column indices of M (indices run from 0 to numcols M - 1)
    -- Find a list q containing p (if possible), of column indices,
    --   s.t. det(M_q) = 1 or -1.
    -- Return q, or null, if one cannot be found.
    S := set toList(0..numcols M - 1);
    others := sort toList(S - set p);
    needed := numrows M - #p;
    if binomial(#others, needed) > opts.Limit then return null;
      -- really: do some random choices of q containing p.
    trythese := subsets(others, needed);
    tried := 0;
    for q1 in trythese do (
        q := join(q1,p);
        tried = tried+1;
        if abs det(M_q) == 1 then (
            << "found suitable set of columns after " << tried << " step" 
            << if tried == 1 then "" else "s" << endl;
            return sort q;
            );
        );
    << "tried all determinants: none had unit determinant" << endl;
    null
    )

fixGLSMDegrees = method()
fixGLSMDegrees NormalToricVariety := List => (V) -> (
    D := transpose degrees ring V;
    M := matrix D;
    p := findFirstUnitVectors M;
    q := findInvertibleSubmatrix(M, p);
    A := M_q;
    {A^-1 * M, q}
    )
*-

gvInputFromToric = method(Options => {
    Heft => null,
    DegreeLimit => infinity,
    Precision => 150
    }
)

gvInputFromToric(NormalToricVariety, List) := opts -> (V, moriGenerators) -> (
    str1 := toString toArray moriGenerators;
    str2 := "[]";
    str3 := toString  toArray if opts.Heft === null then heft ring V else opts.Heft;
    str4 := toString toArrayTable transpose degrees ring V;
    str5 := toString toArray intersectionNumbers(V, X);
    concatenate between("\n", {str1, str2, str3, str4, str5})
    )

-- We need to construct the following data for computeGV...
-- (1) Mori cone generators (Input)
-- (2) empty list [] for now
-- (3) good heft vector (Input, or use one from the ring)
-- (4) GLSM degrees (Take from V)
-- (5) intersection numbers (Take from V)
-- (6) [degree to compute to, precision for RR] (optional inputs).

end--

restart 
load "GV.m2"

topes = kreuzerSkarke(5, Limit => 50)
A = matrix topes_45
P = convexHull A
(V, q) = reflexiveToSimplicialToricVarietyCleanDegrees(P, CoefficientRing => ZZ/32003)
GLSM = transpose matrix degrees ring V

p = findFirstUnitVectors GLSM
q = findInvertibleSubmatrix(GLSM, p)

fixGLSMDegrees V
transpose LLL syz first oo

classGroup V
picardGroup V


topes = kreuzerSkarke(11, 24, Limit => 50)
A = matrix topes_45
A = matrix topes_35
P = convexHull A
(V, q) = reflexiveToSimplicialToricVarietyCleanDegrees(P, CoefficientRing => ZZ/32003)
GLSM = transpose matrix degrees ring V
heft ring V

transpose matrix rays V * matrix degrees ring V == 0
classGroup V
picardGroup V

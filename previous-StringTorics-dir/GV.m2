-- Code to compute the GV invariants of a set of curve classes on a CY 3-fold.
needsPackage "StringTorics"

GVinvariants = method(Options => {
        Precision => 100, -- number of bits of accuracy
        Bound => infinity, -- max Heft bound
        Heft => List -- list of positive weights
        })

gvInputFromToric = method(Options => {
    Heft => null,
    DegreeLimit => infinity,
    Precision => 150
    }
)

-- Can we have the input be:
--  a. rays for the fan
--  b. max simplices
--  c. GLSM charge matrix
--  d. basis indices for GLSM.  
--   assume: this submatrix is the identity submatrix: GLSM_basisIndices

-- Some functions needed:
--   given matrix from KS, get one triangulation (perhaps give input to choose others)
--   given (rays, triang), create good GLSM matrix and bases.
--   given (rays,triang,GLSM,basisIndices), compute:
--     A. topological data
--     B. (just) intersection numbers
--     C. input to computeGV (also needs a list of curve classes).

-- WriteToricData:
--  PolytopeData (contains rays, GLSM charge, basis indices)
--  TopologicalDataOfCY3
-- should have
--  CY3HypersurfaceData
--    contains: PolytopeData, and triangulation.
--    

gvInputFromToric(NormalToricVariety, List) := opts -> (V, moriGenerators) -> (
    str1 := toString moriGenerators;
    str2 := "{}";
    str3 := toString if opts.Heft === null then heft ring V else opts.Heft; -- this might not be correct
    str4 := toString transpose degrees ring V;
    str5 := for x in pairs CY3NonzeroMultiplicities V list append(x#0, x#1);
    str6 := {if opts.DegreeLimit === infinity then -1 else opts.DegreeLimit, opts.Precision}
    concatenate between("\n", {str1, str2, str3, str4, str5, str6})
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

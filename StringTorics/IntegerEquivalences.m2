newPackage(
    "IntegerEquivalences",
    Version => "0.1",
    Date => "13 Jan 2024",
    Headline => "finding invertible integral matrices preserving points, linear forms and ideals",
    Authors => {{ Name => "", Email => "", HomePage => ""}},
    AuxiliaryFiles => false,
    DebuggingMode => true
    )

export {
    "extendToMatrix", -- extendToMatrix(List of integers) ==> Matrix (over ZZ).
    "genericLinearMap", -- genericLinearMap(R).  Constructs two new rings, T, U, a matrix A over T nxn, n = numgens R, and phi = map(U, U, transpose A).
    
    }

importFrom_"LLLBases"{"gcdLLL"};

extendToMatrix = method()
extendToMatrix List := Matrix => L -> (
    if not all(L, a -> instance(a, ZZ))
    then error "expected a list of integers";
    (g, A) := gcdLLL L; -- coming from LLLBases.
    transpose A
    -- TODO: should this insure that the matrix has determinant 1 (not -1)?
    -- I don't really need that...
    )

genericLinearMap = method(Options => {Variable => null})
genericLinearMap Ring := Sequence => opts -> R -> (
    -- R should be a polynomial ring in n >= 1 variables.
    n := numgens R;
    if n == 0 then (
        kk := coefficientRing R; -- TODO: if none, this should give a better error message
        A := map(kk^0, kk^0, {});
        return (A, id_R);
        );
    K := coefficientRing R;
    t := if opts.Variable === null then getSymbol "t" else opts.Variable;
    T := K[t_(1,1)..t_(n,n)];
    U := T [gens R, Join => false];
    A = map(T^n,,transpose genericMatrix(T, T_0, n, n));
    phi := map(U, U, transpose A);
    (A, phi)
    )

TEST ///
-*
  restart
  needsPackage "IntegerEquivalences"
*-
  A = extendToMatrix{10,15,6}
  assert(A * transpose matrix{{10,15,6}} == transpose matrix{{0,0,1}})
  assert(abs det A == 1)

  A = extendToMatrix{10,15,1}
  assert(A * transpose matrix{{10,15,1}} == transpose matrix{{0,0,1}})
  assert(abs det A == 1)

  A = extendToMatrix{10,15,0}
  assert(A * transpose matrix{{10,15,0}} == transpose matrix{{0,0,5}}) -- 5 is gcd!
  assert(abs det A == 1)
///

TEST ///
-*
  restart
  needsPackage "IntegerEquivalences"
*-
  -- trivial case
  R = ZZ[]
  (A, phi) = genericLinearMap R
  assert(numRows A == 0 and numcols A == 0)
  assert(ring A === coefficientRing target phi)
  assert(source phi === target phi)

  R = ZZ[a..d]
  (A, phi) = genericLinearMap R
  assert(numRows A == 4 and numcols A == 4)
  assert(ring A === coefficientRing target phi)
  assert(source phi === target phi)

  (A, phi) = genericLinearMap(R, Variable => symbol s)
  assert(numRows A == 4 and numcols A == 4)
  assert(ring A === coefficientRing target phi)
  assert(source phi === target phi)
  use ring A
  assert(A_(0,0) == s_(1,1))  
  
  RQ = QQ (monoid R)
  (A, phi) = genericLinearMap RQ
  assert(numRows A == 4 and numcols A == 4)
  assert(ring A === coefficientRing target phi)
  assert(source phi === target phi)
  
  R = ZZ[a]  
  (A, phi) = genericLinearMap R
  assert(numRows A == 1 and numcols A == 1)
  assert(ring A === coefficientRing target phi)
  assert(source phi === target phi)

  R = ZZ/101[a..d]
  (A, phi) = genericLinearMap R
  U = target phi
  assert(source phi === U)
  assert(ring A === coefficientRing U)
  for i from 0 to 3 do 
    assert(phi U_i == (A^{i} * (transpose vars U))_(0,0))

  R = ZZ[a..d]
  (A, phi) = genericLinearMap R
  U = target phi
  assert(source phi === U)
  assert(ring A === coefficientRing U)
  for i from 0 to 3 do 
    assert(phi U_i == (A^{i} * (transpose vars U))_(0,0))

  R = QQ[a..e]
  (A, phi) = genericLinearMap R
  U = target phi
  assert(source phi === U)
  assert(ring A === coefficientRing U)
  for i from 0 to numgens R - 1 do 
    assert(phi U_i == (A^{i} * (transpose vars U))_(0,0))
///

end--

-* Documentation section *-
beginDocumentation()

doc ///
Key
  IntegerEquivalences
Headline
Description
  Text
  Tree
  Example
  CannedExample
Acknowledgement
Contributors
References
Caveat
SeeAlso
Subnodes
///

doc ///
Key
Headline
Usage
Inputs
Outputs
Consequences
  Item
Description
  Text
  Example
  CannedExample
  Code
  Pre
ExampleFiles
Contributors
References
Caveat
SeeAlso
///

-* Test section *-
TEST /// -* [insert short title for this test] *-
-- test code and assertions here
-- may have as many TEST sections as needed
///

end--

-* Development section *-
restart
debug needsPackage "IntegerEquivalences"
check "IntegerEquivalences"

uninstallPackage "IntegerEquivalences"
restart
installPackage "IntegerEquivalences"
viewHelp "IntegerEquivalences"

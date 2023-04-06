newPackage(
    "DanilovKhovanskii",
    Version => "0.1",
    Date => "6 April 2023",
    Headline => "Computing Hodge-Deligne polynomials of toric hypersurfaces",
    Authors => {{ Name => "", Email => "", HomePage => ""}},
    PackageExports => {"Polyhedra", "StringTorics"},
    AuxiliaryFiles => false,
    DebuggingMode => true
    )

export {"computeSumqeZ", "eZ2hZ", "computeHodge"}

-* Code section *-
ehrhartNumerator = method();
ehrhartNumerator Polyhedron := P -> (
    d := dim P;
    x := getSymbol "x";
    R := QQ[x];
    f := 1 + sum for i from 1 to d + 1 list (
	a := #latticePoints(i * P);
	a * R_0^i
	);
    g := f * (1 - R_0)^(d + 1);
    for i from 0 to d list (
	lift(coefficient(R_0^i, g), ZZ)
	)
    )

computeSumqeZ = method();--from 4.6 of Danilov and Khovanskii [though that needs a factor of (-1)^d in front of \psi_{d+1}(\Delta)]
computeSumqeZ (Polyhedron, ZZ) := (P, p) -> (
    d := dim P;
    psi := ehrhartNumerator(P);
    --ehrhart(P) computes lattice points in first d dilations to figure out polynomial, so just find lattice points
    (-1)^(d - 1) * ((-1)^p * binomial(d, p + 1) + psi_(p+1))
    )

eZ2hZ = method();
eZ2hZ (Polyhedron, MutableHashTable) := (P, eZ) -> (
    hZ := new MutableHashTable;
    if isSimplicial P then (--not exactly right
	for key in keys eZ do hZ#key = (-1)^(key#0 + key#1) * eZ#key;
	)
    else (
	error "Not enough information. Can only recover Hodge numbers from the Hodge Deligne polynomial for smooth and/or simplicial varieties."
	);
    hZ
    )

computeHodge = method();
computeHodge Polyhedron := P -> (--P := X#"polytope data";
    --Will return both the Hodge numbers and the Hodge-Deligne numbers, in that order.
    d := dim P; print d;
    eZ := new MutableHashTable;
    eZbar := new MutableHashTable;
    --Assume P is full dimensional (X has no torus factors)
    --If dim P != dim X then (eZ = computeHodge(restrict P) * (x*y - 1)^(dim X - dim P)) else 
    
    --Begin by computing eZ of the varieties corresponding to each face of P.
    --This is known by induction.
    fs := flatten for n from 1 to d list facesAsPolyhedra(n,P); print(class(fs));
    eZfaces := for face in fs list (
	h := computeHodge(face);
	h#1); print"faces done";--Hodge-Deligne numbers of the face.
    
    --A couple of Lefschetz-type theorems and Gysin homomorphisms give eZ#(p, q) for p + q > d - 1
    --in terms of eT^d#(p + 1, q + 1)
    --For p + q > d - 1, eZ#(p, q) is 0 for p != q and is (-1)^(d + p + 1) * binomial(d, p + 1) for p == q.
    for p from floor(d / 2) to d - 1 do eZ#(p, p) = (-1)^(d + p + 1) * binomial(d, p + 1); print("p + q > d - 1");
    
    --This gives eZbar for p + q > d - 1.
    --Poincare dualtiy then gives eZbar#(d - 1 - p, d - 1 - q) = eZbar#(p, q).
    --Since eZbar#(p, q) is then known for p + q < d - 1, one can compute obtains eZ#(p, q) for p + q < d - 1.
    --Where to stop?
    for p from 0 to d - 1 do (
	for q from d - p to d - 1 do (print(p,q);
	    eZbar#(p, q) = (if eZ#?(p, q) then eZ#(p, q) else 0) + sum for eZface in eZfaces list (
		if eZface#?(p, q) then eZface#(p, q) else 0
		); print"a";
	    eZbar#(d - 1 - p, d - 1 - q) = eZbar#(p, q); print"b";
	    eZ#(d - 1 - p, d - 1 - q) = eZbar#(d - 1 - p, d - 1 - q) - sum for eZface in eZfaces list (
		if eZface#?(d - 1 - p, d - 1 - q) then eZface#(d - 1 - p, d - 1 - q) else 0); print"c";
	    );
	); print("p + q < d - 1");
    
    --The last remaining number, eZ#(p, d - 1 - p), is then the difference Sum_q eZ#(p, q) - Sum_{q != d - 1 - p} eZ#(p, q).
    --Sum_q eZ#(p, q) can be calculated from the number of lattice points in the interior of each face.
    for p from 0 to d - 1 do (
	eZ#(p, d - 1 - p) = computeSumqeZ(P, p) - sum for q from 0 to d - 1 list (
	    if q == d - 1 - p then continue else (if eZ#?(p, q) then eZ#(p, q) else 0)
	    );
	); print("p + q = d - 1");
    hZ := eZ2hZ(P, eZ);
    (hZ, eZ)
    )

computeHodge CYPolytopeData := P -> computeHodge(convexHull transpose matrix rays P)

computeHodge CYData := X -> computeHodge(convexHull transpose matrix rays X#"polytope data")


-* Documentation section *-
beginDocumentation()

doc ///
Key
  DanilovKhovanskii
Headline
  Computing Hodge-Deligne polynomials of toric hypersurfaces
Description
  Text
References
Caveat
SeeAlso
///

///
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
  P = convexHull transpose matrix {{1,1},{1,-1},{-1,1},{-1,-1}}
  ehrhartNumerator(P)
  R = QQ[x]
  #latticePoints(P)
  #latticePoints(2 * P)
  #latticePoints(3 * P)
  #latticePoints(4 * P)
  (1 + 9*x + 25*x^2 + 49*x^3 + 81*x^4) * (1 - x)^3
  assert (computeSumqeZ(P, 0) == 42)--fix 42
  assert (computeSumqeZ(P, 1) == 42)
  assert (computeSumqeZ(P, 2) == 42)
  P = convexHull transpose matrix {{1,1},{1,-1}}
///

end--

-* Development section *-
restart
path = append(path, "/mnt/c/Users/ccjle/stringtorics/")
debug needsPackage "DanilovKhovanskii"
check "DanilovKhovanskii"

uninstallPackage "DanilovKhovanskii"
restart
installPackage "DanilovKhovanskii"
viewHelp "DanilovKhovanskii"

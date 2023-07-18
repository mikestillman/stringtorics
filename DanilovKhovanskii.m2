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

export {"Cheap",
    "vdot",
    "torusFactor",
    "stdVector",
    "manyMatricesToLargeMatrix",
    "manyPolyhedraToLargeMatrix",
    "manyPolyhedraToLargeOne",
    "ehrhartNumerator",
    "ehrhartNumeratorQuicker",
    "computeSumqeZ",
    "eZ2hZ",
    "computeHodgeDeligne",
    "computeHdogeDeligneAffineAndTorus",
    "computeHodgeDeligneTorusCI",
    "FaceInfo"}

-* Code section *-
stdVector = method();--index from 0
stdVector (ZZ, ZZ) := (n, i) -> (
    for j from 0 to n - 1 list (if j == i then 1 else 0)
    )

manyMatricesToLargeMatrix = method();
manyMatricesToLargeMatrix List := Ms -> (
    r := #Ms;
    transpose matrix prepend(for i from 1 to numrows(Ms#0) + r list 0, flatten for i from 0 to r - 1 list (
	for j from 0 to numcols(Ms#i) - 1 list (
	    entries((Ms#i)_j) | stdVector(r, i)
	    )
	))
    )

manyPolyhedraToLargeMatrix = method();
manyPolyhedraToLargeMatrix List := Ps -> (
    Ms := for P in Ps list vertices P;
    manyMatricesToLargeMatrix(Ms)
    )

manyPolyhedraToLargeOne = method();
manyPolyhedraToLargeOne List := Ps -> (
    c := #Ps;
    convexHull manyPolyhedraToLargeMatrix(Ps)
    )

ehrhartNumerator = method();
ehrhartNumerator Polyhedron := P -> (
    d := dim P;
    t := getSymbol "t";
    R := QQ[t];
    f := 1 + sum for i from 1 to d list (
	a := #latticePoints(i * P);
	a * R_0^i
	);
    g := f * (1 - R_0)^(d + 1);
    for i from 0 to d list (
	lift(coefficient(R_0^i, g), ZZ)
	)
    )

ehrhartNumeratorQuicker = method();
ehrhartNumeratorQuicker Polyhedron := P -> (
    d := dim P;
    l := prepend(1, for i from 1 to ceiling(d / 2) list (
	#latticePoints(i * P)
	));
    lint := for i from 1 to floor(d / 2) list (
	#interiorLatticePoints(i * P)
	);
    prepend(1, for i from 1 to d list (
	if i <= ceiling(d / 2) then (
	    sum for j from max(0, i - d - 1) to i list (-- print(i, j);
		(-1)^(i - j) * binomial(d + 1, i - j) * l#(j)
		)
	    )
	else (
	    sum for j from max(1, i - d - 1) to d + 1 - i list (-- print(i, j);
		(-1)^(d + 1 - i - j) * binomial(d + 1, d + 1 - i - j) * lint#(j - 1)
		)		
	    )
	))
    )

computeSumqeZ = method();--from 4.6 of Danilov and Khovanskii [though that may need a factor of (-1)^d in front of \psi_{d+1}(\Delta)]
computeSumqeZ (Polyhedron, List, ZZ) := (P, psi, p) -> (
    d := dim P;
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

getSparseeZ = method();
getSparseeZ (MutableHashTable, Sequence) := (eZ, pq) -> (
    if eZ#?pq then eZ#pq else 0
    )

getSparseeZ (HashTable, Sequence) := (eZ, pq) -> (
    if eZ#?pq then eZ#pq else 0
    )

vdot = method();
vdot (List, List) := (a, b) -> if #a == #b then (
    sum for i from 0 to #a-1 list a#i*b#i) else (error "Lengths not compatible.")

liftToQQ = method();
liftToQQ Matrix := M -> (
    matrix for i from 0 to numrows(M) - 1 list (
	for j from 0 to numcols(M) - 1 list (
	    lift(M_(i, j), QQ)
	    )
	)
    )

torusFactor = method();
torusFactor (MutableHashTable, ZZ, ZZ) := (eZ, d, D) -> (--d = dimension of polytope; D = dimension of lattice
    dt := D - d;
    eZfull := new MutableHashTable from {};
    if dt == 0 then (print("No torus factors.");
	eZfull = eZ
	)
    else (print("Torus factors." | dt);
	for p from 0 to D - 1 do (
	    for q from 0 to D - 1 do (
		eZfull#(p, q) = sum for i from 0 to dt list (
		    (-1)^(dt - i) * binomial(dt, i) * getSparseeZ(eZ, (p - i, q - i))
		    );
		);
	    );
	);
    eZfull
    )

torusFactor (HashTable, ZZ, ZZ) := (eZ, d, D) -> (--d = dimension of polytope; D = dimension of lattice
    dt := D - d;
    eZfull := new MutableHashTable from {};
    if dt == 0 then (--no torus factors
	eZfull = eZ
	)
    else (
	for p from 0 to D - 1 do (
	    for q from 0 to D - 1 do (
		eZfull#(p, q) = sum for i from 0 to dt list (
		    (-1)^(dt - i) * binomial(dt, i) * getSparseeZ(eZ, (p - i, q - i))
		    );
		);
	    );
	);
    eZfull
    )

computeHodgeDeligne = method(Options => {FaceInfo => {true, new HashTable from {}, -1}});
--FaceInfo: first entry = if this is the full polytope; second = data from lower-dimensional faces;
--third = dimension of ambient variety (will usually, but not always, be the number of rows of the vertex matrix)
computeHodgeDeligne Polyhedron := opts -> P -> (--P := X#"polytope data";
    d := dim P;
    D := opts.FaceInfo#2;
    if D == -1 then (
	D = numrows vertices P
	);
    topdim := opts.FaceInfo#0; print("poly dim = "| d | ", ambient dim = " | D, topdim);
    eZ := new MutableHashTable;
    eZbar := new MutableHashTable;
    --Assume P is full dimensional (X has no torus factors)
    --If dim P != dim X then (eZ = computeHodgeDeligne(restrict P) * (x*y - 1)^(dim X - dim P)) else     
    
    if d == 0 then (--hypersurface is empty
	return (new HashTable from eZ, new HashTable from eZbar, new HashTable from {})
	)
    else if d == 1 then (--print("dim(Z) = 0"); --hypersurface in 1-dimension is a collection of points
	eZ#(0,0) = #latticePoints(P) - 1;
	eZbar#(0,0) = #latticePoints(P) - 1; -- print(eZ#(0,0), eZbar#(0,0));
	eZ = torusFactor(eZ, d, D);
	eZbar = torusFactor(eZbar, d, D);
	return (new HashTable from eZ, new HashTable from eZbar, new HashTable from {})
	);-- print("not 0 or 1");
    
    --Begin by computing eZ of the varieties corresponding to each face of P.
    --This is known by induction.
    eZfaces := new MutableHashTable from opts.FaceInfo#1;
    print(opts.FaceInfo#1);-- print(eZfaces);
    print("eZfaces: " | #eZfaces | " , keys(eZfaces): " | #(keys eZfaces));
    if #(keys eZfaces) == 0 then (--print("no face info");
	verts := entries transpose vertices P;-- print(verts);
	vIndices := new HashTable from for i from 0 to #verts - 1 list i => verts#i;
	oneFaces := facesAsPolyhedra(d - 1, P);
	Pfaces := faces(P);-- print(Pfaces);
	for i from 0 to #oneFaces - 1 do (--print(i, Pfaces#(d - 1)#i#0, oneFaces#i);
	    --eZfaces#(Pfaces#(d - 1)#i#0) = new MutableHashTable from {};
	    --eZfaces#(Pfaces#(d - 1)#i#0)#(0,0) = #latticePoints(oneFaces#i) - 1;
	    eZfaces#(Pfaces#(d - 1)#i#0) = new HashTable from {(0, 0) => #latticePoints(oneFaces#i) - 1};
	    --print(eZfaces#(Pfaces#(d - 1)#i#0));
	    );-- print("done 1");
	for n from 2 to d - 1 do (
	    Fs := facesAsPolyhedra(d - n, P); -- print("n = " | n | ", Fs: " | #Fs);
    	    for i from 0 to #Fs - 1 do (
		eZfaces2 := new HashTable from flatten for j from 1 to dim Fs#i - 1 list (--print(Pfaces#n#i);
		    flatten for k from 1 to n - 1 list (
			for l from 0 to #(Pfaces#(d - k)) - 1 list (-- print(Pfaces#(d - k)#l#0, Pfaces#(d - n)#i#0, isSubset(Pfaces#(d - k)#l#0, Pfaces#(d - n)#i#0));
			    if isSubset(Pfaces#(d - k)#l#0, Pfaces#(d - n)#i#0) and eZfaces#?(Pfaces#(d - k)#l#0) then Pfaces#(d - k)#l#0 => eZfaces#(Pfaces#(d - k)#l#0) else continue
			    )
			)
	            --for k in keys eZfaces list (print(k, Pfaces#(d - n)#i#0);
		    --if isSubset(k, Pfaces#(d - n)#i#0) then (print("yes"); k => eZfaces#k) else continue
		    );-- print("face ready"); print(eZfaces2);
		e := computeHodgeDeligne(Fs#i, FaceInfo => {false, eZfaces2, D});
	    	eZfaces#(Pfaces#(d - n)#i#0) = e#0;-- print(e);
		); print("done " | n);
	    );
	);-- print"faces done";--Hodge-Deligne numbers of the face.
    
    --A couple of Lefschetz-type theorems and Gysin homomorphisms give eZ#(p, q) for p + q > d - 1
    --in terms of eT^d#(p + 1, q + 1)
    --For p + q > d - 1, eZ#(p, q) is 0 for p != q and is (-1)^(d + p + 1) * binomial(d, p + 1) for p == q.
    for p from floor(d / 2) to d - 1 do eZ#(p, p) = (-1)^(d + p + 1) * binomial(d, p + 1); --print("p + q > d - 1");
    
    --This gives eZbar for p + q > d - 1.
    --Poincare dualtiy then gives eZbar#(d - 1 - p, d - 1 - q) = eZbar#(p, q).
    --Since eZbar#(p, q) is then known for p + q < d - 1, one can compute obtains eZ#(p, q) for p + q < d - 1.
    --Where to stop?
    for p from 0 to d - 1 do (
	for q from d - p to d - 1 do (--print(p,q);
	    eZbar#(p, q) = getSparseeZ(eZ, (p, q)) + sum (
		for k in keys eZfaces list getSparseeZ(eZfaces#k, (p, q))
		); --print"a";
	    eZbar#(d - 1 - p, d - 1 - q) = eZbar#(p, q); --print"b";
	    eZ#(d - 1 - p, d - 1 - q) = eZbar#(d - 1 - p, d - 1 - q) - sum (
		for k in keys eZfaces list getSparseeZ(eZfaces#k, (d - 1 - p, d - 1 - q))
		); --print"c";
	    );
	); --print("p + q < d - 1");
    
    --The last remaining number, eZ#(p, d - 1 - p), is then the difference Sum_q eZ#(p, q) - Sum_{q != d - 1 - p} eZ#(p, q).
    --Sum_q eZ#(p, q) can be calculated from the number of lattice points in the interior of each face.
    psi := ehrhartNumeratorQuicker(P);
    for p from 0 to d - 1 do (
	eZ#(p, d - 1 - p) = computeSumqeZ(P, psi, p) - sum (
	    for q from 0 to d - 1 list (
	    	if q == d - 1 - p then continue else getSparseeZ(eZ, (p, q))
	    	)
	    );
	eZbar#(p, d - 1 - p) = eZ#(p, d - 1 - p) + sum (
	    for k in keys eZfaces list getSparseeZ(eZfaces#k, (p, d - 1 - p))
	    );
	); --print("p + q = d - 1");
    --hZ := eZ2hZ(P, eZ);
    if topdim then (print("topdim = true");
	eZ = torusFactor(eZ, d, D);
    	eZbar = torusFactor(eZbar, d, D);
    	for k in keys(eZfaces) do (
	    eZfaces#k = torusFactor(eZfaces#k, d, D);
	    );
	);
    (new HashTable from eZ, new HashTable from eZbar, new HashTable from eZfaces)
    )

--check M versus N lattice
computeHodgeDeligne CYPolytope := P -> computeHodgeDeligne(convexHull transpose matrix rays P)

computeHodgeDeligne CalabiYauInToric := X -> computeHodgeDeligne(convexHull transpose matrix rays X#"polytope data")

--For (a hypersurface in) a toric variety that is a subset of the projective toric variety that has polytope P.
--Not tested.
computeHodgeDeligneInPToric = method();
computeHodgeDeligneInPToric (Polyhedron, List) := (P, whichFaces) -> (
    (eZ, eZbar, eZfaces) := computeHodgeDeligne(P);
    eZtoric := new MutableHashTable from {};
    for k in keys eZfaces do (
	eZtoric#k = eZ#k + sum for f in whichFaces list getSparseeZ(eZfaces#f, k);
	);
    eZtoric
    )

--For a hypersurface in T^n x C^r.
--Cheap option can be (but is not currently) used for complete intersections.
--not Cheap option not finished yet.
computeHodgeDeligneAffineAndTorus = method(Options => {Cheap => false});
computeHodgeDeligneAffineAndTorus (Polyhedron, ZZ, ZZ) := opts -> (P, n, r) -> (
    eZtoric := new MutableHashTable from {};
    subs := subsets(r);
    if opts.Cheap then ( print("Cheap");
	for s in subs do (--lambda_j == 0 <-> j in s
	    vs := vertices P; print(vs);
	    newP := convexHull vs_(for i from 0 to numcols vs - 1 list (
		    if (a := true;
			for j in s do (--column i must have not have std vector e_j
			    a = (a and vs_(n + j, i) == 1)
			    ); print(i, a);
			a
			)
		    then i else continue
		    )
		);
	    D := n + r - #s;--work in T^n x C^(r - #s)
	    (eZ, eZbar, eZfaces) := computeHodgeDeligne(newP, FaceInfo => {false, new HashTable, D}); 
	    for k in keys(eZ) do (print(k);
		eZtoric#k = getSparseeZ(eZtoric, k) + eZ#k
		); print("done subset:" | toString(s));
	    );
	)
    else (
	for s in subs do (
	    vs := vertices P; print(vs);
	    newP := P;--finish...
	    D := n + r - #s;--work in T^n x C^(r - #s)
	    (eZ, eZbar, eZfaces) := computeHodgeDeligne(newP, FaceInfo => {false, new HashTable, D}); 
	    eZ#s = eZ;
	    for k in keys(eZ) do (print(k);
		eZtoric#k = getSparseeZ(eZtoric, k) + eZ#k
		); print("done subset:" | toString(s));	    
	    );
	);
    new HashTable from eZtoric
    )

--For a hypersurface in T^n x C^r. Takes a matrix instead of a polyhedron.
--Cheap option is currently the routine used for complete intersections.
--not Cheap option not finished yet.
computeHodgeDeligneAffineAndTorus (Matrix, ZZ, ZZ) := opts -> (M, n, r) -> (
    eZtoric := new MutableHashTable from {};
    subs := subsets(r);
    if opts.Cheap then (print("Cheap");
	for s in subs do (--lambda_j == 0 <-> j in s
	    newP := convexHull M_(for i from 0 to numcols M - 1 list (
		    if (a := true;
			for j in s do (--column i must have not have std vector e_j
			    a = (a and M_(n + j, i) == 0)
			    ); print(i, a);
			a
			)
		    then i else continue
		    )
		);
	    D := n + r - #s;--work in T^n x C^(r - #s)
	    (eZ, eZbar, eZfaces) := computeHodgeDeligne(newP, FaceInfo => {false, new HashTable, D}); print(eZ);
	    for k in keys(eZ) do (print(k);
		eZtoric#k = getSparseeZ(eZtoric, k) + eZ#k
		); print("done subset:" | toString(s));
	    );
	)
    else (
	for s in subs do (
	    newP := convexHull M;--finish...
	    D := n + r - #s;--work in T^n x C^(r - #s)
	    (eZ, eZbar, eZfaces) := computeHodgeDeligne(newP, FaceInfo => {false, new HashTable, D}); 
	    eZ#s = eZ;
	    for k in keys(eZ) do (print(k);
		eZtoric#k = getSparseeZ(eZtoric, k) + eZ#k
		); print("done subset:" | toString(s));		    
	    );
	);
    new HashTable from eZtoric
    )

--For a complete intersection of hypersurfaces, Y, in a torus, T^n.
computeHodgeDeligneTorusCI = method();
computeHodgeDeligneTorusCI List := Ps -> (
    a := true;
    nr := numrows vertices Ps#0;
    for P in Ps do (
	if not class(P) === Polyhedron then a = false
	else if not numrows vertices P == nr then a = false else continue
	);
    if not a then error("Polyhedra do not sit in the same lattice. Make sure their vertex matrices have the same number of rows.");
    n := nr;
    r := #Ps;
    M := manyPolyhedraToLargeMatrix(Ps); print(M);
    eZtoric := computeHodgeDeligneAffineAndTorus(M, n, r, Cheap => true); print("eZtoric#(0, 0) = " | eZtoric#(0,0));
    eZCI := new MutableHashTable from {};
    for p from 0 to n - r do (
	for q from 0 to n - r do (print("(p, q) = " | toString(p, q));
	        if p == q then (
		    eZCI#(p, q) = (-1)^(n + p) * binomial(n, p) - eZtoric#(p + r - 1, q + r - 1)
		    )
		else (
		     eZCI#(p, q) = -eZtoric#(p + r - 1, q + r - 1)
		     )
	    );
	);
    (new HashTable from eZCI, eZtoric)
    )


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
  R = QQ[x]
  P = convexHull transpose matrix {{1,1},{1,-1},{-1,1},{-1,-1}}
  eNum = ehrhartNumerator(P)
  f = sum for i from 0 to #eNum - 1 list (
      eNum#i * x^i
      )
  h = 1 + sum for i from 1 to 10 list (
      #latticePoints(i * P) * x^i
      )
  assert((h * (1 - x)^3)%x^11 == f)
///

TEST ///
  P = convexHull transpose matrix {{1,1},{1,-1},{-1,1},{-1,-1}}
  psi = ehrhartNumerator(P)
  assert (computeSumqeZ(P, psi, 0) == -8)
  assert (computeSumqeZ(P, psi, 1) == 0)
  (eZ, eZbar, eZfaces) = computeHodgeDeligne(P)
  assert (eZ === new HashTable from {(0,0) => -7, (1,0) => -1, (0,1) => -1, (1,1) => 1})
  assert (eZbar === new HashTable from {(0,0) => 1, (1,0) => -1, (0,1) => -1, (1,1) => 1})
///

TEST ///
  PN = convexHull transpose matrix {{1, 0, 0}, {0, 1, 0}, {-1, -1, -2}, {0, 0, 1}, {0, 0, -1}}
  PM = polar PN
  isReflexive PN
  isReflexive PM
  latticePoints PM
  faces(1,PN)
  psi = ehrhartNumerator(PM)
  assert (computeSumqeZ(PM, psi, 0) == 33)
  assert (computeSumqeZ(PM, psi, 1) == 27)
  assert (computeSumqeZ(PM, psi, 2) == 2)
  (eZ, eZbar, eZfaces) = computeHodgeDeligne(PM)
  -- assert eZ === new HashTable from {(0,0) => 19, (1,0) => 7, (0,1) => 7, (2,0) => 1, (1,1) => 14, (0,2) => 1, (2,2) => 1})
  assert (eZ === new HashTable from {(0,0) => 20, (1,0) => 12, (0,1) => 12, (2,0) => 1, (1,1) => 15, (0,2) => 1, (2,2) => 1})
  assert (eZbar ===  new HashTable from {(0,0) => 1, (1,0) => 0, (0,1) => 0, (2,0) => 1, (1,1) => 20, (0,2) => 1, (2,1) => 0, (1,2) => 0, (2,2) => 1})
///

TEST ///
  topes = kreuzerSkarke 3;
  Q = cyPolytope topes_50
  assert(hh^(1,1) Q == 3)
  assert(hh^(1,2) Q == 73)
  PM = polytope(Q, "M")
  vertices PM
  isReflexive PM
  (eZ, eZbar, eZfaces) = computeHodgeDeligne(PM);
  assert(eZbar === new HashTable from {(0,0) => 1, (1,0) => 0, (0,1) => 0, (1,1) => 3, (2,0) => 0, (0,2) => 0, (3,0) => -1, (2,1) => -73, (0,3) => -1, (1,2) => -73, (1,3) => 0, (2,2) => 3,
      (3,1) => 0, (2,3) => 0, (3,2) => 0, (3,3) => 1})

  Q = cyPolytope topes_70
  assert(hh^(1,1) Q == 3)
  assert(hh^(1,2) Q == 75)
  PM = polytope(Q, "M")
  vertices PM
  isReflexive PM
  (eZ, eZbar, eZfaces) = computeHodgeDeligne(PM);
  assert(eZbar === new HashTable from {(0,0) => 1, (1,0) => 0, (0,1) => 0, (1,1) => 3, (2,0) => 0, (0,2) => 0, (3,0) => -1, (2,1) => -75, (0,3) => -1, (1,2) => -75, (1,3) => 0, (2,2) => 3,
      (3,1) => 0, (2,3) => 0, (3,2) => 0, (3,3) => 1})

  matrix for i from 0 to 3 list (
      for j from 0 to 3 list(
	  (-1)^(i + j) * eZbar#(i,j)
	  )
      )
///

TEST ///
  P = convexHull transpose matrix {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}
  (eZ, eZbar, eZfaces) = computeHodgeDeligne(P)
  assert (eZ === new HashTable from {(0,0) => 5, (1,0) => 0, (0,1) => 0, (2,0) => 1, (1,1) => 0, (0,2) => 1, (2,2) => 1})
  assert (eZbar === new HashTable from {(0,0) => 1, (1,0) => 0, (0,1) => 0, (2,0) => 1, (1,1) => 8, (0,2) => 1, (2,1) => 0, (1,2) => 0, (2,2) => 1})
///

TEST ///
  d = 2
  Q2 = 2 * stdSimplex(d)--not full dimensional
  Q3 = 3 * stdSimplex(d)
  (eZ, eZbar, eZfaces) = computeHodgeDeligne(Q2)
  assert(eZ === new HashTable from {(0,0) => 5, (0,1) => 0, (1,0) => 0, (2,0) => 0, (0,2) => 0, (1,1) => -6, (2,1) => 0, (1,2) => 0, (2,2) => 1})
  (eZ, eZbar, eZfaces) = computeHodgeDeligne(Q3)
  assert(eZ === new HashTable from {(0,0) => 8, (0,1) => 1, (1,0) => 1, (2,0) => 0, (0,2) => 0, (1,1) => -9, (2,1) => -1, (1,2) => -1, (2,2) => 1})
  --computeHodgeDeligneTorusCI({Q2, Q3}, d)
  Q2 = 2 * stdSimplex(5)
  Q3 = 3 * stdSimplex(5)
  --eZCI = computeHodgeDeligneTorusCI({Q2, Q3})
  --assert(eZCI === new HashTable from {(0,0) => 58, (1,0) => 0, (0,1) => 0, (1,1) => 105, (0,2) => 0, (2,0) => 0, (3,0) => 0, (0,3) => 0, (2,1) => 40, (1,2) => 40, (3,1) => 5, (1,3) => 5,
  --    (2,2) => -16, (3,2) => 20, (2,3) => 20, (3,3) => 14})
///

TEST ///
  P2 = convexHull transpose matrix {{0,0},{2,0},{0,2}}
  P3 = convexHull transpose matrix {{0,0},{3,0},{0,3}}
  PP = convexHull transpose matrix {{0,0,0,0},{0,0,1,0},{2,0,1,0},{0,2,1,0},{0,0,0,1},{3,0,0,1},{0,3,0,1}}
  Q = manyPolyhedraToLargeOne({P2, P3})
  Q2 = convexHull (vertices Q)_{0,1,2,3}
  Q3 = convexHull (vertices Q)_{0,4,5,6}
  assert(PP == Q)
  (eZ, eZbar, eZfaces) = computeHodgeDeligne(Q2, FaceInfo => {true, new HashTable, 3})
  assert(eZ === new HashTable from {(0,0) => 6, (1,0) => 0, (0,1) => 0, (2,0) => 0, (1,1) => -3, (0,2) => 0, (2,2) => 1})
  (eZ, eZbar, eZfaces) = computeHodgeDeligne(Q3, FaceInfo => {true, new HashTable, 3})
  assert(eZ === new HashTable from {(0,0) => 9, (1,0) => 1, (0,1) => 1, (2,0) => 0, (1,1) => -3, (0,2) => 0, (2,2) => 1})
  (eZ, eZbar, eZfaces) = computeHodgeDeligne(Q, FaceInfo => {true, new HashTable, 4})
  assert(eZ === new HashTable from {(0,0) => -15, (0,1) => -1, (1,0) => -1, (2,0) => 0, (1,1) => 1, (0,2) => 0, (3,0) => 0, (2,1) => 0, (0,3) => 0, (1,2) => 0,
      (2,2) => -4, (3,3) => 1})
  P = P3
  d = dim P
  for i from 0 to d - 1 do (
      print("codim = " | i);
      fs := faces(i, P);
      Fs := facesAsPolyhedra(i, P);
      for j from 0 to #fs - 1 do (
          print(fs#j#0, ehrhartNumeratorQuicker(Fs#j), ehrhartNumerator(Fs#j));
	  for k from 1 to d - i do (
	      print(k | ": " | #latticePoints(k * Fs#j))
	      );
    	  );
      )
  (eZCI, eZtoric) = computeHodgeDeligneTorusCI({P2, P3})
  assert(eZCI === new HashTable from {(0,0) => 6})
///

TEST ///
  P2 = convexHull transpose matrix {{0,0,0},{1,0,0},{0,1,0}}
  P1 = convexHull transpose matrix {{0,0,0}, {0,0,1}}
  Q1 = P2 + P1
  Q2 = 2*P2 + 3*P1
  (eZCI, eZtoric) = computeHodgeDeligneTorusCI({Q1, Q2})
  assert(eZCI#(0, 1) == -3)
  assert(eZCI#(1, 0) == -3)
  
  V1 = toricProjectiveSpace(1)
  dim(V1)
  V2 = toricProjectiveSpace(2)
  dim(V2)
  V = V1 ** V2
  dim(V)
  rays(V)
  HH^0(V, OO_V(3, 2))
  Q1 = polytope(3 * V_0 + 2 * V_2)
  #latticePoints(Q1)
  Q2 = polytope(V_0 + V_2)
  computeHodgeDeligneTorusCI({Q1, Q2})
  assert(eZCI#(0, 1) == -3)
  assert(eZCI#(1, 0) == -3)
  HH^0(V, OO_V(1, 1))
///


end--

-* Development section *-
restart
debug needsPackage "DanilovKhovanskii"
check "DanilovKhovanskii"

uninstallPackage "DanilovKhovanskii"
restart
installPackage "DanilovKhovanskii"
viewHelp "DanilovKhovanskii"

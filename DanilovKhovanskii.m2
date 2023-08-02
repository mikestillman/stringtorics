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
    "getSparseeZ",
    "eZ2hZ",
    "toeZMatrix",
    "computeHodgeDeligne",
    "computeHodgeDeligneInPToric",
    "computeHodgeDeligneAffineAndTorus",
    "computeHodgeDeligneTorusCI",
    "FaceInfo",
    "EmptyValue"}

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
    Ps := for i from 1 to ceiling(d / 2) list i * P;
    l := prepend(1, for i from 1 to ceiling(d / 2) list (
	#latticePoints(Ps#(i - 1))
	));
    --The coefficients of the Ehrhart numerator of reflexive polytopes
    --is symmetric: c_i = c_{d-i}.
    if isReflexive P then (
	hs1 := prepend(1, for i from 1 to ceiling(d / 2) list (
	    sum for j from max(0, i - d - 1) to i list (-- print(i, j);
	    	(-1)^(i - j) * binomial(d + 1, i - j) * l#(j)
		)
	    )
	    );
	hs2 := for i from ceiling(d / 2) + 1 to d list (
	    hs1#(d - i)
	    );
	return join(hs1, hs2);
	);
    --Instead of using #latticePoints(i * P) for ceiling(d / 2) < i <= d,
    --one can use #interiorLatticePoints(i * P) for 0 < i <= floor (d / 2).
    lint := for i from 1 to floor(d / 2) list (
	#interiorLatticePoints(Ps#(i - 1))
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
    if isSimplicial normalFan P then (--not exactly right
	for key in keys eZ do hZ#key = (-1)^(key#0 + key#1) * eZ#key;
	)
    else (
	error "Not enough information. Can only recover Hodge numbers from the Hodge Deligne polynomial for smooth and/or simplicial varieties."
	);
    hZ
    )

eZ2hZ (Polyhedron, HashTable) := (P, eZ) -> (
    hZ := new MutableHashTable;
    if isSimplicial normalFan P then (--not exactly right
	for key in keys eZ do hZ#key = (-1)^(key#0 + key#1) * eZ#key;
	)
    else (
	error "Not enough information. Can only recover Hodge numbers from the Hodge Deligne polynomial for smooth and/or simplicial varieties."
	);
    new HashTable from hZ
    )

getSparseeZ = method();
getSparseeZ (MutableHashTable, Sequence) := (eZ, pq) -> (
    if eZ#?pq then eZ#pq else 0
    )

getSparseeZ (HashTable, Sequence) := (eZ, pq) -> (
    if eZ#?pq then eZ#pq else 0
    )

toeZMatrix = method(Options => {EmptyValue => 0});
toeZMatrix HashTable := opts -> H -> (
    --assume indexing starts at 0
    nr := 0;
    nc := 0;
    for mn in keys H do (
	if mn#0 > nr then nr = mn#0;
	if mn#1 > nc then nc = mn#1;
	);
    M := matrix for n from 0 to nc list (
	for m from 0 to nr list (
	    if H#?(m, n) then H#(m, n) else opts.EmptyValue
	    )
	);
    M
    )

toeZMatrix MutableHashTable := opts -> H -> (
    --assume indexing starts at 0
    nr := 0;
    nc := 0;
    for mn in keys H do (
	if mn#0 > nr then nr = mn#0;
	if mn#1 > nc then nc = mn#1;
	);
    M := matrix for n from 0 to nc list (
	for m from 0 to nr list (
	    if H#?(m, n) then H#(m, n) else opts.EmptyValue
	    )
	);
    M
    )

vdot = method();
vdot (List, List) := (a, b) -> if #a == #b then (
    sum for i from 0 to #a-1 list a#i*b#i) else (error "Lengths not compatible.")

torusFactor = method();
torusFactor (MutableHashTable, ZZ, ZZ) := (eZ, d, D) -> (--d = dimension of polytope; D = dimension of lattice
    dt := D - d;
    eZproduct := new MutableHashTable from {};
    if dt == 0 then (--print("no torus factors");
	eZproduct = eZ
	)
    else (print("torus factors: " | dt);
	for p from 0 to D - 1 do (
	    for q from 0 to D - 1 do (
		eZproduct#(p, q) = sum for i from 0 to dt list (
		    (-1)^(dt - i) * binomial(dt, i) * getSparseeZ(eZ, (p - i, q - i))
		    );
		);
	    );
	);
    new HashTable from eZproduct
    )

torusFactor (HashTable, ZZ, ZZ) := (eZ, d, D) -> (--d = dimension of polytope; D = dimension of lattice
    dt := D - d;
    eZproduct := new MutableHashTable from {};
    if dt == 0 then (--print("no torus factors");
	eZproduct = eZ
	)
    else (print("torus factors: " | dt);
	for p from 0 to D - 1 do (
	    for q from 0 to D - 1 do (
		eZproduct#(p, q) = sum for i from 0 to dt list (
		    (-1)^(dt - i) * binomial(dt, i) * getSparseeZ(eZ, (p - i, q - i))
		    );
		);
	    );
	);
    new HashTable from eZproduct
    )

faceToCone = method();
--Find the cone of Pfan corresponding to a face of P.
--F is a face of P, Pfan is the normal fan of P.
faceToCone (Polyhedron, Fan) := (F, Pfan) -> (
    rys := rays Pfan;
    vs := vertices F;
    dots := entries (transpose rys * vs);
    Fcone := for r from 0 to #dots - 1 list (
	if all(dots#r, x -> x == -1) then r else continue
	);
    coneFromVData(rys_Fcone)
    )

makeConeTable = method();
--Make a dictionary between faces of P and cones of Pfan.
--Pfan is the normal fan of P.
makeConeTable (Polyhedron, Fan) := (P, Pfan) -> (
    rys := rays Pfan;
    d := dim P;
    new HashTable from for n from 0 to d list (
	Fs := facesAsPolyhedra(d - n, P);
	d - n => for F in Fs list (
	    faceToCone(F, Pfan)
	    )
	)
    )

matchCones = method();
--Find the minimal cone of Pfan containing a given cone of P'fan.
--c1 and c2 are lists of the rays in each cone. clist is a list of the cones in Pfan.
matchCones (List, List, List) := (c1, c2, clist) -> (
    c1 == c2 or (
	Clist := for c in clist list (
	    if isSubset(c1, c) then c else continue
	    );
	isSubset(c1, c2) and not any(clist, c -> c == c1) and all(Clist, c -> isSubset(c2, c))
	)
    )

matchCones (Cone, Cone, List) := (c1, c2, clist) -> (
    c1 == c2 or (
	Clist := for c in clist list (
	    if contains(c, c1) then c else continue
	    );
	contains(c2, c1) and not any(clist, c -> c == c1) and all(Clist, c -> contains(c, c2))
	)
    )

coneToFace = method();
--c is a list of the rays in each cone, Pfanlist is a list of the cones in Pfan,
--Pcones is a HashTable serving as a dictionay between cones in Pfan and faces in P.
coneToFace (List, List, HashTable, Polyhedron) := (c, Pfanlist, Pcones, P) -> (
    c1 := for c2 in Pfanlist do (
	if matchCones(c, c2, Pfanlist) then break c2;
	); print("c1", c1);
    (i, j) := for k in keys Pcones do (
	a := false;
	kl := for l from 0 to #Pcones#k - 1 do (
	    if c1 == Pcones#k#l then (a = true; break (k, l));
	    ); print("kl", kl);
	if a then break kl;
	);
    Fs := facesAsPolyhedra(i, P);
    Fs#j
    )

coneToFace (Cone, List, HashTable, Polyhedron) := (c, Pfanlist, Pcones, P) -> (
    c1 := for c2 in Pfanlist do (
	if matchCones(c, c2, Pfanlist) then break c2;
	); print("c1", rays c1);
    (i, j) := for k in keys Pcones do (
	a := false;
	kl := for l from 0 to #Pcones#k - 1 do (
	    if c1 == Pcones#k#l then (a = true; break (k, l));
	    ); print("kl", kl);
	if a then break kl;
	);
    Fs := facesAsPolyhedra(i, P);
    Fs#j
    )

computeHodgeDeligne = method(Options => {FaceInfo => {true, new HashTable from {}, -1}});
--FaceInfo: first entry = true, if this is the full polytope, and = Sequence containing information about the full poltyope if not; 
--second = data from lower-dimensional faces;
--third = dimension of ambient variety (will usually, but not always, be the number of rows of the vertex matrix)
computeHodgeDeligne Polyhedron := opts -> P -> (
    --determine relevant dimensions (of P and of the ambient space)
    d := dim P;
    D := opts.FaceInfo#2;
    if D == -1 then (
	D = numrows vertices P
	);
    
    --Store information about the normal fan, Pfan, or a simplicial subdivision thereof, P'fan.
    topdim := class opts.FaceInfo#0 === Boolean; print("poly dim = "| d | ", ambient dim = " | D, topdim);
    (Pfan, PfanDim, Pfanfaces, Pfanlist, P'fan, P'fanfaces) := if topdim then (
	Pfan2 := normalFan P;
	PfanDim2 := dim Pfan2;
	Pfanfaces2 := new HashTable from for i from 0 to PfanDim2 list (--not necessarily same order
	    i => facesAsCones(i, Pfan2)
	    );
        Pfanlist2 := flatten for i from 0 to #Pfanfaces2 - 1 list Pfanfaces2#i;
        P'fan2 := if isSimplicial Pfan2 then (print("Pfan is simplicial");
	    Pfan2
    	    )
        else (print("Pfan is not simplicial");
	    VP'2 := if isReflexive P then (
		reflexiveToSimplicialToricVariety P
		)
	    else (
		VP2 := normalToricVariety P;
		makeSimplicial(VP2, Strategy => 1)
		);
	    fan VP'2
	    );
	P'fanfaces2 := new HashTable from for i from 0 to dim P'fan2 list (--not necessarily same order
	    i => facesAsCones(i, P'fan2)
	    );
	(Pfan2, PfanDim2, Pfanfaces2, Pfanlist2, P'fan2, P'fanfaces2)
        )
    else opts.FaceInfo#0;-- print(#Pfan, #Pfanfaces, #Pfanlist, #P'fan, #P'fanfaces);
   
   
    eZ := new MutableHashTable;
    eZbar := new MutableHashTable;
    --For now, assume P is full dimensional (X has no torus factors).
    --If dim P != dim X then (eZ = computeHodgeDeligne(restrict P) * (x*y - 1)^(dim X - dim P)) else     
    --This is implemented at the very end.
    
    --Base cases
    if d == 0 then (--hypersurface is empty
	return (new HashTable from eZ, new HashTable from eZbar, new HashTable from {})
	)
    else if d == 1 then (--print("dim(Z) = 0"); --hypersurface in 1-dimension is a collection of points
	--always simplicial?
	eZ#(0,0) = #latticePoints(P) - 1;
	eZbar#(0,0) = #latticePoints(P) - 1; -- print(eZ#(0,0), eZbar#(0,0));
	eZ = torusFactor(eZ, d, D);
	eZbar = torusFactor(eZbar, d, D);
	return (new HashTable from eZ, new HashTable from eZbar, new HashTable from {})
	);-- print("not 0 or 1");
    
    --Begin by computing eZ of the varieties corresponding to each cone of the [subdivided] normal fan, P'fan.
    --This is known by induction. Build up from lowest dimension, 1.
    eZfaces := new MutableHashTable from opts.FaceInfo#1;
    --print(opts.FaceInfo#1); print(eZfaces);
    --print("eZfaces: " | #eZfaces | " , keys(eZfaces): " | #(keys eZfaces));
    if #(keys eZfaces) == 0 then (--print("no face info");
	Pfaces := faces(P);-- print(Pfaces);
	Pcones := makeConeTable(P, Pfan);--same ordering as Pfaces
	P'fanDim := dim P'fan;
	P'cones := new HashTable from for i from 0 to P'fanDim list (--not necessarily same ordering
	    P'fanDim - i => facesAsCones(i, P'fan)
	    );
	for n from 1 to d - 1 do (
    	    for i from 0 to #(P'cones#(d - n)) - 1 do (print(n, i, P'cones#(d - n)#i);
		eZfaces2 := new HashTable from flatten for k from 1 to n - 1 list (
		    for l from 0 to #(P'cones#(d - k)) - 1 list (-- print(Pfaces#(d - k)#l#0, Pfaces#(d - n)#i#0, isSubset(Pfaces#(d - k)#l#0, Pfaces#(d - n)#i#0));
		        if contains(P'cones#(d - k)#l, P'cones#(d - n)#i) and eZfaces#?(P'cones#(d - k)#l) then P'cones#(d - k)#l => eZfaces#(P'cones#(d - k)#l) else continue
		        )
		    --for k in keys eZfaces list (print(k, Pfaces#(d - n)#i#0);
		    --if isSubset(k, Pfaces#(d - n)#i#0) then (print("yes"); k => eZfaces#k) else continue
		    );-- print("face ready"); print(eZfaces2);
		F := coneToFace(P'cones#(d - n)#i, Pfanlist, Pcones, P);
		e := computeHodgeDeligne(F, FaceInfo => {(Pfan, PfanDim, Pfanfaces, Pfanlist, P'fan, P'fanfaces), eZfaces2, n});
	    	eZfaces#(P'cones#(d - n)#i) = e#0;
		--eZfaces#(Pfaces#(d - n)#i#0) = e#0;--print(e);
		
		--Need to account for torus factors.
		--The dimension of the ambient torus is n. The dimension of the Newton polytope is dim F.
		--eZfaces#(P'cones#(d - n)#i) = torusFactor(e#0, dim F, n); print(dim F, n);
		); print("done " | n);
	    );
	);-- print("faces done");--Hodge-Deligne numbers of the face.
    
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
    eZ = torusFactor(eZ, d, D);
    eZbar = torusFactor(eZbar, d, D);
    if topdim then (print("topdim = true");
        for k in keys(eZfaces) do (
            eZfaces#k = torusFactor(eZfaces#k, d, D);
	    );
        );
    (new HashTable from eZ, new HashTable from eZbar, new HashTable from eZfaces)
    )

computeHodgeDeligne CYPolytope := P -> computeHodgeDeligne(polytope(P, "M"))

computeHodgeDeligne CalabiYauInToric := X -> computeHodgeDeligne(polytope(X, "M"))

--check non-degeneracy
computeHodgeDeligne ToricDivisor := D -> computeHodgeDeligne(polytope(D))

computeHodgeDeligne NormalToricVariety := V -> (
    D := sum for i from 0 to #rays(V) - 1 list (
	V_i
	);
    computeHodgeDeligne(polytope(D))
    )

--For (a hypersurface in) a toric variety that is a subset of the projective normal toric variety that has polytope P.
--Not tested.
computeHodgeDeligneInPToric = method();
computeHodgeDeligneInPToric (Polyhedron, List) := (P, whichFaces) -> (
    (eZ, eZbar, eZfaces) := computeHodgeDeligne(P);
    eZtoric := new MutableHashTable from {};
    for k in keys eZ do (print(k);--maybe not eZ; need to switch to indexing by cones
	eZtoric#k = eZ#k + sum for f in whichFaces list getSparseeZ(eZfaces#f, k);
	);
    new HashTable from eZtoric
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
	    H := {};
	    newP := intersection(P, H);--finish...
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

doc ///
Key
  computeHodgeDeligne
  (computeHodgeDeligne, Polyhedron)
  (computeHodgeDeligne, CYPolytope)
  (computeHodgeDeligne, CalabiYauInToric)
  (computeHodgeDeligne, ToricDivisor)
  (computeHodgeDeligne, NormalToricVariety)
Headline
  compute the Hodge-Deligne polynomial of a hypersurface in a torus.
Usage
  computeHodgeDeligne(P)
Inputs
  P:Polyhedron
    a lattice polytope
Outputs
  :Sequence
    of three HashTables $e_Z$, $e_{\bar{Z}}$, and all of the $e_{Z_\Gamma}$ for $\Gamma \leq P$ a face
Description
  Text
    The Hodge-Deligne polynomial encodes information 
  Example
    P = convexHull matrix {{-1, 4, -1, -1, 0, -1}, {-1, -1, 4, 0, -1, -1}, {-1, -1, -1, 1, 1, 1}}
    latticePoints(P)
    (eZ, eZbar, eZfaces) = computeHodgeDeligne(P)
    eZ
    eZbar
Caveat
SeeAlso
  computeHodgeDeligneInPToric
  computeHodgeDeligneTorusCI
///

doc ///
Key
  stdVector
  (stdVector, ZZ, ZZ)
Headline
  create the i-th standard basis vector in an n-dimensional vector space
Usage
  stdVector(n, i)
Inputs
  n:ZZ
    the dimension of the vector space
  i:ZZ
    the index of the basis vector to create (from 0 to n - 1)
Outputs
  :List
    of 0's in every position except the i-th one, which contains a 1
Description
  Text
    Given a basis of a vector space, any vector in that vector space can be expressed as a linear combination of the basis vectors.
    The coefficients in this linear combination may be arranged into a list, which is also sometimes referred to as the vector.
    When written in this form, the basis vectors themselves take the form of a list with a 1 at the position corresponding to itself and 0's elsewhere.
  Example
    e0 = stdVector(3, 0)
    print transpose matrix {e0}
    e1 = stdVector(3, 1)
    print transpose matrix {e1}
    e2 = stdVector(3, 2)
    print transpose matrix {e2}
Caveat
SeeAlso
  vdot
///

doc ///
Key
  manyMatricesToLargeMatrix
  (manyMatricesToLargeMatrix, List)
Headline
  convert a list of matrices into one large matrix such that the convex hull of the new matrix is the polyhedron desired for computing the Hodge-Deligne numbers of a complete intersection.
Usage
  manyMatricesToLargeMatrix(L)
Inputs
  L:List
    of matrices with the same number of rows
Outputs
  :Matrix
    if there are n matrices, each with m rows, this is the matrix with m + n rows whose first m rows contain the input matrices as adjacent blocks while the last n rows contains blocks of the same size, with the i-th block having 1's in the i-th row and 0's elsewhere.
Description
  Text
    The computation of the Hodge-Deligne numbers of a complete intersection of hypersurfaces in a torus employs an ancilliary hypersurface in a higher-dimensional torus.
    Namely, this hypersurface is defined by a regular function whose Newton polytope can be obtained from the Newton polytopes associated to the original collection of hypersurfaces.
    When the Newton polytopes are expressed as convex hulls of their vertices, and those vertices are put into a matrix as columns, the procedure to construct the Newton polytope of the higher-dimensional hypersurface reduces to the procedure performed by this method.
  Example
    M2 = transpose matrix {{0,0},{2,0},{0,2}}
    M3 = transpose matrix {{0,0},{3,0},{0,3}}
    M = manyMatricesToLargeMatrix({M2, M3})
Caveat
SeeAlso
  manyPolyhedraToLargeMatrix
  manyPolyhedraToLargeOne
  computeHodgeDeligne
  computeHodgeDeligneAffineAndTorus
  computeHodgeDeligneTorusCI
///

doc ///
Key
  manyPolyhedraToLargeMatrix
  (manyPolyhedraToLargeMatrix, List)
Headline
  convert a list of polyhedra into one large matrix such that the convex hull of the matrix is the polyhedron desired for computing the Hodge-Deligne numbers of a complete intersection.
Usage
  manyPolyhedraToLargeMatrix(L)
Inputs
  L:List
    of polyhedra in the same ambient space
Outputs
  :Matrix
    if there are n polyhedra, each sitting in a space of dimension m, this is the matrix with m + n rows whose first m rows contain the matrices of the vertices of the input polyhedra (obtained as vertices(L#i)) as adjacent blocks while the last n rows contains blocks of the same size, with the i-th block having 1's in the i-th row and 0's elsewhere.
Description
  Text
    The computation of the Hodge-Deligne numbers of a complete intersection of hypersurfaces in a torus employs an ancilliary hypersurface in a higher-dimensional torus.
    Namely, this hypersurface is defined by a regular function whose Newton polytope can be obtained from the Newton polytopes associated to the original collection of hypersurfaces.
    Applying convexHull to the output matrix of this method will give the Newton polytope of the higher-dimensional hypersurface.
  Example
    P2 = convexHull transpose matrix {{0,0},{2,0},{0,2}}
    P3 = convexHull transpose matrix {{0,0},{3,0},{0,3}}
    M = manyPolyhedraToLargeMatrix({P2, P3})
Caveat
SeeAlso
  manyMatricesToLargeMatrix
  manyPolyhedraToLargeOne
  computeHodgeDeligne
  computeHodgeDeligneAffineAndTorus
  computeHodgeDeligneTorusCI
///

doc ///
Key
  manyPolyhedraToLargeOne
  (manyPolyhedraToLargeOne, List)
Headline
  convert a list of polyhedra into the higher-dimensional polyhedron desired for computing the Hodge-Deligne numbers of a complete intersection.
Usage
  manyPolyhedraToLargeOne(L)
Inputs
  L:List
    of polyhedra in the same ambient space
Outputs
  :Matrix
    if there are n polyhedra, each sitting in a space of dimension m, this is the matrix with m + n rows whose first m rows contain the matrices of the vertices of the input polyhedra (obtained as vertices(L#i)) as adjacent blocks while the last n rows contains blocks of the same size, with the i-th block having 1's in the i-th row and 0's elsewhere.
Description
  Text
    The computation of the Hodge-Deligne numbers of a complete intersection of hypersurfaces in a torus employs an ancilliary hypersurface in a higher-dimensional torus.
    Namely, this hypersurface is defined by a regular function whose Newton polytope can be obtained from the Newton polytopes associated to the original collection of hypersurfaces.
    Applying convexHull to the output matrix of this method will give the Newton polytope of the higher-dimensional hypersurface.
  Example
    P2 = convexHull transpose matrix {{0,0},{2,0},{0,2}}
    P3 = convexHull transpose matrix {{0,0},{3,0},{0,3}}
    P = manyPolyhedraToLargeOne({P2, P3})
Caveat
SeeAlso
  manyMatricesToLargeMatrix
  manyPolyhedraToLargeMatrix
  computeHodgeDeligne
  computeHodgeDeligneAffineAndTorus
  computeHodgeDeligneTorusCI
///

doc ///
Key
  ehrhartNumerator
  (ehrhartNumerator, Polyhedron)
Headline
  compute the numerator of the rational function expression for the Ehrhart series of a polytope
Usage
  ehrhartNumerator(P)
Inputs
  P:Polyhedron
Outputs
  :List
    the coefficients of the numerator, with the i-th position corresponding to the i-th power of t (beginning from 0)
Description
  Text
    The Ehrhart series of a polytope, P, in which the coefficient of $t^i$ is the number of lattice points in the i-th dilation of P, can be expressed as a rational function with a certain form.
    Namely, one has $Ehr_P(t) = \frac{h^*(t)}{(1-t)^{dim(P)+1}}$, where $h^*(t)$ is a polynomial of degree $dim(P)$.
  Example
    R = QQ[x]
    P = convexHull transpose matrix {{1,1},{1,-1},{-1,1},{-1,-1}}
    eNum = ehrhartNumerator(P)
    h = sum for i from 0 to #eNum - 1 list (
        eNum#i * x^i
        )
    a = 1 + sum for i from 1 to 10 list (
        #latticePoints(i * P) * x^i
        )
    assert((a * (1 - x)^3)%x^11 == h)
Caveat
SeeAlso
  ehrhartNumeratorQuicker
  computeSumqeZ
  computeHodgeDeligne
///

doc ///
Key
  ehrhartNumeratorQuicker
  (ehrhartNumeratorQuicker, Polyhedron)
Headline
  compute the numerator of the rational function expression for the Ehrhart series of a polytope
Usage
  ehrhartNumeratorQuicker(P)
Inputs
  P:Polyhedron
Outputs
  :List
    the coefficients of the numerator, with the i-th position corresponding to the i-th power of t (beginning from 0)
Description
  Text
    The Ehrhart series of a polytope, P, in which the coefficient of $t^i$ is the number of lattice points in the i-th dilation of P, can be expressed as a rational function with a certain form.
    Namely, one has $Ehr_P(t) = \frac{h^*(t)}{(1-t)^{dim(P)+1}}$, where $h^*(t)$ is a polynomial of degree $dim(P)$.
  Example
    R = QQ[x]
    P = convexHull transpose matrix {{1,1},{1,-1},{-1,1},{-1,-1}}
    eNum = ehrhartNumeratorQuicker(P)
    h = sum for i from 0 to #eNum - 1 list (
        eNum#i * x^i
        )
    a = 1 + sum for i from 1 to 10 list (
        #latticePoints(i * P) * x^i
        )
    assert((a * (1 - x)^3)%x^11 == h)
Caveat
SeeAlso
  ehrhartNumerator
  computeSumqeZ
  computeHodgeDeligne
///

doc ///
Key
  computeSumqeZ
  (computeSumqeZ, Polyhedron, List, ZZ)
Headline
  compute $\sum_q e^{p,q}(Z)$ for $Z$ a hypersurface in a torus and any integer $p$
Usage
  computeSumqeZ(P, L, p)
Inputs
  P:Polyhedron
    a lattice polytope
  L:List
    a list representing the Ehrhart numerator of P
  p:ZZ
    the index of the Hodge-Deligne numbers that is held fixed in the sum
Outputs
  :ZZ
    the sum $\sum_q e^{p,q}(Z)$, where $Z$ is a hypersurface defined by a regular function whose polytope is P
Description
  Text
    The value of $\sum_q e^{p,q}(Z)$ for $Z$ a hypersurface in a torus can be computed from Euler-Poincare$\e'$ characteristics of sheaves of differential forms.
    It is used in one of the last step of the Danilov-Khovanskii algorithms to compute $e^{p,q}(Z)$ for $p + q = d - 1$, since for fixed $p$, all other $e^{p,q}(Z)$'s will have been computed and $e^{p,d-1-p}(Z) = \sum_q e^{p,q}(Z) - \sum_{q\neqd-1-p} e^{p,q}(Z)$.
  Example
    P = convexHull transpose matrix {{1,1},{1,-1},{-1,1},{-1,-1}}
    psi = ehrhartNumerator(P)
    computeSumqeZ(P, psi, 0)
    computeSumqeZ(P, psi, 1)
Caveat
SeeAlso
  ehrhartNumerator
  ehrhartNumeratorQuicker
  computeHodgeDeligne
///

doc ///
Key
  eZ2hZ
  (eZ2hZ, Polyhedron, HashTable)
  (eZ2hZ, Polyhedron, MutableHashTable)
Headline
  convert Hodge-Deligne numbers to Hodge numbers for smooth, projective varities
Usage
  eZ2hZ(P, eZ)
Inputs
  P:Polyhedron
  eZ:HashTable
    or @ofClass MutableHashTable@
    the Hodge-Deligne numbers of a variety    
Outputs
  :HashTable
Description
  Text
    If $Z$ is a smooth, projective variety, then $e^{p,q}(Z) = (-1)^{p+q} h^{p,q}(Z)$.
  Example
    PN = convexHull transpose matrix {{1, 0, 0}, {0, 1, 0}, {-1, -1, -2}, {0, 0, 1}, {0, 0, -1}}
    PM = polar PN
    (eZ, eZbar, eZfaces) = computeHodgeDeligne(PM)
    hZ = eZ2hZ(PM, eZ)
    hZbar = eZ2hZ(PM, eZbar)
Caveat
SeeAlso
  toeZMatrix
  computeHodgeDeligne
  computeHodgeDeligneInPToric
  computeHodgeDeligneAffineAndTorus
  computeHodgeDeligneTorusCI
///

doc ///
Key
  getSparseeZ
  (getSparseeZ, HashTable, Sequence)
  (getSparseeZ, MutableHashTable, Sequence)
Headline
  given HashTable(or MutableHashTable) and a key, if the key is in a key-value pair, return the corresponding value; else, return 0
Usage
  getSparseeZ(eZ, pq)
Inputs
  eZ:HashTable
    or @ofClass MutableHashTable@
  pq:Sequence
    a pair of integers
Outputs
  :ZZ
Description
  Text
    Since it is neither possible nor desirable to store the values of $e^{p,q}(Z)$, and all but a finite number of them are 0, this function allows the user to fetch the value of $e^{p,q}(Z)$ for any pair $(p, q)$.
  Example
    P = convexHull transpose matrix {{1,1},{1,-1},{-1,1},{-1,-1}}
    (eZ, eZbar, eZfaces) = computeHodgeDeligne(P)
    getSparseeZ(eZ, (0, 0))
    getSparseeZ(eZ, (1, 1))
    getSparseeZ(eZ, (2, 2))
Caveat
SeeAlso
  computeHodgeDeligne
  computeHodgeDeligneInPToric
  computeHodgeDeligneAffineAndTorus
  computeHodgeDeligneTorusCI
///

doc ///
Key
  toeZMatrix
  (toeZMatrix, HashTable)
  (toeZMatrix, MutableHashTable)
Headline
  convert the eZ HashTable (or MutableHashTable) into a Matrix
Usage
  toeZMatrix(eZ)
Inputs
  eZ:HashTable
    or @ofClass MutableHashTable@
    the Hodge-Deligne numbers of a variety
Outputs
  :Matrix
Description
  Text
    This methods puts the Hodge-Deligne numbers HashTable in a nicer format that resembles a Hodge diamond.
  Example
    PN = convexHull transpose matrix {{1, 0, 0}, {0, 1, 0}, {-1, -1, -2}, {0, 0, 1}, {0, 0, -1}}
    PM = polar PN
    (eZ, eZbar, eZfaces) = computeHodgeDeligne(PM)
    hZ = eZ2hZ(PM, eZ)
    hZbar = eZ2hZ(PM, eZbar)
    toeZMatrix(eZ)
    toeZMatrix(eZbar)
    toeZMatrix(hZ)
    toeZMatrix(hZbar)
Caveat
SeeAlso
  eZ2hZ
  computeHodgeDeligne
  computeHodgeDeligneInPToric
  computeHodgeDeligneAffineAndTorus
  computeHodgeDeligneTorusCI
///

doc ///
Key
  vdot
  (vdot, List, List)
Headline
  compute the vector dot product of two Lists
Usage
  vdot(L1, L2)
Inputs
  L1:List
    a list of elements in a ring, R
  L2:List
    a list of elements in a ring, R
Outputs
  :RingElement
Description
  Text
    The dot product operation takes as input two vectors (stored here as lists) and returns a scalar (in the same ring as the list elements).
  Example
    v1 = {1, 2, 3}
    v2 = {2, 3, 4}
    s = vdot(v1, v2)
Caveat
SeeAlso
  stdVector
///

doc ///
Key
  torusFactor
  (torusFactor, HashTable, ZZ, ZZ)
  (torusFactor, MutableHashTable, ZZ, ZZ)
Headline
  compute the Hodge-Deligne numbers for a variety that is a product of a torus and another variety, whose Hodge-Deligne numbers are given
Usage
  torusFactor(eZ, d, D)
Inputs
  eZ:HashTable
    or @ofClass MutableHashTable@
    the Hodge-Deligne numbers of the variety in the product
  d:ZZ
    the dimension of the polytope
  D:ZZ
    the dimension of the lattice
Outputs
  :HashTable
    or @ofClass MutableHashTable@
    the Hodge-Deligne numbers of the product variety
Description
  Text
    The Hodge-Deligne polynomial of a product of varieties is the product of the Hodge-Deligne polynomial of the varieties.
    In the case of projective normal toric varieties constructed from lattice polytopes, a polytope that is not full-dimensional gives rise to a variety that is a product of a torus and a second variety whose dimension is equal to that of the polytope.
  Example
    P = stdSimplex(2)
    (eZ, eZbar, eZfaces) = computeHodgeDeligne(P)
    Q = convexHull transpose matrix {{0, 0}, {1, 0}, {0, 1}}
    (eZ1, eZbar1, eZfaces1) = computeHodgeDeligne(Q)
    torusFactor(Q, 2, 3)
Caveat
SeeAlso
  computeHodgeDeligne
  computeHodgeDeligneInPToric
  computeHodgeDeligneAffineAndTorus
  computeHodgeDeligneTorusCI
///

doc ///
Key
  computeHodgeDeligneInPToric
  (computeHodgeDeligneInPToric, Polyhedron, List)
Headline
  compute the Hodge-Deligne numbers of a hypersurface in a toric variety that is a subset of the projective normal toric variety that has a given polytope.
Usage
  computeHodgeDeligneInPToric(P, whichFaces)
Inputs
  P:Polyhedron
    a lattice polytope
  whichFaces:List
    the subset of faces of P corresponding to the toric subvariety
Outputs
  :HashTable
    the Hodge-Deligne numbers of the toric subvariety
Description
  Text
    The faces of the polytope P correspond to torus orbits, the union of a a subset of which is a toric variety contained in the toric variety corresponding to P.
    The Hodge-Deligne numbers of this disjoint union is the sum of the Hodge-Deligne numbers of the component tori.
  Example
    P = convexHull transpose matrix {{0, 0}, {1, 0}, {0, 1}}
    eZ = computeHodgeDeligneInPToric(P, {{0, 1}, {1, 2}})
Caveat
SeeAlso
  computeHodgeDeligne
  computeHodgeDeligneAffineAndTorus
  computeHodgeDeligneTorusCI
///

doc ///
Key
  computeHodgeDeligneAffineAndTorus
  (computeHodgeDeligneAffineAndTorus, Polyhedron, ZZ, ZZ)
  (computeHodgeDeligneAffineAndTorus, Matrix, ZZ, ZZ)
Headline
  compute the Hodge-Deligne numbers of a hypersurface in $T^n \times C^r$ defined by a regular function with Newton polytope equal to P
Usage
  computeHodgeDeligneAffineAndTorus(P, n, r)
Inputs
  P:Polyhedron
    a lattice polytope
  n:ZZ
    the dimension of the torus
  r:ZZ
    the dimension of the affine space
Outputs
  :HashTable
    the Hodge-Deligne numbers of the hypersurface
Description
  Text
    $T^n \times C^r$ is a toric variety that can be decomposed into $2^r$ tori.
    The Hodge-Deligne polynomial of a hypersurface in $T^n \times C^r$ is then the sum of the Hodge-Deligne polynomials of the components in each of these tori.
  Example
    P2 = convexHull transpose matrix {{0,0},{2,0},{0,2}}
    P3 = convexHull transpose matrix {{0,0},{3,0},{0,3}}
    M = manyPolyhedraToLargeMatrix({P2, P3})
    n = 2
    r = 2
    computeHodgeDeligneAffineAndTorus(M, n, r, Cheap => true)
Caveat
SeeAlso
  computeHodgeDeligne
  computeHodgeDeligneInPToric
  computeHodgeDeligneTorusCI
///

doc ///
Key
  computeHodgeDeligneTorusCI
  (computeHodgeDeligneTorusCI, List)
Headline
  compute the Hodge-Deligne numbers of a complete intersection of hypersurfaces in a torus
Usage
  computeHodgeDeligneTorusCI(Ps)
Inputs
  Ps:List
    a List of lattice polytopes in the same lattice
Outputs
  :HashTable
    the Hodge-Deligne numbers of a the complete intersection
  :HashTable
    the Hodge-Deligne numbers of a the auxiliary hypersurface used in the calculation
Description
  Text
    The complete intersection of a complete intersection of hypersurfaces in a torus can be computed indirectly from the Hodge-Deligne numbers of a related hypersurfaces in $T^n \times C^r$.
  Example
    P2 = convexHull transpose matrix {{0,0},{2,0},{0,2}}
    P3 = convexHull transpose matrix {{0,0},{3,0},{0,3}}
    (eZCI, eZtoric) = computeHodgeDeligneTorusCI({P2, P3})
Caveat
SeeAlso
  computeHodgeDeligne
  computeHodgeDeligneInPToric
  computeHodgeDeligneAffineAndTorus
///

-* Test section *-
TEST /// -* [insert short title for this test] *-
  R = QQ[x]
  P = convexHull transpose matrix {{1,1},{1,-1},{-1,1},{-1,-1}}
  eNum = ehrhartNumerator(P)
  eNum2 = ehrhartNumeratorQuicker(P)
  h = sum for i from 0 to #eNum - 1 list (
      eNum#i * x^i
      )
  a = 1 + sum for i from 1 to 10 list (
      #latticePoints(i * P) * x^i
      )
  assert((a * (1 - x)^3)%x^11 == h)
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
  P = PP
  d = dim P
  for i from 0 to d - 1 do (
      print("codim = " | i);
      fs := faces(i, P);
      Fs := facesAsPolyhedra(i, P);
      for j from 0 to #fs - 1 do (
          print("vertices = " | toString(fs#j#0), ehrhartNumeratorQuicker(Fs#j), ehrhartNumerator(Fs#j));
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

TEST///
  topes = kreuzerSkarke(3);
  Q = cyPolytope topes_25
  X = makeCY Q
  V = ambient X
  Vfan = fan V
  #rays V
  D = sum for i from 0 to 6 list V_i
  P = polytope D
  isReflexive P
  Pfan = normalFan P
  isSimplicial V
  isSimplicial Vfan
  isSimplicial Pfan
  (eZ, eZbar, eZfaces) = computeHodgeDeligne(P)
  --appears to give h^(1,1) = 2, h^(1,2) = 66
  --compare with topes_50 in previous example with correct result
  
  Pfan = normalFan P;
  PfanDim = dim Pfan;
  Pfanfaces = new HashTable from for i from 0 to PfanDim list (--not necessarily same order
      i => facesAsCones(i, Pfan)
      );
  Pfanlist = flatten for i from 0 to #Pfanfaces - 1 list Pfanfaces#i;
  P'fan = if isSimplicial Pfan then Pfan else (
      VP' = if isReflexive P then reflexiveToSimplicialToricVariety P else (
	  VP = normalToricVariety P;
	  makeSimplicial(VP, Strategy => 1)
          );
      fan VP'
      );
  P'fanfaces = new HashTable from for i from 0 to dim P'fan list (--not necessarily same order
        i => facesAsCones(i, P'fan)
        );
  Pfaces = faces(P);-- print(Pfaces);
  Pcones = makeConeTable(P, Pfan);--same order as Pfaces
  P'fanDim = dim P'fan;
  P'cones = new HashTable from for i from 0 to P'fanDim list (--not necessarily same order
	    P'fanDim - i => facesAsCones(i, P'fan)
	    );

  VP = normalToricVariety P
  --VV = makeSimplicial VP
  VV = reflexiveToSimplicialToricVariety P
  P'fan = fan VV
  Pfanfaces = faces Pfan
  P'fanfaces = faces P'fan
  Pfanlist = flatten for i from 0 to #Pfanfaces - 1 list Pfanfaces#i 
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
  Pfaces := faces(P)
  FacesToCones = new HashTable from flatten for i from 0 to dim P - 1 list (
      for j from 0 to #(Pfaces#i) - 1 list (
	  f := Pfaces#i#j#0;
	  F := (facesAsPolyhedra(i, P))#j;
	  c := faceToCone(F, Pfan); print(f, c);
	  f => c
	  )
      )
  for i from 0 to d - 1 do (
      for j from 0 to #(Pfaces#i) - 1 do (
	  for k from 0 to i - 1 do (
	      for l from 0 to #(Pfaces#k) - 1 do (
		  if isSubset(Pfaces#i#j#0, Pfaces#k#l#0) then (
		      if isSubset(FacesToCones#(Pfaces#k#l#0), FacesToCones#(Pfaces#i#j#0)) then (
			  print(i, j, k, l, Pfaces#i#j#0, Pfaces#k#l#0)
			  )
		      )
		  )
	      )
	  )
      )
  for i from 0 to dim Pfan do (
      for j from 0 to #Pfanfaces#i - 1 do (
	  f := Pfanfaces#i#j;
	  for k from 0 to dim P'fan do (
	      for l from 0 to #P'fanfaces#k - 1 do (
		  f' := P'fanfaces#k#l;
		  if matchCones(f', f, Pfanlist) then print(f', f, d - #f');
		  )
	      )
	  )
      )
  for k in keys Pfanfaces do (
      for l from 0 to #Pfanfaces#k - 1 do (
	  assert(Pfanfaces#k#l == faceToCone(coneToFace(P'fanfaces#k#l, Pfanlist, Pcones, P), Pfan));
	  );
      );
  for k in keys P'fanfaces do (
      for l from 0 to #P'fanfaces#k - 1 do (
	  if not P'fanfaces#k#l == faceToCone(coneToFace(P'fanfaces#k#l, Pfanlist, Pcones, P), Pfan) then print(k, l, P'fanfaces#k#l);
	  );
      );
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

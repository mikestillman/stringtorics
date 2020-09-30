-- also defined: 
--  (1) dim(P,f)
--  (2) genus(P,f)

-- functions from Polyhedra which we use here:
--   latticePoints P -- a different order than 'latticePoints P'
--   vertices P -- a different order than 'vertexList P' !
--   faces(d,P) -- only in faceDimensionHash
--   polar P -- used often
--   dim P -- used often

-- Our plan: stash into a Polyhedron, the information here
protect TCILatticePointList
protect TCIVertexList
protect TCIVertexMatrix
protect TCIFaceDimensionHash
protect TCIInteriorLatticeHash
protect TCILatticePointHash

vertexList = method()
vertexList Polyhedron := (cacheValue symbol TCIVertexList) (P -> (
    verts := vertices P;
    verts = try lift(verts,ZZ) else verts;
    --if liftable(verts,ZZ) then verts = lift(verts,ZZ);
    sort entries transpose verts
    ))

vertexMatrix = method()
vertexMatrix Polyhedron := (cacheValue symbol TCIVertexMatrix) (P -> (
    transpose matrix vertexList P
    ))

faceDimensionHash = method()
if Polyhedra#Options#Version == "1.3" then ( -- version of Polyhedra in M2 <= 1.9.2
  faceDimensionHash Polyhedron := (cacheValue symbol TCIFaceDimensionHash) ((P) -> (
    L := vertexList P;
    vertexHash := hashTable for i from 0 to #L - 1 list (L#i => i);
    hashTable flatten for i from 0 to dim P list for f in faces(dim P-i,P) list (
        -- USING INFO FROM POLYHEDRA
        verts := vertices f;
        --if liftable(verts,ZZ) then verts = lift(verts,ZZ);
        verts = try lift(verts,ZZ) else verts;
        (sort for v in entries transpose verts list vertexHash#v) => i
        )
    ))
) else (
  faceDimensionHash Polyhedron := (cacheValue symbol TCIFaceDimensionHash) ((P) -> (
    L := vertexList P;
    M := vertices P; -- different ordering, possibly, and also a matrix over QQ
    M = try lift(M,ZZ) else M;
    --if liftable(M,ZZ) then M = lift(M,ZZ);
    verticesQ := entries transpose M;
    vertexHash := hashTable for i from 0 to #L - 1 list (L#i => i);
    hashTable flatten for i from 0 to dim P list for f in faces(dim P-i,P) list (
        -- USING INFO FROM POLYHEDRA
        verts := first f;
        if #(last f) > 0 then error "expected a polytope, but received a polyhedron";
        newverts := sort for v in verticesQ_verts list vertexHash#v;
        newverts => i
        )
    ))
)

dim(Polyhedron, List) := (P,f) -> (faceDimensionHash P)#f

faceList = method()
faceList(ZZ,Polyhedron) := (dimF,P) -> (
    H := faceDimensionHash P;
    sort select(keys H, f -> H#f === dimF)
    )
faceList Polyhedron := P -> (
    H := faceDimensionHash P;
    (pairs H)/((k,v) -> (v,k))//sort/last
    )

dualFace = method()
dualFace(Polyhedron, List) := (P,f) -> (
    -- f is a sorted list of integer indices of the vertices, giving a face of P
    -- returns a face of 'polar P', as a list of integer indices of the vertices of 'polar P'
    V1 := vertexMatrix P;
    V2 := transpose vertexMatrix polar P;
    vp := V1_f;
    g := positions(entries (V2 * vp), 
        x -> all(x, x1 -> x1 == -1));
    assert(sort g === g);
    g
    )

-- pts is a list of lattice points in P, or a single lattice point.
-- returns the minimal face of P (as sorted list of vertex indices) containing pts
minimalFace = method()
minimalFace(Polyhedron, List) := (P, pts) -> (
    if all(pts, x -> instance(x,ZZ)) then pts = {pts};
    V2 := transpose vertexMatrix polar P;
    dualfaceverts := positions(entries (V2 * transpose matrix pts), 
        x -> all(x, x1 -> x1 == -1));
    dualFace(polar P, dualfaceverts)
    )

latticePointList = method()
latticePointList Polyhedron := (cacheValue symbol TCILatticePointList) (P -> (
    -- this reorders the lattice points via dimension of smallest face containing them
    Q := P; -- use Q for calling functions in Polyhedra, just for doc...
    lp := latticePoints Q;
    -- the following is our putative list of lattice points.
    L := sort for p in lp list flatten entries lift(p,ZZ);
    -- now we reorder this list, and set:
    --   TCILatticePointList
    --   TCILatticePointHash
    --   TCIInteriorLatticeHash
    L1 := sort for i from 0 to #L-1 list (
        f := minimalFace(P, L#i);
        {dim(P,f), L#i, f}
        );
    --L1 := sort for i from 0 to #L-1 list {faceDim(P,minimalFace(P,L#i)),L#i};
    L2 := L1/(x -> x#1); -- select the actual lattice point
    -- now set the hash tables:
    -- lattice point => index
    H0 := hashTable for i from 0 to #L2-1 list L2#i => i;
    P.cache.TCILatticePointHash = H0;
    -- face => interior lattice points
    H1 := partition(x -> last x, L1); -- partition on minimal face
    H2 := applyPairs(H1, (k,v) -> (k,v/(v1 -> H0#(v1#1))));
    P.cache.TCIInteriorLatticeHash = H2;
    L2
    ))
latticePointList(Polyhedron, List) := (P,f) -> (
    -- f is a sorted list of integer indices of the vertices, giving a face of P
    -- returns the list of indices of lattice point on f.
    V1 := vertexMatrix P;
    LP1 := matrix latticePointList P;
    V2 := vertexMatrix polar P;
    g := dualFace(P,f);
    positions(entries(LP1 * V2_g), x -> all(x, x1 -> x1 == -1))
    )

latticePointHash = method()
latticePointHash Polyhedron := P -> (
    latticePointList P;
    P.cache.TCILatticePointHash
    )

-- f is a sorted list of integer indices of the vertices, giving a face of P
-- returns the list of indices of lattice points in the relative interior of f.
interiorLatticePointList = method()
interiorLatticePointList(Polyhedron, List) := (P,f) -> (
    L := latticePointList P;
    H := P.cache.TCIInteriorLatticeHash;
    if H#?f then H#f else {}
    )

-- number of interior points in the dual face
genus(Polyhedron, List) := (P,f) -> # interiorLatticePointList(polar P, dualFace(P,f))

annotatedFaces = method()
annotatedFaces Polyhedron := (P1) -> (
    P2 := polar P1;
    sort for f in faceList P1 list (
        {dim(P1,f), 
            f, 
            latticePointList(P1,f), 
            # interiorLatticePointList(P1,f), 
            # interiorLatticePointList(P2, dualFace(P1,f))
            }
      )
    )
-- Returns a list for each face of P1 of dimension i:
-- {faceIndices, all lattice pts, #interior lattice pts, #interior lattice pts in dual face of P2}
annotatedFaces(ZZ,Polyhedron) := (i,P1) -> (
    P2 := polar P1;
    sort for f in faceList(i,P1) list (
        {f, 
            latticePointList(P1,f), 
            # interiorLatticePointList(P1,f), 
            # interiorLatticePointList(P2, dualFace(P1,f))
            }
      )
    )


end--


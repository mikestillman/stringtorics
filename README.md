# stringtorics

Here is some functionality we want to include:

Polytope -- preferably just use the class we have?

ReflexivePolytope -- A Polytope, which happens to be reflexive.  Contains
  -- caches for data about the polytope, dual faces, lattice points etc.

CalabiYauToricHypersurface
  -- this is the data of a ReflexivePolytop and a triangulation.
  -- this mostly refers to dimension 3 Calabi-Yaus.  Might contain more later.

  -- functions to include:
  -- a. Picard group
  -- 1. cohomology of line bundles on X, maps between them
  -- 2. intersection ring, "integeral" (uses Picard group).
  -- 3. topological data
  -- 4. GV invariants
  -- 5. Mori cones, effective cones, nef cones ? WANTED: set of isomorphisms coming from birational isomorphisms (i.e. flops)
  -- 6. Flops, topology across a flop
  -- so: extended Kahler cone, ...
  
NefPartition
  -- nef partition of a reflexive polytope
  findNefPartition
  findAllNefPartition(Limit => n)

CalabiYauToricCompleteIntersection
  -- more general than hypersurface.
  -- data: same as hypersurface, together with a nef partition.

CalabiYau3FoldTopology
  -- contains at least 4 pieces of information:
  -- 1. h11, h21
  -- 2. cubic form as a homogeneous polynomial over ZZ
  -- 3. linear form, representing mutl by c2(X), in the same ring.
  -- 4. (not computed, currently)  Includes the torsion part of HH^2(X, ZZ), HH^3(X, ZZ).
  
  -- functionality
  -- 1. determine if 2 topologies are equivalent.  If so, return the change of basis.
  -- 2. given a set of 

-- Eventually: perhaps want other complete intersections, ones that aren't CY.

-- Notes:
-- Would like a general CompleteIntersectionInToric class?
-- We could compute cohomology there, intersectin numbers etc?


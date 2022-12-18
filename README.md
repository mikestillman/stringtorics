# stringtorics

Currently, we have the poorly named types:

  -- CYPolytopeData -- contains essentially the reflexive polytope, but also some computed
    -- data that is meant for ease of creating Calabi Yau hypersurfaces.
    -- Is there a better way to handle this?
  -- CYData
    -- this is a CYPolytopeData, and Triangulation, and also caches other info

These classes are designed to be easily dumped/restored from a string or database file
    (via openDatabase).

TODO: rename these types.  e.g. CYData could be CYHypersurface.
      but I'm not sure about CYPolytopeData.
  Option 1: name it CYPolytope
  Option 2: use Polyhedra, and stash this info into it.  I'm reluctant to do this
    as we want functions that are more easily understood in physics realm.
    e.g. we always have a reflexive polytope, but we always use the dual polytope for
    the fan, and we need the lattice points not interior to facets.
    Reason for reluctance: I don't want Polyhedra to do any computation when I first
    create the polytope (I think).  I'm also concerned about naming conflicts for functions.
      (e.g. allTriangulations).
    So: test whether stashing CYPolytopeData in a Polyhedron is ok?
    Naming conflicts on functions?  e.g. triangulations?

TODO: make sure that dbm files created on apple M1 can be used in linux on intel...

TODO: tests should be more coherrent.
    -- test basics of CYPolytopeData
    -- test basics of CYData
    -- test creation of data bases,use of data bases
    -- test polyhedral functions, including triangulations stuff
    -- test intersection numbers, c2, topology
    -- test GV invariants
    -- code for determining topological equivalence of 2 CYData's.
    -- intersection theory
    -- cohomology of line bundles on V, X.
    -- effective cones, Mori cones, nef cones.
    
Here is some functionality we want to include:

Polytope -- preferably just use the class we have?

ReflexivePolytope -- A Polytope, which happens to be reflexive.  Contains
  -- caches for data about the polytope, dual faces, lattice points etc.

CalabiYauToricHypersurface 
  -- this is the data of a ReflexivePolytope and a triangulation.
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


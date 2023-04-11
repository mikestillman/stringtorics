-- Notes to self:
--  Naomi's code is   /Users/mike/src/stringtorics/naomi-flop-code/Flop Code:
--  Jakob's code is in 

-- Code for constructing potential flops of a CY Hypersurface, or something 
-- constructed from that bvia a sequence of flops.

-- given a curve, gv invariant, topology.  Return the new topology.

-- Determine the types of arguments for:

-- find_nilpotent
-- find_nilpotent_outside_inf
-- is_symmetric_flop
-- find_all_flops

-- what about:
--   toric_curves.compute

--   two_face_triags.all_two_face_triangulations(p)

-- determine if a curve is gv nilpotent
-- determine if a nilpotent ray is "outside the infinity cone"
-- perform a flop

debug needsPackage "StringTorics"


  heftFunction = method()
  heftFunction CalabiYauInToric := X -> (
      mori := hilbertBasisGenerators toricMoriCone(ambient X, basisIndices X);
      sum entries transpose rays dualCone posHull transpose matrix mori
      )

  dot = method()
  dot(List, List) := (v,w) -> (
      if #v =!= #w then error "expected vectors of the same size";
      sum for i from 0 to #v-1 list v#i * w#i
      )

end--
-- Right now, we will do it on an example with h11=3.

-- load this in dir m2-examples.
restart
debug needsPackage "StringTorics"
  -- debug for e.g. hilbertBasisGenerators.
  (Qs, Xs) = readCYDatabase("cys-ntfe-h11-3.dbm");
  #Qs
  #Xs

-- Step 1. Find all gv invariants up to a degree bound.
  X = Xs#(50,0)
  GV = gvInvariants(X, DegreeLimit => 20)
  GVcone = (keys GV)/toList//matrix//transpose//posHull
  rays GVcone
  -- Question: what is the degree limit method?
  

  deglimit = 20
  hf = heftFunction X
  assert(hf == {1,1,3})
  (keys GV)/(v -> dot(toList v, hf))
  allGV = (keys GV)/toList//set
  select(keys allGV, v -> dot(v, hf) <= deglimit // 2)
  allrays = unique for v in keys allGV list (
      c := gcd v;
      for v1 in v list v1//c
      )
  goodray = x -> (
      d := dot(x, hf);
      topval := floor(deglimit/d);
      all(splice{1..topval}, i -> member(i * x, allGV))
      )
  H = partition(goodray, allrays);
  rays posHull transpose matrix H#true
  matrix{hf} * oo
  matrix{hf} * transpose matrix H#false
  for x in oo list (
      d := dot(x, hf);
      topval := floor(deglimit/d);
      if not all(splice{1..topval}, i -> member(i * x, allGV))
        then continue
        else x
      )  
  rays posHull transpose matrix oo

  for k in keys allGV list (
      -- we will select all of the ones
      )

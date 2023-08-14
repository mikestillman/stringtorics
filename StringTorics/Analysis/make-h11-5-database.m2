-- In this file, we create the h^(1,1)=5 database of all of the hypersurfaces 
-- in (simplicial resolutions of) Fano toric 4-folds which are smooth CY3-folds.

  -- h11=5 database construction, 19 June 2023.
  restart
  debug needsPackage "StringTorics"
  DB5 = "../Databases/cys-ntfe-h11-5.dbm"
  topes = kreuzerSkarke(5, Limit => 20000); -- 4990 of these
  assert(#topes == 4990)
  elapsedTime createCYDatabase(DB5, topes)
  -- Let's find which are not favorable, not torsion free.
  -- This took 6300 seconds (a bit less than 2 hours).
  topes2 = drop(topes,4873)
  elapsedTime createCYDatabase(DB5, topes2)
  -- OK made it to here.
    
  elapsedTime Qs = readCYPolytopes DB5; -- 3.6 seconds to read these in.
  assert(sort keys Qs == splice{0..4989})

  -- This takes several hours on heaviside.math.cornell.edu
  elapsedTime for k in keys Qs do (    
      << "---- doing k = " << k << endl;
      Q := Qs#k;
      elapsedTime addToCYDatabase(DB5, Q, NTFE => true);
      );

-- Checking the created database
restart
  debug needsPackage "StringTorics"
  DB5 = "../Databases/cys-ntfe-h11-5.dbm"
  R = ZZ[a,b,c,d,e]
  RZ = R
  RQ = QQ (monoid R);
  elapsedTime (Qs, Xs) = readCYDatabase(DB5, Ring => R);
  assert(#keys Qs == 4990)
  assert(#keys Xs == 13635)

  -- torsion in the class group
  elapsedTime torsionQs = for k in keys Qs list (
      istor := prune coker matrix rays Qs#k;
      if not isFreeModule istor then k else continue
      )
  -- torsion in HH^2.
  elapsedTime nontrivialPis = for k in sort keys Qs list (
      if gens gb matrix transpose rays Qs#k != 1 then k else continue
      )
  select(sort keys Xs, lab -> member(first lab, torsionQs)) === {(1,0),(2,0)}-- one for each polytope
  assert(nontrivialPis === {1,2}) -- matches the torsion in the ambient toric variety class group.

  elapsedTime nonfavorableQs = for k in keys Qs list (
      if not isFavorable Qs#k then k else continue
      )
  assert(#nonfavorableQs == 93) -- 1,2 seem to be both non-favorable and torsion.
  assert(not isFavorable Qs#1)
  assert(not isFavorable Qs#2) 
  nonfavorableQs

-- XXX

  for k in nonfavorables list (Q := Qs#k; findTwoFaceInteriorDivisors Q)

  findBasis = method()  
  findBasis CYPolytope := List => Q -> (
    -- first find 2-face interiors with g>0.
    nonfavAndGenera := findTwoFaceInteriorDivisors Q;
    nonfavs := for f in nonfavAndGenera list f#1#0;
    allindices := splice {0..#rays Q - 1};
    otherindices := sort toList (set allindices - set nonfavs);
    M := transpose matrix rays Q;
    Z := transpose LLL syz M;
    subs := for f in subsets(otherindices, numrows Z - #nonfavs) list (f | nonfavs);
    for g in subs do if abs det(Z_g) == 1 then return g;
    "rats"
    )

  findBasis 
  Q = Qs#429
  nonfavAndGenera = findTwoFaceInteriorDivisors Q
  nonfavs = for f in nonfavAndGenera list (f#1#0#0)
    allindices = splice {0..#rays Q - 1};
    otherindices = sort toList (set allindices - set nonfavs)
    M = transpose matrix rays Q
    Z = transpose LLL syz M
    subs := for f in subsets(otherindices, numrows Z - #nonfavs) list (f | nonfavs);
    for g in subs list if abs det(Z_g) == 1 then g else continue
  
-- One option: is in CYPolytope, compute also the non-favorable divisors keeping: genus, which are on which.
--  one way: "nonfavorable" => {list of 2-face info}, each 2-face info is: {genus, list of divisors interior to that same 2-face}  
--  if we do this, we have to make sure it is written to the database too, I guess, or just recompute it.
--  One option: nonfavorables: index of 2-face (from annotated faces), genus??  Need to figure out most convenient format.
-- basis indices code needed
-- 1. findBasisIndices: also gives the non-favorable divisors as (i, j, g), 
--      maybe instead: (i,j), not the g, as we can collect that from the "nonfavorable" info.
--      i=divisor index, 0 <= j <= g, g = genus of the correesponding 2-face.
-- 2. findToricDivisors: no change
-- 3. toBasisIntersectionNumbers: needs bit more information, e.g. which divisors are on which 2-face.
--      and the genera of these 2-faces.
--      For instance, this function could take the list of non-favorables as an input.
--      to compute {i,j,k}:
--      if i,j,k are not non-favorable, then we have already computed this.
--      if more than one is non-favorable, and they come from different 2-faces, the number is 0.
--      otherwise, use {divisor i, divisor j, divisor k}//(1+g).
-- 4. computeC2
--     needs to be changed.
--     in the loop, if a is a non-favorable divisor, then the intersection number {a,i,j} is {divisor(a),i,j}//(g+1)
--     this might be the only change...
-- 5. in computing basis indices:
--     if possible, we would like to have the non-favorable divisors be part of the toric basis of Cl(V).
-- I think these are possibly the only changes (except reading/dumping the information).

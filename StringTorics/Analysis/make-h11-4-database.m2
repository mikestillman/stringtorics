-- In this file, we create the h^(1,1)=4 database of all of the hypersurfaces 
-- in (simplicial resolutions of) Fano toric 4-folds which are smooth CY3-folds.

  -- h11=4 database construction, 25 Aug 2023.
  restart
  needsPackage "StringTorics"
  DBNAME = "../Databases/cys-ntfe-h11-4.dbm"
  topes = kreuzerSkarke(4, Limit => 5000); -- 1197 of these
  assert(#topes == 1197)
  elapsedTime addToCYDatabase(DBNAME, topes) -- 1485 seconds (includes triangs, autos, moriconecap, all CYs)

----  elapsedTime createCYDatabase(DBNAME, topes) -- 1468 seconds, includes all Xs.

  -- Let's find which are not favorable, not torsion free.
  restart
  needsPackage "StringTorics"
  DBNAME = "../Databases/cys-ntfe-h11-4.dbm"
  (Qs, Xs) = readCYDatabase DBNAME; 

  assert(sort keys Qs == splice{0..1196})
  #sort keys Xs == 1774
  
  nonfavorableQs = for k in keys Qs list (
      if not isFavorable Qs#k then k else continue
      )
  favorablesXs = sort for k in keys Xs list (
      if not isFavorable Qs#(first k) then continue else k
      )
  assert(#nonfavorableQs == 12)
  assert(#nonfavorableXs == 14)
  assert(nonfavorableQs == {796, 800, 803, 1059, 1060, 1064, 1065, 1134, 1135, 1151, 1153, 1155})
  nonfavorableXs == {
      (796, 0), -- potentially same as (798,1)? (which is equiv to 5 others) SAME YES the same!
      (796, 1), -- potentially same as (795,1)? (which is equiv to 5 others) -- seemingly not same
      (800, 0), -- potentially same as (844,0)? (which is equiv to 2 others) SAME 
      (803, 0), -- potentially same as (807,0)? (which is equiv to 5 others) SAME YES the same!
      (1059, 0), (1060, 0), (1064, 0), (1065, 0), -- these 4 have same c2, cubic (exactly same) 
      -- possibly same as (1071,0)? (which is itself equiv to 10 others) SAME
      (1134, 0), -- potentially same as (1136,0)? SAME
      (1135, 0), -- potentially same as (1139,0)? SAME
      (1151, 0), (1153, 0), -- potentially same each other, and as (1156,0)? (which is equiv to a bunch of others).  hessian doesn't factor...
      (1155, 0), -- UNIQUE this one is unique in its invariants bucket. bu hessian does factor.
      (1155, 1)} -- potentially same as (1168,0). hessian doesn't factor

  torsionQs = for k in keys Qs list (
      istor := prune coker matrix rays Qs#k;
      if not isFreeModule istor then k else continue
      )
  torsionXs = sort select(keys Xs, k -> member(first k, torsionQs))
  assert(#torsionQs == 6)
  assert(torsionQs == {0, 3, 4, 5, 12, 15})
  assert(torsionXs == {(0, 0), (3, 0), (4, 0), (5, 0), (12, 0), (15, 0)})

  -- This is with the (fixed) toric mori cone cap's.
  moricones = for lab in sort keys Xs list if not isFavorable Xs#lab then continue else lab => toricMoriConeCap Xs#lab;
  tally (moricones/(x -> #x#1))
  assert(
        (tally for m in moricones list if m === null then continue else #m#1)
        ===
        new Tally from {4 => 1104, 5 => 496, 6 => 151, 7 => 9} -- check this with Andreas.
        )
  netList (moricones/(x -> prepend(x#0, x#1)))

  


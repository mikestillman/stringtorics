-- In this file, we create the h^(1,1)=2 database of all of the hypersurfaces 
-- in (simplicial resolutions of) Fano toric 4-folds which are smooth CY3-folds.

-------------------------
-- The actual creation --
-------------------------
  restart
  needsPackage "StringTorics"
  topes = kreuzerSkarke(2, Limit => 1000);
  assert(#topes == 36)
  DBNAME = "../Databases/cys-ntfe-h11-2.dbm"
  elapsedTime addToCYDatabase(DBNAME, topes) -- 38 seconds

----------------------------------------------------
-- Checking some basic data of examples in the DB --
-- Testing results ---------------------------------
----------------------------------------------------
  restart
  needsPackage "StringTorics"
  DBNAME = "../Databases/cys-ntfe-h11-2.dbm"
  RZ = ZZ[a,b]
  RQ = QQ (monoid RZ);
  (Qs, Xs) = readCYDatabase(DBNAME, Ring => RZ);
  assert(sort keys Qs == splice{0..35})
  assert(# keys Xs == 39)

  -- automorphisms?
  autsizes = tally for lab in sort keys Qs list #automorphisms Qs#lab
  assert (
      autsizes 
      ===
      new Tally from {2 => 5, 4 => 5, 6 => 12, 8 => 1, 12 => 6, 16 => 1, 24 => 3, 36 => 1, 48 => 1, 72 => 1}
      )

  netList for lab in sort keys Xs list toricMoriConeCap Xs#lab

  -- torsions?
  torsionQs = for k in keys Qs list (
      istor := prune coker matrix rays Qs#k;
      if not isFreeModule istor then k else continue
      )
  assert(torsionQs === {0, 1})
  
  nonfavorableQs = for k in keys Qs list (
      if not isFavorable Qs#k then k else continue
      )
  assert(nonfavorableQs === {})

---------------------------------
-- Alternate uses and examples --
---------------------------------  
TEST ///
-- some tests, alternate ways to construct
  DBNAME = "../Databases/test3-cys-ntfe-h11-2.dbm" -- true below is the default, as is NTFE => true.
  elapsedTime addToCYDatabase(DBNAME, topes, "CYs" => true) -- 47 seconds
  F= openDatabase DBNAME  
  sort keys F
///

TEST ///
  -- Here we add in only polytope information, not CY info.
  restart
  needsPackage "StringTorics"
  DBNAME = "../Databases/test4-cys-ntfe-h11-2.dbm"
  topes = kreuzerSkarke(2, Limit => 1000);
  assert(#topes == 36)
  elapsedTime addToCYDatabase(DBNAME, topes, "CYs" => false) -- 42 seconds

  F = openDatabase DBNAME
  sort keys F  
  close F

  Qs = readCYPolytopes DBNAME;
  assert instance(Qs, HashTable)
  assert (sort keys Qs === splice {0..35})

  -- Now create the CY's
  for lab in sort keys Qs do addToCYDatabase(DBNAME, Qs#lab)

  F = openDatabase DBNAME
  sort keys F  
  close F

  RZ = ZZ[a,b]
  (Qs, Xs) = readCYDatabase(DBNAME, Ring => RZ);
///  
  
  
TEST ///
  DBNAME = "../Databases/test5-cys-ntfe-h11-2.dbm"
  elapsedTime createCYDatabase(DBNAME, topes) -- 46 seconds

  DBNAME = "../Databases/test6-cys-ntfe-h11-2.dbm"
  elapsedTime addToCYDatabase(DBNAME, topes_{10..20}) -- xx seconds
  F= openDatabase DBNAME  
  sort keys F  
///

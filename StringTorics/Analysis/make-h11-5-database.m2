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
  elapsedTime (Qs, Xs) = readCYDatabase(DB5, Ring => R); -- 14 seconds on Apple M1 laptop
  assert(#keys Qs == 4990)
  assert(#keys Xs == 13635)

  -- torsion in the class group
  elapsedTime torsionQs = for k in keys Qs list (
      istor := prune coker matrix rays Qs#k;
      if not isFreeModule istor then k else continue
      )  -- 8.6 sec on apple M1 laptop.
  -- torsion in HH^2.
  nontrivialPis = for k in sort keys Qs list (
      if gens gb matrix transpose rays Qs#k != 1 then k else continue
      )
  select(sort keys Xs, lab -> member(first lab, torsionQs)) === {(1,0),(2,0)}-- one for each polytope
  assert(nontrivialPis === {1,2}) -- matches the torsion in the ambient toric variety class group.

  nonfavorableQs = for k in keys Qs list (
      if not isFavorable Qs#k then k else continue
      )
  assert(#nonfavorableQs == 93) -- 1,2 seem to be both non-favorable and torsion.
  assert(not isFavorable Qs#1)
  assert(not isFavorable Qs#2) 
  nonfavorableQs


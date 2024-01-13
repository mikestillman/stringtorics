-- In this file, we create the h^(1,1)=5 database of all of the hypersurfaces 
-- in (simplicial resolutions of) Fano toric 4-folds which are smooth CY3-folds.

  -- h11=5 database construction, 22 Aug 2023.
  -- This crashes because of mutable hashcode overflow.
  -- Total time though is approx. 5-6 hours to construct.
  restart
  debug needsPackage "StringTorics"
  DBNAME = "../Databases/cys-ntfe-h11-5.dbm"
  topes = kreuzerSkarke(5, Limit => 20000); -- 4990 of these
  assert(#topes == 4990)
  elapsedTime addToCYDatabase(DBNAME, topes)
  
  restart
  debug needsPackage "StringTorics"
  DBNAME = "../Databases/cys-ntfe-h11-5.dbm"
  topes = kreuzerSkarke(5, Limit => 20000); -- 4990 of these
  assert(#topes == 4990)
  elapsedTime addToCYDatabase(DBNAME, topes_{2709..#topes-1})


  restart
  debug needsPackage "StringTorics"
  DBNAME = "../Databases/cys-ntfe-h11-5.dbm"
  DBNAME1 = "../Databases/cys-ntfe-h11-5-part1.dbm"
  DBNAME2 = "../Databases/cys-ntfe-h11-5-part2.dbm"
  DBNAME3 = "../Databases/cys-ntfe-h11-5-part3.dbm"
  DBNAME4 = "../Databases/cys-ntfe-h11-5-part4.dbm"
  DBNAME5 = "../Databases/cys-ntfe-h11-5-part5.dbm"
    
  topes = kreuzerSkarke(5, Limit => 20000); -- 4990 of these
  assert(#topes == 4990)
  elapsedTime addToCYDatabase(DBNAME1, topes_{0..1000}) -- running in m2-2
  elapsedTime addToCYDatabase(DBNAME2, topes_{1001..2000}) -- running in m2-3
  elapsedTime addToCYDatabase(DBNAME3, topes_{2001..3000}) -- running in m2-4
  elapsedTime addToCYDatabase(DBNAME4, topes_{3001..4000}) -- running in m2-5
  elapsedTime addToCYDatabase(DBNAME5, topes_{4001..4989}) -- running in m2-6

  combineCYDatabases{DBNAME, DBNAME1, DBNAME2, DBNAME3, DBNAME4, DBNAME5}

  restart
  debug needsPackage "StringTorics"
  DBNAME = "../Databases/cys-ntfe-new-h11-5.dbm"
  DBNAME1 = "../Databases/cys-ntfe-new-h11-5-part1.dbm"
  DBNAME2 = "../Databases/cys-ntfe-new-h11-5-part2.dbm"
  DBNAME3 = "../Databases/cys-ntfe-new-h11-5-part3.dbm"
  DBNAME4 = "../Databases/cys-ntfe-new-h11-5-part4.dbm"
  DBNAME5 = "../Databases/cys-ntfe-new-h11-5-part5.dbm"
    
  topes = kreuzerSkarke(5, Limit => 20000); -- 4990 of these
  assert(#topes == 4990)
  elapsedTime addToCYDatabase(DBNAME1, topes_{0..1000}) -- running in m2-2
  elapsedTime addToCYDatabase(DBNAME2, topes_{1001..2000}) -- running in m2-3
  elapsedTime addToCYDatabase(DBNAME3, topes_{2001..3000}) -- running in m2-4
  elapsedTime addToCYDatabase(DBNAME4, topes_{3001..4000}) -- running in m2-5
  elapsedTime addToCYDatabase(DBNAME5, topes_{4001..4989}) -- running in m2-6

  combineCYDatabases{DBNAME, DBNAME1, DBNAME2, DBNAME3, DBNAME4, DBNAME5}

---------------------------------
-- Check the database somewhat --
---------------------------------
  restart
  debug needsPackage "StringTorics"
  DBNAME = "../Databases/cys-ntfe-h11-5.dbm"
  R = ZZ[a,b,c,d,e]
  RZ = R
  RQ = QQ (monoid R);
  elapsedTime (Qs, Xs) = readCYDatabase(DBNAME, Ring => R); -- 15 seconds on Apple M1 laptop

-- NOTYET  -- This is with the (fixed) toric mori cone cap's.
  moricones = for lab in sort keys Xs list if not isFavorable Xs#lab then continue else lab => toricMoriConeCap Xs#lab;
  tally (moricones/(x -> #x#1))
  netList (moricones/(x -> prepend(x#0, x#1)))
  assert(
        (tally for m in moricones list if m === null then continue else #m)
        ===
        new Tally from {4 => 1104, 5 => 496, 6 => 151, 7 => 9} -- check this with Andreas.
        )

  netList moricones
  
  -- this is the one with the likely faulty setToricMoriConeCap
  new Tally from {5 => 4574, 6 => 3459, 7 => 2178, 8 => 764, 9 => 333, 10 => 128, 11 => 57, 12 => 105, 13 => 40, 14 => 34, 15 => 7, 16 => 18, 17 => 6, 18 => 4, 20 => 3, 21 => 1, 23 => 1, 54 => 1}

-- code to fix the toric mori cones in the database.
debug StringTorics
F = openDatabaseOut DBNAME
for lab in sort keys Xs do (
    X := Xs#lab;
    if not isFavorable X then continue;
    elapsedTime setToricMoriConeCap X;
    print netList prepend(lab, toricMoriConeCap X);
    F#(toString lab) = dump X;
    )
close F


------------------------------------
-- Checking the created database ---
------------------------------------
  restart
  debug needsPackage "StringTorics"
  DBNAME = "../Databases/cys-ntfe-h11-5.dbm"
  R = ZZ[a,b,c,d,e]
  RZ = R
  RQ = QQ (monoid R);
  elapsedTime (Qs, Xs) = readCYDatabase(DBNAME, Ring => R); -- 12 seconds on Apple M1 laptop

  assert(#keys Qs == 4990)
  assert(#keys Xs == 13635)
  assert(#keys Xs == 11847)

  elapsedTime Qs1 = readCYPolytopes DBNAME; -- 4 seconds to read these in.
  assert(sort keys Qs == splice{0..4989})

  -- torsion in the class group
  elapsedTime torsionQs = for k in keys Qs list (
      istor := prune coker matrix rays Qs#k;
      if not isFreeModule istor then k else continue
      )  -- 9 sec on apple M1 laptop.
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

  nonfavorableXs = for k in sort keys Xs list (
      if not isFavorable Qs#(first k) then k else continue
      )
  assert(#nonfavorableXs == 134)

  -- Now check automorphisms
  assert(
    (tally for lab in sort keys Qs list # automorphisms Qs#lab)
    ===
    new Tally from {
        1 => 1808, 
        2 => 2424, 
        4 => 508, 
        6 => 118, 
        8 => 81, 
        12 => 29, 
        16 => 12, 
        24 => 6, 
        48 => 2, 
        72 => 2}
    )

  -- Now check sizes of toric mori cone caps (I've asked Andreas to check these numbers).
  assert(
      (tally for lab in sort keys Xs list if isFavorable Xs#lab then #(toricMoriConeCap Xs#lab) else continue)
      ===
      new Tally from {
          5 => 4574,
          6 => 3459,
          7 => 2178,
          8 => 764,
          9 => 333,
          10 => 128,
          11 => 57,
          12 => 105,
          13 => 40,
          14 => 34,
          15 => 7,
          16 => 18,
          17 => 6,
          18 => 4,
          20 => 3,
          21 => 1,
          23 => 1,
          54 => 1
          }
      )

  heft Xs#(3,0)
  toricMoriConeCap Xs#(3,0)
  rays dualCone posHull transpose matrix oo

-- In this file, we create the h^(1,1)=3 database of all of the hypersurfaces 
-- in (simplicial resolutions of) Fano toric 4-folds which are smooth CY3-folds.

  restart
  needsPackage "StringTorics"
  topes = kreuzerSkarke(3, Limit => 1000);
  assert(#topes == 244)
  DBNAME = "../Databases/cys-ntfe-h11-3.dbm"
  elapsedTime addToCYDatabase(DBNAME, topes) -- 160 seconds, and includes all Xs's, and triangulations, toric mori cones, autos.  Why did it go from 144 to 160??
  
-- Query the results to make sure it seems correct.  
  restart
  needsPackage "StringTorics"
  DBNAME = "../Databases/cys-ntfe-h11-3.dbm"
  R = ZZ[a,b,c]
  RQ = QQ (monoid R);
  (Qs, Xs) = readCYDatabase(DBNAME, Ring => R);

  assert(#Qs == 244)
  assert(#Xs == 275)
  assert(first sort keys Xs === (0,0)) -- make sure at leat one of these has the expected label.

  -- Now we peruse some entries to make sure they have what we expect.  
  -- TODO: make this in to a test, with assertions?
  F = openDatabase DBNAME
  F#"45"
  F#"(45,0)"
  F#"(45,1)"
  F#"(45,2)"
  F#"53"
  F#"(53,0)"
  F#"(53,1)"
  close F

  -- which are favorable (all but 1)
  assert(
    (tally for lab in sort keys Qs list isFavorable Qs#lab)
    === 
    new Tally from {false => 1, true => 243}
    )

  -- check numbers of automorphisms
  assert(
      (tally for lab in sort keys Qs list # automorphisms Qs#lab)
      ===
      new Tally from {1 => 15, 2 => 89, 4 => 52, 6 => 46, 8 => 19, 12 => 15, 16 => 4, 48 => 4}
  )

  assert(
      (tally for lab in sort keys Xs list if not isFavorable Xs#lab then continue else # toricMoriConeCap Xs#lab)
      ===
      new Tally from {3 => 236, 4 => 38} -- TODO: check with Andreas
      --new Tally from {3 => 229, 4 => 44, 6 => 1}  -- TODO: check with Andreas
      )

  moricones = for lab in sort keys Xs list if not isFavorable Xs#lab then continue else prepend(lab, toricMoriConeCap Xs#lab);
  tally for m in moricones list if m === null then continue else #m
  netList moricones
  

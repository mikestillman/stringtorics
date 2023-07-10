-- In this file, we create the h^(1,1)=4 database of all of the hypersurfaces 
-- in (simplicial resolutions of) Fano toric 4-folds which are smooth CY3-folds.

  -- h11=4 database construction, 19 June 2023.
  restart
  needsPackage "StringTorics"
  DB4 = "../Databases/cys-ntfe-h11-4.dbm"
  topes = kreuzerSkarke(4, Limit => 5000); -- 1197 of these
  assert(#topes == 1197)
  elapsedTime createCYDatabase(DB4, topes)
  -- Let's find which are not favorable, not torsion free.
  
  Qs = readCYPolytopes DB4;
  assert(sort keys Qs == splice{0..1196})

  elapsedTime for k in keys Qs do (    
      << "---- doing k = " << k << endl;
      Q := Qs#k;
      elapsedTime addToCYDatabase(DB4, Q, NTFE => true);
      );


  nonfavorables = for k in sort keys Qs list if not isFavorable Qs#k then k else continue
  torsions = for k in keys Qs list (
      istor := prune coker matrix rays Qs#k != ZZ^4;
      if istor then k else continue
      )
  for k in torsions list (
      prune coker matrix rays Qs#k
      )

  assert(
      nonfavorables 
      == {796, 800, 803, 1059, 1060, 
          1064, 1065, 1134, 1135, 1151, 1153, 1155}
  assert(
      torsions
      == 
      {0, 3, 4, 5, 12, 15, 796, 800, 803, 1059, 1060, 
       1064, 1065, 1134, 1135, 1151, 1153, 1155}
      )
  -- There are 6 that are torsion, but still favorable...
  -- Is this correct?  What do I have to change to handle these?
  
  -- Now let's add in all of the NTFE triangulations.
  elapsedTime for k in keys Qs do (    
      if member(k, torsions) then continue;
      << "---- doing k = " << k << endl;
      Q := Qs#k;
      elapsedTime addToCYDatabase(DB4, Q, NTFE => true);
      );


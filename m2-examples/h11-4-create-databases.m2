-- Let's create all of the h11=4 examples.  Should we place them into one file, or into 
-- one per h12?  Currently, I have them in one file per each h12.

-- Let's create one big database file for h11=4
  restart
  needsPackage "StringTorics"
  topes = kreuzerSkarke(4, Limit => 10000);
  assert(#topes == 1197)
  elapsedTime createCYDatabase("foo-cys-ntfe-h11-4.dbm", topes) -- 

  -- we also now want to check that what we have constructed matches
  -- what we expect it to construct.  todo: this should be a test in
  -- the package.  The construction functions need to be documented
  -- too!
  restart
  needsPackage "StringTorics"
  F = openDatabase "foo-cys-ntfe-h11-4.dbm";
  assert(#(keys F) == 1197)
  kF = select((keys F)/value//sort, k -> instance(k, ZZ))
  assert(kF === toList(0..1196))
  F#"10"
  close F
  elapsedTime Qs = readCYPolytopes "foo-cys-ntfe-h11-4.dbm";
  for k from 0 to 1196 do (
      Q := Qs#k;
      assert instance(Q, CYPolytope);
      assert(Q#?"rays");
      assert(Q#?"face dimensions");
      assert(Q.cache#?"favorable");
      assert(Q.cache#?"annotated faces");
      assert(Q.cache#?"basis indices");
      assert(Q.cache#?"glsm");
      assert(Q.cache#?"h11");
      assert(Q.cache#?"h21");
      assert(label Q == k);
      )
  nonfavorables = select(toList(0..1196), i -> not isFavorable Qs#i)
  -- non-favorables appear to be:
  assert(nonfavorables == {796, 800, 803, 1059, 1060, 1064, 1065, 1134, 1135, 1151, 1153, 1155})
  
  h12s = toList(0..1196)/(k -> hh^(1,2) Qs#k)//unique//sort
  assert(#h12s == 87)


  
  -- Now let's add in all the CY's total, including all triangulations
  -- which are distinct when restricted to 2-faces (NTFE => true says don't use all triangulations).
  Qs = readCYPolytopes "foo2-cys-ntfe-h11-3.dbm";
  elapsedTime for Q in values Qs do addToCYDatabase("foo2-cys-ntfe-h11-3.dbm", Q, NTFE => true); -- 11 seconds

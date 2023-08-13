-- In this file, we create the h^(1,1)=2 database of all of the hypersurfaces 
-- in (simplicial resolutions of) Fano toric 4-folds which are smooth CY3-folds.

  restart
  needsPackage "StringTorics"
  topes = kreuzerSkarke(2, Limit => 1000);
  assert(#topes == 36)
  DBNAME = "../Databases/cys-ntfe-h11-2.dbm"
  elapsedTime createCYDatabase(DBNAME, topes) -- 

  -- Now let's add in all the CY's total, including all triangulations.
  Qs = readCYPolytopes DBNAME;
  elapsedTime for Q in values Qs do (
    addToCYDatabase(DBNAME, Q, NTFE => true);
    ) -- 1.6 sec

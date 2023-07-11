-- In this file, we create the h^(1,1)=3 database of all of the hypersurfaces 
-- in (simplicial resolutions of) Fano toric 4-folds which are smooth CY3-folds.

  restart
  needsPackage "StringTorics"
  topes = kreuzerSkarke(3, Limit => 1000);
  assert(#topes == 244)
  elapsedTime createCYDatabase("../Databases/cys-ntfe-h11-3.dbm", topes) -- 104 sec
  
  -- Now let's add in all the CY's total, including all triangulations.
  Qs = readCYPolytopes "../Databases/cys-ntfe-h11-3.dbm";
  elapsedTime for Q in values Qs do (
    addToCYDatabase("../Databases/cys-ntfe-h11-3.dbm", Q, NTFE => true);
    ) -- 10 sec


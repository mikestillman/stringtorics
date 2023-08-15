-- In this file, we create the h^(1,1)=3 database of all of the hypersurfaces 
-- in (simplicial resolutions of) Fano toric 4-folds which are smooth CY3-folds.

  restart
  needsPackage "StringTorics"
  createCYDatabase("../Databases/test-cys-ntfe-h11-3", 3, {0, 10})
  createCYDatabase("../Databases/test-cys-ntfe-h11-3", 3, {10,20})
  createCYDatabase("../Databases/test-cys-ntfe-h11-3", 3, {20,30})
  createCYDatabase("../Databases/test-cys-ntfe-h11-3", 3, {30,36})
  select(readDirectory("../Databases"), s -> match("range", s) and match("h11-3", s))
  (Xs, Qs) = readCYDatabase   "../Databases/test12-cys-ntfe-h11-3.dbm"

  -- TODO: ask Dan about running in different processes a number of these M2's.  What is the best way?
  --       ask Dan: maybe have special character sequence, e.g. @lo@ such that inside a string this is replaced by the string value of (in this case) 'lo'.
  makeDirectory("./Foo"|"3")
  lo = 0;
  hi = 0;
  total = 244;
  nper = 10;
  for i from 0 to 9 do (
      lo = hi;
      num = ceiling(total/nper);
      hi = lo+num;
      if hi > total then hi = total;
      exec := ///M2 -e 'needsPackage "StringTorics"' -e 'prefix = "./Foo3/test-h11-3"' -e 'range={///|lo|///,///|hi|///}' -e 'createCYDatabase(prefix,3,range)' -e 'exit 0' &///;
      print exec;
      run exec;
      )
  -- Separate function perhaps:
  createdFiles = sort for s in select(readDirectory("./Foo3"), s -> match("range", s) and match("h11-3", s)) list ("./Foo3/"|s)
  combineCYDatabases({"../Databases/test14-cys-ntfe-h11-3.dbm"} | createdFiles)


  topes = kreuzerSkarke(3, Limit => 1000);
  assert(#topes == 244)
  DBNAME = "../Databases/test-cys-ntfe-h11-3.dbm"
  elapsedTime addToCYDatabase(DBNAME, topes) -- 224 seconds

  restart
  needsPackage "StringTorics"
  topes = kreuzerSkarke(3, Limit => 1000);
  assert(#topes == 244)
  DBNAME = "../Databases/test-cys-ntfe-h11-3.dbm"
  elapsedTime addToCYDatabase(DBNAME, topes) -- 224 seconds

-- Query the results
  restart
  needsPackage "StringTorics"
  DBNAME = "../Databases/test-cys-ntfe-h11-3.dbm"
  R = ZZ[a,b,c]
  RQ = QQ (monoid R);
  (Qs, Xs) = readCYDatabase(DBNAME, Ring => R);
  

  restart
  needsPackage "StringTorics"
  topes = kreuzerSkarke(3, Limit => 1000);
  assert(#topes == 244)
  elapsedTime createCYDatabase("../Databases/cys-ntfe-h11-3.dbm", topes) -- 104 sec

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

-- A different test: here we want to compute the mori cone cap for each NTFE CY3.
-- and put this info into the database.
  restart
  debug needsPackage "StringTorics"
  topes = kreuzerSkarke(3, Limit => 1000);
  assert(#topes == 244)
  elapsedTime createCYDatabase("../Databases/test2-cys-ntfe-h11-3.dbm", topes) -- 104 sec

  Qs = readCYPolytopes "../Databases/test2-cys-ntfe-h11-3.dbm";
  elapsedTime for Q in values Qs do (
    addToCYDatabase("../Databases/test2-cys-ntfe-h11-3.dbm", Q, NTFE => true);
    ) -- 10 sec

  -- 1. Want to understand nefGenerators, toricMoriCone, in case Cl(V) is not torsion.
  -- 2. in case Cl(V) is torsion, what do we do then?
  -- 3. The map from Pic V --> Cl V is used in nefGenerators.  What do we really want here?
  
  Q = Qs#20
  Xs = findAllCYs Q -- only 1.
  X' = Xs#0
  max X'
  rays X'
  cyPolytope X'
  mC = toricMoriCone X'
  rays dualCone mC

  V = normalToricVariety(rays X, max X)
  nefGenerators V
  fromPicToCl V

  DBR = "~/Dropbox/Collaboration/Physics-Liam/Inequivalent CYs/h11_4Toric/"
  RZ = ZZ[x0,x1,x2,x3]
  X = cyFromCYToolsDB(DBR, (128, 1123, 0), RZ)
  rays X
  max X  
  methods annotatedFaces
code 0  
code 1
code 2
  cyPolytope X
  rays X
  Q = cyPolytope X
  methods cyPolytope
  options cyPolytope -- ID

  -- to create a CalabiYauInToric from CYToolsCY3
  -- 1. create CYPolytope
  -- 2. make a hash table translating CYToolsCY3 rays to CYPolytope rays.
  --    stash this hash table?
  -- 3. grab basis indices from CYTools DB, translate it to new indices.
  
  -- given Q, X
  Xtable = hashTable for i from 0 to #rays X - 1 list (rays X)#i => i+1
  Qtable = hashTable for i from 0 to #rays Q - 1 list (rays Q)#i => i -- indexed starting at 0
  XtoQ = hashTable for p in rays X list if Qtable#?p then Xtable#p => Qtable#p else continue
  for a in basisIndices X list XtoQ#a -- new bases...
  tri = max X
  -- triangulation from X to CalabiYauInToric
  sort for s in max X list (s1 := drop(s, 1); sort for a in s1 list XtoQ#a)

  Qs = readCYPolytopes "../Databases/cys-ntfe-h11-4.dbm";

DB4 = "../Databases/cys-ntfe-h11-4.dbm"
R = ZZ[a,b,c,d]
RQ = QQ (monoid R);
(Qs, Xs) = readCYDatabase(DB4, Ring => R);

max Xs#(1123,0)
max Xs#(1123,1)
rays Xs#(1123,0)
rays Xs#(1123,1)
rays Q

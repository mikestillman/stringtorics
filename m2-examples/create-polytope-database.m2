///
-- Create the h11=3 databases.
-- First we should make sure the 3 files do not exist?
  restart
  needsPackage "StringTorics"
  topes = kreuzerSkarke(3, Limit => 10000);
  assert(#topes == 244)
  -- create the polytopes
  createCYDatabase("polytopes-h11-3.dbm", topes)
  
  -- create the triangulations database (cy's)
  F = openDatabase "polytopes-h11-3.dbm"
  for i in (sort keys F)/value do (
      addToCYDatabase("cys-h11-3.dbm", F, i)
      )

  for i in (sort keys F)/value do (
      addToCYDatabase("cys-ntfe-h11-3.dbm", F, i, NTFE => true)
      )
  close F
///

///
  -- Let's look at these databases a bit more, see what they are like.
  ntfeDB = openDatabase "cys-ntfe-h11-3.dbm"
  cyDB = openDatabase "cys-h11-3.dbm"
  topesDB = openDatabase "polytopes-h11-3.dbm"
  # keys ntfeDB == 306
  # keys cyDB == 526
  -- TODO: the following should be a DB function of some sort?
  elapsedTime Qs = hashTable for k in sort keys topesDB list (value k) => cyPolytopeData(topesDB#k, ID => value k);
  NTFEs = for k in sort keys ntfeDB list cyData(ntfeDB#k, i -> Qs#i);
  R = ZZ[a,b,c]
  for X in NTFEs do (
      elapsedTime T := topologicalData(X, R); -- this fails for some X: TODO: topologicalData should work for all cyData...
      );
  close ntfeDB
  close cyDB
  close topesDB
///

///
  -- Create the h11=4 database.
  restart
  needsPackage "StringTorics"
  topes = kreuzerSkarke(4, Limit => 10000);
  assert(#topes == 1197)
  elapsedTime createCYDatabase("polytopes-h11-4.dbm", topes) -- 730 sec
///

///
  -- Create the h11=5 database.
  restart
  needsPackage "StringTorics"
  topes = kreuzerSkarke(5, Limit => 10000);
  assert(#topes == 4990)
  elapsedTime createCYDatabase("polytopes-h11-5.dbm", topes)
///

///
  -- Create the h11=6 database.
  restart
  needsPackage "StringTorics"
  topes = kreuzerSkarke(6, Limit => 100000);
  assert(#topes == 17101)
  elapsedTime createCYDatabase("polytopes-h11-6.dbm", topes)
///

///
  -- Create the h11=7 databases.  We create one for each (h11, h12) pair.
  -- F#"hodge", F#"count" are the non numeric keys in the data base.
  restart
  needsPackage "StringTorics"
  elapsedTime topes = kreuzerSkarke(7, Limit => 100000);

  -- there were constructed after H below, using:
  -- (keys H)/last//unique//sort
  h12s = {19, 23, 25, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 
      38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 
      52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 
      66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 
      80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 
      94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 
      106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
      117, 119, 120, 121, 122, 123, 124, 125, 127, 130, 131, 
      133, 135, 137, 139, 140, 141, 142, 143, 145, 147, 149, 
      151, 153, 154, 155, 157, 159, 160, 161, 163, 167, 169, 
      171, 173, 175, 178, 179, 181, 183, 185, 187, 192, 193, 
      196, 197, 199, 205, 211, 217, 219, 223, 227, 247, 271, 295}
  
  hodgeNumbers = method()
  hodgeNumbers KSEntry := (ks) -> (
      str := toString ks;
      ans := regex("H:([0-9]+),([0-9]+)", str);
      if #ans != 3 then error "expected 3 matches";
      (value substring(str, ans#1#0, ans#1#1),
      value substring(str, ans#2#0, ans#2#1))
  )

  topes/hodgeNumbers;
  assert(#topes == 50376)
  
  H = partition(hodgeNumbers, topes);
  -- if this doesn't finish (e.g. M2 runs out of hash codes).  Delete the dbm file being
  -- worked on, and start over, replacing the "51" 2 lines below with the new value.
  for h in sort keys H do (
      if h#1 <= 51 then continue;  -- these have been done already
      dbname := "polytopes-h11-"|h#0|"-h12-"|h#1|".dbm";
      << "doing " << h << " placing into " << dbname << endl;
      createCYDatabase(dbname, H#h);
      F := openDatabaseOut dbname;
      F#"hodge" = toString h;
      F#"count" = toString(#H#h);
      close F;
      )
  
///

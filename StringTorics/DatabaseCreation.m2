----------------------------------
-- Code for creating data bases --
----------------------------------

createPolytopeDatabase = method(
    Options => {
        "Hodge" => null, -- TODO: not used yet
        "Count" => null  -- TODO: not used yet
        })

createPolytopeDatabase(String, List) := opts -> (dbfilename, topes) -> (
    -- open data base file
    F := openDatabaseOut dbfilename;
    -- F#"info" = "4990 reflexive polytopes of h11=5"
    -- F#"topes" = toString topes;
    -- loop through topes, create CYPolytopeData, populate it, write it to data base.
    elapsedTime for i from 0 to #topes - 1 do elapsedTime (
        << "computing for polytope " << i << endl;
        V := cyPolytopeData(topes#i, ID => i); -- note that the polytope data is really that of the dual to topes#i.
        -- now fill it with data we want
        basisIndices V; -- compute them
        isFavorable V; -- compute h11, h21, favorability.
        annotatedFaces V; -- compute annotated faces
        -- now write it
        F#(toString i) = dump V;
        );
    --if opts#"Hodge" =!= null then F#"hodge" = toString(opts.Hodge);
    --if opts#"Count" =!= null then F#"count" = #topes;
    close F;
    )

addToCYDatabase = method(Options => {NTFE => false})

addToCYDatabase(String, CYPolytopeData) := opts -> (dbfilename, Q) -> (
    elapsedTime Xs := findAllCYs Q;
    << "  " << #Xs << " triangulations total" << endl;
    if opts.NTFE then (
        elapsedTime H := partition(restrictTriangulation, Xs);
        << "  " << #(keys H) << " NTFE triangulations" << endl;
        Xs = (keys H)/(k -> H#k#0); -- only take one triangulation that matches
        );
    F := openDatabaseOut dbfilename;
    for X in Xs do (
        computeIntersectionNumbers X; -- this should load all of the data we want
        F#(toString label X) = dump X;
        );
    close F;    
    )
addToCYDatabase(String, Database, ZZ) := opts -> (dbfilename, topesDB, i) -> (
    <<  "-- doing polytope " << i << endl;
    Q := cyPolytopeData(topesDB#(toString i), ID => i);
    addToCYDatabase(dbfilename, Q, opts);
    )

///
  -- Example of construction of database for: h11=3, all h12's.
  -- Note, #232 is not favorable.  But we can still construct it
  restart
  needsPackage "StringTorics"
  topes = kreuzerSkarke(3, Limit => 1000);
  assert(#topes == 244)
  elapsedTime createPolytopeDatabase("foo-h11-3.dbm", topes)

  -- Now let's add in all the CY's total, including all triangulations.
  -- Experiment: let's add them to the same database.
  F = openDatabase "foo-h11-3.dbm"
  Qlabels = sort select(keys F, k -> instance(value k, ZZ))
  Qs = for k in Qlabels list cyPolytopeData(F#k, ID => value k);
  close F  

  elapsedTime for Q in Qs do addToCYDatabase("foo-h11-3.dbm", Q);
///

///
  -- Example use of a constructed data base.
  restart
  debug needsPackage "StringTorics"
  RZ = ZZ[a,b,c]
  F = openDatabase "foo-h11-3.dbm"
    Xlabels = sort select(keys F, k -> (a := value k; instance(a, Sequence)))
    Qlabels = sort select(keys F, k -> (a := value k; instance(a, ZZ)))
    assert(#Xlabels == 526)
    assert(#Qlabels == 244)
    elapsedTime Qs = for k in Qlabels list cyPolytopeData(F#k, ID => value k);
    elapsedTime Xs = for k in Xlabels list cyData(F#k, i -> Qs#i, Ring => RZ);
    assert(Xs/label === Xlabels/value)
  close F

  Xs = hashTable for x in Xs list (label x) => x;  
  
  -- these are the ones we will consider.
  torsionfrees = for lab in sort keys Xs list (
      X := Xs#lab;
      V := normalToricVariety(rays X, max X); 
      if classGroup V === ZZ^3 then lab else continue
      )

  -- 6 polytopes are not torsion free (i.e. classGroup is not torsion free)
  nontorsionfrees = for lab in sort keys Xs list (
      X := Xs#lab;
      V := normalToricVariety(rays X, max X); 
      if classGroup V != ZZ^3 then lab else continue
      )

  nonFavorables = select(sort keys Xs, lab -> not isFavorable cyPolytopeData Xs#lab)
  favorables = select(sort keys Xs, lab -> isFavorable cyPolytopeData Xs#lab)

  elapsedTime H1 = partition(lab -> topologicalData Xs#lab, favorables); -- take one from each.
  netList (values H1)
  SAME = hashTable for k in keys H1 list (first H1#k) => drop(H1#k, 1)

  labelYs = sort keys SAME -- these have non-equal topological data 
  #labelYs == 291 -- all are equivalent to one of these, so there are 291 possible different topologies
    -- (but probably less than this).
    
  elapsedTime INV = partition(lab -> invariants Xs#lab, labelYs); -- 98 sec
  #INV == 274
  #INV == 273
  
  -- This selects the different topologies, except it can't tell about (52,0), (53,0).
  (keys INV)/(x -> #INV#x)//tally

  RESULT1 = hashTable for k in sort keys INV list (
      labs := INV#k;
      k => partitionByTopology(labs, Xs, 15)
      )
  
  tally for k in keys RESULT1 list #RESULT1#k -- all 1's, meaning that in each case, all elements of 
  -- INV#k are equivalent.
  
  (keys RESULT1)/(k -> (
          for x in keys RESULT1#k list x => join(RESULT1#k#x, SAME#x)
          ))//flatten//hashTable

  INV2 = hashTable for k in keys INV list if #INV#k === 1 then continue else k => INV#k

  RESULT1 = hashTable for k in sort keys INV2 list (
      labs := INV2#k;
      k => partitionByTopology(labs, Xs, 15)
      )

  F = cubicForm Xs#(52,0)
  G = cubicForm Xs#(53,0)
  see trim saturate(ideal F + ideal jacobian F)
  see trim saturate(ideal G + ideal jacobian G)
--XXX

///

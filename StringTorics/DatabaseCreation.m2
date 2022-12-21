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

---------------------------------------------------------------------
  -- remove below this line?  
  Xs = select(Xs, X -> first label X =!= 232); -- remove (for now) the non-favorable example.
  
  -- Let's do invariants of each h12 separately.
  possibleH12s = {43, 45, 51, 57, 59, 63, 65, 66, 67, 69, 71, 72, 73, 75, 76, 77, 78, 79, 81, 83, 84,
      85, 87, 89, 91, 93, 95, 99, 103, 105, 107, 111, 115, 119, 123, 127, 131, 141, 165,
      195, 231, 243}
  XsByH12 = partition(X -> hh^(1,2) X, Xs);
  sort keys XsByH12 == possibleH12s

  elapsedTime INV = partition(invariants, Xs);
  (sort keys INV)/(k -> #INV#k)  

  GV = for k in sort keys INV list k => (for X in INV#k list (X => partitionGVConeByGV(X, DegreeLimit => 15)))
  
  elapsedTime INV = hashTable for x in sort keys XsByH12 list x => partition(invariants, XsByH12#x);
  hashTable for h12 in keys INV list h12 => (
      netList for k in keys INV#h12 list (INV#h12#k/label)
      )
  
  -- let's try one set.
  netList INV#45#{9, 2, 1, 1, -1, 2, 3, 45}
  cys = INV#45#{9, 2, 1, 1, -1, 2, 3, 45}
  netList oo
  H = partition(topologicalData, cys)
  for k in keys H list k => H#k/label
  k = first select(keys H, k -> #H#k > 1)
  cys = H#k
  cys/label
  (keys H)/(k -> first H#k)

  for k in sort keys INV list for inv in sort keys INV#k list {inv, 
  XGVs = cys/(x -> (x => partitionGVConeByGV x))
  partitionByTopology XGVs
  gv1 = partitionGVConeByGV(cys_0, DegreeLimit => 15)
  gv2 = partitionGVConeByGV(cys_1, DegreeLimit => 15)
  findLinearMaps(gv1, gv2)
  
  Ps = select(Xs, x -> first label x == 27 or first label x == 38)
  Ps/topologicalData//netList
  
  torsionfrees = positions(Xs, X -> (
          V := normalToricVariety(rays X, max X); 
          classGroup V === ZZ^3
          ))

  nontorsionfrees = positions(Xs, X -> (
          V := normalToricVariety(rays X, max X); 
          classGroup V != ZZ^3
          ))
  nontorsionfrees = {0, 14, 15, 16, 17, 18, 103, 110} -- also 232, but that has been removed already
  nontorsionfrees/(t -> label Xs#t)
  elapsedTime INV = partition(invariants, Xs); 
  netList apply(sort keys INV, k -> k => netList (INV#k/label)) 
  (sort keys INV)/(k -> #INV#k)

  -- 14: all 3 triangulations have same data, and don't match any other polytope.
  these = select(Xs, x -> first label x == 14)
  these/topologicalData//unique

  these = select(Xs, x -> first label x == 15 or first label x == 12)
  netList these
  these/topologicalData//unique
  netList oo
  these/cubicForm//unique
  these/c2Form
  
  X1s = Xs_torsionfrees;
  #X1s
  elapsedTime INV = partition(invariants, X1s); 
  netList apply(sort keys INV, k -> k => netList (INV#k/label))
  GV = for k in sort keys INV list k => (for X in INV#k list (X => partitionGVConeByGV(X, DegreeLimit => 15)))

  select(keys INV, k -> #INV#k > 1)
  kINV = sort keys INV;
  
  -- i=26 gives 2 diff tops.
  i = 34
  Ys = INV#(kINV#i)
  netList Ys
  Ys/label
  GVYs = for y in Ys list y => partitionGVConeByGV(y, DegreeLimit => 15)
  partitionByTopology(GVYs)

  ans = {};  
  for i from 0 to #kINV-1 do (
--  for i from 157 to #kINV-1 do (
--  for i from 156 to 156 do (
      Ys := INV#(kINV#i);
      if #Ys === 1 then continue;
      << "doing " << i << endl;
      GVYs := for y in Ys list y => partitionGVConeByGV(y, DegreeLimit => 15);
      ans = append(ans, partitionByTopology(GVYs));
      << "  ans: " << ans#-1 << endl;
      )
  
  hashXs = hashTable for X in Xs list (label X) => X;
  
  RZZp = ZZ/5[gens ring cubicForm hashXs#(227,0)]
  F1 = sub(cubicForm hashXs#(227,0), RZZp)
  F2 = sub(cubicForm hashXs#(229,0), RZZp)
  L1 = sub(c2Form hashXs#(227,0), RZZp)
  L2 = sub(c2Form hashXs#(229,0), RZZp)
  
  use ring F1
  for a from 0 to 4 list for b from 0 to 4 list for c from 0 to 4 list sub(F1, {a => a1, b => b1, 
  tally for i in (0,0,0)..(4,4,4) list sub(F1, {a => i#0, b => i#0, c => i#0});
  tally for i in (0,0,0)..(4,4,4) list sub(F2, {a => i#0, b => i#0, c => i#0})
  
  -- Let's try a different approach: start with a hash table of the X's we want to consider: label => X.
  -- RESULT will be: label => list of {label, matrix}, or just label (if it is exactly equal)
  -- step 1: compute their topologicalData.
  -- step 2: place into clumps based on these.  Keep a table (label of one to keep) => list of the rest.
  -- step 3: for one in each clump, compute its invariants. Partition these labels.
  -- step 4: for each key in this partition, if its value has only one element, place it into the RESULT
  -- step 5: for each key whose value has higher length: run the distinctTopologies code on this.
  -- step 6: 
  
  -- key: invariants
  -- value: hash table of potentially different topologies with those invariants

  -- For now, only do this on the torsion free examples.

  X1s = Xs_torsionfrees;
  X1s = hashTable for X in X1s list (label X) => X;
  labelXs = (keys X1s)//sort;
  elapsedTime H1 = partition(lab -> topologicalData X1s#lab, labelXs); -- take one from each.
  #X1s == 517
  #(keys H1) == 286
  labelYs = sort for top in keys H1 list first H1#top -- these have non-equal topological data 
    -- (but often will still be equivalent topologically).

  -- Now we will partition these according to invariants
  elapsedTime INV = partition(lab -> invariants X1s#lab, labelYs); -- 8.4 seconds
  #(keys INV) == 165 -- this is the minimum number of distinct topologies.
  tally for k in keys INV list #INV#k
  RESULT = new MutableHashTable;
  for k in keys INV do if #INV#k === 1 then RESULT#k = hashTable {INV#k#0 => {}};
  #(keys RESULT)
  
  INV2 = hashTable for k in keys INV list if #INV#k === 1 then continue else k => INV#k
  
  -- now, for each key in INV2, we get gv invariants, 
  #(keys INV2) == 63 -- there are 63 sets of invariants with >= 2 elements that we need to match up
  kINV2 = sort keys INV2

  i = 31
  LGVs = elapsedTime for lab in INV2#(kINV2#i) list X1s#lab => partitionGVConeByGV(X1s#lab, DegreeLimit => 15)
  partitionByTopology(LGVs)
  
  RESULT1 = hashTable for k in sort keys INV2 list (
      LGVs = elapsedTime for lab in INV2#k list X1s#lab => partitionGVConeByGV(X1s#lab, DegreeLimit => 15);
      k => partitionByTopology(LGVs)
      )

  RESULT2 = hashTable for k in sort keys RESULT1 list if # (keys RESULT1#k) == 1 then continue else k => RESULT1#k

  topologicalData(X1s#(227,0))

  -- let's try determine oif (227,0), (229,0) are homeomorphic.
  Xa = X1s#(227,0)
  Xb = X1s#(229,0)
  La = c2Form Xa
  Lb = c2Form Xb
  Fa = cubicForm Xa
  Fb = cubicForm Xb

  -- These are actually different. 
  singb = saturate(ideal Fb + ideal jacobian Fb)
  singa = saturate(ideal Fa + ideal jacobian Fa)

  -- let's try to determine if (225,0), (226,0) are homeomorphic.
  Xa = X1s#(225,0)
  Xb = X1s#(226,0)
  La = c2Form Xa
  Lb = c2Form Xb
  Fa = cubicForm Xa
  Fb = cubicForm Xb

  -- These are actually different. 
  singb = saturate(ideal Fb + ideal jacobian Fb)
  singa = saturate(ideal Fa + ideal jacobian Fa)
  RQ = QQ[a,b,c]
  Ga = sub(Fa, RQ)
  Gb = sub(Fb, RQ)
  saturate(ideal Ga + ideal jacobian Ga)
  saturate(ideal Gb + ideal jacobian Gb)

  Rp = ZZ/5[a,b,c]
  
  pointCount = (F, p) -> (
      R := (ZZ/p) (monoid ring F);
      Fp := sub(F, R);
      # for x in (0,0,0)..(p-1,p-1,p-1) list if sub(Fp, matrix{{x}}) == 0 then x else continue
      )
  pointCount(Fa, 2)
  pointCount(Ga, 2)
  for p in {2,3,5,7,11} list pointCount(Fa, p)
  for p in {2,3,5,7,11} list pointCount(Fb, p)
  
  Xa = X1s#(26,0)
  Xb = X1s#(37,0)
  Fa = cubicForm Xa
  Fb = cubicForm Xb

  Xa = X1s#(53,0)
  Xb = X1s#(53,1)
  Fa = cubicForm Xa
  Fb = cubicForm Xb
  for p in {2,3,5,7,11} list pointCount(Fa, p)
  for p in {2,3,5,7,11} list pointCount(Fb, p)

  Xa = X1s#(56,0)
  Xb = X1s#(68,0)
  Fa = cubicForm Xa
  Fb = cubicForm Xb
  for p in {2,3,5,7,11} list pointCount(Fa, p)
  for p in {2,3,5,7,11} list pointCount(Fb, p)

  Xa = X1s#(63,1)
  Xb = X1s#(72,0)
  Fa = cubicForm Xa
  Fb = cubicForm Xb
  for p in {2,3,5,7,11} list pointCount(Fa, p)
  for p in {2,3,5,7,11} list pointCount(Fb, p)

-----------------
-- Next try -----

  elapsedTime H1 = partition(lab -> topologicalData Xs#lab, torsionfrees); -- take one from each.
  elapsedTime H1 = partition(lab -> topologicalData Xs#lab, sort for x in keys Xs list if first x == 232 then continue else x); -- take one from each.
  #H1 == 286 -- this is the maximum number of different topologies.
  labelYs = sort for top in keys H1 list first H1#top -- these have non-equal topological data 
    -- (but often will still be equivalent topologically).

  -- Now we will partition these according to invariants
  invariants1 Xs#(35,0)
  elapsedTime INV1 = partition(lab -> invariants1 Xs#lab, labelYs); -- 90 sec
  #INV1 == 267 -- so at least this number of different topologies

  elapsedTime INV = partition(lab -> invariants Xs#lab, labelYs); -- 98 sec
  #INV == 269

  tally for k in keys INV list #INV#k
  RESULT = new MutableHashTable;
  for k in keys INV do if #INV#k === 1 then RESULT#k = hashTable {INV#k#0 => {}};
  #(keys RESULT)

  INV2 = hashTable for k in keys INV list if #INV#k === 1 then continue else k => INV#k

  RESULT1 = hashTable for k in sort keys INV2 list (
      labs := INV2#k;
      k => partitionByTopology(labs, Xs, 15)
      )


  INV2 = hashTable for k in keys INV1 list if #INV1#k === 1 then continue else k => INV1#k

  RESULT1 = hashTable for k in sort keys INV2 list (
      labs := INV2#k;
      k => partitionByTopology(labs, Xs, 15)
      )

{(0, 0), (9, 0), (9, 1), (10, 0), (10, 1), (10, 2), (55, 0), (62, 0), (232, 0)}  
///

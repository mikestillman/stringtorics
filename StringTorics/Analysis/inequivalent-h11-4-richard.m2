  -- VERSION #3.  Do the same, but on these examples, conbined with Richard's database (includes flopped CY3s)
  -- How to handle these?
  -- How to compute the GV cone directly?
restart  
debug needsPackage "StringTorics"
-- Here is the code for h11=3.
    RZ = ZZ[a,b,c]
    load "../../m2-examples/richard-example/richard-db.m2"
    DIR = "~/Dropbox/Collaboration/Physics-Liam/Inequivalent CYs/h11_3NontoricNewFlopCode/"
    T2s = flatten flatten for h21 in findH21s DIR list 
      for pol in findPolys(DIR, h21) list
        for cy in findCYs(DIR, h21, pol) list 
          (h21, pol, cy) => toList join(getTopology(DIR, h21, pol, cy, RZ), {3, h21});
    assert(#T2s == 487)
    T2s/first/toList/isToric_DIR//tally -- 220 are not toric, 267 are toric. QUESTION: why 267??
    Ts = T2s;
    hashTs = hashTable Ts;
    keyTs = Ts/first//sort
    #keyTs == 487

    elapsedTime H = partition(lab -> invariantsAll toSequence (hashTs#lab), keyTs); -- 
    #keys H == 277
    (keys H)/(k -> #H#k)//tally

    T1s = for lab in allXs list (
        X := Xs#lab;
        lab => ({c2Form X, cubicForm X, hh^(1,1) X, hh^(1,2) X}, X)
        );

-- Here is the code to get Richard's DB for h11=4
    RZ = ZZ[x0,x1,x2,x3]
    load "../../m2-examples/richard-example/richard-db.m2"
    
    -- These are all of the h11=4 favorable toric examples
    -- modulo 2face-equivalence, and modulo automorphisms.
    DIR = "~/Dropbox/Collaboration/Physics-Liam/Inequivalent CYs/h11_4Toric/"
    T2s = flatten flatten for h21 in findH21s DIR list 
      for pol in findPolys(DIR, h21) list
        for cy in findCYs(DIR, h21, pol) list 
          (h21, pol, cy) => toList join(getTopology(DIR, h21, pol, cy, RZ), {3, h21});
    assert(#T2s == 1760)
    Ts = hashTable T2s
    keyTs = sort keys Ts

    partition(lab -> invariantsAll(Ts#lab);

  allTs = Ts
  --allXs = torsionfrees -- these are the ones we consider
  allT = topologySet(allXs, Xs);
  info allT -- 2014 possibly different topologies

    

restart  
debug needsPackage "StringTorics"
-- Here is the code for h11=3.
    RZ = ZZ[x0,x1,x2,x3]
    DIR = "~/Dropbox/Collaboration/Physics-Liam/Inequivalent CYs/h11_4Toric/"
    Xs = findAllCYToolsCY3s(DIR, RZ);
    # sort keys Xs == 1760

tally for lab in sort keys Xs list isFavorable Xs#lab -- all true
tally for lab in sort keys Xs list isToric Xs#lab -- all true

  allXs = sort keys Xs
  allT = topologySet(allXs, Xs);
  info allT -- 1760 possibly different topologies

  allT1 = combineIfSame(allT, X -> (c2Form X, cubicForm X))
  info allT1
  equivalences allT1

  elapsedTime allT2 = separateIfDifferent(allT, invariantsAll) -- 318 seconds
  info allT2    

  elapsedTime allT3a = combineByGV(allT2, DegreeLimit => 10); --  sec

  elapsedTime allT3 = combineByGV(allT2, DegreeLimit => 15); --  sec

  info allT3
  

-- Example to find all h11=3 id's which 
-- (a) have non-torsion class group
-- (b) are favorable

restart
needsPackage "StringTorics"
getKreuzerSkarke(3, Limit=>1000); -- 244 of these
polys = parseKS getKreuzerSkarke(3, Limit=>1000);

-- for each, we want the simplicial toric variety
-- such that the class group is 
polytopes = for p in polys list (
    convexHull matrixFromString p_1
    );
elapsedTime nonfavs = select(0..#polytopes-1, p -> not isFavorable polytopes#p)
  -- ouch: that line (first time) takes 61 seconds on my mac book pro...
  -- stashes info, still takes 9 seconds to determine information...
  -- only one is not favorable (id#232)

-- Only 5 cannot be handled currently (since they have torsion class group, and so M2 cannot create their Cox rings).
bad = elapsedTime select(#polytopes, i -> (
        p := polytopes#i;
        X := reflexiveToSimplicialToricVariety p;
        r := rank cl X;
        cl X != ZZ^r
        ))
-- {0, 9, 10, 55, 62} are the "bad" ones

restart
debug needsPackage "StringTorics"
kk = ZZ/32003
topes = kreuzerSkarke(4, Access => "wget", Limit => 10000); -- 1197
assert(#topes == 1197)

-*
elapsedTime Vs = topes / (P -> elapsedTime reflexiveToSimplicialToricVariety(convexHull matrix P, CoefficientRing => kk));
torsionFrees = positions(Vs, V -> classGroup V == ZZ^4);
nonTorsionFrees = positions(Vs, V -> classGroup V != ZZ^4);
*-
nonTorsionFrees = {0, 3, 4, 5, 12, 15, 796, 800, 803, 1059, 1060, 1064, 1065, 1134, 1135, 1151, 1153, 1155}
torsionFrees = sort toList(set(0..#topes-1) - set nonTorsionFrees);
assert(#torsionFrees == 1179) -- not 1197!! -- so 18 are torsion...

nonfavs = for i in torsionFrees list elapsedTime if isFavorable convexHull matrix topes_i then continue else (print i; i);
  nonfavs = {} -- there are no non-favorable torsion free polytopes.

-- Let's loop through all examples, write results to a file.
-- For mat of the file:
-*
 -- Polytope 4
 -- Triangulation 0
   h11, h12
   rays
   max simplices
   (cleaned) GLSM degrees
   cubic form (as a polynomial in x,y,z,w)
   linear form (as a polynomial in x,y,z,w)
 -- Triangulation 1
*-


-- Choose one tope, and find all cubic forms, and linear forms for that
-- for topes_10: get 2 triangulations, identical topology (w/o needing change of coordinates.
-- topes_11: 1 triangulation
-- 13: 2 triangulations, not obviously equivalent, but probably so.
kk = ZZ/32003
A = matrix topes_14
RZ = ZZ[x_0..x_3]
V = elapsedTime reflexiveToSimplicialToricVariety(convexHull A, CoefficientRing => kk)
(V, basisIndices) = elapsedTime reflexiveToSimplicialToricVarietyCleanDegrees(convexHull A, CoefficientRing => kk)
GLSM = matrix transpose degrees ring V

P2 = polar convexHull A
LP = select(latticePointList P2, lp -> dim(P2, minimalFace(P2, lp)) <= 2)
mLP = transpose matrix LP
FRSTs = findAllFRSTs mLP
Vs = for frst in FRSTs list normalToricVariety(frst_0, frst_1, WeilToClass => GLSM)
assert(all(Vs, v -> transpose matrix degrees ring v == GLSM))

tops = for v in Vs list (
    t := topologyOfCY3(v, basisIndices);
    {t_0, t_1, sub(t_3, vars RZ), sub(t_2, vars RZ)}
    )
netList tops
netList unique tops

f = tops_0_3
g = tops_1_3
f-g
gens gb ideal jacobian f
gens gb ideal jacobian g
tops_0_2 == tops_1_2
tops_0_3 == tops_1_3

top0 = topologyOfCY3(Vs_0, basisIndices)
top1 = topologyOfCY3(Vs_1, basisIndices)
top0 = {top0_0, top0_1, sub(top0_2, vars RZ), sub(top0_3, vars RZ)}
top1 = {top1_0, top1_1, sub(top1_2, vars RZ), sub(top1_3, vars RZ)}
top0 == top1

RQ = QQ (monoid RZ)
decompose ideal jacobian(fQ = sub(f, RQ))
decompose ideal jacobian(gQ = sub(g, RQ))
inverseSystem fQ
inverseSystem gQ

-- We compute all of the different topologies for h11=4 toric CY3 hypersurfaces
-- Caveats:
--  non-favorables:  These are handled now, except not for GV invariants.
--   I believe the code works for these, just never finds isomorphisms using GV invariants.
--  N examples have Cl V torsion.  Currently our methods do not handle these 
--    examples.  Actually, they do handle these.

restart
debug needsPackage "StringTorics" -- the debug is because some functions are not yet exported.
DBNAME = "../Databases/cys-ntfe-h11-5.dbm"

RZ = ZZ[a,b,c,d,e]
RQ = QQ (monoid RZ);
elapsedTime (Qs, Xs) = readCYDatabase(DBNAME, Ring => RZ);

-- Considering invariants (not coming from GV invariants):
  allXs = sort keys Xs;
  
  nonfavorables = for k in keys Qs list (
      if not isFavorable Qs#k then k else continue
      )
  
  nonfavorableXs = sort select(keys Xs, lab -> member(first lab, nonfavorables))
  #nonfavorables == 93
  #nonfavorableXs == 134

  favorableXs = sort for k in keys Xs list if isFavorable Qs#(first k) then k else continue;
  #favorableXs == 11713

  torsionQs = for k in sort keys Qs list elapsedTime (
      istor := prune coker matrix rays Qs#k;
      if not isFreeModule istor then k else continue
      )
  assert(torsionQs == {1,2})

  torsions = for k in sort keys Xs list (
      istor := prune coker matrix rays Xs#k;
      if not isFreeModule istor then k else continue
      )
  assert(torsions == {(1, 0), (2, 0)})
  -- h12s = sort keys tally for lab in nonfavorableXs list hh^(1,2) Xs#lab
  -- tally for lab in allXs list hh^(1,2) Xs#lab
  
  -- ourXs = select(allXs, lab -> member(hh^(1,2) Xs#lab, h12s));
  -- #ourXs

  -- allT = topologySet(allXs, Xs);
  allT = topologySet(favorableXs, Xs);
  info allT -- (fav) 1 bucket. (11713 possibly different topologies)

  elapsedTime allT = separateIfDifferent(allT, invariantsH11H12) -- 3 sec
  info allT -- (fav) 113 different buckets

  elapsedTime allT = separateIfDifferent(allT, hubschInvariants); -- 127 sec TODO: speed this up!
  info oo -- (fav) 406 different buckets

  elapsedTime allT = separateIfDifferent(allT, hessianInvariants) -- 157 sec
  info oo -- (fav) : 1937 buckets

  elapsedTime PC = pointCounter(RZ, "Primes" => {2,3,5,7,11,13}, Projective => true); -- 8 sec
  elapsedTime allT = separateIfDifferent(allT, pointCounts_PC); -- 3564 sec
  info allT -- (fav) 7685 different.

  elapsedTime allT = separateIfDifferent(allT, X -> elapsedTime polynomialContent det hessian cubicForm X) -- 90 sec
  info oo -- (fav) still 7685

  elapsedTime allT = separateIfDifferent(allT, X -> elapsedTime cubicConductorInvariants X) -- forgot to check total time.
  info oo -- (fav) 7762.

  elapsedTime allT = combineByGV(allT, DegreeLimit => 5); -- 912 sec
  info allT -- (fav) 7762 <= # <= 8687

  elapsedTime allT = combineByGV(allT, DegreeLimit => 10); -- 396  sec
  info allT -- (fav) 7762 <= # <= 8183

  elapsedTime allT = combineByGV(allT, DegreeLimit => 15); -- 2609  sec
  info allT -- (fav) 7762 <= # <= 8079


  elapsedTime allT = separateIfDifferent(allT, X -> elapsedTime cubicLinearConductorInvariants X) -- this is a good one!
  info allT -- (fav) 7814 <= # <= 8079

  -- doing this one now
  elapsedTime PC = pointCounter(RZ, "Primes" => {17,19,23,(2,2),(2,3),(3,2),(2,4)}, Projective => true); -- 75 sec
  elapsedTime allT = separateIfDifferent(allT, pointCounts_PC); -- 1970 sec
  info allT -- (fav) 7818 <= # <= 8079

  "inequiv-reps-h11-5" << representatives allT << close; -- only contains the possibly equivalent sets, and only one rep from each.

--------------------------------------------
-- Now we use this list to find equiv classes or not...
--------------------------------------------
restart
debug needsPackage "StringTorics" -- the debug is because some functions are not yet exported.
DBNAME = "../Databases/cys-ntfe-h11-5.dbm"

RZ = ZZ[a,b,c,d,e]
RQ = QQ (monoid RZ);
needs "../FindEquivalence.m2"
(A,phi) = genericLinearMap RQ

elapsedTime (Qs, Xs) = readCYDatabase(DBNAME, Ring => RZ);

-- Considering invariants (not coming from GV invariants):
  allXs = sort keys Xs;
  
  nonfavorables = for k in keys Qs list (
      if not isFavorable Qs#k then k else continue
      )
  
  nonfavorableXs = sort select(keys Xs, lab -> member(first lab, nonfavorables))
  #nonfavorables == 93
  #nonfavorableXs == 134

  favorableXs = sort for k in keys Xs list if isFavorable Qs#(first k) then k else continue;
  #favorableXs == 11713

  REPS = value get "inequiv-reps-h11-5"

  elapsedTime for x in REPS list (
      elapsedTime for x1 in x list invariantsAll Xs#x1
      );
  INVS = oo/unique/first;

  SETS = partition(x -> {
          INVS#x#"comps sing FQ", 
          drop(INVS#x#"hessian shape", 1)
          }, splice{0..#REPS-1})
  (keys SETS)/(k -> #SETS#k)
  (sort keys SETS)/(k -> {k#0, k#1, #SETS#k, SETS#k})//netList

  -- this is the set where equivalenceBySingularLocus works well.
  SETa = for s in sort keys SETS list (
      t := tally first s;
      if t#?{4,1} and t#{4,1} >= 2 then s else continue
      )
  SET1 = sort toList(set keys SETS - set SETa)

  SETb = for s in SET1 list (
      t := tally ((last s)/first); -- tally the degrees of the generators.
      if t#?1 and t#1 >= 3 then s else continue
      )

  SET2 = sort toList(set SET1 - set SETb)  

  SETc = for s in SET2 list (
      t := tally (first s);
      if t#?{4,2} then s else continue
      )
  
  SET3 = sort toList(set SET2 - set SETc)

  assert(#SETa === 8)
  assert(#SETb === 4)
  assert(#SETc == 4)

  (sort SET2)/(k -> {k#0, k#1, #SETS#k, SETS#k})//netList
  (sort SET3)/(k -> {k#0, k#1, #SETS#k, SETS#k})//netList
  (sort SETb)/(k -> {k#0, k#1, #SETS#k, SETS#k})//netList
  
  SET2/(k -> (SETS#(k)))

  ----------------------------------------------------
  -- SETa --------------------------------------------
  ----------------------------------------------------
  ss = flatten for k in SETa list SETS#k
  ans = for s in ss list (
      << "-- performing set " << s << " " << REPS#s << endl;
      separateAndCombineByFcn(REPS#s, Xs, (A,phi), equivalenceBySingularLocus) -- 
      )
  ans/first//netList
  ntotal = (ss)/(k -> #REPS#k)//sum
  ndifferent = ans/first/length//sum
  nsame = ans/first/(x -> #flatten x - #x)//sum
  nbad = ans/last/length//sum
  -- total:  85 total
  --         81 different
  --          2 redundant
  --          2 bad.  But these are different...
  -- TODO: re-show this.  

  ----------------------------------------------------
  -- SETb --------------------------------------------
  ----------------------------------------------------
  -- this one takes a long time.  Speed it up!
  ss = flatten for k in SETb list SETS#k -- 42 different sets
  elapsedTime ans = for s in ss list (
      << "-- performing set " << s << " " << REPS#s << endl;
      separateAndCombineByFcn(REPS#s, Xs, (A,phi), equivalenceByHessian) --
      ) -- takes 4923 seconds...!
  ans/first//netList
  ntotal = (ss)/(k -> #REPS#k)//sum
  ndifferent = ans/first/length//sum
  nsame = ans/first/(x -> #flatten x - #x)//sum
  nbad = ans/last/length//sum
  -- numbers checked 12 Jan 2023.
  -- total:  96 total
  --         81 different
  --         15 redundant
  --          0 bad.
  -- here is the netList result
  print(ans/first//netList)
+---------------------------------+----------------------+-----------+-----------+
|{(2249, 0)}                      |{(2255, 0)}           |{(2270, 0)}|           |
+---------------------------------+----------------------+-----------+-----------+
|{(2594, 0)}                      |{(2663, 0)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(2863, 0), (2904, 0)}           |                      |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(2906, 2), (3034, 1)}           |                      |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(3229, 0)}                      |{(3251, 2)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(3231, 0)}                      |{(3286, 1)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(4021, 0)}                      |{(4059, 0)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(4021, 2)}                      |{(4059, 1)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(4221, 0)}                      |{(4224, 0)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(3275, 2)}                      |{(3302, 5)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(3688, 1)}                      |{(3746, 0)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(3708, 1)}                      |{(3745, 3)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(4326, 3)}                      |{(4359, 2)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(157, 1), (190, 1)}             |                      |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(693, 3), (736, 0)}             |                      |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(1248, 0)}                      |{(1248, 1)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(1834, 0)}                      |{(1877, 0)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(1837, 2)}                      |{(1882, 0)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(1837, 3)}                      |{(1837, 4)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(2256, 2)}                      |{(2343, 1)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(2612, 0)}                      |{(2674, 0)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(2624, 0), (2677, 0)}           |{(2682, 0), (2715, 0)}|{(2717, 1)}|           |
+---------------------------------+----------------------+-----------+-----------+
|{(2856, 0)}                      |{(2931, 0)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(3248, 0)}                      |{(3321, 3)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(3275, 1)}                      |{(3331, 5)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(3550, 1), (3612, 0)}           |                      |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(3587, 0), (3679, 0)}           |                      |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(3688, 0)}                      |{(3744, 2)}           |{(3746, 1)}|           |
+---------------------------------+----------------------+-----------+-----------+
|{(3701, 2)}                      |{(3702, 2)}           |{(3742, 1)}|{(3747, 2)}|
+---------------------------------+----------------------+-----------+-----------+
|{(4323, 0), (4358, 1)}           |                      |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(4326, 1)}                      |{(4359, 1)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(4343, 0), (4356, 0)}           |                      |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(4360, 1)}                      |{(4360, 3)}           |{(4387, 2)}|           |
+---------------------------------+----------------------+-----------+-----------+
|{(4487, 1)}                      |{(4533, 0)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(4493, 0)}                      |{(4537, 1)}           |{(4548, 3)}|           |
+---------------------------------+----------------------+-----------+-----------+
|{(4619, 1)}                      |{(4636, 1)}           |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(4623, 1), (4625, 2)}           |                      |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(4647, 0), (4652, 0)}           |                      |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(4714, 0), (4715, 0)}           |                      |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(4792, 0), (4794, 0), (4796, 0)}|                      |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(716, 1)}                       |{(807, 0)}            |           |           |
+---------------------------------+----------------------+-----------+-----------+
|{(1183, 1)}                      |{(1207, 0)}           |{(1354, 0)}|{(1370, 2)}|
+---------------------------------+----------------------+-----------+-----------+

  -- Now we are going to fix the FindEquivalence.m2 file to handle
  -- both Hessian factors and components of singular loci.
  -- for hessians: we grab all factors of the same type.
  -- then: in list1: list of all factors of the same type of the X1 side.
  --       in list2: list of all permutations of the corresponding factors on the X2 side.
  -- only then, when constructing the ideal, do we create the cartesian product of all 
  -- of these sets:
  -- list1 becomes the flattening of list1.
  -- list2 becomes the cartesian product of all elements of list2.


  
  separateAndCombineByFcn(REPS#87, Xs, (A,phi), equivalenceByHessian) -- not done
  separateAndCombineByFcn(REPS#132, Xs, (A,phi), equivalenceByHessian) -- good (length 5)
  set4s = positions(REPS, x -> #x == 4)
  set4s/(x -> invariantsAll(Xs#(first REPS#(x))))

  ----------------------------------------------------  
  ss = {{{4,1},{4,1}},{{5,1}}}
  ----------------------------------------------------  
  (SETS#ss)/(k -> #REPS#k)//sum
  ans = for s in SETS#ss list (
      << "-- performing set " << s << " " << REPS#s << endl;
      separateAndCombineByFcn(REPS#s, Xs, (A,phi), equivalenceBySingularLocus) -- 
      )
  assert all(ans, x -> last x === {})
  ans/first//netList
  ntotal = (SETS#ss)/(k -> #REPS#k)//sum
  ndifferent = ans/first/length//sum
  nsame = ans/first/(x -> #flatten x - #x)//sum
  nbad = ans/last/length//sum
  ---------------------------------------------------- 
  -- ss = {{{4,1},{4,1}},{{5,1}}}
  -- total:  40 total
  --          2 redundant
  --         38 different
  --          0 bad
  ----------------------------------------------------  
  
  ----------------------------------------------------    
  ss = {{{4,1},{4,1}},{{1,1},{4,1}}}
  ----------------------------------------------------  
  (SETS#ss)/(k -> #REPS#k)//sum
  ans = for s in SETS#ss list (
      << "-- performing set " << s << " " << REPS#s << endl;
      separateAndCombineByFcn(REPS#s, Xs, (A,phi), equivalenceBySingularLocus) -- 
      )
  assert all(ans, x -> last x === {})
  ans/first//netList
  ntotal = (SETS#ss)/(k -> #REPS#k)//sum
  ndifferent = ans/first/length//sum
  nsame = ans/first/(x -> #flatten x - #x)//sum
  nbad = ans/last/length//sum
  ----------------------------------------------------  
  -- ss = {{{4,1},{4,1}},{{4,1},{4,1}}}
  -- total:   6 total
  --          0 redundant
  --          6 different
  --          0 bad
  ----------------------------------------------------  
  
  ----------------------------------------------------
  ss = {{{4,1},{4,1}}, {{1, 1}, {1, 1}, {1, 1}, {2, 1}}}
  ---------------------------------------------------- 
  (SETS#ss)/(k -> #REPS#k)//sum
  ans = for s in SETS#ss list (
      << "-- performing set " << s << " " << REPS#s << endl;
      separateAndCombineByFcn(REPS#s, Xs, (A,phi), equivalenceBySingularLocus) -- 
      )
  ans/first//netList
  ntotal = (SETS#ss)/(k -> #REPS#k)//sum
  ndifferent = ans/first/length//sum
  nsame = ans/first/(x -> #flatten x - #x)//sum
  nbad = ans/last/length//sum
  -- total:   6 total
  --          5 different
  --          0 redundant
  --          1 bad
  -- bad one: {(4395, 0), (4397, 2)}
  equivalenceBySingularLocus((4395, 0), (4397, 2), Xs, (A,phi))
  (last oo)/see -- every single ideal has no solutino over ZZ: a rational number is forced in all cases.
  -- thus, these are different
  -- total:   6 total
  --          6 different
  --          0 redundant
  --          0 unknown
  ----------------------------------------------------
  
  ----------------------------------------------------
  select(sort keys SETS, x -> first x == {{4,1},{4,1},{4,1}})
  ss = flatten join for k in oo list SETS#k
  ----------------------------------------------------
  ans = for s in ss list (
      << "-- performing set " << s << " " << REPS#s << endl;
      separateAndCombineByFcn(REPS#s, Xs, (A,phi), equivalenceBySingularLocus) -- 
      )
  ans/first//netList
  ntotal = (ss)/(k -> #REPS#k)//sum
  ndifferent = ans/first/length//sum
  nsame = ans/first/(x -> #flatten x - #x)//sum
  nbad = ans/last/length//sum
  -- analyze the bad one
  flatten flatten flatten select(ans, x -> #x#1 > 0)
  equivalenceBySingularLocus((3406, 0), (3448, 3), Xs, (A,phi))
  (last oo)/see -- every single ideal has no solutino over ZZ: a rational number is forced in all cases.
  -- thus, these are different
  -- total:  29 total
  --         29 different
  --          0 redundant
  --          0 unknown
  ----------------------------------------------------
  
  ----------------------------------------------------  
  select(sort keys SETS, x -> first x == {{4,1},{4,1},{4,1},{4,1}})
  ss = flatten join for k in oo list SETS#k
  ----------------------------------------------------
  ans = for s in ss list (
      << "-- performing set " << s << " " << REPS#s << endl;
      separateAndCombineByFcn(REPS#s, Xs, (A,phi), equivalenceBySingularLocus) -- 
      )
  ans/first//netList
  ntotal = (ss)/(k -> #REPS#k)//sum
  ndifferent = ans/first/length//sum
  nsame = ans/first/(x -> #flatten x - #x)//sum
  nbad = ans/last/length//sum
  -- total:   4 total
  --          4 different
  --          0 redundant
  --          0 unknown
  ----------------------------------------------------



  ----------------------------------------------------
  -- SETc --------------------------------------------
  ----------------------------------------------------
  -- this one doesn't work yet.
  ss = flatten for k in SETc list SETS#k -- 42 different sets
  ans = for s in ss list (
      << "-- performing set " << s << " " << REPS#s << endl;
      separateAndCombineByFcn(REPS#s, Xs, (A,phi), equivalenceBySingularLocus) -- XXX doing this now.
      )
  ans/first//netList
  ntotal = (ss)/(k -> #REPS#k)//sum
  ndifferent = ans/first/length//sum
  nsame = ans/first/(x -> #flatten x - #x)//sum
  nbad = ans/last/length//sum


  (lab1, lab2) = ((2340, 3), (2386, 0))
  X1 = Xs#lab1
  X2 = Xs#lab2
  (L1,F1) = (c2Form X1, cubicForm X1)
  (L2,F2) = (c2Form X2, cubicForm X2)
  F1 = sub(F1, RQ)
  F2 = sub(F2, RQ)
  C1 = decompose(ideal F1 + ideal jacobian F1)
  C2 = decompose(ideal F2 + ideal jacobian F2)
  hessianMatches(F1, F2)
  hessianMatches(lab1, lab2, Xs, (A,phi))
  
-- So far: 181 CY3's have been considered.
--   # different ones: 164
--   # that are redundant: 17

netList SET2
------- at end of part I am doing now XXX ----
  (lab1, lab2) = toSequence REPS#150
  (lab1, lab2) = toSequence REPS#167
  (lab1, lab2) = toSequence REPS#176 -- hard, but likely doable: J has an element (var)*(linear) - 2 == 0.
  (lab1, lab2) = toSequence REPS#203

  -- these 4 are all the same.
  (lab1, lab2, lab3, lab4) = toSequence REPS#6
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi))
  equivalenceBySingularLocus(lab1, lab3, Xs, (A,phi))
  equivalenceBySingularLocus(lab1, lab4, Xs, (A,phi))

  -- these 2 are the same
  (lab1, lab2) = toSequence REPS#15
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- same

  --{50, 55, 65, 66, 70, 72, 73, 80, 85, 91, 92, 95, 96, 100, 103, 116, 117, 186, 193}
  -- all 3 of these are different
  (lab1, lab2, lab3) = toSequence REPS#50
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  equivalenceBySingularLocus(lab1, lab3, Xs, (A,phi)) -- diff
  equivalenceBySingularLocus(lab2, lab3, Xs, (A,phi)) -- diff
  -- all 3 of these are different
  (lab1, lab2, lab3) = toSequence REPS#55
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  equivalenceBySingularLocus(lab1, lab3, Xs, (A,phi)) -- diff
  equivalenceBySingularLocus(lab2, lab3, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#65
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#65
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#66
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#70
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#72
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- SAME
  (lab1, lab2) = toSequence REPS#73
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#80
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#85
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#91
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#92
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#95
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#96
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#100
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#103
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#116
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#117
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#186
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#193
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- SAME
  {50, 55, 65, 66, 70, 72, 73, 80, 85, 91, 92, 95, 96, 100, 103, 116, 117, 186, 193}/(x -> REPS#x)
  flatten oo -- 40 potential topologies, but 2 are counted twice, so this set gives 38. (38-19 new tops).

  {1, 2, 19, 20, 26, 64, 67, 97, 121, 123} 
  (lab1, lab2) = toSequence REPS#1
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#2
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- SAME
  (lab1, lab2) = toSequence REPS#19
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#20
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#26
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2, lab3) = toSequence REPS#64 -- all 3 diff
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  equivalenceBySingularLocus(lab1, lab3, Xs, (A,phi)) -- diff
  equivalenceBySingularLocus(lab2, lab3, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#67
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#97
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#121
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  (lab1, lab2) = toSequence REPS#123
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- diff
  {1, 2, 19, 20, 26, 64, 67, 97, 121, 123}/(x -> REPS#x)//flatten -- 21 - 1 = 20 topologies. (20 -10 = 10 new topologies)

  -- HARD ONE?
  (lab1, lab2, lab3) = toSequence REPS#87
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- ??
  equivalenceByHessian(lab1, lab2, Xs, (A,phi)) -- ??

  (lab1, lab2) = toSequence REPS#18
  equivalenceByHessian(lab1, lab2, Xs, (A,phi)) -- diff

  -- HARD ONE (and there are 37 groups in this set!
  (lab1, lab2, lab3, lab4) = toSequence REPS#0
  equivalenceByHessian(lab1, lab2, Xs, (A,phi)) -- diff

  (lab1, lab2) = toSequence REPS#135 -- diff
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- indeterminate!
  J = oo#1#0 -- also try oo#1#1, but it is the same poly here.
  f = J_(numgens J - 1)
  matrix for a from -6 to 6 list for b from -6 to 6 list sub(f, {(first support f) => a, (last support f) => b})
  -- f has no zeros other than (0,0).  Therefore these are NOT equivalent.
  equivalenceByHessian(lab1, lab2, Xs, (A,phi)) -- indeterminate!
  classifyExtremalCurves Xs#lab1
  classifyExtremalCurves Xs#lab2 -- really different, so these *should* be different

  (lab1, lab2) = toSequence REPS#164 -- different
  equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- indeterminate!
  Js = oo#1 -- also try oo#1#1, but it is the same poly here.
  f = Js#0 _ (numgens Js#0 - 1)
  assert(#support f == 2)
  assert(Js#1 _ (numgens Js#1 - 1) == f)
  matrix for a from -6 to 6 list for b from -6 to 6 list sub(f, {(first support f) => a, (last support f) => b})
  for a from -6 to 6 list for b from -6 to 6 list if sub(f, {(first support f) => a, (last support f) => b}) == 0 then (a,b) else continue
  sub(Js#0, {(support f)_0 => 0, (support f)_1 => 1}) -- no ZZ-solution
  sub(Js#0, {(support f)_0 => 1, (support f)_1 => 1}) -- no ZZ-solution

  -- Let's grab all the sets with sing locus containing (4,2), or two (4,1)'s.
  for k in sort keys SETS list if member({4,2}, k#0) or # (for a in k#0 list if a == {4,1} then a else continue) >= 2 then k else continue
  setsForSing = flatten for k in oo list for s in SETS#k list REPS#s
  for s in setsForSing list (
      lab1 = s#0;
      lab2 = s#1;
      (lab1, lab2) => equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)) -- ??
      )

  -- This group has # of labels.
  -- The number of sets:
  -- The number of new topologies on top of that number.
  -- Number of ones that are the same as another.
  hessSets = for k in sort keys SETS list if # (for a in k#1 list if a == {1,1} then a else continue) >= 2 then k else continue
  setsForHess = flatten for k in hessSets list for s in SETS#k list REPS#s
  hessAnswers = for s in setsForHess list (
      lab1 = s#0;
      lab2 = s#1;
      << "doing " << lab1 << " " << lab2 << endl;
      ans = equivalenceByHessian(lab1, lab2, Xs, (A,phi));
      << "  ans = " << ans << endl;
      (lab1, lab2) => ans
      )
  #hessAnswers == 71 -- so far we have just checked the first two for each.
  tally for x in hessAnswers list x#1#0

  -- now let's try all the ones that have larger size:
  hessAnswersLarger = for s in setsForHess list (
      if #s == 2 then continue;
      sall := subsets(s, 2);
      for s1 in sall list (
          lab1 = s1#0;
          lab2 = s1#1;
          << "doing " << lab1 << " " << lab2 << endl;
          ans = equivalenceByHessian(lab1, lab2, Xs, (A,phi));
          << "  ans = " << ans << endl;
          (lab1, lab2) => ans
          )
      )
  
  -- Now consider the rest of the sets
  SETS1 = sort toList(set keys SETS - set hessSets)

  -- Now let's focus on the singular locus (46 of these sets)
  -- This one is quite fast (4 are consistent, 42 are not).
  singSets = for k in sort SETS1 list if member({4,2}, k#0) or # (for a in k#0 list if a == {4,1} then a else continue) >= 2 then k else continue
  setsForSing = flatten for k in singSets list for s in SETS#k list REPS#s
  singAnswer = for s in setsForSing list (
      lab1 = s#0;
      lab2 = s#1;
      << "doing " << lab1 << " " << lab2 << endl;
      ans := ((lab1, lab2) => equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi)));
      << "  ans = " << ans << endl;
      ans
      )
  -- now let's try all the ones that have larger size:
  singAnswersLarger = for s in setsForSing list (
      if #s == 2 then continue;
      sall := subsets(s, 2);
      for s1 in sall list (
          lab1 = s1#0;
          lab2 = s1#1;
          << "doing " << lab1 << " " << lab2 << endl;
          ans = equivalenceBySingularLocus(lab1, lab2, Xs, (A,phi));
          << "  ans = " << ans << endl;
          (lab1, lab2) => ans
          )
      )

  SETS2 = sort toList(set SETS1 - set singSets)
  (sort SETS2)/(k -> {k#0, k#1, #SETS#k, SETS#k})//netList

    
  REPS/length//tally
  -- ones where the hessian completely factors
  -- 39 sets of these
  set1 = positions(INVS, a -> all(a#"hessian shape", x -> first x <= 1))
  set1/(k -> #REPS#k)
  for lab in REPS#(set1#0) list (cubicForm Xs#lab, c2Form Xs#lab)
  for s in set1 list (
      << "--- " << s << " ---- " << REPS#s << " ------- " << endl;
      ans := equivalenceByHessian(REPS#s#0, REPS#s#1, Xs, (A,phi));
      print ans;
      s => ans
      )

  -- 47 sets of these
  set2 = positions(INVS, a -> (
          nlins := # select(a#"hessian shape", x -> first x <= 1);
          nlins >= 3 and nlins <= 4
          ))
  set2/(k -> #REPS#k)//tally
  for s in set2 list (
      << "--- " << s << " ---- " << REPS#s << " ------- " << endl;
      ans := equivalenceByHessian(REPS#s#0, REPS#s#1, Xs, (A,phi));
      print ans;
      s => ans
      )

  tally for i in INVS list i#"comps sing FQ"
  
  -- 13 sets of these (the first pair in each group is INCONSISTENT).
  set3 = positions(INVS, i -> (i)#"comps sing FQ" === {{4, 1}, {4, 1}, {4, 1}})
  set3/(k -> #REPS#k)//tally -- 11:2, 1:3, 1:4.
  for s in set3 list (
      << "--- " << s << " ---- " << REPS#s << " ------- " << endl;
      ans := equivalenceBySingularLocus(REPS#s#0, REPS#s#1, Xs, (A,phi));
      print ans;
      s => ans
      )

  -- 25 sets of these
  set4 = positions(INVS, i -> (i)#"comps sing FQ" === {{4, 1}, {4, 1}})
  set4/(k -> #REPS#k)//tally -- 23:2, 2:3.
  for s in set4 list (
      << "--- " << s << " ---- " << REPS#s << " ------- " << endl;
      ans := equivalenceBySingularLocus(REPS#s#0, REPS#s#1, Xs, (A,phi));
      print ans;
      s => ans
      )
  -- of these, considering first pair in each: all but 2 are INCONSISTENT, 2 are CONSISTENT.

  -- 37 sets of these
  set5 = positions(INVS, i -> i#"comps sing FQ" === {{4, 1}})
  set5/(k -> #REPS#k)//tally -- 31:2, 6:3
  -- some of the following that are excluded we can do, but not all (yet).
  for s in set5 list if member(s, {118, 129, 137, 138, 144, 145, 159, 179, 181, 
          182, 183, 190, 191, 192, 194, 200}) then continue else (
      if s < 200 then continue;
      << "--- " << s << " ---- " << REPS#s << " ------- " << endl;
      ans := equivalenceBySingularLocus(REPS#s#0, REPS#s#1, Xs, (A,phi));
      print ans;
      s => ans
      )

  setBad1 = positions(INVS, i -> (i)#"comps sing FQ" == {{5,1}})
  tally for s in setBad1 list drop((INVS#s)#"hessian shape", 1)
  setBad1A = select(setBad1, i -> drop((INVS#i)#"hessian shape", 1) == {{5,1}})


  -- setBad1A
  {7, 9, 10, 14, 22, 23, 24, 27, 31, 36, 37, 49, 51, 56, 71, 93, 99, 101, 102, 108, 115, 130}
  INVS#7
  REPS#7
  X1 = Xs#(133,0)
  X2 = Xs#(165,0)
  (L1, F1) = (c2Form Xs#(133,0), cubicForm Xs#(133,0))  
  (L2, F2) = (c2Form Xs#(165,0), cubicForm Xs#(165,0))  
  (LQ1, FQ1) = (sub(L1, RQ), sub(F1, RQ))
  (LQ2, FQ2) = (sub(L2, RQ), sub(F2, RQ))
  for lab in {(133,0), (165,0)} list singularContentsQuartic Xs#lab

  A
  phi = map(target phi, target phi, A)

  classifyExtremalCurves X1
  classifyExtremalCurves X2

  -- This checks that A, phi are in sync...
  -- Notice, there is no 'transpose A' in the def of phi.
  I1 = ideal last coefficients(phi sub(LQ1, source phi) - sub(LQ2, source phi))
  I1 = sub(I1, ring A)
  I2 = ideal(A * sub(last coefficients LQ1, ring A) - sub(last coefficients LQ2, ring A))
  I1 == I2
  Ic = sub(ideal last coefficients(phi sub(FQ1, source phi) - sub(FQ2, source phi)), ring A)

  I3 = ideal(A * transpose matrix{{0,0,0,0,-1}} - transpose matrix{{0,0,1,-1,0}})
  I3b = ideal(A * transpose matrix{{0,0,0,0,-1}} - transpose matrix{{0,1,0,-1,0}})

  I4 = ideal(A * transpose matrix{{-1,0,-1,0,0}} - transpose matrix{{-1,0,1,0,2}})
  J = I2 + I3b + I4 + Ic
  J = trim J;
  A % J

  cf1 = last coefficients ((c2Form X1)//(polynomialContent c2Form X1))
  cf2 = last coefficients ((c2Form X2)//(polynomialContent c2Form X1))

  M1 = sub(cf1   | transpose matrix first for kv in pairs classifyExtremalCurves X1 list if kv#0#0 === "FLOP" then kv#1 else continue, QQ)
  possibles = transpose matrix first for kv in pairs classifyExtremalCurves X2 list if kv#0#0 === "FLOP" then kv#1 else continue
  M2s = for p in subsets(7, 4) list sub(cf2 | possibles_p, QQ)
  for m in M2s list m*M1^-1

  B = (ZZ/101) (monoid ring A)
  J = sub(Ic, B) + sub(I1, B)
  gbTrace=3;
  gbJ = ideal gens gb(J, DegreeLimit => 6);
  gbJ_* /(f -> size f)//tally
  deg5s = select(gbJ_*, f -> degree f <= {5})  ;
  (ideal deg5s)_*/size  

  
  Rp = (ZZ/32003) (monoid RZ)
  Lp1 = sub(L1, Rp)
  Lp2 = sub(L2, Rp)
  Fp1 = sub(F1, Rp)
  Fp2 = sub(F2, Rp)
  for a in {0,0,0,0,0}..{1,1,1,1,1} list if sub(F1, matrix{a}) == 0 then a else continue
  for a in {0,0,0,0,0}..{1,1,1,1,1} list if sub(F2, matrix{a}) == 0 then a else continue

  -- from Jakob: h12=41
  --84,-1,73.  84 comes from polytope 41
  --73,-1,84   73 compes from polytopes 36
  -- this doesn't seem right:
  invariantsAll Xs#(36,0)
  invariantsAll Xs#(36,1)
  invariantsAll Xs#(36,2)  
  invariantsAll Xs#(36,3)  
  invariantsAll Xs#(41,0)  
  T = topologySet({(36,0), (36,1), (36,2), (36,3), (41,0)}, Xs)
  separateIfDifferent(T, invariantsAll)
  equivalenceBySingularLocus(lab1, lab2, Xs, phi)  
  -- these are all trivially inequivalent.


  for lab in REPS#(set1#0) list (cubicForm Xs#lab, c2Form Xs#lab)

  mylist = {(38, 1), (43, 0)}
  for lab in mylist list singularContents Xs#lab
  for lab in mylist list singularContentsQuartic Xs#lab
  for lab in mylist list classifyExtremalCurves Xs#lab
  for lab in mylist list invariantsAll Xs#lab

  invariantsAll(Xs#(29,1))
    invariantsAll(Xs#(33,1))
    invariantsAll(Xs#(33,2))
    invariantsAll(Xs#(33,4))
    cubicForm Xs#(29,1)
    cubicForm Xs#(33,1)
    cubicForm Xs#(33,2)

  mylist = {(29,1), (33,1), (33,2), (33,4)}
  for lab in mylist list singularContents Xs#lab
  for lab in mylist list singularContentsQuartic Xs#lab
  for lab in mylist list classifyExtremalCurves Xs#lab

  factor det hessian cubicForm Xs#(mylist#0)
  factor det hessian cubicForm Xs#(mylist#1)
  factor det hessian cubicForm Xs#(mylist#2)
  factor det hessian cubicForm Xs#(mylist#3)
  c2Form Xs#(mylist#0)
  
  (A, phi) = genericLinearMap RQ
  J1 = ideal last coefficients phi (sub(c2Form Xs#(mylist#0), source phi) - sub(c2Form Xs#(mylist#1), source phi))
  J2 = ideal last coefficients phi (sub(a-2*b+c+d, source phi) - sub(c2Form Xs#(mylist#1), source phi))
  J = J1 + J2
  trim J
  J3 = ideal last coefficients phi (sub(cubicForm Xs#(mylist#0), source phi) - sub(cubicForm Xs#(mylist#1), source phi))
  J = trim(J1 + J2 + J3)
  decompose J
------------------------------------------------------------------------------------------
-- Let's do the above, but fix an hh^(1,2) that has a nonfavorable
tally for lab in nonfavorableXs list hh^(1,2) Xs#lab
-- 2 in h12=29.  These are *not* equivalent to any favorable.
-- AND: (3,0), (4,0) are DIFFERENT
h12 = 29
myXs = select(allXs, lab -> hh^(1,2) Xs#lab == h12)
  allT = topologySet(myXs, Xs);
  info allT -- 

  elapsedTime allT = separateIfDifferent(allT, invariantsH11H12) -- 
  info allT -- 

  elapsedTime allT = separateIfDifferent(allT, hubschInvariants); -- 127 sec TODO: speed this up!
  info oo -- 

  elapsedTime allT = separateIfDifferent(allT, hessianInvariants)
  info oo -- 

  PC = pointCounter(RZ, "Primes" => {2,3,5,7,11,13}, Projective => true);
  elapsedTime allT = separateIfDifferent(allT, pointCounts_PC);
  info allT 
  netList allT#"Sets" -- (3,0), (4,0) are different, and not equiv to a favorable.

h12 = 53
myXs = select(allXs, lab -> hh^(1,2) Xs#lab == h12)
  #myXs == 488
  allT = topologySet(myXs, Xs);
  info allT -- 

  elapsedTime allT = separateIfDifferent(allT, invariantsH11H12) -- 
  info allT -- 1 different buckets (fav+nonfav)

  elapsedTime allT = separateIfDifferent(allT, hubschInvariants); -- 127 sec TODO: speed this up!
  info oo -- after H11H12: 2 different buckets

  elapsedTime allT = separateIfDifferent(allT, hessianInvariants)
  info oo -- after hubsch+H11H12: 1956 buckets

  PC = pointCounter(RZ, "Primes" => {2,3,5,7,11,13}, Projective => true);
  elapsedTime allT = separateIfDifferent(allT, pointCounts_PC);
  info allT -- after H11H12+hubsch+hessian+this: >=349 distinct
  netList allT#"Sets" -- 
  select(nonfavorableXs, lab -> hh^(1,2) Xs#lab == 53) -- {(429, 0), (435, 0)} these are in their own class, not equiv to any favorable.
 
  factor det hessian cubicForm  Xs#(429,0) -- 2 linear forms.
  factor det hessian cubicForm  Xs#(435,0) -- 2 linear forms as well.

  -- The following shows that these 2 are isomorphic...  
  (A, phi) = genericLinearMap  RQ
  T = target phi -- and source
  U = coefficientRing T
  X1 = Xs#(429,0)
  X2 = Xs#(435,0)
  I = trim ideal last coefficients( phi matrix{{a+b-3*e, a+b-3*d, 12*a+12*b+32*c-6*d-6*e}} - matrix {{a+c-3*e, a+c-3*d, 12*a+32*b+12*c-6*d-6*e }})
  J = trim ideal last coefficients( phi sub(cubicForm X1, T) - sub(cubicForm X2, T))
  I = sub(I, U)
  I = I + sub(J, U)
  decompose I
  A0 = sub(A % sub(I, U), ZZ) -- is the isomorphism
  isEquivalent(X1, X2, transpose A0)

h12 = 65
myXs = select(allXs, lab -> hh^(1,2) Xs#lab == h12)
  {(1820, 0), (1822, 0), (1822, 1), (1822, 2), (1842, 0), (1842, 1), (1842, 2), (1845, 0)} -- nonfavs
  -- (1820, 0) (1845, 0) -- in one class with 4 favs
  -- (1822, 0) (1842, 1) -- in one bucket with 6 favs 
  -- (1822, 1), (1822, 2), (1842, 0), (1842, 2) in one bucket with a bunch of favs
  #myXs == 664
  allT = topologySet(myXs, Xs);
  info allT -- 

  elapsedTime allT = separateIfDifferent(allT, invariantsH11H12) -- 
  info allT -- 1 different buckets (fav+nonfav)

  elapsedTime allT = separateIfDifferent(allT, hubschInvariants); --
  info oo -- after H11H12: 2 different buckets

  elapsedTime allT = separateIfDifferent(allT, hessianInvariants)
  info oo -- after hubsch+H11H12: 1956 buckets

  PC = pointCounter(RZ, "Primes" => {2,3,5,7,11,13}, Projective => true);
  elapsedTime allT = separateIfDifferent(allT, pointCounts_PC);
  info allT -- after H11H12+hubsch+hessian+this: >= 379 distinct

  netList allT#"Sets" -- 
  select(nonfavorableXs, lab -> hh^(1,2) Xs#lab == h12) -- 
 
  factor det hessian cubicForm  Xs#(1820,0) -- 2 linears
  factor det hessian cubicForm  Xs#(1845,0) -- 2 linears
  factor det hessian cubicForm  Xs#(1856,0) -- 2 linears, fav

  factor det hessian cubicForm  Xs#(1822,0) -- factors completely
  factor det hessian cubicForm  Xs#(1842,1) -- same
  factor det hessian cubicForm  Xs#(1837,0) -- same, fav.

  factor det hessian cubicForm  Xs#(1822,1) -- factors completely
  factor det hessian cubicForm  Xs#(1822,2) -- same
  factor det hessian cubicForm  Xs#(1842,0) -- factors completely
  factor det hessian cubicForm  Xs#(1842,2) -- same
  factor det hessian cubicForm  Xs#(1837,1) -- same.  Favorable.  Probably equivalent.

h12 = 101
myXs = select(allXs, lab -> hh^(1,2) Xs#lab == h12)
  #myXs == 225
  myUnfavs = select(myXs, lab -> not isFavorable Xs#lab)
  {   (4425, 0), (4426, 0), (4429, 0), (4431, 0), 
      (4431, 1), (4437, 0), (4437, 1), (4439, 0), 
      (4441, 0), (4442, 0), (4445, 0), (4449, 0), 
      (4449, 1), (4461, 0), (4463, 0), (4463, 1), 
      (4464, 0), (4468, 0), (4468, 1)}

  for lab in myUnfavs list factor det hessian cubicForm Xs#lab

  allT = topologySet(myXs, Xs);
  info allT -- 

  elapsedTime allT = separateIfDifferent(allT, invariantsH11H12) -- 
  info allT -- 1 different buckets (fav+nonfav)

  elapsedTime allT = separateIfDifferent(allT, hubschInvariants); --
  info oo -- after H11H12+hubsch: 10 different buckets

  elapsedTime allT = separateIfDifferent(allT, hessianInvariants)
  info oo -- after hubsch+H11H12+hessian: 45 buckets

  PC = pointCounter(RZ, "Primes" => {2,3,5,7,11,13}, Projective => true);
  elapsedTime allT = separateIfDifferent(allT, pointCounts_PC);
  info allT -- after H11H12+hubsch+hessian+this: >= 98 distinct

  elapsedTime allT = separateIfDifferent(allT, singularContents)
  info oo -- after hubsch+H11H12+hessian+PC+sing: 98 buckets still

  elapsedTime allT = separateIfDifferent(allT, singularContentsQuartic)
  info oo -- after hubsch+H11H12+hessian:  99 buckets.  3726 seconds for one split...!

  -- let's grab the sets which have a nonfavorable
  for set1 in allT#"Sets" list if any(flatten set1, lab -> not isFavorable Xs#lab) then set1 else continue
  for set1 in allT#"Sets" list for lab 
----------------------------------
-- Try separating given h12 = 35--
----------------------------------
X35s = select(sort keys Xs, lab -> hh^(1,2) Xs#lab == 35)  
  allT = topologySet(X35s, Xs);
  info allT -- 1 bucket. 31 different possible.

  elapsedTime allT = separateIfDifferent(allT, hubschInvariants); -- 
  info oo -- 3 different buckets

  --PC = pointCounter(RZ, "Primes" => {2,3,5,7,11,13}, Projective => true);
  PC = pointCounter(RZ, "Primes" => {2,3,5,7,11,13,(2,2),(3,2),(2,3),(2,4),17}, Projective => true);
  elapsedTime allT1 = separateIfDifferent(allT, pointCounts_PC);
  info allT1 -- 17 different buckets
  netList representatives allT1
  allT = allT1

  elapsedTime allT = separateIfDifferent(allT, singularContents)
  info allT -- still 17 -- no change.

  elapsedTime allT = separateIfDifferent(allT, hessianInvariants)
  info allT -- still 17 -- no change.

  for S in allT#"Sets" list if #S == 1 then continue else (
      partition(lab -> hessianInvariants Xs#lab, S/first)
      )
netList oo

  
  
  netList representatives allT  
  toricMoriConeCap Xs#(9,0)
  toricMoriConeCap Xs#(19,2)
  heft  Xs#(9,0)
  heft  Xs#(19,2)

  flatten flatten select(allT#"Sets", s -> #s == 8)
  set8 = {(8, 0), (11, 0), (12, 0), (12, 1), (12, 2), (21, 0), (21, 1), (21, 2)}
  set8/(lab -> isFavorable Xs#lab)
  set8/(lab -> isFavorable polar Qs#(first lab))
  gvs = hashTable(set8/(lab -> lab => elapsedTime gvInvariants(Xs#lab, DegreeLimit => 16)));
  for lab in set8 list (
      X := Xs#lab;
      degvec := heft X;
      print degvec;
      for c in toricMoriConeCap X list classifyExtremalCurve(gvs#lab, c, 16, degvec)
      )

  transpose matrix toricMoriConeCap Xs#(8,0)
  transpose matrix toricMoriConeCap Xs#(11,0)
  f = (lab) -> (X := Xs#lab; degvec := heft X; for c in toricMoriConeCap X list classifyExtremalCurve(

  C120 = transpose matrix toricMoriConeCap Xs#(12,0)
  C121 = transpose matrix toricMoriConeCap Xs#(12,1)
  A = C120^-1 * C121_{1,0,2,3,4}
  A = C120^-1 * C121_{1,0,3,2,4}
  A = C120^-1 * C121_{0,1,3,2,4}
  phi = map(RZ, RZ, A)
  phi c2Form Xs#(12,0)
  c2Form Xs#(12,1)
  partitionByTopology(set8, Xs, 16)
  
  set8/(lab -> partitionGVConeByGV(Xs#lab, DegreeLimit => 5))
  
-- Here is the group of 66 that might still be equivalent after point counts ---
-- And the second largest group: of 30.
L = value get "allT-h11-5-after-pointCounts";
netList select(L, x -> #x > 9)
for x in L list {
    M := flatten x;
    M = select(M, lab -> isFavorable Xs#lab);
    allLi = topologySet(M, Xs);
    combineByGV(allLi, DegreeLimit => 10);
    print info allLi
    }

set66 = {(1817, 0), (1818, 0), (1821, 0), (1821, 1), (1821, 2), (1821, 3), (1821, 4), (1821, 5), (1821, 6), (1821, 7), (1821, 8), (1821, 9), (1829, 0), (1829, 1), (1829, 2), (1829, 3), (1829, 4), (1838, 0), (1838, 1), (1838, 2), (1838, 3), (1838, 4), (1838, 5), (1838, 6), (1838, 7), (1838, 8), (1838, 9), (1838, 10), (1838, 11), (1838, 12), (1838, 13), (1838, 14), (1839, 0), (1839, 1), (1839, 2), (1839, 3), (1839, 4), (1839, 5), (1839, 6), (1839, 7), (1839, 8), (1839, 9), (1843, 0), (1843, 1), (1843, 2), (1843, 3), (1843, 4), (1843, 5), (1843, 6), (1843, 7), (1843, 8), (1843, 9), (1849, 0), (1849, 1), (1849, 2), (1858, 0), (1858, 1), (1858, 2), (1858, 3), (1923, 0), (1923, 1), (1923, 2), (1923, 3), (1923, 4), (1923, 5), (2009, 0)}  
set30 = {(4425, 0), (4426, 0), (4428, 0), (4429, 0), (4438, 0), (4439, 0), (4440, 0), (4442, 0), (4444, 0), (4445, 0), (4446, 0), (4450, 0), (4454, 0), (4460, 0), (4461, 0), (4464, 0), (4465, 0), (4471, 0), (4476, 0), (4477, 0), (4478, 0), (4479, 0), (4481, 0), (4485, 0), (4488, 0), (4498, 0), (4499, 0), (4503, 0), (4504, 0), (4506, 0)}
partitionByTopology(set66, Xs, 5) -- these 66 are all the same.
partitionByTopology(set30, Xs, 5) -- 19 different classes left

-- 8 of the 30 are not favorable, 22 are.  The 22 are 
all30 = topologySet(set30, Xs);
  elapsedTime all30a = combineByGV(all30, DegreeLimit => 10); -- sec
  elapsedTime all30b = combineByGV(all30a, DegreeLimit => 15); -- sec
  --elapsedTime all30c = combineByGV(all30b, DegreeLimit => 18); -- didn't seem to work?!

  -- 9 different ones left.
  -- can we separate them by various invariants?
  separateIfDifferent(all30c, invariantsAll) -- didn't separate!
  -- let's do more point counts.
  
  elapsedTime PC = pointCounter(RZ, "Primes" => {(2,2),(3,2),(2,3), 17, 19, 23}, Projective => true);
  elapsedTime all30d = separateIfDifferent(all30b, pointCounts_PC);
  info all30d -- still 9 different classes!
  
  the9 = flatten representatives all30d
  toricMoriConeCap(Xs#(the9_0))

  set22 = select(set30, lab -> isFavorable Xs#lab)
  all22 = topologySet(set22, Xs)
  info all22
  elapsedTime all22a = combineByGV(all22, DegreeLimit => 5); -- sec
  info all22a -- <= 11 different.
  elapsedTime all22b = combineByGV(all22a, DegreeLimit => 10); -- sec
  info all22b  -- <= 2 different.
  elapsedTime all22c = combineByGV(all22b, DegreeLimit => 15); -- sec
  info all22c -- 1 DIFFERENT TOP HERE.

-- read in the different sets after point counts, from heaviside.
allS = value get "allT-h11-5-after-pointCounts";
allS = select(allS, x -> #x > 1)
tally (allS/(x -> #x)) -- 1662 pairs of 2 each.
-- Let's attack those.
allS2 = select(allS, x -> #x == 2)
gv90 = gvInvariants(Xs#(9,0), DegreeLimit => 16)
gv192 = gvInvariants(Xs#(19,2), DegreeLimit => 16)

degvec = heft Xs#(9,0)
for c in toricMoriConeCap Xs#(9,0) list dotProduct(c, degvec)
mori = posHull transpose matrix toricMoriConeCap Xs#(9,0)
matrix{hilbertBasis dualCone mori}
(transpose rays mori) * matrix{hilbertBasis dualCone mori}
classifyExtremalCurves(gv90, toricMoriConeCap Xs#(9,0), 16, heft Xs#(9,0))
classifyExtremalCurves(gv192, toricMoriConeCap Xs#(19,2), 16, heft Xs#(19,2))
toricMoriConeCap Xs#(19,2)

--------------- BELOW THIS IS PROBABLY FROM h11=4... ------------------------------
-- We collect the nontorsion, torsion, favorable, nonfavorable's.
   torsions = for k in keys Qs list (
      istor := prune coker matrix rays Qs#k;
      if not isFreeModule istor then k else continue
      )
   nonfavorables = for k in keys Qs list (
       if not isFavorable Qs#k then k else continue
      )
   favorablesXs = sort for k in keys Xs list (
       if not isFavorable Qs#(first k) then continue else k
      )
   #torsions == 6
   #nonfavorables == 12

   #(keys Xs) == 2014 -- this is the maximum number of inequivalent topologies

   nonfavorableQs = {796, 800, 803, 1059, 1060, 1064, 1065, 1134, 1135, 1151, 1153, 1155}
   nonfavorableXs = sort select(keys Xs, k -> member(first k, nonfavorableQs))
   nonfavorableXs = {
       (796, 0), -- potentially same as (798,1)? (which is equiv to 5 others) SAME YES the same!
       (796, 1), -- potentially same as (795,1)? (which is equiv to 5 others) -- seemingly not same
       (800, 0), -- potentially same as (844,0)? (which is equiv to 2 others) SAME 
       (803, 0), -- potentially same as (807,0)? (which is equiv to 5 others) SAME YES the same!
       (1059, 0), (1060, 0), (1064, 0), (1065, 0), -- these 4 have same c2, cubic (exactly same) 
         -- possibly same as (1071,0)? (which is itself equiv to 10 others) SAME

       (1134, 0), -- potentially same as (1136,0)? SAME
       (1135, 0), -- potentially same as (1139,0)? SAME
       (1151, 0), (1153, 0), -- potentially same each other, and as (1156,0)? (which is equiv to a bunch of others).  hessian doesn't factor...
       (1155, 0), -- UNIQUE this one is unique in its invariants bucket. bu hessian does factor.
       (1155, 1)} -- potentially same as (1168,0). hessian doesn't factor

  torsions = {0, 3, 4, 5, 12, 15}
  torsionCYs = sort select(keys Xs, k -> member(first k, torsions))
  assert(torsionCYs == {(0, 0), (3, 0), (4, 0), (5, 0), (12, 0), (15, 0)})

-- There are 6 polytopes whose toric variety has torsion class group.
-- Each of these has 1 NTFE triangulation.
-- These 6 form three topological classes.
-- (0,0) -- is alone, not equivalent to any other in Xs
-- (5,0), (4,0), (3,0) -- have equivalent c2, cubic forms,  not equivalent to any other in Xs
-- (15,0), (12,0) -- have equivalent c2, cubic forms, not equivalent to any other in Xs
-- These are distinct from all others at the invariant level.


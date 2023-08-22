-- I get:
-- 36 different polytopes
-- 39 different NTFE CY3's
-- 29 distinct topologies

restart
debug needsPackage "StringTorics" -- the debug is because some functions are not yet exported.
DBNAME = "../Databases/cys-ntfe-h11-2.dbm"

R = ZZ[a,b]
RQ = QQ (monoid R);
(Qs, Xs) = readCYDatabase(DBNAME, Ring => R);
sort keys Xs
tally for lab in keys Qs list #automorphisms Qs#lab

-- 39 different examples.

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
   assert(#torsions == 2)
   assert(#nonfavorables == 0)

---------------------------------------------------------------
-- Next step: How many of these 39 are distinct topologies? --
---------------------------------------------------------------
  allXs = sort keys Xs
  allT = topologySet(allXs, Xs);
  info allT -- 39 possibly different topologies, so far
  
  allT1 = combineIfSame(allT, X -> (c2Form X, cubicForm X))
  info allT1 -- down to 35 different

  equivalences allT1

  -- Question: are there any torsions or nonfavorables in here?
  -- Torsions: none on this list.
  -- Nonfavorables: (1059, 0), (1060, 0), (1064, 0), (1065, 0) are all the same (on the nose).
  
  elapsedTime allT2 = separateIfDifferent(allT1, invariantsAll) -- 1 second
  info allT2
  equivalences allT2
  representatives allT2
  
  -- Now we take the buckets consisting of more than one set, and try to make equivalences.
  elapsedTime allT3a = combineByGV(allT2, DegreeLimit => 20); -- 
  info allT3a  
  representatives allT3a
  
  -- only two pairs left to check:
  -- {(9, 0), (10, 0)}
  -- {(28, 0), (31, 0)}
  -- below shows that the data is equivalent over QQ for each pair, but not over ZZ.
  
  (A, phi) = genericLinearMap RQ  
  T = target phi

  Ls = {(9, 0), (10, 0)}
  (X1, X2) = toSequence (Ls/(lab -> Xs#lab))
  (L1, F1) = (c2Form X1, cubicForm X1)
  (L2, F2) = (c2Form X2, cubicForm X2)
  -- check X1,X2: SAME
  JL = ideal last coefficients(phi sub(L1, T) - sub(L2, T))
  JF = ideal last coefficients(phi sub(F1, T) - sub(F2, T))
  J = trim(JL + JF)
  trim J
  decompose oo -- not equivalent, no integral points.
  
  Ls = {(28, 0), (31, 0)}
  (X1, X2) = toSequence (Ls/(lab -> Xs#lab))
  (L1, F1) = (c2Form X1, cubicForm X1)
  (L2, F2) = (c2Form X2, cubicForm X2)
  factor det hessian F1
  factor det hessian F2
  -- check X1,X2: SAME
  Ja = ideal last coefficients(phi sub(b, T) % ideal sub(a, T))
  JL = ideal last coefficients(phi sub(L1, T) - sub(L2, T))
  JF = ideal last coefficients(phi sub(F1, T) - sub(F2, T))
  J = trim(Ja + JL + JF)
  decompose J -- no integral points.
  -- somewhat easier:
  JL = ideal last coefficients(phi sub(L1, T) - sub(L2, T))
  JF = ideal last coefficients(phi sub(F1, T) - sub(F2, T))
  J = trim(JL + JF)
  decompose J
    
  

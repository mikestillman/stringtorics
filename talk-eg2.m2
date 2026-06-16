-----------------------
-- A mystery surface --
-----------------------
restart
load "make-mystery-surface.m2"
kk = ZZ/32003
S = kk[a..e]
I = mysterySurface S

-- Let's start as before to get a general picture:
isHomogeneous I
codim I
degree I
isPrime I
decompose I

genus(I + ideal(a))
codim singularLocus I

R = S/I
X = Proj R
HH^0(OO_X)
-- So, X is a smooth connected surface
-- in P^4, of degree 6, and sectional genus 3

-- Where does X fit in the Kodaira-Enriques classfication?
-- Castelnuovo criterion: X rational 
--    iff h^1(OO_X) = h^0(2K) = 0
K = Ext^2(S^1/I,S^{-5}) -- Canonical sheaf
KX = sheaf K
HH^1(OO_X)
HH^0(KX)
HH^0(KX ** KX)

-- So X is a rational surface. Let's write it as the blowup of
-- a minimal rational surface.

-- Let's compute some intersection numbers:
intersectionNum = (M,N) -> euler ring M - euler M - euler N + euler(M ** N)

intersectionNum(K,K) 
  -- K.K = -1
H = (S^1/I) ** (S^{1})
intersectionNum(H, H)
  -- H.H == 6
intersectionNum(K, H)
  --  H.K = -2

-- this is (H+K).(H+K) = 6 - 2*2 -1 == 1  
intersectionNum(H ** K, H ** K)

HH^0((sheaf K)(1)) -- same as:
HH^0(sheaf(H ** K))

-- Our plan: use these 3 sections of H+K to map X ---> P^2.  By the general
-- theory of adjunction, this map is base point free, and nice.
-- The -1 lines on X get blown down, and otherwise it is birational

RS = kk[a,b,c,d,e,x,y,z]
B = intersect(ideal(a,b,c,d,e),ideal(x,y,z))
J = ideal(matrix{{x,y,z}} * sub(presentation K,RS))
Jsat = trim saturate(J,a*x);
betti Jsat
transpose gens Jsat
M = (gens Jsat)_{0..3};
(ms,cfs) = coefficients(M, Variables=>{a..e})
cfs
-- points in P2 (coordinates x,y,z) for which this matrix has rank < 4 
-- have: ideal in P^4 consists of <= 3 linear forms (together with the 
-- equations of X), so the inverse image in P^4 is potentially a line.
Eimage = minors(4,cfs) 
time C = primaryDecomposition Eimage;
netList C -- 10 points rational over kk!
Eimage == radical Eimage -- radical too

-- Let's look at the inverse image of one of these points.
F1 = trim saturate(Jsat + C_0,B)
netList F1_*
E1 = sub(F1,S)
E1 = trim E1
E1module = Hom(E1,S^1/I)
intersection(E1module, E1module)
intersection(E1module, K)

-- Now let's do all of them:
Es = apply(C, P -> trim sub(trim saturate(Jsat + P, B), S));
netList Es
Es/degree
Es/codim
Ms = Es/(Q -> Hom(Q,S^1/I));
matrix table(Ms,Ms,intersection)

-- What did we find out?
-- X is the blowup of P^2 at the 10 points of V(Eimage)
-- X has 10 exceptional lines on it.
-- X is called a Bordiga surface.



-- A cubic hypersurface in P^3 is
-- the blowup of 6 points in P^2.
-- How are these related?
restart
setRandomSeed 42
kk = ZZ/101
R = kk[x,y,z]
pts = random(ZZ^6, ZZ^3)
(vars R) || pts^{0}
trim minors(2, oo) -- ideal of that point
I = intersect for i from 0 to 5 list trim minors(2,(vars R) || pts^{i})
netList I_* -- 4 cubics
netList decompose I

-- Free resolutions are the bread and butter of computing in algebraic geometry
C = freeResolution I
dd^C

gens I
phi = syz gens I
I == minors(3, phi)

-- Now we create the ring map
-- w_i --> I_i
-- its kernel is the defining equation of the image
-- P2 --> P3
-- given by
-- (I_0, I_1, I_2, I_3)
-- it is itself defined by a cubic polynomial, in 4 variables.
S = kk[w_0..w_3]
IX = trim ker map(R, S, gens I)     
F = IX_0

-- How are F, phi related?
RS = kk[gens R, gens S]
wvars = matrix{{w_0..w_3}}
wvars * sub(phi, RS)
coefficients(oo, Variables => {x,y,z})
M = last oo
M = sub(M, S)
G = det M
G = 1/41 * G

F == G -- the same!


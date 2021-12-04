restart
options inverseSystem
S = ZZ/32003[x_0..x_3]
I = inverseSystem product gens S
lin = ()-> random(1,S)
lin()
f = d -> sum for i from 1 to d list (lin())^4
inverseSystem (f(7))
betti res oo
hashTable for n from 4 to 10 list n => betti res inverseSystem (f(n))
#oo
netList oo
g550 = x_0^4+x_1^4+x_2^4+x_3^4+(x_0+x_1+x_2+x_3)^4
g551 = x_0^4+x_1^4+x_2^4+x_3^4+(x_1+x_2+x_3)^4
g562 = x_0^4+x_1^4+x_2^4+x_3^4+(x_2+x_3)^4
g683 = x_0^4+x_1^4+x_2^4+x_3^4
betti res inverseSystem oo

g430 = x_0^4+x_1^4+x_2^4+x_3^4+(x_0+x_1+x_2+x_3)^4 + (x_0+x_1+x_2)^4
--g = x_0^4+x_1^4+x_2^4+x_3^4+(x_0+x_1+x_2+x_3)^4 + (x_0+x_1)^4
g441 = x_0^4+x_1^4+x_2^4+x_3^4+(x_1+x_2+x_3)^4 + (x_1+x_2)^4
--g = x_0^4+x_1^4+x_2^4+x_3^4+(x_1+x_2+x_3)^4 + (x_0+x_1)^4
g420 = x_0^4+x_1^4+x_2^4+x_3^4+(x_0+x_1+x_2+x_3)^4 + (x_0+x_1+2*x_2)^4
g400 = x_0*x_1*x_2*x_3 --- what is it in terms of powers of linear forms 
--g562 = x_0^4+x_1^4+x_2^4+x_3^4+(x_0+x_1)^4 + (x_0+2*x_1)^4
betti res inverseSystem g420


linforms={x_0,x_1,x_2,x_3,x_0+x_1+x_2+x_3,x_0+x_1+x_2}
linforms={x_0,x_1,x_2,x_3,x_0+x_1+x_2+x_3,x_0+x_1+2*x_2}
linforms={x_0,x_1,x_2,x_3,x_0+x_1+x_2,x_0+x_1}
linforms={x_0,x_1,x_2,x_3}
linforms={x_0,x_1,x_2,x_3,x_0+x_1+x_2+x_3}
linforms={x_0,x_1,x_2,x_3,x_0+x_1+x_2}
linforms={x_0,x_1,x_2,x_3,x_0+x_1+x_2+x_3,x_0+2*x_1+3*x_2+4*x_3}
genericTable = L -> (lins=L;
    f = 0;
    for i from 0 to (#lins-1) do f = f + random(ZZ/32003)*lins_i^4;
    betti res inverseSystem f) 
listTable = unique for j from 0 to 100 list genericTable(linforms)

genericFunction = (L, k) -> (lins = L;
    f=0;
    for i from 0 to (#lins-1) do f = f + random(ZZ/32003)*lins_i^k;
    return f)
inverseSystem genericFunction(linforms,4)
inverseSystem ideal(x_0^4,x_1^4,x_2^4,x_3^4)
inverseSystem ideal for elt in linforms list elt^4
linforms
A430 = inverseSystem g430
A683 = inverseSystem g683

listRLF=for i from 1 to 4 list genericFunction(linforms, 1)
f = 1
for l in listRLF do f = f*l
inverseSystem f
betti res oo

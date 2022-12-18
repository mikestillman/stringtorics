-- Let's make our own test for regular triangulations
-- rather: trying to understand the specifics of this
restart
needsPackage "StringTorics"

str = ///3 4  M:5 4 N:35 4 Pic:19 Cor:0
    1    0    0   -1
    0    1    0   -1
    0    0    1   -1
///    
A = matrixFromString last first parseKS str
P = convexHull A
LP = transpose matrix latticePointList P
A
A1 = matrix{{1,1,1,1,1}} || LP
det A1


Aw = matrix{{3,5,-4,1,0}} || LP
Q = convexHull Aw
dim Q
sgn = (x) -> if x > 0 then "+" else if x < 0 then "-" else "0"
for f in ((faces(1,Q))/first) list (
    (sgn det A1_f, sgn det(Aw_f))
    )
Aw

A
Aw
dim convexHull Aw
Aw
A1w = matrix {{1,1,1,1,1}} || Aw
R = QQ[x_1..x_3,y]
det(A1w_{0,1,2,3} | transpose matrix{{1,y,x_1,x_2,x_3}})
det(A1w_{0,1,2,4} | transpose matrix{{1,y,x_1,x_2,x_3}})
det oo


for i from 0 to numColumns A1w -1 list 
  det(A1w_{0,1,2,3} | A1w_{i})

for i from 0 to numColumns A1w -1 list 
  det(A1w_{0,1,2,4} | A1w_{i})

for i from 0 to numColumns A1w -1 list 
  det(A1w_{0,1,3,4} | A1w_{i})

for i from 0 to numColumns A1w -1 list 
  det(A1w_{1,2,3,4} | A1w_{i})

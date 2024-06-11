loadPackage("Tropical")
loadPackage("MatchingFields")
loadPackage("Polyhedra")
loadPackage("ToricDegenerations")

--Computing tuples and Plucker weightvectors of permutation matchingfields
--SIGMA
sigma = matchingFieldFromPermutation(3,6,{6,2,4,3,5,1})
tuplesSigma = getTuples(sigma)
--536
LSigma = grMatchingField(3,6,tuplesSigma)

weightSigma = matrix{getWeightPleucker(LSigma)}



--INTERMEDIATE, aka RHO in thesis
tuplesIntermediate = insert(2, {1,4,3},delete({1,3,4},tuplesSigma))
-- 143 instead of 134
LIntermediate = grMatchingField(3,6,tuplesIntermediate)

weightIntermediate = matrix{getWeightPleucker(LIntermediate)}

--TAU
tau = matchingFieldFromPermutation(3,6,{6,2,5,3,4,1})
tuplesTau = getTuples(tau)
--356
LTau = grMatchingField(3,6,tuplesTau)

weightTau = matrix{getWeightPleucker(LTau)}

-- Computing trop Gr(3,6)
R= QQ[p_(0,1,2), p_(0,1,3), p_(0,1,4), p_(0,1,5), p_(0,2,3), p_(0,2,4), p_(0,2,5), p_(0,3,4), p_(0,3,5), p_(0,4,5), p_(1,2,3), p_(1,2,4), p_(1,2,5), p_(1,3,4), p_(1,3,5), p_(1,4,5), p_(2,3,4), p_(2,3,5), p_(2,4,5), p_(3,4,5), Global=>false]


I = Grassmannian(2,5,R)
T = tropicalVariety(I)

conesT = maxCones(T)
raysT = rays(T)
linealityT = linealitySpace(T)

-- Computing cones containing the weights of the three matching fields

for i from 0 to #conesT-1 do(
    if contains(coneFromVData(submatrix(raysT,conesT#i), linealityT), transpose(weightTau) ) then( print(i); break)
)

for i from 0 to #conesT-1 do(
    if contains(coneFromVData(submatrix(raysT,conesT#i), linealityT), transpose(weightIntermediate) ) then( print(i); break)
)

for i from 0 to #conesT-1 do(
    if contains(coneFromVData(submatrix(raysT,conesT#i), linealityT), transpose(weightSigma) ) then( print(i); break)
)

-- weightSigma in cone 17
-- weightIntermediate in cone 138
-- weightTau in cone 1

--Rays of these cones

conesT#17 -- 12, 46, 53, 54
conesT#138 -- 12, 28, 46, 53
conesT#1 -- 46, 51, 53, 54

--Looking for adjacent non-prime cones
for i from 0 to #conesT -1 do(
    if isSubset( {46,53}, conesT#i) == true then print(i)
)
-- 1 (tau), 13, 17 (sigma), 106 (not prime), 116 (not prime), 124, 138, 142, 194, 206

conesT#13 -- 20, 46, 53, 54
conesT#194 -- 20, 41, 46, 53
conesT#206 -- 12, 41, 46, 53


conesT#116 -- 28, 46, 53, 62 -- NOT PRIME, adjacent to INTERMEDIATE
conesT#124 -- 20, 28, 46, 53
conesT#106 -- 46, 51, 53, 62 -- NOT PRIME, adjacent to 116 and SIGMA
conesT#142 -- 41, 46, 51, 53

-- Find all non-prime cones
for i from 0 to #conesT-1 do(
    weightcone = entries transpose interiorVector(coneFromVData(submatrix(raysT,conesT#i), linealityT));
    Rweight = QQ[p_(0,1,2), p_(0,1,3), p_(0,1,4), p_(0,1,5), p_(0,2,3), p_(0,2,4), p_(0,2,5), p_(0,3,4), p_(0,3,5), p_(0,4,5), p_(1,2,3), p_(1,2,4), p_(1,2,5), p_(1,3,4), p_(1,3,5), p_(1,4,5), p_(2,3,4), p_(2,3,5), p_(2,4,5), p_(3,4,5), Weights=>-weightcone, Global=>false];
    initialideal = ideal leadTerm(1, Grassmannian(2,5,Rweight));
    if isPrime initialideal == false then print(i)
)
-- NON-PRIME CONES: 
nonprime = {105, 106, 107, 108, 109, 110, 111, 112, 113, 115, 116, 117, 118, 119, 120, 121, 122, 135, 136, 178, 179, 180, 181, 182, 183, 184, 185, 1017, 1018, 1019, 1020, 1021, 1022, 1023, 1024, 1025, 1026, 1027, 1028, 1029, 1030, 1031, 1032, 1033, 1034}
adjnonprime={}
for i from 0 to #nonprime -1 do(
    if isSubset( {46, 53}, conesT#(nonprime#i) ) == true then adjnonprime = join(adjnonprime, (nonprime#i, conesT#(nonprime#i)))
)
adjnonprime -- {(106, {46, 51, 53, 62}), (116, {28, 46, 53, 62})}

--These are the only two non-prime cones containing ray 46 and ray 53 (r1 and r7 in thesis)

dim intersection((coneFromVData(submatrix(raysT, conesT#116), linealityT)),(coneFromVData(submatrix(raysT, conesT#106), linealityT)))

--this is one less than the dimension of the maximal cones, hence they are adjacent

--Computing their initial ideals

weight116 = (entries transpose interiorVector(coneFromVData(submatrix(raysT, conesT#116), linealityT)))_0
R116 = QQ[p_(0,1,2), p_(0,1,3), p_(0,1,4), p_(0,1,5), p_(0,2,3), p_(0,2,4), p_(0,2,5), p_(0,3,4), p_(0,3,5), p_(0,4,5), p_(1,2,3), p_(1,2,4), p_(1,2,5), p_(1,3,4), p_(1,3,5), p_(1,4,5), p_(2,3,4), p_(2,3,5), p_(2,4,5), p_(3,4,5), Weights=>-weight116, Global=>false]
initial116 = ideal leadTerm(1, Grassmannian(2,5,R116))


weight106 = (entries transpose interiorVector(coneFromVData(submatrix(raysT, conesT#106), linealityT)))_0
R106 = QQ[p_(0,1,2), p_(0,1,3), p_(0,1,4), p_(0,1,5), p_(0,2,3), p_(0,2,4), p_(0,2,5), p_(0,3,4), p_(0,3,5), p_(0,4,5), p_(1,2,3), p_(1,2,4), p_(1,2,5), p_(1,3,4), p_(1,3,5), p_(1,4,5), p_(2,3,4), p_(2,3,5), p_(2,4,5), p_(3,4,5), Weights=>-weight106, Global=>false]
initial106 = ideal leadTerm(1, Grassmannian(2,5,R106))


-- Manually constructing the modified ideal by applying algorithm to cone 116
I = Grassmannian(2,5,R)
missingBinomial = (findMissingBinomials(initial116))_0
R' = QQ[p_(0,1,2), p_(0,1,3), p_(0,1,4), p_(0,1,5), p_(0,2,3), p_(0,2,4), p_(0,2,5), p_(0,3,4), p_(0,3,5), p_(0,4,5), p_(1,2,3), p_(1,2,4), p_(1,2,5), p_(1,3,4), p_(1,3,5), p_(1,4,5), p_(2,3,4), p_(2,3,5), p_(2,4,5), p_(3,4,5), Y_1, h, Global=>false]
modifiedI = sub(I,R') + ideal(h*Y_1 - p_(0,1,5)*p_(2,3,4)- p_(0,1,3)*p_(2,4,5))

--Computing trop(I'), takes +-45 mins
neT = tropicalVariety(modifiedI)

conesneT = maxCones(neT)
raysneT = rays(neT)
linealityneT = linealitySpace(neT)

-- projecting rays from new embedding down to plucker embedding
n  = #entries (raysneT)_{0}
i=0
projectedraysneT = {}
while  i<numColumns(raysneT) do(
projectedraysneT=append( projectedraysneT, submatrix( transpose((raysneT)_{i}),{0..n-3}));
    i =  i+1)
projectedraysneT

--Turning projected rays into matrix

projectedRaysneTM = transpose(projectedraysneT#0)
for j from 1 to #projectedraysneT -1 do(projectedRaysneTM = projectedRaysneTM|transpose projectedraysneT#j)
projectedRaysneTM

--projecting lineality space

i = 0 
projectedLneT = {}
while  i<numColumns(linealityneT) do(
    projectedLneT =append( projectedLneT, submatrix( transpose((linealityneT)_{i}),{0..n-3}));
    i = i+1)
projectedLneT

--turning projected linealityspace into matrix

projectedLinealityneT = transpose(projectedLneT#0)
for j from 1 to #projectedLneT -1 do(projectedLinealityneT = projectedLinealityneT|transpose projectedLneT#j)
projectedLinealityneT

-- Checking which cones in the new embedding are projected onto C116, C138


C116 = coneFromVData(submatrix(raysT, conesT#116), linealityT)
C138 = coneFromVData(submatrix(raysT, conesT#138), linealityT)
C106 = coneFromVData(submatrix(raysT, conesT#106), linealityT)
C1 = coneFromVData(submatrix(raysT, conesT#1), linealityT)
C17 = coneFromVData(submatrix(raysT, conesT#17), linealityT)

--DONT run this, it quits my macaulay stating too many mutable objects have been created. Instead check only for one or
--two cones. Trying to calculate them all, one by one, consecutively, also causes my macaulay to give up

for Crays in conesneT do(
    C = coneFromVData(submatrix(projectedRaysneTM,Crays), projectedLinealityneT);
    if (dim(intersection(C,C116)) == dim(C116)) then print(position(conesneT,i-> i == Crays), "intersects 116");
    if (dim(intersection(C,C138)) == dim(C138)) then print(position(conesneT,i-> i == Crays), "intersects 138");
    if (dim(intersection(C,C106)) == dim(C106)) then print(position(conesneT,i-> i == Crays), "intersects 106");
    if (dim(intersection(C,C1)) == dim(C1)) then print(position(conesneT,i-> i == Crays), "intersects 1");
    if (dim(intersection(C,C17)) == dim(C17)) then print(position(conesneT,i-> i == Crays), "intersects 17")
)

-- 116: 169, 393
-- 106: 169, 393
-- 138: 183
-- 17: 40
-- 1: 30

Cne169 = coneFromVData(submatrix(raysneT, conesneT#169), linealityneT)
Cne183 = coneFromVData(submatrix(raysneT, conesneT#183), linealityneT)
Cne393 = coneFromVData(submatrix(raysneT, conesneT#393), linealityneT)
Cne40 = coneFromVData(submatrix(raysneT, conesneT#40), linealityneT)
Cne30 = coneFromVData(submatrix(raysneT, conesneT#30), linealityneT)

--Checking adjacency of cones in the new embedding

for i in {Cne169, Cne393} do(
    for j in {Cne183} do(
	F = intersection(i,j);
	if (dim(F) == dim(j)-1) then print(i,j)))

dim(intersection(Cne40,Cne183)) == dim(Cne183) - 1
dim(intersection(Cne169,Cne183)) == dim(Cne183) - 1
dim(intersection(Cne169,Cne30)) == dim(Cne183) - 1

-- so the cones adjacent in T have cones projecting onto each that are adjacent in Tne

-- Computations of the newton okounkov bodies for wall crossing

-- Cone in plucker embedding: cones in new embedding projecting onto it
-- 116: 169, 393
-- 106: 169, 393
-- 138: 183
-- 17: 40
-- 1: 30

conesT#106
conesT#116
conesT#138
conesT#1
conesT#17


-- Computing the newton okounkov bodies for wall crossing

-- Taking the intersection of the cones 169 and 183 and taking the interior vectors to define matrix M
-- The idea for generating interior vectors is to take a fixed interior vector (vtemp) and add it
-- to every ray and to every generator of the lineality space

intersectionCone = intersection(Cne169,Cne183)
vtemp = interiorVector(intersectionCone)
rays(intersectionCone)
linealityIntersection = linealitySpace(intersectionCone)
dim(intersectionCone)

-- In the wall-crossing paper, the top row of the matrix M is taken to be the 1 vector
onevector = transpose( matrix({{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1}}))

linealityIntersection = onevector|submatrix'(linealityIntersection,,{0})

-- sanity check to make sure that this 1 vector is contained in the lineality space
-- and can be taken as generator instead of the first given generator
rank(linealityIntersection)
vsanitycheck = (linealitySpace(intersectionCone))_{0}
for i from 1 to 6 do(
    vsanitycheck = vsanitycheck+(linealitySpace(intersectionCone))_{i})
vsanitycheck

-- To obtain M we now put the generators of the linealityspace and the rays of the intersection
-- cone in a matrix as row vectors, and add vtemp to every row,
-- except the first one since we want it to be the 1-vector

zerovector = transpose( matrix({{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}}))
--indicator = transpose( matrix({{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1}}))
--L2 = linealityIntersection|indicator
--rank(L2)

Mtemp = zerovector
for i from 1 to 9 do(
    Mtemp = Mtemp|vtemp)

Mintersection =transpose((linealityIntersection|rays(intersectionCone))+Mtemp)
rank(Mintersection) 


--Matrix M is of full rank so the chosen interior vectors are indeed independent

--To obtain the matrices of the adjacent cones we add an interior vector as a row to the bottom of M

M169 = Mintersection||transpose(interiorVector(Cne169))
M183 = Mintersection||transpose(interiorVector(Cne183))

rank(M169)
rank(M183)

--Matrices M_1, M_2 in the paper, again of full rank

--Projecting to the original embedding would be forgetting the last two columns of these matrices
--Since these are the coordinates of the vectors corresponding to Y_1,h
--At first i thought the resulting matrices should still be of full rank, but 
-- the homogenization variable h prevents this

rank(submatrix(M169, {0..19}))
rank(submatrix(M183, {0..19}))
submatrix(M169, {0..19}) == submatrix(M183, {0..19})

--They are not of full rank, but they are different
--Checking if the maximal cones are of same dimension in both tropicalisations

dim(C116)
dim(Cne169)

--they are not, resulting in the above matrices never being able to be of full rank.
-- Again: caused by the homogeinisation
--The NO bodies are the convex hull of the columns of these matrices
-- (intersected with x_1 = 1, but our choice of the first row makes this a non-issue)

P169 = convexHull(M169)

dim P169
M169

P183 = convexHull(M183)

vertices(P183)
dim P183
M183

--The NO-bodies obtained using the same method on the projected cones are then

projP169 = convexHull(submatrix(M169, {0..19}))
vertices(projP169)
dim projP169

projP183 = convexHull(submatrix(M183, {0..19}))
vertices(projP183)
dim projP183

--We would like to compare them to the matching field polytopes
--Cone 169 corresponds to the nonprime cone and cone 183 to the intermediate cone

NOintermediate = NOBody(LIntermediate)

fVector(projP183) == fVector(NOintermediate)

--is true, to check they are actually isomorphic, verices were exported and polymake was used

--Constructing non-prime cone polytope from matching fields:
--importing non-lattice polytope vertices
nonprimeM = 1/2*transpose matrix{{2,0,0,0,0,0,0,2,0,0,0,0,0,0,2,0,0,0},
 {2,0,0,0,0,0,0,2,0,0,0,0,0,0,0,2,0,0},
 {2,0,0,0,0,0,0,0,0,2,0,0,0,0,2,0,0,0}, 
 {0,0,2,0,0,0,0,2,0,0,0,0,0,0,0,2,0,0},
 {2,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,2,0},
 {2,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,2,0},
 {0,0,2,0,0,0,0,2,0,0,0,0,0,0,0,0,2,0},
 {2,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,2,0},
 {0,0,0,2,0,0,0,2,0,0,0,0,0,0,0,0,2,0},
 {0,0,2,0,0,0,0,0,0,2,0,0,0,0,0,0,2,0},
 {2,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,2},
 {2,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,2},
 {0,0,2,0,0,0,0,2,0,0,0,0,0,0,0,0,0,2},
 {2,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,2},
 {0,0,0,2,0,0,0,2,0,0,0,0,0,0,0,0,0,2}, 
 {0,0,2,0,0,0,0,0,0,2,0,0,0,0,0,0,0,2},
 {2,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,2},
 {0,0,0,0,2,0,0,2,0,0,0,0,0,0,0,0,0,2}, 
 {0,0,2,0,0,0,0,0,0,0,2,0,0,0,0,0,0,2},
 {0,0,0,0,2,0,0,0,0,2,0,0,0,0,0,0,0,2},
 {1,0,0,0,1,0,0,1,1,0,0,0,0,0,0,1,0,1}}
Pnonprime = convexHull(nonprimeM)
dim Pnonprime

fVector(Pnonprime) == fVector(projP169)
fVector(Pnonprime)

--not true, the polytope from matching fields has one more vertex than the NO-body polytope

--If we remove the non integral vertex however:
nonprimeM' = 1/2*transpose matrix{{2,0,0,0,0,0,0,2,0,0,0,0,0,0,2,0,0,0},
 {2,0,0,0,0,0,0,2,0,0,0,0,0,0,0,2,0,0},
 {2,0,0,0,0,0,0,0,0,2,0,0,0,0,2,0,0,0}, 
 {0,0,2,0,0,0,0,2,0,0,0,0,0,0,0,2,0,0},
 {2,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,2,0},
 {2,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,2,0},
 {0,0,2,0,0,0,0,2,0,0,0,0,0,0,0,0,2,0},
 {2,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,2,0},
 {0,0,0,2,0,0,0,2,0,0,0,0,0,0,0,0,2,0},
 {0,0,2,0,0,0,0,0,0,2,0,0,0,0,0,0,2,0},
 {2,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,2},
 {2,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,2},
 {0,0,2,0,0,0,0,2,0,0,0,0,0,0,0,0,0,2},
 {2,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,2},
 {0,0,0,2,0,0,0,2,0,0,0,0,0,0,0,0,0,2}, 
 {0,0,2,0,0,0,0,0,0,2,0,0,0,0,0,0,0,2},
 {2,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,2},
 {0,0,0,0,2,0,0,2,0,0,0,0,0,0,0,0,0,2}, 
 {0,0,2,0,0,0,0,0,0,0,2,0,0,0,0,0,0,2},
 {0,0,0,0,2,0,0,0,0,2,0,0,0,0,0,0,0,2}}
Pnonprime' = convexHull(nonprimeM')
dim Pnonprime'

vertices P169

vertices P183

fVector(Pnonprime') == fVector(projP169)

--this is true, again polymake can be used to check they are isomorphic

--Eliminating the homogenization vertex by
--adding the last two columns of M169 and taking the polytope

testvertex = submatrix(M169, {20})+submatrix(M169, {21})
testM = submatrix(M169, {0..19})|(1/2*testvertex)
fVector(convexHull(testM))
fVector(convexHull(testM)) == fVector(Pnonprime)

--Again true, and polymake was used to check that they are isomorphic

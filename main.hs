import LinearAlgebra

test_vector = do
  putStrLn "\ntest_vector"
  let v1 = vector [1.0, 0.0, 0.0]
  print $ v1 == [1.0, 0.0, 0.0]

-- TODO vectorInt

test_matrix = do
  putStrLn "\ntest_matrix"
  let ihat = [1.0, 0.0, 0.0]
  let jhat = [0.0, 1.0, 0.0]
  let khat = [0.0, 0.0, 1.0]
  let matrix = [ihat, jhat, khat]
  print $ matrix == [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]

-- test_gjStep = do
--   putStrLn "\ntest_gjStep"
--   let dto = gjmatrix [
--       [1,0,0,0,0]
--     , [0,0,2,0,0]
--     , [0,0,0,1,0]
--     ]
--   let d1 = gjStep dto
--   print d1



main = do
  test_vector
  test_matrix

ihat = Just [1.0, 0.0, 0.0]
jhat = Just [0.0, 1.0, 0.0]
khat = Just [0.0, 0.0, 1.0]

nihat = negative <$> ihat

ih = vector [1.0, 0.0, 0.0]
jh = vector [0.0, 1.0, 0.0]
kh = vector [0.0, 0.0, 1.0]
m = matrix [ih, jh, kh]

v1 = vector [1.0, 2.0, 3.0]
v2 = vector [4.0, 5.0, 6.0]
v3 = vector [7.0, 8.0, 9.0]
v4 = vector [10.0, 11.0, 12.0]
m23 = matrix [v1,v2]
vtest = vector [10,20]

m33 = matrix [v1,v2,v3]
m43 = matrix [v1,v2,v3,v4]

v10 = vector [10.0, 20.0, 30.0]



matrix33 = [
    ['a','b','c']
  , ['d','e','f']
  , ['g','h','i']
  ]

m44 = [
    [1.0, 3.0, 1.0, 4.0]
  , [3.0, 9.0, 5.0, 15.0]
  , [0.0, 2.0, 1.0, 1.0]
  , [0.0, 4.0, 2.0, 3.0]
  ]

m22 = [
    [2,0]
  , [0,2]
  ]

r01 = [
    [0,0]
  , [2,0]
  ]

r02 = [
    [1,0,0]
  , [0,1,0]
  , [0,0,1]
  ]

r03 = [
    [1,0,0,0,0]
  , [0,0,1,0,0]
  , [0,0,0,1,0]
  ]

r04 = [
    [0,0,0,1,0]
  , [0,0,1,0,0]
  , [0,0,0,1,0]
  ]

r05 = [
    [2,6,-2]
  , [1,6,4]
  , [-1,4,9]
  ]

r06 = [
    [1,0,0,0,0]
  , [0,0,2,0,0]
  , [0,0,0,1,0]
  ]


-- main = do
--   putStrLn "Executing"
--   let v1 = Just $ vector [1.0,2.0,3.0]
--   let v2 = Just $ vector [4.0,5.0,6.0]
--   let dp1 = v1 |.| v2
--   putStrLn "Here's a working dot product"
--   print dp1
--   let ov1 = Just $ vector [0.0, 1.0]
--   let ov2 = Just $ vector [0.1, 0.0]
--   let dp2 = ov1 |.| ov2
--   putStrLn "Here's an orthogonal dot product"
--   print dp2
--   let addition1 = vectorAdd v1 v2
--   putStrLn "Vector Addition"
--   print addition1
--   let subtraction1 = vectorSub v2 v1
--   putStrLn "Vector Subtraction"
--   print subtraction1
--   let ihat = Just [1.0, 0.0, 0.0]
--   let jhat = Just [0.0, 1.0, 0.0]
--   let khat = Just [0.0, 0.0, 1.0]
--   let q24a = ihat |*| ihat
--   let q24b = ihat |*| jhat
--   let q24c = ((negative <$> ihat) |*| khat) |+| (jhat |*| ihat)
--   let q24d = (khat |*| jhat) |+| (ihat |*| ihat) |+| (jhat |*| khat) |+| (jhat |*| ihat)
--   putStrLn "Answers"
--   print q24a
--   print q24b
--   print q24c
--   print q24d
  -- let q24b = crossProduct ihat jhat
  -- let q24c = vectorAdd (crossProduct (negative ihat) khat) (crossProduct jhat ihat)
  -- let q24d = vectorAdd (crossProduct khat jhat) (vectorAdd (crossProduct ihat ihat) (vectorAdd (crossProduct jhat khat) (crossProduct jhat ihat)))
  -- putStrLn "Answers to 2.4"

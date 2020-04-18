module LinearAlgebra (
    vector
  , matrix
  , matrixScale
  , vectorJust
  , vectorAdd, (|+|)
  , vectorSub
  , negative
  , vectorScale
  , rows
  , cols
  , matrixVectorMult
  , dotProduct, (|.|), innerProduct
  , outerProduct
  , crossProduct, (|*|)
  , determinant, det
  , norm, magnitude
  , transpose
  , trace
  , display
  , isRREF
  -- , gjw
  , gjWalk
  , gjmatrix
  , gjStep
  -- , gjShowNumber
  , zeroesAtBottom
  , rowsLeadingOneToTheRight
  , findLeadingOnes
) where

import qualified Data.List as DL (transpose, intercalate)
import qualified Data.String as DL (unlines)

type Scalar = Double
type Component = Double
type Vector = [Component]
type Matrix = [Vector]

vector :: [Component] -> Vector
vector cs = cs

matrix :: [Vector] -> Matrix
matrix vs = vs

transpose :: Matrix -> Matrix
transpose matrix = DL.transpose matrix

isSquareMatrix :: Matrix -> Bool
isSquareMatrix matrix = length matrix == length (transpose matrix)

trace :: Matrix -> Maybe Scalar
trace matrix
  | not (isSquareMatrix matrix) = Nothing
  | otherwise = Just $ sum $ map (\n -> getnth n matrix) [0 .. ((length matrix) -1)]
  where getnth n matrix = matrix !! n !! n

matrixScale :: Matrix -> Scalar -> Matrix
matrixScale matrix scalar =
  map (vectorScale scalar) matrix

vectorScale :: Scalar -> Vector -> Vector
vectorScale scalar vector =
  map (* scalar) vector

rows :: Matrix -> Int
rows matrix = length $ matrix !! 0

cols :: Matrix -> Int
cols matrix = length matrix

matrixVectorMult :: Matrix -> Vector -> Maybe Vector
matrixVectorMult matrix vector
  | cols matrix /= length vector = Nothing
  | otherwise = Just $ map sum $ DL.transpose $ zipWith (vectorScale) vector matrix

vectorJust :: [Component] -> Maybe Vector
vectorJust xs = Just xs

vectorAdd :: Maybe Vector -> Maybe Vector -> Maybe Vector
vectorAdd (Just v1) (Just v2)
  | length v1 == length v2 = Just $ zipWith (+) v1 v2
  | otherwise = Nothing
(|+|) = vectorAdd

vectorSub :: Maybe Vector -> Maybe Vector -> Maybe Vector
vectorSub (Just v1) (Just v2)
  | length v1 == length v2 = Just $ zipWith (-) v1 v2
  | otherwise = Nothing

negative :: Vector -> Vector
negative vector = map (* (-1)) vector

dotProduct :: Maybe Vector -> Maybe Vector -> Maybe Scalar
dotProduct (Just v1) (Just v2)
  | length v1 == length v2 = Just $ sum $ zipWith (*) v1 v2
  | otherwise = Nothing
(|.|) = dotProduct
innerProduct = dotProduct

outerProduct :: Vector -> Vector -> Maybe Matrix
outerProduct v1 v2
  | length v1 /= length v2 = Nothing
  | otherwise = Just $ map (\c -> times v1 c) v2
  where times vector c = map (* c) vector

crossProduct :: Maybe Vector -> Maybe Vector -> Maybe Vector
crossProduct (Just (u1:u2:u3:[])) (Just (v1:v2:v3:[])) = Just [(u2*v3 - u3*v2), (u3*v1 - u1*v3), (u1*v2 - u2*v1)]
crossProduct _ _ = Nothing
(|*|) = crossProduct

determinant :: Matrix -> Scalar
determinant ((a:b:[]) : (c:d:[]) : []) = (a * d) - (b * c)
determinant matrix =
  altsum 1 $ headMinorProducts matrix
  where
    altsum _ [] = 0
    altsum 1 (x:xs) = x + (altsum 0 xs)
    altsum 0 (x:xs) = (-x) + (altsum 1 xs)
    index matrix = map (indexify 0) matrix
    indexify _ [] = []
    indexify n (x:xs) = (n,x) : indexify (n+1) xs
    deindex = map deindexify
    deindexify xs = map (\(a,b) -> b) xs
    getMinor matrix i = deindex $ map (filter (\(a,_) -> i /= a)) $ tail $ index matrix
    headMinorProducts matrix = zipWith (*) (head matrix) (map (determinant . getMinor matrix) [0..((length matrix) - 1)])
det = determinant

norm :: Vector -> Scalar
norm vector =
  sqrt $ sum $ map (**2) vector
magnitude = norm

unitVector :: Vector -> Vector
unitVector vector = map (/ (norm vector)) vector

cosTheta :: Maybe Vector -> Maybe Vector -> Maybe Scalar
cosTheta (Just v1) (Just v2) =
  case ((Just v1) |.| (Just v2)) of
    Nothing -> Nothing
    Just dp -> Just $ dp / ((magnitude v1) * (magnitude v2))

display :: Matrix -> String
display matrix =
  unlines $ map (DL.intercalate " ") $ (map . map) set $ transpose matrix

set :: Component -> String
set number =
  let
    shown = show number
    totalLength = 6
    numberLeft = totalLength - (length shown)
  in
    (replicate numberLeft ' ') ++ shown



-------------------------------
--- Gauss Jordan Elimination
-------------------------------

-- rows
-- cols

r06 = [
    [1,0,0,0,0]
  , [0,0,2,0,0]
  , [0,0,0,1,0]
  ]

dto = GJDTO {
    gjMatrix = gjmatrix r06
  , gjRowIndex = 0
  , gjColIndex = 0
  , gjPivotRows = []
  }

rowz = cols -- really dumb need to figure out inversions
colz = rows

data GJDTO = GJDTO {
    gjMatrix :: Matrix
  , gjRowIndex :: Int
  , gjColIndex :: Int
  , gjPivotRows :: [Int]
} deriving (Show)

gjWalk :: [GJDTO] -> [GJDTO]
gjWalk dtos =
  let
    dto = last dtos
    r = gjRowIndex dto
    c = gjColIndex dto
    n = (gjMatrix dto) !! r !! c
    ps = if (n == 1.0) then gjPivotRows dto ++ [r] else gjPivotRows dto
    newDTO = GJDTO {
        gjMatrix = gjMatrix dto
      , gjRowIndex = gjRowIndex $ gjStep dto
      , gjColIndex = gjColIndex $ gjStep dto
      , gjPivotRows = ps
    }
    res = if (n > 1) then [] else gjWalk [newDTO]
  in
    dto:res

gjStep :: GJDTO -> GJDTO
gjStep dto =
  let
    totalRows = rowz $ gjMatrix dto
    r         = gjRowIndex dto
    c         = gjColIndex dto
    pivots    = gjPivotRows dto
    pivotCieling []     = 0
    pivotCieling pivots = (maximum pivots) + 1
    getNewRow dto       = if (r < (totalRows-1)) then r + 1 else pivotCieling pivots
    getNewCol dto       = if (r < (totalRows-1)) then c else c + 1
    -- getNewPivots        = if -- get n then get the pivots
  in
    GJDTO {
        gjMatrix = gjMatrix dto
      , gjRowIndex = getNewRow dto
      , gjColIndex = getNewCol dto
      , gjPivotRows = gjPivotRows dto
    }

gjmatrix :: [[Integer]] -> Matrix
gjmatrix matrix =
  (map . map) fromIntegral matrix

-- gjw :: [[Integer]] -> (Int, Int, Component)
-- gjw matrix = gjWalk (gjmatrix matrix) (0,0) []



-- gjShowNumber :: Matrix -> (Int, Int) -> [Int] -> (Int, Int, Component, [Int])
-- gjShowNumber matrix (c, r) ps =
--   let
--     (fc, fr) = gjStep (rows matrix) (cols matrix) ps (c, r)
--     n = matrix !! fr !! fc
--   in
--     (fc, fr, n, [])

-- nums = [(c,r) | c <- [0..3],  r <- [0..2]]
-- showItBitch = map (\c -> (c, gjStep 2 3 c)) nums

type GJMatrix = [[Integer]]

gj01 = [
    [4,0,-1]
  , [2,-2,3]
  , [7,5,0]
  ]

gj02 = [[1],[2],[3],[4],[5],[6],[7]]

gjswap :: GJMatrix -> Int -> Int -> GJMatrix
gjswap matrix aindex bindex
  | aindex == bindex = matrix
  | aindex > (length matrix) - 1 || aindex < 0 = matrix
  | bindex > (length matrix) - 1 || bindex < 0 = matrix
  | otherwise = start ++ [matrix !! b] ++ middle ++ [matrix !! a] ++ end
    where
      a = min aindex bindex
      b = max aindex bindex
      start = take a matrix
      middle = drop (a + 1) $ take b matrix
      end = drop (b + 1) matrix

gjmult :: GJMatrix -> Int -> Integer -> GJMatrix
gjmult matrix rowIndex scalar =
  start ++ [new] ++ end
  where
    start = take rowIndex matrix
    end = drop (rowIndex + 1) matrix
    new = map (* scalar) $ matrix !! rowIndex

gjadd :: GJMatrix -> Int -> Int -> Integer -> GJMatrix
gjadd matrix originalIndex targetIndex scalar =
  start ++ [new] ++ end
  where
    start = take targetIndex matrix
    end = drop (targetIndex + 1) matrix
    addend = map (* scalar) $ matrix !! originalIndex
    new = zipWith (+) addend (matrix !! targetIndex)


-- http://mathonline.wikidot.com/reduced-row-echelon-form-of-a-matrix-rref

isRREF :: GJMatrix -> Bool
isRREF matrix =
  zeroesAtBottom matrix &&
  rowsHaveLeadingOnes matrix &&
  rowsLeadingOneToTheRight matrix &&
  colsHaveLeadingOnes matrix

-- REF

zeroesAtBottom :: GJMatrix -> Bool
zeroesAtBottom matrix =
  extractedZeroes == endOfMatrix
  where
    extractedZeroes = filter allZeroes matrix
    allZeroes row = row == replicate (length row) 0
    endOfMatrix = take (length extractedZeroes) $ reverse matrix

rowsLeadingOneToTheRight :: GJMatrix -> Bool
rowsLeadingOneToTheRight matrix = rrWalk (head leads) (tail leads)
  where
    leads = findLeadingOnes matrix
rrWalk _ [] = True
rrWalk (a,b) ((c,d):xs)
  | a < c && b < d = rrWalk (c,d) xs
  | otherwise = False

-- RREF

rowsHaveLeadingOnes :: GJMatrix -> Bool
-- I think this is implied in findLeadingOnes
rowsHaveLeadingOnes matrix = True

colsHaveLeadingOnes :: GJMatrix -> Bool
-- I think this is implied in findLeadingOnes
colsHaveLeadingOnes matrix = True

findLeadingOnes :: GJMatrix -> [(Int, Int)]
findLeadingOnes matrix = leadingOnesMatrix matrix 0
leadingOnesMatrix [] _ = []
leadingOnesMatrix (r:rs) n =
  leadingOnes r 0 n ++ leadingOnesMatrix rs (n+1)
leadingOnes [] _ _ = []
leadingOnes (x:xs) n r
  | x == 0 = leadingOnes xs (n+1) r
  | x == 1 = [(r,n)]
  | otherwise = []

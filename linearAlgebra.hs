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
  -- , gjmatrix
  -- , gjStep
  -- , gjShowNumber
  , zeroesAtBottom
  , rowsLeadingOneToTheRight
  , findLeadingOnes
) where

import qualified Data.List as DL (transpose, intercalate)
import qualified Data.String as DL (unlines)
import Data.Ratio

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

r = rMatrix [[1.0, 2.5], [1.5, 3.0]]



displayR :: RMatrix -> IO ()
displayR rMatrix =
  putStrLn $ unlines $ map (DL.intercalate " ") $ (map . map) setR $ rMatrix
  where
    setR num = (replicate (6 - (length $ raw num)) ' ') ++ (raw num)
    raw num =
      if denominator num == 1
        then show (numerator num)
        else show (numerator num) ++ "/" ++ show (denominator num)




-------------------------------
--- Gauss Jordan Elimination
-------------------------------

-- rows
-- cols

r00 = [
    [0,1]
  , [1,0]
  ]

r01 = [
    [1,1]
  , [0,0]
  ]

r02 = [
    [1,1,0]
  , [0,1,0]
  , [0,1,0]
  ]

r03 = [
    [1,1,0]
  , [0,2,0]
  , [0,1,0]
  ]

r04 = [
    [1,1,0]
  , [0,1,0]
  , [0,0,2]
  ]

r05 = [
    [1,0,0,0,0]
  , [0,0,2,0,0]
  , [0,0,0,1,0]
  ]

r06 = [
    [2,3]
  , [4,1]
  ]

type RMatrix = [[Rational]]

data RGJDTO = RGJDTO {
    rgjMatrix :: RMatrix
  , rgjRowIndex :: Int
  , rgjColIndex :: Int
  , rgjPivotRows :: [Int]
} deriving (Show)

rDto :: RMatrix -> RGJDTO
rDto rmatrix = RGJDTO {
    rgjMatrix = rmatrix
  , rgjRowIndex = 0
  , rgjColIndex = 0
  , rgjPivotRows = []
  }

rcolz :: RMatrix -> Int
rcolz matrix = length $ matrix !! 0

rrowz :: RMatrix -> Int
rrowz matrix = length matrix

rMatrix matrix = (map . map) toRational matrix

belowPivot dto = (rgjColIndex dto) `elem` (rgjPivotRows dto)
knockoutWithAdd dto = dto

rgjmatrix :: [[Integer]] -> RMatrix
rgjmatrix matrix =
  (map . map) toRational matrix

gjSolve :: RMatrix -> RMatrix
gjSolve rMatrix =
  let
    currentDTO = last $ gjWalkHist rMatrix
    matrix    = rgjMatrix currentDTO
    rowIndex  = rgjRowIndex currentDTO
    colIndex  = rgjColIndex currentDTO
    pivotRows = rgjPivotRows currentDTO
    scalar = 1.0 / (matrix !! colIndex !! rowIndex)
    ---
    original = firstNonZeroIndex matrix colIndex
    target = rowIndex
    firstNonZeroIndex matrix colIndex = (firstNonZeroRow matrix colIndex) - 1
    firstNonZeroRow [] _ = 0 -- this needs to be corrected
    firstNonZeroRow (r:rs) colIndex =
      if r !! colIndex == 0.0
        then 1 + (firstNonZeroRow rs colIndex)
        else 1

    addScalar = (- 4.0)
  in
    if (colIndex `elem` pivotRows)
      then gjadd matrix original target addScalar
      else gjmult matrix rowIndex scalar

r2 = RGJDTO {
  rgjMatrix = [[1 % 1,3 % 2],[0 % 1,(-5) % 1]],
  rgjRowIndex = 1,
  rgjColIndex = 0,
  rgjPivotRows = [0]
  }

gjWalkHist :: RMatrix -> [RGJDTO]
gjWalkHist rMatrix = gjWalkHist' $ gjMakeDto rMatrix

gjWalkHist' :: RGJDTO -> [RGJDTO]
gjWalkHist' dto =
  dto : if (gjEndDto dto)
    then []
    else gjWalkHist' $ gjNextDto dto

gjMakeDto :: RMatrix -> RGJDTO
gjMakeDto rMatrix = RGJDTO {
  rgjMatrix = rMatrix,
  rgjRowIndex = 0,
  rgjColIndex = 0,
  rgjPivotRows = []
}

gjNextDto :: RGJDTO -> RGJDTO
gjNextDto dto =
  let
    r       = rgjRowIndex dto
    c       = rgjColIndex dto
    n       = gjGetValue dto
    ps      = if (n == 1.0) then pivots ++ [r] else pivots
    tr      = rrowz $ rgjMatrix dto
    pivots  = rgjPivotRows dto
    pivotCieling []     = 0
    pivotCieling pivots = (maximum pivots) + 1
  in RGJDTO {
      rgjMatrix = rgjMatrix dto
    , rgjRowIndex = if (r < (tr-1)) then r + 1 else pivotCieling pivots
    , rgjColIndex = if (r < (tr-1)) then c else c + 1
    , rgjPivotRows = ps
  }

gjGetValue :: RGJDTO -> Rational
gjGetValue dto =
  (rgjMatrix dto) !! (rgjRowIndex dto) !! (rgjColIndex dto)

gjEndDto :: RGJDTO -> Bool
gjEndDto dto =
  let
    n = gjGetValue dto
    outOfBounds =
      rgjRowIndex dto >= (rrowz $ rgjMatrix dto) &&
      rgjColIndex dto >= (rcolz $ rgjMatrix dto)
    letsMakeItAPivot =
      n /= 1.0 && n /= 0.0
  in outOfBounds || letsMakeItAPivot

gjWalk :: [RGJDTO] -> [RGJDTO]
gjWalk dtos =
  let
    dto     = last dtos
    r       = rgjRowIndex dto
    c       = rgjColIndex dto
    n       = (rgjMatrix dto) !! r !! c
    ps      = if (n == 1.0) then pivots ++ [r] else pivots
    tr      = rrowz $ rgjMatrix dto
    tc      = rcolz $ rgjMatrix dto
    pivots  = rgjPivotRows dto
    pivotCieling []     = 0
    pivotCieling pivots = (maximum pivots) + 1
    newDTO  = RGJDTO {
        rgjMatrix = rgjMatrix dto
      , rgjRowIndex = if (r < (tr-1)) then r + 1 else pivotCieling pivots
      , rgjColIndex = if (r < (tr-1)) then c else c + 1
      , rgjPivotRows = ps
    }
    atEnd = r >= tr && c >= tc
    res = if (n /= 1 || atEnd) then [] else gjWalk [newDTO]
    -- res = gjWalk [newDTO]
  in
    dto:res

type GJMatrix = [[Integer]]

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

gjmult :: Num b => [[b]] -> Int -> b -> [[b]]
gjmult rmatrix rowIndex scalar =
  start ++ [new] ++ end
  where
    start = take rowIndex rmatrix
    new = map (* scalar) $ rmatrix !! rowIndex
    end = drop (rowIndex + 1) rmatrix

--gjadd :: GJMatrix -> Int -> Int -> Integer -> GJMatrix
gjadd :: Num b => [[b]] -> Int -> Int -> b -> [[b]]
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


-------------------

-- Ratio shitposting

type Thingy = Ratio
newRational = toRational 3.14
nr1 = toRational 3

type GJMatrix' = [[Rational]]
rr05 = [
    [1,0,0,0,0]
  , [0,0,2,0,0]
  , [0,0,0,1,0]
  ]

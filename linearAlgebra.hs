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
determinant ((a:b:c:[]):(d:e:f:[]):(g:h:i:[]):[]) =
  (a * det [[e,f],[h,i]]) - (b * det [[d,f],[g,i]]) + (c * det [[d,e],[g,h]])
det = determinant

norm :: Vector -> Scalar
norm vector =
  sqrt $ sum $ map (**2) vector
magnitude = norm

unitVector :: Vector -> Vector
unitVector vector = map (/ (norm vector)) vector

-- projection :: Vector -> Vector -> Vector
-- projection direction vector =


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

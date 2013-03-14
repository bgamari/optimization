module Optimization.LineSearch
    ( -- * Line search methods
      LineSearch
    , backtrackingSearch
    , newtonSearch
    , secantSearch
    , constantSearch
    ) where

import Linear

-- | A 'LineSearch' method 'search f df p x' determines a step size
-- in direction 'p' from point 'x' for function 'f' with gradient 'df'
type LineSearch f a = (f a -> a) -> (f a -> f a) -> f a -> f a -> a

-- | Armijo condition
armijo :: (Num a, Additive f, Ord a, Metric f)
       => a -> (f a -> a) -> (f a -> f a) -> f a -> f a -> a -> Bool
armijo c1 f df x p a =
    f (x ^+^ a *^ p) <= f x + c1 * a * (df x `dot` p)

-- | Curvature condition
curvature :: (Num a, Ord a, Additive f, Metric f)
          => a -> (f a -> a) -> (f a -> f a) -> f a -> f a -> a -> Bool
curvature c2 f df x p a =
    df (x ^+^ a *^ p) `dot` p >= c2 * (df x `dot` p)

-- | Backtracking line search algorithm
backtrackingSearch :: (Num a, Ord a, Metric f)
                   => a -> a -> LineSearch f a
backtrackingSearch gamma c f df p x =
    head $ dropWhile (not . armijo c f df x p) $ iterate (*gamma) c

-- | Line search by Newton's method
newtonSearch :: (Num a) => LineSearch f a
newtonSearch f df p x = undefined

-- | Line search by secant method with given tolerance
secantSearch :: (Num a, Fractional a) => a -> LineSearch f a
secantSearch eps f df p x = undefined

-- | Constant line search
constantSearch :: a -> LineSearch f a
constantSearch c f df p x = c

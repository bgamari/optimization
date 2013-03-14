module Optimization.LineSearch.SteepestDescent
    ( -- * Steepest descent method
      steepestDescent
    ) where

import Optimization.LineSearch
import Linear

-- | Steepest descent method
{-# INLINEABLE steepestDescent #-}
steepestDescent :: (Num a, Ord a, Additive f, Metric f)
                => LineSearch f a -> (f a -> a) -> (f a -> f a) -> f a -> [f a]
steepestDescent search f df x0 = iterate go x0
  where go x = let p = negated (df x)
                   a = search f df p x
               in x ^+^ a *^ p

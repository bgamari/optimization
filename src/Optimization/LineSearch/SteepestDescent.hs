module Optimization.LineSearch.SteepestDescent
    ( -- * Steepest descent method
      steepestDescent
    ) where

import Optimization.LineSearch
import Linear

-- | Steepest descent method
--
-- @steepestDescent search f df x0@ optimizes a function @f@ with gradient @df@
-- with step size schedule @search@ starting from initial point @x0@
--
-- The steepest descent method chooses the negative gradient of the function
-- as its step direction.
{-# INLINEABLE steepestDescent #-}
steepestDescent :: (Num a, Ord a, Additive f, Metric f)
                => LineSearch f a -> (f a -> f a) -> f a -> [f a]
steepestDescent search df x0 = iterate go x0
  where go x = let p = negated (df x)
                   a = search df p x
               in x ^+^ a *^ p

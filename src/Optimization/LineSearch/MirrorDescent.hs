module Optimization.LineSearch.MirrorDescent
    ( mirrorDescent ) where

import Optimization.LineSearch
import Linear

-- | Mirror descent method.
--
-- Originally described by Nemirovsky and Yudin and later elucidated
-- by Beck and Teboulle, the mirror descent method is a generalization of
-- the projected subgradient method for convex optimization
-- The mirror descent method requires the gradient of a strongly
-- convex function @psi@ (and its dual) as well as a way to get a
-- subgradient for each point of the objective function @f@.
mirrorDescent :: (Num a, Additive f)
              => LineSearch f a -> (f a -> f a) -> (f a -> f a)
              -> (f a -> f a) -> f a -> [f a]
mirrorDescent search dPsi dPsiStar df = go
  where go y0 = let x0 = dPsiStar y0
                    t0 = search df (df x0) x0
                    y1 = dPsi x0 ^-^ t0 *^ df x0
                    x1 = dPsiStar y1
                in x0 : go y1
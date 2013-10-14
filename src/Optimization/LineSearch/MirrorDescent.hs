module Optimization.LineSearch.MirrorDescent
    ( mirrorDescent ) where

import Optimization.LineSearch
import Linear

-- | Mirror descent method.
--
-- Originally described by Nemirovsky and Yudin and later elucidated
-- by Beck and Teboulle, the mirror descent method is a generalization of
-- the projected subgradient method for convex optimization.
-- Mirror descent requires the gradient of a strongly
-- convex function @psi@ (and its dual) as well as a way to get a
-- subgradient for each point of the objective function @f@.
mirrorDescent :: (Num a, Additive f)
              => LineSearch f a  -- ^ line search method
              -> (f a -> f a)    -- ^ strongly convex function, @psi@
              -> (f a -> f a)    -- ^ dual of @psi@
              -> (f a -> f a)    -- ^ gradient of function
              -> f a             -- ^ starting point
              -> [f a]           -- ^ iterates
mirrorDescent search dPsi dPsiStar df = go
  where go y0 = let x0 = dPsiStar y0
                    t0 = search df (df x0) x0
                    y1 = dPsi x0 ^-^ t0 *^ df x0
                    x1 = dPsiStar y1
                in x1 : go y1
{-# INLINEABLE mirrorDescent #-}


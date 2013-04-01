{-# LANGUAGE ScopedTypeVariables #-}

module Optimization.LineSearch.BFGS (bfgs) where

import Linear
import Optimization.LineSearch
import Data.Distributive
import Data.Foldable

-- | Broyden–Fletcher–Goldfarb–Shanno (BFGS) method
-- @bfgs search df b0 x0@ where @b0@ is the inverse Hessian (the
-- identity is often a good initial value).
bfgs :: (Metric f, Additive f, Distributive f, Foldable f, Num a, Fractional a)
     => LineSearch f a -> (f a -> f a) -> f (f a) -> f a -> [f a]
bfgs search df = go
    where go b0 x0 = let p1 = negated $ b0 !* df x0
                         alpha = search df p1 x0
                         s = alpha *^ p1
                         x1 = x0 ^+^ s
                         y = df x1 ^-^ df x0
                         -- Sherman-Morrison update of inverse Hessian
                         sy = s `dot` y
                         u = ((sy + y `dot` (b0 !* y)) *!! outer s s) !!/ (sy^2)
                         v = (b0 !*! outer y s !+! outer s y !*! b0) !!/ sy
                         b1 = b0 !+! u !-! v
                     in x1 : go b1 x1

(!!/) :: (Functor f, Functor g, Fractional a) => f (g a) -> a -> f (g a)
m !!/ s = m !!* recip s

-- | Outer product
outer :: (Functor f, Functor g, Num a) => f a -> g a -> f (g a)
outer a b = fmap (\x->fmap (*x) b) a
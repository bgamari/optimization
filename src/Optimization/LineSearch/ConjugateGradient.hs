module Optimization.LineSearch.ConjugateGradient
    ( -- * Conjugate gradient methods
      conjGrad
      -- * Beta expressions
    , Beta
    , fletcherReeves
    , polakRibiere
    , hestenesStiefel
    ) where

import Optimization.LineSearch
import Linear

-- | A beta expression 'beta df0 df1 p' is an expression for the
-- conjugate direction contribution given the derivative 'df0' and
-- direction 'p' for iteration 'k', 'df1' for iteration 'k+1'
type Beta f a = f a -> f a -> f a -> a

-- | Conjugate gradient method with given beta and line search method
{-# INLINEABLE conjGrad #-}
conjGrad :: (Num a, RealFloat a, Additive f, Metric f)
         => LineSearch f a -> Beta f a
         -> (f a -> a) -> (f a -> f a) -> f a -> [f a]
conjGrad search beta f df x0 = go (negated $ df x0) x0
  where go p x = let a = search f df p x
                     x' = x ^+^ a *^ p
                     b = beta (df x) (df x') p
                     p' = negated (df x') ^+^ b *^ p
                 in x' : go p' x'

-- | Fletcher-Reeves expression for beta
{-# INLINEABLE fletcherReeves #-}
fletcherReeves :: (Num a, RealFloat a, Metric f) => Beta f a
fletcherReeves df0 df1 p0 = norm df1 / norm df0

-- | Polak-Ribiere expression for beta
{-# INLINEABLE polakRibiere #-}
polakRibiere :: (Num a, RealFloat a, Metric f) => Beta f a
polakRibiere df0 df1 p0 = df1 `dot` (df1 ^-^ df0) / norm df0

-- | Hestenes-Stiefel expression for beta
{-# INLINEABLE hestenesStiefel #-}
hestenesStiefel :: (Num a, RealFloat a, Metric f) => Beta f a
hestenesStiefel df0 df1 p0 =
    - (df1 `dot` (df1 ^-^ df0)) / (p0 `dot` (df1 ^-^ df0))
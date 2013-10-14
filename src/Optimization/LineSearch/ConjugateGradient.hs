module Optimization.LineSearch.ConjugateGradient
    ( -- * Conjugate gradient methods
      conjGrad
      -- * Step size methods
    , module Optimization.LineSearch
      -- * Beta expressions
    , Beta
    , fletcherReeves
    , polakRibiere
    , hestenesStiefel
    ) where

import Optimization.LineSearch
import Linear

-- | A beta expression @beta df0 df1 p@ is an expression for the
-- conjugate direction contribution given the derivative @df0@ and
-- direction @p@ for iteration @k@, @df1@ for iteration @k+1@
type Beta f a = f a -> f a -> f a -> a

-- | Conjugate gradient method with given beta and line search method
--
-- The conjugate gradient method avoids the trouble encountered by the
-- steepest descent method on poorly conditioned problems (e.g. those with
-- a wide range of eigenvalues). It does this by choosing directions which
-- satisfy a condition of @A@ orthogonality, ensuring that steps in the
-- "unstretched" search space are orthogonal.
-- TODO: clarify explanation
conjGrad :: (Num a, RealFloat a, Additive f, Metric f)
         => LineSearch f a  -- ^ line search method
         -> Beta f a        -- ^ beta expression
         -> (f a -> f a)    -- ^ gradient of function
         -> f a             -- ^ starting point
         -> [f a]           -- ^ iterates
conjGrad search beta df x0 = go (negated $ df x0) x0
  where go p x = let a = search df p x
                     x' = x ^+^ a *^ p
                     b = beta (df x) (df x') p
                     p' = negated (df x') ^+^ b *^ p
                 in x' : go p' x'
{-# INLINEABLE conjGrad #-}

-- | Fletcher-Reeves expression for beta
fletcherReeves :: (Num a, RealFloat a, Metric f) => Beta f a
fletcherReeves df0 df1 _ = norm df1 / norm df0
{-# INLINEABLE fletcherReeves #-}

-- | Polak-Ribiere expression for beta
polakRibiere :: (Num a, RealFloat a, Metric f) => Beta f a
polakRibiere df0 df1 _ = df1 `dot` (df1 ^-^ df0) / norm df0
{-# INLINEABLE polakRibiere #-}

-- | Hestenes-Stiefel expression for beta
hestenesStiefel :: (Num a, RealFloat a, Metric f) => Beta f a
hestenesStiefel df0 df1 p0 =
    - (df1 `dot` (df1 ^-^ df0)) / (p0 `dot` (df1 ^-^ df0))
{-# INLINEABLE hestenesStiefel #-}

module Optimization.TrustRegion.Newton
    ( -- * Newton's method
      newton
      -- * Matrix inversion methods
    , bicInv
    , bicInv'
    ) where

import Control.Applicative
import Data.Distributive (Distributive)
import Data.Functor.Bind (Apply)
import Data.Foldable (Foldable)
import Linear

-- | Newton's method
newton :: (Num a, Ord a, Additive f, Metric f, Foldable f)
       => (f a -> f a) -> (f a -> f (f a)) -> f a -> [f a]
newton df ddfInv x0 = iterate go x0
  where go x = x ^-^ ddfInv x !* df x
{-# INLINEABLE newton #-}

-- | Inverse by iterative method of Ben-Israel and Cohen
-- with given starting condition
bicInv' :: (Functor m, Distributive m, Additive m,
            Applicative m, Apply m, Foldable m, Conjugate a)
        => m (m a) -> m (m a) -> [m (m a)]
bicInv' a0 a = iterate go a0
  where go ak = 2 *!! ak !-! ak !*! a !*! ak
{-# INLINEABLE bicInv' #-}

-- | Inverse by iterative method of Ben-Israel and Cohen
-- starting from @alpha A^T@. Alpha should be set such that
-- 0 < alpha < 2/sigma^2 where sigma denotes the largest singular
-- value of A
bicInv :: (Functor m, Distributive m, Additive m,
           Applicative m, Apply m, Foldable m, Conjugate a)
       => a -> m (m a) -> [m (m a)]
bicInv alpha a = bicInv' (alpha *!! adjoint a) a
{-# INLINEABLE bicInv #-}


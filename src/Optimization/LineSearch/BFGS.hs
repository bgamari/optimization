{-# LANGUAGE ScopedTypeVariables #-}

module Optimization.LineSearch.BFGS
    ( -- * Broyden-Fletcher-Goldfarb-Shanna (BFGS) method
      bfgs
      -- * Step size methods
    , module Optimization.LineSearch
    ) where

import Linear
import Optimization.LineSearch
import Control.Applicative
import Data.Traversable
import Data.Foldable

-- | Broyden–Fletcher–Goldfarb–Shanno (BFGS) method
-- @bfgs search df b0 x0@ where @b0@ is the inverse Hessian (the
-- identity is often a good initial value).
bfgs :: ( Metric f, Additive f, Foldable f, Traversable f, Applicative f
        , Fractional a, Epsilon a)
     => LineSearch f a   -- ^ line search method
     -> (f a -> f a)     -- ^ gradient of function
     -> f (f a)          -- ^ inverse Hessian at starting point
     -> f a              -- ^ starting point
     -> [f a]            -- ^ iterates
bfgs search df = go
    where go b0 x0 = let p1 = negated $ b0 !* df x0
                         alpha = search df p1 x0
                         s = alpha *^ p1
                         x1 = x0 ^+^ s
                         y = df x1 ^-^ df x0
                         -- Sherman-Morrison update of inverse Hessian
                         sy = s `dot` y
                         rho = if nearZero sy then 1000 else 1 / sy
                         i = kronecker (pure 1)
                         u = i !-! rho *!! outer y s
                         v = i !-! rho *!! outer s y
                         b1 = u !*! b0 !*! v !+! rho *!! outer s s
                     in x1 : go b1 x1
{-# INLINABLE bfgs #-}    


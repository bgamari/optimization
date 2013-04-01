module Optimization.TrustRegion.Nesterov1983
    ( -- * Nesterov's Optimal Gradient method
      optimalGradient
    ) where

import Linear

-- | Nesterov 1983
-- @optimalGradient kappa l df alpha0 x0@ is Nesterov's optimal
-- gradient method, first described in 1983. This method requires
-- knowledge of the Lipschitz constant @l@ of the gradient, the condition
-- number @kappa@, as well as an initial step size @alpha0@ in @(0,1)@.
{-# INLINEABLE optimalGradient #-}
optimalGradient :: (Additive f, Functor f, Ord a, Floating a, Epsilon a)
                => a -> a -> (f a -> f a) -> a -> f a -> [f a]
optimalGradient kappa l df a0' x0' = go a0' x0' x0'
  where go a0 x0 y0 = let x1 = y0 ^-^ df y0 ^/ l
                          alphas = quadratic 1 (a0^2 - 1/kappa) (-a0^2)
                          a1 = case filter (\x->x >= 0 && x <= 1) alphas of
                                 a:_  -> a
                                 []   -> error "No solution for alpha_{k+1}"
                          b1 = a0 * (1 - a0) / (a0^2 + a1)
                          y1 = x1 ^+^ b1 *^ (x1 ^-^ x0)
                      in x1 : go a0 x1 y1

-- | 'quadratic a b c' is the real solutions to a quadratic equation
-- 'a x^2 + b x + c == 0'
quadratic :: (Ord a, Floating a, Epsilon a)
          => a -> a -> a -> [a]
quadratic a b c
    | discr < 0      = []
    | nearZero discr = [-b / 2 / a]
    | otherwise      = [ (-b + sqrt discr) / 2 / a
                       , (-b - sqrt discr) / 2 / a ]
  where discr = b^2 - 4*a*c

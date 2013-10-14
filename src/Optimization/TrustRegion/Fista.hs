module Optimization.TrustRegion.Fista
    ( -- * Fast Iterative Shrinkage-Thresholding Algorithm
      fista
    ) where

import Linear

-- | Fast Iterative Shrinkage-Thresholding Algorithm (FISTA) with
-- constant step size
fista :: (Additive f, Fractional a, Floating a)
      => a             -- ^ Lipschitz constant, @l@
      -> (f a -> f a)  -- ^ gradient of function
      -> f a           -- ^ starting point
      -> [f a]         -- ^ iterates
fista l df x0' = go x0' x0' 1
  where go x0 y1 t1 = let x1 = y1 ^-^ df y1 ^/ l
                          t2 = (1 + sqrt (1 + 4 * t1^2)) / 2
                          y2 = x1 ^+^ (t1-1) / t2 *^ (x1 ^-^ x0)
                      in x1 : go x1 y2 t2
{-# INLINEABLE fista #-}

module Optimization.LineSearch.BarzilaiBorwein
    ( barzilaiBorwein
    ) where

import Linear

-- | Barzilai-Borwein 1988 is a non-monotonic optimization method
barzilaiBorwein :: (Additive f, Metric f, Functor f, Fractional a)
                => (f a -> a) -> (f a -> f a) -> f a -> f a -> [f a]
barzilaiBorwein f df x0 x1 = go (x0, x1)
  where go (x0,x1) = let s = x1 ^-^ x0
                         z = df x1 ^-^ df x0
                         alpha = (s `dot` z) / (z `dot` z)
                         x2 = x1 ^-^ alpha *^ df x1
                     in x2 : go (x1, x2)
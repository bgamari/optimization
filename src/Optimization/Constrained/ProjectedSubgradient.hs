module Optimization.Constrained.ProjectedSubgradient
    ( -- * Projected subgradient method
      linearProjSubgrad
      -- * Step schedules
    , StepSched
    , optimalStepSched
    , constStepSched
    , sqrtKStepSched
    , invKStepSched
    ) where

import Linear
import Data.Traversable

-- | A step size schedule
type StepSched f a = [f a -> a -> a]

-- | @linearProjSubgrad stepSizes proj a b x0@ minimizes the objective @A
-- x - b@ potentially projecting iterates into a feasible space with
-- @proj@ with the given step-size schedule
linearProjSubgrad :: (Additive f, Traversable f, Metric f, Ord a, Fractional a)
                  => StepSched f a  -- ^ A step size schedule
                  -> (f a -> f a)   -- ^ Function projecting into the feasible space
                  -> f a            -- ^ Coefficient vector @A@ of objective
                  -> a              -- ^ Constant @b@ of objective
                  -> f a            -- ^ Initial solution
                  -> [f a]
linearProjSubgrad stepSizes proj a b = go stepSizes
  where go (alpha:rest) x0 =
            let p = negated $ df x0
                step = alpha a (f x0)
                x1 = proj $ x0 ^+^ step *^ p
            in x1 : go rest x1
        go [] _ = []
        df _ = a
        f x = a `dot` x - b

-- | The optimal step size schedule when the optimal value of the
-- objective is known
optimalStepSched :: (Fractional a, Metric f)
                 => a    -- ^ The optimal value of the objective
                 -> StepSched f a
optimalStepSched fStar =
    repeat $ \gk fk->(fk - fStar) / quadrance gk

-- | Constant step size schedule
constStepSched :: a    -- ^ The step size
               -> StepSched f a
constStepSched gamma =
    repeat $ \_ _ -> gamma

-- | A non-summable step size schedule
sqrtKStepSched :: Floating a
               => a       -- ^ The size of the first step
               -> StepSched f a
sqrtKStepSched gamma =
    map (\k _ _ -> gamma / sqrt (fromIntegral k)) [0..]

-- | A square-summable step size schedule
invKStepSched :: Fractional a
              => a        -- ^ The size of the first step
              -> StepSched f a
invKStepSched gamma =
    map (\k _ _ -> gamma / fromIntegral k) [0..]

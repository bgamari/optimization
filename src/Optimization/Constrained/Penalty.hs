{-# LANGUAGE CPP, DeriveFunctor, DeriveFoldable, DeriveTraversable, DeriveGeneric,
    FlexibleInstances, FlexibleContexts, TypeFamilies,
    KindSignatures, DataKinds, TypeOperators, RankNTypes, ExistentialQuantification #-}

module Optimization.Constrained.Penalty
  ( -- * Building the problem
    Opt
  , FU(..)
  , optimize
  , constrainEQ
  , constrainLT
  , constrainGT
    -- * Optimizing the problem
  , minimize
  , maximize
    -- * Finding the Lagrangian
  , lagrangian
  ) where

#if MIN_VERSION_ad(4,0,0)
import           Numeric.AD.Mode
#else
import           Numeric.AD
#endif
import qualified Data.Vector as V

#if MIN_VERSION_ad(4,0,0)
newtype FU f a = FU { runFU :: forall t. (Mode t, Fractional t, Ord t) => f t -> t}
#else
newtype FU f a = FU { runFU :: forall s. Mode s => f (AD s a) -> AD s a }
#endif

type V = V.Vector

-- | @Opt d f gs hs@ is a Lagrangian optimization problem with objective @f@
-- equality (@g(x) == 0@) constraints @gs@ and less-than (@h(x) < 0@)
-- constraints @hs@
data Opt f a = Opt (FU f a) (V (FU f a)) (V (FU f a))

#if MIN_VERSION_ad(4,0,0)
optimize :: (forall t . Mode t => f t -> t) -> Opt f a
#else
optimize :: (forall s. Mode s => f (AD s a) -> AD s a) -> Opt f a
#endif
optimize f = Opt (FU f) V.empty V.empty

augment :: a -> V a -> V a
augment a xs = V.cons a xs

#if MIN_VERSION_ad(4,0,0)
constrainEQ :: (forall t. Mode t => f t -> t)
#else
constrainEQ :: (forall s. Mode s => f (AD s a) -> AD s a)
#endif
            -> Opt f a -> Opt f a
constrainEQ g (Opt f gs hs) = Opt f (augment (FU g) gs) hs

#if MIN_VERSION_ad(4,0,0)
constrainLT :: (forall t. Mode t => f t -> t)
#else
constrainLT :: (forall s. Mode s => f (AD s a) -> AD s a)
#endif
            -> Opt f a -> Opt f a
constrainLT h (Opt f gs hs) = Opt f gs (augment (FU h) hs)

#if MIN_VERSION_ad(4,0,0)
constrainGT :: (Num a) => (forall t. Mode t => f t -> t)
#else
constrainGT :: (Num a) => (forall s. Mode s => f (AD s a) -> AD s a)
#endif
            -> Opt f a -> Opt f a
constrainGT h (Opt f gs hs) = Opt f gs (augment (FU $ negate . h) hs)

-- | Minimize the given constrained optimization problem
-- This is a basic penalty method approach
minimize :: (Functor f, RealFrac a, Ord a, g ~ V)
         => (FU f a -> f a -> [f a])   -- ^ Primal minimizer
         -> Opt f a                    -- ^ The optimization problem of interest
         -> a                          -- ^ The penalty increase factor
         -> f a                        -- ^ The primal starting value
         -> g a                        -- ^ The dual starting value
         -> [f a]                      -- ^ Optimizing iterates
minimize minX opt alpha = go
  where go x0 l0 = let l1 = fmap (*alpha) l0
                       x1 = head $ drop 100 $ minX (FU $ \x -> augLagrangian opt x (fmap realToFrac l1)) x0
                   in x1 : go x1 l1
{-# INLINEABLE minimize #-}

-- | Maximize the given constrained optimization problem
maximize :: (Functor f, RealFrac a, Ord a, g ~ V)
         => (FU f a -> f a -> [f a])   -- ^ Primal minimizer
         -> Opt f a                    -- ^ The optimization problem of interest
         -> a                          -- ^ The penalty increase factor
         -> f a                        -- ^ The primal starting value
         -> g a                        -- ^ The dual starting value
         -> [f a]                      -- ^ Optimizing iterates
maximize minX (Opt (FU f) gs hs) alpha =
    minimize minX (Opt (FU $ negate . f) gs hs) alpha
{-# INLINEABLE maximize #-}

-- | The Lagrangian for the given constrained optimization
lagrangian :: (Num a) => Opt f a
#if MIN_VERSION_ad(4,0,0)
           -> (forall t. (Mode t, Fractional t, Ord t)=> f t -> V t -> t)
#else
           -> (forall s. Mode s => f (AD s a) -> V (AD s a) -> AD s a)
#endif
lagrangian (Opt (FU f) gs hs) x l =
    f x - V.sum (V.zipWith (\lamb (FU g)->lamb * g x) l gs)
{-# INLINEABLE lagrangian #-}

-- | The augmented Lagrangian for the given constrained optimization
augLagrangian :: (Num a, Ord a) => Opt f a
#if MIN_VERSION_ad(4,0,0)
           -> (forall t. (Mode t, Fractional t, Ord t) => f t -> V t -> t)
#else
           -> (forall s. Mode s => f (AD s a) -> V (AD s a) -> AD s a)
#endif
augLagrangian (Opt (FU f) gs hs) x l =
    f x + V.sum (V.zipWith (*) l $ V.concat [gs', hs'])
  where gs' = V.map (\(FU g) -> (g x)^(2 :: Int)) gs
        hs' = V.map (\(FU h) -> (max 0 $ h x)^(2 :: Int)) hs
{-# INLINEABLE augLagrangian #-}

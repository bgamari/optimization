-- |
-- Module      : Optimization.LineSearch
-- Copyright   : (c) 2012-2013 Ben Gamari
-- License     : BSD-style (see the file LICENSE)
-- Maintainer  : Ben Gamari <bgamari@gmail.com>
-- Stability   : provisional
-- Portability : portable
--
-- Line search algorithms are a class of iterative optimization
-- methods. These methods are distinguished by the characteristic of,
-- starting from a point @x0@, choosing a direction @d@ (by some method)
-- to advance and then finding an optimal distance @a@ (known as the
-- step-size) to advance in this direction.
--
-- Here we provide several methods for determining this optimal
-- distance. These can be used with any of line-search optimization
-- algorithms found in this namespace.
--

module Optimization.LineSearch
    ( -- * Line search methods
      LineSearch
    , backtrackingSearch
      -- * Other line search methods
    , constantSearch
    , newtonSearch
    , secantSearch
      -- * Wolfe conditions
      -- Nocedal gives typical values of 10^-4 for @c1@ and 0.9 for
      -- @c2@
    , armijoSearch
    , wolfeSearch
    ) where

import Prelude hiding (pred)
import Linear

-- | A line search method @search df p x@ determines a step size
-- in direction @p@ from point @x@ for function @f@ with gradient @df@
type LineSearch f a = (f a -> f a)    -- ^ gradient of function
                   -> f a             -- ^ search direction
                   -> f a             -- ^ starting point
                   -> a               -- ^ step size

-- | Armijo condition
--
-- The Armijo condition captures the intuition that step should
-- move far enough from its starting point to change the function enough,
-- as predicted by its gradient. This often finds its place as a criterion
-- for line-search
armijo :: (Num a, Additive f, Ord a, Metric f)
       => a            -- ^ Armijo condition strength
       -> (f a -> a)   -- ^ function value
       -> (f a -> f a) -- ^ gradient of function
       -> f a          -- ^ point to evaulate at
       -> f a          -- ^ search direction
       -> a            -- ^ search step size
       -> Bool         -- ^ is Armijo condition satisfied?
armijo c1 f df x p a =
    f (x ^+^ a *^ p) <= f x + c1 * a * (df x `dot` p)

-- | Curvature condition
curvature :: (Num a, Ord a, Additive f, Metric f)
          => a             -- ^ curvature condition strength c2
          -> (f a -> f a)  -- ^ gradient of function
          -> f a           -- ^ point to evaluate at
          -> f a           -- ^ search direction
          -> a             -- ^ search step size
          -> Bool          -- ^ is curvature condition satisfied
curvature c2 df x p a =
    df (x ^+^ a *^ p) `dot` p >= c2 * (df x `dot` p)

-- | Backtracking line search algorithm
--
-- This is a building block for line search algorithms which reduces
-- its step size until the given condition is satisfied.
-- 
-- @backtrackingSearch gamma alpha pred@ starts with the given step
-- size @alpha@ and reduces it by a factor of @gamma@ until the given
-- condition is satisfied.
backtrackingSearch :: (Num a, Ord a, Metric f)
                   => a            -- ^ step size reduction factor gamma
                   -> a            -- ^ initial step size alpha
                   -> (a -> Bool)  -- ^ search condition
                   -> LineSearch f a
backtrackingSearch gamma alpha pred _ _ _ =
    head $ dropWhile (not . pred) $ nonzero $ iterate (*gamma) alpha
  where nonzero (x:xs) | not $ x > 0 = error "Backtracking search failed: alpha=0" -- FIXME
                       | otherwise   = x : nonzero xs
        nonzero [] = error "Backtracking search failed: no more iterates"
{-# INLINEABLE backtrackingSearch #-}

-- | Armijo backtracking line search algorithm
--
-- @armijoSearch gamma alpha c1@ starts with the given step size @alpha@
-- and reduces it by a factor of @gamma@ until the Armijo condition
-- is satisfied.
armijoSearch :: (Num a, Ord a, Metric f)
             => a                   -- ^ step size reduction factor gamma
             -> a                   -- ^ initial step size alpha
             -> a                   -- ^ Armijo condition strength c1
             -> (f a -> a)          -- ^ function value
             -> LineSearch f a
armijoSearch gamma alpha c1 f df p x =
    backtrackingSearch gamma alpha (armijo c1 f df x p) df p x
{-# INLINEABLE armijoSearch #-}

-- | Wolfe backtracking line search algorithm (satisfies both Armijo and
-- curvature conditions)
--
-- @wolfeSearch gamma alpha c1@ starts with the given step size @alpha@
-- and reduces it by a factor of @gamma@ until both the Armijo and
-- curvature conditions is satisfied.
wolfeSearch :: (Num a, Ord a, Metric f)
             => a                   -- ^ step size reduction factor gamma
             -> a                   -- ^ initial step size alpha
             -> a                   -- ^ Armijo condition strength c1
             -> a                   -- ^ curvature condition strength c2
             -> (f a -> a)          -- ^ function value
             -> LineSearch f a
wolfeSearch gamma alpha c1 c2 f df p x =
    backtrackingSearch gamma alpha wolfe df p x
  where wolfe a = armijo c1 f df p x a && curvature c2 df x p a
{-# INLINEABLE wolfeSearch #-}

-- | Line search by Newton's method
newtonSearch :: (Num a) => LineSearch f a
newtonSearch = undefined

-- | Line search by secant method with given tolerance
secantSearch :: (Num a, Fractional a) => a -> LineSearch f a
secantSearch = undefined

-- | Constant line search
--
-- @constantSearch c@ always chooses a step-size @c@.
constantSearch :: a                 -- ^ step size
               -> LineSearch f a
constantSearch c _ _ _ = c
{-# INLINEABLE constantSearch #-}

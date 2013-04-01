import Linear
import Numeric.AD
import Control.Monad
import Optimization.LineSearch.ConjugateGradient
import Optimization.LineSearch.SteepestDescent
import Optimization.LineSearch.BFGS

-- | Rosenbrock function
rosenbrock :: Num a => V2 a -> a
rosenbrock (V2 x y) = (1-x)^2 + 100*(y-x^2)^2

main = do
    let f = rosenbrock
        df = grad rosenbrock :: V2 Double -> V2 Double
        x0 = V2 2 2
        search = backtrackingSearch 0.1 0.2 f
        beta = fletcherReeves
    forM_ (take 10 $ conjGrad search beta f df x0) $ \x->do print (x, f x)
    forM_ (take 10 $ steepestDescent search f df x0) $ \x->do print (x, f x)
    forM_ (take 10000 $ bfgs search df eye2 x0) $ \x->do print (x, f x)

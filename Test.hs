import Linear
import Numeric.AD
import Control.Monad
import Optimization.LineSearch.ConjugateGradient
import Optimization.LineSearch.SteepestDescent
import Optimization.LineSearch.BFGS
import Optimization.LineSearch.BarzilaiBorwein

-- | Rosenbrock function
rosenbrock :: Num a => V2 a -> a
rosenbrock (V2 x y) = (1-x)^2 + 100*(y-x^2)^2

main = do
    let f = rosenbrock
        df = grad rosenbrock :: V2 Double -> V2 Double
        x0 = V2 2 2
        search = armijoSearch 0.1 0.2 0.2 f
        beta = fletcherReeves

    putStrLn "Conjugate gradient"
    forM_ (take 10 $ conjGrad search beta df x0) $ \x->do print (x, f x)
    putStrLn "Steepest descent"
    forM_ (take 10 $ steepestDescent search df x0) $ \x->do print (x, f x)
    putStrLn "BFGS"
    forM_ (take 10000 $ bfgs search df eye2 x0) $ \x->do print (x, f x)
    putStrLn "Barzilai-Borwein"
    forM_ (take 100 $ barzilaiBorwein df (V2 3 3) x0) $ \x->do print (x, f x)

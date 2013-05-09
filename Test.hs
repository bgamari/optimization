import Linear
import Numeric.AD
import Control.Applicative
import Control.Monad
import Optimization.LineSearch.ConjugateGradient
import Optimization.LineSearch.SteepestDescent
import Optimization.LineSearch.BFGS
import Optimization.LineSearch.BarzilaiBorwein

-- | Rosenbrock's parabolic valley function (Rosenbrock 1960)
rosenbrock :: Num a => V2 a -> a
rosenbrock (V2 x y) = (1-x)^2 + 100*(y-x^2)^2

rosenbrockX = V2 (-1.2) 1

-- | Powell's quartic function (Powell 1962)
powell :: (Num a) => V4 a -> a
powell (V4 x y z t) =
    (x + 10*y)^2 + 5*(z-t)^2 + (y-2*z)^4 + 10*(x-t)^4

powellX = V4 3 (-1) 0 1

-- | Fletcher and Powell's helical valley (Fletcher and Powell 1963)
helical :: (RealFloat a) => V3 a -> a
helical (V3 x y z) =
    100*(z - 10 * atan2 x y / 2 / pi)^2 + (sqrt (x^2 + y^2) - 1)^2 + z^2

helicalX = V3 (-1) 0 0

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
    forM_ (take 100 $ barzilaiBorwein df (pure 3 3) x0) $ \x->do print (x, f x)

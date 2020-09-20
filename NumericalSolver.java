// This numerical solver, attempts to find a value of x such that 
// the function fn(x) = target
// The mothod solveMonotonicFn() assumes that the function is monotonic.
// On the other hand: solve() is more forgiving, but it is not perfect.
// For exponential functions, solve() is slow.
public class NumericalSolver {

    private double             m_target;
    private double             m_startX0; // startX0 and startX1 don't necessarily need to bound the target.
    private double             m_startX1; // i.e. it is not required for: target - fn(startX0) to be the same sign as fn(startX1) - target 
    private double             m_tol;
    private int                m_maxSteps;
    private IfaceFnOneVariable m_fn;
    
    private double             m_adjustment;
    
    ////////////////////////////////////////////////////////
    // The following solver tries to be reasonably forgiving 
    // the startX1 and startX2 don't need to bound the solution. 
    public NumericalSolver(double target, double startX0, double startX1, 
                           double tol,    int maxSteps,     IfaceFnOneVariable fn){
        m_target     = target;
        m_startX0    = startX0;
        m_startX1    = startX1;
        m_tol        = tol;
        m_maxSteps   = maxSteps;
        m_fn         = fn;
        
        m_adjustment = 0.5;
    }
    /////////////////////////////////////////////////////////
    // The following method assums that fn(.) is monotonic, i.e. only increasing or only decreasing.
    public double solveMonotonicFn(){
        
        double[] x = new double[2];
        double[] y = new double[2];
        
        x[0] = Math.min(m_startX0, m_startX1);
        y[0] = m_fn.fn(x[0]);

        if(Math.abs(y[0] - m_target) < m_tol){
            System.out.println("On first iteration, found solution: " + Double.toString(x[0]));
            return (x[0]);
        }
        
        x[1] = Math.max(m_startX0, m_startX1);        
        y[1] = m_fn.fn(x[1]);
        
        if(Math.abs(y[1] - m_target) < m_tol){
                System.out.println("On first iteration, found solution: " + Double.toString(x[1]));
                return (x[1]);
        } else if (!findUpperAndLowerBound(x,y,0)){
            System.out.println("NumericalSolver.solveMonotonicFn() was unable to find a solution.");
            return x[0];
        } else {
            return findBoundedSolution(x,y);
        }
    }
    ////////////////////////////////////////////////////////
    public double findBoundedSolution(double[] x, double[] y) {
        int    counter    = 0;
        double xSolution  = 0.0;
        double ySolution  = 0.0;
        
        boolean finished = false;
        while(!finished){
            counter++;
            
            if( Math.abs(y[0] - m_target) < m_tol){
                xSolution = x[0];
                ySolution = y[0];
                finished  = true;
                    
            } else if( Math.abs(y[1] - m_target) < m_tol){
                xSolution = x[1];
                ySolution = y[1];
                finished  = true;
            } else if (counter > m_maxSteps){
                System.out.println("After " + Integer.toString(counter) + " iterations did not find solution close to the target");
                return x[0];
            } else {
                double nextX =  x[0] +(m_target - y[0]) * ( x[1] - x[0]) / (y[1] - y[0]);
                double nextY =  m_fn.fn(nextX);
                setNextXY(nextX, nextY, x,y);
                
                nextX = 0.5 * ( x[0] + x[1]);
                nextY = m_fn.fn(nextX);                
                setNextXY(nextX, nextY, x,y);
            }
        }
        
        System.out.println("After " + Integer.toString(counter) + " iterations, found solution with x = "
            + Double.toString(xSolution) + " and y = " + Double.toString(ySolution));

        return xSolution;
}
    ////////////////////////////////////////////////////////
    public void setNextXY(double nextX, double nextY, double[] x, double[] y){
        if ( Math.signum( nextY - m_target ) == Math.signum( y[0] - m_target)){
            x[0] = nextX;
            y[0] = nextY;
        } else {
            x[1] = nextX;
            y[1] = nextY;
        }
    }
    ////////////////////////////////////////////////////////
    // The following function is recursive.
    public boolean findUpperAndLowerBound(double[] x, double[] y, int itteration){ 
        
        if ( Math.signum(y[0] - m_target) != Math.signum(y[0] - m_target))
            return true;
        else if( itteration > m_maxSteps){
            System.out.println("Could not bound solution.");
            return false;
        }else if ( y[0] == y[1]) {
            boolean foundDifferentYs = false;
            int     counter          = 0;
            while( ! foundDifferentYs){
                counter++;
                x[0] = adjustDown(x[0]);
                
                y[0] = m_fn.fn(x[0]);

                if(y[0] != y[1]){
                    foundDifferentYs = true;
                } else {
                
                    x[1] = adjustUp  (x[1]);                
                    y[1] = m_fn.fn(x[1]);
                    
                    if(y[0] != y[1]){
                        foundDifferentYs = true;
                    } else if(counter >= m_maxSteps){
                            System.out.println("Could not find distinct function values");
                            return false;
                    } // end of nested if ... else    
                }     // end of if ... else  
            }         // end of while 
        } else {
            double nextX = x[0] +(m_target - y[0]) * ( x[1] - x[0]) / (y[1] - y[0]);
            double nextY = m_fn.fn(nextX);
            
            if( y[1] > y[0]){
                x[0] = x[1];
                x[1] = nextX;
                
                y[0] = y[1];
                y[1] = nextY;
            } else if( y[1] < y[0]){
                x[1] = x[0];
                x[0] = nextX;
                
                y[1] = y[0];
                y[0] = nextY;
            }
            
            if (      Math.signum(y[1] - m_target) != Math.signum(y[0] - m_target))
                return true;
            else if ( Math.abs   (y[0] - m_target)    < Math.abs( y[1] - m_target)){
                x[0] = adjustDown(x[0]);
                y[0] = m_fn.fn(   x[0]);
                if ( Math.signum(y[1] - m_target) != Math.signum(y[0] - m_target))
                    return true;
            } else {
                x[1] = adjustUp(x[1]);
                y[1] = m_fn.fn( x[1]);
                
                if ( Math.signum(y[1] - m_target) != Math.signum(y[0] - m_target))
                    return true;
                else 
                    return findUpperAndLowerBound( x, y, itteration+1);
                
            }            
        }        
        return false;
    }
    ////////////////////////////////////////////////////////
    public double adjustDown(double z){
        if ( z >= 0)
            return ( z * ( 1 - m_adjustment) - m_adjustment);
        else             
            return ( z * ( 1 + m_adjustment) - m_adjustment);
    }
    ////////////////////////////////////////////////////////
    public double adjustUp(double z){
        if ( z >= 0)
            return ( z * ( 1 + m_adjustment) + m_adjustment);
        else             
            return ( z * ( 1 - m_adjustment) + m_adjustment);
    }
    ////////////////////////////////////////////////////////
    public double solve(){ // doesn't work very well with exponentials.
        // The following method does a linear interpolation / extrapolation
        // attempting to find the x that obtains  abs( target - fn(x)) < tol
        
        int     counter  = 0;
        
        double  x0 = m_startX1;
        double  x1 = m_startX0;
        
        if(x0 == x1) {
            // we'll move x2 away from zero.
            x1 = x1 * 1.1 + 0.1 * Math.signum(x1) + 0.01;
        } 

        double  scaleFactor = 1.0;
        if ( Math.abs(m_target) > 1)
            scaleFactor = 1.0 / Math.abs(m_target);
        
        while ( true){
            
            counter++;

            double  y0 = m_fn.fn(x0);
            double  y1 = m_fn.fn(x1);
            if ( Math.abs( y0 - m_target) * scaleFactor < m_tol){
                System.out.println(  "After " + Integer.toString(counter) 
                           + " iterations, for x0 = " + Double.toString(x0) + " found y = " 
                           + Double.toString(y0) + " to be within the tolerance of the target."); 
                 return (x0); 
            } else if ( counter >= m_maxSteps){
                
                System.out.println("After " + Integer.toString(counter) 
                                   + " iterations, with x value: " + Double.toString(x0) 
                                   + " had y value: "  + Double.toString(y0) 
                                   + " which was not within the tolerance.");
                return x0;
            } else if ( y0 == y1){
                x1 = x1 * 1.1 + Math.signum(x1) * 0.1 + 0.01;
            } else {
                // If we use a linear approximation to interpolate or extrapolate, 
                // either way, we have matched slopes.            
                // (m_target - y1)/(nextX - x1) = (y2 - y1)/(x2 - x1)
                double nextX = x0 +(m_target - y0) * ( x1 - x0) / ( y1 - y0);
            
                double nextY = m_fn.fn(nextX);
 

                if (Math.abs(nextY - m_target) * scaleFactor < m_tol) {
                    
                    System.out.println(  "After " + Integer.toString(counter) 
                                       + " iterations, for x = " + Double.toString(nextX) + " found y = " 
                                       + Double.toString(nextY) + " to be within the tolerance of the target."); 
                    return nextX;
                    
                } else if ( Math.abs(y0 - m_target ) < Math.abs(y1 - m_target)){  
                //} else if ( (counter % 2) == 0 ) {
                    x0  = x0 + ( nextX - x0 ) * 0.5;
                    x1 = nextX;
                } else {
                    x0 = nextX;
                    x1  = x1 + ( nextX - x1 ) * 0.5;
                } // end of nested if else
            }     // end of if else
         }         // end of while 
    }             // end of method solve(.)

    ////////////////////////////////////////////////
    public interface IfaceFnOneVariable {
        double fn(double x);
    }    
    //////////////////////////////////////////////////////////////
    public static void testNumericalSolver(){
        sinTest();
        quadraticTest();
        exponentialTest();
        exponentialMonotonicTest();
    }
    //////////////////////////////////////////////////////////////
    public static void sinTest(){
        
        IfaceFnOneVariable sinFn = (new IfaceFnOneVariable()
              {public double fn(double x) {return ( Math.sin(x));}  } );
        
        double target    = sinFn.fn(Math.PI * 0.25);
        double startX0   = 0.0;
        double startX1   = 0.1;
        double tol       = 1e-5;
        int    maxSteps  = 50;
        
        NumericalSolver ns = new NumericalSolver(target, startX0, startX1, tol, maxSteps, sinFn);
        
        double xSoln = ns.solve();
        System.out.println("Found sin(" + Double.toString(xSoln) 
                            + ") is close to target: " + Double.toString(target));      
    }
    //////////////////////////////////////////////////////////////
    public static void quadraticTest(){
        
        IfaceFnOneVariable quadratic = (new IfaceFnOneVariable()
              {public double fn(double x) {return ( 2.0 * x * x + 4 * x + 3);}  } );
        
        double target    = quadratic.fn(4.0);
        double startX0   = 0.0;
        double startX1   = 1.0;
        double tol       = 1e-5;
        int    maxSteps  = 50;

        NumericalSolver ns = new NumericalSolver(target, startX0, startX1, tol, maxSteps, quadratic);
        
        double xSoln = ns.solve();
        System.out.println("Found quadratic(" + Double.toString(xSoln) 
                            + ") is close to target: " + Double.toString(target)); 
    }
    //////////////////////////////////////////////////////////////
    public static void exponentialTest(){
        
        IfaceFnOneVariable exponential = (new IfaceFnOneVariable()
              {public double fn(double x) {return ( Math.exp(x));}  } );
        
        double target    =  exponential.fn(4.0);
        double startX0   = -1.0;
        double startX1   =  1.0;
        double tol       =  1e-5;
        int    maxSteps  =  50;

        NumericalSolver ns = new NumericalSolver(target, startX0, startX1, tol, maxSteps, exponential);
        
        double xSoln = ns.solve();
        System.out.println("Found exponential(" + Double.toString(xSoln) 
                            + ") is close to target: " + Double.toString(target)); 
                
    }
    //////////////////////////////////////////////////////////////
    public static void exponentialMonotonicTest(){
        
        IfaceFnOneVariable exponential = (new IfaceFnOneVariable()
              {public double fn(double x) {return ( Math.exp(x));}  } );
        
        double target    =  exponential.fn(4.0);
        double startX0   = -1.0;
        double startX1   =  1.0;
        double tol       =  1e-5;
        int    maxSteps  =  50;

        NumericalSolver ns = new NumericalSolver(target, startX0, startX1, tol, maxSteps, exponential);
        
        double xSoln = ns.solveMonotonicFn();
        System.out.println("Using monotonic solver, found exponential(" + Double.toString(xSoln) 
                            + ") is close to target: " + Double.toString(target)); 
                
    }
     //////////////////////////////////////////////////////////////
}

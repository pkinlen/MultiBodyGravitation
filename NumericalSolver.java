
public class NumericalSolver {

	private double             m_target;
	private double             m_startX1;
	private double             m_startX2;
	private double             m_tol;
	private int                m_maxSteps;
	private IfaceFnOneVariable m_fn;
	
	////////////////////////////////////////////////////////
	// The following solver tries to be reasonably forgiving 
	// the startX1 and startX2 don't need to bound the solution. 
	public NumericalSolver(double target, double startX1, double startX2, 
			               double tol,    int maxSteps,     IfaceFnOneVariable fn){
		m_target     = target;
		m_startX1    = startX1;
		m_startX2    = startX2;
		m_tol        = tol;
		m_maxSteps   = maxSteps;
		m_fn         = fn;
	}
	////////////////////////////////////////////////////////
    public double solve(){
    	// The following method does a linear interpolation / extrapolation
    	// attempting to find the x that obtains  abs( target - fn(x)) < tol
    	
    	int     counter  = 0;
    	
    	double  x1 = m_startX1;
    	double  x2 = m_startX2;
    	
    	if(x1 == x2) {
    		// we'll move x2 away from zero.
    		x2 = x2 * 1.1 + 0.1 * Math.signum(x2) + 0.01;
    	} 

    	double  scaleFactor = 1.0;
    	if ( Math.abs(m_target) > 1)
    		scaleFactor = 1.0 / Math.abs(m_target);
    	
    	while ( true){
    		
    		counter++;

    		double  y1 = m_fn.fn(x1);
        	double  y2 = m_fn.fn(x2);
        	if ( Math.abs( y1 - m_target) * scaleFactor < m_tol){
        		System.out.println(  "After " + Integer.toString(counter) 
				           + " iterations, for x1 = " + Double.toString(x1) + " found y = " 
	                       + Double.toString(y1) + " to be within the tolerance of the target."); 
          	   return (x1); 
    		} else if ( counter >= m_maxSteps){
    			
    			System.out.println("After " + Integer.toString(counter) 
    					           + " iterations, with x value: " + Double.toString(x1) 
    					           + " had y value: "  + Double.toString(y1) 
    					           + " which was not within the tolerance.");
    			return x1;
    		} else if ( y1 == y2){
        		x2 = x2 * 1.1 + Math.signum(x2) * 0.1 + 0.01;
        	} else {
        		// If we use a linear approximation to interpolate or extrapolate, 
        		// either way, we have matched slopes.        	
        		// (m_target - y1)/(nextX - x1) = (y2 - y1)/(x2 - x1)
        		double nextX = x1 +(m_target - y1) * ( x2 - x1) / ( y2 - y1);
        	
        		double nextY = m_fn.fn(nextX);
 

	    		if (Math.abs(nextY - m_target) * scaleFactor < m_tol) {
	    			
	        		System.out.println(  "After " + Integer.toString(counter) 
	        				           + " iterations, for x = " + Double.toString(nextX) + " found y = " 
	        	                       + Double.toString(nextY) + " to be within the tolerance of the target."); 
	        		return nextX;
	        		
	        	} else if ( Math.abs(y1 - m_target ) < Math.abs(y2 - m_target)){  
	        	//} else if ( (counter % 2) == 0 ) {
	        		x1  = x1 + ( nextX - x1 ) * 0.5;
	        		x2 = nextX;
	        	} else {
	        		x1 = nextX;
	        		x2  = x2 + ( nextX - x2 ) * 0.5;
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
    	quadraticTest();
    	sinTest();
    	exponentialTest();
    }
    //////////////////////////////////////////////////////////////
    public static void sinTest(){
    	double target    = 1.0 / Math.sqrt(2.0);
    	double startX1   = 0.0;
    	double startX2   = 1.0;
    	double tol       = 1e-5;
    	int    maxSteps  = 50;
    	
    	IfaceFnOneVariable sinFn = (new IfaceFnOneVariable()
    	      {public double fn(double x) {return ( Math.sin(x));}  } );
    	
    	NumericalSolver ns = new NumericalSolver(target, startX1, startX2, tol, maxSteps, sinFn);
    	
    	double xSoln = ns.solve();
    	System.out.println("Found sin(" + Double.toString(xSoln) 
    			            + ") is close to target: " + Double.toString(target)); 
    	
    }
    //////////////////////////////////////////////////////////////
    public static void quadraticTest(){
    	double target    = 51.0;
    	double startX1   = 0.0;
    	double startX2   = 1.0;
    	double tol       = 1e-5;
    	int    maxSteps  = 50;
    	
    	IfaceFnOneVariable quadratic = (new IfaceFnOneVariable()
    	      {public double fn(double x) {return ( 2.0 * x * x + 4 * x + 3);}  } );
    	
    	NumericalSolver ns = new NumericalSolver(target, startX1, startX2, tol, maxSteps, quadratic);
    	
    	double xSoln = ns.solve();
    	System.out.println("Found quadratic(" + Double.toString(xSoln) 
    			            + ") is close to target: " + Double.toString(target)); 
    	    	
    }
    //////////////////////////////////////////////////////////////
    public static void exponentialTest(){
    	double target    = Math.exp(4.0);
    	double startX1   = -1.0;
    	double startX2   = 1.0;
    	double tol       = 1e-5;
    	int    maxSteps  = 50;
    	
    	IfaceFnOneVariable quadratic = (new IfaceFnOneVariable()
    	      {public double fn(double x) {return ( Math.exp(x));}  } );
    	
    	NumericalSolver ns = new NumericalSolver(target, startX1, startX2, tol, maxSteps, quadratic);
    	
    	double xSoln = ns.solve();
    	System.out.println("Found quadratic(" + Double.toString(xSoln) 
    			            + ") is close to target: " + Double.toString(target)); 
    	    	
    }
    //////////////////////////////////////////////////////////////
}

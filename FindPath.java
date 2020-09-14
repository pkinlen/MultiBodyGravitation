import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;

public class FindPath {

	private int m_meteoriteIdx;
    ////////////////////////////////////////////////////////////////////////////////
    public Config  getConfigData(){
    	
    	double  flyByTargetDistance = 2e8;      // the closest distance the meteorite gets to the surface.
    	
    	double  timeStepInterval    = 1 * 3600; // seconds, ( model seems to work when when this kept well below shortest_orbit_period / 100
                                                // For the moon, the orbit_period is just less than 30 days.
                                                // Could use 1 hour for m_timeStepInterval i.e. 3600 sec.
    	
    	int     numTimeSteps        = 100000;   // For 8 bodies, with 100,000 time steps, calc time is about 1 second.
    	
    	int     outputInterval      = 10;       // 1 means every time-step will be output, 2 means every second one, etc
                                                // The number of time-steps in the output file will be: 
                                                //      m_numTimeSteps / m_outputInterval
    	boolean displayFinalPos     = false;
    	boolean sendOutputToFile    = false;
    	  
    	String  outputFilePath      = new String ("outputs/movingBodies_" + getTimeStampString() + ".csv");
    	   	
        Config  config              = new Config( flyByTargetDistance, timeStepInterval, numTimeSteps,  outputInterval,
        		                                  displayFinalPos,     sendOutputToFile, outputFilePath                  );
        
        return config;
        		
   }
    /////////////////////////////////////////////////////////////
    public ArrayList<Body> getBodiesDataPolarCoords(double[] thetas){
        // Data source used:
        // https://nssdc.gsfc.nasa.gov/planetary/factsheet/
        // 
        // For motion in a circular orbit round the sun: speed = sqrt(GM/r)
        // Where M is the mass of the sun, 
   	    // and r is the distance between centre of body to centre of sun.
        //
        // For now we don't actually use the radii.  In this model we treat the bodies as point masses.
        
        ArrayList<Body> bodies = new ArrayList<Body>();
        
        double massOfSun = Config.m_massOfSun;
        
        Config config = getConfigData();
        // Here's we're starting with the sun in the centre and each planet starting at an angle theta to the pos0 axis.
        // Will calculate initial velocity assuming circular orbit. 
        
        int p = 0;
        
        //                  Name         Mass(kg,    BodyRadius(m) OrbitRadius(m),    Theta,       Config    )
        bodies.add(new Body("Sun",       massOfSun,  6.96340e8,    0.0,               0.0,         config    )); p++; 
        bodies.add(new Body("Mercury",   3.3e24,     2.44e6,       5.8e10,            thetas[p],   config    )); p++;
        bodies.add(new Body("Venus",     4.87e24,    6.05e6,       1.08e11,           thetas[p],   config    )); p++;
        bodies.add(new Body("Earth",     5.97e24,    6.37e6,       1.5e11,            thetas[p],   config    )); p++;       
        bodies.add(new Body("Mars",      6.42e23,    6.96e6,       2.28e11,           thetas[p],   config    )); p++;
        bodies.add(new Body("Jupiter",   1.89e27,    7.149e7,      7.79e11,           thetas[p],   config    )); p++;
        bodies.add(new Body("Saturn",    5.68e26,    6e7,          1.43e12,           thetas[p],   config    )); p++;
        bodies.add(new Body("Uranus",    8.68e25,    2.55e7,       2.87e12,           thetas[p],   config    )); p++; 
        bodies.add(new Body("Neptune",   1.02e24,    2.47e7,       4.5e12,            thetas[p],   config    )); p++;
        
        m_meteoriteIdx = p;
        bodies.add(new Body("Meteorite", 3e9,        3e2,          5e12,  0.0,  0.0,    -1e5, 0.0, 0.0   )); 

      // To do: put the moon back in        
      //bodies.add(new Body("Test",      10,         1,            1e11,    0.0,    0.0,        0.0,    36434.52363,  0.0    ));         
      
        return bodies;
   }          
	//////////////////////////////////////////////////////////
    public ArrayList<Body> getBodiesDataCartesianCoords(){
        // Data source used:
        // https://nssdc.gsfc.nasa.gov/planetary/factsheet/
        // 
        // For motion in a circular orbit round the sun: speed = sqrt(GM/r)
        // Where M is the mass of the sun, 
   	 // and r is the distance between centre of body to centre of sun.
        //
        // For now we don't actually use the radii.  In this model we treat the bodies as point masses.
        
        ArrayList<Body> bodies = new ArrayList<Body>();
        
        double massOfSun = Config.m_massOfSun;
        
        //                  Name         Mass(kg),   Radius(m)     Initial Pos(m) 0,1,2         Initial Vel (m/s) 0,1,2
        bodies.add(new Body("Sun",       massOfSun,  6.96340e8,    0.0,     0.0,    0.0,        0.0,    0.0,    0.0    )); 
        bodies.add(new Body("Mercury",   3.3e24,     2.44e6,       5.8e10,  0.0,    0.0,        0.0,    4.78e4, 0.0    )); 
        bodies.add(new Body("Venus",     4.87e24,    6.05e6,       1.08e11, 0.0,    0.0,        0.0,    3.5e4,  0.0    )); 
        bodies.add(new Body("Earth",     5.97e24,    6.37e6,       1.5e11,  0.0,    0.0,        0.0,    2.97e4, 0.0    )); 
        bodies.add(new Body("Moon",      7.3e22,     1.737e6,      1.5e11,  4e8,    0.0,        1e3,    2.97e4, 0.0    )); 
        bodies.add(new Body("Mars",      6.42e23,    6.96e6,       2.28e11, 0.0,    0.0,        0.0,    2.4e4,  0.0    )); 
        bodies.add(new Body("Jupiter",   1.89e27,    7.149e7,      7.79e11, 0.0,    0.0,        0.0,    1.3e4,  0.0    )); 
        bodies.add(new Body("Saturn",    5.68e26,    6e7,          1.43e12, 0.0,    0.0,        0.0,    9.63e3, 0.0    )); 
        bodies.add(new Body("Uranus",    8.68e25,    2.55e7,       2.87e12, 0.0,    0.0,        0.0,    6.8e3,  0.0    )); 
        bodies.add(new Body("Neptune",   1.02e24,    2.47e7,       4.5e12,  0.0,    0.0,        0.0,    5.43e3, 0.0    )); 
        bodies.add(new Body("Meteorite", 3e9,        3e2,          4.51e12, 1e10,   0.0,       -1e9,    0,      0.0    )); 

      //bodies.add(new Body("Test",      10,         1,            1e11,    0.0,    0.0,        0.0,    36434.52363,  0.0    ));         
      
        return bodies;
   }          
    ///////////////////////////////////////////////////////////
    public static void main(String[] args) {
    	
         System.out.println("Starting.");
         long startTime = System.nanoTime();

     	 FindPath fp = new FindPath();
    	 
    	 fp.findPositionsForFlyBy();
       //fp.calculateMovementOfBodies();

         double executionTime = 1e-9 * ( System.nanoTime() - startTime);
         System.out.println("Execution time is: " + Double.toString(executionTime) + " seconds.");          
         
         System.out.println("Completed.");
    }
    /* To-do:
     *     convert bodiesMoving calc into a function of one variable
     *     use a solver that will find the initial position that will obtain the desired fly-by distance
     *     
     */
    
    public void findPositionsForFlyBy(){

    	int firstPlanet = 8;
    	int lastPlanet  = 1;
    	
    	Config config = getConfigData();
    	
    	double thetas[] = getInitialPlanetAngles(Math.max(firstPlanet, lastPlanet));
    	reportFlyByDistances(thetas, config);

    	
    	for (int planet = firstPlanet; planet >= lastPlanet; planet--){
    		
    		NumericalSolver solver = getSolver(planet, thetas, config);
    		thetas[planet] = solver.solve();
    	}
    	
    	reportFlyByDistances(thetas, config);
    	
    }
    ////////////////////////////////////////////
	private void reportFlyByDistances(double[] thetas, Config config){
		
		for (int planet = 0; planet <= 8; planet++){
			
			ArrayList<Body> bodies = getBodiesDataPolarCoords(thetas);
	        BodiesMoving movingBodies = new BodiesMoving(bodies,config);
	        double closestDist = movingBodies.calcClosestDist(planet, m_meteoriteIdx);
			
	        Body p = bodies.get(planet);
	        
	        System.out.println("Found closest distance between the meteorite and " 
	                          + p.getName() + " is " + Double.toString(closestDist) 
	                          + " with theta: " + Double.toString(thetas[planet]));
		}
	}
    ////////////////////////////////////////////
	NumericalSolver  getSolver(final int planet, final double[] thetas, final Config config){
		
      	
      	NumericalSolver.IfaceFnOneVariable quadratic = (new NumericalSolver.IfaceFnOneVariable()
      	      {public double fn(double thetaOverRide) 
      	         {return (findClosestPath(config, thetas, planet, thetaOverRide ));}  
      	       } );
      	
      	double flyByTargetDist = config.getFlyByTargetDistance(); 
      	double startTheta1     =  0.0;
      	double startTheta2     =  Math.PI * 0.1;
      	double tol             =  1e-4;
      	int    maxSteps        =  10;
      	

      	NumericalSolver solver = new NumericalSolver(flyByTargetDist, startTheta1, startTheta2, 
      			                                     tol, maxSteps, quadratic);
      	return solver;
	}
    ////////////////////////////////////////////
    private double findClosestPath(final Config config, final double[] thetas, final int planet, double thetaOverRide ){
    	// possibly should make a deep copy of thetas
    	thetas[planet] = thetaOverRide;
    	ArrayList<Body> bodies = getBodiesDataPolarCoords(thetas);
        BodiesMoving movingBodies = new BodiesMoving(bodies,config);
        double closestDist = movingBodies.calcClosestDist(planet, m_meteoriteIdx);

        return closestDist;
    }
	////////////////////////////////////////////
	// could get the output of first solver as the starting point for the next...
    private double[]  getInitialPlanetAngles(int maxPlanet){
    	double[] thetas = new double[maxPlanet+1];
    	return thetas;
    }

    /////////////////////////////////////////////////
    public void calculateMovementOfBodies(){

         BodiesMoving movingBodies = new BodiesMoving(getBodiesDataCartesianCoords(), getConfigData());
         movingBodies.calculatePositions();
    }
   ////////////////////////////////////////////////////////////////////////////////
    public static String getTimeStampString(){
            DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy_MM_dd_HH_mm_ss");  
            LocalDateTime now = LocalDateTime.now();  
            return (dtf.format(now));                   
    }
    ////////////////////////////////////////////////////////////////////////////////    
}

/* Possible extensions and improvements:
bodies could spin
could check for collisions, when have big steps, this is trickier...
could check for overflow errors
could fail more gracefully if the output file cannot be written to.
look into possible division by zero errors, particularly if we don't have collision checks.
could set better initial positions and velocities, could check if distances are to the surface or between centres.
could remove the comma at the end of each line in the output file.

One bug / feature of the model is that if the centre of two bodies are on a direct collision course.
Then the objects can be predicted to accelerate massively, so much that they escape each others gravitational pull.
Ultimately it is caused by using the reciprocal of distance between centres, which can be close to zero.

When large time steps are used relative to the orbit period, the model performance deteriorates.    

*/




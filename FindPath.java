import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;

public class FindPath {

	private int m_meteoriteIdx;
    ////////////////////////////////////////////////////////////////////////////////
    public Config  getConfigData(boolean withOutput){
    	
    	double  flyByTargetDistance = 2e8;      // the closest distance the meteorite gets to the surface.
    	
    	double  timeStepInterval    = 1 * 30;   // seconds, ( model seems to work when when this kept well below shortest_orbit_period / 100
                                                // For the moon, the orbit_period is just less than 30 days.
                                                // Could use 1 hour for m_timeStepInterval i.e. 3600 sec.

    	double  totalTime           = 30 * 24 * 3600;   // total simulation time in seconds.
        int     frames              = 200;               // approximate number of times the positions are output. I.e., number of frames in final animation.
        int     numTimeSteps        = (int)(totalTime /timeStepInterval);        // For 8 bodes, with 100,000 time steps, calc time is about 1 second.
        int     outputInterval      = (int)(totalTime/(timeStepInterval * frames));          // 1 means every time-step will be output, 2 means every second one, etc    	    	
    	    	
    	boolean displayFinalPos     = withOutput;
    	boolean sendOutputToFile    = withOutput;
    	  
    	String  outputFilePath      = new String ("outputs/movingBodies_" + getTimeStampString() + ".csv");
    	   	
        Config  config              = new Config( flyByTargetDistance, timeStepInterval, numTimeSteps,  outputInterval,
        		                                  displayFinalPos,     sendOutputToFile, outputFilePath                  );
        
        return config;
        		
   }
    /////////////////////////////////////////////////////////////
    public ArrayList<Body> getBodiesDataPolarCoords(){
        // Data source used:
        // https://nssdc.gsfc.nasa.gov/planetary/factsheet/
        // 
        // For motion in a circular orbit round the sun: speed = sqrt(GM/r)
        // Where M is the mass of the sun, 
   	    // and r is the distance between centre of body to centre of sun.
        //
        // For now we don't actually use the radii.  In this model we treat the bodies as point masses.
        
        ArrayList<Body> bodies = new ArrayList<Body>();
        
        // Here's we're starting with the sun in the centre and each planet starting at an angle theta to the pos0 axis.
        // Will calculate initial velocity assuming circular orbit. 
        
        double earthOrbitRadius = 1.5e11;
        double earthTheta       = 2.652820628221811;
                
        //                  Name       Mass(kg,              BodyRadius(m) OrbitRadius(m),    Theta,    
        bodies.add(new Body("Sun",     Config.m_massOfSun,   6.96340e8,    0.0,                0.0   )); 
        bodies.add(new Body("Mercury", 3.3e24,               2.44e6,       5.8e10,             2.6026840338805863 )); 
        bodies.add(new Body("Venus",   4.87e24,              6.05e6,       1.08e11,            2.6402093069219816 )); 
        bodies.add(new Body("Earth",   Config.m_massOfEarth, 6.37e6,       earthOrbitRadius,   earthTheta         ));      
        bodies.add(new Body("Mars",    6.42e23,              6.96e6,       2.28e11,            2.6642716373928934 ));
        bodies.add(new Body("Jupiter", 1.89e27,              7.149e7,      7.79e11,            2.681278284481493  )); 
        bodies.add(new Body("Saturn",  5.68e26,              6e7,          1.43e12,            2.6849444632859543  )); 
        bodies.add(new Body("Uranus",  8.68e25,              2.55e7,       2.87e12,            2.6874496533791175 )); 
        bodies.add(new Body("Neptune", 1.02e24,              2.47e7,       4.5e12,             2.688502873555769  )); 

        addMoon(bodies, earthOrbitRadius, earthTheta);

        m_meteoriteIdx = 10;
        //                  Name         Mass(kg),   Radius(m)     Initial Pos(m) 0,1,2    Initial Vel (m/s) 0,1,2
        bodies.add(new Body("Meteorite", 3e9,        3e2,          8e9,  0,  0.0,          -2e6, 1e6, 0.0   )); 
        
      
        return bodies;
   }          
    //////////////////////////////////////////////////////////
    public void addMoon(ArrayList<Body> bodies, double earthOrbitRadius, double earthTheta  ){
    	
    	double moonTheta        = earthTheta - Math.PI * 0.5;
    	double moonOrbitRadius  = 4e8; // distance from earth to moon (m), should be centre to centre, rather than surface to surface
    	
        double moonPos0         = earthOrbitRadius * Math.cos(earthTheta) + moonOrbitRadius * Math.cos(moonTheta);	
        double moonPos1         = earthOrbitRadius * Math.sin(earthTheta) + moonOrbitRadius * Math.sin(moonTheta); 	
        double moonPos2         = 0.0;
        
   	    double sqrtGMoverREarth = Math.sqrt( Config.m_massOfSun   * Config.m_gravitationalConst / earthOrbitRadius);
   	    double sqrtGMoverRMoon  = Math.sqrt( Config.m_massOfEarth * Config.m_gravitationalConst / moonOrbitRadius );
        
        double moonVel0         = - Math.sin(earthTheta) * sqrtGMoverREarth - Math.sin(moonTheta) * sqrtGMoverRMoon;
        double moonVel1         =   Math.cos(earthTheta) * sqrtGMoverREarth + Math.cos(moonTheta) * sqrtGMoverRMoon;
        double moonVel2         = 0.0;
        
    	bodies.add(new Body("Moon",      
    			            7.3e22,                        // mass of moon   (kg)
    			            1.737e6,                       // radius of moon (m)
    			            moonPos0, moonPos1, moonPos2,
    			            moonVel0, moonVel1, moonVel2 ));
    			            

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

         // NumericalSolver.testNumericalSolver();
         FindPath fp = new FindPath();
    	 
          //fp.findPositionsForFlyBy();
          fp.calculateMovementOfBodies();

         double executionTime = 1e-9 * ( System.nanoTime() - startTime);
         System.out.println("Execution time is: " + Double.toString(executionTime) + " seconds.");          
         
         System.out.println("Completed.");
    }
    ///////////////////////////////////////////////
    public void findPositionsForFlyBy(){

    	int firstPlanet = 1;
    	int lastPlanet  = 8;
    	
    	Config config = getConfigData(false);
    	
    	double thetas[] = getInitialPlanetAngles(Math.max(firstPlanet, lastPlanet));
    	reportFlyByDistances(thetas, config);

    	for (int planet = firstPlanet; planet <= lastPlanet; planet++){
    		
    		NumericalSolver solver = getSolver(planet, thetas, config);
    		thetas[planet] = solver.solve();
    	}
    	
    	reportFlyByDistances(thetas, config);    	
    }
    ////////////////////////////////////////////
	private void reportFlyByDistances(double[] thetas, Config config){
		
		for (int planet = 0; planet <= 9; planet++){
			
			ArrayList<Body> bodies = getBodiesDataPolarCoords();
	        BodiesMoving movingBodies = new BodiesMoving(bodies,config);
	        double closestDist = movingBodies.calcClosestDist(planet, m_meteoriteIdx);
			
	        Body p = bodies.get(planet);
	        
	        if ( planet < thetas.length)
	        	System.out.println("Found closest distance between the meteorite and " 
	        						+ p.getName() + " is " + Double.toString(closestDist) 
	        						+ " with theta: " + Double.toString(thetas[planet]));
	        else
	        	System.out.println("Found closest distance between the meteorite and " 
						+ p.getName() + " is " + Double.toString(closestDist));
		}
	}
    ////////////////////////////////////////////
	NumericalSolver  getSolver(final int planet, final double[] thetas, final Config config){
		
      	
      	NumericalSolver.IfaceFnOneVariable fnOneVar = (new NumericalSolver.IfaceFnOneVariable()
      	      {public double fn(double thetaOverRide) 
      	         {return (findClosestPath(config, thetas, planet, thetaOverRide ));}  
      	       } );
      	
      	double flyByTargetDist = config.getFlyByTargetDistance(); 
      	double startTheta1     =  thetas[planet];
      	double startTheta2     =  thetas[planet] + 0.01;
      	double tol             =  1e-4;
      	int    maxSteps        =  40;
      	
      	NumericalSolver solver = new NumericalSolver(flyByTargetDist, startTheta1, startTheta2, 
      			                                     tol, maxSteps, fnOneVar);
      	return solver;
	}
    ////////////////////////////////////////////
    private double findClosestPath(final Config config, final double[] thetas, final int planet, double thetaOverRide ){

    	ArrayList<Body> bodies = getBodiesDataPolarCoords();
    	for(int p = 0; p < thetas.length; p++) {
    		if ( p == planet)
    			bodies.get(p).resetTheta(thetaOverRide);
    		else
    			bodies.get(p).resetTheta(thetas[p]);    	
    	}
    	
        BodiesMoving movingBodies = new BodiesMoving(bodies,config);
        double closestDist = movingBodies.calcClosestDist(planet, m_meteoriteIdx);

        return closestDist;
    }
	////////////////////////////////////////////
	// could get the output of first solver as the starting point for the next...
    private double[]  getInitialPlanetAngles(int maxPlanet){
    	double[] thetas = new double[maxPlanet+1];
    	
    	ArrayList<Body> bodies = getBodiesDataPolarCoords();
    	for(int p = 0; p < thetas.length; p++)
    		thetas[p] = bodies.get(p).getTheta();
    	
    	return thetas;
    }

    /////////////////////////////////////////////////
    public void calculateMovementOfBodies(){

         BodiesMoving movingBodies = new BodiesMoving(getBodiesDataPolarCoords(), getConfigData(true));
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


import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.time.format.DateTimeFormatter;  
import java.time.LocalDateTime;    
import java.util.ArrayList;
import java.util.List;
import java.lang.Math;

// Extensions:
// bodies could spin
// could check for collisions

// This class calculates how bodies move under Newtonian Gravity
public class BodiesMoving {
	private static final double m_gravitationalCons = 0;
	String m_outputFilePath;
	double m_gravitationalConst;
	double m_timeStepInterval; // seconds
	long   m_numTimeSteps;
	long   m_outputInterval;  // 1 means every time-step will be output, 2 means every second one, etc
	Writer  m_fileForOutput;
	ArrayList<Body> m_bodies;
	
	///////////////////////////////////////////////////////////
	public static void main(String[] args) {
		System.out.println("Just starting.");
		
		BodiesMoving movingBodies = new BodiesMoving();
		movingBodies.initialize();
		movingBodies.calculatePositions();
		movingBodies.closeOutputFile();
		System.out.println("Completed.");
	}
	//////////////////////////////////////////////////////////
	public void initialize(){
		setConfigData();
		setBodiesData();
		openFileForOutput();
	}
	//////////////////////////////////////////////////////////
    public void setBodiesData(){
    	// Data source used:
    	// https://nssdc.gsfc.nasa.gov/planetary/factsheet/
    	// Though should used distances to centres and 
    	// could set realistic initial positions and velocities
    	// For motion in a circle: speed = sqrt(GM/r)
    	// Where M is the mass of the sun.
    	
    	m_bodies = new ArrayList<Body>();
    	// -----------------------------------
    	Body sun = new Body("Sun",           // name
    			            1.989e30,        // mass   (kg)
    			            6.96340e8,       // radius (meters)    
    	                    0.0, 0.0, 0.0,   // initial pos (meters)
    	                    0.0, 0.0, 0.0);  // initial vel (meters / sec)
    	m_bodies.add(sun);
    	// -----------------------------------
    	Body mercury = new Body("Mercury",           // name
	            3.3e24,        // mass   (kg)
	            2.44e6,        // radius (meters)    
                5.8e10, 0.0, 0.0,   // initial pos (meters)
                0.0, 4.78e4, 0.0);  // initial vel (meters / sec)
        m_bodies.add(mercury);
    	// -----------------------------------
    	Body venus = new Body("Venus",           // name
	            4.87e24,        // mass   (kg)
	            6.05e6,       // radius (meters)    
                1.08e11, 0.0, 0.0,   // initial pos (meters)
                0.0, 3.5e4, 0.0);  // initial vel (meters / sec)
        m_bodies.add(venus);
    	// -----------------------------------
    	Body earth = new Body("Earth",           // name
	            5.97e24,        // mass   (kg)
	            6.37e6,       // radius (meters)    
                1.5e11, 0.0, 0.0,   // initial pos (meters)
                0.0, 2.97, 0.0);  // initial vel (meters / sec)
        m_bodies.add(earth);
    	// -----------------------------------
    	Body moon = new Body("Moon",  // name
	            7.3e22,               // mass   (kg)
	            1.737e6,              // radius (meters)    
                1.5e11, 3.8e5, 0.0,   // initial pos (meters)
                3.24e4, 2.97, 0.0);       // initial vel (meters / sec)
        m_bodies.add(moon);
    	// -----------------------------------
    	Body mars = new Body("Mars", // name
	            1.989e30,            // mass   (kg)
	            6.96340e8,           // radius (meters)    
                2.28e11, 0.0, 0.0,   // initial pos (meters)
                0.0, 2.4e4, 0.0);      // initial vel (meters / sec)
        m_bodies.add(mars);
    	// -----------------------------------
    	Body jupiter = new Body("Jupiter",           // name
	            1.898e27,            // mass   (kg)
	            7.149e7,             // radius (meters)    
                7.79e11, 0.0, 0.0,   // initial pos (meters)
                0.0, 1.3e4, 0.0);      // initial vel (meters / sec)
        m_bodies.add(jupiter);
    	// -----------------------------------
    	Body meteorite = new Body("Meteorite",           // name
	            3e9,                  // mass   (kg)
	            300,                  // radius (meters)    
	            1.5e11, 1.9e5, 0.0,   // initial pos (meters)
                1e6, 1e4, 1e2);       // initial vel (meters / sec)
        m_bodies.add(meteorite);
        // -----------------------------------
        
        
    }	
    //////////////////////////////////////////////////////////
    public double getDistance(int i, int j, double[] currentPos){
    	
    	double d0 = currentPos[3 * i + 0] - currentPos[3 * j + 0];
    	double d1 = currentPos[3 * i + 1] - currentPos[3 * j + 1];
    	double d2 = currentPos[3 * i + 2] - currentPos[3 * j + 2];
    	
    	return Math.sqrt((d0 * d0) + (d1 * d1) + (d2 * d2));
    	
    }
    //////////////////////////////////////////////////////////
    //                                     input                   outupt
    public void updateMatOfForces(double[] currentPos, double[] matForces){
    
	   // We have a matrix of 1/ dist^2, where dist is the distance between the bodies 
	   double[] matRecipSqDis = new double[ m_bodies.size() * m_bodies.size()];
	   // could check for collisions. compare with matrix of 1 / ( radius_i + radius_j)^2
	   
	   int numBodies = m_bodies.size();
	   
	   for( int i= 0; i < numBodies; i++ ){
		   for (int j = 0; j < numBodies; j++){
			   if ( j < i) {
				   double Gmm      = m_gravitationalCons * m_bodies.get(i).getMass() * m_bodies.get(j).getMass();
				   double dist     = getDistance(i,j, currentPos);
				   
				   // in the following line we divide by distiance to the power of 3.
				   // power of 2 for the gravitational force and one more to normalized vector
				   double GmmOverD3 = Gmm / ( dist * dist * dist);
				   double posDiff0 = currentPos[3*i + 0] - currentPos[3*j + 0];
				   matForces[3*(i * numBodies + j) + 0] = posDiff0 * GmmOverD3;  

				   double posDiff1 = currentPos[3*i + 1] - currentPos[3*j + 1];
				   matForces[3*(i * numBodies + j) + 1] = posDiff1 * GmmOverD3;  

				   double posDiff2 = currentPos[3*i + 2] - currentPos[3*j + 2];
				   matForces[3*(i * numBodies + j) + 2] = posDiff2 * GmmOverD3;  
				   
				   
			   } else if ( j > i){
				   matForces[3*(i * numBodies + j) + 0] = - matForces[3*(j * numBodies + i) + 0];  
				   matForces[3*(i * numBodies + j) + 1] = - matForces[3*(j * numBodies + i) + 1];  
				   matForces[3*(i * numBodies + j) + 2] = - matForces[3*(j * numBodies + i) + 2];  
			   }
			   // else if (i == j) 
			   // then leave the forces at zero
		   }
		   
	   }
    }
    //////////////////////////////////////////////////////////
    public void calculatePositions(){
    	
    	// size of array will be 3 * numBodies. (3 since we have 3 spacial dimensions)
    	double[] currentPos = getInitialPos();
    	storeCurrentPos(0, currentPos);

    	double[] currentVelocity = getInitialVel();
    	
    	double[] matForces     = new double[ m_bodies.size() * m_bodies.size()];
    	    	
    	for( int t = 1; t <= m_numTimeSteps; t++){
    		String posData = new String("");
    		     
    		//                input          output
    	    updateMatOfForces(currentPos, matForces);

    	for( int b = 0; b < m_bodies.size(); b++){	
    	    double[] acceleration = getAcceleration(b, matForces);
    		next_vel = previous_vel + (acceleration * m_timeStepInterval);
    		next_pos = previous_pos + ( previous_vel + next_vel ) * m_timeStepInterval * 0.5;  
    		
    	  } 
    	
    	  storeCurrentPos(t, currentPos);
       }    // next time_step
    }
    //////////////////////////////////////////////////////////    
    public void storeCurrentPos(int t, double[] currentPos){
    	
      if ( (t % m_outputInterval) == 0) {	
    	
        String posDataStr = new String("");
        
        for( int i = 0; i < m_bodies.size(); i++){
        	
        	posDataStr +=  Double.toString(currentPos[i + 0]) + "," 
        			     + Double.toString(currentPos[i + 1]) + ","
                         + Double.toString(currentPos[i + 2]) + ",";
        } // possibly could remove the last ","
        
    	  
	    try {
  	      	 m_fileForOutput.write(posDataStr); // may need to add new line
	  	 } catch (Exception e) {
  		   System.out.println("When trying to write to file got exception: " + e.getMessage());
	  	 }
      }
    }    
    //////////////////////////////////////////////////////////    
    public void openFileForOutput(){
    	try{
    		m_fileForOutput = new PrintWriter(m_outputFilePath, "UTF-8");
    		//m_fileForOutput.println("The first line");
    	} catch(Exception e){
    		System.out.println("When trying to open file for output got exception: " + e.getMessage());
    	}
    }
    ////////////////////////////////////////////////////////////
    public void closeOutputFile(){
    	if ( m_fileForOutput == null)
    		System.out.println("Found that the output file was null.");
    	else {
    		try {m_fileForOutput.close(); }
    		catch(Exception e){ 
    			 System.out.println("When trying to close output file, got exception: " + e.getMessage());
    		}
    	}
    }
	/////////////////////////////////////////////////////////////////////////////
	
    public void setConfigData(){    	//////////////////////////////////////////////////////////
    	
    	m_outputFilePath     = new String ("movingBodies_" + getTimeStampString() + ".csv");
    	m_gravitationalConst = 6.67408e-11; // m^3 kg^-1 s^-2
    	m_timeStepInterval   = 24 * 3600; // seconds
    	m_numTimeSteps       = 1000;
    	m_outputInterval     = 1;  // 1 means every time-step will be output, 2 means every second one, etc
    	
    }
	//////////////////////////////////////////////////////////
    public static String getTimeStampString(){
    	   DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy_MM_dd_HH_mm_ss");  
    	   LocalDateTime now = LocalDateTime.now();  
    	   return (dtf.format(now));     	    	
    }
	//////////////////////////////////////////////////////////
    
    public class Body{
    	String m_name;
    	double m_mass;      // kg
    	double m_radius;    // meters, could be used for checking for collisions
    	double m_initPos1;  // initial position on axis 1 (meters)
    	double m_initPos2;  // initial position on axis 2
    	double m_initPos3;  // initial position on axis 3
    	double m_initVel1;  // initial velocity along axis 1 ( meters/ second)
    	double m_initVel2;
    	double m_initVel3;
    	
    	public Body(String name, double mass, double radius, 
    		        double initPos1, double initPos2, double initPos3,
    		        double initVel1, double initVel2, double initVel3){
    		
    		m_name      = name;
    		m_mass      = mass;
    		m_radius    = radius;
    		m_initPos1  = initPos1;
    		m_initPos2  = initPos2;
    		m_initPos3  = initPos3;
    		m_initVel1  = initVel1;
    		m_initVel2  = initVel2;
    		m_initVel3  = initVel3;
    		
    	}
    	
    	public String   getName()     { return m_name;     }
    	public double   getMass()     { return m_mass;     }
    	public double   getRadius()   { return m_radius;   }
    	public double   getInitPos1() { return m_initPos1; }
    	public double   getInitPos2() { return m_initPos2; }
    	public double   getInitPos3() { return m_initPos3; }
    	public double   getInitVel1() { return m_initVel1; }
    	public double   getInitVel2() { return m_initVel2; }
    	public double   getInitVel3() { return m_initVel3; }
    }
}

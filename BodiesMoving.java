
import java.io.PrintWriter;
import java.io.Writer;
import java.time.format.DateTimeFormatter;  
import java.time.LocalDateTime;    
import java.util.ArrayList;
import java.lang.Math;


// This class calculates how bodies move under Newtonian Gravity

// Notes:
// To get a good estimate of the final velocity, could use 
//    velocity = ( X[T] - X[T-1]) / (m_timeStepsInterval * m_outputInterval)

public class BodiesMoving {
	private String          m_outputFilePath;
	private double          m_timeStepInterval; // seconds
	private long            m_numTimeSteps;
	private long            m_outputInterval;  // 1 means every time-step will be output, 2 means every second one, etc
	private double          m_gravitationalConst;
	private Writer          m_fileForOutput;
	private ArrayList<Body> m_bodies;
	private boolean         m_displayFinalPos;
	private String          m_lineSep;
	/////////////////////////////////////////////////////////////////////////////
    public void setConfigData(){    	
    	
    	m_outputFilePath     = new String ("outputs/movingBodies_" + getTimeStampString() + ".csv");
    	m_timeStepInterval   = 1 * 24 * 3600;   // seconds
    	m_numTimeSteps       = 50;
    	m_outputInterval     = 1;           // 1 means every time-step will be output, 2 means every second one, etc
    	m_displayFinalPos    = true;
    	m_gravitationalConst = 6.67408e-11; // m^3 kg^-1 s^-2

    	
    	m_lineSep            = System.getProperty( "line.separator" );
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
    	
    	//                    Name         Mass(kg),   Radius(m)     Initial Pos(m) 0,1,2         Initial Vel (m/s) 0,1,2
    	m_bodies.add(new Body("Sun",       1.989e30,   6.96340e8,    0.0,     0.0,    0.0,        0.0,    0.0,    0.0    )); 
    	m_bodies.add(new Body("Mercury",   3.3e24,     2.44e6,       5.8e10,  0.0,    0.0,        0.0,    4.78e4, 0.0    )); 
    	m_bodies.add(new Body("Venus",     4.87e24,    6.05e6,       1.08e11, 0.0,    0.0,        0.0,    3.5e4,  0.0    )); 
    	m_bodies.add(new Body("Earth",     5.97e24,    6.37e6,       1.5e11,  0.0,    0.0,        0.0,    2.97e4, 0.0    )); 
    	m_bodies.add(new Body("Moon",      7.3e22,     1.737e6,      1.5e11,  3.8e5,  0.0,        3.24e4, 2.97e4, 0.0    )); 
    	m_bodies.add(new Body("Mars",      6.42e23,    6.96e6,       2.28e11, 0.0,    0.0,        0.0,    2.4e4,  0.0    )); 
    	m_bodies.add(new Body("Jupiter",   1.89e27,    7.149e7,      7.79e11, 0.0,    0.0,        0.0,    1.3e4,  0.0    )); 
    	m_bodies.add(new Body("Meteorite", 3e9,        3e2,          1.5e11,  1.9e5,  0.0,        -1e6,   -1e4,   -1e2    )); 

 //   	m_bodies.add(new Body("Test"",      10,         1,            1.5e11,     0.0,    0.0,     0.0,    3.5e4,    0.0    )); 
   	
    }		
	///////////////////////////////////////////////////////////
	public static void main(String[] args) {
		System.out.println("Starting.");
		
		BodiesMoving movingBodies = new BodiesMoving();
		movingBodies.initialize();
		movingBodies.calculatePositions();
		movingBodies.finalizeOutput();
		
		System.out.println("Completed.");
	}
	//////////////////////////////////////////////////////////
	public void initialize(){
		setConfigData();
		setBodiesData();
		openFileForOutput();
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
    	   
	   int numBodies = m_bodies.size();
	   
	   for( int i= 0; i < numBodies; i++ ){
		   for (int j = 0; j < numBodies; j++){
			   if ( i < j) {
				   double Gmm      = m_gravitationalConst * m_bodies.get(i).getMass() * m_bodies.get(j).getMass();
				   double dist     = getDistance(i,j, currentPos);
				   
				   // in the following line we divide by distiance to the power of 3.
				   // power of 2 for the gravitational force and one more to normalized vector
				   double GmmOverD3 = Gmm / ( dist * dist * dist);
				   double posDiff0 = currentPos[3*i + 0] - currentPos[3*j + 0];
				   matForces[3*(i * numBodies + j) + 0] = -posDiff0 * GmmOverD3;  

				   double posDiff1 = currentPos[3*i + 1] - currentPos[3*j + 1];
				   matForces[3*(i * numBodies + j) + 1] = -posDiff1 * GmmOverD3;  

				   double posDiff2 = currentPos[3*i + 2] - currentPos[3*j + 2];
				   matForces[3*(i * numBodies + j) + 2] = -posDiff2 * GmmOverD3;  
				   
				   
			   } else if ( i > j){
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
    	double[] currentPos = getInitPos();
    	storeCurrentPos(0, currentPos);

    	double[] currentVel = getInitVel();
    	
    	double[] matForces    = new double[ 3 * m_bodies.size() * m_bodies.size()];
    	for( int t = 1; t <= m_numTimeSteps; t++){
    		     
    		//                input       output
    	    updateMatOfForces(currentPos, matForces);

    	for( int b = 0; b < m_bodies.size(); b++){	
    	    double[] acceleration = getAcceleration(b, matForces);
    	    
    	    double nextVel_0 = currentVel[3 * b + 0] + acceleration[0] * m_timeStepInterval;
    	    double nextVel_1 = currentVel[3 * b + 1] + acceleration[1] * m_timeStepInterval;
    	    double nextVel_2 = currentVel[3 * b + 2] + acceleration[2] * m_timeStepInterval;
    	    
    	    currentPos[3 * b + 0] += (nextVel_0 + currentVel[3 * b + 0]) * m_timeStepInterval * 0.5;
    	    currentPos[3 * b + 1] += (nextVel_1 + currentVel[3 * b + 1]) * m_timeStepInterval * 0.5;
    	    currentPos[3 * b + 2] += (nextVel_2 + currentVel[3 * b + 2]) * m_timeStepInterval * 0.5;

    	    currentVel[3 * b + 0] = nextVel_0;
    	    currentVel[3 * b + 1] = nextVel_1;
    	    currentVel[3 * b + 2] = nextVel_2;
    	  } 
    	
    	  storeCurrentPos(t, currentPos);
       }    // next time_step
       if (m_displayFinalPos)
    	   displayPos(currentPos);
    }
    //////////////////////////////////////////////////////////    
    public void displayPos(double[] pos){
    	 System.out.println("Final Positions:");
         for(int b = 0; b < m_bodies.size(); b++){
        	 
            System.out.println("    " + m_bodies.get(b).getName() + ": " 
                               + Double.toString(pos[3 * b + 0]) + ", "
                               + Double.toString(pos[3 * b + 1]) + ", "
                               + Double.toString(pos[3 * b + 2]));                             
         }	
    }
    //////////////////////////////////////////////////////////    
    public double[] getAcceleration(int b, double[] matForces){
    	double[] acceleration = new double[3];
    	int numBodies = m_bodies.size();
    	double oneOverMass = 1.0 / m_bodies.get(b).getMass();
    	for( int B = 0; B < numBodies; B++){
    		acceleration[0] += matForces[ 3*(b * numBodies + B) + 0] * oneOverMass;
    		acceleration[1] += matForces[ 3*(b * numBodies + B) + 1] * oneOverMass;
    		acceleration[2] += matForces[ 3*(b * numBodies + B) + 2] * oneOverMass;
    	}
    	return acceleration;
    }
    //////////////////////////////////////////////////////////    
    public double[] getInitPos(){
    	double[] initPos = new double[m_bodies.size() * 3];
    	
    	for(int b = 0; b < m_bodies.size(); b++) {// could use iterator
    		Body body = m_bodies.get(b);
    		initPos[3*b + 0] = body.getInitPos0();
    		initPos[3*b + 1] = body.getInitPos1();
    		initPos[3*b + 2] = body.getInitPos2();		
    	}  
    	return initPos;    	
    }
    //////////////////////////////////////////////////////////    
    public double[] getInitVel(){
    	double[] initVel = new double[m_bodies.size() * 3];
    	
    	for(int b = 0; b < m_bodies.size(); b++) {// could use iterator
    		Body body = m_bodies.get(b);
    		initVel[3*b + 0] = body.getInitVel0();
    		initVel[3*b + 1] = body.getInitVel1();
    		initVel[3*b + 2] = body.getInitVel2();
    		
    	}  
    	return initVel;    	
    }
    ///////////////////////////////////////////////////////////
    public void storeCurrentPos(int t, double[] currentPos){
    	
      if ( (t % m_outputInterval) == 0) {	
    	
        String posDataStr = new String("");
        
        for( int i = 0; i < m_bodies.size(); i++){
        	
        	posDataStr +=  Double.toString(currentPos[3*i + 0]) + "," 
        			     + Double.toString(currentPos[3*i + 1]) + ","
                         + Double.toString(currentPos[3*i + 2]) + ",";
        } // possibly could remove the last ","
        
    	  
	    try {
  	      	 m_fileForOutput.write(posDataStr + m_lineSep); 
	  	 } catch (Exception e) {
  		   System.out.println("When trying to write to file got exception: " + e.getMessage());
	  	 }
      }
    }    
    //////////////////////////////////////////////////////////    
    public void openFileForOutput(){
    	try{
    		m_fileForOutput = new PrintWriter(m_outputFilePath, "UTF-8");
    		System.out.println("Sending output to file: " + m_outputFilePath);
    		
    		String firstLine = new String("");
    		for(int b = 0; b < m_bodies.size(); b++){
    			Body body = m_bodies.get(b);
    			firstLine += body.getName() + ":" + Double.toString(body.getRadius()) +", ";
    		}
    		m_fileForOutput.write(firstLine + m_lineSep);
    		
    		//m_fileForOutput.println("The first line");
    	} catch(Exception e){
    		System.out.println("When trying to open file for output got exception: " + e.getMessage());
    	}
    }
    ////////////////////////////////////////////////////////////
    public void finalizeOutput(){
    	if ( m_fileForOutput == null)
    		System.out.println("Found that the output file was null.");
    	else {
    		try {m_fileForOutput.close(); }
    		catch(Exception e){ 
    			 System.out.println("When trying to close output file, got exception: " + e.getMessage());
    		}
    	}
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
    	
    	double m_initPos0;  // initial position on axis 0 (meters)
    	double m_initPos1;  // initial position on axis 1
    	double m_initPos2;  // initial position on axis 2
    	
    	double m_initVel0;  // initial velocity along axis 1 ( meters/ second)
    	double m_initVel1;
    	double m_initVel2;
    	//////////////////////////////////////////////////////////////
    	public Body(String name, double mass, double radius, 
    		        double initPos0, double initPos1, double initPos2,
    		        double initVel0, double initVel1, double initVel2){
    		
    		m_name      = name;
    		m_mass      = mass;
    		m_radius    = radius;
    		
    		m_initPos0  = initPos0;
    		m_initPos1  = initPos1;
    		m_initPos2  = initPos2;
    		
    		m_initVel0  = initVel0;
    		m_initVel1  = initVel1;
    		m_initVel2  = initVel2;    		
    	}
    	//////////////////////////////////////////////////////////////
    	public String   getName()     { return m_name;     }
    	public double   getMass()     { return m_mass;     }
    	public double   getRadius()   { return m_radius;   }
    	
    	public double   getInitPos0() { return m_initPos0; }
    	public double   getInitPos1() { return m_initPos1; }
    	public double   getInitPos2() { return m_initPos2; }
    	
    	public double   getInitVel0() { return m_initVel0; }
    	public double   getInitVel1() { return m_initVel1; }
    	public double   getInitVel2() { return m_initVel2; }
    	//////////////////////////////////////////////////////////////
    }
}

/* Possible extensions and improvements:
     bodies could spin
     could check for collisions, when have big steps, this is trickier...
     could check for overflow errors
     could fail more gracefully if the output file cannot be written to.
     look into possible division by zero errors, particularly if we don't have collision checks.
     
 */



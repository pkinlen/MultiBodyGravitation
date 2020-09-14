
import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.lang.Math;

// This class calculates how multiple bodies move, each being pulled by the gravitational force of the others.
// It uses a Newtonian model and does not use Einstein's relativity.
// It will write output to a file. The first line will contain the names of the bodies and their radii.
// The rest of the rows will contain their positions in 3 dimensional cartesian co-ordinates.
// The unit of distance used is meter.

// Notes:
// To get a good estimate of the final velocity from the output file, could use 
//    velocity = ( X[T] - X[T-1]) / (m_timeStepsInterval * m_outputInterval)

public class BodiesMoving {
     private String          m_outputFilePath;
     private double          m_timeStepInterval; // seconds
     private long            m_numTimeSteps;
     private long            m_outputInterval;  // 1 means every time-step will be output, 2 means every second one, etc
     private Writer          m_fileForOutput;
     private ArrayList<Body> m_bodies;
     private boolean         m_displayFinalPos;
     private String          m_lineSep;
     private boolean         m_alreadyReportedNullFile;
     private boolean         m_sendOutputToFile;
     
     private boolean         m_checkFlyBy;
     private double          m_minDistSq;    // the minimum distance squared so far recorded between A and B
     private int             m_flyByBodyA;
     private int             m_flyByBodyB;
     private double[]        m_previousAPos;
     private double[]        m_previousBPos;
     

     private double          m_gravitationalConst; 
     
    /////////////////////////////////////////////////////////////////////////////
     public BodiesMoving(ArrayList<Body> bodies, Config config){
          
          m_outputFilePath     = config.getOutputFilePath();
          m_timeStepInterval   = config.getTimeStepInterval();
          m_numTimeSteps       = config.getNumTimeSteps();
          m_outputInterval     = config.getOutputInterval();
          m_displayFinalPos    = config.displayFinalPos();
          m_lineSep            = config.getLineSep();
          m_sendOutputToFile   = config.sendOutputToFile();
          m_gravitationalConst = Config.m_gravitationalConst;
          
          m_bodies = bodies;
          
          if ( m_sendOutputToFile) 
        	  openFileForOutput();
          
          m_checkFlyBy         = false;    
          m_minDistSq          = -1.0;  // this indicates that the min dist squared has not yet been set.
     }
     //////////////////////////////////////////////////////////
     // The distance squared
     public double getDistanceSq(int threeA, int threeB, double[] currentPos){
         
         double d0 = currentPos[threeA + 0] - currentPos[threeB + 0];
         double d1 = currentPos[threeA + 1] - currentPos[threeB + 1];
         double d2 = currentPos[threeA + 2] - currentPos[threeB + 2];

         return ((d0 * d0) + (d1 * d1) + (d2 * d2));
     }
     /////////////////////////////////////////////////////////
     public double getDistance(int threeA, int threeB, double[] currentPos){
         return Math.sqrt(getDistanceSq(threeA, threeB, currentPos));         
    }
    ///////////////////////////////////////////////////////////
    public double calcClosestDist(int bodyA, int bodyB){
        m_checkFlyBy = true;
        m_flyByBodyA = bodyA;
        m_flyByBodyB = bodyB;
        
        calculatePositions();
        
        double minDist = Math.sqrt(m_minDistSq);
        
        /*
        System.out.println(   "Found the closest distance between "  
                            + m_bodies.get(bodyA).getName() + " and " 
		                    + m_bodies.get(bodyB).getName() 
		                    + " to be : " + Double.toString(minDist));        
         */
        return minDist;
    }     
    //////////////////////////////////////////////////////////
    //                                     input                outupt
    public void updateMatOfForces(double[] currentPos, double[] matForces){
            
        int numBodies = m_bodies.size();
        
        for( int a= 0; a < numBodies; a++ ){
        	 int threeA = 3 * a; // we want to multiply by 3 once rather than numerous times.
        	 int threeANumBodies = threeA * numBodies;
        	         	 
             for (int b = 0; b < numBodies; b++){
            	 int threeB          = 3 * b;
            	 int threeBNumBodies = threeB * numBodies;
                  if ( a < b) {
                       double Gmm      = m_gravitationalConst * m_bodies.get(a).getMass() * m_bodies.get(b).getMass();
                       double dist     = getDistance(threeA, threeB, currentPos);
                       
                       // In the following line we divide by distance to the power of 3.
                       // power of 2 for the gravitational force and one more to normalized the posDiff vector
                       double GmmOverD3 = Gmm / ( dist * dist * dist);
                       
                       double posDiff0 = currentPos[threeB + 0] - currentPos[threeA + 0];
                       double posDiff1 = currentPos[threeB + 1] - currentPos[threeA + 1];
                       double posDiff2 = currentPos[threeB + 2] - currentPos[threeA + 2];

                       matForces[threeANumBodies + threeB + 0] = posDiff0 * GmmOverD3;  
                       matForces[threeANumBodies + threeB + 1] = posDiff1 * GmmOverD3;  
                       matForces[threeANumBodies + threeB + 2] = posDiff2 * GmmOverD3;                         
                       
                  } else if ( a > b){
                       matForces[threeANumBodies + threeB + 0] = - matForces[threeBNumBodies + threeA + 0];  
                       matForces[threeANumBodies + threeB + 1] = - matForces[threeBNumBodies + threeA + 1];  
                       matForces[threeANumBodies + threeB + 2] = - matForces[threeBNumBodies + threeA + 2];  
                  }
                  // else if (a == b) 
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
         double   halfTimeStep = m_timeStepInterval * 0.5;
         
         double[] matForces    = new double[ 3 * m_bodies.size() * m_bodies.size()];
         for( int t = 1; t <= m_numTimeSteps; t++){
                   
             //                input       output
             updateMatOfForces(currentPos, matForces);

             for( int b = 0; b < m_bodies.size(); b++){    
            	int      threeB         = 3 * b;
                double[] acceleration   = getAcceleration(b, matForces);
             
                double nextVel_0        = currentVel[threeB + 0] + acceleration[0] * m_timeStepInterval;
                double nextVel_1        = currentVel[threeB + 1] + acceleration[1] * m_timeStepInterval;
                double nextVel_2        = currentVel[threeB + 2] + acceleration[2] * m_timeStepInterval;
             
                currentPos[threeB + 0] += (nextVel_0 + currentVel[threeB + 0]) * halfTimeStep;
                currentPos[threeB + 1] += (nextVel_1 + currentVel[threeB + 1]) * halfTimeStep;
                currentPos[threeB + 2] += (nextVel_2 + currentVel[threeB + 2]) * halfTimeStep;

                currentVel[threeB + 0]  = nextVel_0;
                currentVel[threeB + 1]  = nextVel_1;
                currentVel[threeB + 2]  = nextVel_2;
             } // next body b 
             
             if( m_checkFlyBy)
                 checkFlyBy(currentPos);
         
             storeCurrentPos(t, currentPos);
       }    // next time_step t
         
       if (m_displayFinalPos)
            displayPos(currentPos);
       
       finalizeOutput();
    }
    //////////////////////////////////////////////////////////
    public void checkFlyBy(double[] pos){
    	
    	double minDistSqOnThisTimeStep = getMinDistSqOnThisTimeStep(pos);
    	
    	if ( (m_minDistSq < 0) || (m_minDistSq > minDistSqOnThisTimeStep) ) 
    		 m_minDistSq = minDistSqOnThisTimeStep;
    }
    //////////////////////////////////////////////////////////
    public double getMinDistSqOnThisTimeStep(double[] pos){
        // rather that just look at the distance at the current position,
    	// we consider over the entire time step.
    	// for simplicity we assume that there is no acceleration.
    	// If we have a(t) = a(0)+ t * u
    	//            b(t) = b(0)+ t * v
    	// where a,b,u,v are all vectors in 3 dimension space.    	
    	// The position of a relative to b will be:
    	// p(t) = a(t) - b(t) = (a(0)-b(0)) + t * (u - v)
    	// let w = u - v
    	// a will be closest to b when p(t) is closest to the origin.
    	// at that time s of shortest distance: p(s) . w = 0
    	// i.e. the cross product is zero because p(s) is perpendicular to w
    	// we introduce a scalar variable alpha(t):
    	// p(t) = p(s) + alpha(t) * w
    	// p(t) . w = alpha(t) * w.w
    	// 
    	// If alpha has changed sign between the beginning and end of the time step,
    	// then the minimum is at p(s)
    	// i.e. sign(p(0) . w) != sign( p(T).w)
    	// p(0) = a(0)-b(0)
    	// p(T) = a(T)-b(T)
    	// w = u-v = a(T) -a(0) + b(0) - b(T)
    	    	
    	if(m_minDistSq < 0)
    		setPreviousPos();
    	
    	double p0_0 = m_previousAPos[0] - m_previousBPos[0]; 
    	double p0_1 = m_previousAPos[1] - m_previousBPos[1]; 
    	double p0_2 = m_previousAPos[2] - m_previousBPos[2]; 
    	
    	double pT_0 = pos[3 * m_flyByBodyA + 0] - pos[3 * m_flyByBodyB + 0]; 
    	double pT_1 = pos[3 * m_flyByBodyA + 1] - pos[3 * m_flyByBodyB + 1]; 
    	double pT_2 = pos[3 * m_flyByBodyA + 2] - pos[3 * m_flyByBodyB + 2];
    	
    	double w0   = pT_0 - p0_0;
    	double w1   = pT_1 - p0_1;
    	double w2   = pT_2 - p0_2;
    	
    	double p0w = (p0_0 * w0) + (p0_1 * w1) + (p0_0 * w2);
    	double pTw = (pT_0 * w0) + (pT_1 * w1) + (pT_0 * w2);
    	
    	double distSq = 0;
    	if( Math.signum(p0w) * Math.signum(pTw) == -1.0){
    		// min dist is p(s) . p(s)
    		// where: p(s) = p(0) - alpha(0) * w
    		// and    alpha(0) = p(0) . w / (w.w)
    		double wDotW  = (w0 * w0) + ( w1 * w1) + (w2 * w2);
    		double alpha0 = ((p0_0 * w0) + (p0_1 * w1) * (p0_2 * w2)) / wDotW;
    		double ps0    = p0_0 - (alpha0 * w0);
    		double ps1    = p0_1 - (alpha0 * w1);
    		double ps2    = p0_2 - (alpha0 * w2);
    		
    		distSq = (ps0 * ps0) + ( ps1 * ps1) + (ps2 * ps2);
    		
    	} else {
    		double dist0Sq =  (p0_0 * p0_0) + (p0_1 * p0_1) + (p0_2 * p0_2);
    		double distTSq =  (pT_0 * pT_0) + (pT_1 * pT_1) + (pT_2 * pT_2);
    		
    		distSq = Math.min(dist0Sq, distTSq);
    	}

        m_previousAPos[0] = pos[3 * m_flyByBodyA + 0];
    	m_previousAPos[1] = pos[3 * m_flyByBodyA + 1];
    	m_previousAPos[2] = pos[3 * m_flyByBodyA + 2];

    	m_previousBPos[0] = pos[3 * m_flyByBodyB + 0]; 
    	m_previousBPos[1] = pos[3 * m_flyByBodyB + 1]; 
    	m_previousBPos[2] = pos[3 * m_flyByBodyB + 2]; 

    	return distSq;
    }
    //////////////////////////////////////////////////////////    
    public void setPreviousPos(){
    	
    	Body bodyA        = m_bodies.get(m_flyByBodyA);
    	Body bodyB        = m_bodies.get(m_flyByBodyB);
    	
    	m_previousAPos    = new double[3];
    	m_previousBPos    = new double[3];
    	
        m_previousAPos[0] = bodyA.getInitPos0();
    	m_previousAPos[1] = bodyA.getInitPos1();
    	m_previousAPos[2] = bodyA.getInitPos2();

    	m_previousBPos[0] = bodyB.getInitPos0(); 
    	m_previousBPos[1] = bodyB.getInitPos1(); 
    	m_previousBPos[2] = bodyB.getInitPos2(); 
    }
    /////////////////////////////////////////////////////////
    public void displayPos(double[] pos){
         System.out.println("Final Positions:");
         for(int b = 0; b < m_bodies.size(); b++){
            int threeB = 3 * b;
            System.out.println("    " 
                               + m_bodies.get(b).getName()        + ": " 
                               + Double.toString(pos[threeB + 0])  + ", "
                               + Double.toString(pos[threeB + 1])  + ", "
                               + Double.toString(pos[threeB + 2])         );                             
         }     
    }
    //////////////////////////////////////////////////////////    
    public double[] getAcceleration(int b, double[] matForces){
         int numBodies = m_bodies.size();
         double[] acceleration = new double[3];
         
         int threeBNumBodies = 3 * b * numBodies; // we do this one off multiplication rather than repeating it in the for loop.
         
         double oneOverMass = m_bodies.get(b).getOneOverMass();
         for( int a = 0; a < numBodies; a++){
        	  int threeA = 3 * a;
        	 
              acceleration[0] += matForces[ threeBNumBodies + threeA + 0] * oneOverMass;
              acceleration[1] += matForces[ threeBNumBodies + threeA + 1] * oneOverMass;
              acceleration[2] += matForces[ threeBNumBodies + threeA + 2] * oneOverMass;
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
    //////////////////////////////////////////////////////////    
    public void openFileForOutput(){
      if ( m_sendOutputToFile) {
         try{
              // If the directory does not exist could create it.
              m_fileForOutput = new PrintWriter(m_outputFilePath, "UTF-8");
              System.out.println("Sending output to file: " + m_outputFilePath);
              
              String firstLine = new String("");
              for(int b = 0; b < m_bodies.size(); b++){
                   Body body = m_bodies.get(b);
                   firstLine += body.getName() + ":" + Double.toString(body.getBodyRadius()) +", ";
              }
              m_fileForOutput.write(firstLine + m_lineSep);
              
              //m_fileForOutput.println("The first line");
         } catch(Exception e){
              System.out.println("When trying to open file for output got exception: ");
              System.out.println(e.getMessage());
         } // end of try...catch
      }    // end of if....
    }      // end of method openFileForOutput
    ///////////////////////////////////////////////////////////
    public void storeCurrentPos(int t, double[] currentPos){
      if ( !m_sendOutputToFile)
    	  return;
      else if( m_fileForOutput == null){
           if ( !m_alreadyReportedNullFile){
                System.out.println("Unable to write to file: " + m_outputFilePath);
                m_alreadyReportedNullFile = true;
           }
      } else if ( (t % m_outputInterval) == 0) {     
         
        String posDataStr = new String("");
        
        for( int b = 0; b < m_bodies.size(); b++){
        	 int threeB = 3 * b;
             
             posDataStr +=  Double.toString(currentPos[threeB + 0]) + "," 
                          + Double.toString(currentPos[threeB + 1]) + ","
                          + Double.toString(currentPos[threeB + 2]) + ",";
        } // possibly could remove the last ","
           
        try {
             m_fileForOutput.write(posDataStr + m_lineSep); 
        } catch (Exception e) {
               System.out.println("When trying to write to file got exception: " + e.getMessage());
        } // end of try catch
      }   // end of if ... else
    }     // end of method
    ////////////////////////////////////////////////////////////
    public void finalizeOutput(){
    	 if ( ! m_sendOutputToFile)
    		 return;
    	 else if ( m_fileForOutput == null)
              System.out.println("Found that the output file was null.");
         else {
              try {m_fileForOutput.close(); }
              catch(Exception e){ 
                    System.out.println("When trying to close output file, got exception: " + e.getMessage());
              } // end of try catch
         }      // end of if ... else
    }           // end of method finalizeOutput
    //////////////////////////////////////////////////////////////    
}   // end of class BodiesMoving


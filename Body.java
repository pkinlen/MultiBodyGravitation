//////////////////////////////////////////////////////////   
// This class mostly just contains the data that is used by the BodiesMoving class.
// It can be initialized with with cartesian or polar coordinates.
// When polar coords are used, a circular orbit round the sun is assumed.
public class Body{
     String m_name;
     double m_mass;          // kg
     double m_bodyRadius;    // meters, could be used for checking for collisions
     double m_orbitRadius;   
     double m_theta;         // initial angle
     
     double m_initPos0;      // initial position on axis 0 (meters)
     double m_initPos1;      // initial position on axis 1
     double m_initPos2;      // initial position on axis 2
     
     double m_initVel0;      // initial velocity along axis 1 ( meters/ second)
     double m_initVel1;
     double m_initVel2;
     
     double m_oneOverMass;
     
     //////////////////////////////////////////////////////////////
     public Body(String name, double mass, double bodyRadius, 
                  double initPos0, double initPos1, double initPos2,
                  double initVel0, double initVel1, double initVel2){

         m_name        = name;
         m_mass        = mass;
         m_bodyRadius  = bodyRadius;
         m_oneOverMass = 1.0 / m_mass;

         initialize(initPos0, initPos1, initPos2, initVel0, initVel1, initVel2);         
     }
     ///////////////////////////////////////////////////////////////////
     public void initialize(double initPos0, double initPos1, double initPos2,
                            double initVel0, double initVel1, double initVel2)  {
         
          m_initPos0    = initPos0;
          m_initPos1    = initPos1;
          m_initPos2    = initPos2;
          
          m_initVel0    = initVel0;
          m_initVel1    = initVel1;
          m_initVel2    = initVel2;    
          
          /*
          System.out.println("Have body:  name: "         +                 name
                                     + ", mass: "         + Double.toString(mass)
                                       + ", body radius: "  + Double.toString(bodyRadius)
                                     + ", pos 0,1,2:   "  + Double.toString(initPos0) + ", " 
                                                            + Double.toString(initPos1) + ", "
                                                            + Double.toString(initPos2) + ", "
                                       + " vel 0,1,2: "     + Double.toString(initVel0) + ", "
                                                            + Double.toString(initVel1) + ", "
                                                            + Double.toString(initVel2) );
          */
     }
     //////////////////////////////////////////////////////////////
     public Body(String name, double mass, double bodyRadius, double orbitRadius, double theta){
         
         m_name        = name;
         m_mass        = mass;
         m_bodyRadius  = bodyRadius;
         m_orbitRadius = orbitRadius;
         m_theta       = theta;
         
         m_oneOverMass = 1.0 / m_mass;
         
         setPosAndVel(orbitRadius, theta);
     }
     //////////////////////////////////////////////////////////////////
     public void setPosAndVel(double orbitRadius, double theta){
         double initPos0    = orbitRadius * Math.cos(theta);
         double initPos1    = orbitRadius * Math.sin(theta);
         double initPos2    = 0.0;
         
         double initVel0, initVel1, initVel2;
         
         /*
          *     F = ma
          *     G M m / r^2 = m * omega^2 * r
          *     omega = sqrt( G M / r^3)
          *     x         = r                (  cos( theta + omega * t), sin( theta + omega * t) )
          *     v = dx/dt = r * omega      * ( -sin( theta + omega *t) , cos( theta + omega * t) )
          *               = sqrt( G M / r) * ( -sin( theta + omega *t) , cos( theta + omega * t) )
          * 
          */
         if( orbitRadius > 0) {
             double sqrtGMoverR = Math.sqrt( Config.m_massOfSun * Config.m_gravitationalConst / orbitRadius);
         
             initVel0    = - Math.sin(theta) * sqrtGMoverR;
             initVel1    =   Math.cos(theta) * sqrtGMoverR;
             initVel2    =   0.0;    
         } else {
             initVel0    = 0.0;
             initVel1    = 0.0;
             initVel2    = 0.0;             
         }
         
         initialize(initPos0, initPos1, initPos2, initVel0, initVel1, initVel2);
          
     }     
     //////////////////////////////////////////////////////////////
     public void resetTheta(double theta){ 
         setPosAndVel(m_orbitRadius, theta);
     }
     //////////////////////////////////////////////////////////////
     public String   getName()        { return m_name;        }
     public double   getMass()        { return m_mass;        }
     public double   getBodyRadius()  { return m_bodyRadius;  }
     public double   getTheta()       { return m_theta;       }
     
     public double   getInitPos0()    { return m_initPos0;    }
     public double   getInitPos1()    { return m_initPos1;    }
     public double   getInitPos2()    { return m_initPos2;    }
     
     public double   getInitVel0()    { return m_initVel0;    }
     public double   getInitVel1()    { return m_initVel1;    }
     public double   getInitVel2()    { return m_initVel2;    }
     
     public double   getOneOverMass() { return m_oneOverMass; }
     //////////////////////////////////////////////////////////////
}    // end of class Body

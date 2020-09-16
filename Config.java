
public class Config {
    private double  m_flyByTargetDistance; // The closest distance we target the meteorite to pass by each planet.
	
    private double  m_timeStepInterval;    // seconds, ( model seems to work when when this kept well below shortest_orbit_period / 100
                                           // For the moon, the orbit_period is just less than 30 days.
                                           // Could use 1 hour for m_timeStepInterval i.e. 3600 sec.
    
    private int     m_numTimeSteps;        // For 11 bodies, with 100,000 time steps, calc time is about 0.7 sec
    
    private int     m_outputInterval;      // 1 means every time-step will be output, 2 means every second one, etc
                                           // The number of time-steps in the output file will be: 
                                           //      m_numTimeSteps / m_outputInterval  
    private boolean m_displayFinalPos;
    
	private boolean m_sendOutputToFile;
    private String  m_outputFilePath;
    
    private String  m_lineSep;
    
    public static final double m_gravitationalConst = 6.67408e-11;   // m^3 kg^-1 s^-2
    public static final double m_massOfSun          = 1.989e30;      // kg
    public static final double m_massOfEarth        = 5.972e24;      // kg

	/////////////////////////////////////////////
	public Config(double  flyByTargetDistance,
			      double  timeStepInterval, int     numTimeSteps,       int    outputInterval, 
			      boolean displayFinalPos,  boolean sendOutputToFile, String outputFilePath  ){
		
		m_flyByTargetDistance = flyByTargetDistance; 
		m_timeStepInterval    = timeStepInterval;
		m_numTimeSteps        = numTimeSteps; 
		m_outputInterval      = outputInterval;
		m_displayFinalPos     = displayFinalPos;
		m_sendOutputToFile    = sendOutputToFile;
		m_outputFilePath      = outputFilePath;
		
		m_lineSep             = System.getProperty( "line.separator" );
	}
	/////////////////////////////////////////////
	public double   getFlyByTargetDistance() { return m_flyByTargetDistance;  }
    public double   getTimeStepInterval()    { return m_timeStepInterval;     }
    public int      getNumTimeSteps()        { return m_numTimeSteps;         }
    public int      getOutputInterval()      { return m_outputInterval;       }
    public boolean  displayFinalPos()        { return m_displayFinalPos;      }
    public boolean  sendOutputToFile()       { return m_sendOutputToFile;     }
    public String   getOutputFilePath()      { return m_outputFilePath;       }
    public String   getLineSep()             { return m_lineSep;              }
	/////////////////////////////////////////////
}

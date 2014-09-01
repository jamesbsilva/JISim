package JISim;

/**
 * 
 *  @(#) MCSimulation
 */  
 
import Backbone.System.SimpleLattice;
import Backbone.System.LatticeMagInt;
import Backbone.Algo.MetropolisMC;
import Backbone.Algo.IsingMC;
import Backbone.System.SimSystem;
import Backbone.Util.DataSaver;
import Backbone.Util.SimProcessParser;
import Backbone.Util.ParameterBank;
import java.util.Random;

/**  
 *   A monte carlo simulation that is meant to be extended and overwritten but
 * not implemented as interface to also be a bare bones example of a simulation
 *  <br>
 * 
 *  @param nRuns - number of runs for determining time of nucleation
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2011-09    
 */



public class MCSimulation {
    public double temperature;
    public double jInteraction;
    public double hField;
    public boolean useLongRange;
    public boolean useDilution;
    public IsingMC mc;
    public int currentSeed;
    public Random Ran = new Random();
    public int currentRun;
    public int nRuns = 50;
    public int instId = 1; // Instance Id - Necessary for running multiple instances

    private SimSystem simSys;
    private DataSaver dSave;
    private SimProcessParser parser;
    private ParameterBank param;
    private int tProcess;
    private int currentTime=0;
    
    public MCSimulation(){}

    public MCSimulation(int runs){nRuns = runs;}

    /**
    *      getSimSystem returns the simulations simSys object
    * 
    *   @see  LatticeMagInt
    */
    public SimSystem getSimSystem(){return simSys;}
    
    /**
    *   getHfield returns the current magnetic field of the simulation
    * 
    *   @return magnetic field value
    */ 
    public double getHfield(){return hField;}
    
    /**
    *      getRuns returns the run number of the current simulation run
    * @return - runs completed
    */    
    public int getRuns(){return currentRun;}
    /**
    *      getParamPostFix returns the postfix of the parameters file for this 
    *   simulation
    */
    public String getParamPostFix(){return "";};
    
    /**
    *      getParams returns the  parameters for this 
    *   simulation
    */
    public ParameterBank getParams(){return param;};
    
    /**
    *      getFixedPostFix returns the postfix of the fixed spin configuration
    *   file for this simulation
    */
    public String getFixedPostFix(){return "";};
    /**
    *      getSimPostFix returns the postfix of the simulation process file for this 
    *   simulation
    */
    public String getSimPostFix(){return "";};
    
    /**
    *         initialize performs all the initialization tasks necessary to
    *    run a simulation run.
    */
    public void initialize(){initialize("");}


    /**
    *      Close updates any configuration or settings before closing program 
    *  instance
    */
    public void Close(){dSave.updateClosedInstance();}

    /**
    *         initialize performs all the initialization tasks necessary to
    *    run a simulation run.
    * 
    *   @param type - post fix for the filename of the parameters file
    */
    public void initialize(String type){
        simSys = new SimpleLattice(param.s,type); 

        //which algo is going to be used
        if(param.mcalgo.equalsIgnoreCase("metropolis")){
                mc = new MetropolisMC(this);
                mc.setRun(currentRun);
        }else{}

        // get a random seed and save it
        currentSeed = Ran.nextInt();
        mc.setSeed(currentSeed);
    }

    
    /**
    *       runSimulation performs the main logic of the whole simulation by 
    *   calling on one simulation after it initializes it
    */
    public void runSimulation(){
        initialize();
        initialRun();
    }


    /**
    *       initialRun performs the first simulation of the system 
    */
    public void initialRun(){
        tProcess = parser.nextSimProcessTime();
        boolean goRun=true;
        while(goRun){
        if(currentTime == tProcess){
                mc = parser.simProcess(mc);
                tProcess = parser.nextSimProcessTime();
        }
        mc.doOneStep();
        currentTime++;
        goRun = mc.nucleated();
            }
    }

    /**
    *       setdsave sets the datasaver object for this simulation .
    * 
    * @see dataSaver
    * @param dSave - input datasaver
    */
    public void setdSave(DataSaver dSave) {
        this.dSave = dSave;
    }

    /**
    *       setdsave gets the datasaver object for this simulation .
    * 
    * @see dataSaver
    * @return  datasaver for simulation
    */
    public DataSaver getdSave() {
        return dSave;
    }

}

package Examples;

/**
 * 
 *  @(#) MCSimulationMovieRuns 1.1  2011-01-17
 */  
 
import Backbone.System.*;
import Backbone.Algo.*;
import Backbone.Util.*;
import AnalysisAndVideoBackend.VideoMakerCV;
import GPUBackend.*;
import JISim.*;
import java.io.*;
import java.util.Random;
/**  
 *   A monte carlo simulation that uses makes movies with different fixed spin
 *  configurations. It does a simulation for a Cross,Square,Box,Circle, and a 
 *  Crack.
 *  <br>
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2011-10    
 */


public final class MCSimulationMovieRuns extends MCSimulation {
    public double temperature;
    public double jInteraction;
    public double hField;
    public boolean useLongRange;
    public boolean useDilution;
    public IsingMC mc;
    private int instId;
    public int currentSeed;
    public Random Ran = new Random();
    public int currentRun=0;
    public int tNucleation;
    public int delPrecision = 3;
    public String fNameData;


    private SimProcessParser parser;
    private ParameterBank param;
    private DataSaver dSave;
    private LatticeMagInt lattice;
    private int tProcess;
    private int currentTime=0;
    private VideoMakerCV  vid;
    private int VideoLength = 0;
    private int FrameRate = 3; // How often to take a capture 
    private boolean makeVideo=true; // Do you want to record initial run
    private boolean useGPU=false; 
    private boolean output=true;
    private int framesPerCapture =1; // each capture is this many frames, make each capture last longer
    private int maxT = 10000;
    private String FixedPostFix = "";


    /**
    *      getSimSystem returns the simulations lattice object
    * 
    *   @see  LatticeMagInt
    */
    public LatticeMagInt getSimSystem(){return lattice;}

    public MCSimulationMovieRuns(){
        // Check for Files 
        (new DirAndFileStructure()).checkForSettingsFiles(this);
        
        param = new ParameterBank(getParamPostFix());
        param.printParameters();

        // Allow for multiple instances
        dSave = new DataSaver();
        instId = dSave.getInstId();
        param = new ParameterBank(getParamPostFix());dSave = new DataSaver(instId,param);
        System.out.println("Instance of program #"+instId);
        System.out.println("______________________________");

        parser = new SimProcessParser(parser,getSimPostFix(),getParamPostFix());
    }

    
    /**
    *      getParamPostFix returns the postfix of the parameters file for this 
    *   simulation
    */
    @Override
    public String getParamPostFix(){return "Movie";};
    /**
    *      getFixedPostFix returns the postfix of the fixed spin configuration
    *   file for this simulation
    */
    @Override
    public String getFixedPostFix(){return FixedPostFix;};
    /**
    *      getSimPostFix returns the postfix of the simulation process file for this 
    *   simulation
    */
    @Override
    public String getSimPostFix(){return "";};
    
    /**
    *      Close updates any configuration or settings before closing program 
    *  instance
    */
    @Override
    public void Close(){dSave.updateClosedInstance();}

    
    /**
    *         initialize performs all the initialization tasks necessary to
    *    run a simulation run.
    */
    @Override
    public void initialize(){
        initialize("");
    }
    /**
    *         initialize performs all the initialization tasks necessary to
    *    run a simulation run.
    * 
    *   @param type - filename of fixed spin configuration
    */
    @Override
    public void initialize(String type){
        lattice = new AtomicLatticeSumSpin(param.s,"Movie",type,instId); 
        parser = new SimProcessParser(parser,getSimPostFix(),getParamPostFix());
        currentTime =0;
        FixedPostFix=type;

        if(type.equalsIgnoreCase("")){
            MakeFixedConfig.clearFixedConfig2D(getFixedPostFix());}
        
        
        //which algo is going to be used
        if(param.mcalgo.equalsIgnoreCase("metropolis")){
            if(useGPU){
                // clear up context on GPU if creating a new context
                if(mc != null){
                    ((CLMetropolisMC)mc).closeOpenCL();
                }else{}
                mc = new CLMetropolisMC(this,output);                    
            }else{
                mc = new MetropolisMC(this,output);
            }
        }else if (param.mcalgo.equalsIgnoreCase("metropolisfix")){
            mc = new MetropolisFixedStrengthMC(this,output);
        }else{}

        // get a random seed and save it
        currentSeed = Math.abs(Ran.nextInt());
        //currentSeed=387483474;
        //currentSeed=98437549;

        mc.setSeed(currentSeed);
        System.out.println("Seed:"+currentSeed);	

        // Make Video
        if(makeVideo){
            vid = new VideoMakerCV(currentSeed,VideoLength,type+"MOVIERUN");
            vid.setFramesPerCapture(framesPerCapture);
        }
    }

    
    /**
    *       runSimulation performs the main logic of the whole simulation by 
    *   calling on one simulation after it initializes it
    */
    @Override
    public void runSimulation(){
            initialize();
            initialRun();
            initialize("square");
            initialRun();
            initialize("circle");
            initialRun();
            initialize("crack");
            initialRun();   
            initialize("cross");
            initialRun();
    }

    
    /**
    *       initialRun performs the first simulation of the system 
    */
    @Override
    public void initialRun(){
        tProcess = parser.nextSimProcessTime();
        System.out.println();
        boolean goRun=true;
        while(goRun){

        if(currentTime == tProcess){
                mc = parser.simProcess(mc);
                tProcess = parser.nextSimProcessTime();
        }
        mc.doOneStep();

        if(makeVideo){
        if(currentTime%FrameRate ==0){

            vid.addLatticeFrame(mc.getSimSystem().getSystemImg());
            makeVideo = vid.isWritten();
            }
        }

        if((currentTime % 10) ==0){System.out.println("t:"+currentTime+"  M:"+mc.getM());}
        currentTime++;


        goRun = !(mc.nucleated());
        if(mc.getM()<0 && currentTime> 300){goRun=false;}
        if(currentTime > maxT){goRun= false;}
        }

        if(makeVideo){
            vid.writeVideo();
        }
    }



    // test the class
    public static void main(String[] args) {
        MCSimulationMovieRuns sim = new MCSimulationMovieRuns();
        sim.runSimulation();
        sim.Close();
        System.out.println("Done!!!");
    }

}

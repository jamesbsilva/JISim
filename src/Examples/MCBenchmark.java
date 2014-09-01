package Examples;

/**
 * 
 *    @(#)   MCBenchmark
 */  
import Backbone.System.*;
import Backbone.Algo.*;
import Backbone.Util.*;
import GPUBackend.*;
import JISim.MCSimulation;
import JISim.MakeFixedConfig;
import java.util.Random;
import java.util.Scanner;

/**
 *      MCBenchmark class runs a simulation to test speed of one algorithm step.
 *  . It takes in a fixed config file post script name as input. Example
 *  input "cross" to use "fixedcross.txt".
 * 
 * <br>
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2012-05
 */
public final class MCBenchmark extends MCSimulation {

    private double temperature;
    private double jInteraction;
    private boolean useLongRange;
    private boolean useDilution;
    private IsingMC mc;
    private int instId;
    private int currentSeed;
    private Random Ran = new Random();
    private String fNameData;
    
    private SimProcessParser parser;
    private ParameterBank param;
    private DataSaver dSave;
    private LatticeMagInt lattice;
    private int tProcess;
    private int currentTime=0;
    private int maxL=2049;
    private int minT=25;
    private int timeN=25;
    private int initL = 512;
    private int initR = 1;
    private int rMax = 84;
    private int fixedSpins = 0;
    private boolean useGPU = true;
    private boolean nearestNeighborOnly = false;
    private boolean fixedL = true;
    private double currentTemp=1.56;
    private double avgTime=0;
    private double sigmaTime=0;
    private int spinFix=-1;
    private String configExplored="crack";

    /**
    *      MCBenchmark constructor
    * 
    * @param conf - fixed spin configuration filename 
    */
    public MCBenchmark(String conf){
        // Allow for multiple instances
        dSave = new DataSaver();
        instId = dSave.getInstId();
        // Check for Files 
        (new DirAndFileStructure()).checkForSettingsFiles(this);
        (new DirAndFileStructure()).createTempSettingsFile(getParamPostFixBase(), getParamPostFix());
        
        System.out.println("Instance of program #"+instId);
        System.out.println("______________________________");
        
        param = new ParameterBank(getParamPostFix());
        dSave = new DataSaver(instId,getParamPostFix());
        
        // enforce zero field
        param.changeField(0, getParamPostFix());
        
        configExplored =conf;
        
        // Check for Files 
        (new DirAndFileStructure()).checkForSettingsFiles(this);
        
        if(fixedSpins==0){
            MakeFixedConfig.clearFixedConfig2D(getFixedPostFix());
        }else{
            MakeFixedConfig.makeFixedConfig2D(
            configExplored,getFixedPostFix(), spinFix,(int)(param.L/2),(int)(param.L/2),param.L,fixedSpins,fixedSpins);
        }
        param.printParameters();
     
        parser = new SimProcessParser(parser,getSimPostFix(),getParamPostFix());
    }
    
    /**
    *      getSimSystem returns the simulations lattice object
    * 
    *   @see  LatticeMagInt
    */
    @Override
    public LatticeMagInt getSimSystem(){return lattice;}

    
    /**
    *      getParamPostFix returns the postfix of the parameters file for this 
    *   simulation.
    */
    @Override
    public String getParamPostFix(){return ""+instId+"-deleteable-";}
    /**
    *      getParamPostFix returns the postfix of the parameters file for this 
    *   simulation.
    */
    public String getParamPostFixBase(){return "";}    
    /**
    *      getFixedPostFix returns the postfix of the fixed spin configuration
    *   file for this simulation.
    */
    @Override
    public String getFixedPostFix(){return (""+configExplored+"-"+instId);};
    /**
    *      getSimPostFix returns the postfix of the simulation process file for this 
    *   simulation.
    */
    @Override
    public String getSimPostFix(){return "";};
    
    
    /**
    *      Close updates any configuration or settings before closing program 
    *  instance.
    */
    @Override
    public void Close(){
        dSave.updateClosedInstance();
        if(useGPU){
            if(param.mcalgo.contains("fully")){
                ((FullyConnectMetroMCCL)mc).closeOpenCL();
            }else{
                ((CLMetropolisMC)mc).closeOpenCL();
            }
        }
    }

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
    *   @param temp - the temperature of this simulation
    *    @param  type - the fixed spin configuration of this simulation
    */
    public void initialize(int lnew,int rnew, String type){

        param.changeLength(lnew,getParamPostFix());
        param.changeRange(rnew,getParamPostFix());
        dSave = new DataSaver(instId,getParamPostFix());

        if(!useGPU && !param.mcalgo.contains("fully")){
            lattice = new AtomicLatticeSumSpin(param.s,getParamPostFix(),getFixedPostFix(),instId);
        } 
        parser = new SimProcessParser(parser,getSimPostFix(),getParamPostFix());
        currentTime =0;

        //which algo is going to be used
        if(param.mcalgo.equalsIgnoreCase("metropolis")){
            if(useGPU){
                // clear up context on GPU if creating a new context
                if(mc != null){
                    ((CLMetropolisMC)mc).closeOpenCL();
                }else{}
                mc = new CLMetropolisMC(this,false);                    
            }else{
                mc = new MetropolisMC(this,false);
            }
        }else if(param.mcalgo.equalsIgnoreCase("metropolisfully")){
            if(useGPU){
                // clear up context on GPU if creating a new context
                if(mc != null){
                    ((FullyConnectMetroMCCL)mc).closeOpenCL();
                }else{}
                mc = new FullyConnectMetroMCCL(this,false,0,spinFix,false, 256);                   
            }else{
                mc = new MetropolisFullyConnectedMC(this,false);
            }
        }else{}
        // turn off trigger
        mc.setTriggerOnOff(false);
        
        // get a random seed and save it
        //currentSeed = Math.abs(Ran.nextInt());
        currentSeed=387483474;
        //currentSeed=98437549;
        
        mc.setSeed(currentSeed);	
    }

    
    /**
    *       runSimulation performs the main logic of the whole simulation by 
    *   calling on one simulation after it initializes it.
    */
    @Override
    public void runSimulation(){
        boolean searching =true;
        if(nearestNeighborOnly){
            initL = 32;
            while(searching){
                initL = 2*initL; 
                System.out.println("********************************");
                System.out.println("Starting Length: "+initL);
                System.out.println("********************************");
                initialize(initL,0,configExplored);
                initialRun();
                saveResult(initL,0);
                searching = continueSearch(initL,0);
            }
        }else if (!nearestNeighborOnly && !fixedL){            
            initL = 32;
            while(searching){
                initL = 2*initL; 
                System.out.println("********************************");
                System.out.println("Starting Length: "+initL);
                System.out.println("********************************");
                initialize(initL,initR,configExplored);
                initialRun();
                saveResult(initL,initR);
                searching = continueSearch(initL,param.R);
            }
        }else if (!nearestNeighborOnly && fixedL){            
            int per =0;
            int range = initR;//0;
            //rMax = (int) (initL*16/100);
            while(searching){
                per = 1+per;
                range *= 2; //(int) (initL*per/100);
                System.out.println("********************************");
                System.out.println("Starting R: "+range);
                System.out.println("********************************");
                initialize(initL,range,configExplored);
                initialRun();
                saveResult(initL,range);
                searching = continueSearch(initL,range);
            }
        }
        printResult(initL);
    }
    
    /**
    *       saveResult  prepares the data to be saved.
    */
    private void saveResult(int lnow,int rnow){	
        double[] meta = new double[8];
        meta[0] = lnow;
        meta[1] = rnow;
        meta[2] = avgTime;
        meta[3] = sigmaTime;
        meta[4] = timeN;
        meta[5] = param.hField;
        meta[6] = currentTemp;
        meta[7] = currentSeed;
        saveData(meta);
    }

    
    /**
    *       saveData saves the data that has been prepared to go out to file.
    * @param tcdata - data to be saved 
    */
    public void saveData(double[] tcdata){
        setDataFileName(configExplored+"-nfix-"+fixedSpins);
        dSave.saveDoubleArrayData(tcdata,fNameData);
    }

    /**
    *       setDataFileName sets the filename of the outgoing data file.
    * @param fix - filename of data file
    */
    public void setDataFileName(String fix){
        if(useGPU){
            fNameData =  "BenchmarkGPUData-"+param.mcalgo+".txt";
        }else{                
            fNameData =  "BenchmarkData-"+param.mcalgo+".txt";
        }
    }
  
    /**
     *      printResult prints the best range to run simulation for given max time
     *  and min time. It does this print out in the console.
     */
    private void printResult(int lnow){
        System.out.println("Avg Step Time:"+avgTime+"      L:"+lnow);
    }
    
    /**
    *       continueSearch performs the logic to determine if the simulation 
    *   should proceed      
    * 
    * @return true if continue simulation 
    */
    public boolean continueSearch(int lnow, int rnow){
        boolean go;
        if((lnow) >= maxL || rnow >= (rMax*lnow/100)){
            go=false;
        }else{go=true;}
        return go;
     }


    
    /**
    *       initialRun performs the first simulation of the system. 
    */
    @Override
    public void initialRun(){
        tProcess = parser.nextSimProcessTime();
        System.out.println();
        boolean goRun=true;
        int maxT = minT+timeN;
        double time = 0;
        double time2 = 0;
        long nanot;
        while(goRun){
            nanot = System.nanoTime();        
            mc.doOneStep();
            nanot = System.nanoTime()-nanot;
            if(currentTime > minT){
                time += (double)(nanot/1000);
                time2 += Math.pow(time, 2);
            }

            if(currentTime > maxT){
                goRun = false;
            }

            if((currentTime % 1) == 0 && currentTime !=0){
                System.out.println("t:"+currentTime+"  M:"+mc.getM());
            }

            currentTime++;
        }
        avgTime = time/timeN;
        sigmaTime = Math.pow((((time2/timeN) - Math.pow(avgTime,2))/timeN),0.5);
    }

// test the class
public static void main(String[] args) {
    boolean ResponseNotDone = true; 
    String response;
    Scanner scan = new Scanner(System.in);

    String fp="random";


    while(ResponseNotDone){

    System.out.println("Default Parameters? (Clear Instances by typing clearInst)");	
    response = scan.next();
    if(response.equalsIgnoreCase("y") || response.equalsIgnoreCase("")){
        ResponseNotDone=false;
    }else if(response.equalsIgnoreCase("n")){
        ResponseNotDone=false;

        System.out.println("Fixed file post?(e for empty input)");
        response = scan.next();
        if(response.equalsIgnoreCase("e")){fp="";}else{
        fp = response;}
    }else if(response.equalsIgnoreCase("clearInst")){
        DataSaver dsave = new DataSaver();
        dsave.clearInstances();
    }else{
        System.out.println("y or n");}
    }
        MCBenchmark sim = new MCBenchmark(fp);
        sim.runSimulation();
        sim.Close();
        System.out.println("*************************************");
        System.out.println("------Done!!!--------");
        System.out.println("*************************************");
    }
}

package Examples;

/**
*    @(#)   MCFieldExplorer
*/  

import Backbone.System.*;
import Backbone.Algo.*;
import Backbone.Util.*;
import JISim.*;
import java.util.Random;
import java.util.Scanner;

/**
*      MCFieldExporer class searches for good hfield range to run simulations
*  on. It takes in a fixed config file post script name as input. Example
*  input "cross" to use "fixedcross.txt".
* 
* <br>
* 
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2012-01
*/
public final class MCFieldExplorer extends MCSimulation {
    private double temperature; private double jInteraction;
    private boolean useLongRange; private boolean useDilution;
    private int instId; private int currentSeed; private boolean singleField = true;
    private IsingMC mc; private boolean useSameConfig= false; private boolean triggerOn = false;
    private int spinFix=-1; private int nfixed = 1000;
    private Random Ran = new Random(); private String fNameData;
    private int tNucleationLast = 0; private int tNucleation = 0;
    private boolean takeSnapshot=true; private int Run = 0;

    private SimProcessParser parser; private ParameterBank param; private DataSaver dSave;
    private int maxT=28000; private int minT=15000; private int fixedSeed;
    private LatticeMagInt lattice; private int runsPerField = 100;
    private int tProcess; private int currentTime=0;
    private double delH = 0.1; private double minDelH =0.01; private double hRight = 1.0;
    private String configExplored; private boolean output =true;
    
    /**
    *      MCFieldExplorer constructor
    * 
    * @param conf - fixed spin configuration filename 
    */
    public MCFieldExplorer(String conf,double xper, double hr, double hd, boolean sf){
        // Allow for multiple instances
        dSave = new DataSaver();instId = dSave.getInstId();
        singleField = sf; hRight = hr; delH = hd;
        // Check for Files 
        (new DirAndFileStructure()).checkForSettingsFiles(this);
        (new DirAndFileStructure()).createTempSettingsFile(getParamPostFixBase(), getParamPostFix());
        param = new ParameterBank(getParamPostFix());
        dSave = new DataSaver(instId,param);
        //param.printParameters();
        configExplored = conf; nfixed = (int)(xper*((double)param.N));
        
        //double temp = ((4.0-((4.0-2.269)*Math.pow(param.R,-2.0)))*4/9)*(1-((double)nfixed/(double)Math.pow(param.L, param.D)));
        double temp = (4.0)*(1-((double)nfixed/(double)Math.pow(param.L, param.D)));
        param.changeTemp(temp, getParamPostFix());
        
        System.out.println("Instance of program #"+instId);
        System.out.println("__________________________________________");
        parser = new SimProcessParser(parser,getSimPostFix(),getParamPostFix(),param);
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
    *   simulation
    */
    @Override
    public String getParamPostFix(){return (getParamPostFixBase()+instId+"-deleteable-");}
    
    /**
    *      getParamPostFix returns the postfix of the parameters file for this 
    *   simulation
    */
    public String getParamPostFixBase(){return "-Explorer";};
    
    /**
    *      getFixedPostFix returns the postfix of the fixed spin configuration
    *   file for this simulation
    */
    @Override
    public String getFixedPostFix(){return (""+configExplored);};
    
    /**
    *      getSimPostFix returns the postfix of the simulation process file for this 
    *   simulation
    */
    @Override
    public String getSimPostFix(){return "Explorer";};
    
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
    *   @param hfield - the magnetic field of this simulation
    *    @param  type - the fixed spin configuration of this simulation
    */
    public void initialize(double hfield,String type, int run){
        param.changeField(hfield,getParamPostFix());
        parser = new SimProcessParser(parser,getSimPostFix(),getParamPostFix(),param);
        currentTime = 0;

        if(!useSameConfig || run == 0){
            fixedSeed = Math.abs(Ran.nextInt());
            if(nfixed==0){
                MakeFixedConfig.clearFixedConfig2D(getFixedPostFix());
            }else{
                if(configExplored.contains("interactionquilt")){   
                    MakeFixedConfig.makeFixedConfig2D(
                    configExplored,getFixedPostFix(), spinFix,(int)(param.L/2),(int)(param.L/2),param.L,param.R,nfixed,fixedSeed);
                }else{
                    MakeFixedConfig.makeFixedConfig2D(
                    configExplored,getFixedPostFix(), spinFix,(int)(param.L/2),(int)(param.L/2),param.L,nfixed,nfixed,fixedSeed);
                }
            }
        }
        
        lattice = new AtomicLatticeSumSpin(param.s,getParamPostFix(),type,instId); 
        //which algo is going to be used
        if(param.mcalgo.equalsIgnoreCase("metropolis")){
            //mc = new CLMetropolisMC(this,output,"GPU");
            mc = new MetropolisMC(this,output);
        }
        if(!triggerOn){mc.setTriggerOnOff(triggerOn);}

        System.out.println("TEMP : "+param.temperature+"    h : "+param.hField);
        output=false;
        // get a random seed and save it
        currentSeed = Math.abs(Ran.nextInt());
        mc.setSeed(currentSeed);	
    }

    /**
    *       runSimulation performs the main logic of the whole simulation by 
    *   calling on one simulation after it initializes it
    */
    @Override
    public void runSimulation(){
        boolean searching =true;
        while(searching){
            for(int u = 0; u < this.runsPerField;u++){
                Run = u;
                initialize(hRight-delH,configExplored,u);
                tNucleationLast = tNucleation;
                System.out.println("+++++++++++++++++++++++++++++++++++++++++++++");
                System.out.println("Run : "+u+"      hField: "+(hRight-delH));
                initialRun();
                System.out.println("_______________________________________");
                System.out.println("tNucleationLeft: "+tNucleation+"     hLeft: "+(hRight-delH)+
                        "      hRight: "+hRight+ "    tNucleationRight: "+tNucleationLast);
                System.out.println("_______________________________________");
            }
            searching = continueSearch();
            if( singleField ){searching = false;}
        }
        printResult();
        (new DirAndFileStructure()).deleteTempSettingsFile(getParamPostFix());        
    }

    /**
    *      printResult prints the best range to run simulation for given max time
    *  and min time. It does this print out in the console.
    */
    private void printResult(){
        System.out.println("Field Range that is good in tNucleation range");
        System.out.println("timeMin:"+minT+"      timeMax:"+maxT);
        System.out.println("FieldL:"+(hRight-delH)+"      FieldR:"+hRight);			
    }
    
    /**
    *       continueSearch performs the logic to determine if the simulation 
    *   should proceed      
    * 
    * @return true if continue simulation 
    */
    public boolean continueSearch(){
        if(tNucleation > minT && tNucleation < maxT){
            return false;
        }else if(delH < minDelH){
            return false;
        }else{
            if(tNucleation < minT){
                hRight = hRight-delH;
            }else{
                delH = (int) (delH/2);
            }
            return true;
        }
    }
    
    /**
    *       initialRun performs the first simulation of the system 
    */
    @Override
    public void initialRun(){
        tProcess = parser.nextSimProcessTime();
        boolean goRun=true;
        while(goRun){
            if(currentTime == tProcess){
                mc = parser.simProcess(mc);
                tProcess = parser.nextSimProcessTime();
            }   
            mc.doOneStep();
            if((currentTime % 1000) == 0 && currentTime !=0){ System.out.println("t:"+currentTime+"  M:"+mc.getM());}
            
            currentTime++; goRun = !(mc.nucleated());
            if(currentTime > maxT){goRun= false;tNucleation = currentTime+1000;}
        }
        if(takeSnapshot){
            dSave.saveLatticeAsImg(maxT, Run ,currentSeed, (LatticeMagInt)mc.getSimSystem(),"Images/Tc/x-"
                +((int)(1000.0*(double)mc.getSimSystem().getNFixed()/(double)param.N))+"/hfield-"
                    +param.hField+"/fixseed-"+fixedSeed);
        }
        if(currentTime > maxT){
            tNucleation = currentTime+1000;
        }else{
            tNucleation = mc.getTriggerTime();
        }
    }

    // test the class
    public static void main(String[] args) {
        boolean ResponseNotDone = true; 
        String response;
        Scanner scan = new Scanner(System.in);
        String fp="Random";
        while(ResponseNotDone){
            System.out.println("Default Parameters? (Clear Instances by typing clearInst)");	
            response = scan.next();
            if(response.equalsIgnoreCase("y") || response.equalsIgnoreCase("")){
                ResponseNotDone=false;
            }else if(response.equalsIgnoreCase("n")){
                System.out.println("Fixed file post?(e for empty input)");
                ResponseNotDone=false; response = scan.next();
                if(response.equalsIgnoreCase("e")){
                    fp="";
                }else{
                    fp = response;
                }
            }else if(response.equalsIgnoreCase("clearInst")){
                DataSaver dsave = new DataSaver(); dsave.clearInstances();
            }else{
                System.out.println("y or n");
            }
        }
        double xper = 0.05; double hdel = 0.01; double hright = 4*xper+hdel;
        MCFieldExplorer sim = new MCFieldExplorer(fp,xper,hright,hdel,true);
        sim.runSimulation();
        sim.Close();
        System.out.println("Done!!!");
    }
}

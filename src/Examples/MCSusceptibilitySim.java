package Examples;

/**
*    @(#)   MCSusceptibilitySim
*/

import Backbone.System.*;
import Backbone.Algo.*;
import Backbone.System.SingleValueLattice;
import Backbone.Util.*;
import GPUBackend.*;
import JISim.MCSimulation;
import JISim.MakeFixedConfig;
import JISim.SetupJISim;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Random;
import java.util.Scanner;

/**
*      MCSusceptibilitySim is a simulation that measures the susceptibility for
*   multiple magnetic field values.
* 
* <br>
* 
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2012-02
*/
public final class MCSusceptibilitySim extends MCSimulation{
    private IsingMC mc;
    private int instId; private int currentSeed;
    private int tMeasure = 28000; private int tMinTime = 25000;
    private int runsPerField= 3;
    private double[] susceptibilityData;
    private double hBound = 1.5;
    private double hMin = 0.5; private double minDelH=0.001;
    private double hOffset = 0.1; private double delH =0.025;
    private double currentField=hMin; private int currentFieldRun=0;
    private int rIndex=0;
    private Random Ran = new Random();
    private int tNucleationLast=0;
    private String fNameData;
    private int spinFix = -1;
    private boolean hIncrease=true; private boolean takeSnapshot=true;
    private boolean useGPU = false; private boolean useSameConfig = false;
    private int currentFixedParam=0;
    private int mStable;

    private SimProcessParser parser; private ParameterBank param; private DataSaver dSave;
    private String susParamPostFix=""; private String susFixedPostFix="";
    private int tProcess; private int currentTime=0;
    private LatticeMagInt lattice;
    private boolean output =true;
    private String configExplored="random";
    
    /**
    *      MCSusceptibilitySim constructor
    * 
    * @param fixPost - fixed configuration file post fix to filename 
    * @param susPost - parameter file post fix to filename
    */
    public MCSusceptibilitySim(String fixPost,String susPost,int nfix,double ho,
            double hd,boolean inc,double bnd){
        // Allow for multiple instances
        dSave = new DataSaver();
        instId = dSave.getInstId();
        
        // Check for Files 
        (new DirAndFileStructure()).checkForSettingsFiles(this);
        (new DirAndFileStructure()).createTempSettingsFile(getParamPostFixBase(), getParamPostFix());
        
        susParamPostFix = susPost; susFixedPostFix=fixPost;
        
        currentFixedParam = nfix;
        hMin = ho; hIncrease = inc;
        delH = hd; hBound = bnd;
        
        // Compensate for change in critical Tc
        param = new ParameterBank(getParamPostFix());dSave = new DataSaver(instId,param);
        
        double TempSim= (4.0-((4.0-2.269)*Math.pow(param.R, -2.0)))*(1.0-((double)nfix/(Math.pow(param.L,2.0))))*4.0/9.0;
        if(param.mcalgo.contains("fully")){
            TempSim= (4.0)*(1.0-((double)nfix/(Math.pow(param.L,1.0))))*4.0/9.0;
        }
        // Asserting Nearest Neighbor case
        if(!param.useLongRange){TempSim = 2.269*4.0/9.0;}
        param.changeTemp(TempSim,getParamPostFix());
        param = new ParameterBank(getParamPostFix());
        param.printParameters();
        //configExplored
        
        System.out.println(configExplored+"    with param: "+currentFixedParam);
        
        if(currentFixedParam==0){
            MakeFixedConfig.clearFixedConfig2D(getFixedPostFix());
        }else if(!param.mcalgo.contains("fully")){
            if(configExplored.contains("interactionquilt")){   
                MakeFixedConfig.makeFixedConfig2D(
                configExplored,getFixedPostFix(), spinFix,1,1,param.L,currentFixedParam,param.R,currentSeed);
            }else{
                MakeFixedConfig.makeFixedConfig2D(
                configExplored,getFixedPostFix(), spinFix,(int)(param.L/2),(int)(param.L/2),param.L,currentFixedParam,currentFixedParam);
            }
        }
        
        susceptibilityData = new double[runsPerField];
        rIndex=0;
        parser = new SimProcessParser(parser,getSimPostFix(),getParamPostFix(),param);
        //
        System.out.println("Only outputting temp and flip changes in first run");
    }

    @Override
    /**
    *      getSimSystem returns the simulations lattice object
    * 
    *   @see  LatticeMagInt
    */
    public LatticeMagInt getSimSystem(){return lattice;}

    
    @Override
    /**
    *      getParamPostFix returns the postfix of the parameters file for this 
    *   simulation
    */
    public String getParamPostFix(){return (getParamPostFixBase()+instId+"-deleteable-");}
    /**
    *      getParamPostFix returns the postfix of the parameters file for this 
    *   simulation
    */
    private String getParamPostFixBase(){return ("-Susceptibility"+susParamPostFix);}
    
    public void deleteTempSettingsFile(){
        (new DirAndFileStructure()).deleteTempSettingsFile(getParamPostFix());
    }
    
    /**
    *      getFixedPostFix returns the postfix of the fixed spin configuration
    *   file for this simulation
    */
    public String getFixedPostFix(){return susFixedPostFix+instId;}
    @Override
    /**
    *      getSimPostFix returns the postfix of the simulation process file for this 
    *   simulation
    */
    public String getSimPostFix(){return "";}
    
    @Override
    /**
    *      Close updates any configuration or settings before closing program 
    *  instance
    */
    public void Close(){dSave.updateClosedInstance();}

    
    /**
    *         initialize performs all the initialization tasks necessary to
    *    run a simulation run.
    */
    @Override
    public void initialize(){
        // Change field for new susceptibility measurement
        param.changeField(currentField, getParamPostFix());
        if(param.mcalgo.contains("ully")){
            lattice = new SingleValueLattice(param.L,param.s,instId);
        }else{
            lattice = new AtomicLatticeSumSpin(param.s,getParamPostFix(),getFixedPostFix(),instId);
        } 
        
        parser = new SimProcessParser(parser,getSimPostFix(),getParamPostFix(),param);
        currentTime =0;

        
        if(!useSameConfig){
            if(currentFixedParam==0){
                MakeFixedConfig.clearFixedConfig2D(getFixedPostFix());
            }else if (!param.mcalgo.contains("fully")){
                if(configExplored.contains("interactionquilt")){   
                    MakeFixedConfig.makeFixedConfig2D(
                    configExplored,getFixedPostFix(), spinFix,(int)(param.L/2),(int)(param.L/2),param.L,currentFixedParam,param.R,Ran.nextInt());
                }else{
                    MakeFixedConfig.makeFixedConfig2D(
                    configExplored,getFixedPostFix(), spinFix,(int)(param.L/2),(int)(param.L/2),param.L,currentFixedParam,currentFixedParam,Ran.nextInt());
                }
            }
        }
        
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
        }else if(param.mcalgo.equalsIgnoreCase("metropolisfully")){
            boolean balanced = false;
            if(configExplored.contains("alanced")){balanced =true;}
            if(useGPU){
                mc = new FullyConnectMetroMCCL(this,output,currentFixedParam,spinFix,balanced, 256);
            }else{
                mc = new MetropolisFullyConnectedMC(this,output,currentFixedParam,param.L,spinFix,balanced);
            }
        }

        mStable = mc.getM();

        output=false;
        // get a random seed and save it
        currentSeed = Math.abs(Ran.nextInt());
        //currentSeed=387483474;
      
        // align field with current state
        mc.alignField();
      
        mc.setSeed(currentSeed);	
    }

    @Override
    public void runSimulation(){
        boolean searching =true;
        currentFieldRun=0;
        currentField = hMin;    
        while(searching){
            System.out.println("_______________________________________________");
            System.out.println("Current Field h: "+currentField+ "      fixed: "+currentFixedParam);
            System.out.println("_______________________________________________");
            for(int i=0;i<runsPerField;i++){
            initialize();
            initialRun(i);
            if(tNucleationLast > tMinTime){
            saveResult(i);}}
            searching = continueSearch();
        }
        (new DirAndFileStructure()).deleteTempSettingsFile(getParamPostFix());
    }


    /**
    *       saveResult saves the resulting susceptibility data to a text
    *   file in the data directory.
    */
    private void saveResult(int i){	
        double[] meta = new double[12];
        double zeros = 0;
        
        meta[0] = param.hField;
        meta[1] = susceptibilityData[i];
        meta[2] = tNucleationLast;
        meta[3] = currentFixedParam;
        meta[4] = (mc.getSimSystem().getNFixed());
        meta[5] = mc.getSimSystem().getN();
        meta[6] = (mc.getSimSystem().getClass().getName().contains("lattice"))? 
                ((LatticeMagInt)mc.getSimSystem()).getRange() : 0 ;
        meta[7] = param.L;
        meta[8] = param.temperature;
        meta[9] = param.jInteraction;
        meta[10] = (param.mcalgo.contains("ully"))? param.L : param.R;
        meta[11] = (configExplored.contains("alanced")) ? 0 : spinFix;
        meta[12] =  (configExplored.contains("alanced")) ? 1 : 0;
        meta[13] =  (configExplored.contains("alanced")) ? 1 : 0;
        setDataFileName(configExplored);
        saveData(meta);
        
    }

    /**
    *       continueSearch performs the logic to determine if the simulation 
    *   should proceed      
    * 
    * @return true if continue simulation 
    */
    public boolean continueSearch(){
        boolean cont = false;
        double del;
        // Change in delH function
        if(currentField< (hBound-hOffset)){
            del=delH;
        }else{
            currentFieldRun++;
            del = (delH/(1.0+currentFieldRun));
            if(del< minDelH){del=minDelH;}
        }
        
        System.out.println("adding "+del+"  to the field . Index for function "+currentFieldRun);
        if(hIncrease){
            currentField = currentField + del;
        }else{
            currentField = currentField - del;
        }
  
        if(hIncrease){
            if(currentField < hBound){cont=true;}
        }else{
            if(currentField > hBound){cont=true;}
        }
        
        if(tNucleationLast < tMinTime){
            System.out.println("Unstable area. Finish search. TRY INCREASING RANGE");
            cont=false;
        }
        return cont;
    }

    
    /**
    *       saveData saves the simulation data as a text file in the data directory. 
    */
    public void saveData(double[] metaTime){
        dSave.saveDoubleArrayData(metaTime,fNameData);
    }

    
    /**
    *       setDataFileName sets the filename of the output data file.
    * 
    * @param fix - fixed spin configuration type. 
    */
    public void setDataFileName(String fix){
        if(useGPU){
            fNameData =  "SusData/SusceptabilityDataGPU-"+mc.getClass().getSimpleName()+
                    "-"+fix+"-L-"+param.L+"-j-"+jInteraction+".txt";
        }else{
            fNameData =  "SusData/SusceptabilityData-"+mc.getClass().getSimpleName()+
                    "-"+fix+"-L-"+param.L+"-j-"+jInteraction+".txt";
        }
    }

    
    /**
    *       initialRun performs the simulation of the system
    */
    public void initialRun(int i){
        tProcess = parser.nextSimProcessTime();
        boolean goRun = true; boolean badTimeout = false;
        currentTime=0;
        
        while(goRun){
            if(currentTime == tProcess){
                //System.out.println("Process time");
                mc = parser.simProcess(mc);
                tProcess = parser.nextSimProcessTime();
            }
            mc.doOneStep();

            if((currentTime % 1000) == 0 && currentTime !=0){System.out.println("t:"+currentTime+"  M:"+mc.getM());}
            currentTime++;

            if(currentTime==tMeasure){
                susceptibilityData[i] = mc.getSusceptibility();
            }
       
            goRun = !(mc.nucleated());
            
            // Time out after taking data and output this fact
            if(currentTime ==tMeasure){
                goRun = false;
                if(mc.getM()< mStable*0.3){
                    System.out.println(" NOT properly timing out with M: "+mc.getM());
                    badTimeout = true;
                }else{
                    System.out.println("Run appears to properly time out with M: "+mc.getM());
                }
            }
        }
        if(!badTimeout && takeSnapshot){
            String f = configExplored+"-seed-"+currentSeed;
            dSave.saveConfigPermanent(mc.getSimSystem(), 0, currentTime, f,currentFixedParam,("Susceptibility-L-"+param.L));
        }
        tNucleationLast = currentTime;
        if(badTimeout){tNucleationLast = 1;} // Dont keep bad data or keep going
    }
    
    // test the class
    public static void main(String[] args) {
        boolean respond=true;
        String responsePost="";
        String responseConfig="";
        
        System.out.println("==========================================================================");
        System.out.println();
        System.out.println("Starting example Susceptibility measurement simulation on fully connected system");
        System.out.println();
        System.out.println("==========================================================================");
   
        Scanner scan = new Scanner(System.in);
        while(respond){
            System.out.println("Default File Parameters ? (y or n)");
            responsePost = scan.next();
            if(responsePost.equalsIgnoreCase("n")){
                System.out.println("post fix for fixed configuration? (e for crack (default))");
                responseConfig = scan.next();
                if(responseConfig.equalsIgnoreCase("e")){responseConfig="crack";}

                System.out.println("post fix for the parameters file? (e for empty (default))");
                responsePost = scan.next();
                if(responsePost.equalsIgnoreCase("e")){responsePost="";}
                respond = false;
            }else if(responsePost.equalsIgnoreCase("y")){
                responsePost = "";
                respond = false;
            }else{System.out.println("y or n");}
        }
        
        String postfix = responsePost;
        
        respond =true;
        double hMin =1.0;
        double delH = 0.015;
        int nfixed = 50;
        boolean inc=true;
        double bnd=1.30;
        
        while(respond){
            System.out.println("Default Sim Parameters ? (y or n)");
            responsePost = scan.next();

            if(responsePost.equalsIgnoreCase("n")){
                System.out.println("Start with which field? (hMin)");
                hMin = scan.nextDouble();
                System.out.println("Change field in intervals of which size? (delH)");
                delH = scan.nextDouble();
                System.out.println("Fix how many spins?");
                nfixed = scan.nextInt();
                System.out.println("Increasing field search? anything "
                        + "other 'n' is increasing confirmation ");
                responsePost = scan.next();
                if(responsePost.equalsIgnoreCase("n")){inc=false;bnd=0.0;}
                respond = false;
            }else if(responsePost.equalsIgnoreCase("y")){
                respond = false;
            }else{System.out.println("y or n");}
        }
            
        MCSusceptibilitySim mc = new MCSusceptibilitySim(responseConfig,postfix,nfixed,hMin,delH,inc,bnd);
        ParameterBank param = new ParameterBank(mc.getParamPostFix());
        mc.deleteTempSettingsFile();
        
        int iter = 8;
        int num = 0;
        if(param.mcalgo.contains("fully")){
            num= param.L;
        }else{
            num= (int)Math.pow(param.L,param.D);//L*L;
        }
        delH = 0.005;
        int perc = (int) (0.01*num);
        for(int i = 1;i<iter;i++){
            nfixed = 0+perc*i;
            // L=115 RB m = 1.2483 b = 1.637
            double p1 = -1.27 - 4.2;//-5.187*10*10*10*10*10; //-3.9187-0.0;
            double p2 = 1.27;//3.974; //1.27;
            bnd = p2 +delH;
            double p3 = 1.27;//1.247;           
            double x = (double)((double)nfixed/(double)(num));
            hMin =(p1*x + p2)- delH*6;
            //hMin = (p1*Math.pow(x,p2)+p3)- delH*5;
            mc = new MCSusceptibilitySim(responseConfig,postfix,nfixed,hMin,delH,inc,bnd);
            mc.runSimulation();
        }
    
        System.out.println("______________________________");
        System.out.println("Done!");
        System.out.println("______________________________");
    }
}

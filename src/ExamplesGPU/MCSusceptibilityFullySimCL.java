package ExamplesGPU;

/**
*    @(#)   MCSusceptibilityFullySimCL
*/

import Examples.*;
import Backbone.System.*;
import Backbone.Algo.*;
import Backbone.Util.*;
import GPUBackend.*;
import JISim.MCSimulation;
import JISim.MakeFixedConfig;
import JISim.SetupJISim;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;


/**
*      MCSusceptibilityFullySimCL is a simulation that measures the susceptibility for
*   multiple magnetic field values.
* 
* <br>
* 
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2012-02
*/
public final class MCSusceptibilityFullySimCL extends MCSimulation{
    private IsingMC mc;  private int instId;
    private int currentSeed;
    private int tMeasure = 28000; private int tMinTime = 23000;
    private int runsPerField= 32;
    private int minAdded = 5;
    private int successAdded = 0;
    private ArrayList<double[]> susceptibilityData;
    private double hBound = 1.3; private double hMin = 0.95;
    private double minDelH=0.001; private double hOffset = 0.1;
    private double delH =0.025;
    private double currentField=hMin;
    private int currentFieldRun=0;
    private int rIndex=0;
    private Random Ran = new Random();
    private int tNucleationLast=0;
    private String fNameData;
    private int spinFix = -1;
    private boolean hIncrease=true;
    private boolean useSameConfig = false;
    private int currentFixedParam=0;
    private int mStable;

    private SimProcessParser parser; private ParameterBank param;
    private DataSaver dSave; private LatticeMagInt lattice;
    private int tProcess; private int currentTime=0;
    private boolean output =true;
    private String configExplored="random";
    private String susParamPostFix=""; private String susFixedPostFix="";

    /**
    *      MCSusceptibilityFullySimCL constructor
    * 
    * @param fixPost - fixed configuration file post fix to filename 
    * @param susPost - parameter file post fix to filename
    */
    public MCSusceptibilityFullySimCL(String fixPost,String susPost,int nfix,double ho,
            double hd,boolean inc,double bnd){
        // Allow for multiple instances
        dSave = new DataSaver("");
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
        
        susceptibilityData = new ArrayList<double[]>();
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
    public String getParamPostFixBase(){return ("-SusceptibilityFully"+susParamPostFix);}
    public void deleteTempSettingsFile(){
        (new DirAndFileStructure()).deleteTempSettingsFile(getParamPostFix());
    }
    
    @Override
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
        susceptibilityData = new ArrayList<double[]>();
        
        // Change field for new susceptibility measurement
        param.changeField(currentField, getParamPostFix());
        lattice = new SingleValueLattice(param.L,param.s,instId);

        parser = new SimProcessParser(parser,getSimPostFix(),getParamPostFix(),param);
        currentTime =0;
        
        //which algo is going to be used
        if(param.mcalgo.equalsIgnoreCase("metropolisfully")){
            // clear up context on GPU if creating a new context
            if(mc != null){
                ((FullyConnectMetroMCCL)mc).closeOpenCL();
            }
            boolean balanced = false;
            if(configExplored.contains("alanced")){balanced =true;}
            mc = new FullyConnectMetroMCCL(this,output,currentFixedParam,spinFix,
                    balanced,runsPerField);
            ((FullyConnectMetroMCCL)mc).setOutputModeOpenCL(false);
            mc.setTriggerOnOff(false);
            ((FullyConnectMetroMCCL)mc).setMeasuringStartTime(tMinTime);
     
        }

        mStable = param.L*param.s;
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
            initialize();
            initialRun();
            saveResult();
            searching = continueSearch();
        }
        (new DirAndFileStructure()).deleteTempSettingsFile(getParamPostFix());
        
    }


    /**
    *       saveResult saves the resulting susceptibility data to a text
    *   file in the data directory.
    */
    private void saveResult(){	
        setDataFileName(configExplored);
        for(int u = 2; u < susceptibilityData.size();u++){
            saveData(susceptibilityData.get(u));
        }
        susceptibilityData.clear();   
    }

    /**
    *       continueSearch performs the logic to determine if the simulation 
    *   should proceed      
    * 
    * @return true if continue simulation 
    */
    public boolean continueSearch(){
        boolean cont = true;
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
            if(currentField < hBound){cont=true;}else{cont = false;}
        }else{
            if(currentField > hBound){cont=true;}else{cont = false;}
        }
        
        if(successAdded < minAdded){cont = false;}
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
        fNameData =  "SusceptCL/SusceptabilityDataCLGPU-"+mc.getClass().getSimpleName()+
               "-"+fix+"-L-"+param.L+".txt";
    }

    
    /**
    *       initialRun performs the simulation of the system
    */
    public void initialRun(){
        tProcess = parser.nextSimProcessTime();
        int tFlip = parser.timeToFlipField();
        boolean goRun=true;
        boolean badTimeout = false;
        currentTime=0;
        
        while(goRun){
            if(currentTime == tProcess){
                    //System.out.println("Process time");
                    mc = parser.simProcess(mc);
                    tProcess = parser.nextSimProcessTime();
            }
            mc.doOneStep();

            //!=0
            if((currentTime % 1000) == 0 && currentTime  > (tMinTime+1+tFlip)){
                System.out.println("t:"+currentTime);

                ArrayList<MeasureIsingSystem> mNow = ((FullyConnectMetroMCCL)mc).getSystemMeasurerMult();
                for(int u = 0; u < 25;u++){
                    System.out.println("M: "+mNow.get(u).getExpectationM(0) + "    <M>: "+
                            mNow.get(u).getExpectationM(1)+"    <M2>: "+mNow.get(u).getExpectationM(2));
                }
            }
            currentTime++;

            if(currentTime==tMeasure){
                addData(mc.getSusceptibilityMult(),
                        ((FullyConnectMetroMCCL)mc).getFieldMult(),
                        ((FullyConnectMetroMCCL)mc).getNfixedMult(),
                        ((FullyConnectMetroMCCL)mc).getMagMult());
            }

            goRun = !(mc.nucleated());
            
             // Time out after taking data and output this fact
            if(currentTime ==tMeasure){
                goRun = false;
                if(mc.getM()< mStable*0.3){
                System.out.println(" NOT properly timing out with M: "+mc.getM());
                badTimeout = true;
                }else{
                System.out.println("Run appears to properly time out with M: "+mc.getM());}
            }
        }
   
        tNucleationLast = currentTime;
        if(badTimeout){tNucleationLast = 1;} // Dont keep bad data or keep going
    }

    private void addData(ArrayList<Double> susIn, ArrayList<Float> hFieldIn, ArrayList<Integer> nfix, ArrayList<Integer> magIn){
        successAdded = 0;
        for(int u = 0; u < runsPerField; u++){
            if( (magIn.get(u)*param.s) > 0 ){
                successAdded++;
                double[] meta = new double[14];     
                meta[0] = susIn.get(u);
                meta[1] = param.hField;
                meta[2] = hFieldIn.get(u);
                meta[3] = currentFixedParam;
                meta[4] = nfix.get(u);
                meta[5] = param.L;
                meta[6] = param.L-nfix.get(u);
                meta[7] = param.L;
                meta[8] = param.temperature;
                meta[9] = param.jInteraction;
                meta[10] = (param.mcalgo.contains("ully"))? param.L : param.R;
                meta[11] = (configExplored.contains("alanced")) ? 0 : spinFix;
                meta[12] =  (configExplored.contains("alanced")) ? 1 : 0;
                meta[13] =  (configExplored.contains("alanced")) ? 1 : 0;
                susceptibilityData.add(meta);
            }
        }
    }

    // test the class
    public static void main(String[] args) {
        boolean respond=true;
        String responsePost="";
        String responseConfig="";

        System.out.println("==========================================================================");
        System.out.println("Starting example Susceptibility measurement simulation on fully connected system with OpenCL");
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
        
        respond =true; int nfixed = 50;
        double hMin =1.0; double delH = 0.015;
        boolean inc=true; double bnd=1.30;
        
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
            
        MCSusceptibilityFullySimCL mc = new MCSusceptibilityFullySimCL(responseConfig,
                    postfix,nfixed,hMin,delH,inc,bnd);
        ParameterBank param = new ParameterBank(mc.getParamPostFix());
        mc.deleteTempSettingsFile();
        
        int iter = 8; int L = 115;int num = param.L;//L*L;
        delH = 0.0075;
        int perc = (int) (0.025*num);
        for(int i = 1;i<iter;i++){
            nfixed = 0+perc*i;
            double p1 = -1.2701-4.2;//-5.187*10*10*10*10*10; //-3.9187-0.0;
            double p2 = 1.2701;//3.974; //1.27;
            bnd = p2 +delH;
            double p3 = 1.27;//1.247;           
            double x = (double)((double)nfixed/(double)(num));
            hMin =(p1*x + p2)- delH*8;
            //hMin = (p1*Math.pow(x,p2)+p3)- delH*5;
            mc = new MCSusceptibilityFullySimCL(responseConfig,
                    postfix,nfixed,hMin,delH,inc,bnd);
            mc.runSimulation();
        }
        
        System.out.println("______________________________");
        System.out.println("Done!");
        System.out.println("______________________________");
    }
}

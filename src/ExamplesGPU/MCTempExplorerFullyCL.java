package ExamplesGPU;

/**
*    @(#)   MCTempExplorerFullyCL
*/ 

import Examples.*;
import Backbone.System.*;
import Backbone.Util.*;
import Backbone.Algo.*;
import GPUBackend.CLMetropolisMC;
import GPUBackend.FullyConnectMetroMCCL;
import JISim.MCSimulation;
import JISim.MakeFixedConfig;
import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;


/**
*      MCTempExplorerFullyCL class searches for the critical temperature of a system
*  using an OpenCL implementation of metropolis algorithm.
*  . It takes in a fixed config file post script name as input. Example
*  input "cross" to use "fixedcross.txt".
* 
* <br>
* 
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2012-05
*/
public final class MCTempExplorerFullyCL extends MCSimulation {
    private double temperature; private double jInteraction;
    private boolean useLongRange; private boolean useDilution;
    private IsingMC mc; private boolean peakOut=true;
    private int instId;
    private int runsPerTemp = 64;
    private boolean takeSnapshot = true;
    private int currentSeed;
    private Random Ran = new Random();
    private String fNameData;
    
    private SimProcessParser parser; private ParameterBank param;
    private DataSaver dSave; private SimSystem simSys;
    private int tProcess;
    private int currentTime=0;
    private int minT=48000; private int maxT=58000;
    private int tMeasureStart=40000;
    private double delT = 0.025;
    private double minTemp = 4.075;//3.845;//15.25; //2.98269;// tc NN = 2.26918531 R=10 3.98269
    private int fixedSpins = 750;
    private int successAdded = 0;
    private int Ninit;
    private ArrayList<double[]> susceptibilityData;
    private int mStable;
    private boolean useGPU = false;
    private boolean useSameConfig = true;
    private boolean useCancelingField = false;
    private boolean reRun;
    private ArrayList<float[]> reRunData;
    private ArrayList<Float> reRunFields;
    private ArrayList<Integer> reRunFixed;
    private int lastTempInd = 0;
    private double currentTemp=minTemp-delT;
    private double maxTemp = 4.0;
    private double susLast; private double sus;
    private double binder;
    private double delSusLast=0;
    private int spinFix = -1;
    private boolean output =true;
    private String configExplored = "random";

    
    /**
    *      MCTempExplorerFullyCL constructor
    * 
    * @param conf - fixed spin configuration filename 
    */
    public MCTempExplorerFullyCL(String conf,double mT,double mxT, int nfix, int sfix,double del){
        this(conf,mT,mxT,nfix,sfix,del,false);
    }
    
    /**
    *      MCTempExplorerFullyCL constructor
    * 
    * @param conf - fixed spin configuration filename 
    */
    public MCTempExplorerFullyCL(String conf,double mT,double mxT, int nfix, int sfix,double del, boolean run){
        reRun = run;
        minTemp = mT; maxTemp = mxT;
        fixedSpins = nfix; spinFix = sfix;
        delT = del; currentTemp = mT;
        // Allow for multiple instances
        dSave = new DataSaver("");
        instId = dSave.getInstId();
        configExplored =conf;
        // Check for Files 
        (new DirAndFileStructure()).checkForSettingsFiles(this);
        (new DirAndFileStructure()).createTempSettingsFile(getParamPostFixBase(), getParamPostFix());
        
        System.out.println("Instance of program #"+instId);
        System.out.println("______________________________");
        
        param = new ParameterBank(getParamPostFix());
        dSave = new DataSaver(instId,param);
        
        // enforce zero fieldc
        param.changeField(0, getParamPostFix());
        param.changeTemp(minTemp, getParamPostFix());
        
        if(reRun){
            reRunData = reRunDataSetNoDilutionGetData();
            param.changeLength((int)((1.0-((double) nfix)/ ((double) param.L))*param.L), getParamPostFix());
        }
        
        param.printParameters();
        parser = new SimProcessParser(parser,getSimPostFix(),getParamPostFix());
    }
    
    /**
    *      getSimSystem returns the simulations simSys object
    * 
    *   @see  Lattice
    */
    @Override
    public SimSystem getSimSystem(){return simSys;}

    /**
    *      getParamPostFix returns the postfix of the parameters file for this 
    *   simulation.
    */
    public String getParamPostFixBase(){return "-TexplorerFully";}
    
    /**
    *      getParamPostFix returns the postfix of the parameters file for this 
    *   simulation
    */
    @Override
    public String getParamPostFix(){return (getParamPostFixBase()+instId+"-deleteable-");}
    public void deleteTempSettingsFile(){
        (new DirAndFileStructure()).deleteTempSettingsFile(getParamPostFix());
    }

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
    public String getSimPostFix(){return "Explorer";};
    
    
    /**
    *      Close updates any configuration or settings before closing program 
    *  instance.
    */
    @Override
    public void Close(){dSave.updateClosedInstance();
        if(useGPU){((CLMetropolisMC)mc).closeOpenCL();}
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
    public void initialize(double temp,String type){
        susceptibilityData = new ArrayList<double[]>();   
        
        if(useCancelingField){
            double antiField = -4.0*((double)spinFix*((double)fixedSpins)/((double)(param.L*param.L)));
            System.out.println("Field: "+antiField);
            param.changeField(antiField,getParamPostFix());
        }
     
        param.changeTemp(temp,getParamPostFix());
        dSave = new DataSaver(instId,param);

        // Initialize simSys for simulation
        simSys = new SingleValueLattice(param.L,param.s,instId);

        parser = new SimProcessParser(parser,getSimPostFix(),getParamPostFix(),param);
        currentTime = 0;
        
        //which algo is going to be used
        if(param.mcalgo.equalsIgnoreCase("metropolisfully")){
            // clear up context on GPU if creating a new context
            if(mc != null){
                ((FullyConnectMetroMCCL)mc).closeOpenCL();
            }
            boolean balanced = false;
            if(configExplored.contains("alanced")){balanced =true;}
            
            mc = (reRun) ? new FullyConnectMetroMCCL(this,output,fixedSpins,spinFix,
                    balanced,runsPerTemp,reRunFields,reRunFixed)
                    : new FullyConnectMetroMCCL(this,output,fixedSpins,spinFix,
                    balanced,runsPerTemp,null,null);
            
            ((FullyConnectMetroMCCL)mc).setOutputModeOpenCL(false);
            mc.setTriggerOnOff(false);
            ((FullyConnectMetroMCCL)mc).setMeasuringStartTime(minT);        
        }
        
        mc.setMeasuringStartTime(tMeasureStart);
        
        mStable = param.L*param.s;
        output=false;
        // get a random seed and save it
        currentSeed = Math.abs(Ran.nextInt());
        //currentSeed=387483474;
        mc.setSeed(currentSeed);	
    }
    
    /**
    *       runSimulation performs the main logic of the whole simulation by 
    *   calling on one simulation after it initializes it.
    */
    @Override
    public void runSimulation(){
        boolean searching =true;
        while(searching){
            currentTemp = getNewTemp(currentTemp);
            if(reRun){
                currentTemp = reRunData.get(lastTempInd)[0];
                setFieldsAtNextTemp(lastTempInd);
                setNfixedAtNextTemp(lastTempInd);
                runsPerTemp = reRunFields.size();
                System.out.println("MCTempExplorerFullyCL | Rerun of data at "
                        + "temp: "+currentTemp+"   for "+reRunFields.size()
                        +"  data points   minT: "+minTemp+"    maxT: "+maxTemp);
                lastTempInd += reRunFields.size();
            }
            System.out.println("********************************");
            System.out.println("Starting Temp: "+currentTemp+"       "
                    + "FixedSites %: "+((double)fixedSpins*100.0/(double)param.L)+"       ConfigType: "+configExplored);
            System.out.println("********************************");
            
            initialize(currentTemp,configExplored);
            initialRun();
            saveResult();
            searching = continueSearch();
            susLast=sus;
        }
        printResult();
        (new DirAndFileStructure()).deleteTempSettingsFile(getParamPostFix());
    }
    
    private double getNewTemp(double currTemp){
        double newTemp = 0.0; 
        if(!peakOut){
            newTemp = currTemp+delT;
        }else{
            double middle = (maxTemp+minTemp)/2;
            if(currTemp == minTemp){
                newTemp = middle;
            }else if(currTemp == middle){
                newTemp = currTemp+delT;
            }else{
                double diff = Math.abs(currTemp - middle);
                if(currTemp > middle){
                    newTemp = middle - diff;
                }else{
                    newTemp = middle+diff+delT;
                }
            }
        }
        return newTemp;
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
    
    private void addData(ArrayList<Double> susIn, ArrayList<Float> hFieldIn, 
        ArrayList<Integer> nfix, ArrayList<Integer> magIn, ArrayList<MeasureIsingSystem> measurers){
        successAdded = 0;
        for(int u = 0; u < runsPerTemp; u++){
            successAdded++;
            double[] meta = new double[50];     
            
            meta[0] = currentTemp;
            meta[1] = susIn.get(u);
            meta[2] = measurers.get(u).calcBinderRatio();
            meta[3] = measurers.get(u).calcFourthCumulantM();
            meta[4] = measurers.get(u).calcFourthCumulantMsymm();
            meta[5] = measurers.get(u).calcFourthCumulantE();
            meta[6] = measurers.get(u).calcFourthCumulantEsymm();
            meta[7] = hField;
            meta[8] = hFieldIn.get(u);
            meta[9] = param.L;
            meta[10] = (param.mcalgo.contains("ully"))? param.L : param.R;
            meta[11] = param.jInteraction;
            meta[12] = currentSeed;
            meta[13] = (configExplored.contains("alanced")) ? 0 : spinFix;
            meta[14] =  (configExplored.contains("alanced")) ? 1 : 0;
            meta[15] =  (configExplored.contains("alanced")) ? 1 : 0;
            meta[16] = nfix.get(u);
            meta[17] = param.L;
            meta[18] = param.L-nfix.get(u);
            meta[19] = measurers.get(u).getExpectationM(1);
            meta[20] = measurers.get(u).getExpectationM(2);
            meta[21] = measurers.get(u).getExpectationM(3);
            meta[22] = measurers.get(u).getExpectationM(4);
            meta[23] = measurers.get(u).getExpectationM(5);
            meta[24] = measurers.get(u).getExpectationM(6);
            meta[25] = measurers.get(u).getExpectationE(1);
            meta[26] = measurers.get(u).getExpectationE(2);
            meta[27] = measurers.get(u).getExpectationE(3);
            meta[28] = measurers.get(u).getExpectationE(4);
            meta[29] = measurers.get(u).getExpectationE(5);
            meta[30] = measurers.get(u).getExpectationE(6);
            meta[31] = measurers.get(u).getExpectationAbsM(1);
            meta[32] = measurers.get(u).getExpectationAbsM(2);
            meta[33] = measurers.get(u).getExpectationAbsM(3);
            meta[34] = measurers.get(u).getExpectationAbsM(4);
            meta[35] = measurers.get(u).getExpectationAbsM(5);
            meta[36] = measurers.get(u).getExpectationAbsM(6);
            meta[37] = measurers.get(u).getExpectationFreeM(1);
            meta[38] = measurers.get(u).getExpectationFreeM(2);
            meta[39] = measurers.get(u).getExpectationFreeM(3);
            meta[40] = measurers.get(u).getExpectationFreeM(4);
            meta[41] = measurers.get(u).getExpectationFreeM(5);
            meta[42] = measurers.get(u).getExpectationFreeM(6);
            meta[43] = measurers.get(u).getExpectationAbsFreeM(1);
            meta[44] = measurers.get(u).getExpectationAbsFreeM(2);
            meta[45] = measurers.get(u).getExpectationAbsFreeM(3);
            meta[46] = measurers.get(u).getExpectationAbsFreeM(4);
            meta[47] = measurers.get(u).getExpectationAbsFreeM(5);
            meta[48] = measurers.get(u).getExpectationAbsFreeM(6);
            meta[49] = tMeasureStart;
 
            susceptibilityData.add(meta);
        }
    }
    
    
    /**
    *       saveData saves the data that has been prepared to go out to file.
    * @param tcdata - data to be saved 
    */
    public void saveData(double[] tcdata){
        if(reRun){
            setDataFileName("-ReRun-"+configExplored+"-nfix-"+fixedSpins);
        }else{
            setDataFileName(configExplored+"-nfix-"+fixedSpins);
        }
        dSave.saveDoubleArrayData(tcdata,fNameData);
    }

    /**
    *       setDataFileName sets the filename of the outgoing data file.
    * @param fix - filename of data file
    */
    public void setDataFileName(String fix){
        if(useCancelingField){
            fNameData =  "TcData/TcCLDataHField-"+fix+"-sfix-"+spinFix+"-L-"+param.L+".txt";
        }else{
            fNameData =  "TcDataCL/TcDataCL-"+fix+"-sfix-"+spinFix+"-L-"+param.L+".txt";
        }
    }
  
    /**
     *      printResult prints the best range to run simulation for given max time
     *  and min time. It does this print out in the console.
     */
    private void printResult(){
        System.out.println("Field Range that is good in tNucleation range");
        System.out.println("timeMin:"+minT+"      timeMax:"+maxT);
    }
    
    /**
    *       continueSearch performs the logic to determine if the simulation 
    *   should proceed      
    * 
    * @return true if continue simulation 
    */
    public boolean continueSearch(){
        boolean go;
        if((currentTemp+delT) > maxTemp){
            go=false;
        }else{go=true;}
        return go;
     }


    /**
    *       initialRun performs the first simulation of the system. 
    */
    @Override
    public void initialRun(){
        currentTime = 0;
        tProcess = parser.nextSimProcessTime();
        System.out.println();
        boolean goRun=true;
        long timer = 0;
        int tdiff;
        while(goRun){
            // NO FIELD FLIPS OR QUENCHING
            if(currentTime == tProcess){
                mc = parser.simProcess(mc);
                tProcess = parser.nextSimProcessTime();
            }
            //System.out.println("Time: "+currentTime);        
            mc.doOneStep();
            
            if(currentTime == maxT){
                 addData(mc.getSusceptibilityMult(),
                        ((FullyConnectMetroMCCL)mc).getFieldMult(),
                        ((FullyConnectMetroMCCL)mc).getNfixedMult(),
                        ((FullyConnectMetroMCCL)mc).getMagMult(),
                        ((FullyConnectMetroMCCL)mc).getSystemMeasurerMult()
                        );
                goRun = false;
            }

            if(((currentTime % 1000) == 0 || (currentTime < 3 && (currentTime%2==0))) 
                    && currentTime !=0){
                    tdiff = (int)(System.nanoTime()/1000) - (int)timer;
                    System.out.println("t:"+currentTime+"  M:"+mc.getM()+"     Temp: "
                        +mc.getTemperature()+"    hField: "+mc.getHfield()+"    computation time elapsed: "+tdiff);
                    timer = System.nanoTime()/1000;
            }
            currentTime++;
        }

    }

    /**
    *      getParams returns the  parameters for this 
    *   simulation
    */
    @Override
    public ParameterBank getParams(){return param;};
   
    private ArrayList<float[]> reRunDataSetNoDilutionGetData(){
        ArrayList<float[]> DataParsed = new ArrayList<float[]>();
        Scanner scan2;
        int fieldIndex = 8; int fixedIndex = 16;
        int nIndex = 9; int tempIndex = 0;
        float curr = 0.0f; float currN = 0;
        float currNfixed = 0; int balancedInd = 15;
        int sfixInd = 13;
        minTemp = 100000000000.0; maxTemp = 0.0;
        
        setDataFileName(configExplored+"-nfix-"+fixedSpins);
        DirAndFileStructure dir = new DirAndFileStructure();
        try {
            scan2 = new Scanner(new File(dir.getDataDirectory()+fNameData));
            Scanner scanner2;
            while(scan2.hasNextLine() && scan2.hasNextFloat()) {
                scanner2 = new Scanner(scan2.nextLine());        
                float[] dataIn = new float[5];
                for(int u = 0; u < fixedIndex;u++){
                    curr = scanner2.nextFloat();
                    if(u==tempIndex){dataIn[0] = curr;}
                    if(u==fieldIndex){dataIn[1] = curr;}
                    if(u==nIndex){currN = curr;}
                    if(u==balancedInd){dataIn[3] = curr;}
                    if(u==sfixInd){dataIn[4] = curr;}                    
                }
                currNfixed = scanner2.nextFloat();
                dataIn[2] = currN-currNfixed; 
                minTemp = (minTemp > dataIn[0]) ? dataIn[0] :minTemp;
                maxTemp = (maxTemp < dataIn[0]) ? dataIn[0] :maxTemp;
                //System.out.println("Data line parsed :  temp: "+dataIn[0]+"   field: "+dataIn[1]+"    n: "+dataIn[2]);
                DataParsed.add(dataIn);	
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(DataSaver.class.getName()).log(Level.SEVERE, null, ex);
        }
        return DataParsed;
    }
    
    private void setFieldsAtNextTemp(int lastInd){
        float currTemp = reRunData.get(lastInd)[0];
        float[] currData;
        reRunFields = new ArrayList<Float>();
        for(int  u=lastInd; u< reRunData.size();u++){
            currData = reRunData.get(u);
            if(currData[0] != currTemp){
                break;
            }else{
                float hIn = 0.0f;
      
                System.out.println("MCTEMPEXPLORER TEST: HIN:::"+currData[1]);
                reRunFields.add(currData[1]);
            }
        }
    }
    
    private void setNfixedAtNextTemp(int lastInd){
        float currTemp = reRunData.get(lastInd)[0];
        float[] currData;
        reRunFixed = new ArrayList<Integer>();
        for(int  u=lastInd; u< reRunData.size();u++){
            currData = reRunData.get(u);
            if(currData[0] != currTemp){break;}else{
            reRunFixed.add((int)currData[2]);}
        }
    }
    
    
    // test the class
    public static void main(String[] args) {
        boolean ResponseNotDone = true; 
        String response;
        Scanner scan = new Scanner(System.in);

        System.out.println("==========================================================================");
        System.out.println("==========================================================================");
        System.out.println();
        System.out.println("Starting example temperature measurement simulation on fully connected system with OpenCL");
        System.out.println();
        System.out.println("==========================================================================");
        System.out.println("==========================================================================");

        String fp="randombalanced";

        while(ResponseNotDone){

        System.out.println("Default Parameters? (Clear Instances by typing clearInst)");	
        response = scan.next();
        if(response.equalsIgnoreCase("y") || response.equalsIgnoreCase("") || args.length > 0){
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
        
        boolean finiteSize = true;
        boolean changingX = false;
        boolean verifyMode = false;
        
        MCTempExplorerFullyCL sim = new MCTempExplorerFullyCL(fp,0,1,0,1,0.0);
        ParameterBank param = new ParameterBank(sim.getParamPostFix());
        sim.deleteTempSettingsFile();
        int N = param.L;
        int iter = 6;
        int initIter = 1;
        int sfix = -1;
        int nfix = 0;
        double mnT = 0;
        double mxT = 10;
        double del = 12;
        double range = 0.085;
        del = 2.0*range/del;
        double tT;double x;
        double perBase = 0.025;
        double per = perBase*N;

        // percentage to run at and how many loops (iter)
        int u = 1;
        int Nnow = N;
        int Nthresh = 30000;
        if(finiteSize){
            iter = 7;
            Nnow = (int)(20000/4);
            param.changeLength(Nnow, sim.getParamPostFixBase());
            N = Nnow;
        }else if (verifyMode){
            iter = u+3;
            //Nnow = (int)(20000);
        }
        
        // gpu time on 550ti for 64 128 256 512 threads at N=5000
        // 26871681 13850465 20976247 34104231
        if(changingX){
            for(int j = 2; j < iter;j++){
                x = per*j;
                nfix = (int) x;
                tT = 4*0.996*(1-perBase*j);
                mnT = tT-range;
                mxT = tT+range;
                System.out.println("MinTemp: "+mnT+"     MaxTemp: "+mxT);
                sim = new MCTempExplorerFullyCL(fp,mnT,mxT,nfix,sfix,del,true);
                sim.runSimulation();
                sim.Close();
            }
        }
        
        if (verifyMode){
            for(int j = u; j < iter;j++){
                x = per*u;
                nfix = (int) x;
                double scalingVal = 0.995;
                if(u == 1){
                    scalingVal =1-364974.2342*(1/(Math.pow(N, 2.0)*3.89));
                    System.out.println("Scaling for "+(perBase*u)+" of N :"+Nnow+"  nfix: "+nfix+"     is "+scalingVal);
                }
                tT = 4*scalingVal*(1-perBase*u);
                mnT = tT-range;
                mxT = tT+range;
                System.out.println("MinTemp: "+mnT+"     MaxTemp: "+mxT);
                if(j==u+1){
                    sim = new MCTempExplorerFullyCL(fp,mnT,mxT,nfix,sfix,del,false);
                }else if (j==u){
                    param.changeLength((int)(Nnow), sim.getParamPostFixBase());
                    tT = 4*scalingVal;
                    mnT = tT-range;
                    mxT = tT+range;
                    sim = new MCTempExplorerFullyCL(fp,mnT,mxT,nfix,sfix,del,true);
                }else if(j ==u+2){
                    tT = 4*scalingVal*(1-perBase*u);
                    mnT = tT-range;
                    mxT = tT+range;
                    sfix = 0;
                    fp = "random";
                    sim = new MCTempExplorerFullyCL(fp,mnT,mxT,nfix,sfix,del,false);
                }
                sim.runSimulation();
                sim.Close();
            }
        }
        //t:9000  M:100     Temp: Temp: 3.7943000000000004    hField: 0.0    computation time elapsed: 158717005
        if(finiteSize){
            // finite size scaling run
            for(int j = 1; j < iter;j++){
                Nnow = N*j*4;
                if(Nnow > Nthresh){Nnow = (int)(Nnow/2);}
                if(j > initIter){
                    param.changeLength(Nnow, sim.getParamPostFixBase());
                    per = perBase*Nnow;
                    x = per*u;
                    nfix = (int) x;
                    double scalingVal = 0.995;
                    if(u ==2){scalingVal = 0.9985;}
                    if(u == 1){
                        scalingVal =1-364974.2342*(1/(Math.pow(Nnow, 2.0)*3.89));
                        System.out.println("Scaling for "+(perBase*u)+"     is "+scalingVal);
                    }
                    tT = 4*scalingVal*(1-perBase*u);
                    tT = tT - 1.0*range;
                    mnT = tT-range;
                    mxT = tT+range;
                    System.out.println("MinTemp: "+mnT+"     MaxTemp: "+mxT);
                    sim = new MCTempExplorerFullyCL(fp,mnT,mxT,nfix,sfix,del,false);
                    sim.runSimulation();
                    sim.Close();
                }
            }
        }
        System.out.println("Done!!!");
    }
}

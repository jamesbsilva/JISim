package Examples.CriticalPoint;

/**
*    @(#)   MCTempExplorer
*/

import AnalysisAndVideoBackend.StructureFactorSlow;
import AnalysisAndVideoBackend.DisplayMakerCV;
import AnalysisAndVideoBackend.StructureFactorPeaks;
import AnalysisAndVideoBackend.VideoMakerCV;
import AnalysisAndVideoBackend.VideoUpdateThreadCV;
import Backbone.System.*;
import Backbone.Util.*;
import Backbone.Algo.*;
import Backbone.SystemInitializer.InitStripes;
import Backbone.Visualization.TriLatticeDiff2D;
import GPUBackend.CLMetropolisMC;
import GPUBackend.FullyConnectMetroMCCL;
import GPUBackend.StructureFactor2DCL;
import JISim.MCSimulation;
import JISim.MCSimulationIntervention;
import JISim.MakeFixedConfig;
import java.util.Random;
import java.util.Scanner;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
*      MCTempExplorer class searches for the critical temperature of a system
*  . It takes in a fixed config file post script name as input. Example
*  input "cross" to use "fixedcross.txt".
* 
* <br>
* 
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2012-05
*/
public final class MCTempExplorer extends MCSimulation {
    private double temperature; private double jInteraction;
    private boolean useLongRange; private boolean useDilution;
    private IsingMC mc; private int instId;
    private Random Ran = new Random(); private boolean peakOut=true;
    private String fNameData; private String structDataFname;
    private int runsPerTemp = 9; private boolean takeSnapshot = true; private int currentSeed=10;
    
    private SimProcessParser parser; private ParameterBank param;
    private DataSaver dSave; private SimSystem simSys;
    private int tProcess; private int currentTime = 0;
    private int minT = 40490; private int maxT = 40500;
    private int tMeasureStart = 25000;
    private double delT = 0.025;
    private double minTemp = 3.1;//15.25; //2.98269;// tc NN = 2.26918531 R=10 3.98269
    private int fixedSpins = 0; private int spinFix = -1;
    private VideoMakerCV  vid; private int VideoLength = 0;
    private int FramesPerCapture = 500; private int FrameRate = 10; // How often to take a capture 
    private int tMovieStart = 40; // How many pre configuration range mcsteps to record
    private int tMovieEnd = 115170;// How many post configuration range mcsteps to record
    private boolean makeVideo = false; private boolean showVideo = false; 
    private int tMovieInitRunMax = -2; 
    private boolean useGPU = false; private boolean useSameConfig = true;
    private boolean useCancelingField = false; private boolean networkSimulation = false;
    private boolean outputMag = false;
    private double currentTemp = minTemp-delT; private double maxTemp = 4.9;
    private double susLast; private double sus;
    private double binder; private double delSusLast=0;
    private boolean output = true; private String configExplored = "crack";
    private MeasureIsingSystem measurer; private StructureFactor2DCL fftslCL;
    
    /**
    *      MCTempExplorer constructor
    * 
    * @param conf - fixed spin configuration filename 
    */
    public MCTempExplorer(String conf,double mT,double mxT, int nfix, int sfix,double del,boolean peakSearch, double hIn){
        minTemp = mT; currentTemp = mT; maxTemp = mxT;
        fixedSpins = nfix; spinFix = sfix;
        delT = del; peakOut = peakSearch;
        
        // Allow for multiple instances
        dSave = new DataSaver(); instId = dSave.getInstId();
        
        // Check for Files 
        (new DirAndFileStructure()).checkForSettingsFiles(this);
        (new DirAndFileStructure()).createTempSettingsFile(getParamPostFixBase(), getParamPostFix());
        System.out.println("Instance of program #"+instId);
        System.out.println("______________________________");
        param = new ParameterBank(getParamPostFix());
        dSave = new DataSaver(instId,param);
        if(param.mcalgo.contains("ully")){
            param.setN(param.L);
        }
        // enforce zero fieldc
        param.changeField(hIn, getParamPostFix()); param.changeTemp(minTemp, getParamPostFix());
        measurer = new MeasureIsingSystem(temperature,param.N);
        measurer.setN(param.N, fixedSpins);  measurer.setFixedVal(spinFix);
        
        if(param.mcalgo.contains("etwork")){networkSimulation =true;}
        configExplored =conf;
        
        param.printParameters();
        parser = new SimProcessParser(parser,getSimPostFix(),getParamPostFix());
        
        temperature = param.temperature;
        // structure factor calculator
        fftslCL = new StructureFactor2DCL(param.L,
                (param.R == 1)? (float)(2.0*param.L) : (float)(1.0*param.L),param.Geo);
    }
    /**
    *      getSimSystem returns the simulations simSys object
    * 
    *   @see  LatticeInt
    */
    @Override
    public SimSystem getSimSystem(){return simSys;}
    /**
    *      getParamPostFix returns the postfix of the parameters file for this 
    *   simulation.
    */
    public String getParamPostFixBase(){return "-Texplorer";}
    /**
    *      getParamPostFix returns the postfix of the parameters file for this 
    *   simulation
    */
    @Override
    public String getParamPostFix(){return (getParamPostFixBase()+instId+"-deleteable-");}
    /**
    *      deleteTempSettingsFile deletes the temporary settings file
    */
    public void deleteTempSettingsFile(){
        (new DirAndFileStructure()).deleteTempSettingsFile(getParamPostFix());
    }
    /**
    *      getFixedPostFix returns the postfix of the fixed spin configuration
    *   file for this simulation.
    */
    @Override
    public String getFixedPostFix(){return (""+configExplored+"-"+instId);}
    /**
    *      getSimPostFix returns the postfix of the simulation process file for this 
    *   simulation.
    */
    @Override
    public String getSimPostFix(){return "Explorer";}
    
    
    /**
    *      Close updates any configuration or settings before closing program 
    *  instance.
    */
    @Override
    public void Close(){
        dSave.updateClosedInstance();
        if(useGPU){((CLMetropolisMC)mc).closeOpenCL();}
        if(!(fftslCL==null)){fftslCL.forceCloseCL();}
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
        //reset Structure factor calculator
        fftslCL.resetAccumulators();   
        if(useCancelingField){
            double antiField = -4.0*((double)spinFix*((double)fixedSpins)/((double)(param.L*param.L)));
            System.out.println("Field: "+antiField);
            param.changeField(antiField,getParamPostFix());
        }
        currentSeed = Math.abs(Ran.nextInt());
        
        if(fixedSpins==0){
            MakeFixedConfig.clearFixedConfig2D(getFixedPostFix());
        }else if(!param.mcalgo.contains("fully")){
            if(configExplored.contains("nteraction")){   
                MakeFixedConfig.makeFixedConfig2D(
                configExplored,getFixedPostFix(), spinFix,1,1,param.L,fixedSpins,param.R,currentSeed);
            }else{
                MakeFixedConfig.clearFixedConfig2D(getFixedPostFix());
                MakeFixedConfig.makeFixedConfig2D(
                configExplored,getFixedPostFix(), spinFix,(int)(param.L/2),(int)(param.L/2),param.L,fixedSpins,fixedSpins,currentSeed);
            }
        }
        param.changeTemp(temp,getParamPostFix());
        measurer.clearAll();
        measurer.changeTemp(temp);
        dSave = new DataSaver(instId,param);

        // Initialize simSys for simulation
        simSys = null;
        if(networkSimulation){
            //simSys = (Network)(new NetworkBA(param.L*param.L,10,1,getFixedPostFix(),param.useHeter)); 
            simSys = (Network)(new NetworkER(param.L*param.L,0.0025,getFixedPostFix(),param.useHeter));
            //((Network)simSys).saveDegreeFrequency();
        }else{
            if(param.mcalgo.contains("ully")){
                simSys = new SingleValueLattice(param.L,param.s,instId);
            }else{
                simSys  = new AtomicLatticeSumSpin(param.s,getParamPostFix(),getFixedPostFix(),instId);
                //InitStripes initializer = new InitStripes();
                //simSys = initializer.initializeSys((LatticeInt)simSys, param.L,param.R, 30, param.s, param.D);
            }
        } 
        parser = new SimProcessParser(parser,getSimPostFix(),getParamPostFix(),param);
        currentTime = 0;

        //which algo is going to be used
        mc = null;
        if(param.mcalgo.equalsIgnoreCase("metropolis")){
            if(useGPU){
                // clear up context on GPU if creating a new context
                if(mc != null){((CLMetropolisMC)mc).closeOpenCL();}
                mc = new CLMetropolisMC(this,output);                    
            }else{
                mc = new MetropolisMC(this,output);
            }
            System.out.println("Starting Sim Run :"+currentRun+"   with fixed spins: "+fixedSpins);
        }else if (param.mcalgo.equalsIgnoreCase("metropolisfix")){
            mc = new MetropolisFixedStrengthMC(this,output);
            System.out.println("Starting Sim Run :"+currentRun+"   with fixed spins: "+fixedSpins);
        }else if (param.mcalgo.equalsIgnoreCase("metropolisfully")){
            boolean balanced = false;
            if(configExplored.contains("alanced")){balanced =true;}
            if(useGPU){
                mc = new FullyConnectMetroMCCL(this,output,fixedSpins,spinFix,balanced, 256);
            }else{
                mc = new MetropolisFullyConnectedMC(this,output,fixedSpins,param.L,spinFix,balanced);
            }
            System.out.println("Starting Sim Run :"+currentRun+"   with fixed spins: "+fixedSpins);
        }else if(param.mcalgo.equalsIgnoreCase("metropolisnetwork")){
            mc = new MetropolisNetworkMC(this,output);
            System.out.println("Starting Sim Run :"+currentRun+"   with fixed spins: "+fixedSpins);
        }else if(param.mcalgo.equalsIgnoreCase("wolff") && hField==0){
            System.out.println("Creating Wolff algo object.");
            mc = new WolffZeroFieldMC(this,output);
        }
        // turn off trigger
        mc.setTriggerOnOff(false);
        mc.setMeasuringStartTime(tMeasureStart);
        output=false;

        // get a random seed and save it
        currentSeed =Math.abs(Ran.nextInt());
        mc.setSeed(currentSeed);
        
        // Make Video
        if(makeVideo){
            vid = new VideoMakerCV(mc.getSeed(),VideoLength,
                    "-fix"+fixedSpins+"-hfield"+param.hField+"-temp-"+param.temperature
                    +"-J-"+mc.getJinteraction()+"-L-"+param.L+"-R-"+param.R+"-Geo-"+param.Geo);
            vid.setFramesPerCapture(FramesPerCapture);
            vid.setDisplayVideo(showVideo);
        }
        if(makeVideo){(mc.getSimSystem()).makeVideo();}
    }

    /**
    *       runSimulation performs the main logic of the whole simulation by 
    *   calling on one simulation after it initializes it.
    */
    @Override
    public void runSimulation(){
        boolean searching = true;
        while(searching){
            currentTemp = getNewTemp(currentTemp);
            System.out.println("********************************");
            System.out.println("Starting Temp: "+currentTemp+
                    "       FixedSites %: "+((double)fixedSpins*100.0/(double)param.N)
                    +"       fixedSpins: "+fixedSpins);
            System.out.println("********************************");
            for(int i=0;i<runsPerTemp;i++){
                initialize(currentTemp,configExplored);
                measurer.clearAll();
                System.out.println("********************************");
                System.out.println("Run: "+i+"  Temp: "+currentTemp
                        +"       FixedSites %: "
                        +((double)mc.getSimSystem().getNFixed()*100.0/(double)mc.getSimSystem().getN())
                        +"     N: "+mc.getSimSystem().getN());
                System.out.println("********************************");
                initialRun();
                saveResult();
                measurer.clearAll();
            }
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
    *       saveResult  prepares the data to be saved.
    */
    private void saveResult(){	
        double[] meta = new double[49];
        meta[0] = currentTemp;
        meta[1] = mc.getSusceptibility();
        meta[2] = measurer.getSusceptibility();
        meta[3] = mc.getSystemMeasurer().calcBinderRatio();
        meta[4] = mc.getSystemMeasurer().calcFourthCumulantM();
        meta[5] = mc.getSystemMeasurer().calcFourthCumulantMsymm();
        meta[6] = mc.getSystemMeasurer().calcFourthCumulantE();
        meta[7] = mc.getSystemMeasurer().calcFourthCumulantEsymm();
        meta[8] = (param.mcalgo.contains("ully"))? param.L : param.R;
        meta[9] = param.hField;
        meta[10] = param.L;
        meta[11] = currentSeed;
        meta[12] = fixedSpins;
        meta[13] = (configExplored.contains("alanced")) ? 0 : spinFix;
        meta[14] = (configExplored.contains("alanced")) ? 1.0 : 0;
        meta[15] = (configExplored.contains("alanced")) ? 1.0 : 0;
        meta[16] = param.jInteraction;
        meta[17] = mc.getSystemMeasurer().getExpectationM(1);
        meta[18] = mc.getSystemMeasurer().getExpectationM(2);
        meta[19] = mc.getSystemMeasurer().getExpectationM(3);
        meta[20] = mc.getSystemMeasurer().getExpectationM(4);
        meta[21] = mc.getSystemMeasurer().getExpectationM(5);
        meta[22] = mc.getSystemMeasurer().getExpectationM(6);
        meta[23] = mc.getSystemMeasurer().getExpectationE(1);
        meta[24] = mc.getSystemMeasurer().getExpectationE(2);
        meta[25] = mc.getSystemMeasurer().getExpectationE(3);
        meta[26] = mc.getSystemMeasurer().getExpectationE(4);
        meta[27] = mc.getSystemMeasurer().getExpectationE(5);
        meta[28] = mc.getSystemMeasurer().getExpectationE(6);
        meta[29] = mc.getSystemMeasurer().getExpectationAbsM(1);
        meta[30] = mc.getSystemMeasurer().getExpectationAbsM(2);
        meta[31] = mc.getSystemMeasurer().getExpectationAbsM(3);
        meta[32] = mc.getSystemMeasurer().getExpectationAbsM(4);
        meta[33] = mc.getSystemMeasurer().getExpectationAbsM(5);
        meta[34] = mc.getSystemMeasurer().getExpectationAbsM(6);
        meta[35] = mc.getSystemMeasurer().getExpectationFreeM(1);
        meta[36] = mc.getSystemMeasurer().getExpectationFreeM(2);
        meta[37] = mc.getSystemMeasurer().getExpectationFreeM(3);
        meta[38] = mc.getSystemMeasurer().getExpectationFreeM(4);
        meta[39] = mc.getSystemMeasurer().getExpectationFreeM(5);
        meta[40] = mc.getSystemMeasurer().getExpectationFreeM(6);
        meta[41] = mc.getSystemMeasurer().getExpectationAbsFreeM(1);
        meta[42] = mc.getSystemMeasurer().getExpectationAbsFreeM(2);
        meta[43] = mc.getSystemMeasurer().getExpectationAbsFreeM(3);
        meta[44] = mc.getSystemMeasurer().getExpectationAbsFreeM(4);
        meta[45] = mc.getSystemMeasurer().getExpectationAbsFreeM(5);
        meta[46] = mc.getSystemMeasurer().getExpectationAbsFreeM(6);
        meta[47] = mc.getSystemMeasurer().getSusceptibilityFreeM();
        meta[48] = tMeasureStart;      
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
        if(useCancelingField){
            if(useGPU){
                fNameData =  "TcData/TcGPUDataHField"+"-sfix-"+spinFix+"-L-"+param.L
                        +"-R-"+param.R+"-j-"+param.jInteraction+"-Geo-"+param.Geo+".txt";
            }else{                
                fNameData =  "TcData/TcDataHField"+"-sfix-"+spinFix+"-L-"+param.L
                        +"-R-"+param.R+"-j-"+param.jInteraction+"-Geo-"+param.Geo+".txt";
            }
        }else{
            if(useGPU){
                fNameData =  "TcData/TcGPUData"+"-sfix-"+spinFix+"-L-"+param.L
                        +"-R-"+param.R+"-j-"+param.jInteraction+"-Geo-"+param.Geo+".txt";
            }else{                
                fNameData =  "TcData/TcData"+"-sfix-"+spinFix+"-L-"+param.L
                        +"-R-"+param.R+"-j-"+param.jInteraction+"-Geo-"+param.Geo+".txt";
            }
        }
        if(networkSimulation){
            fNameData = "NetworkTc/"+simSys.getClass().getName()+"-"+fNameData;
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
        if((currentTemp+delT) > maxTemp || delT == 0){
            go=false;
        }else{
            go=true;
        }
        return go;
    }
    
    /**
    *       initialRun performs the first simulation of the system. 
    */
    @Override
    public void initialRun(){
        currentTime = 0;
        tProcess = parser.nextSimProcessTime();
        System.out.println("First Action at:"+tProcess);
        boolean goRun=true;
        int tFlip = parser.timeToFlipField();
        ExecutorService pool = Executors.newFixedThreadPool(2);

        //StructureFactor fft = new StructureFactor(param.L,param.L,param.Geo);
        StructureFactorSlow fftsl = new StructureFactorSlow(param.L,param.Geo);
        DisplayMakerCV dispSF = new DisplayMakerCV( "Structure Factor CL");
        //DisplayMakerCV dispSF2 = new DisplayMakerCV("Structure Factor Circular Avg");
        //DisplayMakerCV dispSF2 = new DisplayMakerCV("Structure Factor slow");
        DisplayMakerCV dispConf = new DisplayMakerCV( "Configuration" );
        StructureFactorPeaks pkGet = new StructureFactorPeaks(param.Geo,(int)fftslCL.getLq(),param.R);
        //TriLatticeDiff2D diffVis = new TriLatticeDiff2D(param.L);
        //ConfigDifference diffConf = new ConfigDifference((Lattice)simSys);
        //diffVis.initializeImg();
        
        // dSave.saveLatticeAsImg(tFlip, nRuns, instId, (LatticeInt)mc.getSimSystem());
        long cputime = System.nanoTime();
        while(goRun){
            // NO FIELD FLIPS OR QUENCHING
            if(currentTime == tProcess){
                mc = parser.simProcess(mc);
                tProcess = parser.nextSimProcessTime();
                System.out.println("Next Action at:"+tProcess); 
            }
            //System.out.println("Time: "+currentTime);        
            mc.doOneStep();
            
            //dispConf.addDisplayFrame(mc.getSimSystem().getSystemImg());
            /*
             if(currentTime > 4000 ){
                diffConf.update((LatticeInt)simSys);
                diffVis.spinDraw(diffConf);
                dispConf.addDisplayFrame(diffVis.getImageOfVis());
            }
            */
            
            if(outputMag){
                String fname = "TExplMag/"+configExplored+"-hfield-"
                        +mc.getHfield()+"-x-"+((double)fixedSpins/(double)param.N)
                        +"-J0"+param.jInteraction+"-Geo-"+param.Geo+"/"
                        +"Snapshot-temp-"+currentTemp+"-jInter-"+param.jInteraction
                        +"-L-"+param.L+"-R-"+param.R+"-nfix-"+fixedSpins+"-seed-"+mc.getSeed();
                double[] data = new double[2];
                data[0] = currentTime;
                data[1] = mc.getM();
                dSave.saveDoubleArrayData(data, fname);
            }
            
            if(makeVideo){
                if(currentTime%FrameRate ==0){
                    Callable<Boolean> callable = (Callable<Boolean>)(new VideoUpdateThreadCV(mc.getSimSystem(),vid,true));
                    Future<Boolean> future = pool.submit(callable);
                    try {
                        makeVideo = future.get();
                    } catch (InterruptedException ex) {
                        Logger.getLogger(MCSimulationIntervention.class.getName()).log(Level.SEVERE, null, ex);
                    } catch (ExecutionException ex) {
                        Logger.getLogger(MCSimulationIntervention.class.getName()).log(Level.SEVERE, null, ex);
                    }   
                }
            }
            
            if(currentTime == tMeasureStart){
                String fname = "MeasureStartConfigsDiff/"+configExplored+"-hfield-"
                        +mc.getHfield()+"-x-"+((double)fixedSpins/(double)param.N)
                        +"-J0"+param.jInteraction+"-Geo-"+param.Geo+"/"
                        +"Snapshot-temp-"+currentTemp+"-jInter-"+param.jInteraction
                        +"-t-"+currentTime+"-L-"+param.L+"-R-"+param.R+"-nfix-"+fixedSpins+"-seed-"+mc.getSeed();
                //dSave.saveImage(mc.getSimSystem().getSystemImg(), "png", fname);  
            }
            if(currentTime > (minT+tFlip) && currentTime < (maxT)){
                fftslCL.setConfigurationAndCalcSF(mc.getSimSystem());
                float[] temp = fftslCL.getSF();
                pkGet.setData(temp);
                double[] peakVals = pkGet.getPeakAvgData();
                //dispSF.addDisplayFrame(fftslCL.updateAndGetVisImage());    
                
                //fft.takeFT(mc.getSimSystem());
                //dispSF.addDisplayFrame(fft.getVisImg());
                //fftsl.calc2DSF(mc.getSimSystem());
                //dispSF2.addDisplayFrame(fftsl.getVisImg());
                String outdir = "TestSFCL/Geo-"+param.Geo+"L"+param.L+"R"+param.R+
                        "/dilution-"+((int)(100*fixedSpins/param.N))+"-h-"+param.hField+"-Seed-"+currentSeed;
                dSave.saveDoubleArrayData(peakVals, outdir+"/peakValues.txt");
                
                if(currentTime % 25 ==0){
                    dSave.saveImage(fftslCL.getVisImg(),"png",
                            "sf-geo-"+param.Geo+"-t-"+currentTime+"-seed-"+currentSeed+"-h-"+param.hField+"-temp-"+param.temperature,"Data/"+outdir+"/SFImg/");
                }
                
                if(currentTime == (maxT-1)){
                    String outdir2 = "TestSFCL/Geo-"+param.Geo+"L"+param.L+"R"+param.R+"PhaseData";
                    double[] phase = pkGet.getPhaseDataOnly();
                    phase[2] = ((double)fixedSpins/(double)param.N);
                    phase[3] = param.temperature;
                    phase[4] = param.hField;
                    dSave.saveDoubleArrayData(phase, outdir2+"/phaseState.txt");
                    
                    dSave.saveDoubleArrayListData(fftslCL.getRadVals(), 
                           outdir+"/radVals-geo-"+param.Geo+"-t-"+currentTime+"-seed-"+currentSeed+"-h-"+param.hField+"-temp-"+param.temperature+".txt");
                    dSave.saveDoubleArrayData(fftslCL.getSFcirAvg(), 
                            outdir+"/sfCirAvg-geo-"+param.Geo+"-t-"+currentTime+"-seed-"+currentSeed+"-h-"+param.hField+"-temp-"+param.temperature+".txt"); 
                    dSave.saveDoubleArrayData(fftslCL.getSFradAvg(), 
                            outdir+"/sfRadAvg-geo-"+param.Geo+"-t-"+currentTime+"-seed-"+currentSeed+"-h-"+param.hField+"-temp-"+param.temperature+".txt");
                    dSave.saveDoubleArrayData(fftslCL.getSFavg(), 
                            outdir+"/sfAvg-geo-"+param.Geo+"-t-"+currentTime+"-seed-"+currentSeed+"-h-"+param.hField+"-temp-"+param.temperature+".txt"); 
                    dSave.saveImage(fftslCL.getVisImg(),"png",
                            "sf-geo-"+param.Geo+"-t-"+currentTime+"-seed-"+currentSeed+"-h-"+param.hField+"-temp-"+param.temperature,"Data/"+outdir);
                    int[] peakLoc = pkGet.getPeakLoc();
                    fftslCL.setPeaksLoc(peakLoc);
                    dSave.saveImage(fftslCL.getVisImg(),"png",
                            "sfNonTrackPeaks-geo-"+param.Geo+"-t-"+currentTime+"-seed-"+currentSeed+"-h-"+param.hField+"-temp-"+param.temperature,"Data/"+outdir);
                    dSave.saveImage(mc.getSimSystem().getSystemImg(),"png",
                            "ConfSF-geo-"+param.Geo+"-t-"+currentTime+"-seed-"+currentSeed+"-h-"+param.hField+"-temp-"+param.temperature,"Data/"+outdir);
                }
            }
            
            if(currentTime > (minT+tFlip)){
                measurer.updateM(mc.getM(), true); measurer.updateE(mc.getEnergy(), true);
            }
            if(currentTime > maxT){
                goRun = false;
            }
            if(((currentTime % 50) == 0 || (currentTime < 3 && (currentTime%2==0))) 
                    && currentTime !=0){
                System.out.println("t:"+currentTime+"  M:"+mc.getM()
                        +"   cputime:  "+(System.nanoTime()-cputime)/(1000000) 
                        +"   Mstag: "+mc.getMStag()+"     Temp: "
                        +mc.getTemperature()+"    hField: "+mc.getHfield()
                        +"    jInteraction: "+mc.getJinteraction()+"    L: "+param.L);
                cputime = System.nanoTime();
            }
            currentTime++;
            if(makeVideo && currentTime == tMovieInitRunMax ){vid.writeVideo();makeVideo=false;}
        }
        if(takeSnapshot){
            String f = configExplored+"-seed-"+currentSeed;
            dSave.saveConfigPermanent(mc.getSimSystem(), 0, currentTime, 
                    f,fixedSpins,("CriticalT-L-"+param.L));
        }
        if(makeVideo && VideoLength==0){
            vid.writeVideo();makeVideo=false;
        }
    }

    /**
    *      getParams returns the  parameters for this simulation
    */
    @Override
    public ParameterBank getParams(){return param;};
    
    // test the class
    public static void main(String[] args) {
        boolean ResponseNotDone = true; 
        String response;
        Scanner scan = new Scanner(System.in);
        String fp="random";

        while(ResponseNotDone){
            System.out.println("Default Parameters? (Clear Instances by typing clearInst)");	
            response = scan.next();
            if(response.equalsIgnoreCase("y") || response.equalsIgnoreCase("") || args.length > 0){
                ResponseNotDone=false;
            }else if(response.equalsIgnoreCase("n")){
                ResponseNotDone=false;
                System.out.println("Fixed file post?(e for empty input)");
                response = scan.next();
                fp = (response.equalsIgnoreCase("e")) ? "" : response;
            }else if(response.equalsIgnoreCase("clearInst")){
                DataSaver dsave = new DataSaver();
                dsave.clearInstances();
            }else{
                System.out.println("y or n");
            }
        }
        
        MCTempExplorer sim = new MCTempExplorer(fp,0,1,0,1,0.0,true,0);
        ParameterBank param = new ParameterBank(sim.getParamPostFixBase());
        sim.deleteTempSettingsFile();
        
        int N = param.N;
        int sfix = 0; int nfix = 0;
        double mnT = 0; double mxT = 10;
        double del = 30; double range = 0.0;
        del = 2.0*range/del;
        double tT; double x;
        double maxPer = 1.0; double perBase = 0.03;
        int tempSteps = 30; double tempDel = 0.025;
        int iter = 10;//(int)(maxPer/perBase);
        double per = perBase*N;
        int maxLiter = 2;
        int lnow = 75;
        double hdel = 0.15;
        int uoffset = 0;
        
        for(int ln = 1; ln < maxLiter;ln++){
            for(int tempInd = 0;  tempInd < tempSteps ;tempInd++){
                lnow = 64*(ln+1)+0; //(int)Math.pow(2, 6+ln)+100;
                param.changeLength(lnow, sim.getParamPostFixBase());                           
                N = (int) Math.pow(lnow, param.D);
                per = perBase*N;
                
                uoffset = (int) ((-1.2792*(tempInd*tempDel) + 0.7639)/perBase) - (int)((double)iter/2.0);
                if(uoffset < 0){uoffset = 0;}
                
                //System.out.println("   per: "+per+"    N: "+N);
                for(int u = uoffset; u < (uoffset+iter);u++){
                    for(int fieldInd = 0; fieldInd < 1; fieldInd++){
                        int hInd = fieldInd;
                        double hnew = fieldInd*(1.1-0.72*(0.45-perBase*u));
                        // L=80 R=8 hfield
                        //y = p1*x + p2   p1 = -0.72674p2 = 1.0356
                        for(int k = 0; k < 1;k++){
                            if(k == 0){
                                param.changeJinteraction(-1.0, param.postFix);
                            }else if(k == 1){
                                param.changeJinteraction(1.0, param.postFix);
                            }
                            if(tempInd == 0){   }
                        
                            x = per*u;
                            //x = (0.45*(lnow*lnow))-per*u;
                            nfix = (int) x;
                            if(param.R > 1){
                                if(param.Geo == 3){
                                    tT = (0.35)*0.99;//*(1-perBase*u);
                                }else{
                                    tT = 3.7;//(tempInd*tempDel);//*0.99 ;//*(1-perBase*u);//*(1-perBase*u);
                                }
                            }else{
                                if(param.Geo == 3){
                                    tT = (1.7)*0.99*(1-perBase*u);
                                }else{
                                    tT = (0.000001);//*0.99*(1-perBase*u);
                                }
                            }
                            mnT = tT-range; mxT = tT+range;
                            sim = new MCTempExplorer(fp,mnT,mxT,nfix,sfix,del,true,hnew);
                            sim.runSimulation();
                            sim.Close();
                        }
                }}
            }
        }
        System.out.println("Done!!!");
    }
}

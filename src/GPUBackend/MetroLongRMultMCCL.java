package GPUBackend;
/**
 *   @(#) MetroLongRMultMCCL 
 *
 */

import Backbone.Algo.IsingMC;
import Backbone.System.*;
import Backbone.Util.*;
import JISim.MCSimulation;
import JISim.MakeFixedConfig;
import Triggers.*;
import java.io.File;
import java.io.FileNotFoundException;
import java.nio.IntBuffer;
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;
import scikit.opencl.CLHelper;
import scikit.opencl.RandomNumberCL;

/**
 *  MetroLongRMultMCCL 
 *
 *
 *   <br>
 *
 * @author James B. Silva <jbsilva@bu.edu>
 * @since May 2012
 */
public final class MetroLongRMultMCCL implements IsingMC {
    private CLHelper clhandler;
    private ParameterBank param;
    private String metropolisKernel = "ising2d_longmultmetro";
    private String updateFlKernel = "update_fl_2buffer";
    private LatticeMagInt lattice;
    private ArrayList<Float> jInteractionS;
    private ArrayList<Float> temperatureS;
    private ArrayList<Float> hFieldS;
    private ArrayList<Integer> spinS;
    private ArrayList<Integer> spinSumsS;
    private ArrayList<Integer> Ms;
    private ArrayList<Integer> MstagS;
    private ArrayList<Float> energyS;
    private ArrayList<Integer> fixedLocs;    
    private ArrayList<Integer> nfixedS;    
    private double temperature;
    private double jInteraction;
    private int nfixed = 0;
    private int spinFix = 1;
    private double hField;
    private boolean balancedConfig = false;
    private boolean useLongRange;
    private boolean useDilution;
    private boolean useHeter;
    private Trigger trigger;
    private int tTrigger =0;
    private int tAccumulate =5000; // begin accumulator
    private boolean suppressFieldMessages=false;
    private boolean triggerOn=true;
    private int numSystems = 256;
    private int s0=1;
    private int nThreads=1;
    private int nThreadRuns=1;
    private int maxThreads=101;
    private int nThreadGroup=1;
    private int minDivision =10;
    
    private int N;
    private int R;
    private int run;
    private int instId=0;
    private boolean output = false;
    private boolean triggerReady=false;
    private Random Ran = new Random();
    private ArrayList<MeasureIsingSystem> measurers = new ArrayList<MeasureIsingSystem>();
    private String SimPost="";
    private String ParamPost="";  
    private String FixedPost="";  
    private int magnetization; 
    private int magStaggered;
    private RandomNumberCL RNG;
    private int RNGKey;
    private DirAndFileStructure dir = new DirAndFileStructure();
    private int tFlip;
    public double energy;  
    public int currentSeed;
    public int mcs = 0; // number of MC moves per spin
    public int currentTime;
    
    public MetroLongRMultMCCL(String device){
        param = new ParameterBank("");
        R = param.R;
        hField = param.hField;
        temperature = param.temperature;
        useLongRange = param.useLongRange;
        useDilution = param.useDilution;
        useHeter = param.useHeter;
        N = param.L*param.L;
        System.out.println("MetroLongRMultMCCL |  Evolving "+numSystems+"  systems in OpenCL");
        
        createSystemBuffers();
        
        for(int u = 0; u < numSystems;u++){
            measurers.add(new MeasureIsingSystem(temperature,N));
        }
        
        triggerReady=false;
        triggerOn = false;
                
        clhandler = new CLHelper();
        clhandler.setOutputMode(true);
        clhandler.initializeOpenCL(device);
        for(int u = 0; u < numSystems;u++){
            measurers.add(new MeasureIsingSystem(temperature,N));
        }
        
        maxThreads = (int)clhandler.getCurrentDevice1DMaxWorkItems()/2;  
        //Ran.setSeed(983745);
        // User imposed thread limit
        int maxConcurThreads = 128;
        if(maxConcurThreads>0){maxThreads = maxConcurThreads;}
        
        clhandler.createKernel("", metropolisKernel);
        clhandler.createKernel("", updateFlKernel);

        divideWork();
        
        int nRandom = (2*N*(numSystems+1)+5);
       
        // float parameters
        // jInteraction for systems
        clhandler.createFloatBuffer(metropolisKernel, 0,jInteractionS.size(), jInteractionS,"r",false);
        
        // hField for systems
        clhandler.createFloatBuffer(metropolisKernel, 1,hFieldS.size(), hFieldS,"r",false);
        // copy to updateFL kernel to allow for updating outside of Metropolis kernel
        clhandler.copyFlBufferAcrossKernel(metropolisKernel, 1, updateFlKernel, 0);
        
        // temperatures for system
        clhandler.createFloatBuffer(metropolisKernel, 2,temperatureS.size(), temperatureS,"r",false);
        // copy to updateFL kernel to allow for updating outside of Metropolis kernel
        clhandler.copyFlBufferAcrossKernel(metropolisKernel, 2, updateFlKernel, 1);
        
        // random
        clhandler.createFloatBuffer(metropolisKernel, 3,nRandom,0,"rw");
        
        // energy for systems
        clhandler.createFloatBuffer(metropolisKernel, 4,nRandom,0,"rw");
       
        // error buffer
        clhandler.createFloatBuffer(metropolisKernel, 5,400, 0,"rw");
        
        clhandler.setKernelArg(updateFlKernel);
        
        // integer buffers
        // spin values
        clhandler.createIntBuffer(metropolisKernel,0,spinS.size(), spinS,"rw",false);
        // fixed spin values
        clhandler.createIntBuffer(metropolisKernel,1,fixedLocs.size(), fixedLocs,"r",false);
        // sum of spins in interaction ranges
        clhandler.createIntBuffer(metropolisKernel,2,spinSumsS.size(), spinSumsS,"rw",false);
        // magnetization of systems
        clhandler.createIntBuffer(metropolisKernel,3,Ms.size(), Ms,"rw",false);
        
        // integer arguments
        // Number of systems
        clhandler.setIntArg(metropolisKernel, 0, numSystems);
        clhandler.setIntArg(metropolisKernel, 1, N);
        clhandler.setIntArg(metropolisKernel, 2, param.L);
        
        clhandler.setKernelArg(metropolisKernel);
   
        currentSeed = Ran.nextInt();
        RNG = new RandomNumberCL(currentSeed,clhandler);
        // Get the keyname for the RNG to buffer 4 the randome number buffer
        RNGKey = RNG.addRandomBuffer(metropolisKernel, 4, nRandom);
        
        // Initialize trigger
        setTrigger();
        
        for(int u = 0; u < numSystems;u++){
            measurers.get(u).setN(N, nfixedS.get(u));
        }
        if(output){printParameters();}
        currentTime=0;
    }
    
    public MetroLongRMultMCCL(MCSimulation mc, boolean out,int nfix,int sfix,boolean bal){
        this(mc,out,"GPU",nfix,sfix,bal,256,-1);
    }
    public MetroLongRMultMCCL(MCSimulation mc, boolean out,int nfix,int sfix,boolean bal, int systems){
        this(mc,out,"GPU",nfix,sfix,bal,systems,-1);
    }
    public MetroLongRMultMCCL(MCSimulation mc, boolean out,String device,int nfix,
            int sfix,boolean bal, int systems,int maxConcurThreads){
        param = mc.getParams();
        if(param==null){param = new ParameterBank(mc.getParamPostFix());}
        SimPost = mc.getSimPostFix();
        ParamPost = mc.getParamPostFix();
        FixedPost = mc.getFixedPostFix();
        currentSeed = mc.currentSeed;
        R = param.R;
        numSystems = systems;
        nfixed = nfix;
        spinFix = sfix;
        balancedConfig = bal;
        instId = mc.instId;
        
        // determine if outputting
        output = out;
        if(output){
            System.out.println("MetroLongRMultMCCL | ###########################");
            System.out.println("MetroLongRMultMCCL | Run:"+mc.getRuns()+"     ");
            System.out.println("MetroLongRMultMCCL | ###########################");
        }
        hField = param.hField;
        temperature = param.temperature;
        useLongRange = param.useLongRange;
        useDilution = param.useDilution;
        useHeter = param.useHeter;
        s0 = param.s;
        N = param.N;
        param.setN(N);
        triggerReady = false;
        
        createSystemBuffers();
        
        clhandler = new CLHelper();
        clhandler.setOutputMode(output);
        clhandler.initializeOpenCL(device);
        for(int u = 0; u < numSystems;u++){
            measurers.add(new MeasureIsingSystem(temperature,N));
            measurers.get(u).setFixedVal(spinFix);
        }
        
        maxThreads = (int)clhandler.getCurrentDevice1DMaxWorkItems()/2;  
        //Ran.setSeed(983745);
        // User imposed thread limit
        if(maxConcurThreads > 0){maxThreads = maxConcurThreads;}
        
        clhandler.createKernel("", metropolisKernel);
        clhandler.createKernel("", updateFlKernel);

        divideWork();
        
        int nRandom = (2*N*(numSystems+1)+5);
       
       
        // float parameters
        // jInteraction for systems
        clhandler.createFloatBuffer(metropolisKernel, 0,jInteractionS.size(), jInteractionS,"r",false);
        
        // hField for systems
        clhandler.createFloatBuffer(metropolisKernel, 1,hFieldS.size(), hFieldS,"r",false);
        // copy to updateFL kernel to allow for updating outside of Metropolis kernel
        clhandler.copyFlBufferAcrossKernel(metropolisKernel, 1, updateFlKernel, 0);
        
        // temperatures for system
        clhandler.createFloatBuffer(metropolisKernel, 2,temperatureS.size(), temperatureS,"r",false);
        // copy to updateFL kernel to allow for updating outside of Metropolis kernel
        clhandler.copyFlBufferAcrossKernel(metropolisKernel, 2, updateFlKernel, 1);
        
        // random
        clhandler.createFloatBuffer(metropolisKernel, 3,nRandom,0,"rw");
        
        // energy for systems
        clhandler.createFloatBuffer(metropolisKernel, 4,nRandom,0,"rw");
       
        // error buffer
        clhandler.createFloatBuffer(metropolisKernel, 5,400, 0,"rw");
        
        clhandler.setKernelArg(updateFlKernel);
        
        // integer buffers
        // spin values
        clhandler.createIntBuffer(metropolisKernel,0,spinS.size(), spinS,"rw",false);
        // fixed spin values
        clhandler.createIntBuffer(metropolisKernel,1,fixedLocs.size(), fixedLocs,"r",false);
        // sum of spins in interaction ranges
        clhandler.createIntBuffer(metropolisKernel,2,spinSumsS.size(), spinSumsS,"rw",false);
        // magnetization of systems
        clhandler.createIntBuffer(metropolisKernel,3,Ms.size(), Ms,"rw",false);
        
        // integer arguments
        // Number of systems
        clhandler.setIntArg(metropolisKernel, 0, numSystems);
        clhandler.setIntArg(metropolisKernel, 1, N);
        clhandler.setIntArg(metropolisKernel, 2, param.L);
        
        clhandler.setKernelArg(metropolisKernel);
   
        currentSeed = Ran.nextInt();
        RNG = new RandomNumberCL(currentSeed,clhandler);
        // Get the keyname for the RNG
        RNGKey = RNG.addRandomBuffer(metropolisKernel, 4, nRandom);
        
        System.out.println("MetroLongRMultMCCL |  Filling Lattice Buffer");

        // Initialize trigger
        setTrigger();
        
        for(int u = 0; u < numSystems;u++){
            measurers.get(u).setN(N, nfixedS.get(u));
            measurers.get(u).changeTemp(temperature);
        }
        if(output){printParameters();}
        currentTime=0;
    }
    
    public void printParameters(){
        System.out.println("MetroLongRMultMCCL | ***********************************************");
        System.out.println("MetroLongRMultMCCL | OpenCL System Parameters");
        System.out.println("MetroLongRMultMCCL | Current M: "+magnetization);
        System.out.println("MetroLongRMultMCCL | Device MB Used: "+clhandler.getDeviceUsedMB(metropolisKernel));
        System.out.println("MetroLongRMultMCCL | Max Work Items: "+maxThreads);
        
        System.out.println("MetroLongRMultMCCL | ***********************************************");
        System.out.println("MetroLongRMultMCCL | Work Division Parameters");
        System.out.println("MetroLongRMultMCCL | N Site Visits divided into "+ nThreadRuns
                +" runs of thread groups of "+nThreads+
                "   with local work groups of "+nThreadGroup+"  threads");
        System.out.println("MetroLongRMultMCCL | ***********************************************");
        //printW();
    }
    
    private void divideWork(){   
        nThreads = numSystems;
        boolean divisible = false;
        int div = 2;
        int remainder = 1;
        // If too many threads than divide threads into local workgroups
        while(nThreadGroup> maxThreads && div < nThreads){
            while(!divisible){
                remainder = nThreads % div;
                divisible = (remainder == 0) ? true : false;
                if(!divisible){div++;}
            }
            nThreadGroup = (int)(nThreads/div);
            if(nThreadGroup> maxThreads){div++;}
        }
        nThreadGroup = (int)(nThreads/div);
        System.out.println("MetroLongRMultMCCL | Number of threads total : "+nThreads+"    worked on in groups of "+nThreadGroup);
        
    }
    
    /**
    *         setTrigger should initialize the trigger to be used in the simulation
    */ 
    @Override
    public void setTrigger(){
        trigger = param.setProperTrigger(param, trigger, SimPost, output);
        
        SimProcessParser parser = new SimProcessParser();
        tFlip = parser.timeToFlipField();
    }

    @Override
    public void doOneStep() {
        
        //System.out.println("MetroLongRMultMCCL | Getting Random");
        RNG.fillBufferWithRandom(RNGKey);
        //System.out.println("MetroLongRMultMCCL | Running Metropolis");
        clhandler.runKernel(metropolisKernel, nThreads, nThreadGroup);
        //System.out.println("MetroLongRMultMCCL | Done with Kernel Run");
        
        //clhandler.getBufferIntAsArrayList(metropolisKernel, 0, 5, true);
        //clhandler.getBufferIntAsArrayList(metropolisKernel, 3, 5, true);
        //clhandler.getBufferFlAsArrayList(metropolisKernel, 2, 5, true);
        //error buffer
        //clhandler.getBufferFlAsArrayList(metropolisKernel, 5, 15, true);
        
        currentTime++;
        
        if(triggerOn){triggerReady = trigger.triggerNow(magnetization);}
        if(triggerReady){tTrigger = trigger.getTriggerTime();}

        // Wait till after the flip to measure the susceptibility
        if(currentTime>(tFlip+tAccumulate)){
            updateAccumulators();
        }
    }

    /**
    *         setSeed should set the random number seed to the value given
    * 
    *  @param seed - random number seed
    */ 
    @Override
    public void setSeed(int seed){
        currentSeed = seed;
        Ran.setSeed(seed);
        RNG.setSeed(seed);
    }
   
    /**
    *         resetData resets all the accumulators of data and number of mcsteps made.
    */
    public void resetData() {
        hField = param.hField;
        temperature = param.temperature;
        mcs = 0;
        for(int u = 0; u < numSystems;u++){
            measurers.get(u).clearAll();
        }
    }
    
    /**
    *         resetSimulation should reset all simulation parameters.
    */
    @Override
    public void resetSimulation(){
        currentTime = 0;
        resetData();
        trigger.reset();
        //lattice.initialize(param.s);
        //magnetization = lattice.getMagnetization();
    }
  
   /**
    *         resetSimulation should reset all simulation parameters but 
    *  set seed to the given value.
    * 
    *  @param seed - random number seed
    */ 
    @Override
    public void resetSimulation(int seed){
        currentTime = 0;
        resetData();
        trigger.reset();
        //lattice.initialize(param.s);
        //if(param.useHeter){lattice.setFixedLatticeValues();}
        //magnetization = lattice.getMagnetization();
        setSeed(seed);
    }

   /**
    *         resetSimulation should reset all simulation parameters, set
    *   all lattice values to lattice at the time given, and set the 
    *   seed to the given value
    * 
    *  @param t - time to set lattice to
    *  @param seed - random number seed
    */ 
    @Override
    public void resetSimulation(int t, int seed){
        currentTime = t;
        resetData();
        //change to allow more triggers
        trigger.reset2();//reset but keep std 
        //if(t==0){lattice.initialize(param.s);
        //if(param.useHeter){lattice.setFixedLatticeValues();}
        //}else{
        //lattice.setInitialConfig(t,run,"");}
        //magnetization = lattice.getMagnetization();
        setSeed(seed);
    }

    /**
    *   getSeed should return the current random number seed
    * 
    */
    @Override
    public int getSeed(){return currentSeed;}

    /**
    *       updateAccumulators updates accumulators used to calculate susceptibility
    */
    private void updateAccumulators(){
        ArrayList<Float> energySys = clhandler.getBufferFlAsArrayList(metropolisKernel, 4,numSystems, false);
        ArrayList<Integer> magSys = clhandler.getBufferIntAsArrayList(metropolisKernel, 3,numSystems, false);
        for(int u = 0; u < numSystems;u++){
            (measurers.get(u)).updateE((double)energySys.get(u));
            (measurers.get(u)).updateM(magSys.get(u));         
        }
        mcs++;
    }
    /**
    *   getSpecificHeat should return the current specific heat.
    * 
    */
    @Override
    public double getSpecificHeat() {
        return(measurers.get(0).getSpecificHeat());
    }
    /**
    *   getSusceptibility should return the current susceptibility
    * 
    */
    @Override
    public double getSusceptibility() {
        return(measurers.get(0).getSusceptibility());
    }
    /**
    *   getSusceptibilityMult should return the current susceptibility multiple values if possible
    * 
    */
    @Override
    public ArrayList<Double> getSusceptibilityMult() {
        ArrayList<Double> sus = new ArrayList<Double>();
        for(int u = 0; u < numSystems;u++){
            sus.add(measurers.get(u).getSusceptibility());
        }
        return (sus);
    }
    
    public ArrayList<Double> getSpecificHeatMult() {
        ArrayList<Double> s = new ArrayList<Double>();
        for(int u = 0; u < numSystems;u++){
            s.add(measurers.get(u).getSpecificHeat());
        }
        return (s);
    }
    
    public ArrayList<Float> getFieldMult() {
        return (hFieldS);
    }

    public ArrayList<Integer> getNfixedMult() {
        return nfixedS;
    }
    
    public ArrayList<Integer> getMagMult() {
        ArrayList<Integer> magSys = clhandler.getBufferIntAsArrayList(metropolisKernel, 3,numSystems, false);
        return (magSys);
    }
    
    /**
    *   flipField should flip the magnetic field
    * 
    */
    @Override
    public void flipField(){
        updateField( -1.0*hField, hField);
        hField = -hField;
        if(!suppressFieldMessages){System.out.println("MetroLongRMultMCCL | Field Flipped at t:"+currentTime+" now hField:"+hField);}
    }
    /**
    *   alignField should align the magnetic field with current state
    * 
    */
    @Override
    public void alignField(){
        if(Math.signum(hField) != Math.signum(magnetization)){flipField();}
    }
    /**
    *         changeT should set the temperature to the given value.
    * 
    *  @param temp - temperature to set simulation to
    */ 
    @Override
    public void changeT(double temp){
        updateTemp(temp, temperature);
        updateMeasurerTemp(temp);
        temperature = temp;
    }

    /**
    *         changetH should set the magnetic field to the given value.
    * 
    *  @param hnow - new magnetic field
    */ 
    @Override
    public void changeH(double hnow){
        updateField(hnow, hField);
        hField = hnow;
    }

    /**
    *         changeTandH should set the temperature to the given value and
    *   the magnetic field.
    * 
    * 
    *  @param temp - new temperature
    *  @param hnow - new field
    */ 
    @Override
    public void changeTandH(double temp,double hnow){
        updateTempField(temp, temperature, hnow, hField);
        updateMeasurerTemp(temp);
        hField = hnow;
        temperature = temp;
    }
    
    private void updateMeasurerTemp(double newTemp){
        for(int u =0;u<measurers.size();u++){
            measurers.get(u).changeTemp(newTemp);
        }
    }
    
    /**
    *         changeTFlipField should flip the field and set the temperature to
    *   the given value
    * 
    *  @param temp - new temperature
    */ 
    @Override
    public void changeTFlipField(double temp){
        updateTempField(temp, temperature, -1.0*hField, hField);
        updateMeasurerTemp(temp);
        temperature = temp;
        hField = -hField; 
    }

    /**
    *   nucleated should return true if nucleation has occurred
    * 
    */
    @Override
    public boolean nucleated(){
        return triggerReady;
    }

    /**
    *   getSimSystem should return the lattice object
    * 
    */
    @Override
    public LatticeMagInt getSimSystem(){return lattice;}
    /**
    *   getM should return the current magnetization.
    * 
    */
    @Override
    public int getM(){
        int averagingNum = 10;
        ArrayList<Integer> magSys = clhandler.getBufferIntAsArrayList(metropolisKernel, 3,averagingNum, false);
        long magsum = 0;
        for(int u = 0; u < averagingNum;u++){
            //System.out.println("System "+u+"    m: "+magSys.get(u));
            magsum += magSys.get(u);
        }  
        magnetization = (int)(magsum/averagingNum);
        /*System.out.println("FullyConnectMetrMCCL | hField");
        clhandler.getBufferFlAsArrayList(metropolisKernel, 1, 5, true);
        System.out.println("FullyConnectMetrMCCL | Temp");
        clhandler.getBufferFlAsArrayList(metropolisKernel,  2, 5, true);
        System.out.println("FullyConnectMetrMCCL | Random");
        clhandler.getBufferFlAsArrayList(metropolisKernel, 4, 5, true);
        */
        return magnetization;
    }
    /**
    *   getMStag should return the current regular magnetization for fully connected model.
    * 
    */
    public int getMStag(){return getM();}
    /**
    *   getEnergy should return the current energy.
    * 
    */
    @Override
    public double getEnergy(){
        return energy;
    }
    /**
    *   getRun should return the current run number
    * 
    */
    @Override
    public int getRun(){return run;}
    /**
    *         setRun should set the number of the current run to the value given.
    *
    *  @param cr - current run value 
    */
    @Override
    public void setRun(int cr){run = cr;}
    /**
    *   getTriggerTime should return the time the trigger went off.
    * 
    */
    @Override
    public int getTriggerTime(){return tTrigger;}
    /**
    *         setTriggerOnOff allows turning triggers off
    * 
    * @param tr - on is true
    */ 
    @Override
    public void setTriggerOnOff(boolean tr){
        triggerOn=tr;
    }
    
    /**
    *   getJinteraction should return the current strength of interaction.
    * 
    */
    @Override
    public double getJinteraction(){
        return jInteraction;
    }
    
    /**
    *   getHfield should return the magnetic field value.
    * 
    */
    @Override
    public double getHfield(){return hField;}    
    /**
    *   getTemperature should return the temperature value.
    * 
    */
    @Override
    public double getTemperature(){return temperature;}
    /**
    *         getConfigRange should save the configurations ins the given range
    * 
    *  @param tInitial - initial time in range
    *  @param tFinal - final time in range
    */ 
    @Override
    public void getConfigRange(int tInitial, int tFinal){
        //suppressFieldMessages = true;

        resetSimulation(currentSeed);
        SimProcessParser parser = new SimProcessParser();
        DataSaver save = new DataSaver(lattice.getInstId());
        int tProcess = parser.nextSimProcessTime();


        //assert field is right
        if(hField == param.hField){}else{flipField();}

        for(int t = 0; t<tFinal;t++){
        //Need to flip at step tflip not tflip+1
        if(currentTime == tProcess){
            parser.simProcess(this);
            tProcess = parser.nextSimProcessTime();
        }
        doOneStep();
        if(t>tInitial){save.saveConfig(lattice,run,t);}
        }
    }

    /**
    *       getSystemMeasurer should return the system measurer which contains system data
    *   and accumulators.
    * 
    */
    @Override
    public MeasureIsingSystem getSystemMeasurer(){
        return measurers.get(0);
    }
    
    public ArrayList<MeasureIsingSystem> getSystemMeasurerMult(){
        return measurers;
    }
    
    private void updateTempField(double temp,double tempOld,double field,double fieldOld){
        float diff = (float)(field-fieldOld);
        float difft = (float)(temp-tempOld);
        for(int u = 0; u < numSystems; u++){
            hFieldS.set(u, hFieldS.get(u)+diff);
            temperatureS.set(u, temperatureS.get(u)+difft);
        }
        clhandler.setFloatBuffer(metropolisKernel, 1, hFieldS.size(),hFieldS,true);
        clhandler.setFloatBuffer(metropolisKernel, 2, temperatureS.size(),temperatureS,true);
        clhandler.runKernel(updateFlKernel, hFieldS.size(), 1);
        //System.out.println("New Fields");
        //clhandler.getBufferFlAsArrayList(metropolisKernel, 1, 10, true);
        //System.out.println("New Temp");
        //clhandler.getBufferFlAsArrayList(metropolisKernel,2, 10, true);
    }
    
    private void updateField(double field, double fieldOld){
        float diff = (float)(field-fieldOld);
        for(int u = 0; u < numSystems; u++){
            hFieldS.set(u, hFieldS.get(u)+diff);
        }
        clhandler.setFloatBuffer(metropolisKernel, 1, hFieldS.size(),hFieldS,true);
        clhandler.runKernel(updateFlKernel, hFieldS.size(), 1);
        //System.out.println("New Fields");
        //clhandler.getBufferFlAsArrayList(metropolisKernel,  1, 10, true);
        
    }

    private void updateTemp(double temp,double tempOld){
        float diff = (float)(temp-tempOld);
        for(int u = 0; u < numSystems; u++){
            temperatureS.set(u, temperatureS.get(u)+diff);
        }
        clhandler.setFloatBuffer(metropolisKernel, 2, temperatureS.size(),temperatureS,true);
        clhandler.runKernel(updateFlKernel, temperatureS.size(), 1);
        //System.out.println("New Temp");
        //clhandler.getBufferFlAsArrayList(metropolisKernel, 2, 10, true);
    }
    
    
    public void createSystemBuffers(){
        for(int u = 0; u < numSystems; u++){
            insertAsystem(nfixed, (balancedConfig) ? "randombalanced": "random", FixedPost, u, 0);
        }
    }
    
    public void insertAsystem( int fixedSpins,String configExplored, String fixPost, int sysNum,int t){
        int[][] fix =  null;
        if(fixedSpins==0){
                MakeFixedConfig.clearFixedConfig2D(fixPost);
        }else if(!param.mcalgo.contains("fully")){
            if(configExplored.contains("nteraction")){   
                fix = MakeFixedConfig.makeFixedConfig2D(
                configExplored,fixPost, spinFix,1,1,param.L,fixedSpins,param.R,currentSeed);
            }else{
                fix = MakeFixedConfig.makeFixedConfig2D(
                configExplored,fixPost, spinFix,(int)(param.L/2),(int)(param.L/2),param.L,fixedSpins,fixedSpins);
            }
        }
        ArrayList<Integer> fixIn = new ArrayList<Integer>();
        ArrayList<Integer> spinsIn  = setInitialConfig(t, run, fixPost);
        spinsIn = adjustForFixed(spinsIn, fix);
        fixIn = adjustForFixed(fixIn, fix);
        int sysStart = N*sysNum;
        ArrayList<Integer> sumIn  = calcSum(spinsIn);
        boolean update = ((sysStart+1) > spinS.size()) ? false : true;
        for(int u = 0; u < param.N;u++){
            if(update){
                spinS.set(sysStart+u,spinsIn.get(u));
                fixedLocs.set(sysStart+u,fixIn.get(u));
                spinSumsS.set(sysStart+u,sumIn.get(u));
            }else{
                spinS.add(spinsIn.get(u));
                fixedLocs.add(fixIn.get(u));
                spinSumsS.add(sysStart+u,sumIn.get(u));
            }
        }
        if(update){
            hFieldS.set(sysNum,(float)(hField));
            Ms.set(sysNum,Ms.set(sysNum, magnetization+spinsIn.get(N)));
            MstagS.set(sysNum, magStaggered+spinsIn.get(N+1));
            jInteractionS.set(sysNum, (float)((param.jInteraction < 0) ? -4.0/(Math.pow((2*R+1),2.0)-1) : 4.0/(Math.pow((2*R+1),2.0)-1)));
            temperatureS.set(sysNum, (float)(temperature));
            energyS.set(sysNum,(float)calcIsingEnergy(spinsIn,
                    (param.jInteraction < 0) ? -4.0/(Math.pow((2*R+1),2.0)-1) : 4.0/(Math.pow((2*R+1),2.0)-1),sumIn));
        }else{
            hFieldS.add(sysNum,(float)(hField));
            Ms.add(sysNum, magnetization+spinsIn.get(N));
            MstagS.add(sysNum, magStaggered+spinsIn.get(N+1));
            jInteractionS.add( (float)((param.jInteraction < 0) ? -4.0/(Math.pow((2*R+1),2.0)-1) : 4.0/(Math.pow((2*R+1),2.0)-1)));
            temperatureS.add(sysNum, (float)(temperature));
            energyS.add(sysNum,(float)calcIsingEnergy(spinsIn,
                (param.jInteraction < 0) ? -4.0/(Math.pow((2*R+1),2.0)-1) : 4.0/(Math.pow((2*R+1),2.0)-1),sumIn));
        }
        
    }
    
    public ArrayList<Integer> adjustForFixed(ArrayList<Integer> spins,int[][] fix){
        int x; int y; int s;
        int sum = 0;int sumStag = 0;
        for(int u = 0; u < fix.length;u++){
            x = fix[u][0];
            y = fix[u][1];
            s = fix[u][2];
            if(spins.get(+y*param.L) == 0){
                sum += (spins.get(x+y*param.L) != s) ? s : 0;
                sumStag += ((x+y*param.L)%2 == 0 && spins.get(x+y*param.L) != s) ? s : -1*s;
            }else{
                sum += (spins.get(x+y*param.L) != s) ? 2*s : 0;
                sumStag += ((x+y*param.L)%2 == 0 && spins.get(x+y*param.L) != s) ? 2*s : -2*s;
            }
            spins.set(x+y*param.L,s);
        }
        spins.add(sum);
        spins.add(sumStag);
        return spins;
    }
    
    
    public ArrayList<Integer> adjustFixed(ArrayList<Integer> fixList,int[][] fix){
        int x; int y; int s;
        for(int u = 0; u < fix.length;u++){
            x = fix[u][0];
            y = fix[u][1];
            fixList.set(x+y*param.L,1);
        }
        return fixList;
    }
    
    
    
    public ArrayList<Integer> calcSum(ArrayList<Integer> spins){
        ArrayList<Integer> sums = new ArrayList<Integer>();
        for(int u = 0; u < N; u++){
            sums.add(calcSumNeig(spins, u,param.L,param.D));
        }        
        return sums;
    }    
    
    private int calcSumNeig(ArrayList<Integer> spins, int i,  int L, int D){
        int j = ((int)i/L)%L;
        int k = ((int)(i/(L*L)))%L;
        i = i%L;
        int u;
        int v;
        int z;
        int sum = 0;
        for (int m = 0; m < (2 * R + 1); m++) {
            for (int n = 0; n < (2 * R + 1); n++) {
                u = ((i - R + m + L) % L);
                v = ((j - R + n + L) % L);
                if(D ==3){
                    for (int p = 0; p < (2 * R + 1); p++) {
                        z = ((k - R + p + L) % L);
                        if (!(u == i && v == j && z == k)) {
                            sum += spins.get(u+v*L+z*L*L);
                        }
                    }
                }else{
                    if (!(u == i && v == j)) {
                        sum += spins.get(u+v*L);
                    }
                }
            }
        }
        return sum;
    }
    
    /**
    *       calcEnergy calculates the energy of the lattice in an Ising Model 
    * 
    * @return energy in lattice for Ising model
    */
    public double calcIsingEnergy(ArrayList<Integer> lat, double jInt , ArrayList<Integer> sumIn){
        double e=0.0;
        
        for (int i = 0; i < param.L; i++)
        for (int j = 0; j < param.L; j++){{
            if(param.D==3){
                for (int k = 0; k < param.L; k++){			    		  
                    e += -1.0*jInt*lat.get(i+j*param.L+k*param.L*param.L)*sumIn.get(i+j*param.L+k*param.L*param.L)/2
                        -1.0*hField*lat.get(i+j*param.L+k*param.L*param.L);   
                }
            // 2d case    
            }else{
                    // divide the bonds term by 2 since they will be double counted
                    e += -1.0*jInt*lat.get(i+j*param.L)*sumIn.get(i+j*param.L)/2
                        -1.0*hField*lat.get(i+j*param.L);
            }
        }}
        return e;
    }


    /**
    *         setInitialConfig should set all spins in the lattice to the values
    *   given by a file in the Config directory.
    * 
    *   @param t - time of config to set lattice to
    *   @param run - run of time to set lattice to
    *   @param  post - postfix or filename of lattice file
    */
    public synchronized ArrayList<Integer> setInitialConfig(int t, int run, String post) {	 
        int L  = param.L;
        ArrayList<Integer> spins  = new ArrayList<Integer>();
        // If Heterogenous- Fix some spins in the fixed lattice and update the magnetization
        int sum=0; int sumStag=0;   
        if(t==0){
            if(param.D==3){
                for(int x=0;x<L;x++){
                for(int y=0;y<L;y++){
                for(int z=0;z<L;z++){   
                    spins.add(s0);
                    sum+=s0;
                    sumStag += ((x+y*L+z*L*L)%2==0) ? s0: -1*s0; 
                }}}
            }else{
                for(int x=0;x<L;x++){
                for(int y=0;y<L;y++){   
                    int z=0; 
                    spins.add(s0);
                    sum+=s0;
                    sumStag += ((x+y*L+z*L*L)%2==0) ? s0: -1*s0;
                }}
            }
        }else{
            
            // Open file and set spins to values in file
            try {
                // using data saver to get file location
                String fname = dir.getLatticeConfig(instId, t, run, post);
                Scanner scanner = new Scanner(new File(fname));
                sum = 0;
                while(scanner.hasNextInt()) {
                    int fSpin;
                    int x = scanner.nextInt();
                    int y = scanner.nextInt();
                    int z = 0;
                    if (param.D==2){
                        fSpin = scanner.nextInt();
                    }else{				
                        z = scanner.nextInt();
                        fSpin = scanner.nextInt();
                    }
                    sum = fSpin + sum;
                    sumStag += ((x+y*L+z*L*L)%2==0) ? fSpin: -1*fSpin;
                    spins.add(fSpin);
                }
            //System.out.println("value " + s)
            } catch (FileNotFoundException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
        magnetization = sum;
        magStaggered = sumStag;
        return spins;
    }
    
    
    /**
    *       calcEnergy calculates the energy of the lattice in an Ising Model 
    * 
    * @return energy in lattice for Ising model
    */
    private float calcIsingEnergy(ArrayList<Integer> sys){
        double e=0.0;
        
        return ((float)e);
    }
    
    /**
    *         setMeasuringStartTime sets the time to begin all measurement of all variables
    *   by a system measurer.
    * 
    *  @param tin - initial measurement time of all data
    */ 
    @Override
    public void setMeasuringStartTime(int tin){
        tAccumulate = tin;
    }
    
    
    
    public void closeOpenCL(){
        clhandler.closeOpenCL();
    }
    
    public void setOutputModeOpenCL(boolean md){
        clhandler.setOutputMode(md);
    }
    
    // test the class
    public static void main(String[] args) {
        MetroLongRMultMCCL mc = new MetroLongRMultMCCL("GPU");
        long time;
        
        mc.changeH(0.005);
        mc.doOneStep();
        mc.changeH(0.003);
        mc.doOneStep();
        mc.changeT(25.008);
        mc.doOneStep();
        
        System.out.println("||||||  Done with parameter changes.");
        
        for(int j = 0; j< 50;j++){
            time = System.nanoTime();
            mc.doOneStep();
            time = System.nanoTime()-time;
            System.out.println("MetroLongRMultMCCL | M: "+ mc.getM());
            System.out.println("MetroLongRMultMCCL | time/system (us) for step: "+(time/(1000*mc.numSystems)));
        }

        mc.changeH(-0.005);
        
        for(int j = 0; j< 400;j++){
            time = System.nanoTime();
            mc.doOneStep();
            time = System.nanoTime()-time;
            System.out.println("MetroLongRMultMCCL | M: "+mc.getM());
            System.out.println("MetroLongRMultMCCL | time/system (us) for step: "+(time/(1000*mc.numSystems)));
        }
        
        mc.closeOpenCL();
        System.out.println("MetroLongRMultMCCL | DONE!");
    }
    
}
package GPUBackend;

/**
 *   @(#) CLMetropolisMC 
 *
 */

import Backbone.Util.MeasureIsingSystem;
import Backbone.System.LatticeMagInt;
import Backbone.Algo.IsingMC;
import Backbone.Util.DataSaver;
import Backbone.Util.SimProcessParser;
import Backbone.Util.ParameterBank;
import JISim.MCSimulation;
import Triggers.*;
import java.util.ArrayList;
import java.util.Random;
import scikit.graphics.dim2.GridCL;
import scikit.jobs.SimulationCL;
import scikit.opencl.CLHelper;
import scikit.opencl.RandomNumberCL;

/**
 *      CLMetropolisMC is an interface to the OpenCL kernels for OpenCL kernels
 *  implementing metropolis algorithms for ising model. 
 *          
 *
 *   <br>
 *
 * @author James B. Silva <jbsilva@bu.edu>
 * @since May 2012
 */
public final class CLMetropolisMC implements IsingMC {
    private CLHelper clhandler;
    private ArrayList<Float> w;
    private ArrayList<Float> sysParams;
    private ParameterBank param;
    private String metropolisBoostKernel = "ising2d_longmetro";
    private String metropolisKernel = metropolisBoostKernel;//"ising2d_longmetro";
    private String metropolisNNKernel = "ising2d_nnmetro";
    private String updateFlKernel = "update_fl_buffer";
    private LatticeMagInt lattice;
    private double temperature;
    private double jInteraction;
    private double hField;
    private boolean useLongRange;
    private boolean useDilution;
    private boolean useHeter;
    private Trigger trigger;
    private int tTrigger =0;
    private int tAccumulate =1000; // begin accumulator
    private boolean suppressFieldMessages=false;
    private boolean triggerOn=true;
    private int Q;
    private int L;
    private int subL=L;
    private int nThreads=1;
    private int nThreadRuns=1;
    private int maxThreads=101;
    private int nThreadGroup=1;
    private int minDivision =10;
    private GridCL grid = null;
    private int N;
    private int D;
    private int R;
    private int run;
    private boolean output = true;
    private boolean triggerReady=false;
    private Random Ran = new Random();
    private MeasureIsingSystem measurer;
    private String SimPost;
    private String ParamPost;  
    private int magnetization;  
    private int magStaggered = 0;
    private RandomNumberCL RNG;
    private int RNGKey;
    private int buffNumSysParams = 0;
    
    private int tFlip;
    public double energy;  
    public int currentSeed;
    public int mcs = 0; // number of MC moves per spin
    public int currentTime;
    
    /**
    *      CLMetropolisMC constructor.
    * 
    * @param device - type of device to run simulation on ("GPU/CPU"
    */
    public CLMetropolisMC(String device){
        param = new ParameterBank("");
        w = new ArrayList<Float>();
        sysParams = new ArrayList<Float>();
        jInteraction = param.jInteraction;
        hField = param.hField;
        temperature = param.temperature;
        useLongRange = param.useLongRange;
        useDilution = param.useDilution;
        useHeter = param.useHeter;
        L =  param.L;
        D = param.D;
        R = param.R;
        
        N = (int)Math.pow(L,D);
        measurer = new MeasureIsingSystem(temperature,N);
        triggerReady=false;
        triggerOn = false;
        clhandler = new CLHelper();
        
        grid = new GridCL("CLMetropolisMC | ",clhandler,param.L,param.L,param.Geo);
        if(grid != null){
            grid.getCLScheduler().waitForInitiatedCLGL();
            clhandler = grid.getCLHelper();
        }else{
            clhandler.initializeOpenCL(device);
        }
        maxThreads = (int)clhandler.getCurrentDevice1DMaxWorkItems()/2;
        //Ran.setSeed(983745);
        
        if(!useLongRange){R=1;metropolisKernel = metropolisNNKernel;}
        clhandler.createKernel("", metropolisKernel);
        clhandler.createKernel("", updateFlKernel);

        divideWork();
        initSysParams();
        
        int nRandom = (nThreads*(2*nThreadRuns)+1);
        
        // sysParams
        clhandler.createFloatBuffer(metropolisKernel, 0,sysParams.size(), sysParams,"r",true);
        clhandler.copyFlBufferAcrossKernel(metropolisKernel, 0, updateFlKernel, buffNumSysParams);
        
        clhandler.setKernelArg(updateFlKernel);
        
        //Random
        clhandler.createFloatBuffer(metropolisKernel, 1,nRandom, 0,"rw");
        //dE
        clhandler.createFloatBuffer(metropolisKernel, 2,nThreads, 0,"rw");
        
        //error buffer
        clhandler.createFloatBuffer(metropolisKernel, 3,400, 0,"rw");
        
        // Spin system and fixed buffer
        clhandler.createIntBuffer(metropolisKernel, 0,N, 1,"rw");
        clhandler.createIntBuffer(metropolisKernel, 1,N, 0,"r");
        // dM
        clhandler.createIntBuffer(metropolisKernel, 2,nThreads, 0,"rw");
        
        if(metropolisKernel.contains("boost")){
            ArrayList<Integer> sumNeigh = new ArrayList<Integer>();
            for(int u = 0;u < N;u++){
                sumNeigh.add((2*param.R+1)*(2*param.R+1)-1);
            }
            clhandler.createIntBuffer(metropolisKernel, 3,N, sumNeigh,"rw",false);
        }
        
        // 
        clhandler.setIntArg(metropolisKernel, 0, nThreads);
        if(useLongRange){      
            clhandler.setIntArg(metropolisKernel, 1, nThreadRuns);
        }
   
        if(grid != null){
            // Coloring kernel
            String colorKernel = "colorLatticeInt";
            clhandler.createKernel("", colorKernel);
            clhandler.copyIntBufferAcrossKernel(metropolisKernel, 0, colorKernel, 0);
            grid.getCLScheduler().schedCopyColSharedBuffer(colorKernel, 0);
            System.out.println("Number of Threads: "+N);
            grid.getCLScheduler().setColorKernel1D(colorKernel, N, nThreadGroup);
            SciKitSimInterface scikit = new SciKitSimInterface(this);
            System.out.println("SCIKIT INTERFACE | "+scikit);
            grid.getCLScheduler().setSimulation((SimulationCL)scikit);
            grid.getCLScheduler().setKernelChunkRunsOff(true);
        }
        
        currentSeed = Ran.nextInt();
        RNG = new RandomNumberCL(currentSeed,clhandler);
        RNGKey = RNG.addRandomBuffer(metropolisKernel, 1, nRandom);
       
        lattice =  new LatticeCL(param.s,"",
                "",0,clhandler,metropolisKernel,0,1);
        
        magnetization = ((LatticeCL)lattice).calcCorrectionMag(param.s);
        energy = ((LatticeCL)lattice).calcIsingEnergy();
        printParameters();
        measurer.setN(N, lattice.getNFixed());
        currentTime=0;
    }
    
    /**
    *      CLMetropolisMC constructor.
    * 
    * @param mc - monte carlo simulation to be updated based on this algo
    * @param out - suppress output if true
    */
    public CLMetropolisMC(MCSimulation mc,boolean out){
        this(mc,out,"GPU",0);
    }
    /**
    *      CLMetropolisMC constructor.
    * 
    * @param mc - monte carlo simulation to be updated based on this algo
    * @param out - suppress output if true
    * @param maxConcurThreads - maximum work-items/threads to run
    */
    public CLMetropolisMC(MCSimulation mc,boolean out, int max){
        this(mc,out,"GPU",max);
    }
    /**
    *      CLMetropolisMC constructor.
    * 
    * @param mc - monte carlo simulation to be updated based on this algo
    * @param out - suppress output if true
    * @param device - type of device to run simulation on ("GPU/CPU"
    */
    public CLMetropolisMC(MCSimulation mc,boolean out, String device){
        this(mc,out,device,0);
    }    
    /**
    *      CLMetropolisMC constructor.
    * 
    * @param mc - monte carlo simulation to be updated based on this algo
    * @param out - suppress output if true
    * @param device - type of device to run simulation on ("GPU/CPU"
    * @param maxConcurThreads - maximum work-items/threads to run
    */
    public CLMetropolisMC(MCSimulation mc,boolean out, String device, int maxConcurThreads){
        param = mc.getParams();
        if(param==null){param = new ParameterBank(mc.getParamPostFix());}
        SimPost = mc.getSimPostFix();
        ParamPost = mc.getParamPostFix();
        
        // determine if outputting
        output = out;
        if(output){
            System.out.println("CLMetropolisMC | ###########################");
            System.out.println("CLMetropolisMC | Run:"+mc.getRuns()+"     ");
            System.out.println("CLMetropolisMC | ###########################");
        }
       
        w = new ArrayList<Float>();
        sysParams = new ArrayList<Float>();
        jInteraction = param.jInteraction;
        hField = param.hField;
        temperature = param.temperature;
        useLongRange = param.useLongRange;
        useDilution = param.useDilution;
        useHeter = param.useHeter;
        L = param.L;
        D = param.D;
        R = param.R;
        N = (int)Math.pow(L,D);
        triggerReady=false;
        clhandler = new CLHelper();
        clhandler.initializeOpenCL(device);
        measurer = new MeasureIsingSystem(temperature,N);
        
        maxThreads = (int)clhandler.getCurrentDevice1DMaxWorkItems()/2;  
        //Ran.setSeed(983745);
        // User imposed thread limit
        if(maxConcurThreads>0){maxThreads = maxConcurThreads;}
        
        if(!useLongRange){R=1;metropolisKernel = metropolisNNKernel;}
        clhandler.createKernel("", metropolisKernel);
        clhandler.createKernel("", updateFlKernel);

        divideWork();
        initSysParams();
     
        int nRandom = (2*nThreadRuns*(nThreads)+1);
        
        // sysParams
        clhandler.createFloatBuffer(metropolisKernel, 0,sysParams.size(), sysParams,"r",true);
        clhandler.copyFlBufferAcrossKernel(metropolisKernel, 0, updateFlKernel, buffNumSysParams);
        
        clhandler.setKernelArg(updateFlKernel);
        
        //Random
        clhandler.createFloatBuffer(metropolisKernel, 1,nRandom, 0,"rw");
        //dE
        clhandler.createFloatBuffer(metropolisKernel, 2,nThreads, 0,"rw");
        //error buffer
        clhandler.createFloatBuffer(metropolisKernel, 3,400, 0,"rw");
        
        // Spin system and fixed buffer
        clhandler.createIntBuffer(metropolisKernel, 0,N, param.s,"rw");
        clhandler.createIntBuffer(metropolisKernel, 1,N, 0,"r");
        // dM
        clhandler.createIntBuffer(metropolisKernel, 2,nThreads, 0,"rw");
        if(metropolisKernel.contains("boost")){
            ArrayList<Integer> sumNeigh = new ArrayList<Integer>();
            for(int u = 0;u < N;u++){
                sumNeigh.add((2*param.R+1)*(2*param.R+1)-1);
            }
            clhandler.createIntBuffer(metropolisKernel, 3,N, sumNeigh,"rw",false);
        }
        
        // 
        clhandler.setIntArg(metropolisKernel, 0, nThreads);
        if(useLongRange){      
            clhandler.setIntArg(metropolisKernel, 1, nThreadRuns);
        }
        
        currentSeed = Ran.nextInt();
        RNG = new RandomNumberCL(currentSeed,clhandler);
        RNGKey = RNG.addRandomBuffer(metropolisKernel, 1, nRandom);
        
        System.out.println("CLMetropolisMC |  Filling Lattice Buffer");
        lattice =  new LatticeCL(param.s,mc.getParamPostFix(),
                mc.getFixedPostFix(),mc.instId,clhandler,metropolisKernel,0,1);
        
        System.out.println("CLMetropolisMC |  Lattice Buffer Filled");
        magnetization = ((LatticeCL)lattice).calcCorrectionMag(param.s);
        energy = ((LatticeCL)lattice).calcIsingEnergy();
        System.out.println("CLMetropolisMC |  Energy Calculated");
        
        // Initialize trigger
        setTrigger();
        
        measurer.setN(N, lattice.getNFixed());
        if(output){printParameters();}
        currentTime=0;
    }
    // rescale interaction
    private void rescaleJ(){
        if(param.Geo == 4 || param.Geo == 2){
            jInteraction = 4.0/((2.0*R+1)*(2.0*R+1)-1.0);
        }
    }
    /**
    *   printParameters prints the parameters of this metropolis algorithm
    */
    public void printParameters(){
        System.out.println("CLMetropolisMC | ***********************************************");
        System.out.println("CLMetropolisMC | System Parameters");
        System.out.println("CLMetropolisMC | N : "+N);
        System.out.println("CLMetropolisMC | L : "+L);
        System.out.println("CLMetropolisMC | R : "+R);
        System.out.println("CLMetropolisMC | jInteraction : "+jInteraction);
        System.out.println("CLMetropolisMC | hField : "+hField);
        System.out.println("CLMetropolisMC | temp : "+temperature);
        System.out.println("CLMetropolisMC | Current M: "+magnetization);
        System.out.println("CLMetropolisMC | Device MB Used: "+clhandler.getDeviceUsedMB(metropolisKernel));
        System.out.println("CLMetropolisMC | Max Work Items: "+maxThreads);
        
        System.out.println("CLMetropolisMC | ***********************************************");
        System.out.println("CLMetropolisMC | Work Division Parameters");
        System.out.println("CLMetropolisMC | N Site Visits divided into "+ nThreadRuns
                +" runs of thread groups of "+nThreads+
                "   with local work groups of "+nThreadGroup+"  threads");
        System.out.println("CLMetropolisMC | divided L : "+subL);
        System.out.println("CLMetropolisMC | ***********************************************");
        //printW();
    }
    // divides work for kernel worksize
    private void divideWork(){
        int remainder = L;
        ArrayList<Integer> subLengths = new ArrayList<Integer>();

        subL = 2*R+1;

        for(int u = 0 ; u < 4; u++){
        // Find optimal division of lattice
        while(remainder !=0 && subL<L){
            subL++;
            remainder = L%subL;
        }
        remainder =L;
        subLengths.add(subL);
        }
        int divLengths = 0;
        subL = subLengths.get(0);
        while(subL < minDivision && divLengths < subLengths.size()){
            subL = subLengths.get(divLengths);
            divLengths++;
        }

        // define amount of threads by this division
        nThreads = ((int)L/subL)*((int)L/subL);
        //subL = (int)(L/subL);

        if(metropolisKernel.equalsIgnoreCase("ising2d_nnmetro_checker")){
            nThreads = (int)(N/2);
            subL=L;
        }

        nThreadGroup = nThreads;

        int i=1;
         // If too many threads than divide threads into local workgroups
        while(nThreadGroup> maxThreads && i < nThreads){
            while(remainder !=0 && i < nThreads){
                i++;
                remainder = nThreads % i; 
            }
            nThreadGroup = (int) (nThreads/i);
            if(nThreadGroup> maxThreads){
                i++;
                remainder= nThreads % i; 
                nThreadGroup = (int) (nThreads/i);}
        }
        nThreadGroup = (int) (nThreads/i);
        nThreadRuns = (int)(N/nThreads);
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
    /**
    *      doOneStep performs one step of metropolis algorithm 
    */
    @Override
    public void doOneStep() {
        int dM;
        double dE;
        
        //System.out.println("CLMetropolisMC | Getting Random");
        RNG.fillBufferWithRandom(RNGKey);
        //System.out.println("CLMetropolisMC | Running Metropolis");
        clhandler.runKernel(metropolisKernel, nThreads, nThreadGroup);
        //System.out.println("CLMetropolisMC | Done with Kernel Run");
        
        //Gather total dM and dE
        dM =(int) clhandler.getBufferSumAsDouble(metropolisKernel,"int", 2, nThreads);
        dE = clhandler.getBufferSumAsDouble(metropolisKernel,"float", 2, nThreads);
       
        //clhandler.getFloatBufferAsArray(metropolisKernel, 1, 6, true);
        //clhandler.getIntBufferAsArray(metropolisKernel, 0, 4, true);
        //clhandler.getFloatBufferAsArray(metropolisKernel, 4, (8), true);
        //System.out.println("CLMetropolisMC | Change in M: "+dM);
        magnetization += dM;
        energy += dE;
        
        currentTime++;
        
        if(triggerOn){triggerReady = trigger.triggerNow(magnetization,energy);}
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
        mcs = 0;
        measurer.clearAll();
    }    
    /**
    *         resetSimulation should reset all simulation parameters.
    */
    @Override
    public void resetSimulation(){
        currentTime = 0;
        resetData();
        trigger.reset();
        lattice.initialize(param.s);
        magnetization = lattice.getMagnetization();
        magStaggered = lattice.getMagStaggered();
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
        lattice.initialize(param.s);
        if(param.useHeter){lattice.setFixedLatticeValues();}
        magnetization = lattice.getMagnetization();
        magStaggered = lattice.getMagStaggered();
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
        if(t==0){lattice.initialize(param.s);
        if(param.useHeter){lattice.setFixedLatticeValues();}
        }else{
        lattice.setInitialConfig(t,run,"");}
        magnetization = lattice.getMagnetization();
        magStaggered = lattice.getMagStaggered();
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
        measurer.updateE(energy, true);
        // use staggered magnetization for antiferromagnet
        measurer.updateM((param.jInteraction > 0) ? magnetization:magStaggered, true);
        mcs++;
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
    *   getSpecificHeat should return the current specific heat.
    * 
    */
    @Override
    public double getSpecificHeat() {
        return(measurer.getSpecificHeat());
    }
    /**
    *   getSusceptibility should return the current susceptibility
    * 
    */
    @Override
    public double getSusceptibility() {
        return(measurer.getSusceptibility());
    }
    /**
    *   getSusceptibilityMult should return the current susceptibility multiple values if possible
    * 
    */
    @Override
    public ArrayList<Double> getSusceptibilityMult() {
        ArrayList<Double> sus = new ArrayList<Double>();
        sus.add(measurer.getSusceptibility());
        return (sus);
    }
    /**
    *   flipField should flip the magnetic field
    * 
    */
    @Override
    public void flipField(){
        hField = -hField;
        if(!suppressFieldMessages){
            System.out.println("CLMetropolisMC | Field Flipped at t:"+currentTime+" now hField:"+hField);
        }
        updateSysParams();
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
        measurer.changeTemp(temp);
        temperature = temp;
        updateSysParams();
    }
    /**
    *         changeH should set the magnetic field to the given value.
    * 
    *  @param hnow - new magnetic field
    */ 
    @Override
    public void changeH(double hnow){
        hField = hnow;
        updateSysParams();
    }
    /**
    *         changetTandH should set the temperature to the given value and
    *   the magnetic field.
    * 
    * 
    *  @param temp - new temperature
    *  @param hnow - new field
    */ 
    @Override
    public void changeTandH(double temp,double hnow){
        hField = hnow;
        temperature = temp;
        updateSysParams();
    }
    /**
    *         changeTFlipField should flip the field and set the temperature to
    *   the given value
    * 
    *  @param temp - new temperature
    */ 
    @Override
    public void changeTFlipField(double temp){
        temperature = temp;
        measurer.changeTemp(temp);
        hField = -hField;
        updateSysParams();
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
    public int getM(){return magnetization;}
    /**
    *   getMStag should return the current staggered magnetization.
    * 
    */
    public int getMStag(){return magStaggered;}
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
        return measurer;
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
    
    
    private void initSysParams(){
        rescaleJ();
        sysParams.add((float)jInteraction);
        sysParams.add((float)hField);
        sysParams.add((float)L);
        sysParams.add((float)subL);
        if(!useLongRange){
            sysParams.add((float)temperature);  
        }else{  
            sysParams.add((float)R);  
            sysParams.add((float)temperature);
        }
    }
    
    private void updateSysParams(){
        sysParams.set(0,(float)jInteraction);
        sysParams.set(1,(float)hField);
        sysParams.set(2, (float)L);
        sysParams.set(3,(float)subL);
        if(!useLongRange){
            sysParams.set(4,(float)temperature);
        }else{
            sysParams.set(4,(float)R);
            sysParams.set(5,(float)temperature);
        }
        clhandler.setFloatBuffer(metropolisKernel, buffNumSysParams, sysParams.size(),sysParams,true);
    }
    
    public void closeOpenCL(){
        clhandler.closeOpenCL();
    }
    
    // test the class
    public static void main(String[] args) {
        CLMetropolisMC mc = new CLMetropolisMC("GPU");
        // only true if not using grid
        if(false){
            long time;
            mc.changeH(0.5);
            mc.doOneStep();
            mc.changeH(0.3);
            mc.doOneStep();
            mc.changeT(1.008);
            mc.doOneStep();

            for(int j = 0; j<10;j++){
                time = System.nanoTime();
                mc.doOneStep();
                time = System.nanoTime()-time;
                System.out.println("CLMetropolisMC | M: "+mc.getM());
                System.out.println("CLMetropolisMC | time for step: "+time);}

            mc.changeH(-0.5);

            for(int j = 0; j<20;j++){
                time = System.nanoTime();
                mc.doOneStep();
                time = System.nanoTime()-time;
                System.out.println("CLMetropolisMC | M: "+mc.getM());
                System.out.println("CLMetropolisMC | time for step: "+time);
            }
            mc.closeOpenCL();
        }
        System.out.println("CLMetropolisMC | DONE!");
    }
    
}
package GPUBackend;
/**
 *   @(#) FullyConnectMetroMCCL 
 *
 */

import Backbone.Algo.IsingMC;
import Backbone.System.LatticeMagInt;
import Backbone.Util.DataSaver;
import Backbone.Util.MeasureIsingSystem;
import Backbone.Util.ParameterBank;
import Backbone.Util.SimProcessParser;
import JISim.MCSimulation;
import Triggers.*;
import java.util.ArrayList;
import java.util.Random;
import scikit.opencl.CLHelper;
import scikit.opencl.RandomNumberCL;

/**
 *  FullyConnectMetroMCCL 
 *
 *
 *   <br>
 *
 * @author James B. Silva <jbsilva@bu.edu>
 * @since May 2012
 */
public final class FullyConnectMetroMCCL implements IsingMC {
    private CLHelper clhandler;
    private ParameterBank param;
    private String metropolisKernel = "ising_fullymetro";
    private String updateFlKernel = "update_fl_2buffer";
    private LatticeMagInt lattice;
    private ArrayList<Float> jInteractions;
    private ArrayList<Float> temperatures;
    private ArrayList<Float> hFields;
    private ArrayList<Integer> Ms;
    private ArrayList<Float> energies;
    private ArrayList<Integer> NupS;
    private ArrayList<Integer> NdownS;
    private ArrayList<Integer> NupFixS;
    private ArrayList<Integer> NdownFixS;    
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
    private int run;
    private boolean output = false;
    private boolean triggerReady=false;
    private Random Ran = new Random();
    private ArrayList<MeasureIsingSystem> measurers = new ArrayList<MeasureIsingSystem>();
    private String SimPost="";
    private String ParamPost="";  
    private int magnetization;
    private RandomNumberCL RNG;
    private int RNGKey;
    
    private int tFlip;
    public double energy;  
    public int currentSeed;
    public int mcs = 0; // number of MC moves per spin
    public int currentTime;
    
    public FullyConnectMetroMCCL(String device){
        param = new ParameterBank("");
        
        hField = param.hField;
        temperature = param.temperature;
        useLongRange = param.useLongRange;
        useDilution = param.useDilution;
        useHeter = param.useHeter;
        N = param.L*param.L;
        jInteraction = 4.0*param.jInteraction/N;
        System.out.println("FullyConnectMetroMCCL |  Evolving "+numSystems+"  systems in OpenCL");
        
        createSystemBuffers(null,null);
        
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
        
        int nRandom = (2*N*(numSystems)+1);
       
        // float parameters
        // energy 
        clhandler.createFloatBuffer(metropolisKernel, 0,energies.size(), energies,"rw",false);
        
        // hField for systems
        clhandler.createFloatBuffer(metropolisKernel, 1,hFields.size(), hFields,"r",false);
        // copy to updateFL kernel to allow for updating outside of Metropolis kernel
        clhandler.copyFlBufferAcrossKernel(metropolisKernel, 1, updateFlKernel, 0);
        
        // temperatures for system
        clhandler.createFloatBuffer(metropolisKernel, 2,temperatures.size(), temperatures,"r",false);
        // copy to updateFL kernel to allow for updating outside of Metropolis kernel
        clhandler.copyFlBufferAcrossKernel(metropolisKernel, 2, updateFlKernel, 1);
        
        // jInteractions for system
        clhandler.createFloatBuffer(metropolisKernel, 3,jInteractions.size(), jInteractions,"r",false);
        
        // random
        clhandler.createFloatBuffer(metropolisKernel, 4,nRandom,0,"rw");
        
        // error buffer
        clhandler.createFloatBuffer(metropolisKernel, 5,400, 0,"rw");
        
        clhandler.setKernelArg(updateFlKernel);
        
        // integer buffers
        // spin values
        clhandler.createIntBuffer(metropolisKernel,0,Ms.size(), Ms,"rw",false);
        // N up
        clhandler.createIntBuffer(metropolisKernel,1,NupS.size(), NupS,"rw",false);
        // N down
        clhandler.createIntBuffer(metropolisKernel,2,NdownS.size(), NdownS,"rw",false);
        // N up fixed
        clhandler.createIntBuffer(metropolisKernel,3,NupFixS.size(), NupFixS,"r",false);
        // N down fixed
        clhandler.createIntBuffer(metropolisKernel,4,NdownFixS.size(), NdownFixS,"r",false);
        
        // integer arguments
        // Number of systems
        clhandler.setIntArg(metropolisKernel, 0, numSystems);
        clhandler.setIntArg(metropolisKernel, 1, N);
        
        clhandler.setKernelArg(metropolisKernel);
   
        currentSeed = Ran.nextInt();
        RNG = new RandomNumberCL(currentSeed,clhandler);
        // Get the keyname for the RNG to buffer 4 the randome number buffer
        RNGKey = RNG.addRandomBuffer(metropolisKernel, 4, nRandom);
        
        //System.out.println("FullyConnectMetroMCCL |  Filling LatticeMagInt Buffer");        
        //lattice =  new LatticeCL(param.s,"",
          //      "",0,clhandler,metropolisKernel,0,1);

        //System.out.println("FullyConnectMetroMCCL |  LatticeMagInt Buffer Filled");
        //System.out.println("FullyConnectMetroMCCL |  Energy Calculated");
        
        // Initialize trigger
        setTrigger();
        
        for(int u = 0; u < numSystems;u++){
            measurers.get(u).setN(N, (NdownFixS.get(u)+NupFixS.get(u)));
        }
        if(output){printParameters();}
        currentTime=0;
    }
    
    public FullyConnectMetroMCCL(MCSimulation mc, boolean out,int nfix,int sfix,boolean bal){
        this(mc,out,"GPU",nfix,sfix,bal,256,-1);
    }
    
    public FullyConnectMetroMCCL(MCSimulation mc, boolean out,int nfix,int sfix,
            boolean bal, ArrayList<Float> fieldsIn, ArrayList<Integer> Nin){
        this(mc,out,"GPU",nfix,sfix,bal,256,-1,fieldsIn,Nin);
    }
    
    public FullyConnectMetroMCCL(MCSimulation mc, boolean out,int nfix,int sfix,boolean bal, int systems){
        this(mc,out,"GPU",nfix,sfix,bal,systems,-1);
    }
    
    public FullyConnectMetroMCCL(MCSimulation mc, boolean out,int nfix,int sfix,boolean bal, 
            int systems, ArrayList<Float> fieldsIn,ArrayList<Integer> Nin){
        this(mc,out,"GPU",nfix,sfix,bal,systems,-1,fieldsIn,Nin);
    }
    
    public FullyConnectMetroMCCL(MCSimulation mc, boolean out,String device,int nfix,
            int sfix,boolean bal, int systems,int maxConcurThreads){
        this(mc,out,device,nfix,sfix,bal,systems,maxConcurThreads,null,null);
    }
    
    public FullyConnectMetroMCCL(MCSimulation mc, boolean out,String device,int nfix,
            int sfix,boolean bal, int systems,int maxConcurThreads, 
            ArrayList<Float> fieldsIn,ArrayList<Integer> Nin){
        param = mc.getParams();
        if(param==null){param = new ParameterBank(mc.getParamPostFix());}
        SimPost = mc.getSimPostFix();
        ParamPost = mc.getParamPostFix();
        currentSeed = mc.currentSeed;
        numSystems = systems;
        nfixed = nfix;
        spinFix = sfix;
        balancedConfig = bal;
        
        // determine if outputting
        output = out;
        if(output){
            System.out.println("FullyConnectMetroMCCL | ###########################");
            System.out.println("FullyConnectMetroMCCL | Run:"+mc.getRuns()+"     ");
            System.out.println("FullyConnectMetroMCCL | ###########################");
        }
        hField = param.hField;
        temperature = param.temperature;
        useLongRange = param.useLongRange;
        useDilution = param.useDilution;
        useHeter = param.useHeter;
        s0 = param.s;
        N = param.L;
        param.setN(N);
        jInteraction = 4.0*param.jInteraction/N;
        triggerReady=false;
        
        createSystemBuffers(fieldsIn,Nin);
        
        clhandler = new CLHelper();
        clhandler.setOutputMode(output);
        clhandler.initializeOpenCL(device);
        for(int u = 0; u < numSystems;u++){
            measurers.add(new MeasureIsingSystem(temperature,N));
            measurers.get(u).setFixedVal((balancedConfig)? 0 :spinFix);
        }
        
        maxThreads = (int)clhandler.getCurrentDevice1DMaxWorkItems()/2;  
        //Ran.setSeed(983745);
        // User imposed thread limit
        if(maxConcurThreads > 0){maxThreads = maxConcurThreads;}
        
        clhandler.createKernel("", metropolisKernel);
        clhandler.createKernel("", updateFlKernel);

        divideWork();
        
        int nRandom = (2*N*(numSystems)+1);
       
        // float parameters
        // energy 
        clhandler.createFloatBuffer(metropolisKernel, 0,energies.size(), energies,"rw",false);
        
        // hField for systems
        clhandler.createFloatBuffer(metropolisKernel, 1,hFields.size(), hFields,"r",false);
        // copy to updateFL kernel to allow for updating outside of Metropolis kernel
        clhandler.copyFlBufferAcrossKernel(metropolisKernel, 1, updateFlKernel, 0);
        
        // temperatures for system
        clhandler.createFloatBuffer(metropolisKernel, 2,temperatures.size(), temperatures,"r",false);
        // copy to updateFL kernel to allow for updating outside of Metropolis kernel
        clhandler.copyFlBufferAcrossKernel(metropolisKernel, 2, updateFlKernel, 1);
        
        // jInteractions for system
        clhandler.createFloatBuffer(metropolisKernel, 3,jInteractions.size(), jInteractions,"r",false);
        
        // random
        clhandler.createFloatBuffer(metropolisKernel, 4,nRandom,0,"rw");
        
        // error buffer
        clhandler.createFloatBuffer(metropolisKernel, 5,400, 0,"rw");
        
        clhandler.setKernelArg(updateFlKernel);
        
        // integer buffers
        // spin values
        clhandler.createIntBuffer(metropolisKernel,0,Ms.size(), Ms,"rw",false);
        // N up
        clhandler.createIntBuffer(metropolisKernel,1,NupS.size(), NupS,"rw",false);
        // N down
        clhandler.createIntBuffer(metropolisKernel,2,NdownS.size(), NdownS,"rw",false);
        // N up fixed
        clhandler.createIntBuffer(metropolisKernel,3,NupFixS.size(), NupFixS,"r",false);
        // N down fixed
        clhandler.createIntBuffer(metropolisKernel,4,NdownFixS.size(), NdownFixS,"r",false);
        
        // integer arguments
        // Number of systems
        clhandler.setIntArg(metropolisKernel, 0, numSystems);
        clhandler.setIntArg(metropolisKernel, 1, N);
        
        clhandler.setKernelArg(metropolisKernel);
   
        currentSeed = Ran.nextInt();
        RNG = new RandomNumberCL(currentSeed,clhandler);
        // Get the keyname for the RNG
        RNGKey = RNG.addRandomBuffer(metropolisKernel, 4, nRandom);
        
        System.out.println("FullyConnectMetroMCCL |  Filling Lattice Buffer");
        //lattice =  new LatticeCL(param.s,mc.getParamPostFix(),
        //        mc.getFixedPostFix(),mc.instId,clhandler,metropolisKernel,0,1);
        
        //System.out.println("FullyConnectMetroMCCL |  LatticeMagInt Buffer Filled");
        //magnetization = ((LatticeCL)lattice).calcCorrectionMag(param.s);
        //energy = ((LatticeCL)lattice).calcIsingEnergy();
        //System.out.println("FullyConnectMetroMCCL |  Energy Calculated");
        
        // Initialize trigger
        setTrigger();
        
        for(int u = 0; u < numSystems;u++){
            int Ndil = N - (NupS.get(u)+NdownS.get(u)+NupFixS.get(u)+NdownFixS.get(u));
            measurers.get(u).setN(N, (NdownFixS.get(u)+NupFixS.get(u)+Ndil));
            measurers.get(u).changeTemp(temperature);
        }
        if(output){printParameters();}
        currentTime=0;
    }
    
    public void printParameters(){
        System.out.println("FullyConnectMetroMCCL | ***********************************************");
        System.out.println("FullyConnectMetroMCCL | System Parameters");
        //System.out.println("FullyConnectMetroMCCL | N : "+N);
        //System.out.println("FullyConnectMetroMCCL | jInteraction : "+jInteraction);
        //System.out.println("FullyConnectMetroMCCL | hField : "+hField);
        //System.out.println("FullyConnectMetroMCCL | temp : "+temperature);
        System.out.println("FullyConnectMetroMCCL | Current M: "+magnetization);
        System.out.println("FullyConnectMetroMCCL | Device MB Used: "+clhandler.getDeviceUsedMB(metropolisKernel));
        System.out.println("FullyConnectMetroMCCL | Max Work Items: "+maxThreads);
        
        System.out.println("FullyConnectMetroMCCL | ***********************************************");
        System.out.println("FullyConnectMetroMCCL | Work Division Parameters");
        System.out.println("FullyConnectMetroMCCL | N Site Visits divided into "+ nThreadRuns
                +" runs of thread groups of "+nThreads+
                "   with local work groups of "+nThreadGroup+"  threads");
        System.out.println("FullyConnectMetroMCCL | ***********************************************");
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
        System.out.println("FullyConnectMetroMCCL | Number of threads total : "+nThreads+"    worked on in groups of "+nThreadGroup);
        
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
        
        //System.out.println("FullyConnectMetroMCCL | Getting Random");
        RNG.fillBufferWithRandom(RNGKey);
        //System.out.println("FullyConnectMetroMCCL | Running Metropolis");
        clhandler.runKernel(metropolisKernel, nThreads, nThreadGroup);
        //System.out.println("FullyConnectMetroMCCL | Done with Kernel Run");
        
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
        ArrayList<Float> energySys = clhandler.getBufferFlAsArrayList(metropolisKernel, 0,numSystems, false);
        ArrayList<Integer> magSys = clhandler.getBufferIntAsArrayList(metropolisKernel, 0,numSystems, false);
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
        return (hFields);
    }

    public ArrayList<Integer> getNfixedMult() {
        ArrayList<Integer> fixSys = new ArrayList<Integer>();
        for(int u = 0; u < numSystems;u++){
            fixSys.add(NdownFixS.get(u)+NupFixS.get(u));
        }
        return (fixSys);
    }
    
    public ArrayList<Integer> getMagMult() {
        ArrayList<Integer> magSys = clhandler.getBufferIntAsArrayList(metropolisKernel, 0,numSystems, false);
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
        if(!suppressFieldMessages){System.out.println("FullyConnectMetroMCCL | Field Flipped at t:"+currentTime+" now hField:"+hField);}
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
    *   getJinteraction should return the current strength of interaction.
    * 
    */
    @Override
    public double getJinteraction(){
        return jInteraction;
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
        ArrayList<Integer> magSys = clhandler.getBufferIntAsArrayList(metropolisKernel, 0,averagingNum, false);
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
            hFields.set(u, hFields.get(u)+diff);
            temperatures.set(u, temperatures.get(u)+difft);
        }
        clhandler.setFloatBuffer(metropolisKernel, 1, hFields.size(),hFields,true);
        clhandler.setFloatBuffer(metropolisKernel, 2, temperatures.size(),temperatures,true);
        clhandler.runKernel(updateFlKernel, hFields.size(), 1);
        //System.out.println("New Fields");
        //clhandler.getBufferFlAsArrayList(metropolisKernel, 1, 10, true);
        //System.out.println("New Temp");
        //clhandler.getBufferFlAsArrayList(metropolisKernel,2, 10, true);
    }
    
    private void updateField(double field, double fieldOld){
        float diff = (float)(field-fieldOld);
        for(int u = 0; u < numSystems; u++){
            hFields.set(u, hFields.get(u)+diff);
        }
        clhandler.setFloatBuffer(metropolisKernel, 1, hFields.size(),hFields,true);
        clhandler.runKernel(updateFlKernel, hFields.size(), 1);
        //System.out.println("New Fields");
        //clhandler.getBufferFlAsArrayList(metropolisKernel,  1, 10, true);
        
    }

    private void updateTemp(double temp,double tempOld){
        float diff = (float)(temp-tempOld);
        for(int u = 0; u < numSystems; u++){
            temperatures.set(u, temperatures.get(u)+diff);
        }
        clhandler.setFloatBuffer(metropolisKernel, 2, temperatures.size(),temperatures,true);
        clhandler.runKernel(updateFlKernel, temperatures.size(), 1);
        //System.out.println("New Temp");
        //clhandler.getBufferFlAsArrayList(metropolisKernel, 2, 10, true);
    }
    
    
    public void createSystemBuffers(ArrayList<Float> fieldsIn,ArrayList<Integer> Nin){
        ArrayList<Integer> ms = new ArrayList<Integer>();
        ArrayList<Integer> NuS = new ArrayList<Integer>();
        ArrayList<Integer> NdS = new ArrayList<Integer>();
        ArrayList<Integer> NuSF = new ArrayList<Integer>();
        ArrayList<Integer> NdSF = new ArrayList<Integer>();
        ArrayList<Float> Es = new ArrayList<Float>();
        ArrayList<Float> hs = new ArrayList<Float>();
        ArrayList<Float> ts = new ArrayList<Float>();
        ArrayList<Float> js = new ArrayList<Float>();
        int avgMag = 0; int mag0=0;
        
        // if inputting field then treat as a 0 fixed spin system
        boolean reRun = (Nin != null)? true: false;
        if(reRun){
            nfixed = 0;balancedConfig = false;
            System.out.println("FullyConnectMetroMCCL | Rerun data in pure system");
        }
        for(int u = 0; u < numSystems; u++){
            ArrayList<Integer> temp = makeMag(s0);
            mag0 = temp.get(2);
            if(reRun){N = Nin.get(u);} 
            ArrayList<Integer> temp2 = createFixed(temp.get(0), temp.get(1),nfixed, N, spinFix, balancedConfig);
            avgMag = temp2.get(2)-temp2.get(3);
            NuS.add(temp2.get(0));
            NdS.add(temp2.get(1));      
            NuSF.add(temp2.get(2));
            NdSF.add(temp2.get(3));
            ms.add(temp2.get(0)+temp2.get(2)-temp2.get(1)-temp2.get(3));
            float heff = 4.0f*((float)(temp2.get(2)-temp2.get(3)))/((float)N);
            hs.add(heff+((float)hField));
            ts.add(((float)temperature));
            jInteraction = (Nin != null)? (Math.signum(jInteraction)*4/Nin.get(u)): jInteraction;
            js.add((float)jInteraction);
            Es.add(calcIsingEnergy(temp2.get(0)+temp2.get(2)-temp2.get(1)-temp2.get(3)));
        }
        hFields = (reRun)? fieldsIn : hs;
        temperatures = ts;
        jInteractions = js;
        Ms = ms;
        energies = Es;
        NupS = NuS;
        NdownS = NdS;
        NupFixS = NuSF;
        NdownFixS = NdSF;
        
        //for(int u = 0;u < 25; u++){
        //    System.out.println("TESTTESTTEST++++++++++==== Mag IN: "+Ms.get(u));
        //}
        
        magnetization = ((int)((double)mag0/(double)N)+(int)((double)avgMag/(double)numSystems))*N;
        
    }
    
    private ArrayList<Integer> makeMag(int s){
        int Nup = 0; int Ndown = 0;
        if(s == 1){
            Nup = N; Ndown = 0;
        }else if(s == (-1)){
            Nup = 0; Ndown = N;
        }else if(s == 0){
            Nup = (N/2);
            Ndown = Nup;
        }else{
            Nup = (int)(N * Ran.nextDouble());
            Ndown = N = Nup;
        }
        magnetization = Nup-Ndown;
        ArrayList<Integer> sys = new ArrayList<Integer>();
        sys.add(Nup);
        sys.add(Ndown);
        sys.add((Nup-Ndown));
        return sys;
    }
    
    public ArrayList<Integer> createFixed(int Nup, int Ndown,int fix, int num, int sFix, boolean balanced){
        if( num > N){
            int dif = num-N;
            Nup += dif/2;
            Ndown += dif/2;
        }
        int NfixUp = 0;
        int NfixDown = 0;
        double fixProb = (double)(fix)/((double)num);
        // balanced spin fix zero is same as random diluted in intent
        if(balanced && sFix == 0){balanced = false;}
        for(int u = 0; u < num;u++){
            if(Ran.nextDouble() <= fixProb){
                // choose fix 
                if(!balanced){
                    if(sFix == 1){NfixUp++;}
                    if(sFix == (-1)){NfixDown++;}
                }else{
                    if(Ran.nextDouble() <= 0.5){
                        NfixUp++;
                    }else{
                        NfixDown++;
                    }
                }
                // Choose fix point in lattice to fix
                if(Nup == 0 ){
                    Ndown--;
                }else if(Ndown == 0 ){
                    Nup--;
                }else{
                    if(Ran.nextDouble() < 0.5){
                            Nup--;
                    }else{
                            Ndown--;
                    }
                }
            }
        }
        ArrayList<Integer> sys = new ArrayList<Integer>();
        sys.add(Nup);
        sys.add(Ndown);
        sys.add(NfixUp);
        sys.add(NfixDown);        
        return sys;
    }
    
        /**
    *       calcEnergy calculates the energy of the lattice in an Ising Model 
    * 
    * @return energy in lattice for Ising model
    */
    private float calcIsingEnergy(int mag){
        double e=0.0;
        // Simple version 
        double m = (double)mag;
        double num = (double) N;
        e = -((0.5*jInteraction*m*m)+(hField*m));
        // ammar
        //e =  (jInteraction)*(0.5*(N - magnetization*magnetization)) - (magnetization*hField);
        // Kang
        //e = -jInteraction*(mag)+jInteraction+hField;
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
        FullyConnectMetroMCCL mc = new FullyConnectMetroMCCL("GPU");
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
            System.out.println("FullyConnectMetroMCCL | M: "+ mc.getM());
            System.out.println("FullyConnectMetroMCCL | time/system (us) for step: "+(time/(1000*mc.numSystems)));
        }

        mc.changeH(-0.005);
        
        for(int j = 0; j< 400;j++){
            time = System.nanoTime();
            mc.doOneStep();
            time = System.nanoTime()-time;
            System.out.println("FullyConnectMetroMCCL | M: "+mc.getM());
            System.out.println("FullyConnectMetroMCCL | time/system (us) for step: "+(time/(1000*mc.numSystems)));
        }
        
        
            mc.closeOpenCL();
            System.out.println("FullyConnectMetroMCCL | DONE!");
    }
    
}
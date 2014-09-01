package Backbone.Algo;

/**
 * 
 *    @(#)   MetropolisMCThreaded
 */  

import Backbone.System.AtomicLattice;
import Backbone.System.LatticeMagInt;
import Backbone.Util.MeasureIsingSystem;
import Backbone.Util.DataSaver;
import Backbone.Util.SimProcessParser;
import Backbone.Util.ParameterBank;
import JISim.MCSimulation;
import Triggers.DeviationTrigger;
import Triggers.SimpleValFitTrigger;
import Triggers.SimpleValueTrigger;
import Triggers.Trigger;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.*;
import java.util.logging.Level;
import java.util.logging.Logger;


 /** 
 *   Threaded Implementation of Metropolis Monte Carlo Simulation.
 *  Many threaded implementation of simulation based on amount of available cores.
 *  Synchronization seems to drag down class.
 *  <br>
 *  @param mc - MCSimulation - Monte Carlo simulation class hosting 
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2012-01    
 */


public final class MetropolisMCThreaded implements IsingMC  {
	private ConcurrentHashMap<Integer,Double> w ;
        private ParameterBank param;
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
	private int Q;
	private int L;
	private int N;
	private int D;
	private int run;
        private int threadChunk =1;
	private boolean triggerReady=false;
        private MeasureIsingSystem measurer;
        private Random Ran = new Random();
        private int magnetization = 0;
        private int magStaggered = 0;
        private boolean triggerOn = true;
        private String SimPost="";
        
        public double energy;  
        public int currentSeed;
        public int mcs = 0; // number of MC moves per spin
        public int currentTime;
        public int tFlip=0;
        
    public MetropolisMCThreaded(){
        param = new ParameterBank("");

        jInteraction = param.jInteraction;
        hField = param.hField;
        temperature = param.temperature;
        useLongRange = param.useLongRange;
        useDilution = param.useDilution;
        useHeter = param.useHeter;
        L = param.L;
        N = param.N;
        D = param.D;
        triggerReady=false;
        measurer = new MeasureIsingSystem(temperature,N);
        threadChunk = divideStepWork();
        w = new ConcurrentHashMap<Integer,Double>();
        lattice = new AtomicLattice(1);
        measurer.setN(N, lattice.getNFixed());
    }

    /**  
     * @param mc - MCSimulation - Monte Carlo simulation class hosting 
    */
    public MetropolisMCThreaded(MCSimulation mc){
        this(mc,"");
    }
    
    /** 
    *     @param mc - MCSimulation - Monte Carlo simulation class hosting 
    *     @param postfixFname - postfix for file.
    */
    public MetropolisMCThreaded(MCSimulation mc, String postfixfname){
        param = new ParameterBank(postfixfname);
        SimPost = mc.getSimPostFix();
        jInteraction = param.jInteraction;
        hField = param.hField;
        temperature = param.temperature;
        useLongRange = param.useLongRange;
        useDilution = param.useDilution;
        useHeter = param.useHeter;
        L = param.L;
        N = param.N;
        D = param.D;
        triggerReady=false;
        threadChunk = divideStepWork();
        measurer = new MeasureIsingSystem(temperature,N);
        //DELETE
        System.out.print("Run:"+mc.getRuns()+"     ");
        run = mc.getRuns();  
        w = new ConcurrentHashMap<Integer,Double>();
        if(mc.getSimSystem().getClass().getName().contains("attice")){
            lattice = (LatticeMagInt)mc.getSimSystem();
        }else{
            System.err.println("THIS CLASSS (METROPOLISMCTHREADED NEEDS A LATTICE AS SYSTEM.");
        }
        magnetization = lattice.getMagnetization();
        magStaggered = lattice.getMagStaggered();
        measurer.setN(N, lattice.getNFixed());
        SimProcessParser parser = new SimProcessParser(null,SimPost,postfixfname);
        tFlip = parser.timeToFlipField();
        
        calcW();
        setTrigger();
        currentTime=0;
    }

    /**
    *         calcW calculates the transition matrix that is necessary
    *   for the metropolis to get the statistics with the right temperature
    *   and field. Supports Long Range and Nearest Neighbor
    * 
    */ 
    public void calcW(){

        if (useDilution == false && useLongRange == false){
        if(suppressFieldMessages){}else{

        System.out.println("Calc W with j:"+jInteraction+"      h:"+hField+"    temp:"+temperature);}


        w.put(1+1*4+4*4*2, Math.exp((-8*jInteraction-2*hField)/temperature));
        w.put(0+1*4+4*4*2, Math.exp((8*jInteraction+2*hField)/temperature));
        w.put(1+1*4+2*4*2, Math.exp((4*jInteraction-2*hField)/temperature));
        w.put(0+1*4+2*4*2, Math.exp((4*jInteraction+2*hField)/temperature));
        w.put(1+1*4, Math.exp((-2*hField)/temperature));
        w.put(0+1*4, Math.exp((2*hField)/temperature));
        w.put(1+2*4*2, Math.exp((4*jInteraction-2*hField)/temperature));
        w.put(0+2*4*2, Math.exp((-4*jInteraction+2*hField)/temperature));
        w.put(1+4*4*2, Math.exp((8*jInteraction-2*hField)/temperature));
        w.put(4*4*2, Math.exp((-8*jInteraction+2*hField)/temperature));
        }
        else if(useDilution == true && useLongRange == false){
        w.put(1+4*1+4*2*4, Math.exp((-8*jInteraction-2*hField)/temperature)); 
        w.put(0+4*1+4*2*4, Math.exp((8*jInteraction+2*hField)/temperature));
        w.put(1+4*1+4*2*3, Math.exp((-6*jInteraction-2*hField)/temperature));
        w.put(0+4*1+4*2*3, Math.exp((6*jInteraction+2*hField)/temperature));
        w.put(1+4*1+4*2*2, Math.exp((-4*jInteraction-2*hField)/temperature));
        w.put(0+4*1+4*2*2, Math.exp((4*jInteraction+2*hField)/temperature));
        w.put(1+4*1+4*2*1, Math.exp((-2*jInteraction-2*hField)/temperature));
        w.put(0+4*1+4*2*1, Math.exp((2*jInteraction+2*hField)/temperature));
        w.put(1+4*1+4*2*0, Math.exp((-2*hField)/temperature));
        w.put(0+4*1+4*2*0, Math.exp((2*hField)/temperature));
        w.put(1+4*0+4*2*1, Math.exp((2*jInteraction-2*hField)/temperature));
        w.put(0+4*0+4*2*1, Math.exp((-2*jInteraction+2*hField)/temperature));
        w.put(1+4*0+4*2*2, Math.exp((4*jInteraction-2*hField)/temperature));
        w.put(0+4*0+4*2*2, Math.exp((-4*jInteraction+2*hField)/temperature));
        w.put(1+4*0+4*2*3, Math.exp((6*jInteraction-2*hField)/temperature));
        w.put(0+4*0+4*2*3, Math.exp((-6*jInteraction+2*hField)/temperature));
        w.put(1+4*0+4*2*4, Math.exp((8*jInteraction-2*hField)/temperature));
        w.put(0+4*0+4*2*4, Math.exp((-8*jInteraction+2*hField)/temperature));
        // 2 and 3 are both current spin 0 changed to - and + respectively
        w.put(3+4*0+4*2*4, Math.exp((-4*jInteraction-hField)/temperature)); 
        w.put(2+4*1+4*2*4, Math.exp((4*jInteraction+hField)/temperature));
        w.put(3+4*0+4*2*3, Math.exp((-3*jInteraction-hField)/temperature));
        w.put(2+4*1+4*2*3, Math.exp((3*jInteraction+hField)/temperature));
        w.put(3+4*0+4*2*2, Math.exp((-2*jInteraction-hField)/temperature));
        w.put(2+4*1+4*2*2, Math.exp((2*jInteraction+hField)/temperature));
        w.put(3+4*0+4*2*1, Math.exp((-jInteraction-hField)/temperature));
        w.put(2+4*1+4*2*1, Math.exp((jInteraction+hField)/temperature));
        w.put(3+4*0+4*2*0, Math.exp((-hField)/temperature));
        w.put(2+4*1+4*2*0, Math.exp((hField)/temperature));
        }
        else{

        Q =  lattice.getNinRange();

        jInteraction = (Math.signum(jInteraction) > 0)?4.0/Q:-4.0/Q;

        if(suppressFieldMessages){}else{
            System.out.println("Long R Calc W with j:"+jInteraction+
            "      h:"+hField+"    temp:"+temperature
            + "   neighborsInRange:"+Q);
        }

        for(int k = 0;k<Q+1;k++){    
        for(int j = 0;j<2;j++){
            for(int i = 0;i<4;i++){
            int jpre,hpre;
            if((i+j) % 2 == 0){jpre = -1;}else{jpre=1;}
            if((i) % 2 == 0){hpre = 1;}else{hpre=-1;}
            double ener = 2*(k*jpre*jInteraction+hpre*hField);

            w.put(i+j*4+k*4*2, Math.exp((ener)/temperature)); 	
        }}}

        }
    }
    
    /**
    *        doOneStep should do a single monte carlo step.
    */ 
    @Override
    public void doOneStep(){

        currentTime++;

        // threadChunk = divideStepWork();
        ExecutorService pool = Executors.newFixedThreadPool((int)(N/threadChunk));
        Set<Future<Double>> set = new HashSet<>();

        energy = 0;


        // Pass same seed so that random movement is synchronized across amount of threads
        int stemp = Ran.nextInt();
        // Create threads that split the MCstep work
        for(int i=0; i< ((int)(N/threadChunk));i++){
            //MetropolisStepThread(int t, LatticeMagInt lat,double w0[][][], 
            //int loop,int r0,double jit,
            //double hfiel, boolean useHet,boolean useDil, int see)
            Callable<Double> callable = new MetropolisStepThread(currentTime,lattice,w,
                    threadChunk,param.R,jInteraction,hField,
                    param.useHeter,param.useDilution,stemp,i);
            Future<Double> future = pool.submit(callable);
            set.add(future);

            try {
            //System.out.println("Result"+future.get());
                energy += future.get();
            } catch (    InterruptedException | ExecutionException ex) {
            Logger.getLogger(MetropolisMCThreaded.class.getName()).log(Level.SEVERE, null, ex);
            }
        }

        /*// Sum up change in energy from the threads
        for (Future<Double> future : set) {
            try {
                double en = future.get().doubleValue();
                //System.out.println("Adding segment de ="+en);
                //energy += en;
            } catch (InterruptedException ex) {
                    Logger.getLogger(MetropolisMCThreaded.class.getName()).log(Level.SEVERE, null, ex);
            } catch (ExecutionException ex) {
                    Logger.getLogger(MetropolisMCThreaded.class.getName()).log(Level.SEVERE, null, ex);
            }
        }*/

        //System.out.println((int)(N/threadChunk)+"  threads and "+threadChunk+"  size chunk");


        // Wait till after the flip to measure the susceptibility
        if(currentTime>(tFlip+tAccumulate)){
            updateAccumulators();
        }
    
        if(triggerOn){triggerReady = trigger.triggerNow(magnetization);}
        if(triggerReady){tTrigger = trigger.getTriggerTime();}
    }

    /**
    *    divideStepWork divides the N amount of work into segments determined
    *   by the amount of available processor cores
    */
    private int divideStepWork(){
        int profactor =1;
        // largest number up to processor that evenly divides work;    

        int nmax = Runtime.getRuntime().availableProcessors();

        for(int i =1; i<(nmax+1);i++){
        if ((N % i) ==0){profactor=i;}
        }

        System.out.println(nmax+ "  processors available but going to use only "+profactor);
        System.out.println("_____________________________________________________________");
        profactor = (int) (N/profactor);

        return profactor;
    }
    /**
    *        isPostitive returns true if number is positive
    * 
    *   @param number - input number
    */ 
    public static int isPositive(int number){
        if (number >= 0){return 1;}else{return 0;}
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
        measurer.updateM(magnetization, true);
        mcs++;
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
    *   flipField should flip the magnetic field
    * 
    */
    @Override
    public void flipField(){
        hField = -hField;
        if(suppressFieldMessages){}else{
            System.out.println("Field Flipped at t:"+currentTime+" now hField:"+hField);}
        calcW();
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
        calcW();
    }

    /**
    *         changeH should set the magnetic field to the given value.
    * 
    *  @param hnow - new magnetic field
    */ 
    @Override
    public void changeH(double hnow){
        hField = hnow;
        calcW();
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
        hField = hnow;
        temperature = temp;
        measurer.changeTemp(temp);
        calcW();
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
        calcW();
    }

    
    /**
    *         setTrigger should initialize the trigger to be used in the simulation
    */ 
    @Override
    public void setTrigger(){
        trigger = param.setProperTrigger(param, trigger, SimPost,false);
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
    *         getW returns the value of the transition matrix w at given index
    * 
    *   @param k - index
    */ 
    public double getW(int k){return (double) w.get(k);} 

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
    *   getJinteraction should return the current strength of interaction.
    * 
    */
    @Override
    public double getJinteraction(){
        return jInteraction;
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
        
    // test the class
    public static void main(String[] args) {
        MetropolisMCThreaded mc = new MetropolisMCThreaded();
        mc.calcW();
        mc.setSeed(1928494);
        mc.setTrigger();
        long t = System.nanoTime();
        System.out.println(mc.getW(8));
        System.out.println("Magnetization Initial:"+mc.getM());

        mc.doOneStep();
        t = System.nanoTime()-t;
        System.out.println("t:"+(double)(t/1000000)+"     Magnetization Final:"+mc.getM());
        t = System.nanoTime();
        mc.doOneStep();
        t = System.nanoTime()-t;
        System.out.println("t:"+(double)(t/1000000)+"     Magnetization Final2:"+mc.getM());

        int Lx = 4;
        int Ly = 2;
        int Lz = 8;

        /*for(int k = 0;k<Lz;k++){
        for(int j = 0;j<Ly;j++){            
        for(int i = 0;i<Lx;i++){
            if(i==(Lx-1)){System.out.println("  "+(i+j*Lx+k*Lx*Ly));}else{
            System.out.print("  "+(i+j*Lx+k*Lx*Ly));}       
        }}}*/
        
        System.out.println("Done!");
    }
}



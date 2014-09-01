package Backbone.Algo;

/**
 * 
 *      @(#)  WangLandauMC
 */  

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
import java.util.Random;

/**  
 *   Basic implementation of Wang-Landau Monte Carlo sampling simulation.
 *  Single thread implementation of simulation
 *  <br>
 *  @param mc - MCSimulation - Monte Carlo simulation class hosting 
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2012-08    
 */
public final class WangLandauMC  implements IsingMC {
    private ParameterBank param;
    private LatticeMagInt lattice;
    private Trigger trigger;
    private double temperature;
    private double jInteraction;
    private double hField;
    private boolean useLongRange;
    private boolean useDilution;
    private boolean useHeter;
    private int tTrigger =0;
    private int tAccumulate =1000; // begin accumulator
    private int Q;
    private int L;
    private int N;
    private int D;
    private int run;
    private boolean output = true;
    private boolean noFieldMessages = true;
    private boolean triggerReady=false;
    private boolean triggerOn=false;
    private Random Ran = new Random();
    private String SimPost="";
    private String ParamPost="";  
    
    private double energy;  
    private int magStaggered = 0;
    private int magnetization = 0;
    public int currentSeed;
    public int mcs = 0; // number of MC moves per spin
    public int tFlip=0;
    public int currentTime;
  
    double[] g;    // logarithm of the density of states (energy argument// translated by 2N)
    int[] H; // histogram (reduce f when it is "flat")
    int binN;
    double eRangeMultiplier=4;
    double binSize;
    double f=Math.exp(1.0);      // multiplicative modification factor to g
    double eMax;
    double eMin;
    long intMax = 2^31-1;
    
    /**  @param mc - MCSimulation - Monte Carlo simulation class hosting 
    */
    public WangLandauMC(){
        hField =0.95;
        useLongRange=false;
        jInteraction = 1.0;//0.009090909;
        temperature=1.7700844444444443;
        updateParams();
    }
    
    /**  @param mc - MCSimulation - Monte Carlo simulation class hosting 
    */
    public WangLandauMC(MCSimulation mc){
            this(mc,true);
    }
    
    /**  @param mc - MCSimulation - Monte Carlo simulation class hosting 
    *     @param out - boolean to determine if outputting text confirmations into console
    */
    public WangLandauMC(MCSimulation mc, boolean out){
        triggerReady=false;
        initParams(mc);
        
        // determine if outputting
        output = out;
        if(output){System.out.print("Run:"+mc.getRuns()+"     ");}
        run = mc.getRuns();  
        if(mc.getSimSystem().getClass().getName().contains("attice")){
            lattice = (LatticeMagInt)mc.getSimSystem();
        }else{
            System.err.println("THIS CLASSS (WANGLANDAUMC NEEDS A LATTICE AS SYSTEM.");
        }
        magnetization = lattice.getMagnetization();
        magStaggered = lattice.getMagStaggered();
        
        makeEnergyBins();
        g = new double[binN+1];
        H = new int[binN+1];
        
        for(int e = 0; e<=binN; e++) {
            H[e] = 0;
        }
        
        // initialize accumulators and energy
        resetData();
        energy=calcIsingEnergy();
        SimProcessParser parser = new SimProcessParser(null,SimPost,ParamPost);
        tFlip = parser.timeToFlipField();
        
        updateParams();
        if(triggerOn)setTrigger();
        currentTime=0;
    }
    
    private void makeEnergyBins(){
        eMax = Math.abs(2*jInteraction*N);
        eMin = -2*Math.abs(jInteraction*N)-Math.abs(hField*N);
        double eRange = 4*Math.abs(jInteraction*N)+Math.abs(hField*N);
        // Only need integer divisions for zero field 
        if(hField==0){eRangeMultiplier=1;}
        long binNlong = (long) (eRange*eRangeMultiplier);
        
        if(binNlong > intMax){
            System.err.print("More bins than integer max. "
                + "WILL HAVE INDEX RETRIEVAL ERRORS. CHOOSE SMALLER MULTIPLIER");
        }
        binN = (int) binNlong;
        binSize = (double)(eRange/binN);       
    }
    
    private int getEnergyBin(double e){
        double diff = e-eMin;
        return (int) (diff/binSize);
    }
    
    private void initParams(MCSimulation mc){
        param = new ParameterBank(mc.getParamPostFix());
        SimPost = mc.getSimPostFix();
        ParamPost = mc.getParamPostFix();
        jInteraction = param.jInteraction;
        if(useLongRange){
            jInteraction 
                    = (Math.signum(jInteraction) > 0) ? jInteraction*4/((2*param.R+1)*(2*param.R+1)-1): -1*jInteraction*4/((2*param.R+1)*(2*param.R+1)-1);
        }
        hField = param.hField;
        temperature = param.temperature;
        useLongRange = param.useLongRange;
        useDilution = param.useDilution;
        useHeter = param.useHeter;
        L = param.L;
        D = param.D;
        N = (int)Math.pow(L, D);
    }
    
    /**
    *         setTriggerOnOff allows turning triggers off
    * 
    * @param tr - on is true
    */ 
    @Override
    public void setTriggerOnOff(boolean tr){
        triggerOn=tr;
    }

    private void updateParams(){}

    /**
    *        doOneStep should do a single monte carlo step.
    */ 
    @Override
    public void doOneStep(){
        currentTime++;
        int mcsMax = mcs+Math.max(100000/N, 1);
        for(; mcs<mcsMax; mcs++) {
            flipSpins();
        }
        double z = 0;
        for(int e = 0; e<=binN; e++) {
            if(g[e]>0) {
                z += Math.exp(g[e]-g[0])*Math.exp(-(1/temperature)*(eMin+e*binSize));
            }
        }
        if(isFlat()) {
            f = Math.sqrt(f);
            
            for(int e = 0; e<=binN; e++) {
                H[e] = 0;
            }
        }
        magnetization = lattice.getMagnetization();    
        magStaggered = lattice.getMagStaggered();    
        // Wait till after the flip to measure the susceptibility

        if(triggerOn){
            triggerReady = trigger.triggerNow(lattice.getMagnetization(),energy);}
        if(triggerReady){tTrigger = trigger.getTriggerTime();}
    }
    /**
    *         setSeed should set the random number seed to the value given
    * 
    *  @param seed - random number seed
    */ 
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
    *   getSpecificHeat should return the current specific heat.
    * 
    */
    @Override
    public double getSpecificHeat() {
        return(getSpecificHeat(2.2));
    }
    /**
    *   getSpecificHeat should return the current specific heat.
    * 
    */
    public double getSpecificHeat(double temp) {
        return (heatCapacity(N,g,1/temp)/N);
    }
    /**
    *   getSusceptibility should return the current susceptibility
    * 
    */
    @Override
    public double getSusceptibility() {
        return(0);
    }
    /**
    *   getSusceptibilityMult should return the current susceptibility multiple values if possible
    * 
    */
    @Override
    public ArrayList<Double> getSusceptibilityMult() {
        ArrayList<Double> sus = new ArrayList<Double>();
        sus.add(0.0);
        return (sus);
    }
    /**
    *         resetData resets all the accumulators of data and number of mcsteps made.
    */
    public void resetData() {
        mcs = 0;
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
        energy=calcIsingEnergy();
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
        energy=calcIsingEnergy();
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
        trigger.reset2();//reset but keep std or go into intervention mode
        if(t==0){lattice.initialize(param.s);
        if(param.useHeter){lattice.setFixedLatticeValues();}
        }else{
        lattice.setInitialConfig(t,run,"");}
        magnetization = lattice.getMagnetization();
        magStaggered = lattice.getMagStaggered();
        setSeed(seed);
        energy=calcIsingEnergy();
    }
    /**
    *   flipField should flip the magnetic field
    * 
    */
    @Override
    public void flipField(){
        hField = -hField;
        if(noFieldMessages){}else{
            if(output){
                System.out.println("Field Flipped at t:"+currentTime+" now hField:"+hField);
            }}
        updateParams();
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
        temperature = temp;
        updateParams();
    }
    /**
    *         changeH should set the magnetic field to the given value.
    * 
    *  @param hnow - new magnetic field
    */ 
    @Override
    public void changeH(double hnow){
        hField = hnow;
        updateParams();
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
        updateParams();
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
        hField = -hField;
        updateParams();
    }


     /**
    *         setTrigger should initialize the trigger to be used in the simulation
    */ 
    @Override
    public void setTrigger(){
        trigger = param.setProperTrigger(param, trigger, SimPost, output);
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
        SimProcessParser parser = new SimProcessParser(null,SimPost,ParamPost);
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
    *       calcEnergy calculates the energy of the lattice in an Ising Model 
    * 
    * @return energy in lattice for Ising model
    */
    public double calcIsingEnergy(){
        double e=0.0;
        for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++){{
            if(D==3){
                for (int k = 0; k < L; k++){			    		  
                    e += -1.0*param.jInteraction*lattice.getValue(i,j,k)*lattice.getNeighSum(i, j, k)/2
                        -1.0*param.hField*lattice.getValue(i, j, k);   
                }
            // 2d case    
            }else{
                    // divide the bonds term by 2 since they will be double counted
                    e += -1.0*param.jInteraction*lattice.getValue(i,j)*lattice.getNeighSum(i,j)/2
                        -1.0*param.hField*lattice.getValue(i,j);
            }
        }}
        return e;
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
    *   getEnergy should return the current energy.
    * 
    */
    @Override
    public double getEnergy(){
        return energy;
    }
    
    /**
    *       getSystemMeasurer should return the system measurer which contains system data
    *   and accumulators.
    * 
    */
    @Override
    public MeasureIsingSystem getSystemMeasurer(){
        return null;
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
    
    public double getMultiplicativeFactor(){return f;}
    
    private boolean isFlat() {
        int netH = 0;
        double numEnergies = 0;
        int hLast =0;
        int hLargest = 0;
        
        for(int e = 0; e<=binN; e++) {
            if(H[e]>0) {
                netH += H[e];
                hLast = H[e];
                numEnergies++;
                if(H[e]> hLargest){hLargest=H[e];}
            }
        }
        System.out.println("WangLandauMC | f: "+
                getMultiplicativeFactor()+"    Nenergy: "+numEnergies+
                "      <H(E)>: "+(netH/numEnergies)+"   H(0): "+hLast+"    hMax: "+hLargest);            
        for(int e = 0; e<=binN; e++) {
            if((0<H[e])&&(H[e]<0.8*netH/numEnergies)) {
                return false;
            }
        }
        return true;
    }
    
    private void flipSpins() {
        for(int steps = 0; steps<N; steps++) {
            int i = (int) (Ran.nextDouble()*L);
            int j = (int) (Ran.nextDouble()*L);
            int k;
                if(D==3){
                    k = (int) (Ran.nextDouble()*L);
                }else{
                    k=0;
                }
            int currSpin = lattice.getValue(i, j,k);
            int nearestSum = lattice.getNeighSum(i, j,k);
            int delE = (int)(nearestSum*currSpin*jInteraction+currSpin*hField);
            if(Math.random()<Math.exp(g[getEnergyBin(energy)]-g[getEnergyBin(energy+delE)])) {
                lattice.setValue(i, j,k,(-1*currSpin), 0);
                energy += delE;
            }
                g[getEnergyBin(energy)] += Math.log(f);
                H[getEnergyBin(energy)] += 1;
            
        }
    }
    
    
    private double logZ(int N, double[] g, double beta) {
        double m = 0;
        for(int E = 0; E<=binN; E++) {
            m = Math.max(m, g[E]-beta*(eMin+E*binSize));
        }
        double s = 0;
        for(int E = 0; E<=binN; E++) {
            s += Math.exp(g[E]-beta*(eMin+E*binSize)-m);
        }
        return Math.log(s)+m;
    }

    private double heatCapacity(int N, double[] g, double beta) {
        double logZ = logZ(N, g, beta);
        double E_avg = 0;
        double E2_avg = 0;
        for(int E = 0; E<=binN; E++) {
            if(g[E]==0) {
                continue;
            }
            E_avg += E*Math.exp(g[E]-beta*(eMin+E*binSize)-logZ);
            E2_avg += E*E*Math.exp(g[E]-beta*(eMin+E*binSize)-logZ);
        }
        return(E2_avg-E_avg*E_avg)*beta*beta;
    }

    
    // test the class
    public static void main(String[] args) {
        WangLandauMC mc = new WangLandauMC();   
    }
}

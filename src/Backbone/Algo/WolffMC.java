package Backbone.Algo;

/**
 * 
 *   @(#) WolffMC
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
import java.util.HashMap;
import java.util.Random;

/**  
 *   Basic Implementation of Wolff Cluster Algorithm Monte Carlo Simulation.
 *  Single thread implementation of simulation
 *  <br>
 *  @param mc - MCSimulation - Monte Carlo simulation class hosting 
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2012-08    
 */
public class WolffMC implements IsingMC{
    
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
    private int R;
    private int run;
    private boolean output = true;
    private boolean triggerReady=false;
    private boolean triggerOn=true;
    private Random Ran = new Random();
    private MeasureIsingSystem measurer;
    private String SimPost;
    private String ParamPost;  
    private ArrayList<Integer> currCluster;
    private int clustSpinVal;
    private ArrayList<Integer> perimeterSpins;
    private HashMap<Integer,Boolean> visitedSpins;
    private int magnetization = 0;  
    private int magStaggered = 0;
    private double bondProbability;
    private boolean fixedAllowedInCluster = true;
    
    public double energy;  
    public int currentSeed;
    public int mcs = 0; // number of MC moves per spin
    public int tFlip=0;
    public int currentTime;
    
    /**  @param mc - MCSimulation - Monte Carlo simulation class hosting 
    *     @param out - boolean to determine if outputting text confirmations into console
    */
    public WolffMC(MCSimulation mc, boolean out){
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
        N = param.N;
        D = param.D;
        measurer = new MeasureIsingSystem(temperature,N);
        triggerReady=false;
        
        // determine if outputting
        output = out;
        if(output){System.out.print("Run:"+mc.getRuns()+"     ");}
        run = mc.getRuns();  
        if(mc.getSimSystem().getClass().getName().contains("attice")){
            lattice = (LatticeMagInt)mc.getSimSystem();
        }else{
            System.err.println("THIS CLASSS (WOLFFMC NEEDS A LATTICE AS SYSTEM.");
        }
        measurer.setN(N, lattice.getNFixed());
        magnetization = lattice.getMagnetization();
        magStaggered = lattice.getMagStaggered();
        
        // initialize accumulators and energy
        resetData();
        energy=calcIsingEnergy();
        SimProcessParser parser = new SimProcessParser(null,SimPost,ParamPost);
        tFlip = parser.timeToFlipField();
        
        calcBondProb();
        setTrigger();
        currentTime=0;
        
        ArrayList<Integer> currCluster = new ArrayList<Integer>();
        ArrayList<Integer> perimeterSpins = new ArrayList<Integer>();
    }
    
    private void calcBondProb(){
        bondProbability = 1-Math.exp(-2*(jInteraction/(temperature)));
    }
    
    private void seedSpin(){
        boolean attemptAgain = true;
        int x=0,y=0;
        
        // Dont start from fixed spin seed
        while(attemptAgain){
            x = (int)(Ran.nextDouble()*L);
            y = (int)(Ran.nextDouble()*L);
            if(!lattice.isThisFixed(x, y, 0)){attemptAgain=false;}
        }
        
        clustSpinVal = lattice.getValue(x, y);
        currCluster.add((x+y*L));
        visitedSpins.put((x+y*L), true);
        addPerimeterSpins(x,y);
    }
    
    private void addPerimeterSpins(int x, int y){
        if(R==0){
            tryPerimeterAdd(((x+1)%L),y);
            tryPerimeterAdd(((x+L-1)%L),y);
            tryPerimeterAdd(x,((y+1)%L));
            tryPerimeterAdd(x,((y+L-1)%L));
        }
    }
    
    private void tryPerimeterAdd(int x, int y){
        // dont revisit site
        if(previousVisited(x,y)){return;}else{
            visitedSpins.put((x+y*L), true);
        }
        // only add if free and parallel
        if(lattice.getValue(x, y)==clustSpinVal){
            perimeterSpins.add(x+y*L);           
        }
    }
    
    private boolean previousVisited(int x, int y){
        boolean visited = false;
        if(visitedSpins.containsKey(x+y*L)){
            if(visitedSpins.get(x+y*L)){visited = true;}
        }   
        return visited;
    }
    
    private int popFromPerimeter(){
        int head = perimeterSpins.get(0);
        perimeterSpins.set(0,perimeterSpins.get(perimeterSpins.size()));        
        perimeterSpins.remove(perimeterSpins.size());
        return head;
    }
    
    private void bondSpin(int i){
        if(Ran.nextDouble() <= bondProbability){
            int x = i%L;
            int y = (int)(i/L);
            y = y%L;
            currCluster.add((x+y*L));
            visitedSpins.put((x+y*L), true);
            addPerimeterSpins(x, y);
        }
    }
    
    private void flipClusterSpins(){
        int newSpin = -1*clustSpinVal;
        int i0,x,y;
        int nearestSum;
        double dE = 0.0;
        
        for(int u = 0;u < currCluster.size();u++){
            i0 = currCluster.get(u);
            x = i0%L;
            y = ((int)(i0/L))%L;
            // DO NOT FLIP FIXED
            if(!lattice.isThisFixed(x,y, 0)){
                lattice.setValue(x, y, 0, newSpin, 0);
                nearestSum = lattice.getNeighSum(x,y);   
                dE += 2*clustSpinVal*(nearestSum*jInteraction+hField);
            }
        }
        energy += dE;
    }
    
    /**
    *        doOneStep should do a single monte carlo step.
    */ 
    public void doOneStep(){
        currentTime++;
        
        visitedSpins = new HashMap<Integer,Boolean>();
        
        seedSpin();
        int i0;
        while(perimeterSpins.size()>0){
            i0 = popFromPerimeter();
            bondSpin(i0);
        }
        flipClusterSpins();
        
        // clear clusterInfo
        currCluster.clear();
        perimeterSpins.clear();
        visitedSpins.clear();
        
        // Wait till after the flip to measure the susceptibility
        if(currentTime>(tFlip+tAccumulate)){
            updateAccumulators();
        }

        if(triggerOn){
            triggerReady = trigger.triggerNow(lattice.getMagnetization(),energy);}
        if(triggerReady){tTrigger = trigger.getTriggerTime();}
    }
    
    /**
    *       updateAccumulators updates accumulators used to calculate susceptibility
    */
    private void updateAccumulators(){
        measurer.updateE(energy, true);
        // use staggered Magnetization for antiferromagnet
        measurer.updateM((param.jInteraction < 0)?magnetization:magStaggered, true);
        mcs++;
    }
    /**
    *   flipField should flip the magnetic field
    * 
    */
    @Override
    public void flipField(){
        hField = -hField;
        if(suppressFieldMessages){}else{
            if(output){
                System.out.println("Field Flipped at t:"+currentTime+" now hField:"+hField);
            }}
        calcBondProb();
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
        measurer.changeTemp(temp);
        calcBondProb();
    }
    /**
    *         changeH should set the magnetic field to the given value.
    * 
    *  @param hnow - new magnetic field
    */ 
    @Override
    public void changeH(double hnow){
        hField = hnow;
        calcBondProb();
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
        calcBondProb();
    }
    /**
    *         changetTFlipField should flip the field and set the temperature to
    *   the given value
    * 
    *  @param temp - new temperature
    */ 
    @Override
    public void changeTFlipField(double temp){
        temperature = temp;
        measurer.changeTemp(temp);
        hField = -hField;
        calcBondProb();
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
    *         setSeed should set the random number seed to the value given
    * 
    *  @param seed - random number seed
    */ 
    public void setSeed(int seed){
        currentSeed = seed;
        Ran.setSeed(seed);
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
     
    /**
    *   getJinteraction should return the current strength of interaction.
    * 
    */
    @Override
    public double getJinteraction(){
        return jInteraction;
    }
    
    /**
    *   getSeed should return the current random number seed
    * 
    */
    @Override
    public int getSeed(){return currentSeed;}
    /**
    *   getEnergy should return the current energy.
    * 
    */
    @Override
    public double getEnergy(){
        return energy;
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
    *         setTriggerOnOff allows turning triggers off
    * 
    * @param tr - on is true
    */ 
    @Override
    public void setTriggerOnOff(boolean tr){
        triggerOn=tr;
    }   
}

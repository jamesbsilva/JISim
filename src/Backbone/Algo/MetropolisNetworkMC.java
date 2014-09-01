package Backbone.Algo;

/**
 * 
 *      @(#)  MetropolisNetworkMC
 */  

import Backbone.System.Network;
import Backbone.Util.MeasureIsingSystem;
import Backbone.Util.DataSaver;
import Backbone.Util.SimProcessParser;
import Backbone.Util.ParameterBank;
import JISim.MCSimulation;
import Triggers.*;
import java.util.ArrayList;
import java.util.Random;

/**  
 *   Basic Implementation of Metropolis Monte Carlo Simulation.
 *  Single thread implementation of simulation
 *  <br>
 *  @param mc - MCSimulation - Monte Carlo simulation class hosting 
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2012-11    
 */
public final class MetropolisNetworkMC  implements IsingMC {
    private ParameterBank param;
    private Network network;
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
    private boolean output = true;
    private boolean triggerReady=false;
    private boolean triggerOn = true;
    private Random Ran = new Random();
    private MeasureIsingSystem measurer;
    private String SimPost;
    private String ParamPost;  
    
    public double energy;  
    private int magnetization = 0;
    public int currentSeed;
    public int mcs = 0; // number of MC moves per spin
    public int tFlip=0;
    public int currentTime;

    
    /**  @param mc - MCSimulation - Monte Carlo simulation class hosting 
    */
    public MetropolisNetworkMC(){
        hField =0.95;
        useLongRange=false;
        jInteraction = 1.0;//0.009090909;
        temperature=1.7700844444444443;
    }
    
    /**  @param mc - MCSimulation - Monte Carlo simulation class hosting 
    */
    public MetropolisNetworkMC(MCSimulation mc){
        this(mc,true);
    }
    
    /**  @param mc - MCSimulation - Monte Carlo simulation class hosting 
    *     @param out - boolean to determine if outputting text confirmations into console
    */
    public MetropolisNetworkMC(MCSimulation mc, boolean out){
        param = mc.getParams();
        if(param==null){param = new ParameterBank(mc.getParamPostFix());}
        SimPost = mc.getSimPostFix();
        ParamPost = mc.getParamPostFix();        
        hField = param.hField;
        temperature = param.temperature;
        useLongRange = param.useLongRange;
        useDilution = param.useDilution;
        useHeter = param.useHeter;
        L = param.L;
        N = L*L;
        D = param.D;
        
        jInteraction = param.jInteraction;
        if(useLongRange){
            jInteraction 
                    = (Math.signum(jInteraction) > 0) ? jInteraction*4/((2*param.R+1)*(2*param.R+1)-1): -1*jInteraction*4/((2*param.R+1)*(2*param.R+1)-1);
        }
        measurer = new MeasureIsingSystem(temperature,N);
        triggerReady=false;
        
        // determine if outputting
        output = out;
        if(output){System.out.print("Run:"+mc.getRuns()+"     ");}
        run = mc.getRuns();  
        if(mc.getSimSystem().getClass().getName().contains("etwork")){
            network = (Network)mc.getSimSystem();
        }else{
            System.err.println("THIS CLASSS (METROPOLISNETWORKMC NEEDS A NETWORK AS SYSTEM.");
        }
        measurer.setN(N, network.getNFixed());
        magnetization = network.getMagnetization();
        
        // initialize accumulators and energy
        resetData();
        
        energy=calcIsingEnergy();
        SimProcessParser parser = new SimProcessParser(null,SimPost,ParamPost);
        tFlip = parser.timeToFlipField();
        
        System.out.println("Initial M: "+network.getMagnetization()+"    E:"+(energy/(L*L)));
        
        if(triggerOn){setTrigger();}
        currentTime=0;
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
    /**
    *        doOneStep should do a single monte carlo step.
    */ 
    @Override
    public void doOneStep(){
        currentTime++;
        for(int f = 0;f<N;f++) {
            int i = (int) (Ran.nextDouble()*N);
            
            // Check if fixed. If fixed then its been visited but skip
            if (useHeter == true){
                if(network.isThisFixed(i)){continue;}
            }

            int nearestSum;
            nearestSum = network.getNeighSum(i);

            // diluted version
            if(useDilution ==true){
                dilutedStep(i,nearestSum);
            }else{
                // regular Step
                regularStep(i,nearestSum);
            }
        }
        magnetization = network.getMagnetization();
        
        // Wait till after the flip to measure the susceptibility
        if(currentTime>(tFlip+tAccumulate)){
            updateAccumulators();
        }

        if(triggerOn){
            triggerReady = trigger.triggerNow(network.getMagnetization());}
        if(triggerReady){tTrigger = trigger.getTriggerTime();}
    }
    /**
    *       updateAccumulators updates accumulators used to calculate susceptibility
    */
    private void updateAccumulators(){
        measurer.updateE(energy, true);
        measurer.updateM(magnetization, true);
        mcs++;
    }
    /**
    *        dilutedStep performs the main part of the metropolis algorithm step
    *   for a diluted system step.
    * 
    *   @param i - i coordinate
    *   @param nearestSum - the sum of the spins in interaction range 
    */ 
    private void dilutedStep(int i, int nearestSum){
        int currentSpin = network.getValue(i);
        double acceptance;

        int r = (int) (Ran.nextDouble()*2);    
        if(currentSpin == 0){
        // change to +1
        if(r == 0){
            double dE = -currentSpin*(nearestSum*jInteraction+hField);
            acceptance = Math.exp(-1.0*dE/temperature);
            if((dE<=0)||(acceptance>Ran.nextDouble())) {
            network.setValue(i,1,currentTime);
            //---------------add 3D

            energy += dE;
            }}
            // change to -1
            else{
            double dE = currentSpin*(nearestSum*jInteraction+hField);
            acceptance = Math.exp(-1.0*dE/temperature);
            if((dE<=0)||(acceptance>Ran.nextDouble())) {
                network.setValue(i,-1,currentTime);
            //---------------add 3D

            energy += dE;
            }
            }
        }
        else if(currentSpin == 1){
        // change to - 1
            if(r == 0){
            double dE = 2*currentSpin*(nearestSum*jInteraction+hField);
            acceptance = Math.exp(-1.0*dE/temperature);
            if((dE<=0)||(acceptance>Ran.nextDouble())) {
                network.setValue(i,-1,currentTime);
                //---------------add 3D


                energy += dE;
            }}
            // change to 0
            else{
            double dE = currentSpin*(nearestSum*jInteraction+hField);
            acceptance = Math.exp(-1.0*dE/temperature);
            if((dE<=0)||(acceptance>Ran.nextDouble())) {
                network.setValue(i, 0,currentTime);
            //---------------add 3D

            energy += dE;
                }
            }  
        }
        // current Spin is -1
        else{
        // change to + 1
        if(r == 0){
            double dE = 2*(nearestSum*currentSpin*jInteraction+currentSpin*hField);
            acceptance = Math.exp(-1.0*dE/temperature);
            if((dE<=0)||(acceptance>Ran.nextDouble())) {
                network.setValue(i, 1,currentTime);
                
                energy += dE;
            }}
            // change to 0
            else{
            double dE = (nearestSum*currentSpin*jInteraction+currentSpin*hField);
            acceptance = Math.exp(-1.0*dE/temperature);
            if((dE<=0)||(acceptance>Ran.nextDouble())) {
                network.setValue(i, 0,currentTime);
            //---------------add 3D

            energy += dE;}}  
        }
    }
    /**
    *        regularStep performs the main part of the metropolis algorithm step
    *   for a regular non-diluted system step.
    * 
    *   @param i - i coordinate
    *   @param nearestSum - the sum of the spins in interaction range 
    */ 
    private void regularStep(int i, int nearestSum){
        int currentSpin = network.getValue(i);
        double acceptance;
        double dE = 2*currentSpin*(nearestSum*jInteraction+hField);
        acceptance = Math.exp(-1.0*dE/temperature);
        double toss = Ran.nextDouble();
        if( (dE <= 0) || (acceptance > toss) ){
            int newSpin = -currentSpin;
            if(useHeter){
                if(network.isThisFixed(i)){System.out.println("MetropolisNetworkMC | Changing Fixed");}
            }
            network.setValue(i,newSpin,currentTime);
            energy = energy+dE;
        }
    }
    /**
    *        isPostitive returns true if number is positive
    * 
    *   @param number - input number
    */ 
    public int isPositive(int number){
        if (number >= 0){return 1;}else{return 0;}
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
        network.setNetworkSeed(network.getNetworkSeed());
        network.initialize(param.s);
        magnetization = network.getMagnetization();
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
        network.setNetworkSeed(network.getNetworkSeed());
        network.initialize(param.s);
        if(param.useHeter){network.setFixedNetworkValues();}
        magnetization = network.getMagnetization();
        setSeed(seed);
        System.out.println("Reset Simulation Seed - Current Seed : "+currentSeed 
                + "     M: "+magnetization+"       links: "+network.getTotalLinks());
        energy=calcIsingEnergy();
    }
    /**
    *         resetSimulation should reset all simulation parameters, set
    *   all network values to network at the time given, and set the 
    *   seed to the given value
    * 
    *  @param t - time to set network to
    *  @param seed - random number seed
    */ 
    @Override
    public void resetSimulation(int t, int seed){
        currentTime = t;
        resetData();
        //change to allow more triggers
        trigger.reset2();//reset but keep std or go into intervention mode
        if(t==0){
            network.initialize(param.s);
            if(param.useHeter){network.setFixedNetworkValues();}
        }else{
            network.setInitialConfig(t,run,"");
            if(param.useHeter){network.setFixedNetworkValues();}
        }
        magnetization = network.getMagnetization();
        setSeed(seed);
        System.out.println("Reset Simulation Seed+ T- Current Seed : "+currentSeed 
        + "     M: "+magnetization+"       links: "+network.getTotalLinks());
        energy=calcIsingEnergy();
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
                System.out.println("MetropolisNetworkMC | Field Flipped at t:"+currentTime+" now hField:"+hField);
            }}
        
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
    }
    /**
    *         changeH should set the magnetic field to the given value.
    * 
    *  @param hnow - new magnetic field
    */ 
    @Override
    public void changeH(double hnow){
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
        hField = hnow;
        temperature = temp;
        measurer.changeTemp(temp);
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
        DataSaver save = new DataSaver(network.getInstId());
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
            if(t>tInitial){save.saveConfig(network,run,t);}
        }
    }

    /**
    *   getSimSystem should return the network object
    * 
    */
    @Override
    public Network getSimSystem(){return network;}
    /**
    *   getM should return the current magnetization.
    * 
    */
    @Override
    public int getM(){return magnetization;}
    /**
    *   getMStag should return the current regular magnetization for network.
    * 
    */
    public int getMStag(){return magnetization;}
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
    *       calcEnergy calculates the energy of the network in an Ising Model 
    * 
    * @return energy in network for Ising model
    */
    public double calcIsingEnergy(){
        double e=0.0;
        
        for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++){{
            if(D==3){
                for (int k = 0; k < L; k++){			    		  
                    e += -1.0*jInteraction*network.getValue(i)*network.getNeighSum(i)/2
                        -1.0*hField*network.getValue(i);   
                }
            // 2d case    
            }else{
                    // divide the bonds term by 2 since they will be double counted
                    e += -1.0*jInteraction*network.getValue(i)*network.getNeighSum(i)/2
                        -1.0*hField*network.getValue(i);
            }
        }}
        return e;
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
    
    // test the class
    public static void main(String[] args) {
        MetropolisNetworkMC mc = new MetropolisNetworkMC();
    }
}

package Backbone.Algo;

/**
*      @(#)  MetropolisMC
*/  

import Backbone.System.LatticeMagInt;
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
* 
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2011-09    
*/
public final class MetropolisMC  implements IsingMC {
    private LatticeMagInt lattice;
    private double[][][] w;
    private ParameterBank param;
    private double temperature; private double jInteraction;
    private boolean useLongRange; private boolean useDilution;
    private double hField;
    private boolean useHeter;
    private Trigger trigger;
    private int tTrigger =0;
    private int tAccumulate =1000; // begin accumulator
    private boolean suppressFieldMessages=false;
    private int L; private int N;
    private int D; private int Q;
    private int run;
    private boolean output = true;
    private boolean triggerReady=false;
    private boolean triggerOn = true;
    private Random Ran = new Random();
    private MeasureIsingSystem measurer;
    private String SimPost;
    private String ParamPost;  
    
    public double energy;  
    private int magStaggered = 0; private int magnetization = 0;
    public int currentSeed;
    public int mcs = 0; // number of MC moves per spin
    public int tFlip=0;
    public int currentTime;

    
    public MetropolisMC(){
            hField =0.95;
            useLongRange=false;
            jInteraction = 1.0;//0.009090909;
            temperature=1.7700844444444443;
            calcW();
            printW();
    }
    
    /** 
    * @param mc - MCSimulation - Monte Carlo simulation class hosting 
    */
    public MetropolisMC(MCSimulation mc){
            this(mc,true);
    }
    
    /**  
    *    @param mc - MCSimulation - Monte Carlo simulation class hosting 
    *    @param out - boolean to determine if outputting text confirmations into console
    */
    public MetropolisMC(MCSimulation mc, boolean out){
        param = mc.getParams();
        if(param==null){param = new ParameterBank(mc.getParamPostFix());}
        SimPost = mc.getSimPostFix();
        ParamPost = mc.getParamPostFix();        
        hField = param.hField; temperature = param.temperature;
        useLongRange = param.useLongRange; useDilution = param.useDilution;
        useHeter = param.useHeter;
        L = param.L; N = param.N; D = param.D;
        
        jInteraction = param.jInteraction;
        if(param.R > 1){
            jInteraction 
                    = (Math.signum(jInteraction) > 0) ? jInteraction*4/((2*param.R+1)*(2*param.R+1)-1): -1*jInteraction*4/((2*param.R+1)*(2*param.R+1)-1);
        }
        measurer = new MeasureIsingSystem(temperature,N);
        triggerReady=false;
        
        // determine if outputting
        output = out;
        if(output){System.out.print("Run:"+mc.getRuns()+"     ");}
        run = mc.getRuns();
        if(mc.getSimSystem().getClass().getName().contains("attice")){
            lattice = (LatticeMagInt)mc.getSimSystem();
        }else{
            System.err.println("THIS CLASSS (METROPOLISMC NEEDS A LATTICE AS SYSTEM.");
        }
        measurer.setN(N, lattice.getNFixed());
        magnetization = lattice.getMagnetization(); magStaggered = lattice.getMagStaggered();
        
        // initialize accumulators and energy
        resetData();
        calcW();
        energy=calcIsingEnergy();
        SimProcessParser parser = new SimProcessParser(null,SimPost,ParamPost);
        tFlip = parser.timeToFlipField();
        
        System.out.println("Initial M: "+((jInteraction > 0)?lattice.getMagnetization(): lattice.getMagStaggered())+"    E:"+(energy/(L*L)));
        
        if(triggerOn){setTrigger();}
        currentTime=0;
    }

    /**
    *         setTriggerOnOff allows turning triggers off
    * 
    * @param tr - on is true
    */ 
    @Override
    public void setTriggerOnOff(boolean tr){ triggerOn=tr; }

    /**
    *         calcW calculates the transition matrix that is necessary
    *   for the metropolis to get the statistics with the right temperature
    *   and field. Supports Long Range and Nearest Neighbor
    * 
    */ 
    public void calcW(){
        jInteraction = param.jInteraction;
        if (useDilution == false && useLongRange == false){
            if(suppressFieldMessages){}else{
                if(output){
                System.out.println("MetropolisMC | Calc W (Nearest Neighbor) with j:"+jInteraction
                    +"      h:"+hField+"    temp:"+temperature);}
            }

            // First index 4 due to 0,1 for spin sign and 2 and 3 for current spin = 0;
            w = new double[4][2][5]; 
            Q=4;
            w[1][1][4] = Math.exp((-8*jInteraction-2*hField)/temperature); 
            w[0][1][4] = Math.exp((8*jInteraction+2*hField)/temperature);
            w[1][1][3] = Math.exp((-6*jInteraction-2*hField)/temperature);
            w[0][1][3] = Math.exp((6*jInteraction+2*hField)/temperature);
            w[1][1][2] = Math.exp((-4*jInteraction-2*hField)/temperature);
            w[0][1][2] = Math.exp((4*jInteraction+2*hField)/temperature);
            w[1][1][1] = Math.exp((-2*jInteraction-2*hField)/temperature);
            w[0][1][1] = Math.exp((2*jInteraction+2*hField)/temperature);
            w[1][1][0] = Math.exp((-2*hField)/temperature);
            w[0][1][0] = Math.exp((2*hField)/temperature);
            w[1][0][1] = Math.exp((2*jInteraction-2*hField)/temperature);
            w[0][0][1] = Math.exp((-2*jInteraction+2*hField)/temperature);
            w[1][0][2] = Math.exp((4*jInteraction-2*hField)/temperature);
            w[0][0][2] = Math.exp((-4*jInteraction+2*hField)/temperature);
            w[1][0][3] = Math.exp((6*jInteraction-2*hField)/temperature);
            w[0][0][3] = Math.exp((-6*jInteraction+2*hField)/temperature);
            w[1][0][4] = Math.exp((8*jInteraction-2*hField)/temperature);
            w[0][0][4] = Math.exp((-8*jInteraction+2*hField)/temperature);
        }else if(useDilution == true && useLongRange == false){
            w = new double[4][2][5]; 
            Q=4;
            w[1][1][4] = Math.exp((-8*jInteraction-2*hField)/temperature); 
            w[0][1][4] = Math.exp((8*jInteraction+2*hField)/temperature);
            w[1][1][3] = Math.exp((-6*jInteraction-2*hField)/temperature);
            w[0][1][3] = Math.exp((6*jInteraction+2*hField)/temperature);
            w[1][1][2] = Math.exp((-4*jInteraction-2*hField)/temperature);
            w[0][1][2] = Math.exp((4*jInteraction+2*hField)/temperature);
            w[1][1][1] = Math.exp((-2*jInteraction-2*hField)/temperature);
            w[0][1][1] = Math.exp((2*jInteraction+2*hField)/temperature);
            w[1][1][0] = Math.exp((-2*hField)/temperature);
            w[0][1][0] = Math.exp((2*hField)/temperature);
            w[1][0][1] = Math.exp((2*jInteraction-2*hField)/temperature);
            w[0][0][1] = Math.exp((-2*jInteraction+2*hField)/temperature);
            w[1][0][2] = Math.exp((4*jInteraction-2*hField)/temperature);
            w[0][0][2] = Math.exp((-4*jInteraction+2*hField)/temperature);
            w[1][0][3] = Math.exp((6*jInteraction-2*hField)/temperature);
            w[0][0][3] = Math.exp((-6*jInteraction+2*hField)/temperature);
            w[1][0][4] = Math.exp((8*jInteraction-2*hField)/temperature);
            w[0][0][4] = Math.exp((-8*jInteraction+2*hField)/temperature);
            // 2 and 3 are both current spin 0 changed to - and + respectively
            w[3][0][4] = Math.exp((-4*jInteraction-hField)/temperature); 
            w[2][1][4] = Math.exp((4*jInteraction+hField)/temperature);
            w[3][0][3] = Math.exp((-3*jInteraction-hField)/temperature);
            w[2][1][3] = Math.exp((3*jInteraction+hField)/temperature);
            w[3][0][2] = Math.exp((-2*jInteraction-hField)/temperature);
            w[2][1][2] = Math.exp((2*jInteraction+hField)/temperature);
            w[3][0][1] = Math.exp((-jInteraction-hField)/temperature);
            w[2][1][1] = Math.exp((jInteraction+hField)/temperature);
            w[3][0][0] = Math.exp((-hField)/temperature);
            w[2][1][0] = Math.exp((hField)/temperature);
        }else{
            Q = lattice.getNinRange();
            jInteraction = (Math.signum(param.jInteraction) > 0) ? 4.0/Q : -4.0/Q;

            if(!suppressFieldMessages){
                if(output){System.out.println("MetropolisMC | Long R Calc W with j:"+jInteraction+
                "      h:"+hField+"    temp:"+temperature
                + "   neighborsInRange:"+Q);}
            }

            w = new double[4][2][Q+2];

            for(int i = 0;i<4;i++){ for(int j = 0;j<2;j++){
                for(int k = 0;k<Q+1;k++){
                    int jpre,hpre;
                    if((i+j) % 2 == 0){jpre = -1;}else{jpre=1;}
                    if((i) % 2 == 0){hpre = 1;}else{hpre=-1;}
                    double en = 2*(k*jpre*jInteraction+hpre*hField);
                    w[i][j][k] = Math.exp((en)/temperature); 	
                }
            }}
        }
    }

    /**
    *        doOneStep should do a single monte carlo step.
    */ 
    @Override
    public void doOneStep(){
        currentTime++;
        //if(currentTime%100==0){System.out.print("\r time:"+currentTime+"       jInt: "+jInteraction);}
        for(int f = 0;f<N;f++) {
            int i = (int) (Ran.nextDouble()*L); int j = (int) (Ran.nextDouble()*L);
            int k = ( D == 3 ) ? (int) (Ran.nextDouble()*L) : 0;

            // Check if fixed. If fixed then its been visited but skip
            if (useHeter == true){ if(lattice.isThisFixed(i,j,k)){continue;} }

            int nearestSum = (D == 3) ? lattice.getNeighSum(i,j,k) : lattice.getNeighSum(i, j) ;
            
            // diluted version
            if(useDilution == true){
                if(D==3){dilutedStep(i,j,k,nearestSum);}else{dilutedStep(i,j,0,nearestSum);}
            }else{
                // regular Step
                if(D==3){regularStep(i,j,k,nearestSum);}else{regularStep(i,j,0,nearestSum);}
            }
        }
        magnetization = lattice.getMagnetization(); magStaggered = lattice.getMagStaggered();
        
        // Wait till after the flip to measure the susceptibility
        if(currentTime>(tFlip+tAccumulate)){ updateAccumulators(); }
        if(triggerOn){ triggerReady = trigger.triggerNow(lattice.getMagnetization(),energy);}
        if(triggerReady){tTrigger = trigger.getTriggerTime();}
    }

    /**
    *       updateAccumulators updates accumulators used to calculate susceptibility
    */
    private void updateAccumulators(){
        measurer.updateE(energy, true);
        // use staggered magnetization for antiferromagnet
        measurer.updateM((param.jInteraction > 0) ? magnetization : magStaggered, true);
        mcs++;
    }

    /**
    *        dilutedStep performs the main part of the metropolis algorithm step
    *   for a diluted system step.
    * 
    *   @param i - i coordinate
    *   @param j - j coordinate
    *   @param k - k coordinate
    *   @param nearestSum - the sum of the spins in interaction range 
    */ 
    private void dilutedStep(int i, int j, int k, int nearestSum){
        int currentSpin = lattice.getValue(i,j);
        double acceptance;
        int r = (int) (Ran.nextDouble()*2);    
        if(currentSpin == 0){
            // change to +1
            if(r == 0){
                acceptance = w[3][isPositive(nearestSum)][Math.abs(nearestSum)];
                double dE = -currentSpin*(nearestSum*jInteraction+hField);
                if((dE<=0)||(acceptance>Ran.nextDouble())) {
                    lattice.setValue(i, j,k, 1,currentTime);
                    //---------------add 3D
                    energy += dE;
                }
            // change to -1
            }else{
                acceptance = w[2][isPositive(nearestSum)][Math.abs(nearestSum)];
                double dE = currentSpin*(nearestSum*jInteraction+hField);
                if((dE<=0)||(acceptance>Ran.nextDouble())) {
                    lattice.setValue(i, j,k,-1,currentTime);
                    //---------------add 3D
                    energy += dE;
                }
            }
        }else if(currentSpin == 1){
            // change to - 1
            if(r == 0){
                acceptance = w[isPositive(currentSpin)][isPositive(nearestSum)][Math.abs(nearestSum)];
                double dE = 2*currentSpin*(nearestSum*jInteraction+hField);
                if((dE<=0)||(acceptance>Ran.nextDouble())) {
                    lattice.setValue(i,j,k,-1,currentTime);
                    //---------------add 3D
                    energy += dE;
                }
            // change to 0
            }else{
                acceptance = w[isPositive(currentSpin)][isPositive(nearestSum)][Math.abs(nearestSum)];
                double dE = currentSpin*(nearestSum*jInteraction+hField);
                if((dE<=0)||(acceptance>Ran.nextDouble())) {
                    lattice.setValue(i, j,k, 0,currentTime);
                    //---------------add 3D
                    energy += dE;
                }
            }  
        }
        // current Spin is -1
        else{
            // change to + 1
            if(r == 0){
                acceptance = w[isPositive(currentSpin)][isPositive(nearestSum)][Math.abs(nearestSum)];
                double dE = 2*(nearestSum*currentSpin*jInteraction+currentSpin*hField);
                if((dE<=0)||(acceptance>Ran.nextDouble())) {
                    lattice.setValue(i, j,k, 1,currentTime);
                    energy += dE;
                }
            // change to 0
            }else{
                acceptance = w[isPositive(currentSpin)][isPositive(nearestSum)][Math.abs(nearestSum)];
                double dE = (nearestSum*currentSpin*jInteraction+currentSpin*hField);
                if((dE<=0)||(acceptance>Ran.nextDouble())) {
                    lattice.setValue(i, j,k, 0,currentTime);
                    //---------------add 3D
                    energy += dE;
                }
            }  
        }
    }

    /**
    *        regularStep performs the main part of the metropolis algorithm step
    *   for a regular non-diluted system step.
    * 
    *   @param i - i coordinate
    *   @param j - j coordinate
    *   @param k - k coordinate
    *   @param nearestSum - the sum of the spins in interaction range 
    */ 
    private void regularStep(int i, int j, int k, int nearestSum){
        int currentSpin = lattice.getValue(i,j);
        double acceptance;
        double dE = 2*currentSpin*(nearestSum*jInteraction+hField);
        //if(currentTime > 600)System.out.println("Nearest Sum : "+nearestSum);
        if(param.Geo != 2 || param.Geo != 4 ){
            //System.out.println("non-square geo"+"--hfield: "+hField);
            acceptance = Math.exp(-dE/temperature);
        }else{
            acceptance = w[isPositive(currentSpin)][isPositive(nearestSum)][Math.abs(nearestSum)];
        }
        double toss = Ran.nextDouble();
        
        if( (dE<=0) || (acceptance>toss) ){
            int newSpin = -currentSpin;
            if(useHeter){
                if(lattice.isThisFixed(i,j,k)){System.out.println("MetropolisMC | Changing Fixed");}
            }
            lattice.setValue(i,j,k,newSpin,currentTime);
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
    *       printW prints the transition matrix
    */
    public void printW(){
        for(int k = 0;k<Q+1;k++){
            for(int j = 0;j<2;j++){
                for(int i = 0;i<4;i++){
                    if(w[i][j][k]!=0)System.out.println("MetropolisMC | W["+(i+j*4+k*8)+"] = "+w[i][j][k]);
                }   
            }   
        }
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
        currentTime = 0; resetData();
        trigger.reset();
        lattice.initialize(param.s);
        magnetization = lattice.getMagnetization(); magStaggered = lattice.getMagStaggered();
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
        setSeed(seed);
        currentTime = 0; resetData();
        trigger.reset(); lattice.initialize(param.s);
        if(param.useHeter){lattice.setFixedLatticeValues();}
        magnetization = lattice.getMagnetization(); magStaggered = lattice.getMagStaggered();
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
        setSeed(seed); currentTime = t; resetData();
        //change to allow more triggers
        trigger.reset2();//reset but keep std or go into intervention mode
        if( t == 0 ){
            lattice.initialize(param.s);
            if( param.useHeter ){ lattice.setFixedLatticeValues(); }
        }else{
            lattice.setInitialConfig(t,run,"");
        }
        magnetization = lattice.getMagnetization(); magStaggered = lattice.getMagStaggered();
        energy=calcIsingEnergy();
    }
    /**
    *   flipField should flip the magnetic field
    * 
    */
    @Override
    public void flipField(){
        hField = -hField;
        if( !suppressFieldMessages ){
            if(output){ System.out.println("MetropolisMC | Field Flipped at t:"+currentTime+" now hField:"+hField); }
        }
        calcW();
    }
    /**
    *   alignField should align the magnetic field with current state 
    */
    @Override
    public void alignField(){
        if( Math.signum(hField) != Math.signum(magnetization) ){ flipField(); }
    }
    
    /**
    *         changeT should set the temperature to the given value.
    * 
    *  @param temp - temperature to set simulation to
    */ 
    @Override
    public void changeT(double temp){
        measurer.changeTemp(temp);
        temperature = temp; calcW();
    }
    /**
    *         changeH should set the magnetic field to the given value.
    * 
    *  @param hnow - new magnetic field
    */ 
    @Override
    public void changeH(double hnow){
        hField = hnow; calcW();
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
        hField = hnow; temperature = temp;
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
        temperature = temp; hField = -hField;
        measurer.changeTemp(temp);
        calcW();
    }


    /**
    *         setTrigger should initialize the trigger to be used in the simulation
    */ 
    @Override
    public void setTrigger(){ trigger = param.setProperTrigger(param, trigger, SimPost, output); }


    /**
    *   nucleated should return true if nucleation has occurred
    * 
    */
    @Override
    public boolean nucleated(){ return triggerReady; }

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
        for(int t = 0; t < tFinal; t++){
            //Need to flip at step tflip not tflip+1
            if(currentTime == tProcess){
                parser.simProcess(this);
                tProcess = parser.nextSimProcessTime();
            }
            doOneStep();
            if( t > tInitial ){ save.saveConfig(lattice,run,t); }
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
        for (int i = 0; i < L; i++){ for (int j = 0; j < L; j++){
            if(D==3){
                for (int k = 0; k < L; k++){			    		  
                    e += -1.0*jInteraction*lattice.getValue(i,j,k)*lattice.getNeighSum(i, j, k)/2
                        -1.0*hField*lattice.getValue(i, j, k);   
                }
            // 2d case    
            }else{
                // divide the bonds term by 2 since they will be double counted
                e += -1.0*jInteraction*lattice.getValue(i,j)*lattice.getNeighSum(i,j)/2
                    -1.0*hField*lattice.getValue(i,j);
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
        MetropolisMC mc = new MetropolisMC();
    }
}

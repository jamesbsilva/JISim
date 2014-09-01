package Backbone.Algo;

/**
 * 
 *      @(#)  MetropolisFixedStrengthMC
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
 *   Basic Implementation of Metropolis Monte Carlo Simulation with variable fixed 
 *  spin strength interaction.
 * 
 *  Single thread implementation of simulation
 *  <br>
 *  @param mc - MCSimulation - Monte Carlo simulation class hosting 
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2011-09    
 */
public final class MetropolisFixedStrengthMC  implements IsingMC {
    private double[][][] w;
    private ParameterBank param;
    private LatticeMagInt lattice;
    private double temperature;
    private double jInteraction;
    private double jInteractionFixedSpin = 1.5;
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
    private int[] fixedSum;
    private boolean output = true;
    private boolean triggerReady=false;
    private boolean triggerOn=true;
    private Random Ran = new Random();
    private MeasureIsingSystem measurer;
    private String SimPost;
    private String ParamPost;  
    
    public double energy;  
    private int magnetization = 0;
    private int magStaggered = 0;
    public int currentSeed;
    public int mcs = 0; // number of MC moves per spin
    public int tFlip=0;
    public int currentTime;

    
    /**  @param mc - MCSimulation - Monte Carlo simulation class hosting 
    */
    public MetropolisFixedStrengthMC(){
        hField =0.95;
        useLongRange=false;
        jInteraction = 1.0;//0.009090909;
        temperature=1.7700844444444443;
        printW();
    }
    
    /**  @param mc - MCSimulation - Monte Carlo simulation class hosting 
    */
    public MetropolisFixedStrengthMC(MCSimulation mc){
            this(mc,true);
    }
    
    /**  @param mc - MCSimulation - Monte Carlo simulation class hosting 
    *     @param out - boolean to determine if outputting text confirmations into console
    */
    public MetropolisFixedStrengthMC(MCSimulation mc, boolean out){
        param = mc.getParams();
        if(param==null){param = new ParameterBank(mc.getParamPostFix());}
        SimPost = mc.getSimPostFix();
        ParamPost = mc.getParamPostFix();
        jInteraction = param.jInteraction;
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
            System.err.println("THIS CLASSS (METROPOLISMCFIXEDSTRENGTH NEEDS A LATTICE AS SYSTEM.");
        }        
        measurer.setN(N, lattice.getNFixed());
        magnetization = lattice.getMagnetization();
        magStaggered = lattice.getMagStaggered();
        
        // deal with interaction scaling  
        if(param.R != 0){
            createFixedSumArrayBoxR();
        }else{
            createFixedSumArrayNN();
        }
        properlyNormalizeJ();
        
        // initialize accumulators and energy
        resetData();
        energy=calcIsingEnergy();
        SimProcessParser parser = new SimProcessParser(null,SimPost,ParamPost);
        tFlip = parser.timeToFlipField();
        System.out.println("Initial M: "+lattice.getMagnetization()+"    E:"+(energy/(L*L)));
        
        setTrigger();
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

    private void createFixedSumArrayBoxR(){
        int x;int y;int sum=0;
        fixedSum = new int[L*L];
        for(int i=0;i<L;i++){for(int j=0;j<L;j++){            
            sum =0;
            for(int u=0;u<(2*param.R+1);u++){for(int v=0;v<(2*param.R+1);v++){
                x = (i-param.R+u+L)%L;
                y = (j-param.R+v+L)%L;
                if(lattice.isThisFixed(x, y, 0)){
                    sum += lattice.getValue(x, y);
                }
            }}
            fixedSum[i+j*L] = sum;
        }}
    }
    
    
    private void createFixedSumArrayNN(){
        int x;int y;int sum=0;
        fixedSum = new int[L*L];
        for(int i=0;i<L;i++){for(int j=0;j<L;j++){            
            sum =0;
            x = (i+1)%L;
            y = j;
            if(lattice.isThisFixed(x, y, 0)){
                sum += lattice.getValue(x, y);
            }
            x = (i-1+L)%L;
            y = j;
            if(lattice.isThisFixed(x, y, 0)){
                sum += lattice.getValue(x, y);
            }
            x = i;
            y = (j+1)%L;
            if(lattice.isThisFixed(x, y, 0)){
                sum += lattice.getValue(x, y);
            }
            x = i;
            y = (j-1+L)%L;
            if(lattice.isThisFixed(x, y, 0)){
                sum += lattice.getValue(x, y);
            }
            fixedSum[i+j*L] = sum;
        }}
    }
    
    private void properlyNormalizeJ(){
        if(param.R == 0){Q=4;}else{Q = (int)Math.pow((2*param.R+1),2.0)-1;}
        double sum = 0;
        for(int i=0;i<(L*L);i++){
            sum += (double)(fixedSum[i]/Q);  
        }
        sum = sum/(L*L);
        double r =sum;
        double jIntOld = jInteraction;
        jInteraction = 4*jInteraction/((jInteraction+(jInteractionFixedSpin-jInteraction)*r)*Q);
        jInteractionFixedSpin = jInteractionFixedSpin*jInteraction/jIntOld;
        System.out.println("j: "+jInteraction+"         Fixed Spin j: "+jInteractionFixedSpin);
    }
    
    private int getFixSum(int i , int j){
        return fixedSum[i+j*L];
    }
    

    /**
    *        doOneStep should do a single monte carlo step.
    */ 
    @Override
    public void doOneStep(){

        currentTime++;

        //if(currentTime%100==0){System.out.print("\r time:"+currentTime+"       ");}

        for(int f = 0;f<N;f++) {
            int i = (int) (Ran.nextDouble()*L);
            int j = (int) (Ran.nextDouble()*L);
            int k;
            if(D==3){k=(int) (Ran.nextDouble()*L);}else{k=0;}

            // Check if fixed. If fixed then its been visited but skip
            if (useHeter == true){
                if(lattice.isThisFixed(i,j,k)){continue;}
            }

            int nearestSum;

            if (D==3){
                
                nearestSum = lattice.getNeighSum(i,j,k);
            }else{

                nearestSum = lattice.getNeighSum(i, j);}   


            // diluted version
            if(useDilution ==true){
                if(D==3){dilutedStep(i,j,k,nearestSum);}else{dilutedStep(i,j,0,nearestSum);}

            }else{
                // regular Step
                if(D==3){regularStep(i,j,k,nearestSum);}else{regularStep(i,j,0,nearestSum);}
            }
        }
        magnetization = lattice.getMagnetization();
        magStaggered = lattice.getMagStaggered();
        
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
        measurer.updateM((param.jInteraction < 0)?magnetization:magStaggered, true);
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
            }}
            // change to -1
            else{
            acceptance = w[2][isPositive(nearestSum)][Math.abs(nearestSum)];
            double dE = currentSpin*(nearestSum*jInteraction+hField);
            if((dE<=0)||(acceptance>Ran.nextDouble())) {
                lattice.setValue(i, j,k,-1,currentTime);
            //---------------add 3D

            energy += dE;
            }
            }
        }
        else if(currentSpin == 1){
        // change to - 1
            if(r == 0){
            acceptance = w[isPositive(currentSpin)][isPositive(nearestSum)][Math.abs(nearestSum)];
            double dE = 2*currentSpin*(nearestSum*jInteraction+hField);
            if((dE<=0)||(acceptance>Ran.nextDouble())) {
                lattice.setValue(i,j,k,-1,currentTime);
                //---------------add 3D


                energy += dE;
            }}
            // change to 0
            else{
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
            }}
            // change to 0
            else{
            acceptance = w[isPositive(currentSpin)][isPositive(nearestSum)][Math.abs(nearestSum)];
            double dE = (nearestSum*currentSpin*jInteraction+currentSpin*hField);
            if((dE<=0)||(acceptance>Ran.nextDouble())) {
                lattice.setValue(i, j,k, 0,currentTime);
            //---------------add 3D

            energy += dE;}}  
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
        double dE = 2*currentSpin*(nearestSum*jInteraction+hField+(jInteractionFixedSpin-jInteraction)*getFixSum(i, j));
        acceptance = Math.exp(-1*dE/param.temperature);
        
        double toss = Ran.nextDouble();
        
        if( (dE<=0) || (acceptance>toss) ){
            int newSpin = -currentSpin;
            if(useHeter){
                if(lattice.isThisFixed(i,j,k)){System.out.println("MetropolisFixedStrengthCL | Changing Fixed");}
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
             if(w[i][j][k]!=0)System.out.println("MetropolisFixedStrengthCL | W["+(i+j*4+k*8)+"] = "+w[i][j][k]);
        }}}
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
        if(suppressFieldMessages){}else{
            if(output){
                System.out.println("MetropolisFixedStrengthCL | Field Flipped at t:"+currentTime+" now hField:"+hField);
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
        if(suppressFieldMessages){}else{
            if(output){
                System.out.println("MetropolisFixedStrengthCL | Change T at t:"+currentTime+" now temp:"+temp);
            }}
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
        if(suppressFieldMessages){}else{
            if(output){
                System.out.println("MetropolisFixedStrengthCL | "
                        + "Change H at t:"+currentTime+"     h: "+hField);
            }
        }
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
        if(suppressFieldMessages){}else{
            if(output){
                System.out.println("MetropolisFixedStrengthCL | "
                        + "Change T and H at t:"+currentTime
                        +" now temp:"+temp+"     h: "+hField);
            }
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
        temperature = temp;
        measurer.changeTemp(temp);
        hField = -hField;
        if(suppressFieldMessages){}else{
            if(output){
                System.out.println("MetropolisFixedStrengthCL | "
                        + "Change T and Flip Field at t:"+currentTime
                        +" now temp:"+temp+"     h: "+hField);
            }
        }
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
                    e += -1.0*jInteraction*lattice.getValue(i,j,k)*lattice.getNeighSum(i, j, k)/2
                        -1.0*hField*lattice.getValue(i, j, k)
                            -(1.0/2.0)*(jInteractionFixedSpin-jInteraction)*getFixSum(i, j)*lattice.getValue(i, j);   
                }
            // 2d case    
            }else{
                    // divide the bonds term by 2 since they will be double counted
                    e += -1.0*jInteraction*lattice.getValue(i,j)*lattice.getNeighSum(i,j)/2.0
                        -1.0*hField*lattice.getValue(i,j)
                            -(1.0/2.0)*(jInteractionFixedSpin-jInteraction)*getFixSum(i, j)*lattice.getValue(i, j);
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
        MetropolisFixedStrengthMC mc = new MetropolisFixedStrengthMC();
    }
}

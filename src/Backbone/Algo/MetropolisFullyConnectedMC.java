/**
 * 
 *   @(#) MetropolisFullyConnectedMC
 */  

package Backbone.Algo;

import Backbone.System.*;
import Backbone.Util.*;
import JISim.MCSimulation;
import Triggers.*;
import java.util.ArrayList;
import java.util.Random;

 /**  
 *      Ising monte carlo simulation interface to be implemented 
 *   depending on which algorithm one decides to use.
 *  <br>
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2013-01    
 */


public final class MetropolisFullyConnectedMC  implements IsingMC {
    private ParameterBank param;
    private LatticeMagInt lattice;
    private double temperature;
    private double jInteraction;
    private double hField;
    private boolean useDilution;
    private boolean useHeter;
    private Trigger trigger;
    private int tTrigger =0;
    private int tAccumulate =1000; // begin accumulator
    private boolean suppressFieldMessages=false;
    private int N;
    private int Nfixed;
    private int Nfixed0;
    private int s0;
    private int sFix0;
    private double hField0;
    private boolean balancedConf;
    private int run;
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
    public double Nup;         // Number of Up spins
    public double Ndown;	   // Number of down spins	 
    public double Pup;                 // probability of choosing an up spin.
    public double chi;                         // defined = (1/K*T)*[<M^2> - <M>^2]
    public double k = 1;                       // (1.381)*Math.pow(10, -23);
    //public double x = Math.sqrt(4*beta);	   // Input parameter for arcCosh. 
    //public double arcCoshsq4beta = Math.log(x+Math.sqrt(x*x-1));
    //public double H_s = -((1/beta)*arcCoshsq4beta - 4*Math.sqrt((4*beta-1)/(4*beta)));       // spinodal field; the minus changes to the sign to
    public MetropolisFullyConnectedMC(MCSimulation mc, boolean out){
        param = mc.getParams();
        if(param==null){param = new ParameterBank(mc.getParamPostFix());}
        SimPost = mc.getSimPostFix();
        ParamPost = mc.getParamPostFix();        
        hField = param.hField;
        hField0 = hField;
        temperature = param.temperature;
        useDilution = param.useDilution;
        useHeter = param.useHeter;
        N = param.L;param.setN(param.L);
        s0 = param.s;
        jInteraction = param.jInteraction;
        
        jInteraction = 4.0/((double)N);
        measurer = new MeasureIsingSystem(temperature,N);
        triggerReady=false;
    
        lattice = new SingleValueLattice(N,param.s,mc.instId);
        magnetization = lattice.getMagnetization();
        magStaggered = lattice.getMagStaggered();
        measurer.setN(N, lattice.getNFixed());
        // initialize accumulators and energy
        resetData();
        energy=calcIsingEnergy();
        SimProcessParser parser = new SimProcessParser(null,SimPost,ParamPost);
        tFlip = parser.timeToFlipField();
        
        systemOut("Initial M: "+lattice.getMagnetization()+"    E:"+(energy/((double)N)));
        
        if(triggerOn){setTrigger();}
        currentTime=0;
    }

    public MetropolisFullyConnectedMC(MCSimulation mc, boolean out,int nfix,int num,int sfix,boolean bal){
        this(mc,out);
        Nfixed = nfix;
        Nfixed0 = nfix;
        sFix0 = sfix;
        balancedConf = bal;
        setFixed(Nfixed0, N, sFix0, balancedConf);
        measurer.setN(N, lattice.getNFixed());
    }
    
    private void initialize(){
        N = ((SingleValueLattice)lattice).getN();
        Nup = ((SingleValueLattice)lattice).getUpN();
        Ndown = ((SingleValueLattice)lattice).getDownN();	 
        Pup = ((double)Nup)/((double)((SingleValueLattice)lattice).getNDynamic());
    }
    
    /**
    *        doOneStep should do a single monte carlo step.
    */ 
    public void doOneStep(){
        currentTime++;
        
        for (int i =0;i<N;i++){                  
            double rand = Ran.nextDouble();
            // current state of system
            N = ((SingleValueLattice)lattice).getN();
            Nup = ((SingleValueLattice)lattice).getUpN();         // Number of Up spins
            Pup = ((double)Nup)/((double)((SingleValueLattice)lattice).getNDynamic());
            int currMag = lattice.getMagnetization();
            // figure out flip direction
            boolean flipDown = false;
            if(rand < Pup){flipDown = true;}
            int newSpin = (flipDown) ? (-1) : 1;
            
            // Metropolis
            int newMag = lattice.getMagnetization()+2*newSpin;
            double acceptance;
            double dE = calcIsingEnergy(newMag)- calcIsingEnergy(currMag);
            acceptance = Math.exp(-1.0*dE/temperature);
            double toss = Ran.nextDouble();
                if( (dE<=0) || (acceptance>toss) ){
                    ((SingleValueLattice)lattice).setValue(newSpin,currentTime);
                    energy += dE;
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
        measurer.updateM(magnetization, true);
        mcs++;
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
        if(param.useHeter){
            setFixed(Nfixed0, N, sFix0, balancedConf);
        }
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
            if(param.useHeter){setFixed(Nfixed0, N, sFix0, balancedConf);}
        }else{
            lattice.setInitialConfig(t,run,"");
        }
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
                System.out.println("MetropolisFullyConnectedMC |  Field Flipped at t:"+currentTime+" now hField:"+hField);
            }
        }
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
    *         changetT should set the temperature to the given value.
    * 
    *  @param temp - temperature to set simulation to
    */ 
    @Override
    public void changeT(double temp){
        measurer.changeTemp(temp);
        temperature = temp;
    }
    /**
    *         changetH should set the magnetic field to the given value.
    * 
    *  @param hnow - new magnetic field
    */ 
    @Override
    public void changeH(double hnow){
        hField = hnow;
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
        measurer.changeTemp(temp);
        hField = hnow;
        temperature = temp;
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
        double mag = (double)magnetization;
        double num = (double) N;
        e = -((0.5*jInteraction*mag*mag)+(hField*mag));
        return e;
    }

    /**
    *       calcEnergy calculates the energy of the lattice in an Ising Model 
    * 
    * @return energy in lattice for Ising model
    */
    private double calcIsingEnergy(int mag){
        double e=0.0;
        // Simple version 
        double m = (double)mag;
        double num = (double) N;
        N = lattice.getN();
        e = -((0.5*jInteraction*m*m)+(hField*m));
        // ammar
        //e =  (jInteraction)*(0.5*(N - magnetization*magnetization)) - (magnetization*hField);
        // Kang
        //e = -jInteraction*(mag)+jInteraction+hField;
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
    *         setSeed should set the random number seed to the value given
    * 
    *  @param seed1 - random number seed
    */ 
    public void setSeed(int seed1){
        currentSeed = seed1;
    }

    @Override
    public int getSeed() {
        return currentSeed;
    }

    @Override
    public double getSpecificHeat() {        
        return(measurer.getSpecificHeat());
    }

    @Override
    public void setTriggerOnOff(boolean tr) {
        triggerOn = tr;
    }
    
    public void setFixed(int fix, int num, int sfix,boolean balanced){
        ((SingleValueLattice)lattice).setFixed(fix, num, sfix, balanced);
        adjustField(((SingleValueLattice)lattice).getNFixedUp(),((SingleValueLattice)lattice).getNFixedDown());
    }
    
    private void adjustField(int nfixup,int nfixdown){
        double hEff = 4.0*((double)(nfixup-nfixdown))/((double)N);
        systemOut("Adding effective Field : "+hEff);
        hField = hField0+hEff;
    }
    /**
    *   getJinteraction should return the current strength of interaction.
    * 
    */
    @Override
    public double getJinteraction(){
        return jInteraction;
    }
    
    
    private void systemOut(String out){
        System.out.println("MetropolisFullyConnectedMC | "+out);
    }
}

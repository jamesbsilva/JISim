package Backbone.Algo;

/**
 * 
 *   @(#) InvasionPercLat
 */  

import AnalysisAndVideoBackend.DisplayMakerCV;
import Backbone.System.AtomicBondsSiteLattice;
import Backbone.System.SimSystem;
import Backbone.Util.DataSaver;
import Backbone.Util.DirAndFileStructure;
import Backbone.Util.ParameterBank;
import JISim.MCSimulation;
import Triggers.*;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;

/**  
 *   Basic Implementation of InvasionPercLat Algorithm Monte Carlo Simulation.
 *  Single thread implementation of simulation
 *  <br>
 *  @param mc - MCSimulation - Monte Carlo simulation class hosting 
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2013-04    
 */
public class InvasionPercLat implements MCAlgo{
    private class bondVal implements Comparable<bondVal>{
        private int bond;
        private double value;
        public bondVal(int b, double v){
            bond = b;
            value = v;
        }
        public double getValue(){return value;}
        public int getInd(){return bond;}
        public void setValue(double v){value =  v;}
        @Override
        public int compareTo(bondVal t) {
            double old = t.getValue();
            if( old == this.value ){
                return 0;
            }else if( old > this.value ){
                return 1;
            }else{
                return -1;
            }
        }
    }
    private ParameterBank param;
    private AtomicBondsSiteLattice lattice;
    private boolean useLongRange;
    private boolean useDilution;
    private boolean useHeter;
    private boolean crossEdge = false;private boolean savedCrossEdge = false;
    private Trigger trigger;
    private int tTrigger = 0;
    private int tAccumulate = 1000; // begin accumulator
    private boolean suppressFieldMessages=false;
    private ArrayList<bondVal> fringeFront;
    private HashMap<Integer,Boolean> fringe;
    private HashMap<Integer,Double> fringeVal;
    private int fringeSize = 1000;
    private int fringeAll = 0;
    private double maxMin = 0;
    private int coordination;
    private double fringeBound = 2.0;
    private int initLocation;
    
    private int Q;
    private int L;private int N;
    private int D;private int R;
    private int minX;private int maxX;
    private int minY;private int maxY;
    private int run;
    private boolean output = true;
    private Random Ran = new Random();
    private String SimPost;
    private String ParamPost;  
    private int percSize = 0;  
    public int currentSeed;
    public int mcs = 0; // number of MC moves per spin
    public int tFlip = 0;
    public int currentTime;
    
    /**  @param mc - MCSimulation - Monte Carlo simulation class hosting 
    *     @param out - boolean to determine if outputting text confirmations into console
    */
    public InvasionPercLat(MCSimulation mc, boolean out){
        param = new ParameterBank(mc.getParamPostFix());
        SimPost = mc.getSimPostFix();
        ParamPost = mc.getParamPostFix();
        useLongRange = param.useLongRange; useDilution = param.useDilution; useHeter = param.useHeter;
        L = param.L; N = param.N; D = param.D; R = param.R;
        initLocation = (int)((double)(N+L)/2.0);
        
        // determine if outputting
        output = out;
        if(output){System.out.print("Run:"+mc.getRuns()+"     ");}
        run = mc.getRuns();  
        
        if( mc.getSimSystem().getClass().getName().contains("AtomicBondsSiteLattice") ){
            lattice = (AtomicBondsSiteLattice)mc.getSimSystem();
        }else{
            System.err.println("THIS CLASS (InvasionPercLat NEEDS A ATOMICBONDSSITELATTICE AS SYSTEM.");
        }
        
        percSize = lattice.getMagnetization();
    
        fringeFront = new ArrayList<bondVal>();
        fringe = new HashMap<Integer,Boolean>();
        fringeVal = new HashMap<Integer,Double>();
        coordination = lattice.getBonds(0, 0, 0).length;
        for(int u = 0; u < (coordination*N);u++){
            fringe.put(u, false);
            fringeVal.put(u, 0.0);
        }
        
        // initialize accumulators and energy
        resetData();    
        currentTime=0;
    }
    
    /**  
    *     @param out - boolean to determine if outputting text confirmations into console
    */
    public InvasionPercLat(String par, String sim,int runIn, boolean out){
        param = new ParameterBank(par);
        SimPost = sim; ParamPost = par;
        useLongRange = param.useLongRange; useDilution = param.useDilution; useHeter = param.useHeter;
        L = param.L; N = param.N; D = param.D; R = param.R;
        minX = (int)((double)(L)/2.0);
        maxX = minX;maxY = minX; minY = minX;
        initLocation = (int)((double)(N+L)/2.0);
        
        // determine if outputting
        output = out;
        if(output){System.out.print("Run:"+runIn+"     ");}
        run = runIn;  
        
        lattice = new AtomicBondsSiteLattice(0,par);
        
        percSize = lattice.getMagnetization();
    
        fringeFront = new ArrayList<bondVal>();
        fringe = new HashMap<Integer,Boolean>();
        fringeVal = new HashMap<Integer,Double>();
        coordination = lattice.getBonds(0, 0, 0).length;
        for(int u = 0; u < (coordination*N);u++){
            fringe.put(u, false);
            fringeVal.put(u, 0.0);
        }
        
        // initialize accumulators and energy
        resetData();    
        currentTime=0;
    }
    
    
    /**+
    *        doOneStep should do a single monte carlo step.
    */ 
    public void doOneStep(){
        currentTime++;
    
        if(currentTime  == 1){
            System.out.println("Init In | "+initLocation+"    Total Sites | "+N);
            addToFringe( getX(initLocation),getY(initLocation),getZ(initLocation) );
            lattice.setValue(initLocation, 1, currentTime);       
        }
        
        int minBondInd = getMin();
        growToSite(minBondInd);
        addToFringe( getX(getSiteOutFromBond(minBondInd)), 
                getY(getSiteOutFromBond(minBondInd)), 
                getZ(getSiteOutFromBond(minBondInd)) );
        percSize = lattice.getMagnetization();
        //System.out.println("Bond | "+minBondInd+"     Index | "+(int)Math.floor((double)minBondInd/(double)coordination));
    }

    private boolean doneCrossingEdge(){
        return crossEdge;
    }
    
    private void refillFringeFront(){
        int ind;
        double val;
        for(int u = 0; u < fringe.size();u++){
            if(fringe.get(u)){
                if( fringeFront.size() < fringeSize ){
                    fringeFront.add(new bondVal(u,fringeVal.get(u)));                      
                    Collections.sort(fringeFront);
                    fringeBound = fringeFront.get(0).getValue();
                }else{
                    if(fringeVal.get(u) < fringeBound){
                        fringeFront.set(0, new bondVal(u,fringeVal.get(u)));
                        Collections.sort(fringeFront);
                    }
                }
            }
        }
    }
    
    
    
    private void growToSite(int bond){
        int oldVal = lattice.getValue(getSiteOutFromBond(bond));
        //System.out.println("Grow to | "+getSiteOutFromBond(bond)+"     Bond | "+bond
        //        +"     t | "+currentTime+"    clusterSize | "+lattice.getMagnetization()
        //        +"    fringeFront | "+fringeFront.size()+
        //        "    old | "+lattice.getValue(getSiteOutFromBond(bond)));
        if( oldVal == 1 )printFringe();
        if( oldVal == 1 )System.out.println("------------------------------------------");
        int site = getSiteOutFromBond(bond);
        lattice.setValue(site, 1, currentTime);
        if(getX(site) > maxX ){maxX = getX(site);}if(getX(site) < minX ){minX = getX(site);}
        if(getY(site) > maxY ){maxY = getY(site);}if(getY(site) < minY ){minY = getY(site);}
        crossedEdge();
        removeMin();
        if( oldVal == 1 )printFringe();    
    }
    
    private void crossedEdge(){
        if(savedCrossEdge){
        }else{
            savedCrossEdge = false;
            if(minX <= 0)savedCrossEdge = true;
            if(maxX >= (L-1))savedCrossEdge = true;
            if(minY <= 0)savedCrossEdge = true;
            if(maxY >= (L-1))savedCrossEdge = true;
            if(savedCrossEdge)saveClusterData();
        }
    }
    
    public void saveClusterData(){
        crossEdge = true;
        String nameOut = "InvasionClust2D";
        String fname = nameOut+"cluster-L"+L+"-size-"+percSize+".dat";
        DirAndFileStructure dir = new DirAndFileStructure();
        fname = dir.getDataDirectory(nameOut)+fname;
        System.out.println("Saving cluster | "+fname);
        //read file first    
        /*
        try{
            PrintStream out = new PrintStream(new FileOutputStream(
                fname,false));
                ArrayList<Integer> temp = clhandler.getBufferIntAsArrayList(growKernel, 0, N, false);
                for(int u = 0; u < temp.size();u++){
                    out.println(temp.get(u));
                }
                out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }*/
        fname= "percSize"+"cluster-L"+L+".dat";
        fname = dir.getDataDirectory(nameOut)+fname;
        //read file first    
        try{
            PrintStream out = new PrintStream(new FileOutputStream(
                fname,true));
                out.println(percSize);
                out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
    

    private void printFringe(){
        System.out.println("++++++++++++++++++++++++++++++++++++++" );
        for(int u = 0; u < fringeFront.size();u++){
            System.out.println("Index | "+u+
                    "    siteInd | "+getSiteOutFromBond(fringeFront.get(u).getInd())+
                    "   Bond | "+fringeFront.get(u).getInd()+"     Val | "+fringeFront.get(u).getValue() );
        }
    }
    
    private void removeMin(){
        int newLoc = getSiteOutFromBond(fringeFront.get(fringeFront.size()-1).getInd()); 
        fringe.put(fringeFront.get(fringeFront.size()-1).getInd(), false);
        fringeAll--;
        fringeVal.put(fringeFront.get(fringeFront.size()-1).getInd(), 0.0);
        fringeFront.remove(fringeFront.size()-1);
        boolean searching = true;
        int u = 0;
        while(searching){
            if( getSiteOutFromBond(fringeFront.get(u).getInd()) == newLoc ){
                fringeVal.put(fringeFront.get(u).getInd(), 0.0);
                fringe.put(fringeFront.get(u).getInd(), false);
                fringeFront.remove(u);
                fringeAll--;
                if(fringeFront.size() == 0 ){
                    refillFringeFront();
                }
            }else{
                u++;
            }
            if(u >= (fringeFront.size())){searching = false;}
        }
    }
    
    private int getMin(){
        //System.out.println("Fringe Size | "+fringeAll+"  percSize | "+percSize+"  Bond Value | "+fringeFront.get(fringeFront.size()-1).getValue()+"    diffToNext | "
        //        +(fringeFront.get(fringeFront.size()-2).getValue()-fringeFront.get(fringeFront.size()-1).getValue())
        //        +"  maxMin | "+maxMin);
        if(maxMin < fringeFront.get(fringeFront.size()-1).getValue()){
            maxMin = fringeFront.get(fringeFront.size()-1).getValue();
        }
        return fringeFront.get(fringeFront.size()-1).getInd();
    }
    
    
    private void addToFringe(int i, int j, int k){
        double[] siteBonds = lattice.getBonds(i, j, k);
        double val;
        int ind; int indOut;
        for(int u = 0; u < siteBonds.length;u++){
            val = siteBonds[u];
            ind = siteBonds.length*(i+j*L+k*L*L)+u;
            indOut = getSiteOutFromBond(ind);
            if( lattice.getValue(indOut) != 1 && !lattice.isThisFixed(indOut) ){
                if(!fringe.get(ind)){
                    fringe.put(ind, true);
                    fringeAll++;
                    fringeVal.put(ind, val);
                    // determine if putting in fringe front
                    if( fringeFront.size() < fringeSize ){
                        fringeFront.add(new bondVal(ind,val));                      
                        Collections.sort(fringeFront);
                        fringeBound = fringeFront.get(0).getValue();
                    }else{
                        if(val < fringeBound){
                            fringeFront.set(0, new bondVal(ind,val));
                            Collections.sort(fringeFront);
                        }
                    }
                }
            }
        }
    }
    
    private int getX(int ind){
        return (ind%L);
    }
    private int getY(int ind){
        if(D < 2){return 0;}
        return ((int)((double)ind/(double)L)%L);
    }
    private int getZ(int ind){
        if(D < 3){return 0;}
        return ((int)((double)ind/(double)(L*L))%L);
    }
    
    private int getSiteInFromBond(int bind){
        return (int) Math.floor((double)bind/(double)coordination);
    }
    
    private int getSiteOutFromBond(int bind){
        int site = (int) Math.floor((double)bind/(double)coordination);
        int ind = bind%coordination; 
        int i = getX(site);int j = getY(site);
        int k = 0;
        if(D == 3){
             k = getZ(site);
        }
        switch( ind ){
            case 0: ind = ((i+1)%L)+j*L+k*(L*L);
                     break;
            case 1: ind = ((i+L-1)%L)+j*L+k*(L*L);
                     break;
            case 2: ind = i+((j+1)%L)*L+k*(L*L);
                     break;
            case 3: ind = i+((j+L-1)%L)*L+k*(L*L);
                     break;
            case 4: ind = i+j*L+((k+1)%L)*(L*L);
                     break;
            case 5: ind = i+j*L+((k+L-1)%L)*(L*L);
                     break;
            default:ind = site;
                     break;
        }
        return ind;
    }
    
    /**
    *         setTrigger should initialize the trigger to be used in the simulation
    */ 
    public void setTrigger(){
        trigger = param.setProperTrigger(param, trigger, SimPost, output);
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
        DataSaver save = new DataSaver(lattice.getInstId());
        //assert field is right
        for(int t = 0; t<tFinal;t++){
            doOneStep();
            if(t>tInitial){save.saveConfig(lattice,run,t);}
        }
    }
    
    /**
    *   getSimSystem should return the lattice object
    * 
    */
    @Override
    public SimSystem getSimSystem(){return (SimSystem)lattice;}
    
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
        percSize = lattice.getMagnetization();
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
        percSize = lattice.getMagnetization();
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
        trigger.reset2();//reset but keep std or go into intervention mode
        if(t==0){lattice.initialize(param.s);
        if(param.useHeter){lattice.setFixedLatticeValues();}
        }else{
        lattice.setInitialConfig(t,run,"");}
        percSize = lattice.getMagnetization();
        setSeed(seed);
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
    *   getSeed should return the current random number seed
    * 
    */
    @Override
    public int getSeed(){return currentSeed;}

    // test the class
    public static void main(String[] args) {
        DisplayMakerCV dispConf = new DisplayMakerCV( "Configuration" );
        for(int u = 0; u < 2500;u++){
            InvasionPercLat inv = new InvasionPercLat("","",0,true);
            while(!inv.doneCrossingEdge()){
                inv.doOneStep();
                //dispConf.addDisplayFrame(inv.lattice.getSystemImg());
            }
            inv.doOneStep();
        }
        System.out.println("Done!!");
    }
}


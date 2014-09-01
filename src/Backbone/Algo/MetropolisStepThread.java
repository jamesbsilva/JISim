package Backbone.Algo;

/**
 * 
 *      @(#)  MetropolisStepThread
 */  
 
import Backbone.System.LatticeMagInt;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;

/** 
 *   Thread for a metropolis monte carlo step. Work is divided by amount of
 *  cores ie for N size lattice and C cores then each step thread gets N/Cmax
 *  where Cmax is the closest number to C that divides N.
 *  <br>
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2012-01    
 */


//------------------------------------------------------------ 
 //------------------------------------------------------------          
 // Information needs to come back Energy occurs
 // Coming In-  W,lattice,              
  public class MetropolisStepThread implements Callable {
    private LatticeMagInt lattice;
    private double jInteraction;
    private double hField;
    private ConcurrentHashMap<Integer,Double> w;
    private boolean useDilution;
    private boolean useHeter;
    private int L;
    private int range;
    private int N;
    private int D;
    private double energy;
    private int currentTime;
    private Random Ran = new Random();
    private int loops;
    private int tId;
    private int xoff=0;
    private int yoff=0;
    private int divisions; 
    

 /** 
 *  @param t - currentTime for mc step
 *  @param lat - the lattice for step
 *  @param w0 - w matrix for mc steps 
 *  @param loop - Chunk assigned to thread
 *  @param r0 - Range for interaction
 *  @param jit - interaction constant
 *  @param hfiel - magnetic field
 *  @param useHet - boolean for if this is a heter simulation
 *  @param useDil - boolean for if this is a simulation as dilution
 *  @param see - seed for step
 *  @param threadNumber - thread id
 */ 
    public MetropolisStepThread(int t, LatticeMagInt lat,ConcurrentHashMap w0, 
           int loop,int r0,double jit,
            double hfiel, boolean useHet,boolean useDil, int see,int threadNumber) {
      this.lattice = lat;
      this.currentTime = t;
      this.jInteraction = jit;
      this.hField = hfiel;
      this.w = w0;
      this.range = r0;
      this.loops = loop;
      this.useDilution = useDil;
      this.useHeter = useHet;
      this.L = lat.getLength();
      this.N = lat.getN();
      this.D = lat.getDimension();
      Ran.setSeed(see);
      this.tId = threadNumber;
    }
    
    /**
     *      call performs the main code of this callable metropolis step.
     * 
     * @return - energy change
     */
    @Override
    public Double call() {
        Double updateLat;
        doOneStepPiece();
        //output step movement
        updateLat = energy;
        return updateLat;
    }
    
    /**
     *      doOneStepPiece performs a single random site visit of the 
     *  metropolis algo.
     */
    public  void doOneStepPiece(){

        //if(currentTime%100==0){System.out.print("\r time:"+currentTime+"       ");}
        energy=0;
        divisions =(int)(N/loops);

        //if not a perfect square just offset in i
        if((Math.sqrt(divisions)-((int)(Math.sqrt(divisions))) != 0)){
        xoff = (int) (L*tId/(divisions+1));
        }else{ 
            divisions = (int)Math.sqrt(divisions);
            xoff= ((int)(L/(divisions)))*(tId%divisions);
            yoff = ((int)(L/(divisions)))*((int)((tId-xoff)/divisions));
        }

        for(int f = 0;f<loops;f++) {
        int i= ((int) (Ran.nextDouble()*L)+xoff)%L ;
        int j = ((int) (Ran.nextDouble()*L)+yoff)%L ;
        int k;
        if(D==3){k=(int) (Ran.nextDouble()*L);}else{k=0;}

        // Check if fixed. If fixed then its been visited but skip
        if (useHeter == true){
            if(lattice.isThisFixed(i,j,k)){continue;}
        }

        int nearestSum;

        if (D==3){
            nearestSum =lattice.getNeighSum(i, j);
        }else{
            // Synchronized
            nearestSum = lattice.getNeighSum(i, j);}     

        // diluted version
        if(useDilution ==true){
            if(D==3){dilutedStep(i,j,k,nearestSum);}else{dilutedStep(i,j,0,nearestSum);}
        }else{
            // regular Step
            if(D==3){regularStep(i,j,k,nearestSum);}else{regularStep(i,j,0,nearestSum);}
        }
        }
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
                acceptance = (double) w.get(
                        3+isPositive(nearestSum)*4+Math.abs(nearestSum)*4*2);
                double dE = -currentSpin*(nearestSum*jInteraction+hField);
                if((dE<=0)||(acceptance>Ran.nextDouble())) {
                    lattice.setValue(i, j,k, 1,currentTime);
                    //---------------add 3D

                    energy += dE;
                    }}
            // change to -1
            else{
                acceptance = (double) w.get(
                        2+isPositive(nearestSum)*4+Math.abs(nearestSum)*4*2);
                double dE = currentSpin*(nearestSum*jInteraction+hField);
                if((dE<=0)||(acceptance>Ran.nextDouble())) {
                    lattice.setValue(i, j,k,-1,currentTime);
                //---------------add 3D

                energy += dE;

                //Removing since setValue is doing now
                //    	            magnetization += -1;
                }
            }
        }
        else if(currentSpin == 1){
            // change to - 1
            if(r == 0){
            acceptance = (double) w.get(
                    isPositive(currentSpin)+isPositive(nearestSum)*4+Math.abs(nearestSum)*4*2);
            double dE = 2*currentSpin*(nearestSum*jInteraction+hField);
                if((dE<=0)||(acceptance>Ran.nextDouble())) {
                lattice.setValue(i,j,k,-1,currentTime);
                //---------------add 3D

                energy += dE;
                //Removing since setValue is doing now
                //magnetization += -2;
            }}
            // change to 0
            else{
            acceptance = (double) w.get(
                    isPositive(currentSpin)+isPositive(nearestSum)*4+Math.abs(nearestSum)*4*2);
            double dE = currentSpin*(nearestSum*jInteraction+hField);
            if((dE<=0)||(acceptance>Ran.nextDouble())) {
                lattice.setValue(i, j,k, 0,currentTime);
            //---------------add 3D

            energy += dE;
            //Removing since setValue is doing now
            // magnetization += -1;
            }
            }  
        }
        // current Spin is -1
        else{
        // change to + 1
        if(r == 0){
            acceptance = (double) w.get(
                    isPositive(currentSpin)+isPositive(nearestSum)*4+Math.abs(nearestSum)*4*2);
            double dE = 2*(nearestSum*currentSpin*jInteraction+currentSpin*hField);
            if((dE<=0)||(acceptance>Ran.nextDouble())) {
                lattice.setValue(i, j,k, 1,currentTime);

                energy += dE;
                //Removing since setValue is doing now
                //        magnetization += 2;
        }}
        // change to 0
        else{
            acceptance = (double) w.get(
                    isPositive(currentSpin)+isPositive(nearestSum)*4+Math.abs(nearestSum)*4*2);
            double dE = (nearestSum*currentSpin*jInteraction+currentSpin*hField);
        if((dE<=0)||(acceptance>Ran.nextDouble())) {
            lattice.setValue(i, j,k, 0,currentTime);
        //---------------add 3D

        energy += dE;
        //Removing since setValue is doing now
        //magnetization += 1;
        }}  
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

        acceptance = (double) w.get(isPositive(currentSpin)+
                isPositive(nearestSum)*4+Math.abs(nearestSum)*4*2);
        double toss = Ran.nextDouble();

        if((dE<=0)||(acceptance>toss)){
        int newSpin = -currentSpin;
        if(useHeter){
        if(lattice.isThisFixed(i,j,k)){System.out.println("Changing Fixed");}}
        //System.out.print("Mi="+lattice.getMagnetization());
        lattice.setValue(i,j,k,newSpin,currentTime);
        //System.out.println("      Mf="+lattice.getMagnetization());
        energy += dE;
        }
    }


    
    /**
    *        isPostitive returns true if number is positive
    * 
    *   @param number - input number
    */ 
    public static int isPositive(int number){
        if (number >= 0){return 1;}else{return 0;}
    }

}	      

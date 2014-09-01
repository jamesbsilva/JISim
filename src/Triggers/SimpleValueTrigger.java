package Triggers;

import Backbone.Util.SimProcessParser;

/**
* 
*    @(#) SimpleValueTrigger 
*/  
/**  
*      Trigger based on the value of the magnetization. If the magnetization is 
*  larger than a given percentage of the stable value.
*  <br>
* 
*  @param s - Initialized value of lattice
*  @param thresh - What percent of stable state to trigger off
*  @param N - Size of lattice
* 
* @author      James B Silva <jbsilva @ bu.edu>                 
* @since       2011-11   
*/
public final class SimpleValueTrigger implements Trigger{
    private int stable;
    private double valThreshold=0.3;
    private int t=0;
    private int tFlip=0;
    private int tTrigger =0;
    private int tOffset =0;

    public SimpleValueTrigger(int stab){
        stable = -stab;
    }
	
    /**  @param s - Initialized value of lattice
    *  @param thresh - What percent of stable state to trigger off
    *  @param N - Size of lattice
    *  @param paramPost - postfix for parameter file
    */ 
    public SimpleValueTrigger(int s,double thresh,int N, String paramPost){
        this( s, thresh, N, paramPost, paramPost);
    }

    // MCSimulation Going though here

    /**  @param s - Initialized value of lattice
    *  @param thresh - What percent of stable state to trigger off
    *  @param N - Size of lattice
    *  @param paramPost - postfix for parameter file
    *  @param fname - postfix for simulation process file
    */ 
    public SimpleValueTrigger(int s,double thresh,int N, String fname, String paramPost){
        stable = -s; 

        if(thresh> 1){System.out.println("Supposed to input as a percentage");thresh = thresh/100;}
        valThreshold = thresh;
        stable  = stable*N;

        SimProcessParser sim = new SimProcessParser(null,fname,paramPost);
        tFlip = sim.timeToFlipField();

        printTriggerSettings();		
    }

    @Override
    public int getTriggerTime(){return (tTrigger-tOffset);}

    /**
    *    printTriggerSettings outputs the settings for this instance of a deviation trigger.
    */
    public void printTriggerSettings(){
        System.out.println("---------------------------------------------");
        System.out.println("New run starting.  Stable state :"+(-stable));
        System.out.println("Trigger if m  is "+(valThreshold*100)+"  percent of m stable after t ="+tFlip);
        System.out.println("---------------------------------------------");
    }

    /**
    * 	triggerNow determines if trigger is ready to go off.
    * 
    * @param m - Magnetization of lattice
    * 
    */
    public boolean triggerNow(double m,double e){
        return triggerNow(m);
    }
        
    @Override
    public boolean triggerNow(){return false;}

    /**
    * 	triggerNow determines if trigger is ready to go off.
    * 
    * @param m - Magnetization of lattice
    * 
    */
    @Override
    public boolean triggerNow(double m){
        t++;
        boolean trigger=false;

        // Assert trigger working right
        if(t>tFlip){
        if(Math.abs(m)>Math.abs(valThreshold*stable) && Math.signum(m)==Math.signum(stable)){
                System.out.println("Triggered at t= "+t+"     M: "+m);trigger=true;}
        }

        if(trigger){tTrigger = t;}

        return trigger;
    }

    /**
    *    reset() resets the trigger including the current standard deviation value.
    */
    @Override
    public void reset(){
        t=0;
    }

    /**
    *    reset() resets the trigger including the current standard deviation value.
    *    
    *    @param s - new stable state 
    */
    public void reset(int s){
    }

    @Override
    public void resetTime(){t=0;}
    @Override
    public int getT(){return t;}
    @Override
    public void setT(int s){}
    @Override
    public void reset2(){};

    @Override
    public int getTcutoff() {
        throw new UnsupportedOperationException("Not supported yet.");
    }  
        
    // test the class
    public static void main(String[] args) {
    }
}

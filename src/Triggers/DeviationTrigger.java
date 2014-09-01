package Triggers;

import Backbone.Util.ParameterBank;
import Backbone.Util.SimProcessParser;

/**
* 
*  @(#) DeviationTrigger 
*/  
/** 		Trigger based on change in magnetization and how large it is 
*  relative to the standard deviation. Each the change is bigger than thresh 
*  (the variable) times the standard deviation the trigger counts a hit.
*  Once 5 hits are reached the trigger goes off. The amount of hits 
*  before trigger goes off can be changed.
*  <br>
* 
*  @param s - Initialized value of lattice
*  @param thresh - How many standard deviation before trigger goes off as a hit
*  @param tcut- Time after which stop calculating the standard deviation and begin looking for trigger to go off
*  
* @author      James B Silva <jbsilva @ bu.edu>                 
* @since       2011-10   
*/
public final class DeviationTrigger implements Trigger{
    private double stdThreshold=5;
    private double std;
    private double mean;
    private double mSquaredSum=0;
    private double mSum=0;
    private int size;
    private int stable;
    private double devThreshold=0.3;
    private int tCutoff; 
    private double mLast;
    private double mdel;
    private int t=0;
    private int tFlip =500;
    private int hitThreshold =0;
    private int hits = 5;
    private int tTrigger =0;
    private int tOffset =100;

    // Testing variables
    /*private boolean saveStd=false;
    private int tDebug = 2300;
    private LinkedList<Double> delData;
    private LinkedList<Double> stdData;*/

    public DeviationTrigger(int stab){
        std =0;size=0;stable = -stab;
    }
    
    public DeviationTrigger(int s,double stab){
        std =0;size=0; stable = -s; stdThreshold = stab;
    }


    // MCSimulation Going though here
    /**  @param s - Initialized value of lattice
    *  @param thresh - How many standard deviation before trigger goes off as a hit
    *  @param tcut - Time after which stop calculating the standard deviation and begin looking for trigger to go off
    */  
    public DeviationTrigger(int s,double thresh, double tcut, String fname){
        std =0;size=0;stable = -s; stdThreshold = thresh;tCutoff = (int) tcut;
        printTriggerSettings();
        ParameterBank param = new ParameterBank(fname);
        devThreshold = devThreshold*param.N;
        param=null;
        SimProcessParser sim = new SimProcessParser(fname);
        tFlip = sim.timeToFlipField();
        System.out.println("time flip:"+tFlip);
    }
    public DeviationTrigger(double standard, int stab){
        std= standard;
        stable = -stab;
    }

    public DeviationTrigger(double standard,double standThresh, int unstab,String fname){
        std= standard;
        stable = -unstab;
        stdThreshold=standThresh;
        SimProcessParser sim = new SimProcessParser(fname);
        tFlip = sim.timeToFlipField();
        System.out.println("time flip:"+tFlip);

        ParameterBank param = new ParameterBank();
        devThreshold = devThreshold*param.N;
    }

    @Override
    public int getTriggerTime(){return (tTrigger-tOffset);}
    @Override
    public int getT(){return t;}
    public double getStd(){return std;}
    public double getMdel(){return mdel;}
    
    /**
    * 
    * 	setHits sets the amount of hits before trigger goes off.
    * 
    * @param k - Amount of hits before trigger goes off
    */
    public void setHits(int k){hits=k;}
    
    public void setStd(double m){
        size++;
        mSum+=m;
        mSquaredSum+=m*m;

        mean = mSum/size;

        std = Math.sqrt(Math.abs(mSquaredSum-size*mean*mean)/(size-1));
    }
    /**
    *    printTriggerSettings outputs the settings for this instance of a deviation trigger.
    */
    public void printTriggerSettings(){
        System.out.println("Starting with stable state :"+(-stable));
        System.out.println("Trigger if m del is "+stdThreshold+"  times sigma m after "+tCutoff+ "  seconds.");
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

        if(t>tCutoff){
            if(t==(tCutoff+3)){System.out.println("Std:"+std+ "   std threshold:"+(std*stdThreshold)+ "    mean:"+mean+"       trig:"+devThreshold);}
            mdel = m-mLast;

            if ((stable == 1)  && (mdel>0) && (Math.abs(mdel)> stdThreshold*std)){
                hitThreshold++;if(hitThreshold>hits){trigger = true;}
            }else if((stable == -1)  && (mdel<0) && (Math.abs(mdel)> stdThreshold*std)){
                hitThreshold++;if(hitThreshold>hits){trigger = true;}
            }else if((stable == -2)  && (Math.abs(mdel)> stdThreshold*std)){
                hitThreshold++;if(hitThreshold>hits){trigger = true;}
            }else{
                trigger = false;
            }

            if(Math.abs(m-mean)> devThreshold){
                hitThreshold++;
                if(hitThreshold>hits){trigger = true;}
            }
        }else{
            if(t>tFlip){setStd(m);trigger=false;}
        }

        // Assert trigger working right
        if(Math.abs(stable)==1){
        if(Math.abs(m)>Math.abs((0.99*stable)) && (m/Math.abs(m))==stable){System.out.println("Trigger Problem");}}

        mLast= m;

        if(trigger){tTrigger = t;}

        return trigger;
    }


    /**
    *    reset() resets the trigger including the current standard deviation value.
    */
    public void reset(){
        t=0;std=0;size=0;mSum=0;mSquaredSum=0;mLast=stable;hitThreshold =0;mean = 0;
    }

    /**
    *    reset2() resets the trigger but preserves the current standard deviation value.
    */
    @Override
    public void reset2(){
        size=0;mSum=0;mSquaredSum=0;hitThreshold =0;
    }

    /**
    *    reset() resets the trigger including the current standard deviation value.
    *    
    *    @param s - new stable state 
    */
    public void reset(int s){
        stable =s;
        t=0;std=0;size=0;mSum=0;mSquaredSum=0;mLast=stable;hitThreshold =0;
    }

    @Override
    public void resetTime(){t=0;}

    @Override
    public int getTcutoff() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean triggerNow() {
        return false;
    }

    @Override
    public void setT(int t) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    // test the class
    public static void main(String[] args) {
        DeviationTrigger trig = new DeviationTrigger(-1,2.0,8.0,"");

        System.out.print(trig.triggerNow(-2));System.out.println("std :"+trig.std);
        System.out.print(trig.triggerNow(-4));System.out.println("std :"+trig.std);
        System.out.print(trig.triggerNow(-4));System.out.println("std :"+trig.std);
        System.out.print(trig.triggerNow(-4));System.out.println("std :"+trig.std);
        System.out.print(trig.triggerNow(-5));System.out.println("std :"+trig.std);
        System.out.print(trig.triggerNow(-5));System.out.println("std :"+trig.std);
        System.out.print(trig.triggerNow(-7));System.out.println("std :"+trig.std);
        System.out.print(trig.triggerNow(-9));System.out.println("std :"+trig.std);
        System.out.println("Standard Dev:"+trig.std);
        System.out.print(trig.triggerNow(5));System.out.println("std :"+trig.std);
        System.out.print(trig.triggerNow(10));System.out.println("std :"+trig.std);
        System.out.println("______________________________________________________________");

        trig = new DeviationTrigger(1,2.0,8.0,"");

        System.out.println(trig.triggerNow(2));
        System.out.println(trig.triggerNow(4.0));
        System.out.println(trig.triggerNow(4));
        System.out.println(trig.triggerNow(4.0));
        System.out.println(trig.triggerNow(5));
        System.out.println(trig.triggerNow(5));
        System.out.println(trig.triggerNow(7.0));
        System.out.println(trig.triggerNow(9));
        System.out.println(trig.triggerNow(-8));
        System.out.println(trig.triggerNow(-12));

        System.out.println("Done!");
    }
}

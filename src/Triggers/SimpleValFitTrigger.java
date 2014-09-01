package Triggers;

/**
 * 
 *    @(#) SimpleValFitTrigger 
 */  

import AnalysisAndVideoBackend.LinearRegression;
import Backbone.Util.DirAndFileStructure;
import Backbone.Util.SimProcessParser;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

 /**  
 *      Trigger based on the value of the magnetization and linear regression. If the magnetization is 
 *  larger than a given percentage of the stable value it triggers initially but
 *  a second time is used as a second final trigger (m = -.2). Linear regression
 *  is then used to extrapolate a trigger time. The trigger also has an 
 *  intervention mode where it just checks for m greater than second parameter.
 *  A negative value for the first parameter will put the trigger in calibration
 *  mode which will output m and parameter value for the dev and simple triggers.
 *  <br>
 * 
 *  @param s - Initialized value of lattice
 *  @param thresh - What percent of stable state to trigger off
 *  @param N - Size of lattice
 * 
 * @author      James B Silva <jbsilva @ bu.edu>                 
 * @since       2011-11   
 */

public class SimpleValFitTrigger implements Trigger{
    private int stable;
    private double valThreshold=0.3;
    private double val2Threshold = 0.2;
    private double valReset2Mode = 0.8;
    private boolean Reset2mode=false;
    private boolean calibrationMode=false;
    private boolean tFirstTrigger =false;
    private boolean tSecondTrigger =false;
    private int t=0;
    private int tFlip=0;
    private int tTrigger =0;
    private int tOffset =8; // offset for trigger time
    private int tFirst =0;
    private int tSecond=0;
    private int tTimeout = 200;
    private double mLast=0;
    private List<Double> dataX;
    private List<Double> dataY;
    private double std;
    private double mean;
    private double mSquaredSum=0;
    private double mSum=0;
    private double m0;
    private int size;
    private int samples=200;
    private double meanM=0;
    private int N;
    private boolean output=false;

    public SimpleValFitTrigger(int stab){
        stable = -stab;
        dataX = new ArrayList<Double>();
        dataY = new ArrayList<Double>();
    }

    public SimpleValFitTrigger(int s,double thresh,int num, String fname, 
            double res2){
        this(s,thresh,num,fname,res2,false);
    }
        
    // MCSimulation Going though here
    /**  @param s - Initialized value of lattice
    *  @param thresh - What percent of stable state to trigger off initially
    *  @param res2 - what percent to trigger of in Mode 2
    *  @param fname  - SimProcess Postfix for file
    *  @param num - Size of lattice
    */ 
    public SimpleValFitTrigger(int s,double thresh,int num, String fname, double res2, boolean output){
        stable = -s; 
        dataX = new ArrayList<Double>();
        dataY = new ArrayList<Double>();

        if(thresh> 1){
            if(output){
            System.out.println("Supposed to input as a percentage");}
            thresh = thresh/100;} 
        if(Math.abs(res2) > 1){
            if(output){
            System.out.println("Supposed to input as a percentage");}
            valReset2Mode = res2/100;
        }else{valReset2Mode = res2;}

        if(thresh<0){calibrationMode=true;thresh=-thresh;}

        valThreshold = thresh;
        stable  = stable*num;
        N=num;

        SimProcessParser sim = new SimProcessParser(fname);
        tFlip = sim.timeToFlipField();
    }

    @Override
    public int getTriggerTime(){return (tTrigger);}

    /**
    *    printTriggerSettings outputs the settings for this instance of a deviation trigger.
    */
    public void printTriggerSettings(){
        System.out.println("___________________________________________________");
        System.out.println("Starting with stable state :"+(-stable));
        System.out.println("Trigger initially if m  is "+(valThreshold*100)+" "
                + " percent of m initial or m="+(valThreshold*meanM*stable)+" after t ="+tFlip);
        System.out.println("Trigger again if m  is "+(val2Threshold*stable));
        System.out.println("Then perform linear regression to obtain trigger value."
                + "In intervention mode it is just a simple value trigger");
        System.out.println("___________________________________________________");
    }

        
    public void setStd(double m){
        size++;
        mSum+=m;
        mSquaredSum+=m*m;

        mean = mSum/size;

        std = Math.sqrt(Math.abs(mSquaredSum-size*mean*mean)/(size-1));
    }
        
    @Override
    public boolean triggerNow(){return false;}
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

        if(Reset2mode){
              if(Math.abs(m)>Math.abs(valReset2Mode*stable) 
                      && Math.signum(m)==Math.signum(stable)){
                //System.out.println("Triggered off m: "+m+ "  at t: "+t);
                trigger=true;}
        }else{
            if(t>tFlip){
                if(t<(tFlip+samples)){
                    meanM += m;
                }
                if(t==(tFlip+samples)){
                    meanM = meanM/samples;
                    meanM = Math.abs(meanM/stable);
                    if(output){printTriggerSettings();}
                }   

                // First time trigger for user defined part. This is same magnitude as unstable
                if(Math.abs(m)<Math.abs(valThreshold*meanM*stable) && t>(tFlip+samples)
                        && Math.signum(m)==Math.signum(-stable) && !(tFirstTrigger)){
                    tFirst = t;tFirstTrigger =true;
                    System.out.println("Triggered1 off m: "+m+ "  at t: "+t);
                }
                // calculate standard dev for first 300 time steps
                if(t<(tFlip+325) && (t>(tFlip+25)) ){setStd(m);}

                if(tFirstTrigger){
                    //input the data for regression
                    dataX.add((double)t);
                    dataY.add(Math.abs(m-mLast)/std);
                }

                // Second time trigger always done with value 
                if(Math.abs(m)>Math.abs(val2Threshold*stable) && 
                        Math.signum(m)==Math.signum(stable)){
                    tSecond = t;tSecondTrigger =true;
                    System.out.println("Triggered off m: "+m+ "  at t: "+t);
                    trigger=true;
                }

                mLast = m;
            }
            if(trigger && !Reset2mode){
                // do regression and calculate tTrigger
                LinearRegression reg = new LinearRegression();
                reg.setDataX(dataX);
                reg.setDataY(dataY);
                reg.performRegression();
                if(calibrationMode){saveSettings(m);}
                tTrigger = (int)Math.round(reg.calculateTriggerTime()) - tOffset;
                if(tTrigger<0){
                    tTrigger=tFlip;
                    if(output){System.out.println("UNSTABLE EVOLUTION!!!. Lower the field.");}
                }
            }
        }
        return trigger;
    }

        
    public void saveSettings(double m){
        DirAndFileStructure dir = new DirAndFileStructure();
        String fName = dir.getRootDirectory()+"calibration.txt";

        try {
            PrintStream out = new PrintStream(new FileOutputStream(
                fName,true));
            double del = Math.abs(m-mLast);
            double delNorm = del/std;
            System.out.println("m:"+m+"    del:"+del+ "   delNorm:"
                    +delNorm+ "     delta t:"+(tSecond-tFirst));
            out.close();
        } catch (FileNotFoundException e) {
          e.printStackTrace();
        }        
    }
        

    /**
    *    reset() resets the trigger including the current standard deviation value.
    */
    @Override
    public void reset(){
        t=0;tFirstTrigger=false;tSecondTrigger=false;
        std=0;size=0;mSum=0;mSquaredSum=0;mean = 0;
    }

    /**
    *    reset() resets the trigger including the current standard deviation value.
    *    
    *    @param s - new stable state 
    */
    public void reset(int s){
        t=0;std=0;size=0;mSum=0;mSquaredSum=0;stable = s;mean = 0;
        tFirstTrigger=false;tSecondTrigger=false;
    }

    @Override
    public void resetTime(){t=0;}
    @Override
    public int getT(){return t;}
    @Override
    public void setT(int s){}
        
    /**
    *  reset2() sets the trigger into an intervention mode which just 
    * checks intervention happens very by using the second parameter as 
    * the new trigger 
    *
    */
    @Override
    public void reset2(){
        Reset2mode=true;t=0;std=0;
        //System.out.println("Going into second mode with thresh:"+this.valReset2Mode);
        size=0;mSum=0;mSquaredSum=0;mean = 0;
    }
    
    @Override
    public int getTcutoff() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
        
    // test the class
    public static void main(String[] args) {
    }
}

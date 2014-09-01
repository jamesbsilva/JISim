package Backbone.Util;
/**
 * 
 *      @(#)  SystemAverager
 */  

import java.util.ArrayList;

/**  
 *   Class for managing how measurements of system variables and observables are made to keep average of values.
 *  
 *  <br>
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2013-07    
 */
public class SystemAverager {
    private ArrayList<ArrayList<Double>> Observables;
    private ArrayList<Long> updatesObs;
    private short largestOrder = 7;
    
    public SystemAverager(){
        Observables = new ArrayList<ArrayList<Double>>();
        updatesObs = new ArrayList<Long>();
    }
    
    public int addObservable(){
        ArrayList<Double> temp = new ArrayList<Double>();
        for(int i=0;i <= largestOrder;i++){
            temp.add(i,0.0);
        }
        Observables.add(temp);
        updatesObs.add(0l);
        return Observables.size();
    }
    
    
    private double calcFourthCumulantObs(int u){
        double cum;
        ArrayList<Double> obs = Observables.get(u);
        double numerator = getExpectationObs(u,4)+getExpectationObs(u,1,4)
                -8*(getExpectationObs(u,3)*getExpectationObs(u,1))
                +6*(getExpectationObs(u,2)*getExpectationObs(u,1,2)) ;
        double denominator = 3*(
                getExpectationObs(u,2,2)+getExpectationObs(u,1,4)
                -2*(getExpectationObs(u,1,2)*getExpectationObs(u,2)));
       cum = 1-(numerator/denominator);
       return cum;
    }

    public void updateObs(int u, double val){
        updateObs(u,val,true);
    }
    
    public void updateObs(int u, double val , boolean update){
        if(update){updatesObs.set(u,updatesObs.get(u)+1);}
        ArrayList<Double> obs = Observables.get(u);
        obs.set(0,val);
        for(short i=1; i<=largestOrder;i++){
            obs.set(i, (Math.pow(val,i)+obs.get(i)));
        }
        Observables.set(u, obs);
    }
    
    public void clearObs(int u){
        updatesObs.set(u,0l);
        ArrayList<Double> obs = Observables.get(u);   
        obs.set(0,0.0);
        for(short i=1; i<=largestOrder;i++){
            obs.set(i, 0.0);
        }
        Observables.set(u, obs);
    }

    
    public void clearAll(){
        if(!Observables.isEmpty()){
            for(int u = 0;u<Observables.size();u++){
                clearObs(u);
            }
        }
    }
    
    public double getStdDevObs(int i){
        return getExpectationObs(i,2, 1)-getExpectationObs(i,1, 2);
    }
    
    public double getExpectationObs(int i, int u){
        return getExpectationObs(i,u, 1);
    }
    
    public double getExpectationObs(int i, int u, int power){
        ArrayList<Double> obs = Observables.get(i);
        double expect = obs.get(u)/updatesObs.get(i);
        if(u==0){expect*=updatesObs.get(i);}
        if(power > 1){
            expect = Math.pow(expect, power);
        }
        return expect;
    }   
}

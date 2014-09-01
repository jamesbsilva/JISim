package Backbone.Util;
/**
 * 
 *      @(#)  MeasureIsingSystem
 */  

import java.util.ArrayList;

/**  
 *   Class for managing how measurements of system variables and observables are made.
 *  
 *  <br>
 * 
 * @author      James Silva <jbsilva @ bu.edu>                 
 * @since       2012-08    
 */
public class MeasureIsingSystem {
    private ArrayList<Double> eAccum;
    private ArrayList<Long> mAccum;
    private ArrayList<Long> mAccumFree;
    private ArrayList<Long> mAccumAbsFree;
    private ArrayList<Long> mAccumAbs;
    private ArrayList<ArrayList<Double>> Observables;
    private Long updatesE;
    private Long updatesM;
    private ArrayList<Long> updatesObs;
    private int N;
    private int fixedVal=-1;
    private short largestOrder = 7;
    private int Nfixed;
    private double temperature;
    
    public MeasureIsingSystem(double temp, int n){
        eAccum = new ArrayList<Double>();
        mAccum = new ArrayList<Long>();
        mAccumFree = new ArrayList<Long>();
        mAccumAbsFree = new ArrayList<Long>();
        mAccumAbs = new ArrayList<Long>();
        Observables = new ArrayList<ArrayList<Double>>();
        updatesObs = new ArrayList<Long>();
        updatesE = 0l;
        updatesM = 0l;
        for(int i=0;i <= largestOrder;i++){
            eAccum.add(i,0.0);
            mAccum.add(i,0l);
            mAccumAbs.add(i,0l);
            mAccumFree.add(i,0l);
            mAccumAbsFree.add(i,0l);
        }
        temperature = temp;
        N = n;
    }
    
    public void changeTemp(double temp){
        temperature = temp;
    }
    
    public void setN(int n, int nfix){
        Nfixed = nfix;
        N = n-nfix;
    }
    
    public void setFixedVal( int nfix){
        fixedVal = nfix;
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

    /**
    *   calcFourthCumulantM should return the fourth order cumulant of the magnetization.
    * 
    */
    public double calcFourthCumulantM(){
        double cum;
        double numerator = (double) (getExpectationM(4)+getExpectationM(1,4)
                -8*(getExpectationM(3)*getExpectationM(1))+6*(getExpectationM(2)*getExpectationM(1,2)) );
        double denominator = (double) ( 3*(
                getExpectationM(2,2)+getExpectationM(1,4)
                -2*(getExpectationM(1,2)*getExpectationM(2))) );
       cum = 1-(numerator/denominator);
       return cum;
    }
    
    /**
    *   calcFourthCumulantM should return the fourth order cumulant of the magnetization.
    * 
    */
    public double calcFourthCumulantMsymm(){
        double cum;
        double numerator = (double) (getExpectationM(4));
        double denominator = (double) ( 3*(getExpectationM(2,2)));
        cum = 1-(numerator/denominator);
        return cum;
    }
    
    /**
    *   calcFourthCumulantE should return the fourth order cumulant of the energy.
    * 
    */
    public double calcFourthCumulantE(){
        double cum;
        double numerator = getExpectationE(4)+getExpectationE(1,4)
                -8*(getExpectationE(3)*getExpectationE(1))+6*(getExpectationE(2)*getExpectationE(1,2)) ;
        double denominator = 3*(
                getExpectationE(2,2)+getExpectationE(1,4)
                -2*(getExpectationE(1,2)*getExpectationE(2)));
       cum = 1-(numerator/denominator);
       return cum;
    }

    /**
    *   calcFourthCumulantE should return the fourth order cumulant of the energy.
    * 
    */
    public double calcFourthCumulantEsymm(){
        double cum;
        double numerator = getExpectationE(4);
        double denominator = 3*(Math.pow(getExpectationE(2),2));
       cum = 1-(numerator/denominator);
       return cum;
    }

    /**
    *   calcBinderRatio should return the binder ratio.
    * 
    */
    public double calcBinderRatio(){
       return (mAccum.get(2) /updatesM)/(Math.pow((mAccumAbs.get(1) /updatesM),2.0));
    }
    /**
    *   getSpecificHeat should return the current specific heat.
    * 
    */
    public double getSpecificHeat() {
        double heatCapacity = getExpectationE(2)-getExpectationE(1,2);
        heatCapacity = heatCapacity/(temperature*temperature);
        return(heatCapacity/N);
    }
    
    /**
    *   getSusceptibility should return the current susceptibility per site.
    * 
    */
    public double getSusceptibility() {
        return (getExpectationM(2)-getExpectationM(1,2))/(temperature*N);
    }

    /**
    *   getSusceptibility should return the current susceptibility.
    * 
    */
    public double getSusceptibilityFreeM() {
        return (getExpectationFreeM(2)-getExpectationFreeM(1,2))/(temperature*(N-Nfixed));
    }
    
    public void updateM(int m){
        updateM(m,true);
    }
    
    public void updateM(int m, boolean update){
        if(update){updatesM++;}
        int mFree = m-fixedVal*Nfixed;
        mAccum.set(0,(long)m);
        mAccumAbs.set(0,Math.abs((long)m));
        mAccumAbsFree.set(0,Math.abs((long)mFree));
        mAccumFree.set(0,(long)mFree);    
        for(short i=1; i<=largestOrder;i++){
            mAccum.set(i, ((long)Math.pow((double)m,i)+mAccum.get(i)));
            mAccumAbs.set(i, Math.abs((long)Math.pow((double)m,i))+mAccumAbs.get(i));
            mAccumFree.set(i, (long)Math.pow((double)(mFree),i)+mAccumFree.get(i));
            mAccumAbsFree.set(i, Math.abs((long)Math.pow((double)(mFree),i))+mAccumAbsFree.get(i));
        }
    }

    public void updateE(double e){
        updateE(e,true);
    }
    
    public void updateE(double e, boolean update){
        if(update){updatesE++;}
        if(Math.abs(e) > 4){e = e/N;}
        eAccum.set(0,e);
        for(short i=1; i<=largestOrder;i++){
            eAccum.set(i, (Math.pow(e,i)+eAccum.get(i)));
        } 
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

    public void clearE(){
        updatesE = 0l;
        eAccum.set(0,0.0);
        for(short i=1; i<=largestOrder;i++){
            eAccum.set(i, 0.0);
        }
    }
    
    public void clearM(){
        updatesM = 0l;
        mAccum.set(0,0l);
        mAccumAbs.set(0,0l);
        mAccumFree.set(0,0l);
        mAccumAbsFree.set(0,0l);
        for(short i=1; i<=largestOrder;i++){
            mAccum.set(i, 0l);
            mAccumAbs.set(i, 0l);
            mAccumFree.set(i, 0l);
            mAccumAbsFree.set(i, 0l);
        }
    }
    
    public void clearAll(){
        if(!Observables.isEmpty()){
            for(int u = 0;u<Observables.size();u++){
                clearObs(u);
            }
        }
        clearM();
        clearE();
    }
    
    public double getExpectationM(int u){
        return getExpectationM(u, 1);
    }
    
    public double getExpectationM(int u, int power){
        double expect = ((double)mAccum.get(u))/((double)updatesM);
        if(u==0){expect*=updatesM;}
        if(power > 1){
            expect = Math.pow(expect, power);
        }
        //System.out.println("MEASURESYSTEM MExpect:"+expect+"   powerOfExpect: "+u+"    powerToMult: "+power);
        return expect;
    }

    public double getExpectationFreeM(int u){
        return getExpectationFreeM(u, 1);
    }
    
    public double getExpectationFreeM(int u, int power){
        double expect = ((double)mAccumFree.get(u))/((double)updatesM);
        if(u==0){expect*=updatesM;}
        if(power > 1){
            expect = Math.pow(expect, power);
        }
        expect = expect/Math.pow(N-Nfixed,power+u);
        return expect;
    }
    
    public double getExpectationAbsFreeM(int u){
        return getExpectationAbsFreeM(u, 1);
    }
    
    public double getExpectationAbsFreeM(int u, int power){
        double expect = ((double)mAccumAbsFree.get(u))/((double)updatesM);
        if(u==0){expect*=updatesM;}
        if(power > 1){
            expect = Math.pow(expect, power+u);
        }
        expect = expect/Math.pow(N-Nfixed,power+u);
        return expect;
    }
    
    public double getExpectationAbsM(int u){
        return getExpectationAbsM(u, 1);
    }
    
    public double getExpectationAbsM(int u, int power){
        double expect = mAccumAbs.get(u)/updatesM;
        if(u==0){expect*=updatesM;}
        if(power > 1){
            expect = Math.pow(expect, power);
        }
        expect = expect/Math.pow(N,power+u);
        return expect;
    }
    
    public double getExpectationE(int u){
        return getExpectationE(u, 1);
    }
    
    public double getExpectationE(int u, int power){
        double expect = eAccum.get(u)/updatesE;
        if(u==0){expect*=updatesE;}
        if(power > 1){
            expect = Math.pow(expect, power);
        }
        return expect;
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

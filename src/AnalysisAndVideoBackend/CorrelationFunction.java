package AnalysisAndVideoBackend;

import Backbone.System.LatticeMagInt;
import static java.lang.Math.PI;
import java.util.ArrayList;
import java.util.Collections;

public class CorrelationFunction{
    ArrayList<Integer> corrValXY;
    ArrayList<Integer> corrValX;
    ArrayList<Double> corrValXYt;
    ArrayList<Double> corrValXt;
    ArrayList<Double> corrValYt;
    ArrayList<Double> corrFunc;
    ArrayList<Double> corrFuncT;
    ArrayList<Integer> valTime;
    ArrayList<Integer> corrValY;
    ArrayList<Double> corrLen;
    ArrayList<Integer> corrNspacial;
    ArrayList<Integer> corrNspacialT;
    int NinTime = 0;
    int L;                 
    int R;
    
    public CorrelationFunction(int len, int range) {
        R = range;
        L = len;
        initializeArrays();
    }

    public void initializeArrays(){
        corrLen = createCorrArr(L);
        corrValXY =  new ArrayList<Integer>(corrLen.size());
        corrValX =  new ArrayList<Integer>(corrLen.size());
        corrValY =  new ArrayList<Integer>(corrLen.size());
        corrNspacial =  new ArrayList<Integer>(corrLen.size());
        corrNspacialT =  new ArrayList<Integer>(corrLen.size());
        corrValXYt =  new ArrayList<Double>(corrLen.size());
        corrValXt =  new ArrayList<Double>(corrLen.size());
        corrValYt =  new ArrayList<Double>(corrLen.size());
        valTime =  new ArrayList<Integer>(L*L);
        corrFunc =  new ArrayList<Double>(corrLen.size());
        corrFuncT =  new ArrayList<Double>(corrLen.size());
        NinTime = 0;
    }
    
    public ArrayList<Double> createCorrArr(int len){
        ArrayList<Double> carr = new ArrayList<Double>();
        double r;
        boolean newVal = true;
        int ind = 0;
        for(int u = 0; u < (int)(L/2)+1;u++){
            for(int v = 0; v < (int)(L/2);v++){
                r = Math.sqrt(u*u+v*v);
                while(ind < carr.size() && newVal){
                    if(r == carr.get(ind)){
                        newVal = false;
                    }
                    ind++;
                }    
                if(newVal){
                    carr.add(r);
                }
            }
        }
        //sort the list
        Collections.sort(carr);
        return carr;
    }
    
    public void processLatticeSpacial(LatticeMagInt lat){
        int val = 0;
        for(int v = 0; v < L;v++){
        for(int u = 0; u < L;u++){
            val  = lat.getValue(u, v);
            for(int y = 0; y < L;y++){
            for(int x = 0; x < L;x++){
                if(x == u && y == v){}else{
                    insertCorrVal(u, v, x, y, val, lat.getValue(x, y));
                }
            }}        
        }}
        for(int u = 0; u < corrFunc.size();u++){
            corrFunc.set(u,(corrValXY.get(u)/corrNspacial.get(u) - 
                    (corrValX.get(u)*corrValY.get(u)/Math.pow(corrNspacial.get(u),2) )));       
        }
    }
    
    public void processLatticeTime(){
        int val = 0;
        for(int v = 0; v < L;v++){
        for(int u = 0; u < L;u++){
            val  = valTime.get(u+v*L)/NinTime;
            for(int y = 0; y < L;y++){
            for(int x = 0; x < L;x++){
                if(x == u && y == v){}else{
                    insertCorrVal(u, v, x, y, val, valTime.get(x+y*L)/NinTime);
                }
            }}        
        }}
        for(int u = 0; u < corrFunc.size();u++){
            corrFuncT.set(u,(corrValXYt.get(u)/corrNspacialT.get(u) - 
                    (corrValXt.get(u)*corrValYt.get(u)/Math.pow(corrNspacialT.get(u),2) )));       
        }
    }
    
    public void addLatticeData(LatticeMagInt lat){
        int val = 0;
        for(int v = 0; v < L;v++){
        for(int u = 0; u < L;u++){
            val  = lat.getValue(u, v);
            valTime.set(u+v*L,val+valTime.get(u+v*L));    
        }}
        NinTime++;
    }
    
    public void insertCorrVal(int u, int v, int x, int y, double x1, double x2){
        double r = Math.sqrt((x-u)*(x-u)+(y-v)*(y-v));
        int ind = 0;
        boolean valNotFound = true;
        while(ind < corrLen.size() && valNotFound){
            ind++;
        }
        corrValXt.set(ind,x1+corrValX.get(ind));
        corrValYt.set(ind,x2+corrValY.get(ind));
        corrValXYt.set(ind,x1*x2+corrValXY.get(ind));
        corrNspacialT.set(ind,1+corrNspacialT.get(ind));    
    }
    
    public void insertCorrVal(int u, int v, int x, int y, int x1, int x2){
        double r = Math.sqrt((x-u)*(x-u)+(y-v)*(y-v));
        int ind = 0;
        boolean valNotFound = true;
        while(ind < corrLen.size() && valNotFound){
            ind++;
        }
        corrValX.set(ind,x1+corrValX.get(ind));
        corrValY.set(ind,x2+corrValY.get(ind));
        corrValXY.set(ind,x1*x2+corrValXY.get(ind));
        corrNspacial.set(ind,1+corrNspacial.get(ind));
    }
    
    
    
}
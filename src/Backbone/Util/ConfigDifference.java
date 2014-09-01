package Backbone.Util;

/**
*    @(#)   ConfigDifference
*/

import Backbone.System.Lattice;
import Backbone.System.LatticeDbl;
import Backbone.System.LatticeMagInt;
import java.util.ArrayList;

/**
 *      ConfigDifference determines the difference between simulation steps.
 * 
 * @author James B. Silva <jbsilva@bu.edu>
 * @since July 2013
 */
public class ConfigDifference {
    private ArrayList<Double> lastConf;
    private ArrayList<Double> diffConf;
    private ArrayList<Double> diffSumConf;
    private int sumN = 0;
    private double sumMax = -10000;
    private double sumMin;
    private double minDiff = 1000000;
    private double maxDiff = -10000;
    private boolean binaryAvg = true;
    private int N;
    private boolean initConfSet = false;
    
    public ConfigDifference(Lattice s){
        N = s.getN();
        lastConf = new ArrayList<Double>();
        diffConf = new ArrayList<Double>();
        diffSumConf = new ArrayList<Double>();
        for(int u = 0; u < N;u++){
            lastConf.add(0.0);
            diffConf.add(0.0);
            diffSumConf.add(0.0);      
        }
        sumMin =  100000;
        sumN =  0;
    }
    
    public void update(LatticeMagInt s){
        double currVal;
        for(int u = 0; u < N; u++){
            currVal = (double)s.getValue(u);
            diffConf.set(u, currVal - lastConf.get(u));
            lastConf.set(u, currVal);
            double sumVal = diffConf.get(u);
            if(binaryAvg){
                sumVal = (sumVal == 0) ? 0 : 1;
            }
            diffSumConf.set(u, sumVal+diffSumConf.get(u));
            //System.out.println(" sum now | "+diffSumConf.get(u)+"     max| "+sumMax+"    diff| "+sumVal);
            if( diffSumConf.get(u) > sumMax ){sumMax = diffSumConf.get(u);}
            if( diffSumConf.get(u) < sumMin ){sumMin = diffSumConf.get(u);}
            if( diffConf.get(u) > maxDiff ){maxDiff = diffConf.get(u);}
            if( diffConf.get(u) < minDiff ){minDiff = diffConf.get(u);}
        }
        
        sumN++;    
    }
    
    public void update(LatticeDbl s){
        for(int u = 0; u < N; u++){
            diffConf.set(u, s.getValue(u) - lastConf.get(u));
            lastConf.set(u, s.getValue(u));
            double sumVal = diffConf.get(u);
            if(binaryAvg){
                sumVal = (sumVal == 0) ? 0 : 1;
            }
            diffSumConf.set(u, sumVal+diffSumConf.get(u));
            if( diffSumConf.get(u) > sumMax ){sumMax = diffSumConf.get(u);}
            if( diffSumConf.get(u) < sumMin ){sumMin = diffSumConf.get(u);}
            if( diffConf.get(u) > maxDiff ){maxDiff = diffConf.get(u);}
            if( diffConf.get(u) < minDiff ){minDiff = diffConf.get(u);}
        }
        sumN++;
    }
    
    public ArrayList<Double> getConfigDiff(){
        return diffConf;
    }
    public double getMaxAvgSum(){
        return sumMax;
    }
    
    public double getMinAvgSum(){
        return sumMin;
    }
    public double getMinDiff(){
        return minDiff;
    }
    public double getMaxDiff(){
        return minDiff;
    }
    
    public double getConfigN(){
        return sumN;
    }
    
    public void setBinaryModeOnOff(boolean bl){
        binaryAvg = bl;
    }
    
    public ArrayList<Double> getConfigDiffAvgSum(){
        return diffSumConf;
    }
}

package AnalysisAndVideoBackend;

import Backbone.System.LatticeMagInt;
import Backbone.System.SimSystem;
import Backbone.Visualization.SquareDoubleLattice2D;
import java.awt.image.BufferedImage;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.atan;
import static java.lang.Math.floor;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import org.apache.commons.math3.complex.Complex;


public class StructureFactorSlow{
    private double accH;
    private double accV;
    private double accTriBase1;
    private double accTriBase2;
    private int L;
    private int Lq=L;
    private int Geo=4;
    private int accN=0;
    private double lastSFv=0;
    private double lastSFh=0;
    private double lastSFtri=0;
    private double lastSFtri2=0;
    private double qfactor;
    private double maxSF = 0;
    private ArrayList<Double> basisVecs;
    private boolean visualize = true;
    private boolean initializedImg = false;
    private SquareDoubleLattice2D vis;
    private SquareDoubleLattice2D visCir;
    private int circularBins = 128;
    private int radialBins = 128;
    private double[] circularSF;
    private double[] radialSF;
    private double[] sFactor;
    private double[] circularSFacc;
    private double[] radialSFacc;
    private double[] sFactorAcc;
    private int nMeasured;
    private ArrayList<Double> radialVals;
    private HashMap<Double,Integer> radialMap;
    
    public StructureFactorSlow(int lin, int g){
        L = lin;
        Lq = L;
        Geo = g;
        accN = 0;
        accTriBase1 = 0;
        accTriBase2 = 0;
        accH = 0;
        accV = 0;
        sFactor = new double[L*Lq];
        circularBins = L*Lq/16;
        circularSF = new double[circularBins];
        basisVecs = new ArrayList<Double>();
        qfactor = Math.PI*2/L;
        if(Geo ==6){
            addBasisVecs(1, 0, 0, 1);
            //addBasisVecs(0, 1, 0.5, Math.sqrt(3)/2.0);
            //addBasisVecs(1, 0, 0.5, Math.sqrt(3)/2.0);
            //addBasisVecs(0, 1, Math.sqrt(3)/2.0,-1.0);
        }else{
            addBasisVecs(1, 0, 0, 1);
        }
        vis = new  SquareDoubleLattice2D(Lq);
        vis.initializeImg();
        visCir = new  SquareDoubleLattice2D(Lq/4);
        visCir.initializeImg();
        makeSkeletonRadialSF(Lq);
        sFactorAcc = new double[L*Lq];
        circularSFacc = new double[circularBins];
        radialSFacc = new double[radialBins];
        nMeasured = 0;
        System.out.println("FourierTransform|     L: "+L+"      Geo: "+Geo);
    }
    
    public void makeSkeletonRadialSF(int Lin){
        int max = 1+(int)(Lin/2.0);
        radialVals = new ArrayList<Double>();
        for(int i = 0; i < Lq;i++){for(int j = 0; j < Lq;j++){
            boolean newVal = false;
            newVal = notInList(radialVals,Math.sqrt(getXrCentered(i,j)*getXrCentered(i,j)+getYrCentered(i,j)*getYrCentered(i,j)));
            if(newVal){radialVals.add(Math.sqrt(getXrCentered(i,j)*getXrCentered(i,j)+getYrCentered(i,j)*getYrCentered(i,j)));}
        }}
        // sort Possible
        Collections.sort(radialVals);
        radialMap = new HashMap<Double,Integer>();
        for(int i=0;i < radialVals.size();i++){
            radialMap.put(radialVals.get(i), i);
        }
        radialBins = radialVals.size();
    }
    
    public boolean notInList(ArrayList<Double> list, double val){
        int ind = 0;
        boolean searching = true;
        boolean found = false;
        while(ind < list.size() && searching){
            if(list.get(ind) == val){found = true;}
            ind++;
        }
        return !found;
    }
    
    public void calc2DSF(int[] spins){
    double argX = 0;double argY = 0;double arg = 0;
        int a = 0;        
        int offset = 0;//(int)(L/2.0);
        int lhalf = (int)((double)Lq/(double)(2.0));
        for(int i = (-lhalf); i < (lhalf);i++){for(int j = (-lhalf); j < (lhalf);j++){
            Complex sum = new Complex(0,0);
            Complex temp;
            // iterate over the lattice
            for(int u = 0; u < L;u++){for(int v = 0; v < L;v++){
                // calculate components for dot product
                // G1*nx+G2*nx where G1 and G2 are the basis vectors 
                argX = i*basisVecs.get(a)+j*basisVecs.get(a+2);
                argY = i*basisVecs.get(a+1)+j*basisVecs.get(a+3);
                // calculate distance in functions getXr and getYr to allow for more complicated lattices
                arg = dotproduct(argX,argY,getXrCentered(u,v),getYrCentered(u,v));
                // qfactor is 2 pi/L
                arg = arg*qfactor;
                // i * q* r
                temp = new Complex(0,arg);
                //exp(iqr)*s
                temp = temp.exp();
                temp = temp.multiply(spins[u+v*L]);
                sum = sum.add(temp);
            }}
            // take |exp(iqr)s|^2/N
            arg = Math.pow(getAbsoluteComplex(sum),2.0)/spins.length;
            //sFactor[(((i+offset+(int)(1.0*Lq/2))%L)+((j+offset)%L)*Lq)] = arg;
            //sFactorAcc[(((i+offset)%L)+((j+offset)%L)*Lq)] += arg;
            //sFactor[(((i+offset+(int)(1.0*Lq/2))%L)+((j+offset+(int)(1.0*Lq/2))%L)*Lq)] = arg;
            sFactor[(((i+offset+(int)(1.0*Lq/2))%L)+((j+offset+(int)(1.0*Lq/2))%L)*Lq)] = arg;
            sFactorAcc[(((i+offset+(int)(1.0*Lq/2))%L)+((j+offset+(int)(1.0*Lq/2))%L)*Lq)] += arg;
            if(arg > maxSF){maxSF = arg;
                System.out.println("sFactorSlow|    i:"+i+"   j: "+j+"    sf: "+arg+"   N: "+spins.length);
            }           
        }}
        vis.spinDraw(sFactor,(maxSF != 0) ? maxSF : 1);
        initializedImg = true;
        calcSFCircular();
        calcSFradial();
        nMeasured++;
    }
    public void calc2DSF(SimSystem simsys){
        calc2DSF(simsys,0);
    }
    
    public void calc2DSF(SimSystem simsys,int time){
        if(simsys.getClass().getName().contains("Lattice")){
            calc2DSF((LatticeMagInt)simsys,time);
        }else{
            System.err.println("SlowFourierTranform |  System needs to be a lattice  : "+simsys.getClass().getName());
        }
    }
    
    public void calc2DSF(LatticeMagInt spins, int time){
        double argX = 0;double argY = 0;double arg = 0;
        int a = 0;        
        int offset = 0;//(int)(L/2.0);
        int lhalf = (int)((double)Lq/(double)(2.0));
        for(int i = (-lhalf); i < (lhalf);i++){for(int j = (-lhalf); j < (lhalf);j++){
            Complex sum = new Complex(0,0);
            Complex temp;
            // iterate over the lattice
            for(int u = 0; u < L;u++){for(int v = 0; v < L;v++){
                // calculate components for dot product
                // G1*nx+G2*nx where G1 and G2 are the basis vectors 
                argX = i*basisVecs.get(a)+j*basisVecs.get(a+2);
                argY = i*basisVecs.get(a+1)+j*basisVecs.get(a+3);
                // calculate distance in functions getXr and getYr to allow for more complicated lattices
                arg = dotproduct(argX,argY,getXrCentered(u,v),getYrCentered(u,v));
                // qfactor is 2 pi/L
                arg = arg*qfactor;
                // i * q* r
                temp = new Complex(0,arg);
                //exp(iqr)*s
                temp = temp.exp();
                temp = temp.multiply(spins.getValue(u, v));
                sum = sum.add(temp);
            }}
            // take |exp(iqr)s|^2/N
            arg = Math.pow(getAbsoluteComplex(sum),2.0)/spins.getN();
            //sFactor[(((i+offset+(int)(1.0*Lq/2))%L)+((j+offset)%L)*Lq)] = arg;
            //sFactorAcc[(((i+offset)%L)+((j+offset)%L)*Lq)] += arg;
            //sFactor[(((i+offset+(int)(1.0*Lq/2))%L)+((j+offset+(int)(1.0*Lq/2))%L)*Lq)] = arg;
            sFactor[(((i+offset+(int)(1.0*Lq/2))%L)+((j+offset+(int)(1.0*Lq/2))%L)*Lq)] = arg;
            sFactorAcc[(((i+offset+(int)(1.0*Lq/2))%L)+((j+offset+(int)(1.0*Lq/2))%L)*Lq)] += arg;
            if(arg > maxSF){maxSF = arg;
                System.out.println("sFactorSlow|    i:"+i+"   j: "+j+"    sf: "+arg+"   N: "+spins.getN());
            }           
        }}
        vis.spinDraw(sFactor,(maxSF != 0) ? maxSF : 1);
        initializedImg = true;
        calcSFCircular();
        calcSFradial();
        nMeasured++;
    }
    
    public BufferedImage getVisImg(){
        return vis.getImageOfVis();
    }
    
    
    public BufferedImage getVisCirImg(){
        return visCir.getImageOfVis();
    }
    
    public double standardizeAngle(double x, double y,double angle){
        if (x <= 0 && y >=0)
                angle = PI -angle;
        else if (x <=0 && y<=0)
                angle += PI;
        else if (x >= 0 && y <= 0)
                angle = 2*PI - angle;
        if (angle == 2*PI)
                angle = 0;
        return angle;
    }
    
    // Radial Average
    public void calcSFCircular(){
        double maxCir = 0;
        // create array with [0,2pi] binned into circularBins amount of bins
        circularSF = new double[circularBins];
        double[] cN = new double[circularBins];
        for(int i = 0; i < Lq;i++){for(int j = 0; j < Lq;j++){
            // Get angle relative to center of lattice
            double angle = atan((double)abs(getXrCentered(i,j))/(double)abs(getYrCentered(i,j)));        
            // Convert angle to range [0,2pi]
            angle = standardizeAngle(getXrCentered(i,j),getYrCentered(i,j), angle);
            // Convert angle to proper bin
            int ind = (int)floor(angle*circularBins/(2*PI));
            // Add contribution from angle
            circularSF[ind] += sFactor[i+j*L];
            cN[ind]++;
            // keep track of largest value
            if(circularSF[ind] > maxCir){maxCir = circularSF[ind];}
        }}
        for(int i = 0; i < circularBins;i++){
            if(cN[i] != 0){
                circularSF[i] = circularSF[i]/cN[i];
            }else{
                circularSF[i] = 0;
            }
            circularSFacc[i] += circularSF[i];
        }
        visCir.spinDraw(circularSF,(maxCir != 0) ? maxCir : 1);
    }
    
    //Circular Average
    public void calcSFradial(){
        double maxRad = 0;double rad;
        radialSF = new double[radialBins];
        double[] rN = new double[radialBins];
        
        for(int i = 0; i < Lq;i++){for(int j = 0; j < Lq;j++){
            // Get distance from center
            rad= Math.sqrt(getXrCentered(i,j)*getXrCentered(i,j)+getYrCentered(i,j)*getYrCentered(i,j));
            // Convert distance into index in lattice using precomputed hashmap
            int ind = radialMap.get(rad);
            // Add contribution for this r
            radialSF[ind] += sFactor[i+j*L];
            rN[ind]++;
            // keep track of largest value
            if(radialSF[ind] > maxRad){maxRad = radialSF[ind];}
        }}
        for(int i = 0; i < radialBins;i++){
            if(rN[i] != 0){
                radialSF[i] = radialSF[i]/rN[i];
            }else{
                radialSF[i] = 0;
            }
            radialSFacc[i] += radialSF[i];
        }
        // visualization
        //visRad.spinDraw(radialSF,(maxRad != 0) ? maxRad : 1);
    }
    
    public void addBasisVecs(double a, double b,double c, double d){
        basisVecs.add(a);
        basisVecs.add(b);
        basisVecs.add(c);
        basisVecs.add(d);
    }
    
    public double dotproduct(double a,double b, double c, double d){
        return (a*c+b*d);
    }
    
    private double getXrCentered(int u,int v){
        int cenInd = (int)(L/2.0);
        double xcen = getXr(cenInd,cenInd);
        return (getXr(u,v)-xcen);
    }
    
    private double getYrCentered(int u,int v){
        int cenInd = (int)(L/2.0);
        double ycen = getYr(cenInd,cenInd);
        return (getYr(u,v)-ycen);
    }
    
    private double getXr(int u,int v){
        if(Geo ==  6){
            return getXcoordTri(u,v); 
        }else if(Geo ==  3){
            return getXcoordHoney(u,v); 
        }else{
            return u;
        }
    }
    
    private double getYr(int u,int v){
        if(Geo ==  6){
            return getYcoordTri(u,v); 
        }else if(Geo ==  3){
            return getYcoordHoney(u,v); 
        }else{
            return v;
        }
    }
    
    public double getXcoordTri(int i, int j){
        return (j%2 == 0) ? (double)i : (double)i + 0.5;
    }
    
    public double getYcoordTri(int i, int j){
        return ((double)j)*Math.sqrt(3)/2.0;
    }

    public double getXcoordHoney(int i, int j){
        double scaleX = 1.0;
        double offset = (j%2 == 0) ? Math.floor(i/2)*scaleX : (Math.floor(i/2)+1.5)*scaleX;
        return (((double)i)*scaleX)+offset;
    }
    
    public double getYcoordHoney(int i, int j){
        double scaleY = 1.0;
        return ((double)j)*Math.sqrt(3)*scaleY/2.0;
    }
    
    public double[] getSF(){
        return sFactor;
    }
 
    public double[] avgArray(double[] arr){
        if(nMeasured == 0){return arr;}
        double[] avg = new double[arr.length];
        for(int u = 0;u < arr.length;u++){
            avg[u] = arr[u]/nMeasured;
        }
        return avg;
    }
    public double[] getSFavg(){
        return avgArray(sFactor);
    }
    
    public double[] getSFcirAvg(){
        return avgArray(circularSF);
    }
    public double[] getSFradAvg(){
        return avgArray(radialSF);
    }
    
    public double[] getSFcir(){
        return circularSF;
    }
    public double[] getSFrad(){
        return radialSF;
    }
    public ArrayList<Double> getRadVals(){
        return radialVals;
    }
    public double getAbsoluteComplex(Complex cnumber){
        return ((cnumber.multiply(cnumber.conjugate())).sqrt()).getReal();
    }
    
    // test the class
    public static void main(String[] args) {
        Complex test = new Complex(3,5); 
        System.out.println("test re: "+test.getReal()+"     comp: "+test.getImaginary());
        Complex test2 = new Complex(1,2); 
        test = test.add(test2);
        System.out.println("test re: "+test.getReal()+"     comp: "+test.getImaginary());
        test = test.multiply(2);
        System.out.println("test re: "+test.getReal()+"     comp: "+test.getImaginary());
        test = (test.multiply(test.conjugate())).sqrt();
        System.out.println("test re: "+test.getReal()+"     comp: "+test.getImaginary());
        Complex test3 = new Complex(0,4); 
        test = test3.exp();
        System.out.println("test re: "+test.getReal()+"     comp: "+test.getImaginary());
    }
    
}
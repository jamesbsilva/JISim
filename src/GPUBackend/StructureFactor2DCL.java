package GPUBackend;
/**
 * @(#)  StructureFactor2DCL
 */

// imports
import Backbone.System.LatticeMagInt;
import Backbone.System.SimSystem;
import Backbone.Visualization.SquareDoubleLattice2D;
import java.awt.Color;
import java.awt.image.BufferedImage;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.atan;
import static java.lang.Math.floor;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import scikit.opencl.*;

/**
*      StructureFactor2DCL  
*
*   <br>
* 
* @author James B. Silva <jbsilva@bu.edu>
* @date June 28, 2013
*/
public class StructureFactor2DCL {
    private CLHelper clhandler;
    private String SFKernel = "structurefactor2d";
    private String intUpdateKernel = "update_int_buffer";
    private boolean visualize = true;
    private boolean initializedImg = false;
    private SquareDoubleLattice2D vis;
    private HashMap<Integer,KernelArgInfo> kernelsInfo;
    private HashMap<KernelArgInfo,Integer> indexInfo;
    private HashMap<Integer,int[]> workgroupInfo; 
    private ArrayList<Float> testSet;
    private int currentBuffer=0;    
    private int nBuffers=0;
    private int L;
    private float Lq;
    private double[] sFactorAvg;
    private int nMeasurements = 0;
    private int geo;
    private int maxThreads=129;
    private int circularBins = 128;
    private int radialBins = 128;
    private int[] peaksLoc = null;
    private double[] circularSF;
    private double[] radialSF;
    private double[] circularSFacc;
    private double[] radialSFacc;
    private ArrayList<Double> radialVals;
    private HashMap<Double,Integer> radialMap;
    
    public StructureFactor2DCL(int lin, float lqin, int ge){
        this(null,lin,lqin,ge);
    }
    
    public StructureFactor2DCL(CLHelper cl, int lin, float lqin, int ge){
        if(cl == null){
            clhandler = new CLHelper();
            clhandler.initializeOpenCL("GPU");
            cl = clhandler;
        }
        L = lin;
        geo =ge;
        Lq = lqin;
        
        vis = new  SquareDoubleLattice2D((int)(Lq));
        vis.initializeImg();
    
        kernelsInfo = new HashMap<Integer,KernelArgInfo>();
        indexInfo = new HashMap<KernelArgInfo,Integer>();
        workgroupInfo = new HashMap<Integer,int[]>();
        maxThreads = cl.getCurrentDevice1DMaxWorkItems();
        
        divideWork((int)(Lq*Lq));
        
        clhandler = cl;
        clhandler.createKernel("", SFKernel);
        // spins buffer
        clhandler.createIntBuffer(SFKernel, 0, L*L, 0, "rw");
        // structure factor buffer
        clhandler.createFloatBuffer(SFKernel, 0, (int)(Lq*Lq), 0);
        //error buffer
        clhandler.createFloatBuffer(SFKernel, 1,400, 0,"rw");
        
        setSimpleArgs();
        clhandler.setKernelArg(SFKernel);
        initializeUpdateKernel();
    
        makeSkeletonRadialSF((int)Lq);
        circularBins = (int)(L*Lq/16);
        
        sFactorAvg = new double[(int)(Lq*Lq)];
        circularBins = (int)(L*Lq/16);
        circularSF = new double[circularBins];
        circularSFacc = new double[circularBins];
        radialSFacc = new double[radialBins];

    }

    public void resetAccumulators(){
        sFactorAvg = new double[(int)(Lq*Lq)];
        circularSF = new double[circularBins];
        circularSFacc = new double[circularBins];
        radialSFacc = new double[radialBins];
        nMeasurements = 0;
    }
    
    private void initializeUpdateKernel(){
        clhandler.createKernel("", intUpdateKernel);
        clhandler.copyIntBufferAcrossKernel(SFKernel, 0, intUpdateKernel, 0);
    }
    
    private void setSimpleArgs(){
        clhandler.setIntArg(SFKernel, 0, (int)Lq);
        clhandler.setIntArg(SFKernel, 1, L);
        clhandler.setIntArg(SFKernel, 2, geo);
        clhandler.setFloatArg(SFKernel, 0, (float)Math.PI);
    }
    
    public void setConfigurationAndCalcSF(ArrayList<Integer> spins){
        setConfiguration(spins);
        calcSFexistingConf();
    }
    
    public void setConfigurationAndCalcSF(SimSystem sys){
        setConfigurationAndCalcSF((LatticeMagInt)sys);
    }
    public void setConfigurationAndCalcSF(LatticeMagInt spins){
        ArrayList<Integer> arr =  new ArrayList<Integer>();
        int v;int L = spins.getLength();
        for(int u = 0; u < spins.getN();u++){
            v =  ((int)((double)u/((double)L)))%L;
            arr.add(spins.getValue(u%L,v));
        }
        setConfiguration(arr);
        calcSFexistingConf();
    }
    
    public void setConfiguration(ArrayList<Integer> spins){
        int[] workInfo = workgroupInfo.get(0);
        // fill with test set and update
        clhandler.setIntBuffer(intUpdateKernel, 0, workInfo[0], spins, true);
       // clhandler.setKernelArg(intUpdateKernel, true); 
        clhandler.runKernel(intUpdateKernel,workInfo[0],workInfo[1]);
    }
    
    
    public BufferedImage getVisImg(){
        return vis.getImageOfVis();
    }
    
    public void calcSFexistingConf(){
        int[] workInfo = workgroupInfo.get(0);
        //clhandler.setKernelArg(SFKernel, true); 
        clhandler.runKernel(SFKernel,workInfo[0],workInfo[1]);
    }
    
    public float[] getSF(){
        float[] sf = clhandler.getFloatBufferAsArray(SFKernel, 0, (int)(Lq*Lq), false);
        float maxSF = (float) clhandler.quickGetLastBufferMax();;
        for(int u = 0; u < (int)(Lq*Lq);u++){
            sFactorAvg[u] += sf[u];
        }
        nMeasurements++;
        calcSFradial(sf);
        calcSFCircular(sf);
        if(visualize){
            double[] temp = sFactorAvg; double max = 0;
            for(int u = 0; u < (int)(Lq*Lq);u++){
                temp[u] = (double) sf[u];
                //if(temp[u] > 100){
                //    System.out.println("MAX AT : "+u+"    VAL: "+temp[u]);
                //    }
                if(temp[u] > max){
                    max = temp[u];
                    //System.out.println("sFactorCL | MAXSF: "+max);
                }
            }
            vis.spinDraw(temp,(max != 0) ? max : 1);
        }
        return sf;
    }
    
    
    public void drawPeaks(){
        for(int u = 0; u < peaksLoc.length;u++){
            vis.updateImg(peaksLoc[u], Color.GREEN);
        }
    }
    
    public void setPeaksLoc(int[]peakLoc){
        peaksLoc = peakLoc;
        if(peaksLoc != null){drawPeaks();}
    }
    
    public double[] getSFavg(){
        if(nMeasurements ==0){return sFactorAvg;}
        double[] temp = sFactorAvg;
        for(int u = 0; u < (int)(Lq*Lq);u++){
            temp[u] /= nMeasurements;
        }
        return temp;
    }
    
    //Circular Average
    public void calcSFradial(float[] sFactor){
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
    
    
    // Radial Average
    public void calcSFCircular(float[] sFactor){
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
    
    public BufferedImage updateAndGetVisImage(){
        float[] sf = clhandler.getFloatBufferAsArray(SFKernel, 0, (int)(Lq*Lq), false);
        float maxSF = (float) clhandler.quickGetLastBufferMax();;
        if(visualize){vis.spinDraw(sf,(maxSF != 0) ? maxSF : 1);}
        if(peaksLoc != null){drawPeaks();}
        return vis.getImageOfVis();
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
        if(geo ==  6){
            return getXcoordTri(u,v); 
        }else if(geo ==  3){
            return getXcoordHoney(u,v); 
        }else{
            return u;
        }
    }
    
    private double getYr(int u,int v){
        if(geo ==  6){
            return getYcoordTri(u,v); 
        }else if(geo ==  3){
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
    
    public double[] avgArray(double[] arr){
        if(nMeasurements == 0){return arr;}
        double[] avg = new double[arr.length];
        for(int u = 0;u < arr.length;u++){
            avg[u] = arr[u]/nMeasurements;
        }
        return avg;
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

    public void setMaxThreads(int u){
        maxThreads=u;
    }

    public void divideWorkWithMaxThreads(int u, int size){
        maxThreads=u;
        divideWork(size);
    }
    
    private void divideWork(int size){
        int i=1;
        int remainder = size;
        
        int local = size;

        // If too many threads than divide threads into local workgroups
        while(local> maxThreads && i < size){
            while(remainder !=0 && i < size){
                i++;
                remainder = size % i; 
            }
            local = (int) (size/i);
            if(local> maxThreads && remainder ==0){
                i++;
                remainder= size % i; 
                local = (int) (size/i);}
        }

        local = (int) (size/i);        
        int[] temp = new int[2];
        temp[0] = size;
        temp[1] = local;
        workgroupInfo.put(nBuffers, temp);
    }
    
    private class KernelArgInfo{
        private String KernelName;
        private int Argument;
        public KernelArgInfo(String k,int u){
            KernelName= k;
            Argument = u;
        }
        public void setInfo(String k,int u){
            KernelName= k;
            Argument = u;
        }
        public int getArg(){
            return Argument;
        }
        public String getKernelName(){
            return KernelName;
        }
    }
    
    public double getLq(){
        return Lq;
    }
    
    public void forceCloseCL(){
        clhandler.closeOpenCL();
    }
    
    // test the class
    public static void main(String[] args) {
        CLHelper clhandler = new CLHelper();
        clhandler.initializeOpenCL("GPU");  
        clhandler.closeOpenCL();
    }
}


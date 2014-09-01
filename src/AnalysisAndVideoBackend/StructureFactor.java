package AnalysisAndVideoBackend;

// Modified version of Kang Liu's structure factor code to work with these classes
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
import scikit.numerics.fft.managed.ComplexDouble2DFFT;

public class StructureFactor{
    ComplexDouble2DFFT fft;
    double[] fftData;       // Fourier transform data
    double maxfft;
    private double sFactor [];
    private boolean initialized = false;
    private boolean visualize = true;
    private boolean initializedImg = false;
    private SquareDoubleLattice2D vis;
    int scale = 1;
    int geo =4;
    int Lp;                 // # elements per side
    double L;               // the actual system length, L = Lp*dx, where dx is lattice spacing
    static double squarePeakValue = 4.4934092;
    static double circlePeakValue = 5.135622302;
    private int circularBins = 128;
    private int radialBins = 128;
    private double[] circularSF;
    private double[] radialSF;
    private double[] circularSFacc;
    private double[] radialSFacc;
    private double[] sFactorAcc;
    private int nMeasured;
    private ArrayList<Double> radialVals;
    private HashMap<Double,Integer> radialMap;
    
    public StructureFactor(int Lp, double L, int gin) {
            this.Lp = Lp;
            this.L = L;
            geo = gin;
            sFactor = new double [Lp*Lp];
            circularBins = (int) L*Lp/16;
            circularSF = new double[circularBins];
            fft = new ComplexDouble2DFFT(Lp, Lp);
            
            fftData = new double[2*Lp*Lp*scale*scale];
            vis = new  SquareDoubleLattice2D((int)(Lp));
            vis.initializeImg();
            makeSkeletonRadialSF((int)Lp);
            sFactorAcc = new double[(int)L*Lp];
            circularSFacc = new double[circularBins];
            radialSFacc = new double[radialBins];
            nMeasured = 0;
    
    }

    
    public void makeSkeletonRadialSF(int Lin){
        int max = 1+(int)(Lin/2.0);
        radialVals = new ArrayList<Double>();
        for(int i = 0; i < (int)Lp;i++){for(int j = 0; j < (int)Lp;j++){
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
    
    public void takeFT(SimSystem lat){
        takeFT((LatticeMagInt)lat);
    }
    
    public void takeFT(LatticeMagInt lat){
            double dx = (L/Lp);
            for (int j = 0; j < Lp; j++){
                for (int i = 0; i < Lp; i++){
                        fftData[2*(i+j*(int)Lp)] = lat.getValue(i, j);
                        fftData[2*(i+j*(int)Lp)+1] = 0.0;
                }
            }
            fft.transform(fftData);
            fftData = fft.toWraparoundOrder(fftData);
            maxfft = 0;
            for (int i=0; i < Lp*Lp; i++){
                double re = fftData[2*i];
                double im = fftData[2*i+1];
                sFactor[i] = (re*re + im*im)/(L*L);
                if(sFactor[i] > maxfft){maxfft = sFactor[i];}
            }
            shiftSFactor();
    }
    
    public void takeFT(double[] data){
            double dx = (L/Lp);
            for (int i = 0; i<Lp*Lp; i++){
                    fftData[2*i] = data[i]*dx*dx;
                    fftData[2*i+1] = 0;
            }
            fft.transform(fftData);
            fftData = fft.toWraparoundOrder(fftData);
            for (int i=0; i < Lp*Lp; i++){
                    double re = fftData[2*i];
                    double im = fftData[2*i+1];
                    sFactor[i] = (re*re + im*im)/(L*L);
            }
            shiftSFactor();
    }

    
    public BufferedImage getVisImg(){
        vis.spinDraw(sFactor, maxfft);
        return vis.getImageOfVis();
    }
    public void takeFTUp(SimSystem lat){
        takeFTUp((LatticeMagInt)lat);
    }
    
    public void takeFTUp(LatticeMagInt lat){
            double dx = (L/Lp);
            int u; int v;
            for (int j = 0; j < Lp; j++){
                for (int i = 0; i < Lp; i++){
                        fftData[2*(i+j*(int)Lp)] = (lat.getValue(i, j) > 0) ? 1*dx*dx: 0;
                        fftData[2*(i+j*(int)Lp)+1] = 0.0;
                }
            }
            fft.transform(fftData);
            fftData = fft.toWraparoundOrder(fftData);
            maxfft = 0;
            for (int i=0; i < Lp*Lp; i++){
                double re = fftData[2*i];
                double im = fftData[2*i+1];
                sFactor[i] = (re*re + im*im)/(L*L);
                if(sFactor[i] > maxfft){maxfft = sFactor[i];}
            }
            shiftSFactor();
    }

    public double getX(int i, int j){
        if(geo == 6){
            return getXcoordTri(i,j);
        }else if(geo == 3){
            return getXcoordHoney(i,j);
        }else{
            return i;
        }
    }
    
    public double getY(int i, int j){
        if(geo == 6){
            return getYcoordTri(i,j);
        }else if(geo == 3){
            return getYcoordHoney(i,j);
        }else{
            return j;
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
    
    private double getXrCentered(int u,int v){
        int cenInd = (int)(L/2.0);
        double xcen = getX(cenInd,cenInd);
        return (getX(u,v)-xcen);
    }
    
    private double getYrCentered(int u,int v){
        int cenInd = (int)(L/2.0);
        double ycen = getY(cenInd,cenInd);
        return (getY(u,v)-ycen);
    }
    
    public void takeFTDilution(LatticeMagInt lat){
            double dx = (L/Lp);
            for (int j = 0; j < Lp; j++){
                for (int i = 0; i < Lp; i++){ 
                        fftData[2*(i+j*(int)L)] = (lat.getValue(i, j) == 0) ? 1*dx*dx: 0;
                        fftData[2*(i+j*(int)L)+1] = 0.0;
                }
            }
            fft.transform(fftData);
            fftData = fft.toWraparoundOrder(fftData);
            for (int i=0; i < Lp*Lp; i++){
                    double re = fftData[2*i];
                    double im = fftData[2*i+1];
                    sFactor[i] = (re*re + im*im)/(L*L);
            }
            shiftSFactor();
    }

    public void takeFTDown(LatticeMagInt lat){
            double dx = (L/Lp);
            for (int j = 0; j < Lp; j++){
                for (int i = 0; i < Lp; i++){ 
                        fftData[2*(i+j*(int)L)] = (lat.getValue(i, j) < 0) ? 1*dx*dx: 0;
                        fftData[2*(i+j*(int)L)+1] = 0.0;
                }
            }
            fft.transform(fftData);
            fftData = fft.toWraparoundOrder(fftData);
            for (int i=0; i < Lp*Lp; i++){
                    double re = fftData[2*i];
                    double im = fftData[2*i+1];
                    sFactor[i] = (re*re + im*im)/(L*L);
            }
            shiftSFactor();
    }

    public void shiftSFactor(){
            double [] temp = new double [Lp*Lp];
            //sFactor[0]=0;
            for (int i = 0; i<Lp*Lp; i++){
                    int x = i%Lp;
                    int y = i/Lp;
                    x += Lp/2; y += Lp/2;
                    x = x%Lp; y = y%Lp;
                    int j = Lp*((y+Lp)%Lp) + (x+Lp)%Lp;
                    temp[j] = sFactor[i];
            }
            for(int i = 0; i<Lp*Lp; i++)
                    sFactor[i] = temp[i];
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
        for(int i = 0; i < Lp;i++){for(int j = 0; j < Lp;j++){
            // Get angle relative to center of lattice
            double angle = atan((double)abs(getXrCentered(i,j))/(double)abs(getYrCentered(i,j)));        
            // Convert angle to range [0,2pi]
            angle = standardizeAngle(getXrCentered(i,j),getYrCentered(i,j), angle);
            // Convert angle to proper bin
            int ind = (int)floor(angle*circularBins/(2*PI));
            // Add contribution from angle
            circularSF[ind] += sFactor[(int)(i+j*L)];
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
    
    //Circular Average
    public void calcSFradial(){
        double maxRad = 0;double rad;
        radialSF = new double[radialBins];
        double[] rN = new double[radialBins];
        
        for(int i = 0; i < Lp;i++){for(int j = 0; j < Lp;j++){
            // Get distance from center
            rad= Math.sqrt(getXrCentered(i,j)*getXrCentered(i,j)+getYrCentered(i,j)*getYrCentered(i,j));
            // Convert distance into index in lattice using precomputed hashmap
            int ind = radialMap.get(rad);
            // Add contribution for this r
            radialSF[ind] += sFactor[(int)(i+j*L)];
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
    
    
    public void SpinSFUp(LatticeMagInt lat){
    	takeFTUp(lat);
    }
    public void SpinSFDown(LatticeMagInt lat){
    	takeFTDown(lat);
    }
    
    //sfactor[Lp/2*Lp+Lp/2]=0;
    public double getSF(int i){
        return sFactor[i];
    }

    public double[] getSFArr(){
        return sFactor;
    }
    
}
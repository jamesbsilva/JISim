package AnalysisAndVideoBackend;

/**
 * 
 *   @(#) StructureFactorPeaks
 */ 

import Backbone.Util.SystemAverager;

/**  An analysis class for structure factor data.
 *  <br>
 * 
 *  @param g - geometry  of lattice (6 - triangular / 3 honey comb / else square}
 * 
 * @author      James B. Silva <jbsilva @ bu.edu>                 
 * @since       2013-07    
 */
public class StructureFactorPeaks{
    private double[] maxPeaks;
    private int geo;private int q;
    private int L;private int R;
    private int goodPeaks = 0;
    private int peaks = 0;
    private double angle;
    private int drPeak = 4;
    private int center;
    private int pkDataN;
    private int[] peaksLoc;
    private double[] avgPeakData;
    private double[] peaksStd;
    private double[] peaksMean;
    private double[] phaseState;
    private int maxPeakIndex;
    private double backgroundMean = 0;
    private double backgroundStd = 0;
    private int backgroundSamples = 500;
    private double[] backgroundData = new double[backgroundSamples];
    private boolean stripesPhase = false;
    private boolean singleStripesPhase = false;
    private boolean clumpsPhase = false; private boolean uniformPhase = false;
    private SystemAverager peaksAvg = new SystemAverager();
         
    public StructureFactorPeaks(int g, int lin, int rin){
        geo = g;L = lin;R = rin;
        angle = 60*Math.PI/180;
        peaksAvg.addObservable();peaksAvg.addObservable();
        peaksAvg.addObservable();peaksAvg.addObservable();
        phaseState = new double[7];
        if(g == 6 || geo == 3){
            maxPeaks = new double[6];peaksLoc = new int[6];
            peaksStd = new double[6];peaksMean = new double[6];
            center = (int)((L*(L+1))/2.0);
            peaksAvg.addObservable();peaksAvg.addObservable();
            if(g == 6){setQTri();}
        }else{
            peaksLoc = new int[4];maxPeaks = new double[4];
            peaksStd = new double[4];peaksMean = new double[4];          
        }
        pkDataN = 3*maxPeaks.length+3;
        avgPeakData = new double[pkDataN+getSizePeakData()];
        //System.out.println("SF-Peaks | base | "+pkDataN+"     rest| "+getSizePeakData());
    }
    public int getNumberOfMaxPeaks(){
        return goodPeaks;
    }
    public int getNumberPeaksOverNoise(){
        return peaks;
    }
    public boolean isSingleStripes(){
        return singleStripesPhase;
    }
    public boolean isStripes(){
        return stripesPhase;
    }
    public boolean isClumps(){
        return clumpsPhase;
    }
    public boolean isUniform(){
        return uniformPhase;
    }
    
    public double[] setData(float[] din){
        double maxVal = 0;
        maxPeaks[0] = getMax( din,getPossibleLeftPeak(din), 0 );
        maxPeaks[1] = getMax( din,getPossibleRightPeak(din),1 );
        if(geo == 3 || geo  == 6){
            maxPeaks[2] = getMax( din,getPossibleULeftPeak(din), 2 );
            maxPeaks[3] = getMax( din,getPossibleURightPeak(din),3 );
            maxPeaks[4] = getMax( din,getPossibleDLeftPeak(din), 4 );
            maxPeaks[5] = getMax( din,getPossibleDRightPeak(din),5 );
        }else if(geo == 4){
            maxPeaks[2] = getMax( din,getPossibleUpPeak(din), 2 );
            maxPeaks[3] = getMax( din,getPossibleDownPeak(din),3 );
        }
        // determine maxPeak
        for(int u = 0; u < maxPeaks.length;u++){
            if( maxPeaks[u] > maxVal ){
                maxPeakIndex = u; maxVal = maxPeaks[u];
            }
        }
        updatePeakAvg(maxPeaks); updateBackground(din); setPhase(maxPeaks);
        //System.out.println("Qmeasured: "+((peaksLoc[1]-peaksLoc[0])/2.0));
        return maxPeaks;
    }
    
    
    private void setPhase(double[] maxPeak){
        int goodPeak = 0;int peaksOverNoise = 0;
        double delTolerance = 1.30;
        double absMax = maxPeak[maxPeakIndex];double absStd = peaksStd[maxPeakIndex];
        for( int u = 0; u < maxPeak.length;u++ ){
            //System.out.println("Max Peak | "+maxPeak[u]+"  absStd | "+absStd+"   backgroundMean | "+backgroundMean+"      backStd | "+backgroundStd);
            if( (maxPeak[u] > (backgroundMean+10*delTolerance*backgroundStd)) ){
                peaksOverNoise++;
                if( (maxPeak[u] > (absMax-delTolerance*absStd)) ){goodPeak += 1;}
            }
        }
        if( goodPeak == 2 ){
            stripesPhase = true;
            if( peaksOverNoise == 2 ){singleStripesPhase = true;}
        }else{
            stripesPhase = false;
        }
        clumpsPhase = (goodPeak == maxPeaks.length) ? true : false;
        uniformPhase = (peaksOverNoise == 0) ? true : false ;
        //if(clumpsPhase){System.out.println("THIS IS CLUMPS");}
        //if(stripesPhase){System.out.println("THIS IS STRIPES");}
        //if(goodPeak < 1){System.out.println("THIS IS RANDOM");}
        phaseState[3] = goodPeak;goodPeaks = goodPeak;peaks = peaksOverNoise;
    }
    
    public double[] getPeakAvgData(){
        return avgPeakData;
    }

    public void updateBackground(float[]din){
        for(int u = 0; u < backgroundSamples;u++){
            backgroundData[u] = din[(int)(Math.random()*din.length)];
        }
        backgroundMean = getMean(backgroundData);
        backgroundStd = getStdDev(backgroundData,backgroundMean);
    }
    
    public int[] getPeakLoc(){
        return peaksLoc;
    }
    
    public double[] getBackgroundData(){
        double[] temp  = new double[2];
        temp[0] = backgroundMean;
        temp[1] = backgroundStd;
        return temp;
    }
    
    public double[] getPhaseDataOnly(){
        return phaseState;
    }
    public int getSizePeakData(){
        int sizeIn = 0;
        if(geo == 4) {
            sizeIn += 2*((1+2*drPeak)*(1+2*drPeak))+2*(1+2*drPeak);
        }else if(geo == 6 || geo == 3) {
            sizeIn += 6*((1+2*drPeak)*(1+2*drPeak));
        }
        return sizeIn;
    }
    public double getStdDev(double[] din, double mean){
        double sum = 0;
        for(int u = 0; u < din.length; u++){
            sum += Math.pow(din[u]-mean , 2.0);
        }
        sum /= din.length; sum = Math.sqrt(sum);
        return sum;
    }
    public double getStdDev(double[] din){
        double sum = 0;double mean = getMean(din);
        for(int u = 0; u < din.length; u++){
            sum += Math.pow(din[u]-mean , 2.0);
        }
        sum /= (double)din.length;sum = Math.sqrt(sum);
        return sum;
    }   
    public double getMean(double[] din){
        double sum = 0;
        for(int u = 0; u < din.length; u++){
            sum += din[u];
        }
        sum /= (double)din.length;
        return sum;
    }
    public double getMax(float[] din, int[] locs , int peak){
        double max = 0;
        for(int u = 0; u < locs.length;u++){
            if(locs[u] < din.length){
                if(din[locs[u]] >max){
                    max = din[locs[u]];
                    peaksLoc[peak] = locs[u];
                }
            }
        }
        return max;
    }
    public void setQTri(){
        q = (int) ((L*0.86757-2.2604)/(0.2441+R));
        if(R == 1){
            q = (int)(q/2.0);
        }
        if(q < 10){drPeak = drPeak/2;}if(q < 3){drPeak = 1;}
    }
    public void updatePeakAvg(double[] din){
        int chunk = 3;
        for(int u = 0; u < din.length;u++){
            if(u == 0){
                peaksAvg.updateObs(0, din[u], true);
                avgPeakData[chunk*u+1] = peaksAvg.getStdDevObs(0);
            }else if(u == 1){
                peaksAvg.updateObs(1, din[u], true);
                avgPeakData[chunk*u+1] = peaksAvg.getStdDevObs(1);
            }else if(u == 2){
                peaksAvg.updateObs(2, din[u], true);
                avgPeakData[chunk*u+1] = peaksAvg.getStdDevObs(2);
            }else if(u == 3){
                peaksAvg.updateObs(3, din[u], true);
                avgPeakData[chunk*u+1] = peaksAvg.getStdDevObs(3);
            }else if(u == 4){
                peaksAvg.updateObs(4, din[u], true);
                avgPeakData[chunk*u+1] = peaksAvg.getStdDevObs(4);
            }else if(u == 5){
                peaksAvg.updateObs(5, din[u], true);
                avgPeakData[chunk*u+1] = peaksAvg.getStdDevObs(5);
            }
            avgPeakData[chunk*u] = din[u];
            avgPeakData[chunk*u+2] = peaksLoc[u];
        }        
        avgPeakData[3*maxPeaks.length] = maxPeakIndex;
        avgPeakData[3*maxPeaks.length+1] = (stripesPhase) ? 1.0 : 0.0 ;
        avgPeakData[3*maxPeaks.length+2] = (clumpsPhase) ? 1.0 : 0.0 ;
        phaseState[0] = (stripesPhase) ? 1.0 : 0.0 ;
        phaseState[1] = (clumpsPhase) ? 1.0 : 0.0 ;
        phaseState[2] = (singleStripesPhase) ? 1.0 : 0.0 ;
    }
    public int[] getPossibleLeftPeak(float[] din){
        int i = getPeakLocationLeft(q);
        int[] peaks = new int[(1+2*drPeak)*(1+2*drPeak)];
        double[] peaksVal = new double[(1+2*drPeak)*(1+2*drPeak)];
        double sum = 0.0;double val;
        int x = getXcoord(i); int y = getYcoord(i);
        int ind = pkDataN+0;
        int m; int n;
        for(int u = 0; u < (2*drPeak+1);u++){for(int v = 0; v < (2*drPeak+1); v++){
            m = (x-drPeak+u+L)%L; n = (y-drPeak+v+L)%L;
            peaks[u+(2*drPeak+1)*v] = (m+n*L);
            val = din[(m+n*L)]; peaksVal[u+(2*drPeak+1)*v] = val;
            avgPeakData[u+v*(2*drPeak+1)+ind] = val;
            sum += val;
        }}
        peaksStd[0] = getStdDev(peaksVal, sum/(double)peaksVal.length);
        peaksMean[0] = sum/(double)peaksVal.length;
        return peaks;
    }
    public int[] getPossibleUpPeak(float[] din){
        int i = getPeakLocationUp(q);
        int[] peaks = new int[(1+2*drPeak)];
        double[] peaksVal = new double[(1+2*drPeak)];
        double sum = 0.0;double val;
        int x = getXcoord(i);
        int y = getYcoord(i);
        int ind = pkDataN+2*((2*drPeak+1)*(2*drPeak+1));
        int m;
        for(int u = 0; u < (2*drPeak+1);u++){
            m = (y-drPeak+u+L)%L;val = din[(x+m*L)];
            peaks[u] = (x+m*L); peaksVal[u] = val;
            avgPeakData[u+ind] = (x+m*L);
            sum += val;
        }
        peaksStd[2] = getStdDev( peaksVal, sum/(double)peaksVal.length );
        peaksMean[2] = sum/(double)peaksVal.length;
        return peaks;
    }
    public int[] getPossibleDownPeak(float[] din){
        int i = getPeakLocationDown(q);
        int[] peaks = new int[(1+2*drPeak)]; double[] peaksVal = new double[(1+2*drPeak)];
        double sum = 0.0; double val;
        int x = getXcoord(i); int y = getYcoord(i); int m;
        int ind = pkDataN+3*((2*drPeak+1)*(2*drPeak+1));
        for(int u = 0; u < (2*drPeak+1);u++){
            m = (y-drPeak+u+L)%L;
            val = din[(x+m*L)];
            peaks[u] = (x+m*L); peaksVal[u] = val;
            avgPeakData[u+ind] = (x+m*L);
            sum += val;
        }
        peaksStd[3] = getStdDev(peaksVal, sum/(double)peaksVal.length);
        peaksMean[3] = sum/(double)peaksVal.length;
        return peaks;
    }
    public int[] getPossibleRightPeak(float[] din){
        int i = getPeakLocationRight(q);
        int[] peaks = new int[(1+2*drPeak)*(1+2*drPeak)];
        double[] peaksVal = new double[(1+2*drPeak)*(1+2*drPeak)];
        double sum = 0.0;double val;
        int x = getXcoord(i);
        int y = getYcoord(i);
        int ind = pkDataN+((2*drPeak+1)*(2*drPeak+1));
        int m; int n;
        for(int u = 0; u < (2*drPeak+1);u++){for(int v = 0; v < (2*drPeak+1); v++){
            m = (x-drPeak+u+L)%L; n = (y-drPeak+v+L)%L;
            val = din[(m+n*L)];
            peaks[u+(2*drPeak+1)*v] = (m+n*L); peaksVal[u+(2*drPeak+1)*v] = val;
            avgPeakData[u+v*(2*drPeak+1)+ind] = (m+n*L);
            sum += val;
        }}
        peaksStd[1] = getStdDev(peaksVal, sum/(double)peaksVal.length);
        peaksMean[1] = sum/(double)peaksVal.length;
        return peaks;
    }
    public int[] getPossibleDLeftPeak(float[] din){
        int i = getPeakLocationDL(q);
        int[] peaks = new int[(1+2*drPeak)*(1+2*drPeak)];
        double[] peaksVal = new double[(1+2*drPeak)*(1+2*drPeak)];
        double sum = 0.0;double val;
        int x = getXcoord(i); int y = getYcoord(i);
        int ind = pkDataN+2*((2*drPeak+1)*(2*drPeak+1));
        int m; int n;
        for(int u = 0; u < (2*drPeak+1);u++){for(int v = 0; v < (2*drPeak+1); v++){
            m = (x-drPeak+u+L)%L; n = (y-drPeak+v+L)%L;
            val = din[(m+n*L)];
            peaks[u+(2*drPeak+1)*v] = (m+n*L); peaksVal[u+(2*drPeak+1)*v] = val;
            avgPeakData[u+v*(2*drPeak+1)+ind] = (m+n*L);
            sum += val;
        }}
        peaksStd[4] = getStdDev(peaksVal, sum/(double)peaksVal.length);
        peaksMean[4] = sum/(double)peaksVal.length;
        return peaks;
    }
    public int[] getPossibleDRightPeak(float[] din){
        int i = getPeakLocationDR(q);
        int[] peaks = new int[(1+2*drPeak)*(1+2*drPeak)];
        double[] peaksVal = new double[(1+2*drPeak)*(1+2*drPeak)];
        double sum = 0.0;double val;
        int x = getXcoord(i); int y = getYcoord(i);
        int ind = pkDataN+3*((2*drPeak+1)*(2*drPeak+1));
        int m; int n;
        for(int u = 0; u < (2*drPeak+1);u++){for(int v = 0; v < (2*drPeak+1); v++){
            m = (x-drPeak+u+L)%L; n = (y-drPeak+v+L)%L;
            val = din[(m+n*L)];
            peaks[u+(2*drPeak+1)*v] = (m+n*L); peaksVal[u+(2*drPeak+1)*v] = val;
            avgPeakData[u+v*(2*drPeak+1)+ind] = (m+n*L);
            sum += val;
        }}
        peaksStd[5] = getStdDev(peaksVal, sum/(double)peaksVal.length);
        peaksMean[5] = sum/(double)peaksVal.length;
        return peaks;
    }
    public int[] getPossibleURightPeak(float[] din){
        int i = getPeakLocationUR(q);
        int[] peaks = new int[(1+2*drPeak)*(1+2*drPeak)];
        double[] peaksVal = new double[(1+2*drPeak)*(1+2*drPeak)];
        double sum = 0.0; double val;
        int x = getXcoord(i); int y = getYcoord(i);
        int ind = pkDataN+4*((2*drPeak+1)*(2*drPeak+1));
        int m; int n;
        for(int u = 0; u < (2*drPeak+1);u++){for(int v = 0; v < (2*drPeak+1); v++){
            m = (x-drPeak+u+L)%L; n = (y-drPeak+v+L)%L;
            val = din[(m+n*L)];
            peaks[u+(2*drPeak+1)*v] = (m+n*L); peaksVal[u+(2*drPeak+1)*v] = val;
            avgPeakData[u+v*(2*drPeak+1)+ind] = (m+n*L);
            sum += val;
        }}
        peaksStd[3] = getStdDev(peaksVal, sum/(double)peaksVal.length);
        peaksMean[3] = sum/(double)peaksVal.length;
        return peaks;
    }
    
    public int[] getPossibleULeftPeak(float[] din){
        int i = getPeakLocationUL(q);
        int[] peaks = new int[(1+2*drPeak)*(1+2*drPeak)];
        double[] peaksVal = new double[(1+2*drPeak)*(1+2*drPeak)];
        double sum = 0.0;double val;
        int x = getXcoord(i); int y = getYcoord(i);
        int ind = pkDataN+5*((2*drPeak+1)*(2*drPeak+1));
        int m; int n;
        for(int u = 0; u < (2*drPeak+1);u++){for(int v = 0; v < (2*drPeak+1); v++){
            m = (x-drPeak+u+L)%L; n = (y-drPeak+v+L)%L;
            val = din[(m+n*L)];
            peaks[u+(2*drPeak+1)*v] = (m+n*L); peaksVal[u+(2*drPeak+1)*v] = val;
            avgPeakData[u+v*(2*drPeak+1)+ind] = (m+n*L);
            sum += val;
        }}
        peaksStd[2] = getStdDev(peaksVal, sum/(double)peaksVal.length);
        peaksMean[2] = sum/(double)peaksVal.length;
        return peaks;
    }
    public int getXcoord(int i){
        return (i%L);
    }
    public int getYcoord(int i){
        return (((int)((double)i/(double)L))%L);
    }    
    public int getPeakLocationUp(int q){
        // u peak
        int x = getXcoord(center); int y = getYcoord(center);
        int peakU = (x+(y+L-q)%L*L);
        return peakU;
    }
    public int getPeakLocationDown(int q){
        // u peak
        int x = getXcoord(center); int y = getYcoord(center);
        int peakD = (x+(y+L+q)%L*L);
        return peakD;
    }
    public int getPeakLocationLeft(int q){
        // l peak
        int peakL = center - q;
        return peakL;
    }
    public int getPeakLocationRight(int q){
        // r peak
        int peakR = center + q;
        return peakR;
    }
    public int getPeakLocationDL(int q){
        int y = (int)(q*Math.sin(angle)); int x = (int)(q*Math.cos(angle));
        int peakDL = center + (int)(-x+y*L);
        return peakDL;
    }
    public int getPeakLocationUR(int q){
        int y = (int)(q*Math.sin(angle)); int x = (int)(q*Math.cos(angle));
        int peakUR = center + (int)(x-y*L);
        return peakUR;
    }
    public int getPeakLocationUL(int q){
        int y = (int)(q*Math.sin(angle)); int x = (int)(q*Math.cos(angle));
        int peakUL = center + (int)(-x-y*L)-1;
        return peakUL;
    }
    public int getPeakLocationDR(int q){
        // DR
        int y = (int)(q*Math.sin(angle)); int x = (int)(q*Math.cos(angle));
        int peakDR = center + (int)(x+y*L)+1;
        return peakDR;
    }
}
package Backbone.System;

/**
* 
*   @(#) DblValThreshRange
*/  
import Backbone.Visualization.BarsDouble1D;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Random;

/**
*      DblValThreshRange is a rangeDbl class that is composed of a fully connected system
*  which is represented by a single value. 
* 
* 
* <br>
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2013-07    
*/
public final class DblValThreshRange implements LatticeDbl{
    private ArrayList<Double> rangeDbl;
    private ArrayList<Double> rangeDblThresh;
    private double baseThresh = 1;
    private int initVal = 100;
    private int N;
    private int R;
    private int L;
    private int Nfixed;
    private int magnetization;
    private int currentTime=0;
    private int instId;
    private Random Ran = new Random();
    private BarsDouble1D vis;
    private boolean initialized=false;
    private boolean useHeter;
    private boolean makingVideo=false;
    private double maxValue = 1;
    private double asperPer = 0.0;
    private double asperRatio = 1.0;
    
    public DblValThreshRange(int num, double s){
        N = num;
        makeLatticeBasic();
        vis = new  BarsDouble1D(N); 
        vis.initializeImg();
        initialize(s,(int)(0.1*N));
    }
    public DblValThreshRange(int num, double s, int inst){
        N = num;
        makeLatticeBasic();
        instId = inst;
        vis = new  BarsDouble1D(N);
        vis.initializeImg();
        initialize(s,(int)(0.1*N));
    }
    public DblValThreshRange(int num,int loc, double s, int inst){
        N = num;
        makeLatticeBasic();
        instId = inst;
        vis = new  BarsDouble1D(N); 
        vis.initializeImg();
        initialize(s,loc);
    }   
    
    public DblValThreshRange(int num,int loc, double s, int inst,double  bthresh){
        N = num;
        baseThresh = bthresh;
        makeLatticeBasic();
        instId = inst;
        vis = new  BarsDouble1D(N); 
        vis.initializeImg();
        initialize(s,loc);
    }
    
    public void initialize(double s){
        maxValue = baseThresh;
        makeLatticeBasic();
        initialize(s,(int)(0.1*N));
        initialized = true;
    }
    
    public void makeLatticeBasic(){
        rangeDbl = new ArrayList<Double>();
        rangeDblThresh = new ArrayList<Double>();
        for(int u = 0; u < (10*N);u++){
            rangeDbl.add(0.0);
            rangeDblThresh.add(baseThresh);
        }
        N = 10*N;
    }
    
    public void reshapeLat(int newMax){
        //System.out.println("Reshape from | "+rangeDbl.size()+"    to  | "+newMax);
        if( newMax < rangeDbl.size() ){
            return;
        }
        int diff = newMax-rangeDbl.size()+100;
        for(int u = 0; u < (2*diff);u++){
            rangeDbl.add(0.0);
            if(Ran.nextDouble() < asperPer){
                rangeDblThresh.add(baseThresh*asperRatio);
            }else{
                rangeDblThresh.add(baseThresh);
            }
        }
        N = rangeDbl.size();
        vis.spinDraw(rangeDbl,maxValue);
    }
    
    /**
    *         setValue which updates the spin value in the rangeDbl and the magnetization
    *       with update of time in image layer .
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    *  @param k - k coordinate
    *  @param s - new spin value
    *  @param t - time 
    */ 
    public void setValue(int u, double s, int t){
        currentTime = t;
        if(u >= rangeDbl.size()){
            reshapeLat(u);
        }
        rangeDbl.set(u , s);
        //System.out.println("N: "+N+"    Nup: "+Nup+"    Ndown: "+Ndown);
        if(initialized && makingVideo){
            if(s > maxValue){maxValue = s;}
            updateImg(u,s,t);
        }
    }
    
    public void setFixedRandom(int nfix, double largeVal){
        asperPer = (double)nfix/rangeDbl.size();
        asperRatio = largeVal/baseThresh;
        
        for(int u = 0;u < nfix;u++){
            rangeDblThresh.set((int)(Ran.nextDouble()*rangeDblThresh.size()), largeVal);
        }
        maxValue = (maxValue < largeVal ) ? largeVal : maxValue;
    }
    
    public void setValue(int i,int j, int k,double s, int t){
        setValue(i+j*L,s,t);
    }
    /**
    *         getValue which gets the spin value in the rangeDbl. 3d case
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    *  @param k - k coordinate
    */ 
    public double getValue(int i,int j, int k){
        return rangeDbl.get(i+j*L+k*L*L);
    }
    /**
    * getValue which gets the spin value in the rangeDbl.
    *
    * @param i - coordinate
    */
    @Override
    public double getValue(int i) {
        return rangeDbl.get(i);
    }
    
    public double getThresh(int i){
        return (rangeDblThresh.get(i));
    }
    
    public double getDiffValThresh(int i){
        return (rangeDblThresh.get(i)-getValue(i));
    }
    
    /**
    *         getLength which gets the length of the rangeDbl
    * 
    */ 
    public int getLength(){return rangeDbl.size();}
    /**
    *         getDimension gives the dimensionality of the rangeDbl
    */ 
    public int getDimension(){return 1;}
    /**
    *         getGeo gives the geometry of the rangeDbl
    */ 
    public int getGeo(){return 1;}
    /**
    *         getN gives the size of the rangeDbl
    */ 
    public int getN(){return rangeDbl.size();}
    /**
    *         getRange gives the interaction range of the rangeDbl
    */ 
    public int getRange(){return rangeDbl.size();}
    /**
    *         getMagnetization gives the magnetization of the rangeDbl
    */ 
    public int getMagnetization(){return magnetization;}
    /**
    *         setFixedLatticeValues should set all fixed values in the rangeDbl.
    */ 
    public void setFixedLatticeValues(){N = (1-(Nfixed/N));}
    
    /**
    *         initialize sets the rangeDbl values to their initial values which is given as input.
    *
    *
    *  @param s - Initial rangeDbl value; -2 for random values
    */ 
    public void initialize(double s, int loc){
        if(s != 0.0){
            // draw
            maxValue = baseThresh;
            initialized = true;
            return;
        }
        if(s < 1){s = s*baseThresh;}
        if(loc > 0 && loc < N){
            for(int u = 0; u < N; u++){
                rangeDbl.set(u, (u == loc) ? s : s*Ran.nextDouble()*0.9);
            }
        }
        if(loc < 0){
            for(int u = 0; u < N;u++){
                rangeDbl.set(u, baseThresh*Ran.nextDouble());
            }
        }
        // draw
        maxValue = baseThresh;
        vis.spinDraw(rangeDbl,maxValue);
        initialized = true;
    }
    /**
    *         isThisFixed should return true if the rangeDbl coordinate is fixed.
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    *  @param k - k coordinate 
    */ 
    public boolean isThisFixed(int i,int j,int k){
        return  (rangeDblThresh.get(i+j*L) > baseThresh) ? true: false ;
    }
    
    /**
    *         isThisFixed should return true if the rangeDbl coordinate is fixed.
    * 
    *  @param i - i coordinate
    */ 
    public boolean isThisFixed(int i){
        if(rangeDblThresh.get(i) > baseThresh){
            //System.out.println("Base | "+baseThresh+"    latThresh | "+rangeDblThresh.get(i));
            return true;
        }else{
            return false;
        }
    }
    
    /**
    *         getBaseThresh returns the base threshold.
    * 
    */ 
    public double getBaseThresh(){
        return baseThresh;
    }
    
    /**
    *         getNeighSum should return the sum of the spins within the interaction
    *   range centered by the coordinates given. 3d case
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    *  @param k - k coordinate 
    */ 
    public int getNeighSum(int i,int j,int k){
        return getMagnetization();
    }
    /**
    *         getNeighSum should return the sum of the spins within the interaction
    *   range centered by the coordinates given. 2d case
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    */ 
    public int getNeighSum(int i,int j){
        return getMagnetization();
    }
    /**
    *         getNinRnage should return the amount of spins within the interaction 
    *   range. Useful for metropolis algorithm in long range case.
    */ 
    public int getNinRange(){return R;}
    /**
    *         setInitialConfig should set all spins in the rangeDbl to the values
    *   given by a file in the Config directory.
    *   
    *   @param t - time of config to set rangeDbl to
    *   @param run - run of time to set rangeDbl to
    *   @param  post - postfix or filename of rangeDbl file
    */
    public void setInitialConfig(int t,int run,String post){}
    /**
    *         getInstId returns the identification number of the instance of this
    *  program running. Useful for allowing each instance a directory to work in.
    */ 
    public int getInstId(){return instId;}
    /**
    *         spinDraw completely redraws the image layer based on the current 
    *   configuration of the rangeDbl.
    */
    private void spinDraw(){
        vis.spinDraw(rangeDbl, maxValue);
    }
    
    /**
    *         getSystemImg should return an image of the rangeDbl
    */ 
    public BufferedImage getSystemImg(){
        return vis.getImageOfVis();
    }    
    /**
    *         updateImg just updates the image layer instead of redrawing the layer.
    *   Change time text of image
    * 
    *   @param i - i coordinate
    *   @param j - j coordinate
    *   @param spin  - new spin value
    *   @param t  - new time in image
    * 
    */
    private void updateImg(int loc,double val,int t){
        vis.updateImg(loc%N, val, t, maxValue);
    }
    // No update of time text
    /**
    *         updateImg just updates the image layer instead of redrawing the layer.
    *   Does not change time text of image
    * 
    *   @param i - i coordinate
    *   @param j - j coordinate
    *   @param spin  - new spin value
    * 
    */
    private void updateImg(int i, double spin){
        updateImg(i,spin,0);
    }

    /**
    * makeVideo sets flag notifying images returned should be updated
    */
    @Override
    public void makeVideo(){makingVideo=true;spinDraw();}
    /**
    *         getMagStaggered gives the regular magnetization for single value rangeDbl
    */ 
    public int getMagStaggered(){
        return magnetization;
    }    

    @Override
    public int getFixSpinVal() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public int getNFixed() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public double getValue(int i, int j) {
        throw new UnsupportedOperationException("Not supported for 1D object."); //To change body of generated methods, choose Tools | Templates.
    }
}



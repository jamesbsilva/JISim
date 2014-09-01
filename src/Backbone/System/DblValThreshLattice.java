package Backbone.System;

/**
* 
*   @(#) DblValThreshLattice
*/  
import Backbone.Visualization.BarsDouble1D;
import Backbone.Visualization.SquareDoubleLattice2D;
import Backbone.Visualization.SysVisualizationDbl;
import Backbone.Visualization.SysVisualizationDbl2D;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.util.Random;

/**
*      DblValThreshLattice is a lattice class that is composed of a double value threshold
*  which is represented by a single value. 
* 
* 
* <br>
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2011-09    
*/
public final class DblValThreshLattice implements LatticeDbl{
    private double[] lattice;
    private double[] latticeThresh;
    private double baseThresh = 1.0;
    private int N;
    private int R;
    private int L;
    private int Nfixed;
    private int magnetization;
    private int currentTime=0;
    private int instId;
    private Random Ran = new Random();
    private SysVisualizationDbl2D vis2D;
    private SysVisualizationDbl vis1D;
    private int ImageTextOffset = 90;
    private Graphics2D g;
    private int scale = 5;
    private int maxImgSize = 1600;
    private int minImgSize = 600;
    private boolean initialized=false;
    private int D = 2;
    private boolean useHeter;
    private boolean makingVideo=true;
    private double maxValue = 1;
    
    public DblValThreshLattice(int num, double s){
        N = num;
        lattice = new double[N];
        latticeThresh = new double[N];
        initVis();
        initialize(s,(int)(0.1*N));
    }
    
    public DblValThreshLattice(int num, double s, int inst){
        N = num;
        lattice = new double[N];
        latticeThresh = new double[N];
        instId = inst;
        initVis();
        initialize(s,(int)(0.1*N));
    }
    public DblValThreshLattice(int num,int loc, double s, int inst){
        N = num;
        lattice = new double[N];
        latticeThresh = new double[N];
        instId = inst;
        initVis();
        initialize(s,loc);
    }   
    
    public DblValThreshLattice(int num,int loc, double s, int inst,double  bthresh){
        this(num,loc,s,inst,bthresh,2);
    }
    
    public DblValThreshLattice(int num,int loc, double s, int inst,double  bthresh,int dim){
        D = dim;
        N = num;
        lattice = new double[N];
        baseThresh = bthresh;
        instId = inst;
        initVis();
        initialize(s,loc);
    }
    
    private void initVis(){
        if(D == 2){
            vis2D = new  SquareDoubleLattice2D((int)Math.sqrt(N));
            vis2D.initializeImg();
        }else if(D == 1){
            vis1D =  new BarsDouble1D(N);
            vis1D.initializeImg();
        }
    }
    
    private void visDraw(double[] lat, double bthr){
        if(D == 2){
            vis2D.spinDraw(lat, bthr);
        }else if(D == 1){
            vis1D.spinDraw(lat, bthr);
        }
    }
    
    public void initialize(double s){
        maxValue = baseThresh;
        initialize(s,(int)(0.1*N));
        visDraw(lattice, baseThresh);
        initialized = true;
    }
    
    /**
    *         setValue which updates the spin value in the lattice and the magnetization
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
        lattice[u] = s;
        //System.out.println("N: "+N+"    Nup: "+Nup+"    Ndown: "+Ndown);
        if(initialized && makingVideo){
            updateImg(u,s,t);
        }
    }
    
    public void setLatThresh(double[] lin){latticeThresh = lin;}

    public void setFixed(double[] fix){
        latticeThresh = fix;
        double max = -1.0;
        for(int u = 0;u < fix.length;u++){
            if(fix[u] > max){
                max =  fix[u];
            }
        }
        maxValue = max;
        visDraw(lattice,max);
    }
    
    public void setValue(int i,int j, int k,double s, int t){
        setValue(i+j*L,s,t);
    }
    /**
    *         getValue which gets the spin value in the lattice. 3d case
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    *  @param k - k coordinate
    */ 
    public double getValue(int i,int j, int k){
        return lattice[i+j*L+k*L*L];
    }
    /**
    * getValue which gets the spin value in the lattice.
    *
    * @param i - coordinate
    */
    @Override
    public double getValue(int i) {
        return lattice[i];
    }
    /**
    *         getValue which gets the spin value in the lattice. 2d case
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    */ 
    public double getValue(int i,int j){
        return lattice[i+j*L];
    }
    
    public double getThresh(int i){
        return (latticeThresh[i]);
    }
    
    public double getDiffValThresh(int i,int j){
        return (latticeThresh[i+j*L]-getValue(i,j));
    }
    
    public double getDiffValThresh(int i){
        return (latticeThresh[i]-getValue(i,0));
    }
    
    /**
    *         getLength which gets the length of the lattice
    * 
    */ 
    public int getLength(){return N;}
    /**
    *         getDimension gives the dimensionality of the lattice
    */ 
    public int getDimension(){return 1;}
    /**
    *         getGeo gives the geometry of the lattice
    */ 
    public int getGeo(){return 1;}
    /**
    *         getN gives the size of the lattice
    */ 
    public int getN(){return N;}
    /**
    *         getRange gives the interaction range of the lattice
    */ 
    public int getRange(){return N;}
    /**
    *         getMagnetization gives the magnetization of the lattice
    */ 
    public int getMagnetization(){return magnetization;}
    /**
    *         setFixedLatticeValues should set all fixed values in the lattice.
    */ 
    public void setFixedLatticeValues(){N = (1-(Nfixed/N));}
    
    /**
    *         initialize sets the lattice values to their initial values which is given as input.
    *
    *
    *  @param s - Initial lattice value; -2 for random values
    *
    */ 
    public void initialize(double s, int loc){
        if(s < 1){s = s*baseThresh;}
        if(loc > 0 && loc < N){
            for(int u = 0; u < N; u++){
                lattice[u] = (u == loc) ? s : s*Ran.nextDouble()*0.9;
            }
        }
        if(loc < 0){
            for(int u = 0; u < N;u++){
                lattice[u] = baseThresh*Ran.nextDouble();
            }
        }
        initialized = true;
    }
    /**
    *         isThisFixed should return true if the lattice coordinate is fixed.
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    *  @param k - k coordinate 
    */ 
    public boolean isThisFixed(int i,int j,int k){
        return  (latticeThresh[i+j*L] > baseThresh) ? true: false ;
    }
    
    /**
    *         isThisFixed should return true if the lattice coordinate is fixed.
    * 
    *  @param i - i coordinate
    */ 
    public boolean isThisFixed(int i){
        if(latticeThresh[i] > baseThresh){
            //System.out.println("Base | "+baseThresh+"    latThresh | "+latticeThresh[i]);
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
    *         setInitialConfig should set all spins in the lattice to the values
    *   given by a file in the Config directory.
    *   
    *   @param t - time of config to set lattice to
    *   @param run - run of time to set lattice to
    *   @param  post - postfix or filename of lattice file
    */
    public void setInitialConfig(int t,int run,String post){}
    /**
    *         getInstId returns the identification number of the instance of this
    *  program running. Useful for allowing each instance a directory to work in.
    */ 
    public int getInstId(){return instId;}
    /**
    *         spinDraw completely redraws the image layer based on the current 
    *   configuration of the lattice.
    */
    private void spinDraw(){
        double val;
        for(int u = 0;u<N;u++){
            val = getValue(u);
            updateImg(u,val,currentTime);
        }
    }
    /**
    *         getSystemImg should return an image of the lattice
    */ 
    public BufferedImage getSystemImg(){
        return (D == 1) ? vis1D.getImageOfVis() : vis2D.getImageOfVis();
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
        if(D == 2){
            int L = (int) Math.sqrt(N);
            vis2D.updateImg(loc%L, (int)(((double)loc/(double)L)%L) , val, t, maxValue);
        }else if (D == 1){
            vis1D.updateImg(loc%N, val, t, maxValue);
        }
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
    private void updateImg(int i, int j, double spin){
        updateImg((i+j*L),spin,0);
    }

    /**
    * makeVideo sets flag notifying images returned should be updated
    */
    @Override
    public void makeVideo(){makingVideo=true;spinDraw();}
    /**
    *         getMagStaggered gives the regular magnetization for single value lattice
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
}



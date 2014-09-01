package Backbone.System;

/**
* 
*   @(#) SingleValueLattice
*/  
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.util.Random;

/**
*      SingleValueLattice is a lattice class that is composed of a fully connected system
*  which is represented by a single value. 
* 
* 
* <br>
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2011-09    
*/
public final class SingleValueLattice implements LatticeMagInt{
    private int Nup;
    private int Ndown;
    private int NfixUp = 0;
    private int NfixDown = 0;
    private int N;
    private int Nfixed;
    private int magnetization;
    private int currentTime;
    private int instId;
    private Random Ran = new Random();
    private BufferedImage latticeImg;
    private int ImageTextOffset = 90;
    private Graphics2D g;
    private int scale = 5;
    private int maxImgSize = 1600;
    private int minImgSize = 600;
    private boolean initialized=false;
    private boolean useHeter;
    private boolean makingVideo=false;
    
    public SingleValueLattice(int num, int s){
        N = num;
        initialize(s);
    }
    public SingleValueLattice(int num, int s, int inst){
        N = num;
        instId = inst;
        initialize(s);
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
    public void setValue(int s, int t){
        currentTime = t;
        if(s == (+1)){
            Nup++;Ndown--;
        }else{
            Ndown++;Nup--;
        }
        magnetization = Nup-Ndown;//+NfixUp-NfixDown;
        //System.out.println("N: "+N+"    Nup: "+Nup+"    Ndown: "+Ndown);
        if(initialized && makingVideo){updateImg(s,t);}
    }
    
    /**
    *         setValue updates the spin value in the lattice and the magnetization
    *       with update of time in image layer .
    * 
    *  @param i - coordinate
    *  @param s - new spin value
    *  @param t - time 
    */ 
    @Override
    public void setValue(int i, int s, int t) {
        setValue(s,t);
    }
    
    public void setValue(int s){
        setValue(s,0);
    }
    public void setValue(int i,int j, int k,int s, int t){
        setValue(s,t);
    }
    /**
    *         getValue which gets the spin value in the lattice. 3d case
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    *  @param k - k coordinate
    */ 
    public int getValue(int i,int j, int k){
        return getMagnetization();
    }
    /**
    *         getValue which gets the spin value in the lattice. 3d case
    * 
    *  @param i - i coordinate
    */ 
    public int getValue(int i){
        return getMagnetization();
    }
    /**
    *         getValue which gets the spin value in the lattice. 2d case
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    */ 
    public int getValue(int i,int j){
        return getMagnetization();
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
    public void initialize(int s){
        if(s == 1){
            Nup = N; Ndown = 0;
        }else if(s == (-1)){
            Nup = 0; Ndown = N;
        }else if(s == 0){
            Nup = (N/2);
            Ndown = Nup;
        }else{
            Nup = (int)(N * Ran.nextDouble());
            Ndown = N = Nup;
        }
        magnetization = Nup-Ndown+NfixUp-NfixDown;
             
        // Create image layer
        initializeImg();
    }
    /**
    *         isThisFixed should return true if the lattice coordinate is fixed.
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    *  @param k - k coordinate 
    */ 
    public boolean isThisFixed(int i,int j,int k){
        return  false;
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
    public int getNinRange(){return N;}
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
    *         initializeImg creates the image layer of the lattice
    * 
    */  
    private void initializeImg(){
        scale = 5;
        int sizeImage = (int)(scale*Math.sqrt(N)*Math.sqrt(N));
        while(sizeImage < minImgSize){scale = scale*2;sizeImage = (int)(scale*Math.sqrt(N));}
        while(sizeImage > maxImgSize){scale = scale/2;sizeImage = (int)(scale*Math.sqrt(N));}

        latticeImg = new BufferedImage(sizeImage, sizeImage+ImageTextOffset, BufferedImage.TYPE_INT_RGB);  
        g = latticeImg.createGraphics();  
    }
    /**
    *         spinDraw completely redraws the image layer based on the current 
    *   configuration of the lattice.
    */
    private void spinDraw(){
        int spin;
        float baseCol = ((float)magnetization)/((float)N);
        Color col = new Color(baseCol,baseCol,baseCol);
        g.fillRect(0, 0+ImageTextOffset, (int)(scale*Math.sqrt(N)), (int)(scale*Math.sqrt(N)));		 
    }
    /**
    *         getSystemImg should return an image of the lattice
    */ 
    public BufferedImage getSystemImg(){
        return latticeImg;
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
    private void updateImg( int spin,int t){
        //erase old text
        g.setColor(Color.BLACK);
        g.fillRect(((int)(Math.sqrt(N))*scale*6/10),0,((int)(Math.sqrt(N))*scale*4/10),ImageTextOffset);
        g.setColor(Color.WHITE);
        if(t!=0){
        Font font = new Font("Courier New", Font.BOLD, 72);
        g.setFont(font);
        // 70 is the x coordinate of text
        g.drawString("t = "+ t, ((int) (Math.sqrt(N))*scale*6/10), ((int) ImageTextOffset*7/10));}

        float baseCol = ((float)magnetization)/((float)N);
        Color col = new Color(baseCol,baseCol,baseCol);
        g.fillRect(0, 0+ImageTextOffset, (int)(scale*Math.sqrt(N)), (int)(scale*Math.sqrt(N)));		 
       
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
    private void updateImg(int i, int j, int spin){
        updateImg(spin,0);
    }

    public void setFixed(int fix, int num, int sFix, boolean balanced){
        if(!balanced){
            System.out.println("FIXING "+fix+"  SPINS WITH S : "+sFix);
        }else{
            System.out.println("FIXING "+fix+"  SPINS WITH SUM EQUAL TO NEARLY OR AT ZERO ");
        }
        if( num > N){
            int dif = num-N;
            Nup += dif/2;
            Ndown += dif/2;
        }
        double fixProb = (double)(fix)/((double)num);
        for(int u = 0; u < num;u++){
            if(Ran.nextDouble() <= fixProb){
                // choose fix 
                if(!balanced){
                    if(sFix == 1){NfixUp++;}
                    if(sFix == (-1)){NfixDown++;}
                }else{
                    if(Ran.nextDouble() <= 0.5){
                        NfixUp++;
                    }else{
                        NfixDown++;
                    }
                }
                // Choose fix point in lattice to fix
                if(Ran.nextDouble() < 0.5 && Nup > 0){
                        Nup--;
                }else{
                        Ndown--;
                }
            }
        }
        Nfixed = NfixUp+NfixDown;
        //System.out.println("NBEFORE: "+N);
        N = Nup+Ndown+NfixUp+NfixDown;
        //System.out.println("NAFTER: "+N+ "    NupFix: "+ NfixUp+"     NdownFix: "+NfixDown);
        // using whole lattice as M;
        magnetization = Nup-Ndown+NfixUp-NfixDown;        
    }
    /**
    * makeVideo sets flag notifying images returned should be updated
    */
    @Override
    public void makeVideo(){makingVideo=true;spinDraw();}
    /**
    *         getNDynamic gives the amount of dynamic spins in the lattice
    */ 
    public int getNDynamic(){return N-NfixDown-NfixUp;}
    /**
    *         getNFixed gives the amount of fixed spins in the lattice
    */ 
    public int getNFixed(){return Nfixed;}
    /**
    *         getNFixedUp gives the amount of fixed spins up in the lattice
    */ 
    public int getNFixedUp(){return NfixUp;}
    /**
    *         getNFixedDown gives the amount of fixed spins down in the lattice
    */ 
    public int getNFixedDown(){return NfixDown;}
    /**
    *         getUpN gives the amount of spins up in the lattice
    */ 
    public int getUpN(){return Nup;}
    /**
    *         getDownN gives the amount of spins down in the lattice
    */ 
    public int getDownN(){return Ndown;}
    @Override
    public boolean highFixedAndStableMatched() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public int getHighMagIndex() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public int getHighFixIndex() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public double getLastAvgMagSubLat() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public int getSubLatSum(int i) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public int getFixedSubLatSum(int i) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    @Override
    public int getFixSpinVal() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    /**
    *         getMagStaggered gives the regular magnetization for single value lattice
    */ 
    public int getMagStaggered(){
        return magnetization;
    }    
}



package GPUBackend;

/**
 *   @(#) LatticeCL 
 *
 */

 // imports   
import Backbone.Util.DirAndFileStructure;
import Backbone.System.LatticeMagInt;
import Backbone.Util.ParameterBank;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.nio.IntBuffer;
import java.util.*;

import scikit.opencl.CLHelper;
import scikit.opencl.RandomNumberCL;

/**
 *  LatticeCL - An interface to lattice which is actually on OpenCL device
 *
 *
 *   <br>
 *
 * @author James B. Silva <jbsilva@bu.edu>
 * @since Apr 28, 2012
 */
public final class SingleValueLatticeCL implements LatticeMagInt {
    private int L = 256;
        private int N;
        private int D = 2;          // dimension only up to 3
        private int Geo = 2;       // geometry - 2 is rectangular
        private int fixSpinVal;
        private boolean allowDilution = false;
        private SortedMap<Integer,Integer> fixedSorted;
        private HashMap<Integer,Boolean> fixedLookup;
        private int instId=1;
        private int magnetization;
        private int sumNeigh;
        private DirAndFileStructure dir;
        private ParameterBank param=null;
        private BufferedImage latticeImg;
        private int ImageTextOffset = 90;
        private Graphics2D g;
        private int scale;
        private int nfixed=0;
        private int stable=-1;
        private String fixedFname;
        private int R = 1;
        private boolean initialized=false;
        private boolean useHeter;
        private boolean makingVideo=false;											//lattice else only sum is maintained
        private CLHelper clhandler;
        private String metropolisKernel = "ising2d_longmetro";
        private int s0;
        private String updateIntKernel="update_int_2buffer";
        
    /**
    *         LatticeCL constructor
    *
    *
    *  @param length - Length of lattice
    *  @param dim - Dimension of lattice
    *  @param g - Geometry of lattice; 2 is square lattice
    *  @param s - Initial spin of lattice; -2 is for random initialization
    */
    public SingleValueLatticeCL(int length,int rad, int dim,int g,int s, CLHelper cl){
        
        clhandler = cl;
        clhandler.createKernel("", updateIntKernel);
        clhandler.copyIntBufferAcrossKernel(metropolisKernel, 0, updateIntKernel, 0);
        clhandler.copyIntBufferAcrossKernel(metropolisKernel, 1, updateIntKernel, 1);
        clhandler.setKernelArg(updateIntKernel);
        
        instId = 0;
        R=1;
        L=length;
        D=dim;
        Geo=g;
        s0 = s;
        N = 0;
    //Rectangular lattice. Would implement single index but mod 
        //would cost computationally
        if(Geo==2){N=(int)Math.pow(L, D);}

        // Triangular/Honeycomb lattice
        if(Geo==3){
            System.err.println("LatticeCL | TRIANGULAR HONEYCOMB LATTICE NOT SUPPORTED: DEFAULTING TO SQUARE");
        }


        // If Heterogenous- Fix some spins
        useHeter = false; 
        fixedFname = "fixed.txt";
        
        // Initialize spins to spin state s
        initialize(s);
        initializeImg();
        
        initialized = true;
    }
        
    /**
    *         LatticeCL constructor.
    * 
    *  @param s - spin to initialize the lattice to
    */
    public SingleValueLatticeCL(int s, CLHelper clhandler){
        this(s,"","",1,clhandler,"",0,-1);
    }
    
    /**
    *         LatticeCL constructor.
    * 
    *  @param s - spin to initialize the lattice to
    *  @param id - instance id of currently running instance of program
    */ 
    public SingleValueLatticeCL(int s, int id, CLHelper clhandler){
        this(s,"","",id,clhandler,"",0,-1);
    }

    /**
    *         LatticeCL constructor.
    * 
    *  @param s - spin to initialize the lattice to
    *  @param postfix - postfix for the parameter file name
    */ 
    public SingleValueLatticeCL(int s,String postfix, CLHelper clhandler){
            this(s,postfix,"",1,clhandler,"",0,-1);
    }
    /**
    *         LatticeCL constructor.
    * 
    *  @param s - spin to initialize the lattice to
    *  @param postfix - postfix for the parameter file name
    *  @param id - instance id of currently running instance of program
    */ 
    public SingleValueLatticeCL(int s,String postfix,int id, CLHelper clhandler){
            this(s,postfix,"",id,clhandler,"",0,-1);
    }

    /**
    *         LatticeCL constructor.
    * 
    *  @param s - spin to initialize the lattice to
    *  @param postfix - postfix for the parameter file name
    *  @param fixedpost - postfix for the fixed file name
    *  @param id - instance id of currently running instance of program
    */ 
    public SingleValueLatticeCL(int s,String postfix,String fixedpost,int id, CLHelper clhandler){
            this(s,postfix,fixedpost,id,clhandler,"",0,-1);
    }
    
    /**
    *         LatticeCL constructor.
    * 
    *  @param s - spin to initialize the lattice to
    *  @param postfix - postfix for the parameter file name
    *  @param fnamepost - postfix for the fixed spin configuration file
    *  @param id - instance id of currently running instance of program
    */ 
    public SingleValueLatticeCL(int s,String postfix,String fnamepost,int id, CLHelper cl,String kernel, int argn, int arg2){
        
        clhandler = cl;
        
        if(kernel.equals("")){argn=0;}else{
        metropolisKernel =kernel;}
        if(arg2<0){arg2= argn+1;}
        
        clhandler.createKernel("", updateIntKernel);
        clhandler.copyIntBufferAcrossKernel(metropolisKernel, argn, updateIntKernel, 0);
        clhandler.copyIntBufferAcrossKernel(metropolisKernel, arg2, updateIntKernel, 1);
        clhandler.setKernelArg(updateIntKernel);
        
        param = new ParameterBank(postfix);
        instId = id;
        dir = new DirAndFileStructure();
        L = param.L;
        D= param.D;
        s0=s;
        Geo = param.Geo;
        N = 0;
        R = param.R;

        
        //Rectangular lattice. Would implement single index but mod 
        //would cost computationally
        if(Geo==2){N=(int)Math.pow(L, D);}

        // Triangular/Honeycomb lattice
        if(Geo==3){
            System.err.println("LatticeCL | TRIANGULAR HONEYCOMB LATTICE NOT SUPPORTED: DEFAULTING TO SQUARE");
        }


        // If Heterogenous- Fix some spins
        useHeter = param.useHeter; 

        fixedFname = "fixed"+fnamepost+".txt";
        
        // Initialize spins to spin state s
        initialize(s);
        initializeImg();
    
        // Update the lattice on device
        clhandler.runKernel(updateIntKernel, N, 1);
    
        initialized = true;
    }

    /**
    *         LatticeCL constructor.
    * 
    *  @param id - instance id of currently running instance of program
    */ 
    public SingleValueLatticeCL(int id, CLHelper cl,String kernel, int argn){
        this("","",id,cl,kernel,argn);
    }
    
    /**
    *         LatticeCL constructor.
    * 
    *  @param postfix - postfix for the parameter file name
    *  @param fnamepost - postfix for the fixed spin configuration file
    *  @param id - instance id of currently running instance of program
    */ 
    public SingleValueLatticeCL(String postfix,String fnamepost,int id, CLHelper cl,String kernel, int argn){
        this((new ParameterBank(postfix)).s,postfix,fnamepost,id,cl,kernel,argn,-1);
    }
    
    /**
    *         setValue updates the spin value in the lattice and the magnetization
    *       with update of time in image layer .
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    *  @param k - k coordinate
    *  @param s - new spin value
    *  @param t - time 
    */ 
    @Override
    public void setValue(int i, int j, int k, int s, int t) {
            IntBuffer spinBuff = (IntBuffer) clhandler.getArgHandler().getDirectIntBuffer(metropolisKernel, 0);
            spinBuff.put(i+j*L+k*L*L,s);
            if(t!=0){updateImgTime(t);}
            spinBuff.rewind();
            System.out.println("LatticeCL | Set value to :"+spinBuff.get(i+j*L+k*L*L)+"    at "+(i+j*L+k*L*L));
            spinBuff.rewind();
            clhandler.runKernel(updateIntKernel, N, 1);
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
        setValue(i%L,(i/L)%L,(i/(L*L))%L,s,t);
    }
    
    /**
    *         getValue which gets the spin value in the lattice. 3d case
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    *  @param k - k coordinate
    */ 
    @Override
    public int getValue(int i, int j, int k) {
            IntBuffer spinBuff = (IntBuffer) clhandler.getArgHandler().getDirectIntBuffer(metropolisKernel, 0);
            int val = spinBuff.get(i+j*L+k*L*L);
            spinBuff.rewind();
            
            return val;
    }

    /**
    *         getValue which gets the spin value in the lattice. 2d case
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    */ 
    @Override
    public int getValue(int i, int j) {
        return getValue(i,j,0);
    }

    /**
    *         getValue which gets the spin value in the lattice. 2d case
    * 
    *  @param i - i coordinate
    */ 
    @Override
    public int getValue(int i) {
        return getValue(i,0,0);
    }
    
    @Override
    /**
    *         getLength which gets the length of the lattice
    * 
    */ 
    public int getLength() {
        return N;
    }

    @Override
    /**
    *         getDimension gives the dimensionality of the lattice
    */ 
    public int getDimension() {
        return 1;
    }

    @Override
    /**
    *         getGeo gives the geometry of the lattice
    */ 
    public int getGeo() {
        return Geo;
    }

    @Override
    /**
    *         getN gives the size of the lattice
    */ 
    public int getN() {
        return N;
    }
    
    /**
    *         getNFixed gives the amount of fixed spins in the lattice
    */ 
    @Override
    public int getNFixed(){return nfixed;}

    @Override
    /**
    *         getRange gives the interaction range of the lattice
    */
    public int getRange() {
     return N;
    }

    @Override
    /**
    *         getMagnetization gives the magnetization of the lattice
    */ 
    public int getMagnetization() {
        return magnetization;
    }
    /**
    *         getMagStaggered gives the regular magnetization of the lattice for single value lattice.
    */ 
    public int getMagStaggered() {
        return magnetization;
    }
    
    @Override
    /**
    *         setFixedLatticeValues should set all fixed values in the lattice.
    */ 
    public synchronized void setFixedLatticeValues() {
    }

    private int calcMagnetization(){
        // assert starting fresh
        IntBuffer spinBuff = (IntBuffer) clhandler.getArgHandler().getDirectIntBuffer(metropolisKernel, 0);
        
        int sum = 0;
        for(int u = 0;u < (L*L);u++){
            sum  += spinBuff.get(u);
        }
        // Rewind buffer back to input state
        spinBuff.rewind();
        
        return sum;
    }

    /**
    *         initialize sets the lattice values to their initial values which is given as input.
    *
    *  @param s - Initial lattice value
    *
    */ 
    @Override
    public void initialize(int s) {
        s0=s;
        setInitialConfig(0, 0, "");
        clhandler.runKernel(updateIntKernel, N, 1);
    }

    @Override
    public boolean isThisFixed(int i, int j, int k) {
        boolean fixStatus = false;
        if(fixedLookup.containsKey((i+j*L+k*L*L))){
            if(fixedLookup.get((i+j*L+k*L*L))){fixStatus=true;}
        }
        return fixStatus;
    }

    @Override
    public int getNeighSum(int i, int j, int k) {
        return magnetization;
    }

    @Override
    public int getNeighSum(int i, int j) {
        return magnetization;
    }

    /**
    *         getNinRange gets the amount of neighbors for the simulation interaction range.
    *      2R+1 since this implementation only supports this type of range
    *
    *
    */ 
    @Override
    public int getNinRange(){
        return N;
    }

   @Override
    /**
    *         setInitialConfig should set all spins in the lattice to the values
    *   given by a file in the Config directory.
    * 
    *   @param t - time of config to set lattice to
    *   @param run - run of time to set lattice to
    *   @param  post - postfix or filename of lattice file
    */
    public synchronized void setInitialConfig(int t, int run, String post) {	 
       
    }

    @Override
    /**
    *         getInstId returns the identification number of the instance of this
    *  program running. Useful for allowing each instance a directory to work in.
    */ 
    public int getInstId() {
        return instId;
    }

    /**
    *         getSystemImg should return an image of the lattice
    */
    @Override
    public BufferedImage getSystemImg() {
            spinDraw();
            return latticeImg;
    }

    /**
    * makeVideo sets flag notifying images returned should be updated
    * 
    */
    @Override
    public void makeVideo() {
        makingVideo=true;
    }

    /**
    *         spinDraw completely redraws the image layer based on the current 
    *   configuration of the lattice.
    */
    private void spinDraw(){
        int spin;
        IntBuffer spinBuff = (IntBuffer) clhandler.getArgHandler().getDirectIntBuffer(metropolisKernel, 0);
        IntBuffer fixBuff = (IntBuffer) clhandler.getArgHandler().getDirectIntBuffer(metropolisKernel, 1);
            
        int currFixed,currSpin;
        int i,j;
        
        for(int u = 0;u < (L*L);u++){
            spin = spinBuff.get(u);
            currFixed = fixBuff.get(u);
            
            i= u%L;
            j= ((int)(u/L))%L;
            
            if(useHeter && currFixed ==1 ){
                    // Paint different color for fixed spins
                    // Cyan is like a blueish
                    // Magenta is purpleish
                    if(spin==1){g.setColor(Color.GRAY);} 
                    if(spin==(-1)){g.setColor(Color.RED);}
                    if(spin==(0)){g.setColor(Color.MAGENTA);}

            }else{
                if(spin==1){g.setColor(Color.YELLOW);}
                    if(spin==(-1)){g.setColor(Color.CYAN);}
                            if(spin==(0)){g.setColor(Color.WHITE);}}

            g.fillRect(i*scale, j*scale+ImageTextOffset, scale, scale);		 
        }
        
        // rewind buffers to default position
        fixBuff.rewind();
        spinBuff.rewind();
    }

    /**
    *         calcIsingEnergy calculates the energy in this configuration in the ising model.
    */
    public double calcIsingEnergy(){
        IntBuffer spinBuff = (IntBuffer) clhandler.getArgHandler().getDirectIntBuffer(metropolisKernel, 0);
            
        double energy = 0.0;
        
        ArrayList<Integer> spins = new ArrayList<Integer>();
        for(int u = 0;u < (L*L);u++){
            spins.add(spinBuff.get(u));
        }
        
        for(int i = 0;i < (L*L);i++){
            int x0 = (i)%(L);
            int y0 = (int) (i/L);
            y0=y0%L;
            int i0 = x0+y0*L;
            
            int leftX = x0+L-R;
            int leftY = y0+L-R;
            int index;
            
            if(R==0){
                index = ((leftX)%L)+(y0)*L;
                sumNeigh = sumNeigh +spins.get(index);

                index = ((leftX+2)%L)+(y0)*L;
                sumNeigh = sumNeigh +spins.get(index);

                index = (x0)+((leftY)%L)*L;
                sumNeigh = sumNeigh +spins.get(index);

                index = (x0)+((leftY+2)%L)*L;
                sumNeigh = sumNeigh +spins.get(index);

            }else{
                for(int u=0;u<(2*R+1);u++){
                    for(int v=0;v<(2*R+1);v++){
                        index = ((leftX+u)%L)+((leftY+v)%L)*L;
                        if(i0==index){}else{
                        sumNeigh = sumNeigh + spins.get(index);}
                    }
                }
            }
            
            // divide the bonds term by 2 since they will be double counted
            energy += -1.0*param.jInteraction*spins.get(i0)*sumNeigh/2
                    -1.0*param.hField*spins.get(i0);
      
        } 
        
        //Return buffer to input state
        spinBuff.rewind();
        
        return energy;
    }

    /**
    *         initializeImg creates the image layer of the lattice
    * 
    */  
    private void initializeImg(){

        // Standard scale of images
        scale = 5;
        int sizeImage = scale*L;
        while(sizeImage < 600){scale = scale*2;sizeImage = scale*L;}

        // Make the image
        latticeImg = new BufferedImage(sizeImage, sizeImage+ImageTextOffset, BufferedImage.TYPE_INT_RGB);  
        g = latticeImg.createGraphics(); 

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
    private void updateImgTime(int t){
        //erase old text
        g.setColor(Color.BLACK);
        g.fillRect(((int)L*scale*6/10),0,((int)L*scale*4/10),ImageTextOffset);
        g.setColor(Color.WHITE);
        if(t!=0){
        Font font = new Font("Courier New", Font.BOLD, 72);
        g.setFont(font);
        // 70 is the x coordinate of text
        g.drawString("t = "+ t, ((int) L*scale*6/10), ((int) ImageTextOffset*7/10));}
    }
    

    public int calcCorrectionMag(int s){
        String fileName;
       	fileName=fixedFname;
       	magnetization = L*L*s;
        fileName = dir.getFixedDirectory() + fileName;
        //Scanner scanner = new Scanner(fileName);
        Scanner scanner;
        int sum = 0;
            
        try {
            scanner = new Scanner(new File(fileName));

	    while(scanner.hasNextInt()) {
	    int x = scanner.nextInt();
	    int y = scanner.nextInt();
	    int z=0;
            
	    if(D==3){z = scanner.nextInt();}
	    
            int fSpin = scanner.nextInt();
	    nfixed++;
	        
	    	// sum used to fix magnetization
	        if(s==(-1*fSpin)){
                    sum += 2*fSpin;
	        }else if(s==0){sum += fSpin;}

            }
            
        }catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        magnetization=magnetization+sum;
    
        return magnetization;
    }

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

    public void setKernelName(String kernel){
        metropolisKernel = kernel;
    }
    
    
    // test the class
    public static void main(String[] args) {
    CLHelper clhandler = new CLHelper();
    clhandler.initializeOpenCL("CPU");
    
    int elementCount = 16;                                  // Length of arrays to process
    int lsize = 1;  // Local work size dimensions
    int gsize = 16;   // rounded up to the nearest multiple of the localWorkSize

    String kernel = "ising2d_longmetro";
    
    clhandler.createKernel("", kernel);

    clhandler.createIntBuffer(kernel, 0, 16, 0, 0);
    clhandler.createIntBuffer(kernel, 1, 16, 0, 0);
    clhandler.createIntBuffer(kernel, 2, 16, 0, 0);
    clhandler.createFloatBuffer(kernel, 0, 16, 1.0f, 0);
    clhandler.createFloatBuffer(kernel, 1, 16, 1.0f, 0);
    clhandler.createFloatBuffer(kernel, 2, 16, 1.0f, 0);    
    clhandler.createFloatBuffer(kernel, 3, 16, 1.0f, 0);    
    clhandler.createFloatBuffer(kernel, 4, 16, 0.0f, 0);    

    clhandler.setIntArg(kernel, 0,0);
    clhandler.setIntArg(kernel, 1,2);
    clhandler.setIntArg(kernel, 2,2);

    clhandler.setKernelArg(kernel);
    clhandler.runKernel(kernel, gsize, lsize);

    clhandler.getIntBufferAsArray(kernel, 0, 1, true);
    
    SingleValueLatticeCL lat = new SingleValueLatticeCL(4,1,2,2,1,clhandler);
    lat.setKernelName(kernel);
    
    System.out.println("i = 0 value is : "+lat.getValue(0, 0));
    System.out.println("i = 1 value is : "+lat.getValue(1, 0));
    lat.setValue(0,0,0,9,0);
    System.out.println("i = 0 value is : "+lat.getValue(0, 0));
    System.out.println("i = 1 value is : "+lat.getValue(1, 0));
    lat.setValue(0,0,0,4,0);
    System.out.println("i = 0 value is : "+lat.getValue(0, 0));
    System.out.println("i = 1 value is : "+lat.getValue(1, 0));
    System.out.println("i = 0 value is : "+lat.getValue(0, 0));
    System.out.println("i = 1 value is : "+lat.getValue(1, 0));
    
    
    clhandler.closeOpenCL();
    }
}


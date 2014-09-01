package Backbone.System;

/**
* 
*  @(#) SimpleLattice
* 
*/  
import Backbone.Util.DataSaver;
import Backbone.Util.DirAndFileStructure;
import Backbone.Util.ParameterBank;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.Vector;

/**        
*      Basic LatticeMagInt Class implemented as an int array. 
*  <br>
* 
*  @param N - System Size  (256 default)
*  @param D - Dimension (2 default)
*  @param Geo - Geometry,  2 is rectangular
*  @param Neighbors 
* 
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2011-09    
*/
public final class SimpleLattice implements LatticeMagInt{
    private int L = 256;
    private int N;
    private int D = 2;          // dimension only up to 3
    private int Geo = 2;       // geometry - 2 is rectangular
    private boolean allowDilution = false;
    private int instId=1;
    private Neighbors neigh;
    private int sumNeigh;
    private int [][][] lattice;   // lattice want to implement using just one index
    private DirAndFileStructure dir;
    private ParameterBank param=null;
    private BufferedImage latticeImg;
    private int ImageTextOffset = 90;
    private Graphics2D g;
    private int scale = 5;
    private int nfixed=0;
    private String fixedFname;
    private int R = 1;
    private int granularity=0; // divide lattice into sublattices
    private Vector<int[][][]> sublattices; // for sublattices 
    private int[] subLatticeSum;
    private boolean initialized=false;
    private boolean useHeter;
    private boolean makingVideo=false;
    private boolean completeGranularity = false; // true if need whole new structure for 
                                                //lattice else only sum is maintained
    private int magnetization;
    private int magStaggered;
    public int[][][] fixed;

    /**
    *         LatticeMagInt constructor
    *
    *
    *  @param length - Length of lattice
    *  @param rad - Interaction range
    *  @param dim - Dimension of lattice
    *  @param g - Geometry of lattice; 2 is square lattice
    *  @param s - Initial spin of lattice; -2 is for random initialization
    *
    */
    public SimpleLattice(int length,int rad, int dim,int g,int s){
        this(length,rad,dim,g,s,0,false);  			  
    }

    
    /**
    *         LatticeMagInt constructor
    *
    *
    *  @param length - Length of lattice
    *  @param rad - Interaction range
    *  @param dim - Dimension of lattice
    *  @param g - Geometry of lattice; 2 is square lattice
    *  @param s - Initial spin of lattice; -2 is for random initialization
    *  @param grandSub - granularity ie divide the lattice length into this amount 
    *  @param gran - true if keep track of lattice values in sublattices
    */
    public SimpleLattice(int length,int rad, int dim,int g,int s,int grandSub,boolean gran){
        granularity = grandSub;
        completeGranularity = gran;

        dir = new DirAndFileStructure();

        L = length;
        D= dim;
        Geo = g;
        lattice = null;
        N = 0;
        R = rad;

        //Rectangular lattice. Would implement single index but mod 
        //would cost computationally
        if(Geo==2){
            N=(int)Math.pow(L, D);

            if(D==2){
                if(granularity>0 && completeGranularity){
                    int setSize = (int) L/granularity;
                    setSize = setSize*setSize;
                    setSize = L*L/setSize;
                    int subL = (int) L/granularity;
                    sublattices = new Vector<int[][][]>();
                    for(int u=0;u < setSize;u++){
                        sublattices.add(new int[subL][subL][1]);
                    }
                    subLatticeSum = new int[subL*subL];

                    for(int u = 0;u < (subL*subL);u++){
                        subLatticeSum[u]=0;
                    }
                }else{
                    if(granularity > 0){ 
                        int subL = (int) L/granularity;subLatticeSum = new int[subL*subL];
                        for(int u = 0;u < (subL*subL);u++){
                            subLatticeSum[u]=0;
                        }
                    }
                    lattice = new int [L][L][1];
                }
            }
            
            if(D==3){
                if(granularity>0 && completeGranularity){
                    int setSize = (int) L/granularity;
                    setSize = setSize*setSize*setSize;
                    setSize = L*L*L/setSize;
                    sublattices = new Vector<int[][][]>();
                    for(int u=0;u < setSize;u++){
                            int subL = (int) L/granularity;
                            sublattices.add(new int[subL][subL][subL]);}
                }else{
                    if(granularity > 0){ 
                        int subL = (int) L/granularity;subLatticeSum = new int[subL*subL];
                        for(int u = 0;u < (subL*subL);u++){
                            subLatticeSum[u]=0;
                        }  
                    }
                    lattice = new int [L][L][L];
                }
            }
        }
        // Triangular/Honeycomb lattice
        if(Geo==3){}
        //Initialize spins to spin state s
        initialize(s);

        neigh = new Neighbors(this);
        initialized = true;			  
    }

    /**
    *         SimpleLattice constructor.
    * 
    *  @param s - spin to initialize the lattice to
    */
    public SimpleLattice(int s){
        this(s,"","",1);
    }
    
    /**
    *         SimpleLattice constructor.
    * 
    *  @param s - spin to initialize the lattice to
    *  @param id - instance id of currently running instance of program
    */ 
    public SimpleLattice(int s, int id){
        this(s,"","",id);
    }

    /**
    *         SimpleLattice constructor.
    * 
    *  @param s - spin to initialize the lattice to
    *  @param postfix - postfix for the parameter file name
    */ 
    public SimpleLattice(int s,String postfix){
        this(s,postfix,"",1);
    }
    /**
    *         SimpleLattice constructor.
    * 
    *  @param s - spin to initialize the lattice to
    *  @param postfix - postfix for the parameter file name
    *  @param id - instance id of currently running instance of program
    */ 
    public SimpleLattice(int s,String postfix,int id){
        this(s,postfix,"",id);
    }

    /**
    *         SimpleLattice constructor.
    * 
    *  @param s - spin to initialize the lattice to
    *  @param postfix - postfix for the parameter file name
    *  @param fnamepost - postfix for the fixed spin configuration file
    *  @param id - instance id of currently running instance of program
    */ 
    public SimpleLattice(int s,String postfix,String fnamepost,int id){
        param = new ParameterBank(postfix);
        instId = id;
        dir = new DirAndFileStructure();
        L = param.L;
        D= param.D;
        Geo = param.Geo;
        granularity = param.Granularity;
        completeGranularity = param.CompleteInfoGranular;
        lattice = null;
        N = 0;
        R = param.R;

        //Rectangular lattice. Would implement single index but mod 
        //would cost computationally
        if(Geo==2){
            N=(int)Math.pow(L, D);
            if(D==2){
                if(granularity>0 && completeGranularity){
                    int setSize = (int) L/granularity;
                    setSize = setSize*setSize;
                    setSize = L*L/setSize;
                    sublattices = new Vector<int[][][]>();
                    for(int u=0;u < setSize;u++){
                        int subL = (int) L/granularity;
                        sublattices.add(new int[subL][subL][1]);
                    }
                }else{
                    if(granularity > 0){ 
                        int subL = (int) L/granularity;subLatticeSum = new int[subL*subL];
                        for(int u = 0;u < (subL*subL);u++){
                            subLatticeSum[u]=0;
                        }
                    }
                    lattice = new int [L][L][1];}
            }
            if(D==3){
                if(granularity>0 && completeGranularity){
                    int setSize = (int) L/granularity;
                    setSize = setSize*setSize*setSize;
                    setSize = L*L*L/setSize;
                    sublattices = new Vector<int[][][]>();
                    for(int u=0;u < setSize;u++){
                        int subL = (int) L/granularity;
                        sublattices.add(new int[subL][subL][subL]);
                    }
                }else{
                    if(granularity > 0){ 
                        int subL = (int) L/granularity;subLatticeSum = new int[subL*subL];
                        for(int u = 0;u < (subL*subL);u++){
                            subLatticeSum[u]=0;
                        }  
                    }
                lattice = new int [L][L][L];}
            }
        }
        // Triangular/Honeycomb lattice
        if(Geo==3){}

        // Initialize spins to spin state s
        initialize(s);

        // If Heterogenous- Fix some spins
        useHeter = param.useHeter;
        if(useHeter){readFixedFile(fnamepost);} 

        spinDraw();
        neigh = new Neighbors(this);
        initialized = true;
    }

    /**
    * getValue which gets the spin value in the lattice.
    *
    * @param i - coordinate
    */
    public int getValue(int i) {
        int x = i%L;
        int y = (int)((int)(double)i/(double)L)%L;
        int z = (int)((int)(double)i/(double)(L*L))%L;
        return lattice[x][y][z];
    }
    
    @Override
    /**
    *         getValue which gets the spin value in the lattice. 3d case
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    *  @param k - k coordinate
    */ 
    public int getValue(int i,int j, int k){
        int retValue;
        if(granularity>0 && completeGranularity){
            retValue = getSublatticeValue(i,j,0);
        }else{
            retValue = lattice[i][j][k];
        }    
        return retValue;
    }

    @Override
    /**
    *         getValue which gets the spin value in the lattice. 2d case
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    */ 
    public  int getValue(int i,int j){
        int retValue;
        if(granularity>0 && completeGranularity){
            retValue = getSublatticeValue(i,j,0);			   
        }else{
            retValue = lattice[i][j][0];
        }
        
        return retValue;
    }


    @Override
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
    public synchronized void setValue(int i,int j, int k, int s,int t){
        if(getValue(i, j) ==s){return;}
        if(granularity>0 && completeGranularity){
            if(getValue(i, j) !=0){
                magnetization += 2*s;
            }else{
                magnetization +=s;
            }
            if((i+j+k)%2 == 0) { 
                magStaggered+= 2*s;
            }else{ 
                magStaggered -= 2*s;
            }
            setSublatticeValue(i,j,k,s);
            updateSublatticeSum(i,j,k,s);
        }else{
            if(granularity>0){updateSublatticeSum(i,j,k,s);}  
            if(getValue(i, j) !=0){
                magnetization += 2*s;
            }else{
                magnetization +=s;
            }
            if((i+j+k)%2 == 0) { 
                magStaggered+= 2*s;
            }else{ 
                magStaggered -= 2*s;
            }
            lattice[i][j][k]=s;}
        if(initialized && makingVideo){updateImg(i,j,s,t);}
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
        setValue(i%L,(int)((double)i/(double)L)%L,(int)((double)i/(double)(L*L))%L,s,t);
    }
    
    /**
    *         setValue updates the spin value in the lattice and the magnetization
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    *  @param k - k coordinate
    *  @param s - new spin value
    */ 
    public void setValue(int i,int j, int k, int s){
        if(getValue(i,j)==s){return;}  
        if(granularity>0 && completeGranularity){
            if(getValue(i, j) !=0){
                magnetization += 2*s;
            }else{
                magnetization +=s;
            }
            if((i+j+k)%2 == 0) { 
                magStaggered+= 2*s;
            }else{ 
                magStaggered -= 2*s;
            }
            setSublatticeValue(i,j,k,s);
            updateSublatticeSum(i,j,k,s);

        }else{
            if(granularity>0){updateSublatticeSum(i,j,k,s);}
                if(getValue(i, j) !=0){
                magnetization += 2*s;
            }else{
                magnetization +=s;
            }
            if((i+j+k)%2 == 0) { 
                magStaggered+= 2*s;
            }else{ 
                magStaggered -= 2*s;
            }
            lattice[i][j][k]=s;
        }
        if(initialized && makingVideo){updateImg(i,j,s);}
    }

    @Override
    /**
    *         getLength which gets the length of the lattice
    * 
    */ 
    public int getLength() {
        return L;
    }

    @Override
    /**
    *         getDimension gives the dimensionality of the lattice
    */ 
    public int getDimension() {
        return D;
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
        return R;
    }

    @Override
    /**
    *         getMagnetization gives the magnetization of the lattice
    */ 
    public int getMagnetization() {
        return magnetization;
    }

    /**
    *         getMagStaggered gives the staggered magnetization of the lattice
    */ 
    public int getMagStaggered(){return magStaggered;}
    /**
    *       getLastNeighSum gets the sum of spin values obtained in last call 
    *   of get Neighbors.
    * 
    * @return sum of last query of getting neighbors
    */
    public synchronized int getLastNeighSum(){return neigh.getLastNeighborSum();}
    
    
    @Override
    /**
    *         getNeighSum should return the sum of the spins within the interaction
    *   range centered by the coordinates given. 2d case
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    */ 
    public synchronized int getNeighSum(int i, int j){
        neigh.getNeighbors(i, j, 0);
        return neigh.getLastNeighborSum();
    }
       
    @Override
    /**
    *         getNeighSum should return the sum of the spins within the interaction
    *   range centered by the coordinates given. 3d case
    * 
    *  @param i - i coordinate
    *  @param j - j coordinate
    *  @param k - k coordinate 
    */ 
    public synchronized int getNeighSum(int i, int j,int k){
        neigh.getNeighbors(i, j, k);
        return neigh.getLastNeighborSum();
    }
    
    /**
    *         getNinRange gets the amount of neighbors for the simulation interaction range.
    *
    */ 
    @Override
    public int getNinRange(){
        return neigh.getNeighbors(1,0, 0).length;
    }

    /**
    *         getNeighbors gets the lattice values within the interaction range at the location given by input.
    *
    * returns all lattice values in range in the form of an integer array
    *
    *
    *  @param i - i coordinate
    *  @param j - j coordinate
    *  @param k - k coordinate
    *  
    *
    */ 
    public int[] getNeighbors(int i ,int j, int k){return neigh.getNeighbors(i, j, k);}
   
    /**
    *       getRandom returns random integer in given range
    * 
    * @param min - minimum integer in range
    * @param max - maximum integer in range
    * @return random number in range
    */
    public int getRandom(int min,int max){return (min + (int)(Math.random() * ((max - min) + 1)));}

    /**
    *      getSublatticeValue returns the value of the spin when using a complete
    *  information of lattice in sublattices. Only use with completegranularity is 
    *  true
    * 
    * 
    * @param i - i coordinate
    * @param j - j coordinate
    * @param k - k coordinate
    * @return spin value
    */
    private int getSublatticeValue(int i, int j, int k){
        int subI = (int) i/granularity;
        int subJ = (int) j/granularity;
        int subL = (int) L/granularity;
        int subK =0;

        if(D==3){subK= (int) k/granularity;}

        int subIndex= subI+subJ*subL+subK*subL*subL;

        int[][][] lat = sublattices.get(subIndex);

        return lat[i-subI*granularity][j-subJ*granularity][k-subK*granularity];
    }

    /**
    *          getSublatticeSum returns the value of the sum of the spins in the 
    *   sublattice that the given coordinates is located in.
    * 
    * 
    * @param i - i coordinate
    * @param j - j coordinate
    * @param k - k coordinate
    * @return sum of sublattice of given index location
    */
    public int getSublatticeSum(int i, int j, int k){
        int subI = (int) i/granularity;
        int subJ = (int) j/granularity;
        int subL = (int) L/granularity;
        int subK =0;
                        if(D==3){subK= (int) k/granularity;}

        int subIndex= subI+subJ*subL+subK*subL*subL;

        return  subLatticeSum[subIndex];
    }


    /**
    *       setSublatticeValue updates the sublattice collection that holds
    *   each spin value in a sublattice.
    * 
    * @param i - i coordinate
    * @param j - j coordinate
    * @param k - k coordinate
    * @param s - new spin value 
    */
    private void setSublatticeValue(int i, int j, int k,int s){
        int subI = (int) i/granularity;
        int subJ = (int) j/granularity;
        int subL = (int) L/granularity;
        int subK =0;
        if(D==3){subK= (int) k/granularity;}

        int subIndex= subI+subJ*subL+subK*subL*subL;

        int[][][] lat = sublattices.get(subIndex);

        lat[i-subI*granularity][j-subJ*granularity][k-subK*granularity] = s;

        sublattices.setElementAt(lat, subIndex);
    }

    /**
    *       updateSublatticeSum updates the sublattice collection that holds
    *   the sum of the spins in each sublattice.
    * 
    * @param i - i coordinate
    * @param j - j coordinate
    * @param k - k coordinate
    * @param s - new spin value 
    */
    private void updateSublatticeSum(int i, int j, int k, int s){
        int subI = (int) i/granularity;
        int subJ = (int) j/granularity;
        int subL = (int) L/granularity;
        int subK =0;
        if(D==3){subK= (int) k/granularity;}

        int subIndex= subI+subJ*subL+subK*subL*subL;

        if(!initialized){
                subLatticeSum[subIndex] = subLatticeSum[subIndex]+s;
        }else{subLatticeSum[subIndex] = subLatticeSum[subIndex]+2*s;}
    }

    /**
    *         initialize sets the lattice values to their initial values which is given as input.
    *
    *
    *  @param s - Initial lattice value; -2 for random values
    *
    */ 
    @Override
    public void  initialize(int s){
        int spin = s;int sum=0;int sumStag = 0;makingVideo=false;
        // Square lattice initialize all lattice elements to value s
        // unless s is -2 which is a flag for random initialization of spins
        if(Geo==2){
            for(int i = 0;i < L;i++){ 
            for(int j = 0;j < L;j++){
                if(s==(-2)){
                    if(allowDilution){
                            spin = getRandom(-1,1);
                    }else{
                        spin = getRandom(0,1);if(spin ==0){spin=-1;}
                    }
                }	 
                if(D==3){
                    for(int k = 0;k < L;k++){setValue(i,j,k,spin);sum +=spin;}
                }else{
                    setValue(i,j,0,spin);sum +=spin;
                    sumStag += ((i+j)%2==0) ? spin: -1*spin;
                }
            }}
        }
        // Add for non square geometries
        // Create image layer
        initializeImg();

        magnetization = sum;
        magStaggered = sumStag;
    }


    /**
    *         initializeImg creates the image layer of the lattice
    * 
    */  
    private void initializeImg(){
        scale = 5;
        int sizeImage = scale*L;
        while(sizeImage < 600){scale = scale*2;sizeImage = scale*L;}

        latticeImg = new BufferedImage(sizeImage, sizeImage+ImageTextOffset, BufferedImage.TYPE_INT_RGB);   
        g = latticeImg.createGraphics();  
    }

    /**
    *         spinDraw completely redraws the image layer based on the current 
    *   configuration of the lattice.
    */
    private void spinDraw(){
        int spin;
        for(int i = 0;i < L;i++){ 
        for(int j = 0;j < L;j++){
            if(this.isThisFixed(i, j, 0) ){
                // Paint different color for fixed spins
                spin = getValue(i, j);
                if(spin==1){g.setColor(Color.BLUE);} 
                if(spin==(-1)){g.setColor(Color.RED);}
                if(spin==(0)){g.setColor(Color.MAGENTA);}
            }else{
                spin = getValue(i, j);
                if(spin==1){g.setColor(Color.WHITE);}
                    if(spin==(-1)){g.setColor(Color.BLACK);}
                            if(spin==(0)){g.setColor(Color.YELLOW);}
            }
            g.fillRect(i*scale, j*scale+ImageTextOffset, scale, scale);		 
        }}
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
    private void updateImg(int i, int j, int spin,int t){
        //erase old text
        g.setColor(Color.BLACK);
        g.fillRect(((int)L*scale*6/10),0,((int)L*scale*4/10),ImageTextOffset);
        g.setColor(Color.WHITE);
        if(t!=0){
        Font font = new Font("Courier New", Font.BOLD, 72);
        g.setFont(font);
        // 70 is the x coordinate of text
        g.drawString("t = "+ t, ((int) L*scale*6/10), ((int) ImageTextOffset*7/10));}
        spin = getValue(i, j);
        if(spin==1){g.setColor(Color.WHITE);}else 
            if(spin==(-1)){g.setColor(Color.BLACK);}else
                    if(spin==(0)){g.setColor(Color.YELLOW);}
        g.fillRect(i*scale, j*scale+ImageTextOffset, scale, scale);
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
        updateImg(i,j,spin,0);
    }

    /**
    *         getSystemImg should return an image of the lattice
    */ 
    public BufferedImage getSystemImg(){
        return latticeImg;
    }

    public void saveLatticeImage(String name){
        DataSaver dat = new DataSaver(instId);
        dat.saveImage(latticeImg, "png",name);
    }

    /**
    * printConfiguration                           (1)
    *
    * Prints out the current configuration
    *
    */ 
    public void printConfiguration(){
        for (int i = 0; i < L; i++)
        for (int j = 0; j < L; j++){
            if(D==3){
                for (int k = 0; k < L; k++){			    		  
                    System.out.println(i + "   " + j + "  "+ k+ "   "+  getValue(i,j,0));
                }
            }else{
                System.out.println(i + "   " + j + "  "+   getValue(i,j,0));
            }
        }
    }

    /**
    *         isThisFixed determines if input lattice coordinates are for a fixed lattice element.
    *
    * returns indexed location
    *
    *
    *  @param i - i coordinate
    *  @param j - j coordinate
    *  @param k - k coordinate
    *  
    *
    */ 
    @Override
    public boolean isThisFixed(int i, int j, int k){
        boolean yesItIs = false;
        if(useHeter){if(fixed[i][j][k] == 1){yesItIs = true;}}
        return yesItIs;
    }

    /**
    *         isPositive determines if number is positive.
    *
    *  @param number 
    *  
    *
    */ 
    public int isPositive(int number){if (number >= 0){return 1;}else{return 0;}}


     
    /**
    *         setInitialConfig should set all spins in the lattice to the values
    *   given by a file in the Config directory.
    * 
    *   @param t - time of config to set lattice to
    *   @param run - run of time to set lattice to
    */
    public void setInitialConfig(int t,int run){
        setInitialConfig(t,run,"");
    }

    /**
    *         getInstId returns the identification number of the instance of this
    *  program running. Useful for allowing each instance a directory to work in.
    */ 
    @Override
    public int getInstId(){return instId;}

    //TODO ::   Read and Empty fixed file

    /**
    *         setInitialConfig sets the lattice to the desired configuration given by the run and time in the simulation.
    *
    *
    *  @param t - Time in simulation to initialize lattice to.
    *  @param run - Run which simulation should initialize to time t for. 
    *  @param appendStr - (optional) string post fix for filename
    *  
    *
    */
    @Override
    public void setInitialConfig(int t,int run,String appendStr){
        if(t==0){
            initialize(param.s);
        }else{
            try {
                String fname = dir.getLatticeConfig(instId, t, run, appendStr);
                Scanner scanner = new Scanner(new File(fname));
                int sum = 0;int sumStag = 0;

                while(scanner.hasNextInt()) {
                    if (D==2){
                        int x = scanner.nextInt();
                        int y = scanner.nextInt();
                        int fSpin = scanner.nextInt();
                        setValue(x,y,0,fSpin);
                        sum = fSpin + sum;
                        sumStag += ((x+y)%2==0) ? fSpin: -1*fSpin;
                    }else{				
                        int x = scanner.nextInt();
                        int y = scanner.nextInt();
                        int z = scanner.nextInt();
                        int fSpin = scanner.nextInt();
                        setValue(x,y,z,fSpin);
                        sum = fSpin + sum; 
                        sumStag += ((x+y+z)%2==0) ? fSpin: -1*fSpin;
                    }
                }
                magnetization = sum;
                magStaggered = sumStag;
                //System.out.println("SimpleLattice | value " + s)
            } catch (FileNotFoundException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
            }
        }
        // If Heterogenous- Fix some spins in the fixed lattice and update the magnetization
        if(useHeter){setFixedLatticeValues();}
    }

    /**
    *         readFixedFile reads the file "fixed.txt" which contains the fixed lattice values and updates the lattice with this information.
    *
    *
    */
    public void readFixedFile(String fname){
        fixedFname = "fixed"+fname+".txt";
        // Setting some spins fixed using File
        if(D == 3){
            fixed = new int[L][L][L];
        }else if(D==2){
            fixed = new int[L][L][1];
        }
        initializeFixedArray();
        setFixedLatticeValues();
    }

    /**
    *         setFixedLatticeValues sets the fixed lattice values in the lattice given by the "fixed.txt" file.
    *
    *
    */
    @Override
    public void setFixedLatticeValues(){
	String fileName;
	fileName=fixedFname;
        fileName = dir.getFixedDirectory() + fileName;
        //Scanner scanner = new Scanner(fileName);
        Scanner scanner;
        int sum = 0;
        nfixed = 0;
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
                if(getValue(x,y,z)==(-1*fSpin)){
                        sum += 2*fSpin;
                }else if(getValue(x,y,z)==0){
                    sum += fSpin;
                }
                setValue(x,y,z,fSpin);
                fixed[x][y][z] = 1;	
	    }
	} catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
	}
        magnetization = calcMagnetization()[0];
        magStaggered = calcMagnetization()[1];
    }

    private int[] calcMagnetization(){
        // assert starting fresh
        int[] sum = new int[2];
        for(int u=0;u<L;u++){
	for(int v=0;v<L;v++){
            sum[0] += getValue(u,v);
            sum[1] += ((u+v)%2==0) ? getValue(u,v): -1*getValue(u,v);
        }}
        return sum;
    }

    /**
    *         initializeFixedArray sets up the array of fixed lattice values.
    *
    */
    private void initializeFixedArray(){
        for(int i=0;i<L;i++){
        for(int j=0;j<L;j++){
            if(D==3){
                for(int k=0;k<L;k++){
                    fixed[i][j][k] = 0;
                }
            }else{
                fixed[i][j][0] = 0;
            }
        }}		    
    }

    /**
    * makeVideo sets flag notifying images returned should be updated
    * 
    */
    @Override
    public void makeVideo(){makingVideo=true;spinDraw();}
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
        // TODO Auto-generated method stub
        return 0;
    }
    @Override
    public int getSubLatSum(int i) {
        // TODO Auto-generated method stub
        return 0;
    }
    @Override
    public int getFixedSubLatSum(int i) {
        // TODO Auto-generated method stub
        return 0;
    }
    @Override
    public int getFixSpinVal() {
        // TODO Auto-generated method stub
        return 0;
    }

    // test the class
    public static void main(String[] args) {
        SimpleLattice lat = new SimpleLattice(4,1,2,2,1,2,true);


        System.out.println("SimpleLattice | sum:"+lat.getSublatticeSum(2,2,0));

        lat.printConfiguration();
        System.out.println("SimpleLattice | __________________________________");
        System.out.println(lat.magnetization);
        lat = new SimpleLattice(1);
        lat.getNeighbors(3, 4, 0);
        System.out.println("SimpleLattice |  Sum of neighbors"+lat.getLastNeighSum());
        System.out.println("SimpleLattice |  Sum of neighbors"+lat.getNeighSum(3,4));

        System.out.println("SimpleLattice | __________________________________");
        System.out.println(lat.magnetization);
        lat.saveLatticeImage("");
        lat.setValue(3, 4, 0, -1);
        lat.setValue(3, 5, 0, -1);
        lat.setValue(3, 6, 0, -1);
        lat.setValue(4, 4, 0, -1);
        lat.getNeighbors(3, 4, 0);
        System.out.println("SimpleLattice |  Sum of neighbors"+lat.getLastNeighSum());
        System.out.println("SimpleLattice | N in range: "+lat.getNinRange()+"  with R : "+lat.R);
        lat.saveLatticeImage("post");

        System.out.println("SimpleLattice | Done!");
    }
}

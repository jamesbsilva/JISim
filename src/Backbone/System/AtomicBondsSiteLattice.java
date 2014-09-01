package Backbone.System;

/**
 *
 * @(#) AtomicBondsSiteLattice
 *
 */
import Backbone.Util.DataSaver;
import Backbone.Util.DirAndFileStructure;
import Backbone.Util.ParameterBank;
import Backbone.Visualization.SquareLattice2D;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Random;
import java.util.Scanner;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicIntegerArray;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * AtomicBondsSiteLattice class implements lattice as AtomicIntegerArray. It also
 * implements functions to determine sector of highest stable spin density and
 * for fixed spin systems the highest fixed spin density. <br>
 *
 * @param N - System Size (256 default)
 * @param D - Dimension (2 default)
 * @param Geo - Geometry, 2 is rectangular
 * @param Neighbors
 *
 * @author James Silva <jbsilva @ bu.edu>
 * @since 2012-01
 */
public final class AtomicBondsSiteLattice implements LatticeInt {
    private int L = 256;
    private int N;
    private int D = 2;          // dimension only up to 3
    private int Geo = 2;       // geometry - 2 is rectangular
    private int fixSpinVal;
    private boolean allowDilution = false;
    private int instId = 1;
    private Random Ran = new Random();
    private AtomicIntegerArray lattice;   // lattice want to implement using just one index
    private ArrayList<Double> latticeBonds;   // lattice want to implement using just one index
    private AtomicInteger percSize;
    private AtomicIntegerArray fixed;
    private DirAndFileStructure dir;
    private ParameterBank param = null;
    private SquareLattice2D vis;
    private Graphics2D g;
    private int scale;
    private int nfixed = 0;
    private int stable = -1;
    private String fixedFname;
    private int R = 1;
    private boolean initialized = false;
    private boolean useHeter;
    private boolean makingVideo = true;											//lattice else only sum is maintained
    private double[] bondsSite;
    

    
    /**
    * AtomicLattice constructor.
    *
    * @param s - spin to initialize the lattice to
    */
    public AtomicBondsSiteLattice(int s) {
        this(s, "", "", 1);
    }

    /**
    * AtomicLattice constructor.
    *
    * @param s - spin to initialize the lattice to
    * @param id - instance id of currently running instance of program
    */
    public AtomicBondsSiteLattice(int s, int id) {
        this(s, "", "", id);
    }

    /**
     * AtomicLattice constructor.
     *
     * @param s - spin to initialize the lattice to
     * @param postfix - postfix for the parameter file name
     */
    public AtomicBondsSiteLattice(int s, String postfix) {
        this(s, postfix, "", 1);
    }

    /**
    * AtomicLattice constructor.
    *
    * @param s - spin to initialize the lattice to
    * @param postfix - postfix for the parameter file name
    * @param id - instance id of currently running instance of program
    */
    public AtomicBondsSiteLattice(int s, String postfix, int id) {
        this(s, postfix, "", id);
    }

    /**
    * AtomicLattice constructor.
    *
    * @param s - spin to initialize the lattice to
    * @param postfix - postfix for the parameter file name
    * @param fnamepost - postfix for the fixed spin configuration file
    * @param id - instance id of currently running instance of program
    */
    public AtomicBondsSiteLattice(int s, String postfix, String fnamepost, int id) {
        param = new ParameterBank(postfix);
        instId = id;
        dir = new DirAndFileStructure();
        L = param.L;
        D = param.D;
        // stable opposite of initialized, assumption
        stable = -s;
        Geo = param.Geo;
        lattice = null;
        N = 0;
        R = param.R;
        if (param.useLongRange == false) {
            R = 0;
        }
        percSize = new AtomicInteger(0);
        
        //Rectangular lattice. Would implement single index but mod 
        //would cost computationally
        if (Geo == 2 || Geo == 4) {
            N = (int) Math.pow(L, D);
            int coordination = 1;
            if (D == 2) {
                lattice = new AtomicIntegerArray(L * L);
                fixed = new AtomicIntegerArray(L * L);
                latticeBonds = new ArrayList<Double>();
                coordination = 4;
            }else if (D == 3) {
                lattice = new AtomicIntegerArray(L * L * L);
                fixed = new AtomicIntegerArray(L * L * L);
                latticeBonds = new ArrayList<Double>();
                coordination = 6;
            }
            for(int u = 0; u < (coordination*lattice.length());u++){
                latticeBonds.add(2.0);
            }
            bondsSite = new double[coordination];
            // Triangular/Honeycomb lattice
        }else if (Geo == 3) {
        }

        // Initialize spins to spin state s
        initialize(s);

        // If Heterogenous- Fix some spins
        useHeter = param.useHeter;
        if (useHeter) {
            readFixedFile(fnamepost);
        }
        
        // Initial draw of lattice image
        vis = new SquareLattice2D(L,R,1,true,true);
        vis.initializeImg();
        vis.spinDraw(this);
        
        // throw initial bonds
        throwBonds();
        
        // done with initialization
        initialized = true;
    }

    @Override
    /**
    * getValue which gets the spin value in the lattice. 3d case
    *
    * @param i - i coordinate
    * @param j - j coordinate
    * @param k - k coordinate
    */
    public int getValue(int i, int j, int k) {
        return lattice.get(i + j * L + k * L * L);
    }

    /**
    * getValue which gets the spin value in the lattice.
    *
    * @param i - coordinate
    */
    public int getValue(int i) {
        return lattice.get(i);
    }
    
    @Override
    /**
    * setValue updates the spin value in the lattice and the magnetization with
    * update of time in image layer .
    *
    * @param i - i coordinate
    * @param j - j coordinate
    * @param k - k coordinate
    * @param s - new spin value
    * @param t - time
    */
    public void setValue(int i, int j, int k, int s, int t) {
        if ( getValue(i, j) == s ) {
            return;
        }
        // Magnetization changes differently in diluted site
        if ( getValue(i,j,k) == 0 ) {
            if(s == 1){percSize.set(percSize.addAndGet(1));}
        }else if( getValue(i,j,k) == (-1)){
            if(s == 0){percSize.set(percSize.addAndGet(-1));}
        } 

        //System.out.println("Setting | "+getValue(i, j)+"    TO | "+percSize.get()+"     makeVid | "+makingVideo+"    init | "+initialized);
        
        // Change the spin at site    
        lattice.set(i + j * L + k * L * L, s);
        
        // If doing video update the image of lattice including the time
        if (initialized && makingVideo) {
            updateImg(i, j, s, t);
        }
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
    public void setValue(int i, int s, int t){
        if(D == 2){
            setValue(i%L,(int)((double)i/(double)L)%L,0,s,t);
        }else{
            setValue(i%L,(int)((double)i/(double)L)%L,(int)((double)i/(double)(L*L))%L,s,t);
        }
    }
    
    /**
    * setValue which updates the spin value in the lattice and the
    * magnetization
    *
    *
    * @param i - i coordinate
    * @param j - j coordinate
    * @param k - k coordinate
    * @param s - new spin value
    */
    public void setValue(int i, int j, int k, int s) {
        setValue(i,j,k,s,0);
    }

    @Override
    /**
    * getValue which gets the spin value in the lattice. 2d case
    *
    * @param i - i coordinate
    * @param j - j coordinate
    */
    public int getValue(int i, int j) {
        return lattice.get(i + j * L);
    }

    public int getRandom(int min, int max) {
        return (min + (int) (Math.random() * ((max - min) + 1)));
    }

    @Override
    /**
    * getLength which gets the length of the lattice
    *
    */
    public int getLength() {
        return L;
    }

    @Override
    /**
    * getDimension gives the dimensionality of the lattice
    */
    public int getDimension() {
        return D;
    }

    @Override
    /**
    * getGeo gives the geometry of the lattice
    */
    public int getGeo() {
        return Geo;
    }

    @Override
    /**
    * getN gives the size of the lattice
    */
    public int getN() {
        return N;
    }

    /**
    * getNFixed gives the amount of fixed spins in the lattice
    */
    @Override
    public int getNFixed() {
        return nfixed;
    }

    @Override
    /**
    * getRange gives the interaction range of the lattice
    */
    public int getRange() {
        return R;
    }

    @Override
    /**
    *  getMagnetization gives the magnetization of the lattice
    */
    public int getMagnetization() {
        return percSize.get();
    }

    public double[] getBonds(int i, int j , int k){
        for(int u = 0; u < bondsSite.length;u++){
            bondsSite[u] = latticeBonds.get(bondsSite.length*(i+j*L+k*L*L)+u);
        }
        return bondsSite;
    }
    
    /**
    * getSystemImg should return an image of the lattice
    */
    @Override
    public BufferedImage getSystemImg() {
        return vis.getImageOfVis();
    }

    /**
    * makeVideo should set the lattice to begin making video which should
    * require updating of the lattice image.
    */
    @Override
    public void makeVideo() {
        makingVideo = true;
    }

    /**
    * saveLatticeImage saves the image of the lattice with the given name
    *
    * @param name - filename
    */
    public void saveLatticeImage(String name) {
        DataSaver dat = new DataSaver(instId);
        dat.saveImage(vis.getImageOfVis(), "png", name);
    }

    @Override
    /**
    * setFixedLatticeValues should set all fixed values in the lattice.
    */
    @SuppressWarnings("CallToThreadDumpStack")
    public synchronized void setFixedLatticeValues() {
        String fileName;
        fileName = fixedFname;
        int fsum = 0;
        fileName = dir.getFixedDirectory() + fileName;
        //Scanner scanner = new Scanner(fileName);
        Scanner scanner;
        int sum = 0;
        nfixed = 0;
        try {
            scanner = new Scanner(new File(fileName));
            while (scanner.hasNextInt()) {
                int x = scanner.nextInt();
                int y = scanner.nextInt();
                int z = 0;
                if (D == 3) {
                    z = scanner.nextInt();
                }
                int fSpin = scanner.nextInt();
                fsum += fSpin;
                nfixed++;

                // sum used to fix magnetization
                if ( (getValue(x, y, z) == 0) && (fSpin == 1) ) {
                    sum += 1;
                } else if ( (getValue(x, y, z) == 1) && (fSpin == 0)  ) {
                    sum -= 1;
                }
         
                setValue(x, y, z, fSpin);
         
                fixed.set(x + y * L + z * L * L, 1);        
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(AtomicBondsSiteLattice.class.getName()).log(Level.SEVERE, null, ex);
        }
        if (nfixed > 0) {
            fixSpinVal = (int) (fsum / nfixed);
        }

        // Assert that magnetization is recalculated fresh
        percSize.addAndGet(sum);
    }

    /*
    *   throwBonds throws the random bond values.
    */
    public void throwBonds(){
        for(int u = 0; u < latticeBonds.size();u++){
            latticeBonds.set(u,Ran.nextDouble());
        }
    }
    
    @Override
    /**
    * initialize sets the lattice values to their initial values which is given
    * as input.
    *
    *
    * @param s - Initial lattice value; -2 for random values
    *
    */
    public synchronized void initialize(int s) {
        int spin = s;
        int sum = 0;

        // Square lattice initialize all lattice elements to value s
        // unless s is -2 which is a flag for random initialization of spins
        if (Geo == 2) {
            for (int i = 0; i < L; i++) {
                for (int j = 0; j < L; j++) {
                    if (s == (-2)) {
                        if (allowDilution) {
                            spin = getRandom(-1, 1);
                        } else {
                            spin = getRandom(0, 1);
                            if (spin == 0) {
                                spin = -1;
                            }
                        }
                    }
                    if (D == 3) {
                        for (int k = 0; k < L; k++) {
                            setValue(i, j, k, spin);
                            sum += spin;
                        }
                    } else {
                        setValue(i, j, 0, spin);
                        sum += spin;
                    }
                }
            }
        }
        // Add for non square geometries

    }


    /**
    * updateImg just updates the image layer instead of redrawing the layer.
    * Change time text of image
    *
    * @param i - i coordinate
    * @param j - j coordinate
    * @param spin - new spin value
    * @param t - new time in image
    *
    */
    private synchronized void updateImg(int i, int j, int spin, int t) {
        vis.updateImg(i, j, spin, t);
    }

    @Override
    /**
    * isThisFixed should return true if the lattice coordinate is fixed.
    *
    * @param i - i coordinate
    * @param j - j coordinate
    * @param k - k coordinate
    */
    public boolean isThisFixed(int i, int j, int k) {
        return isThisFixed(i + j * L + k * L * L);
    }

    /**
    * isThisFixed should return true if the lattice coordinate is fixed.
    *
    * @param i - i coordinate
    */
    public boolean isThisFixed(int i) {
        if (fixed.get(i) == 1) {
            return true;
        } else {
            return false;
        }
    }


    @Override
    /**
    * getNinRnage should return the amount of spins within the interaction
    * range. Useful for metropolis algorithm in long range case.
    */
    public int getNinRange() {
        if (param.useLongRange) {
            return ((2 * R + 1) * (2 * R + 1) - 1);
        } else {
            return 4;
        }
    }

    @Override
    /**
    * setInitialConfig should set all spins in the lattice to the values given
    * by a file in the Config directory.
    *
    * @param t - time of config to set lattice to
    * @param run - run of time to set lattice to
    * @param post - postfix or filename of lattice file
    */
    @SuppressWarnings("CallToThreadDumpStack")
    public synchronized void setInitialConfig(int t, int run, String post) {
        if (t == 0) {
            initialize(param.s);
        } else {
            // Open file and set spins to values in file
            try {
                // using data saver to get file location
                String fname = dir.getLatticeConfig(instId, t, run, post);
                Scanner scanner = new Scanner(new File(fname));
                int sum = 0;
                while (scanner.hasNextInt()) {
                    int x = scanner.nextInt();
                    int y = scanner.nextInt();
                    int z = 0;
                    if (D == 3) {
                        z = scanner.nextInt();
                    }
                    int fSpin = scanner.nextInt();
                    setValue(x, y, z, fSpin);
                    sum = fSpin + sum;
                    
                }
                //System.out.println("AtomicLattice | value " + s)
            } catch (FileNotFoundException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
        // If Heterogenous- Fix some spins in the fixed lattice and update the magnetization
        if (param.useHeter) {
            setFixedLatticeValues();
        }

    }

    @Override
    /**
    * getInstId returns the identification number of the instance of this
    * program running. Useful for allowing each instance a directory to work
    * in.
    */
    public int getInstId() {
        return instId;
    }

    /**
    * isPositive determines if number is positive.
    *
    * @param number
    *
    */
    public int isPositive(int number) {
        if (number >= 0) {
            return 1;
        } else {
            return 0;
        }
    }

    /**
    * readFixedFile reads the file "fixed.txt" which contains the fixed lattice
    * values and updates the lattice with this information.
    *
    *
    */
    public void readFixedFile(String fname) {
        fixedFname = "fixed" + fname + ".txt";

        // Setting some spins fixed using File
        if (D == 3) {
            fixed = new AtomicIntegerArray(L * L * L);
        } else if (D == 2) {
            fixed = new AtomicIntegerArray(L * L);
        }
        initializeFixedArray();
        setFixedLatticeValues();
    }

    /**
    * initializeFixedArray sets up the array of fixed lattice values.
    *
    */
    private void initializeFixedArray() {
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                if (D == 3) {
                    for (int k = 0; k < L; k++) {
                        fixed.set(i + j * L + k * L * L, 0);
                    }
                } else {
                    fixed.set(i + j * L, 0);
                }
            }
        }
    }

    /**
    * printConfiguration (1)
    *
    * Prints out the current configuration
    *
    */
    public void printConfiguration() {
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {
                if (D == 3) {
                    for (int k = 0; k < L; k++) {
                        System.out.println(i + "   " + j + "  " + k + "   " + getValue(i, j, 0));
                    }
                } else {
                    System.out.println(i + "   " + j + "  " + getValue(i, j, 0));
                }
            }
        }
    }

    @Override
    public int getFixSpinVal() {
        return fixSpinVal;
    }

    @Override
    public int getMagStaggered() {
        return percSize.get();
    }

    // test the class
    public static void main(String[] args) {
    }
}

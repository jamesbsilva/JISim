package Backbone.System;

/**
*
* @(#) AtomicLattice
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
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicIntegerArray;

/**
* Atomic LatticeMagInt class implements lattice as AtomicIntegerArray. It also
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
public final class AtomicLattice implements LatticeMagInt {

    private int L = 256;
    private int N;
    private int D = 2;          // dimension only up to 3
    private int Geo = 2;       // geometry - 2 is rectangular
    private int fixSpinVal;
    private boolean allowDilution = false;
    private int instId = 1;
    private Neighbors neigh;
    private int sumNeigh;
    private AtomicIntegerArray lattice;   // lattice want to implement using just one index
    private AtomicInteger magnetization;
    private AtomicInteger magStaggered;
    private AtomicInteger maxFixedSum;
    private AtomicIntegerArray fixed;
    private DirAndFileStructure dir;
    private ParameterBank param = null;
    private BufferedImage latticeImg;
    private int ImageTextOffset = 90;
    private Graphics2D g;
    private int scale;
    private int nfixed = 0;
    private int stable = -1;
    private String fixedFname;
    private int R = 1;
    private int granularity = 0; // divide lattice into sublattices
    private AtomicIntegerArray subLatticeSum;
    private AtomicIntegerArray fixedSubLatticeSum;
    private double lastMagSubAvg = 0;
    private boolean staggeredSpinDraw = true;											
    private boolean initialized = false;
    private boolean useHeter;
    private boolean makingVideo = false;											//lattice else only sum is maintained

    /**
    * AtomicLattice constructor.
    *
    * @param s - spin to initialize the lattice to
    */
    public AtomicLattice(int s) {
        this(s, "", "", 1);
    }

    /**
    * AtomicLattice constructor.
    *
    * @param s - spin to initialize the lattice to
    * @param id - instance id of currently running instance of program
    */
    public AtomicLattice(int s, int id) {
        this(s, "", "", id);

    }

    /**
    * AtomicLattice constructor.
    *
    * @param s - spin to initialize the lattice to
    * @param postfix - postfix for the parameter file name
    */
    public AtomicLattice(int s, String postfix) {
        this(s, postfix, "", 1);
    }

    /**
    * AtomicLattice constructor.
    *
    * @param s - spin to initialize the lattice to
    * @param postfix - postfix for the parameter file name
    * @param id - instance id of currently running instance of program
    */
    public AtomicLattice(int s, String postfix, int id) {
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
    public AtomicLattice(int s, String postfix, String fnamepost, int id) {
        param = new ParameterBank(postfix);
        instId = id;
        dir = new DirAndFileStructure();
        L = param.L;
        D = param.D;
        // stable opposite of initialized, assumption
        stable = -s;
        Geo = param.Geo;
        granularity = param.Granularity;
        if (granularity == L) {
            granularity = 1;
        }
        lattice = null;
        N = 0;
        R = param.R;
        if (param.useLongRange == false) {
            R = 0;
        }
        magnetization = new AtomicInteger(0);
        magStaggered = new AtomicInteger(0);

        //Rectangular lattice. Would implement single index but mod 
        //would cost computationally
        if (Geo == 2) {
            N = (int) Math.pow(L, D);

            if (D == 2) {

                // Initialize sublattices if needed
                if (granularity > 1) {
                    int subL = (int) L / granularity;
                    subLatticeSum = new AtomicIntegerArray(subL * subL);
                    fixedSubLatticeSum = new AtomicIntegerArray(subL * subL);
                    maxFixedSum = new AtomicInteger(0);
                    for (int u = 0; u < (subL * subL); u++) {
                        subLatticeSum.set(u, 0);
                    }
                }
                lattice = new AtomicIntegerArray(L * L);
            }
            if (D == 3) {

                // Initialize sublattices if needed
                if (granularity > 1) {
                    int subL = (int) L / granularity;
                    subLatticeSum = new AtomicIntegerArray(subL * subL);
                    fixedSubLatticeSum = new AtomicIntegerArray(subL * subL);
                    maxFixedSum = new AtomicInteger(0);
                    for (int u = 0; u < (subL * subL); u++) {
                        subLatticeSum.set(u, 0);
                    }
                }
                lattice = new AtomicIntegerArray(L * L * L);
            }
        }

        // Triangular/Honeycomb lattice
        if (Geo == 3) {
        }

        // Initialize spins to spin state s
        initialize(s);

        // If Heterogenous- Fix some spins
        useHeter = param.useHeter;
        if (useHeter) {
            readFixedFile(fnamepost);
        }

        // Assert that magnetization is recalculated fresh
        magnetization.set(calcMagnetization2D()[0]);
        magStaggered.set(calcMagnetization2D()[1]);
    
        // Initial draw of lattice image
        spinDraw();

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
        if (this.getValue(i, j) == s) {
            return;
        }
        // Magnetization changes differently in diluted site
        if (granularity > 1) {
            updateSublatticeSum(i, j, k, s);
        }
        // Magnetization changes differently in diluted site
        if (getValue(i,j,k) != 0) {
            magnetization.set(magnetization.get() + 2 * s);
            if((i+j)%2 == 0) { 
                magStaggered.set(magStaggered.get()+2*s);
            }else{ 
                magStaggered.set(magStaggered.get()-2*s);
            }
        } else {
            magnetization.set(magnetization.get() + s);
            if((i+j)%2 == 0) { 
                magStaggered.set(magStaggered.get()+s);
            }else{ 
                magStaggered.set(magStaggered.get()-s);
            }
        }

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
        setValue(i%L,(int)((double)i/(double)L)%L,(int)((double)i/(double)(L*L))%L,s,t);
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
        if (getValue(i, j) == s) {
            return;
        }

        if (granularity > 1) {
            updateSublatticeSum(i, j, k, s);
        }
        // Magnetization changes differently in diluted site
        if (getValue(i,j,k) != 0) {
            magnetization.set(magnetization.get() + 2 * s);
            if((i+j+k)%2 == 0) { 
                magStaggered.set(magStaggered.get()+2*s);
            }else{ 
                magStaggered.set(magStaggered.get()-2*s);
            }
        } else {
            magnetization.set(magnetization.get() + s);
            if((i+j+k)%2 == 0) { 
                magStaggered.set(magStaggered.get()+s);
            }else{ 
                magStaggered.set(magStaggered.get()-s);
            }
        }
        
        // Change the spin at site    
        lattice.set(i + j * L + k * L * L, s);

        // If doing video update the image of lattice, kill time output.
        if (initialized && makingVideo) {
            updateImg(i, j, s, 0);
        }
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
        return magnetization.get();
    }

    /**
    *  getMagStaggered gives the staggered magnetization of the lattice
    */
    public int getMagStaggered() {
        return magStaggered.get();
    }
    
    /**
    * getSystemImg should return an image of the lattice
    */
    @Override
    public BufferedImage getSystemImg() {
        if (makingVideo == false) {
            spinDraw();
        }
        return latticeImg;
    }

    /**
    * makeVideo should set the lattice to begin making video which should
    * require updating of the lattice image.
    */
    @Override
    public void makeVideo() {
        makingVideo = true;
        spinDraw();
    }

    /**
    * saveLatticeImage saves the image of the lattice with the given name
    *
    * @param name - filename
    */
    public void saveLatticeImage(String name) {
        DataSaver dat = new DataSaver(instId);
        dat.saveImage(latticeImg, "png", name);
    }

    /**
    * updateSublatticeSum updates the sublattices with the new spin value
    */
    private void updateSublatticeSum(int i, int j, int k, int s) {
        int subI = (int) i / granularity;
        int subJ = (int) j / granularity;
        int subL = (int) L / granularity;
        int subK = 0;
        if (D == 3) {
            subK = (int) k / granularity;
        }

        int subIndex = subI + subJ * subL + subK * subL * subL;

        // Have to initialize the sublattice sums before adding in 
        if (!initialized) {
            subLatticeSum.set(subIndex, subLatticeSum.get(subIndex) + s);
        } else {
            // Magnetization and sum change differently for diluted case
            if (this.getValue(i, j) != 0) {
                subLatticeSum.set(subIndex, subLatticeSum.get(subIndex) + 2 * s);
            } else {
                subLatticeSum.set(subIndex, subLatticeSum.get(subIndex) + s);
            }
        }
    }

    /**
    * initializeFixedSubSum initializes the sublattices to zero
    */
    private void initializeFixedSubSum() {
        int subL = (int) L / granularity;

        for (int i = 0; i < (subL * subL); i++) {
            fixedSubLatticeSum.set(i, 0);
        }
    }

    @Override
    /**
    * setFixedLatticeValues should set all fixed values in the lattice.
    */
    @SuppressWarnings("CallToThreadDumpStack")
    public synchronized void setFixedLatticeValues() {
        if (granularity > 1) {
            initializeFixedSubSum();
        }
        String fileName;
        fileName = fixedFname;
        int fsum = 0;

        fileName = dir.getFixedDirectory() + fileName;
        //Scanner scanner = new Scanner(fileName);
        Scanner scanner;
        int sum = 0;
        nfixed = 0;

        int subI = 0;
        int subJ = 0;
        int subL = (int) L / granularity;
        int subK = 0;
        int subIndex;
        int maxSum = 0;
        int curr;

        try {
            scanner = new Scanner(new File(fileName));

            while (scanner.hasNextInt()) {
                int x = scanner.nextInt();
                int y = scanner.nextInt();
                int z = 0;

                if (granularity > 1) {
                    //fixed sub sum
                    subI = (int) x / granularity;
                    subJ = (int) y / granularity;
                    subL = (int) L / granularity;
                }

                if (D == 3) {
                    z = scanner.nextInt();
                    if (granularity > 1) {
                        subK = (int) z / granularity;
                    }
                }

                int fSpin = scanner.nextInt();
                fsum += fSpin;
                nfixed++;

                // sum used to fix magnetization
                if (getValue(x, y, z) == (-1 * fSpin)) {
                    sum += 2 * fSpin;
                } else if (getValue(x, y, z) == 0) {
                    sum += fSpin;
                }

                setValue(x, y, z, fSpin);

                if (granularity > 1) {
                    subIndex = subI + subJ * subL + subK * subL * subL;
                    curr = fixedSubLatticeSum.get(subIndex);
                    fixedSubLatticeSum.set(subIndex, (curr + 1));
                    curr = fixedSubLatticeSum.get(subIndex);
                    if (curr > maxSum) {
                        maxSum = curr;
                    }
                }
                fixed.set(x + y * L + z * L * L, 1);
            }

        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

        if (nfixed > 0) {
            fixSpinVal = (int) (fsum / nfixed);
        }

        if (granularity > 1) {
            maxFixedSum.set(maxSum);
        }
    }

    private int[] calcMagnetization2D() {
        // assert starting fresh
        int[] sum = new int[2];
        for (int u = 0; u < L; u++) {
            for (int v = 0; v < L; v++) {
                sum[0] += getValue(u, v);
                if((u+v)%2==0) { 
                    sum[1] += getValue(u, v);
                }else{ 
                    sum[1] -= getValue(u, v);
                }
            }
        }
        return sum;
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
        int sumStag = 0;
        makingVideo = false;

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
                            if((i+j+k)%2 == 0){
                                sumStag += spin;
                            }else{
                                sumStag -= spin;
                            }
                        }
                    } else {
                        setValue(i, j, 0, spin);
                        sum += spin;
                        if((i+j)%2 == 0){
                            sumStag += spin;
                        }else{
                            sumStag -= spin;
                        }
                    }
                }
            }
        }
        // Add for non square geometries

        // Create image layer
        initializeImg();

        magnetization.set(sum);
        magStaggered.set(sumStag);

    }

    /**
    * initializeImg creates the image layer of the lattice
    *
    */
    private void initializeImg() {

        // Standard scale of images
        scale = 5;
        int sizeImage = scale * L;
        while (sizeImage < 600) {
            scale = scale * 2;
            sizeImage = scale * L;
        }

        // Make the image
        latticeImg = new BufferedImage(sizeImage, sizeImage + ImageTextOffset, BufferedImage.TYPE_INT_RGB);
        g = latticeImg.createGraphics();
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
        //erase old text
        g.setColor(Color.BLACK);
        g.fillRect(((int) L * scale * 6 / 10), 0, ((int) L * scale * 4 / 10), ImageTextOffset);
        g.setColor(Color.WHITE);

        if (t != 0) {
            Font font = new Font("Courier New", Font.BOLD, 72);
            g.setFont(font);
        }
        // 70 is the x coordinate of text
        g.drawString("t = " + t, ((int) L * scale * 6 / 10), ((int) ImageTextOffset * 7 / 10));

        // Draw the update
        spin = getValue(i, j); 
        if(staggeredSpinDraw && (param.jInteraction < 0) && R > 0 && (((i+j)%2) == 0) ){
            spin *= -1;
        }
        if (spin == 1) {
            g.setColor(Color.WHITE);
        } else if (spin == (-1)) {
            g.setColor(Color.BLACK);
        } else if (spin == (0)) {
            g.setColor(Color.YELLOW);
        }
        g.fillRect(i * scale, j * scale + ImageTextOffset, scale, scale);
    }

    /**
    * spinDraw completely redraws the image layer based on the current
    * configuration of the lattice.
    */
    private void spinDraw() {
        int spin;
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < L; j++) {

                if (useHeter && isThisFixed(i, j, 0)) {
                    // Paint different color for fixed spins
                    // Cyan is like a blueish
                    // Magenta is purpleish
                    spin = getValue(i, j);
                    if(staggeredSpinDraw && (param.jInteraction < 0) && (((i+j)%2) == 0) ){
                        spin *= -1;
                    }
                    if (spin == 1) {
                        g.setColor(Color.BLUE);
                    }
                    if (spin == (-1)) {
                        g.setColor(Color.RED);
                    }
                    if (spin == (0)) {
                        g.setColor(Color.MAGENTA);
                    }
                } else {
                    spin = getValue(i, j); 
                    if(staggeredSpinDraw && (param.jInteraction < 0) && (((i+j)%2) == 0) ){
                        spin *= -1;
                    }
                    if (spin == 1) {
                        g.setColor(Color.WHITE);
                    }
                    if (spin == (-1)) {
                        g.setColor(Color.BLACK);
                    }
                    if (spin == (0)) {
                        g.setColor(Color.YELLOW);
                    }
                }

                g.fillRect(i * scale, j * scale + ImageTextOffset, scale, scale);
            }
        }
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
        if (fixed.get(i + j * L + k * L * L) == 1) {
            return true;
        } else {
            return false;
        }
    }

    @Override
    /**
    * getNeighSum should return the sum of the spins within the interaction
    * range centered by the coordinates given. 3d case
    *
    * @param i - i coordinate
    * @param j - j coordinate
    * @param k - k coordinate
    */
    public int getNeighSum(int i, int j, int k) {
        int neighSum = 0;

        if (R < 1) {

            // Get neighbors deal with periodic conditions, Calculate Sum while getting values
            neighSum += lattice.get(((i + 1) % L) + j * L + k * L * L);
            neighSum += lattice.get((i + ((j + 1) % L) * L + k * L * L));
            neighSum += lattice.get(((i + L - 1) % L) + j * L + k * L * L);
            neighSum += lattice.get(i + ((j + L - 1) % L) * L + k * L * L);
            // 3d neighbors
            if (D == 3) {
                neighSum += lattice.get(i + j * L + ((k + 1) % L) * L * L);
                neighSum += lattice.get(i + j * L + ((k + L - 1) % L) * L * L);
            }

            // 3D Case of square range
        } else if (D == 3) {

            int u;
            int v;
            int z;

            for (int m = 0; m < (2 * R + 1); m++) {
                for (int n = 0; n < (2 * R + 1); n++) {
                    for (int p = 0; p < (2 * R + 1); p++) {
                        u = ((i - R + m + L) % L);
                        v = ((j - R + n + L) % L);
                        z = ((k - R + p + L) % L);

                        if (u == i && v == j && z == k) {
                        } else {
                            neighSum += lattice.get(u + v * L + z * L * L);
                        }
                    }
                }
            }

            // 2D Case of square long range	
        } else {

            int u;
            int v;

            for (int m = 0; m < (2 * R + 1); m++) {
                for (int n = 0; n < (2 * R + 1); n++) {
                    u = ((i - R + m + L) % L);
                    v = ((j - R + n + L) % L);

                    if (u == i && v == j) {
                    } else {
                        neighSum += lattice.get(u + v * L);
                    }
                }
            }
        }

        return neighSum;
    }

    @Override
    /**
    * getNeighSum should return the sum of the spins within the interaction
    * range centered by the coordinates given. 2d case
    *
    * @param i - i coordinate
    * @param j - j coordinate
    */
    public int getNeighSum(int i, int j) {
        return getNeighSum(i, j, 0);
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
                int sumStag = 0;
                while (scanner.hasNextInt()) {
                    if (D == 2) {
                        int x = scanner.nextInt();
                        int y = scanner.nextInt();
                        int fSpin = scanner.nextInt();
                        setValue(x, y, 0, fSpin);
                        sum = fSpin + sum;
                        if((x+y)%2 == 0){
                            sumStag += fSpin;
                        }else{
                            sumStag -= fSpin;
                        }
                    } else {
                        int x = scanner.nextInt();
                        int y = scanner.nextInt();
                        int z = scanner.nextInt();
                        int fSpin = scanner.nextInt();
                        setValue(x, y, z, fSpin);
                        sum = fSpin + sum;
                        if((x+y+z)%2 == 0){
                            sumStag += fSpin;
                        }else{
                            sumStag -= fSpin;
                        }
                    }

                }
                magnetization.set(sum);
                magStaggered.set(sumStag);
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
    * getFixedSubLatSum should return the sum of the fixed spin amounts in
    * sublattice in the index given
    *
    * @param j - index of sublattice
    */
    @Override
    public int getFixedSubLatSum(int j) {
        return fixedSubLatticeSum.get(j);
    }

    /**
    * getFixSpinVal should return the average value of the fixed spins
    *
    */
    @Override
    public int getFixSpinVal() {
        return fixSpinVal;
    }

    /**
    * getHighFixIndex returns the index of the highest density of fixed spins
    */
    @Override
    public int getHighFixIndex() {

        // Find max fixed sum area
        int iMax = 0;
        int maxM = 0;
        int subL = (int) L / granularity;

        for (int i = 0; i < (subL * subL); i++) {
            if (fixedSubLatticeSum.get(i) > maxM) {
                iMax = i;
            }
        }

        /*
         * System.out.println("AtomicLattice | Highest possible index : "+(subL*subL)); int subX
         * = iMax%subL; int subY = (int)((iMax-subX)/subL);
         * System.out.println("AtomicLattice | i: "+subX+" j: "+subY);
         */
        return iMax;
    }

    @Override
    /**
    * highFixedAndStableMatched should return if the sublattice index for the
    * highest density stable spins and fixed spins match.
    */
    public boolean highFixedAndStableMatched() {
        // Do not bother if not granular
        if (granularity == 1) {
            return false;
        } else {
            // Find max stable sum area
            int iMax = 0;
            int maxM = 0;
            int subL = (int) L / granularity;

            for (int i = 0; i < (subL * subL); i++) {
                if ((subLatticeSum.get(i) * stable) > maxM) {
                    iMax = i;
                }
            }

            //iMax =iMax-subL;
            //System.out.println("AtomicLattice | index maxSumM: "+iMax+"   index fixMax: "+highFixedIndex());

            // if area of max density is equal to max stable density return true
            if (fixedSubLatticeSum.get(iMax) == maxFixedSum.get()) {
                return true;
            }
        }
        return false;
    }

    /**
    * getHighMagIndex returns the index of the highest density of stable spins
    */
    @Override
    public int getHighMagIndex() {
        // Find max stable sum area
        int iMax = 0;
        double maxM = 0;
        double sumM = 0.0;
        double effM;
        int subL = (int) L / granularity;

        for (int i = 0; i < (subL * subL); i++) {
            // Calculate effective magnetization ie dont count the fixed spins
            effM = (subLatticeSum.get(i) - getFixedSubLatSum(i) * getFixSpinVal()) / (subL * subL - getFixedSubLatSum(i));
            if ((effM * stable) > maxM) {
                iMax = i;
            }
            sumM += effM;
        }
        sumM = sumM / (subL * subL);
        lastMagSubAvg = sumM;
        return iMax;
    }

    /**
    * getLastAvgMagSubLat should return the average magnetization of the
    * sublattices.
    */
    @Override
    public double getLastAvgMagSubLat() {
        return lastMagSubAvg;
    }

    /**
    * getSubLatSum should return the sum of the sublattice in the index given
    *
    * @param u - index of sublattice
    */
    @Override
    public int getSubLatSum(int u) {
        return subLatticeSum.get(u);
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

    // test the class
    public static void main(String[] args) {
        AtomicLattice lat = new AtomicLattice(1);
        System.out.println("AtomicLattice |  Sum of neighbors " + lat.getNeighSum(3, 4));
        //System.out.println("AtomicLattice | index high fixed: " + lat.getHighFixIndex());
        System.out.println("AtomicLattice | __________________________________");
        System.out.println(lat.magnetization);
        lat.saveLatticeImage("");
        lat.setValue(3, 4, 0, -1);
        lat.setValue(3, 5, 0, -1);
        lat.setValue(3, 6, 0, -1);
        lat.setValue(4, 4, 0, -1);
        System.out.println("AtomicLattice |  Sum of neighbors " + lat.getNeighSum(0, 0));
        lat.saveLatticeImage("post");
        System.out.println("AtomicLattice | N in range: " + lat.getNinRange() + "  with R : " + lat.R);
        System.out.println("AtomicLattice | Done!");
    }
}

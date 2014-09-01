package Backbone.System;

/**
*       @(#)     AtomicLatticeSumSpin
*/

import Backbone.Util.*;
import Backbone.Visualization.*;
import com.googlecode.javacv.CanvasFrame;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicIntegerArray;

/**
*       Atomic LatticeMagInt class implements lattice as AtomicIntegerArray. It also
* implements functions to determine sector of highest stable spin density and
* for fixed spin systems the highest fixed spin density. <br>
*
*
* @author James Silva <jbsilva @ bu.edu>
* @since 2012-01
*/
public final class AtomicLatticeSumSpin implements LatticeMagInt {
    private int fixSpinVal; private String fixedFname;
    private int L = 256;private int N;
    private int R = 1;private int D = 2;          // dimension only up to 3
    private int Geo = 2;       // geometry - 2 is rectangular
    private boolean allowDilution = false;
    private int instId = 1;
    private int sumNeigh;
    private AtomicIntegerArray lattice;   // lattice want to implement using just one index
    private AtomicIntegerArray latticeSum; 
    private AtomicInteger magnetization; private AtomicInteger magStaggered;
    private AtomicInteger maxFixedSum; private AtomicIntegerArray fixed;
    private DirAndFileStructure dir; private ParameterBank param = null;
    private int nfixed = 0; private int nInteract = 4;
    private int stable = -1;
    private TriangularLattice2D triangLatVis; private SquareLattice2D sqLatVis;
    private HoneycombLattice2D honeyLatVis;
    private int granularity = 0; // divide lattice into sublattices
    private AtomicIntegerArray subLatticeSum; private AtomicIntegerArray fixedSubLatticeSum;
    private AtomicIntegerArray neighborOffsets; private AtomicIntegerArray neighborOffsets2;
    private AtomicIntegerArray neighbors;
    private boolean useHeter;
    private double lastMagSubAvg = 0;
    private boolean initialized = false;
    private boolean staggeredSpinDraw = true;											
    private boolean makingVideo = false;
    private boolean circularRange = false;
    private CanvasFrame canvasFrame ;
    
    //lattice else only sum is maintained

    /**
    *       AtomicLatticeSumSpin constructor.
    *
    * @param s - spin to initialize the lattice to
    */
    public AtomicLatticeSumSpin(int s){ this(s, "", "", 1); }

    /**
    *       AtomicLatticeSumSpin constructor.
    *
    * @param s - spin to initialize the lattice to
    * @param id - instance id of currently running instance of program
    */
    public AtomicLatticeSumSpin(int s, int id) {this(s, "", "", id); }

    /**
    *       AtomicLatticeSumSpin constructor.
    *
    * @param s - spin to initialize the lattice to
    * @param postfix - postfix for the parameter file name
    */
    public AtomicLatticeSumSpin(int s, String postfix) { this(s, postfix, "", 1); }

    /**
    *       AtomicLatticeSumSpin constructor.
    *
    * @param s - spin to initialize the lattice to
    * @param postfix - postfix for the parameter file name
    * @param id - instance id of currently running instance of program
    */
    public AtomicLatticeSumSpin(int s, String postfix, int id) { this(s, postfix, "", id); }

    /**
    *       AtomicLatticeSumSpin constructor.
    *
    * @param s - spin to initialize the lattice to
    * @param postfix - postfix for the parameter file name
    * @param fnamepost - postfix for the fixed spin configuration file
    * @param id - instance id of currently running instance of program
    */
    public AtomicLatticeSumSpin(int s, String postfix, String fnamepost, int id) {
        param = new ParameterBank(postfix); dir = new DirAndFileStructure();
        instId = id;
        L = param.L; D = param.D;
        // stable opposite of initialized, assumption
        stable = -s;
        Geo = param.Geo;
        useHeter = param.useHeter;
        granularity = param.Granularity;
        if( granularity == L ){ granularity = 1; }
        lattice = null;
        N = 0; R = param.R;
        if( param.useLongRange == false ){ R = 0; }
        magnetization = new AtomicInteger(0); magStaggered = new AtomicInteger(0);
        
        //Rectangular lattice. Would implement single index but mod 
        //would cost computationally
        if (Geo == 2 || Geo == 4) {
            N = (int) Math.pow(L, D);
            if( D == 2 ){
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
                lattice = new AtomicIntegerArray(L*L);
                latticeSum = new AtomicIntegerArray(L*L);
            }
            if( D == 3 ){
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
                latticeSum = new AtomicIntegerArray(L*L*L);
            }
        }else if( Geo == 3 || Geo == 6 ){
            // Triangular/Honeycomb lattice
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
                lattice = new AtomicIntegerArray(L*L);
                latticeSum = new AtomicIntegerArray(L*L);
            }
        }else{
            System.err.println("AtomicLatticeSumSpin | ERROR : "+Geo+" IS NOT A COMPATIBLE GEOMETRY");
        }

        // Initial draw of lattice image
        if( Geo == 2 || Geo  == 4 ){
            sqLatVis = new SquareLattice2D(L,R,param.jInteraction,useHeter,staggeredSpinDraw);
        }else if(Geo == 6){
            triangLatVis = new TriangularLattice2D(L,R,param.jInteraction,useHeter,staggeredSpinDraw);
            initNeighborOffsets();
        }else if(Geo == 3){
            honeyLatVis  = new HoneycombLattice2D(L,R,param.jInteraction,useHeter,staggeredSpinDraw);
            initNeighborOffsets();
        }
        
        // Initialize spins to spin state s
        initialize(s);

        // If Heterogenous- Fix some spins
        useHeter = param.useHeter;
        if( useHeter ){ readFixedFile(fnamepost); }

        spinDraw();
        // initialize spin sum
        initNeighSum();
        
        // Assert that magnetization is recalculated fresh
        magnetization.set(calcMagnetization2D()[0]); magStaggered.set(calcMagnetization2D()[1]);
        
        // done with initialization
        initialized = true;
    }

    /**
    *       getValue which gets the spin value in the lattice. 3d case
    *
    * @param i - i coordinate
    * @param j - j coordinate
    * @param k - k coordinate
    */
    @Override
    public int getValue(int i, int j, int k) { return lattice.get(i + j * L + k * L * L); }
    
    /**
    * getValue which gets the spin value in the lattice.
    *
    * @param i - coordinate
    */
    public int getValue(int i) { return lattice.get(i); }
    
    /**
    *       setValue updates the spin value in the lattice and the magnetization with
    * update of time in image layer .
    *
    * @param i - i coordinate
    * @param j - j coordinate
    * @param k - k coordinate
    * @param s - new spin value
    * @param t - time
    */
    @Override
    public void setValue(int i, int j, int k, int s, int t) {
        if (getValue(i, j) == s) {
            return;
        }
        // Magnetization changes differently in diluted site
        if (granularity > 1) {
            updateSublatticeSum(i, j, k, s);
        }
        // Magnetization changes differently in diluted site
        if (getValue(i,j,k) != 0) {
            magnetization.set(magnetization.get()+2*s);
            if((Geo == 3 && (i%2 == 0)) || (Geo !=3 && (i+j)%2 == 0)){
                magStaggered.set(magStaggered.get()+2*s);
            }else{ 
                magStaggered.set(magStaggered.get()-2*s);
            }
        } else {
            magnetization.set(magnetization.get()+s);
            if((Geo == 3 && (i%2 == 0)) || (Geo !=3 && (i+j)%2 == 0)){
                magStaggered.set(magStaggered.get()+s);
            }else{ 
                magStaggered.set(magStaggered.get()-s);
            }
        }

        if(Geo == 2 || Geo  == 4 ){
            updateNeighSumSqLat(i, j, k, s, (getValue(i,j,k) == 0) ? true:false);
        }else if(Geo == 3){
            updateNeighSumHoneyLat(i, j, s, (getValue(i,j,k) == 0) ? true:false,R,-1);
        }else if(Geo == 6){
            updateNeighSumTriangleLat(i, j, s, (getValue(i,j,k) == 0) ? true:false,R,-1);
        }
        
        // Change the spin at site    
        lattice.set(i + j * L + k * L * L, s);
        // If doing video update the image of lattice including the time
        if (initialized && makingVideo) {
            updateImg(i, j, s, t);
        }
    }

    public void updateImg(int i, int j, int s, int t){
        if(Geo == 2 || Geo  == 4){
            sqLatVis.updateImg(i, j, s, t);
        }else if(Geo == 6){
            triangLatVis.updateImg(i, j, s, t);
        }else if(Geo == 3){
            honeyLatVis.updateImg(i, j, s, t);
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
    public void setValue(int i, int s, int t) {
        setValue(i%L,(int)((double)i/(double)L)%L,(int)((double)i/(double)(L*L))%L,s,t);
    }
    
    /**
    *       setValue which updates the spin value in the lattice and the
    * magnetization
    *
    *
    * @param i - i coordinate
    * @param j - j coordinate
    * @param k - k coordinate
    * @param s - new spin value
    */
    public void setValue(int i, int j, int k, int s) {
        if( getValue(i, j) == s ){ return; }
        if( granularity > 1 ){ updateSublatticeSum(i, j, k, s); }
        // Magnetization changes differently in diluted site
        if( getValue(i,j,k) != 0 ){
            magnetization.set(magnetization.get() + 2 * s);
            if((Geo == 3 && (i%2 == 0)) || (Geo !=3 && (i+j)%2 == 0)){
                magStaggered.set(magStaggered.get()+2*s);
            }else{ 
                magStaggered.set(magStaggered.get()-2*s);
            }
        } else {
            magnetization.set(magnetization.get() + s);
            if((Geo == 3 && (i%2 == 0)) || (Geo !=3 && (i+j)%2 == 0)){
                magStaggered.set(magStaggered.get()+s);
            }else{ 
                magStaggered.set(magStaggered.get()-s);
            }
        }
        
        updateNeighSum(i,j,k,s);
    
        // Change the spin at site    
        lattice.set(i + j * L + k * L * L, s);

        // If doing video update the image of lattice, kill time output.
        if( initialized && makingVideo ){
            updateImg(i, j, s, 0);
        }
    }
    
    private void updateNeighSum(int i, int j, int k, int s){
        if( Geo == 2 || Geo == 4 ){
            updateNeighSumSqLat(i, j, k, s, (getValue(i,j,k) == 0) ? true:false);
        }else if( Geo == 3 ){
            updateNeighSumHoneyLat(i, j, s, (getValue(i,j,k) == 0) ? true:false,R,-1);
        }else if( Geo == 6 ){
            updateNeighSumTriangleLat(i, j, s, (getValue(i,j,k) == 0) ? true:false,R,-1);
        }
    }
    
    /**
    *       getValue which gets the spin value in the lattice. 2d case
    *
    * @param i - i coordinate
    * @param j - j coordinate
    */
    @Override
    public int getValue(int i, int j) { return lattice.get(i + j * L); }

    public int getRandom(int min, int max) { return (min + (int) (Math.random() * ((max - min) + 1))); }

    /**
    *       getLength which gets the length of the lattice
    */
    @Override
    public int getLength() { return L; }

    /**
    *       getDimension gives the dimensionality of the lattice
    */
    @Override
    public int getDimension() { return D; }

    /**
    *       getGeo gives the geometry of the lattice
    */
    @Override
    public int getGeo() {return Geo; }

    /**
    *       getN gives the size of the lattice
    */
    @Override
    public int getN() { return N; }

    /**
    *       getNFixed gives the amount of fixed spins in the lattice
    */
    @Override
    public int getNFixed() { return nfixed; }

    /**
    *       getRange gives the interaction range of the lattice
    */
    @Override
    public int getRange() { return R; }

    /**
    *       fixAsite randomly fixes a site
    */
    public void fixAsite(int spinVal){
        boolean oldSite = true;
        int site;
        while(oldSite){
            site = (int)(Math.random()*N);
            int i = site%L;
            int j = (int)((double)site/(double)L)%L;
            //System.out.println("Fixing spin | " +isThisFixed(i, j, 0)+"   i| "+i+"  j| "+j);
            if(!isThisFixed(i, j, 0)){
                oldSite = false;
                setValue(site , spinVal,0);
                fixed.set(site , 1);
                if(Geo == 6){
                    triangLatVis.fixAspin(i, j, spinVal);
                }else if(Geo == 3){
                    honeyLatVis.fixAspin(i, j, spinVal);
                }else{
                    sqLatVis.fixAspin(i, j, spinVal);
                }
            }
        }
        nfixed++;
    }
    
    /**
    *       getMagnetization gives the magnetization of the lattice
    */
    @Override
    public int getMagnetization() { return magnetization.get(); }

    /**
    *       getMagnetization gives the magnetization of the lattice
    */
    public int getMagStaggered() { return magStaggered.get(); }
    
    /**
    *       getSystemImg should return an image of the lattice
    */
    @Override
    public BufferedImage getSystemImg() {
        if(makingVideo == false){ spinDraw(); }
        return getVisImg();
    }
    
    private void spinDraw(){
        if( Geo == 2 || Geo  == 4 ){
            sqLatVis.spinDraw(this);
        }else if( Geo == 6 ){
            triangLatVis.spinDraw(this);
        }else if( Geo == 3 ){
            honeyLatVis.spinDraw(this);
        }
    }
    
    private BufferedImage getVisImg(){
        if( Geo == 2 || Geo  == 4 ){
            return sqLatVis.getImageOfVis();
        }else if( Geo == 6 ){
            return triangLatVis.getImageOfVis();
        }else if( Geo == 3 ){
            return honeyLatVis.getImageOfVis();
        }
        return null;
    }
    
    /**
    *       makeVideo should set the lattice to begin making video which should
    * require updating of the lattice image.
    */
    @Override
    public void makeVideo() { makingVideo = true; spinDraw(); }

    /**
    *       saveLatticeImage saves the image of the lattice with the given name
    *
    * @param name - filename
    */
    public void saveLatticeImage(String name) {
        DataSaver dat = new DataSaver(instId);
        dat.saveImage(getVisImg(), "png", name);
    }

    /**
    *     updateSublatticeSum updates the sublattices with the new spin value
    */
    private void updateSublatticeSum(int i, int j, int k, int s) {
        int subI = (int) i / granularity; int subJ = (int) j / granularity;
        int subL = (int) L / granularity; int subK = 0;
        if (D == 3) { subK = (int) k / granularity; }

        int subIndex = subI + subJ * subL + subK * subL * subL;

        // Have to initialize the sublattice sums before adding in 
        if (!initialized) {
            subLatticeSum.set(subIndex, subLatticeSum.get(subIndex) + s);
        } else {
            // Magnetization and sum change differently for diluted case
            if (getValue(i, j) != 0) {
                subLatticeSum.set(subIndex, subLatticeSum.get(subIndex) + 2 * s);
            } else {
                subLatticeSum.set(subIndex, subLatticeSum.get(subIndex) + s);
            }
        }
    }

    /**
    *       initializeFixedSubSum initializes the sublattices to zero
    */
    private void initializeFixedSubSum() {
        int subL = (int) L / granularity;
        for (int i = 0; i < (subL * subL); i++) { fixedSubLatticeSum.set(i, 0); }
    }

    @Override
    /**
    *       setFixedLatticeValues should set all fixed values in the lattice.
    */
    @SuppressWarnings("CallToThreadDumpStack")
    public synchronized void setFixedLatticeValues() {
        if (granularity > 1) { initializeFixedSubSum(); }
        String fileName; fileName = fixedFname; int fsum = 0;

        fileName = dir.getFixedDirectory() + fileName;
        Scanner scanner; int sum = 0; nfixed = 0;

        int subI = 0; int subJ = 0; int subL = (int) L / granularity;
        int subK = 0; int subIndex; int maxSum = 0; int curr;

        if( (new File(fileName)).exists() ){
            try {
                scanner = new Scanner( new File(fileName) );
                while( scanner.hasNextInt() ) {
                    int x = scanner.nextInt(); int y = scanner.nextInt(); int z = 0;
                    if ( granularity > 1 ) {
                        //fixed sub sum
                        subI = (int) x / granularity;
                        subJ = (int) y / granularity;
                        subL = (int) L / granularity;
                    }
                    if( D == 3 ){
                        z = scanner.nextInt();
                        if( granularity > 1 ){ subK = (int) z / granularity; }
                    }

                    int fSpin = scanner.nextInt();
                    fsum += fSpin;
                    nfixed++;

                    // sum used to fix magnetization
                    if ( getValue(x, y, z) == (-1 * fSpin) ) {
                        sum += 2 * fSpin;
                    } else if ( getValue(x, y, z) == 0 ) {
                        sum += fSpin;
                    }
                    
                    setValue(x, y, z, fSpin);
                    
                    if ( granularity > 1 ) {
                        subIndex = subI + subJ * subL + subK * subL * subL;
                        curr = fixedSubLatticeSum.get(subIndex);
                        fixedSubLatticeSum.set(subIndex, (curr + 1));
                        curr = fixedSubLatticeSum.get(subIndex);
                        if( curr > maxSum ){ maxSum = curr; }
                    }
                    fixed.set(x + y * L + z * L * L, 1);
                }
            } catch (FileNotFoundException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
        if (nfixed > 0) { fixSpinVal = (int) (fsum / nfixed); }
        // Assert that magnetization is recalculated fresh
        magnetization.set(calcMagnetization2D()[0]); magStaggered.set(calcMagnetization2D()[1]);
        if( granularity > 1 ){ maxFixedSum.set(maxSum); }
    }

    private int[] calcMagnetization2D() {
        // assert starting fresh
        int[] sum = new int[2];
        for (int u = 0; u < L; u++) {
            for (int v = 0; v < L; v++) {
                sum[0] += getValue(u, v);
                if((Geo == 3 && (u%2 == 0)) || (Geo !=3 && (u+v)%2 == 0)){
                    sum[1] += getValue(u, v);
                }else{ 
                    sum[1] -= getValue(u, v);
                }
            }
        }
        return sum;
    }

    /**
    *       initialize sets the lattice values to their initial values which is given
    * as input.
    *
    * @param s - Initial lattice value; -2 for random values
    */
    @Override
    public synchronized void initialize(int s) {
        int spin = s; int sum = 0;
        int sumStag = 0; makingVideo = false;

        // Square/Triangular lattice initialize all lattice elements to value s
        // unless s is -2 which is a flag for random initialization of spins
        if (Geo == 2 || Geo == 3 || Geo == 4 || Geo  == 6) {
            for (int i = 0; i < L; i++) {
                for (int j = 0; j < L; j++) {
                    if (s == (-2)) {
                        if (allowDilution) {
                            spin = getRandom(-1, 1);
                        } else {
                            spin = getRandom(0, 1);
                            if (spin == 0) { spin = -1; }
                        }
                    }
                    if (D == 3) {
                        for (int k = 0; k < L; k++) {
                            setValue(i, j, k, spin);
                            sum += spin;
                            if((Geo == 3 && (i%2 == 0)) || (Geo !=3 && (i+j)%2 == 0)){
                                sumStag += spin;
                            }else{
                                sumStag -= spin;
                            }
                        }
                    } else {
                        setValue(i, j, 0, spin);
                        sum += spin;
                        if((Geo == 3 && (i%2 == 0)) || (Geo !=3 && (i+j)%2 == 0)){
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
        magnetization.set(sum); magStaggered.set(sumStag);
    }

    private void initializeImg(){
        if( Geo == 2 || Geo  == 4 ){
            sqLatVis.initializeImg();
        }else if( Geo == 6 ){
            triangLatVis.initializeImg();
        }else if( Geo == 3 ){
            honeyLatVis.initializeImg();
        }
    }
    
    /**
    *       isThisFixed should return true if the lattice coordinate is fixed.
    *
    * @param i - i coordinate
    * @param j - j coordinate
    * @param k - k coordinate
    */
    @Override
    public boolean isThisFixed(int i, int j, int k) {
        return (fixed.get(i + j * L + k * L * L) == 1);
    }

    private void updateNeighSumSqLat(int i,int j,int k, int newSpin,boolean diluted){
        if (R < 1) {
            int curr = ((i + 1) % L) + j * L + k * L * L;
            // Get neighbors deal with periodic conditions, Calculate Sum while getting values
            if(!diluted){
                latticeSum.set(curr,latticeSum.get(curr)+2*newSpin);
            }else{
                latticeSum.set(curr,latticeSum.get(curr)+newSpin);
            }
            curr = (i + ((j + 1) % L) * L + k * L * L);
            if(!diluted){
                latticeSum.set(curr,latticeSum.get(curr)+2*newSpin);
            }else{
                latticeSum.set(curr,latticeSum.get(curr)+newSpin);
            }
            curr = ((i + L - 1) % L) + j * L + k * L * L;
            if(!diluted){
                latticeSum.set(curr,latticeSum.get(curr)+2*newSpin);
            }else{
                latticeSum.set(curr,latticeSum.get(curr)+newSpin);
            }
            curr = (i + ((j + L - 1) % L) * L + k * L * L);
            if(!diluted){
                latticeSum.set(curr,latticeSum.get(curr)+2*newSpin);
            }else{
                latticeSum.set(curr,latticeSum.get(curr)+newSpin);
            }
            // 3d neighbors
            if (D == 3) {
                curr = (i + j * L + ((k + 1) % L) * L * L);
                if(!diluted){
                    latticeSum.set(curr,latticeSum.get(curr)+2*newSpin);
                }else{
                    latticeSum.set(curr,latticeSum.get(curr)+newSpin);
                }
                curr = (i + j * L + ((k + L - 1) % L) * L * L);
                if(!diluted){
                    latticeSum.set(curr,latticeSum.get(curr)+2*newSpin);
                }else{
                    latticeSum.set(curr,latticeSum.get(curr)+newSpin);
                }
            }
        }else{
            int u; int v; int z;
            for (int m = 0; m < (2 * R + 1); m++) {
                for (int n = 0; n < (2 * R + 1); n++) {
                    u = ((i - R + m + L) % L); v = ((j - R + n + L) % L);
                    if(D ==3){
                        for (int p = 0; p < (2 * R + 1); p++) {
                            z = ((k - R + p + L) % L);
                            if (!(u == i && v == j && z == k)) {
                                if(!circularRange ||
                                        (circularRange && (Math.sqrt((R-m)*(R-m)+(R-n)*(R-n)+(R-p)*(R-p))) < R)  ){
                                    if(!diluted){
                                        latticeSum.set(u+v*L+z*L*L,latticeSum.get(u+v*L+z*L*L)+2*newSpin);
                                    }else{
                                        latticeSum.set(u+v*L+z*L*L,latticeSum.get(u+v*L+z*L*L)+newSpin);
                                    }
                                }
                            }
                        }
                    }else{
                        if (!(u == i && v == j)) {
                            if(!circularRange ||
                                        (circularRange && (Math.sqrt((R-m)*(R-m)+(R-n)*(R-n))) < R)  ){
                                if(!diluted){
                                    latticeSum.set(u+v*L,latticeSum.get(u+v*L)+2*newSpin);
                                }else{
                                    latticeSum.set(u+v*L,latticeSum.get(u+v*L)+newSpin);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    private void initNeighSum(){
        for(int x = 0;x <L;x++){ for(int y = 0;y <L;y++){
            if(Geo == 2 || Geo  == 4){
                if(D==3){
                    for(int z = 0;z <L;z++){
                        latticeSum.set((x+y*L+z*L*L), neighSumSqLat(x,y,z));
                    }
                }else{
                    latticeSum.set((x+y*L), neighSumSqLat(x,y,0));
                }
            }else if(Geo == 3){
                latticeSum.set((x+y*L), neighSumHoneyLat(x,y));
            }else if(Geo == 6){
                latticeSum.set((x+y*L), neighSumTriangleLat(x,y));
            }
        }}
    }
    
    private int neighSumSqLat(int i, int j, int k){
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
            int u; int v; int z;
            for (int m = 0; m < (2 * R + 1); m++) {
                for (int n = 0; n < (2 * R + 1); n++) {
                    for (int p = 0; p < (2 * R + 1); p++) {
                        u = ((i - R + m + L) % L); v = ((j - R + n + L) % L);
                        z = ((k - R + p + L) % L);
                        if ( !(u == i && v == j && z == k) ) {
                            if(!circularRange ||
                                        (circularRange && (Math.sqrt((R-m)*(R-m)+(R-n)*(R-m)+(R-p)*(R-p))) < R)  ){
                                neighSum += lattice.get(u + v * L + z * L * L);
                            }
                        }
                    }
                }
            }
            // 2D Case of square long range	
        } else {
            int u; int v;
            for (int m = 0; m < (2 * R + 1); m++) {
                for (int n = 0; n < (2 * R + 1); n++) {
                    u = ((i - R + m + L) % L); v = ((j - R + n + L) % L);
                    if(!(u == i && v == j) ){
                        if(!circularRange ||
                                        (circularRange && (Math.sqrt((R-m)*(R-m)+(R-n)*(R-n)) < R) )){    
                            neighSum += lattice.get(u + v * L);
                        }
                    }
                }
            }
        }
        return neighSum;
    }
    
    private int neighSumHoneyLat( int i, int j ){
        int neighSum = 0; int  u = 0; int v = 0;
        setHoneyNeighborInRange(i, j);
        for(int k = 0; k < neighbors.length();k++){
            int ind = neighbors.get(k);
            neighSum += (ind == (i+j*L)) ? 0 :lattice.get(ind);
        }
        return neighSum;
    }
        
    private int neighSumTriangleLat(int i, int j){
        int  u = 0; int v = 0; int neighSum = 0;
        setTriangleNeighborInRange(i, j);
        for(int k = 0; k < neighbors.length();k++){
            int ind = neighbors.get(k);
            neighSum += (ind == (i+j*L)) ? 0 :lattice.get(ind);
        }
        return neighSum;
    }
    
    private void updateNeighSumHoneyLat(int i, int j, int newSpin,boolean diluted, int range, int callingNeigh){
        setHoneyNeighborInRange(i, j); nInteract = 0;
        for(int k = 0; k < neighbors.length();k++){
            int ind = neighbors.get(k);
            if(ind != (i+j*L)){
                nInteract++;
                if(!diluted){
                    latticeSum.set(ind,latticeSum.get(ind)+2*newSpin);
                }else{
                    latticeSum.set(ind,latticeSum.get(ind)+newSpin);
                }
            }
        }
    }
    
    private void updateNeighSumTriangleLat(int i, int j, int newSpin,boolean diluted, int range, int callingNeigh){
        setTriangleNeighborInRange(i, j); nInteract = 0;
        for(int k = 0; k < neighbors.length();k++){
            int ind = neighbors.get(k);
            if(ind != (i+j*L)){
                nInteract++;
                if(!diluted){
                    latticeSum.set(ind,latticeSum.get(ind)+2*newSpin);
                }else{
                    latticeSum.set(ind,latticeSum.get(ind)+newSpin);
                }
            }
        }        
    }

    
    private synchronized void setHoneyNeighborInRange(int i, int j){
        int ind = 0; int x; int y;
        for(int u = 0; u < (neighborOffsets.length()/2); u++){
            if( j % 2 != 0 ){
                if( i % 2 == 0 ){
                    x = (((j+neighborOffsets.get(2*u+1)+L)%L) % 2 == 0)?
                        (i+2+neighborOffsets.get(2*u)+L)%L : (i+neighborOffsets.get(2*u)+L)%L;
                }else{
                    // L+1
                    x = (((j+neighborOffsets2.get(2*u+1)+L)%L) % 2 == 0)?
                        (i+neighborOffsets2.get(2*u)+L)%L : (i+neighborOffsets2.get(2*u)+L)%L;
                }
            }else{
                if( i % 2 == 0 ){
                    x = (((j+neighborOffsets.get(2*u+1)+L)%L) % 2 == 0)?
                        (i+neighborOffsets.get(2*u)+L)%L : (i+neighborOffsets.get(2*u)+L)%L;
                }else{
                    // L = +1 patch
                    x = (((j+neighborOffsets2.get(2*u+1)+L)%L) % 2 == 0)?
                        (i+neighborOffsets2.get(2*u)+L)%L : (i-2+neighborOffsets2.get(2*u)+L)%L;
                }
            }
            y = ((j+neighborOffsets.get(2*u+1)+L)%L);
            ind = x+y*L;
            neighbors.set(u,ind);
        }
    }
    
    private synchronized void setTriangleNeighborInRange(int i, int j){
        int ind = 0; int x; int y;
        for(int u = 0; u < (neighborOffsets.length()/2); u++){
            if( j % 2 != 0 ){
                x = (((j+neighborOffsets.get(2*u+1)+L)%L) % 2 == 0)?
                    (i+1+neighborOffsets.get(2*u)+L)%L : (i+neighborOffsets.get(2*u)+L)%L;
            }else{
                x = (i+neighborOffsets.get(2*u)+L)%L;
            }
            y = ((j+neighborOffsets.get(2*u+1)+L)%L);
            ind = x+y*L;
            neighbors.set(u,ind);
        }
    }
    
    public void verifyTriangleNeigh(int i, int j){
        ArrayList<Integer> neighReg = getTriangleNeighborInRange(i, j, R, -1);
        ArrayList<Integer> neighNew = new ArrayList<Integer>();
        int mutualVals = 0;
        for( int u = 0; u < neighbors.length(); u++ ){
            int ind = 0;boolean found = false;
            while( ind < neighReg.size() && !found ){
                if( neighReg.get(ind) == neighbors.get(u) ){
                    found = true;
                    mutualVals++;
                    neighReg.remove(ind);
                }
                ind++;
            }
        }
        System.out.println("Verify| i: "+i+"   j: "+j+"   Neighbors: "+neighbors.length()+"     mutual: "+mutualVals);        
    }
    
    private synchronized ArrayList<Integer> getTriangleNeighborInRange(int i, int j){
        ArrayList<Integer> neigh = new ArrayList<Integer>((neighborOffsets.length()/2));
        int ind = 0;int x;int y;   
        for(int u = 0; u < (neighborOffsets.length()/2); u++){
            if( j % 2 != 0 ){
                x = (((j+neighborOffsets.get(2*u+1)+L)%L) % 2 == 0)?
                    (i+1+neighborOffsets.get(2*u)+L)%L : (i+neighborOffsets.get(2*u)+L)%L;
            }else{
                x = (i+neighborOffsets.get(2*u)+L)%L;
            }
            y = ((j+neighborOffsets.get(2*u+1)+L)%L);
            ind = x+y*L;
            neigh.add(ind);
        }
        return neigh;
    }
 
    private synchronized ArrayList<Integer> getHoneyNeighborInRange(int i, int j){
        ArrayList<Integer> neigh = new ArrayList<Integer>((neighborOffsets.length()/2));
        int ind = 0;int x;int y;
        for(int u = 0; u < (neighborOffsets.length()/2); u++){
            if( j % 2 != 0 ){
                if( i % 2 == 0 ){
                    x = (((j+neighborOffsets.get(2*u+1)+L)%L) % 2 == 0)?
                        (i+2+neighborOffsets.get(2*u)+L)%L : (i+neighborOffsets.get(2*u)+L)%L;
                }else{
                    // L+1
                    x = (((j+neighborOffsets2.get(2*u+1)+L)%L) % 2 == 0)?
                        (i+neighborOffsets2.get(2*u)+L)%L : (i+neighborOffsets2.get(2*u)+L)%L;
                }
            }else{
                if( i % 2 == 0 ){
                    x = (((j+neighborOffsets.get(2*u+1)+L)%L) % 2 == 0)?
                        (i+neighborOffsets.get(2*u)+L)%L : (i+neighborOffsets.get(2*u)+L)%L;
                }else{
                    // L=+1 patch
                    x = (((j+neighborOffsets2.get(2*u+1)+L)%L) % 2 == 0)?
                        (i+neighborOffsets2.get(2*u)+L)%L : (i-2+neighborOffsets2.get(2*u)+L)%L;
                }
            }
            y = ((j+neighborOffsets.get(2*u+1)+L)%L);
            ind = x+y*L; neigh.add(ind);
        }
        return neigh;
    }
    
    private void initNeighborOffsets(){
        int center = (int)(L/2); ArrayList<Integer> neighs;
        if( Geo == 6 ){
            neighs = getTriangleNeighborInRange(center,center,R,-1);
        }else if( Geo == 3 ){
            center = (int)(L/2)+1;
            neighs = getHoneyNeighborInRange(center,center,R,-1);
            nInteract = (R < 2) ? neighs.size():(neighs.size()-1);            
            neighborOffsets2 = new AtomicIntegerArray(2*(nInteract));
            int x = 0;int y = 0; int i = 0;int ind = 0;
            for( int u = 0; u < neighs.size(); u++ ){
                i = neighs.get(u);
                x = i%L;y = ((int)((double)i/(double)L))%L;
                if( i != (center+center*L) ){
                    neighborOffsets2.set(2*ind, x-center);
                    neighborOffsets2.set(2*ind+1, y-center);
                    ind++;
                }
            }
            center = (int)(L/2);
            neighs = getHoneyNeighborInRange(center,center,R,-1);         
        }else{
            neighs = null;
        }
        
        nInteract = (R < 2) ? neighs.size():(neighs.size()-1);            
        neighborOffsets = new AtomicIntegerArray(2*(nInteract));
        int x = 0; int y = 0; int i = 0; int ind = 0;
        for( int u = 0; u < neighs.size(); u++ ){
            i = neighs.get(u);
            x = i%L;y = ((int)((double)i/(double)L))%L;
            if( i != (center+center*L) ){
                neighborOffsets.set(2*ind, x-center);
                neighborOffsets.set(2*ind+1, y-center);
                ind++;
            }
        }
        neighbors = new AtomicIntegerArray(nInteract);
    }
    
    private synchronized ArrayList<Integer> getTriangleNeighborInRange(int i, int j, int range, int callingNeigh){
        int  u = 0; int v = 0;
        ArrayList<Integer> neighLat = new ArrayList<Integer>();        
        if( range < 2 ){
            ArrayList<Integer> neighNew = addTriangleNeighbors(neighLat, i, j);                
            for(int m = 0; m < neighNew.size(); m++){
                neighLat.add(neighNew.get(m));
            }
        }else{
            ArrayList<Integer> neighAtRange = addTriangleNeighbors(neighLat, i, j);
            for(int m = 0; m < neighAtRange.size(); m++){
                neighLat.add(neighAtRange.get(m));
            }
            for(int k = 0; k < (range-1); k++){
                int toProc = neighAtRange.size();
                for(int w = 0; w < toProc; w++){
                    int ind =neighAtRange.get(0);
                    u = ind%L; v = ((int)((double)ind/(double)L))%L;
                    ArrayList<Integer> neighNew = addTriangleNeighbors(neighLat, u, v);                
                    if( neighAtRange.size() > 0 ){ neighAtRange.remove(0); }
                    for(int m = 0; m < neighNew.size(); m++){
                        neighAtRange.add(neighNew.get(m));
                        neighLat.add(neighNew.get(m));
                    }
                }
            }
        }
        return neighLat;
    }
    
    private ArrayList<Integer> getHoneyNeighborInRange(int i, int j, int range, int callingNeigh){
        int  u = 0; int v = 0;
        ArrayList<Integer> neighLat = new ArrayList<Integer>();        
        if( range < 2 ){
            ArrayList<Integer> neighNew = addHoneyNeighbors(neighLat, i, j);                
            for(int m = 0; m < neighNew.size(); m++){ neighLat.add(neighNew.get(m)); }
        }else{
            ArrayList<Integer> neighAtRange = addHoneyNeighbors(neighLat, i, j);
            for(int m = 0; m < neighAtRange.size(); m++){ neighLat.add(neighAtRange.get(m)); }
            for(int k = 0; k < (range-1); k++){
                int toProc = neighAtRange.size();
                for(int w = 0; w < toProc; w++){
                    int ind =neighAtRange.get(0);
                    u = ind%L; v = ((int)((double)ind/(double)L))%L;
                    ArrayList<Integer> neighNew = addHoneyNeighbors(neighLat, u, v);                
                    if(neighAtRange.size() > 0){ neighAtRange.remove(0); }
                    for(int m = 0; m < neighNew.size(); m++){
                        neighAtRange.add(neighNew.get(m));
                        neighLat.add(neighNew.get(m));
                    }
                }
            }
        }
        return neighLat;
    }
    
    private ArrayList<Integer> addHoneyNeighbors(ArrayList<Integer> neigh , int sx, int sy){
        int u = 0; int v = 0;
        ArrayList<Integer> neighLat = new ArrayList<Integer>();        
        u = ((sx % 2) == 0) ? (sx+1)%L : (sx+L-1)%L ; v = sy;
        if(neighborsVisited(neigh,u+v*L)){ neighLat.add(u+v*L); }
        u = ((sy % 2) == 0) ? (sx+L-1)%L : (sx+1)%L ; v = (sy+1)%L;
        if(neighborsVisited(neigh,u+v*L)){ neighLat.add(u+v*L); }
        u =  ((sy % 2) == 0) ? (sx+L-1)%L : (sx+1)%L ; v = (sy+L-1)%L;
        if(neighborsVisited(neigh,u+v*L)){ neighLat.add(u+v*L); }
        return neighLat;
    }
    
    private synchronized  boolean neighborsVisited(ArrayList<Integer> neighbors , int neigh){
        int ind = 0; boolean newNeigh = true;
        while( ind < neighbors.size() && newNeigh ){
            if(neighbors.get(ind) == neigh){newNeigh = false;}
            ind++;
        }
        return newNeigh;
    }
    
    private synchronized  ArrayList<Integer> addTriangleNeighbors(ArrayList<Integer> neigh , int sx, int sy){
        int u = 0; int v = 0;
        ArrayList<Integer> neighLat = new ArrayList<Integer>();        
        
        u = (sx+L-1)%L ; v = sy;
        if(neighborsVisited(neigh,u+v*L)){ neighLat.add(u+v*L); }
        u = (sx+1)%L; v = sy;
        if(neighborsVisited(neigh,u+v*L)){ neighLat.add(u+v*L); }
        u = sx; v = (sy+L-1)%L;
        if(neighborsVisited(neigh,u+v*L)){ neighLat.add(u+v*L); }
        u = sx; v = (sy+1)%L;
        if(neighborsVisited(neigh,u+v*L)){ neighLat.add(u+v*L); }
        
        if((sy % 2) == 0){
            u =  (sx+L-1)%L ; v = (sy+1)%L;
            if(neighborsVisited(neigh,u+v*L)){ neighLat.add(u+v*L); }
            u = (sx+L-1)%L; v = (sy+L-1)%L;
            if(neighborsVisited(neigh,u+v*L)){ neighLat.add(u+v*L); }
        }else{
            u =  (sx+1)%L ; v = (sy+L-1)%L;
            if(neighborsVisited(neigh,u+v*L)){ neighLat.add(u+v*L); }
            u =  (sx+1)%L ; v = (sy+1)%L;
            if(neighborsVisited(neigh,u+v*L)){ neighLat.add(u+v*L); }
        }       
        return neighLat;
    }
    
    public void showNeighbors(int i, int j, int range, int callingNeigh){
        if(Geo == 3){
            honeyLatVis.showNeighbors(getHoneyNeighborInRange(i, j)
                    ,i,j);
        }else if(Geo == 6){
            triangLatVis.showNeighbors(getTriangleNeighborInRange(i, j)
                    ,i,j);
        }else if(Geo == 2 || Geo  == 4){
            ArrayList<Integer> neighLat = new ArrayList<Integer>();
            int x; int y;
            for(int u = 0; u < (2 * range + 1);u++){
                for(int v = 0; v < (2 * range + 1);v++){
                    x = ((i - range + u + L) % L);
                    y = ((j - range + v + L) % L);
                    if(!(u == i && v == j)){
                        if(!circularRange 
                                || (circularRange && Math.sqrt((v-range)*(v-range)+(u-range)*(u-range)) <= range)){
                            neighLat.add(x+y*L);
                        }
                    }
                }
            }
            sqLatVis.showNeighbors(neighLat,i,j);
        }
        updateImage();
    }
    
    /**
    *       getNeighSum should return the sum of the spins within the interaction
    * range centered by the coordinates given. 3d case
    *
    * @param i - i coordinate
    * @param j - j coordinate
    * @param k - k coordinate
    */
    @Override
    public int getNeighSum(int i, int j, int k) { return latticeSum.get(i+j*L+k*L*L); }

    /**
    *       getNeighSum should return the sum of the spins within the interaction
    * range centered by the coordinates given. 2d case
    *
    * @param i - i coordinate
    * @param j - j coordinate
    */
    @Override
    public int getNeighSum(int i, int j) { return getNeighSum(i, j, 0); }

    /**
    *       getNinRnage should return the amount of spins within the interaction
    * range. Useful for metropolis algorithm in long range case.
    */
    @Override
    public int getNinRange() {
        if(Geo == 2 || Geo  == 4){
            if (param.useLongRange) {
                return ((2 * R + 1) * (2 * R + 1) - 1);
            } else {
                return 4;
            }
        }else if(Geo == 3 || Geo  == 6){
            return (param.useLongRange) ? nInteract : Geo ;
        }
        return 0;
    }

    /**
    *       setInitialConfig should set all spins in the lattice to the values given
    * by a file in the Config directory.
    *
    * @param t - time of config to set lattice to
    * @param run - run of time to set lattice to
    * @param post - postfix or filename of lattice file
    */
    @Override
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
                int sum = 0; int sumStag = 0;
                while ( scanner.hasNextInt() ) {
                    if ( D == 2 ) {
                        int x = scanner.nextInt(); int y = scanner.nextInt();
                        int fSpin = scanner.nextInt();
                        setValue(x, y, 0, fSpin);
                        sum = fSpin + sum;
                        if((Geo == 3 && (x%2 == 0)) || (Geo !=3 && (x+y)%2 == 0)){
                                sumStag += fSpin;
                        }else{
                            sumStag -= fSpin;
                        }
                    } else {
                        int x = scanner.nextInt(); int y = scanner.nextInt();
                        int z = scanner.nextInt();
                        int fSpin = scanner.nextInt();
                        setValue(x, y, z, fSpin);
                        sum = fSpin + sum;
                        if((Geo == 3 && (x%2 == 0)) || (Geo !=3 && (x+y+z)%2 == 0)){
                                sumStag += fSpin;
                        }else{
                            sumStag -= fSpin;
                        }
                    }

                }
                magnetization.set(sum); magStaggered.set(sumStag);
                //System.out.println("AtomicLattice | value " + s)
            } catch (FileNotFoundException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }

        // If Heterogenous- Fix some spins in the fixed lattice and update the magnetization
        if( param.useHeter ){ setFixedLatticeValues(); }
    }

    /**
    *   getInstId returns the identification number of the instance of this
    * program running. Useful for allowing each instance a directory to work
    * in.
    */
    @Override
    public int getInstId() { return instId; }

    /**
    *   isPositive determines if number is positive.
    *
    * @param number
    *
    */
    public int isPositive(int number) {
        return (number >= 0) ? 1 : 0 ;
    }

    /**
    *       readFixedFile reads the file "fixed.txt" which contains the fixed lattice
    * values and updates the lattice with this information.
    *
    *
    */
    public void readFixedFile(String fname) {
        fixedFname = "fixed" + fname + ".txt";
        // Setting some spins fixed using File
        if( D == 3 ){
            fixed = new AtomicIntegerArray(L * L * L);
        } else if( D == 2 ){
            fixed = new AtomicIntegerArray(L * L);
        }
        initializeFixedArray();
        setFixedLatticeValues();
    }

    /**
    *       getFixedSubLatSum should return the sum of the fixed spin amounts in
    * sublattice in the index given
    *
    * @param j - index of sublattice
    */
    @Override
    public int getFixedSubLatSum(int j) { return fixedSubLatticeSum.get(j); }

    /**
    *   getFixSpinVal should return the average value of the fixed spins
    */
    @Override
    public int getFixSpinVal() { return fixSpinVal; }

    /**
    *       getHighFixIndex returns the index of the highest density of fixed spins
    */
    @Override
    public int getHighFixIndex() {
        // Find max fixed sum area
        int iMax = 0; int maxM = 0; int subL = (int) L / granularity;
        for( int i = 0; i < (subL * subL); i++ ){
            if (fixedSubLatticeSum.get(i) > maxM) { iMax = i; }
        }
        return iMax;
    }

    /**
    *       highFixedAndStableMatched should return if the sublattice index for the
    *   highest density stable spins and fixed spins match.
    */
    @Override
    public boolean highFixedAndStableMatched() {
        // Do not bother if not granular
        if (granularity == 1) {
            return false;
        } else {
            // Find max stable sum area
            int iMax = 0; int maxM = 0;
            int subL = (int) L / granularity;
            for (int i = 0; i < (subL * subL); i++) {
                if ((subLatticeSum.get(i) * stable) > maxM) { iMax = i; }
            }
            // if area of max density is equal to max stable density return true
            if (fixedSubLatticeSum.get(iMax) == maxFixedSum.get()) { return true; }
        }
        return false;
    }

    /**
    * getHighMagIndex returns the index of the highest density of stable spins
    */
    @Override
    public int getHighMagIndex() {
        // Find max stable sum area
        int iMax = 0; double maxM = 0;
        double sumM = 0.0; double effM;
        int subL = (int) L / granularity;

        for (int i = 0; i < (subL * subL); i++) {
            // Calculate effective magnetization ie dont count the fixed spins
            effM = (subLatticeSum.get(i) - getFixedSubLatSum(i) * getFixSpinVal()) / (subL * subL - getFixedSubLatSum(i));
            if( (effM * stable) > maxM ){ iMax = i; }
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
    public double getLastAvgMagSubLat() { return lastMagSubAvg; }

    /**
    * getSubLatSum should return the sum of the sublattice in the index given
    *
    * @param u - index of sublattice
    */
    @Override
    public int getSubLatSum(int u) { return subLatticeSum.get(u); }

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
    
    public void initShowImage(){
        canvasFrame = new CanvasFrame("Lattice");
        canvasFrame.setCanvasSize(getVisImg().getWidth(),getVisImg().getHeight());
    }
    
    public void updateImage(){ canvasFrame.showImage(getVisImg()); }
    
    public void endShowImage(){ canvasFrame.dispose(); }
    
    // test the class
    public static void main(String[] args) {
        AtomicLatticeSumSpin lat = new AtomicLatticeSumSpin(1);
        System.out.println("AtomicLattice |  Sum of neighbors " + lat.getNeighSum(3, 1));
        System.out.println("AtomicLattice |  length: " + lat.L+"    range: "+lat.R+"   Geo: "+lat.getGeo());
        //System.out.println("AtomicLattice | index high fixed: " + lat.getHighFixIndex());
        System.out.println("AtomicLattice | __________________________________");
        System.out.println(lat.magnetization);
        lat.saveLatticeImage("");
        lat.makeVideo();
        lat.spinDraw();
        lat.initShowImage();
        lat.updateImage();
        //lat.setValue(3, 4, 0, -1);
        //lat.setValue(3, 5, 0, -1);
        //lat.setValue(3, 6, 0, -1);
        //lat.setValue(4, 4, 0, -1);
        lat.showNeighbors(32, 34, 1, -3);
        //lat.endShowImage();
        System.out.println("AtomicLattice |  Sum of neighbors " + lat.getNeighSum(0, 0));
        lat.saveLatticeImage("post");
        System.out.println("AtomicLattice | N in range: " + lat.getNinRange() + "  with R : " + lat.R);
        System.out.println("AtomicLattice | Done!");
    }
}

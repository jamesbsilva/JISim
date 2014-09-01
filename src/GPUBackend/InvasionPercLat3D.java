package GPUBackend;

/**
*   @(#) InvasionPercLat3D
*/  

import AnalysisAndVideoBackend.DisplayMakerCV;
import Backbone.Algo.MCAlgo;
import Backbone.System.AtomicBondsSiteLattice;
import Backbone.System.SimSystem;
import Backbone.Util.DataSaver;
import Backbone.Util.DirAndFileStructure;
import Backbone.Util.ParameterBank;
import JISim.MCSimulation;
import Triggers.*;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;
import scikit.graphics.dim3.Grid3DCL;
import scikit.jobs.SimulationCL;
import scikit.opencl.CLHelper;
import scikit.opencl.RandomNumberCL;

/**  
*   Basic Implementation of InvasionPercLat3D Algorithm Monte Carlo Simulation  
*  <br>
*  @param mc - MCSimulation - Monte Carlo simulation class hosting 
* 
* @author      James Silva <jbsilva @ bu.edu>                 
* @since       2013-07    
*/
public class InvasionPercLat3D implements MCAlgo{
    private CLHelper clhandler;
    private ParameterBank param;private Trigger trigger;
    private AtomicBondsSiteLattice lattice;
    private boolean useLongRange;private boolean useDilution;private boolean useHeter;
    private int tTrigger = 0; private int tAccumulate = 1000; // begin accumulator
    private boolean suppressFieldMessages=false;
    private int coordination = 6;
    private int initLocation;
    
    private int L;private int R;
    private int N;private int Q;
    private int BondsN;
    private int run;private int maxThreads;
    private boolean output = true;
    private Random Ran = new Random();private RandomNumberCL RNG;
    private String SimPost;private String ParamPost;  
    private int minX;private int maxX;
    private int minY;private int maxY;
    private int minZ;private int maxZ;  
    private int percSize = 0;  
    private int maxSubL = 6; 
    private int Chunk = 8;private int ChunkL = 1;private int ChunkNum = 8; 
    private int ChunkXrun = 0;private int ChunkYrun = 0;private int ChunkZrun = 0;
    private int ChunkXoffset = 0;private int ChunkYoffset = 0;private int ChunkZoffset = 0; 
    private int minXlsize = 0;private int minYlsize = 0;private int minZlsize = 0;
    private int lsizeN = 1;private int lsizeBondsN = 1;private int lsizeChunkN = 1; 
    private int subL = L;  
    private int boxCountSubL = L/2;private int boxCountDiv = 4; private int boxDataSize = L/2;
    private int consolidateN = 1;  
    public int currentSeed;
    public int mcs = 0; // number of MC moves per spin
    public int tFlip = 0;public int currentTime;
    public int lastBond;
    private String growKernel = "invperc3Dgrow";private String minKernel = "invperc3Dmin";
    private String reduceKernel = "invperc3Dreduce";private String consKernel = "invperc3Dconsolidate";
    private String minPrepKernel = "invperc3Dminprep";private String RNGKernel = "random123";
    private String sortKernel = "localSort3D";private String boxKernel = "boxCounting";
    private String initKernel = "makeOrderIntBuffer";
    private Grid3DCL grid = null;
    private boolean isSpanning=false; private boolean savedSpanning =false;
    private boolean crossedEdge=false; private boolean savedCrossEdge =false;
    
    
    /**  @param mc - MCSimulation - Monte Carlo simulation class hosting 
    *     @param out - boolean to determine if outputting text confirmations into console
    */
    public InvasionPercLat3D(MCSimulation mc, boolean out,String device){
        param = new ParameterBank(mc.getParamPostFix());
        SimPost = mc.getSimPostFix();ParamPost = mc.getParamPostFix();
        useLongRange = param.useLongRange;useDilution = param.useDilution;
        useHeter = param.useHeter;
        L = param.L;N = param.N;R = param.R;
        initLocation = (int)((double)N/2.0);
        
        // determine if outputting
        output = out;
        if(output){System.out.print("Run:"+mc.getRuns()+"     ");}
        run = mc.getRuns();  
        
        if( mc.getSimSystem().getClass().getName().contains("AtomicBondsSiteLattice") ){
            lattice = (AtomicBondsSiteLattice)mc.getSimSystem();
        }else{
            System.err.println("THIS CLASS (InvasionPercLat NEEDS A ATOMICBONDSSITELATTICE AS SYSTEM.");
        }
        
        clhandler = new CLHelper();
        
        grid = new Grid3DCL("InvasionPercLat3D | ",clhandler,param.L,param.L,param.L);
        if(grid != null){
            grid.getCLScheduler().waitForInitiatedCLGL();
            clhandler = grid.getCLHelper();
        }else{
            clhandler.initializeOpenCL(device);
        }
        
        percSize = lattice.getMagnetization();
        
        // initialize accumulators and energy
        resetData();    
        currentTime=0;
    }
    
    /**  
    *     @param out - boolean to determine if outputting text confirmations into console
    */
    public InvasionPercLat3D(String par, String sim,int runIn, boolean out){
        param = new ParameterBank(par);
        SimPost = sim;ParamPost = par;
        useLongRange = param.useLongRange;useDilution = param.useDilution;
        useHeter = param.useHeter;
        L = param.L;N = param.N;R = param.R;
        BondsN = coordination*N;
        int initX = (int)((double)L/2.0);int initY = (int)((double)L/2.0);int initZ = (int)((double)L/2.0);
        initLocation = initX+initY*L+initZ*(L*L);
        boxCountSubL = initX;boxCountDiv = 5;
        
        // determine if outputting
        output = out;
        if(output){System.out.print("Run:"+runIn+"     ");}
        run = runIn;  
        
        String device = "GPU";
        clhandler = new CLHelper();
        
        grid = new Grid3DCL("InvasionPercLat3D | ",clhandler,param.L,param.L,param.L,0.2);
        if(grid != null){
            grid.getCLScheduler().waitForInitiatedCLGL();
            clhandler = grid.getCLHelper();
        }else{
            clhandler.initializeOpenCL(device);
        }
        clhandler.setPrintMode(false);
        maxThreads = (int)clhandler.getCurrentDevice1DMaxWorkItems()/2;
        lsizeN = clhandler.maxLocalSize1D(N);
        lsizeBondsN = clhandler.maxLocalSize1D(BondsN);
        
        
        // build necessary kernels
        clhandler.createKernel(sortKernel);clhandler.createKernel(initKernel);
        clhandler.createKernel(growKernel);clhandler.createKernel(consKernel);
        clhandler.createKernel(minPrepKernel);clhandler.createKernel(minKernel);
        clhandler.createKernel(reduceKernel);clhandler.createKernel(boxKernel);

        subL = getSubL(L);
        Chunk = (int) Math.pow(subL,3);
        ChunkNum = (int)((double)N/(double)Chunk);
        ChunkL = (int)((double)L/(double)subL);
        lsizeChunkN = clhandler.maxLocalSize1D(ChunkNum);
        
        System.out.println(" Chunks Total | "+ChunkNum + " Chunk Size | "+Chunk + "   N | "+N );
        
        // Initialize bond order buffer 
        clhandler.createIntBuffer(initKernel, 0, BondsN, 0, "rw");
        clhandler.setIntArg(initKernel, 0, BondsN);
        clhandler.setIntArg(initKernel, 1, BondsN);
        
        // sort setup 
        clhandler.createFloatBuffer(sortKernel, 0, BondsN, 0.0f, "rw");
        clhandler.copyIntBufferAcrossKernel(initKernel, 0, sortKernel, 0);
        clhandler.setIntArg(sortKernel, 0, subL);
        clhandler.setIntArg(sortKernel, 1, L);
        clhandler.setIntArg(sortKernel, 2, coordination*Chunk);
        clhandler.setIntArg(sortKernel, 3, ChunkL);
        
        // initialize bond order array
        clhandler.runKernel(initKernel, BondsN, lsizeBondsN);
        
        // throw random bonds
        RandomNumberCL RNG = new RandomNumberCL(Ran.nextInt(),clhandler);
        int RNGKey = RNG.addRandomBuffer(sortKernel, 0, BondsN);
        RNG.fillBufferWithRandom(RNGKey);
        int temp = 1;
        
        // pre sort
        clhandler.runKernel(sortKernel, ChunkNum, lsizeChunkN);
 
        //post sort
        int[] paramInt = new int[10];
        paramInt[0] = L;
        paramInt[1] = subL;
        paramInt[2] =  (Chunk*coordination);
        paramInt[3] = ChunkL;
        
        minX = ((initX - subL) > 0) ? initX - subL - 1: 0;
        maxX = ((initX + subL) < L) ? initX + subL + 1: L;
        minY = ((initY - subL) > 0) ? initY - subL - 1: 0;
        maxY = ((initY + subL) < L) ? initY + subL + 1: L;
        minZ = ((initZ - subL) > 0) ? initZ - subL - 1: 0;
        maxZ = ((initZ + subL) < L) ? initZ + subL + 1: L;
        
        //maxX = L;maxZ = L;maxY = L;
        //minX = 0;minY = 0;minZ = 0;
        
        paramInt[4] = minX;paramInt[5] = maxX;
        paramInt[6] = minY;paramInt[7] = maxY;
        paramInt[8] = minZ;paramInt[9] = maxZ;
        
        System.out.println("System of size | "+N+"    | Using initial bond | "
                +((initLocation-1)*coordination+1)+"     initialLoc | "
                +initLocation+"    subL : "+subL+"    L: "+L+"   chunkL | "+ChunkL);
        
        // grow kernel 
        clhandler.createIntBuffer(growKernel, 0, N, 0, "rw");
        clhandler.createIntBuffer(growKernel, 1, BondsN, 0, "rw");
        clhandler.createIntBuffer(growKernel, 2, paramInt.length,paramInt,0,true);
        clhandler.setIntArg(growKernel, 0,(initLocation-1)*coordination+1);
        clhandler.setIntArg(growKernel, 1,(initLocation));
        clhandler.setIntArg(growKernel, 2,coordination);
        
        // error buffer
        clhandler.createFloatBuffer(growKernel, 0, 20, 0.0f, "rw");
        
        clhandler.runKernel(growKernel, coordination, coordination);
        
        // set init location to unreachable value once initialized
        clhandler.setIntArg(growKernel, 1,-1);
        
        // minimum bond in fringe kernel
        clhandler.copyIntBufferAcrossKernel(growKernel, 0, minKernel, 0);
        clhandler.copyFlBufferAcrossKernel(sortKernel, 0, minKernel, 0);
        clhandler.copyIntBufferAcrossKernel(sortKernel, 0, minKernel, 1);
        clhandler.copyIntBufferAcrossKernel(growKernel, 1, minKernel, 2);
        clhandler.createIntBuffer(minKernel, 3, ChunkNum, -1, "rw");
        clhandler.copyIntBufferAcrossKernel(growKernel, 2, minKernel, 4);
        clhandler.setIntArg(minKernel, 0,0);
        clhandler.setIntArg(minKernel, 1,0);
        clhandler.setIntArg(minKernel, 2,0);
        clhandler.setIntArg(minKernel, 3,ChunkNum);
        //error buffer
        clhandler.copyFlBufferAcrossKernel(growKernel, 0, minKernel, 1);
        
        // reduction kernel
        clhandler.copyIntBufferAcrossKernel(growKernel, 0, reduceKernel, 0);
        clhandler.copyFlBufferAcrossKernel(sortKernel, 0, reduceKernel, 0);
        clhandler.copyIntBufferAcrossKernel(minKernel, 3, reduceKernel, 1);
        clhandler.createFloatBuffer(reduceKernel, 1, ChunkNum, 0.0f, "rw");
        clhandler.setIntArg(reduceKernel, 0,1);
        clhandler.setIntArg(reduceKernel, 1,ChunkNum);
        //error buffer
        clhandler.copyFlBufferAcrossKernel(growKernel, 0, reduceKernel, 2);
        
        // consolidate kernel
        clhandler.copyIntBufferAcrossKernel(reduceKernel, 1, consKernel, 0);
        clhandler.copyFlBufferAcrossKernel(reduceKernel, 1, consKernel, 0);
        clhandler.setIntArg(consKernel, 0,consolidateN);
        clhandler.setIntArg(consKernel, 1,ChunkNum);
        //error buffer
        clhandler.copyFlBufferAcrossKernel(growKernel, 0, consKernel, 1);
       
        // consolidate kernel
        clhandler.copyIntBufferAcrossKernel(reduceKernel, 1, minPrepKernel, 0);
        clhandler.copyFlBufferAcrossKernel(reduceKernel, 1, minPrepKernel, 0);
        clhandler.setIntArg(minPrepKernel, 0,ChunkNum);
        //error buffer
        clhandler.copyFlBufferAcrossKernel(growKernel, 0, minPrepKernel, 1);

        // box counting analysis
        boxDataSize = (int)(Math.pow(boxCountSubL,3.0)*(1-Math.pow(0.5, boxCountDiv+1)/(1-0.5)));
        clhandler.copyIntBufferAcrossKernel(growKernel, 0, boxKernel, 0);
        clhandler.createIntBuffer(boxKernel, 1, boxDataSize, 0, "rw");
        clhandler.setIntArg(boxKernel, 0,L);
        clhandler.setIntArg(boxKernel, 1,N);
        clhandler.setIntArg(boxKernel, 2,boxCountSubL);
        clhandler.setIntArg(boxKernel, 3,3);
        clhandler.setIntArg(boxKernel, 4,boxCountDiv);
        
        
        //checkChunk();
        //clhandler.getIntBufferAsArray(minKernel, 3, ChunkNum, true);
        //clhandler.getFloatBufferAsArray(growKernel, 0, 10, true);
  
        // initialize accumulators and energy
        resetData();    
        currentTime=0;
        updateChunkBounds();
        
        if(grid != null){
            // Coloring kernel
            String colorKernel = "colorLatticeInt2";
            clhandler.createKernel(colorKernel);
            clhandler.copyIntBufferAcrossKernel(minKernel, 0, colorKernel, 0);
            grid.getCLScheduler().schedCopyColSharedBuffer(colorKernel, 0);
            System.out.println("Number of Threads: "+N);
            grid.getCLScheduler().setColorKernel1D(colorKernel, N, subL*subL*subL);
            SciKitSimInterface scikit = new SciKitSimInterface(this);
            System.out.println("SCIKIT INTERFACE | "+scikit);
            grid.getCLScheduler().setSimulation((SimulationCL)scikit);
            grid.getCLScheduler().setKernelChunkRunsOff(true);
        }
    }
    
    public void findPossibleMinimumBond(){
        //printMinKernelRunParam();
        clhandler.runKernel3D(minKernel, ChunkXrun, minXlsize,ChunkYrun, minYlsize, ChunkZrun, minZlsize);
        //clhandler.getIntBufferAsArray(minKernel, 3, ChunkNum, true);
        //System.out.println("++++++++++++++++++++++++++++");
        //clhandler.getFloatBufferAsArray(growKernel, 0, 6, true);
        //System.out.println("++++++++++++++++++++++++++++");
        //clhandler.getIntBufferAsArray(minKernel, 0, N, true);
        //System.out.println("++++++++++++++++++++++++++++");
    }
    
    public void printMinKernelRunParam(){
        System.out.println("ChunkRunX | "+ChunkXrun+"   from : "+ChunkXoffset+
                "    ChunkRunY | "+ChunkYrun+"   from : "+ChunkYoffset+
                "     ChunkRunZ | "+ChunkZrun+"   from : "+ChunkZoffset
                +"    offsetFromKernelX | "+clhandler.getIntArg(minKernel, 0)
                +"    offsetFromKernelY | "+clhandler.getIntArg(minKernel, 1)
                ); 
    }
    
    public int getNextMinBond(){
        reduceAndConsolidatePossibleMin();
        int redNum = ChunkNum/(1+consolidateN);
        ArrayList<Integer> possibleMinBond = clhandler.getBufferIntAsArrayList(consKernel, 0, redNum, false);
        ArrayList<Float> possibleMinVal = clhandler.getBufferFlAsArrayList(consKernel, 0, redNum, false);
        float minVal = 1.0f;
        int minBond = 0;
        for(int u = 0; u < redNum;u++){
            if(possibleMinVal.get(u) < minVal){
                minVal = possibleMinVal.get(u);
                minBond = possibleMinBond.get(u);
            }
        }
        return minBond;
    }
    
    public void growToMinBondLoc(){
        int minBond = getNextMinBond();
        //System.out.println("Growing to | "+minBond);
        getSiteOutFromBond(minBond);
        updateChunkBounds();
        //if(lastBond == minBond){this.turnOffoffset();}    
        //turnOffoffset();    
        clhandler.setIntArg(growKernel, 0,minBond);
        clhandler.setKernelArg(growKernel,true);
        clhandler.runKernel(growKernel, coordination, coordination);
        percSize++;
        lastBond = minBond;
        
        int redNum = ChunkNum/(1+consolidateN);
        clhandler.runKernel(minPrepKernel, ChunkNum, lsizeChunkN);
        //clhandler.getIntBufferAsArray(growKernel, 1, BondsN, true);
        //System.out.println("++++++++++++++++++++++++++++");
        //clhandler.getFloatBufferAsArray(growKernel, 0, 10, true);
    }
    
    public void reduceAndConsolidatePossibleMin(){
        int redNum = ChunkNum/(1+consolidateN);
        for(int u = 0; u < (consolidateN);u++){
            redNum = ChunkNum/(1+u);
            clhandler.setIntArg(reduceKernel, 0,u+1);
            clhandler.runKernel(reduceKernel, redNum, (u+1));
        }
        redNum = ChunkNum/(1+consolidateN);
        clhandler.runKernel(consKernel, redNum, consolidateN);
    }
    
    public int getSubL(int L){
        int subL = 1;
        int ind = 1;
        while(ind < maxSubL){
            ind++;
            if(L % ind == 0){subL = ind;}
        }
        return subL;
    }
    
    public void checkChunk(){
        int chunkX = 0; int chunkY = 0; int chunkZ =0;
        int x; int y; int z;int ind; int indSite = 0;
        int[] order = clhandler.getIntBufferAsArray(sortKernel, 0, BondsN , false);
        float[] values = clhandler.getFloatBufferAsArray(sortKernel, 0, BondsN , false);
        for(int u = 0; u < (Chunk);u++){
            int subX = u%subL;
            int subY = (int)((double)u/(double)subL)%subL;
            int subZ = (int)((double)u/(double)(subL*subL))%subL;
            x = subX + chunkX*subL;
            y = subY + chunkY*subL;
            z = subZ + chunkZ*subL;
            indSite = (x+y*L+ z*(L*L))*coordination;
            for(int v = 0; v < coordination;v++){
                ind = indSite + v;       
                System.out.println("Position Ordered | "+(ind) 
                        + "   PositionActual | "+ (order[ind])+ "    Value | "+values[order[ind]]
                        +"   mult | "+v+"   x | "+x+"   y | "+y+"   z | "+z );
            }
        }
    }
    
    
    /**+
    *        doOneStep should do a single monte carlo step.
    */ 
    public void doOneStep(){
        currentTime++;
        findPossibleMinimumBond();
        growToMinBondLoc();
        if(currentTime == (-1)){
            int div = 2;
            clhandler.runKernel(boxKernel, 
                    (int)Math.pow((L/subL),3), 
                    clhandler.maxLocalSize1D((int)Math.pow((L/subL),3)/20));
            clhandler.getIntBufferAsArray(boxKernel, 1, boxDataSize, true);
            System.out.println("Done with box count");
        }
        if(isSpanning){
            //if(!savedSpanning)saveClusterData();
            savedSpanning = true;
            
        }
        if(crossedEdge){
            if(!savedCrossEdge)saveClusterData();
            savedCrossEdge = true;
        }
        //System.out.println("Min/Max X | "+minX+"  "+maxX+"   Y | "+minY+"  "+maxY+"    Z | "+minZ+"    "+maxZ);
    }
    
    public void saveClusterData(){
        String nameOut = "InvasionClust3D";
        String fname = nameOut+"cluster-L"+L+"-size-"+percSize+".dat";
        DirAndFileStructure dir = new DirAndFileStructure();
        fname = dir.getDataDirectory(nameOut)+fname;
        System.out.println("Saving cluster | "+fname);
        //read file first    
        /*
        try{
            PrintStream out = new PrintStream(new FileOutputStream(
                fname,false));
                ArrayList<Integer> temp = clhandler.getBufferIntAsArrayList(growKernel, 0, N, false);
                for(int u = 0; u < temp.size();u++){
                    out.println(temp.get(u));
                }
                out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }*/
        fname= "percSize"+"cluster-L"+L+".dat";
        fname = dir.getDataDirectory(nameOut)+fname;
        //read file first    
        try{
            PrintStream out = new PrintStream(new FileOutputStream(
                fname,true));
                out.println(percSize);
                out.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
    
    private int getX(int ind){
        return (ind%L);
    }
    private int getY(int ind){
        return ((int)((double)ind/(double)L)%L);
    }
    private int getZ(int ind){
        return ((int)((double)ind/(double)(L*L))%L);
    }
    
    private int getSiteInFromBond(int bind){
        return (int) Math.floor((double)bind/(double)coordination);
    }
    
    private int getSiteOutFromBond(int bind){
        int site = (int) Math.floor((double)bind/(double)coordination);
        int ind = bind%coordination; 
        int i = getX(site);int j = getY(site);
        int k = getZ(site);
        switch( ind ){
            case 0: j = ((j+1)%L);
                    if((j+subL) > L){maxY = L;}else{
                        if((j+subL) > maxY){
                            maxY = (j+subL);
                        }
                    } 
                    break;
            case 1: i = ((i+1)%L);
                    if((i+subL) > L){maxX = L;}else{
                        if((i+subL) > maxX){
                            maxX = (i+subL);
                        }
                    } 
                    break;
            case 2: j = ((j+L-1)%L);
                    if((j-subL) < 0){minY = 0;}else{
                        if((j-subL) < minY){
                            minY = (j-subL);
                        }
                    }
                    break;
            case 3: i = ((i+L-1)%L);
                    if((i-subL) < 0){minX = 0;}else{
                        if((i-subL) < minX){
                            minX = (i-subL);
                        }
                    } 
                    break;
            case 4: k = ((k+1)%L);
                    if((k+subL) > L){maxZ = L;}else{
                        if((k+subL) > maxZ){
                            maxZ = (k+subL);
                        }
                    }
                    break;
            case 5: k = ((k+L-1)%L);
                    if((k-subL) < 0){minZ = 0;}else{
                        if((k-subL) < minZ){
                            minZ = (k-subL);
                        }
                    } 
                    break;
            default:ind = site;
                     break;
        }
        ind = i+j*L+k*(L*L);
        //System.out.println("Bond | "+bind+"  going to ind | "+ind+"     x: "+i+"   y: "+j+"   z : "+k);
        return ind;
    }
    
    public void updateChunkBounds(){
        int chunkXmax = (int)Math.ceil((double)maxX/(double)subL);
        int chunkYmax = (int)Math.ceil((double)maxY/(double)subL);
        int chunkZmax = (int)Math.ceil((double)maxZ/(double)subL);
        int chunkXmin = (int)Math.floor((double)minX/(double)subL);
        int chunkYmin = (int)Math.floor((double)minY/(double)subL);
        int chunkZmin = (int)Math.floor((double)minZ/(double)subL);
        int chunkDiffX = ((maxX-minX) > 0) ? chunkXmax - chunkXmin : ChunkL;
        int chunkDiffY = ((maxY-minY) > 0) ? chunkYmax - chunkYmin : ChunkL;
        int chunkDiffZ = ((maxZ-minZ) > 0) ? chunkZmax - chunkZmin : ChunkL;
        /*System.out.println("maxX : "+maxX+"   minX : "+minX);
        System.out.println("maxY : "+maxY+"   minY : "+minY);
        System.out.println("maxZ : "+maxZ+"   minZ : "+minZ);
        System.out.println("chunkmaxX : "+chunkXmax+"   chunkminX : "+chunkXmin);
        System.out.println("chunkmaxY : "+chunkYmax+"   chunkminY : "+chunkYmin);
        System.out.println("chunkmaxZ : "+chunkZmax+"   chunkminZ : "+chunkZmin);
        */
        if(chunkDiffX == 0){
            chunkDiffX = 2;
            chunkXmin--;
        }
        if(chunkDiffY == 0){
            chunkDiffY = 2;
            chunkYmin--;
        }
        if(chunkDiffZ == 0){
            chunkDiffZ = 2;
            chunkZmin--;
        }
        
        if(ChunkXrun < chunkDiffX){
            ChunkXrun = chunkDiffX;
            ChunkXoffset = chunkXmin;
            minXlsize = findMaxLsize3D(ChunkXrun);
            clhandler.setIntArg(minKernel,0,chunkXmin);
            clhandler.setKernelArg(minKernel);
        }else if(ChunkXoffset > chunkXmin){
            ChunkXrun = ChunkXrun +ChunkXoffset - chunkXmin;
            ChunkXoffset = chunkXmin;
            minXlsize = findMaxLsize3D(ChunkXrun);
            clhandler.setIntArg(minKernel,0,chunkXmin);
            clhandler.setKernelArg(minKernel);
        } 
        if(ChunkYrun < chunkDiffY){
            ChunkYrun = chunkDiffY;
            ChunkYoffset = chunkYmin;
            minYlsize = findMaxLsize3D(ChunkYrun);
            clhandler.setIntArg(minKernel,1,chunkYmin);
            clhandler.setKernelArg(minKernel);
        }else if(ChunkYoffset > chunkYmin){
            ChunkYrun = ChunkYrun +ChunkYoffset - chunkYmin;
            ChunkYoffset = chunkYmin;
            minYlsize = findMaxLsize3D(ChunkYrun);
            clhandler.setIntArg(minKernel,1,chunkYmin);
            clhandler.setKernelArg(minKernel);
        }
        if(ChunkZrun < chunkDiffZ){
            ChunkZrun = chunkDiffZ;
            ChunkZoffset = chunkZmin;
            minZlsize = findMaxLsize3D(ChunkZrun);
            clhandler.setIntArg(minKernel,2,chunkZmin);
            clhandler.setKernelArg(minKernel);
        }else if(ChunkZoffset > chunkZmin){
            ChunkZrun = ChunkZrun +ChunkZoffset - chunkZmin;
            ChunkZoffset = chunkZmin;
            minZlsize = findMaxLsize3D(ChunkZrun);
            clhandler.setIntArg(minKernel,2,chunkZmin);
            clhandler.setKernelArg(minKernel);
        }
        checkIfSpanning();
        //turnOffoffset();
    }
    
    private void checkIfSpanning(){
        boolean spanning = false;
        double span = 0;
        //check X
        span = ( (maxX - minX )> 0) ? maxX - minX : minX-maxX;
        if( span > (L-2) ){
            spanning = true;
        }
        //check Y
        span = ( (maxY - minY )> 0) ? maxY - minY : minY-maxY;
        if( span > (L-2) ){
            spanning = true;
        }
        //check X
        span = ( (maxZ - minZ )> 0) ? maxZ - minZ : minZ-maxZ;
        if( span > (L-2) ){
            spanning = true;
        }
        isSpanning = spanning;
        if( maxX >= L || maxY >= L || maxZ >= L ){crossedEdge = true;}
        if( minX <= 0 || minY <= 0 || minZ <= 0 ){crossedEdge = true;}
    }
    
    public int findMaxLsize3D(int xchunk){
        int max = 1;int temp;
        for(int u = 1; u < 7;u++){
            if(xchunk % u == 0){
                max = u;
            }
        }
        return max;
    }
    
    public void turnOffoffset(){
        ChunkZrun = ChunkL;ChunkXrun = ChunkL;ChunkYrun = ChunkL;
        minXlsize = clhandler.maxLocalSize1D(ChunkL);
        minYlsize = minXlsize;
        minZlsize = minXlsize;
        ChunkXoffset = 0;ChunkYoffset = 0;ChunkZoffset = 0;
        clhandler.setIntArg(minKernel,0,0);
        clhandler.setIntArg(minKernel,1,0);
        clhandler.setIntArg(minKernel,2,0);   
    }
    
    
    /**
    *         setTrigger should initialize the trigger to be used in the simulation
    */ 
    public void setTrigger(){
        trigger = param.setProperTrigger(param, trigger, SimPost, output);
    }

    /**
    *         getConfigRange should save the configurations ins the given range
    * 
    *  @param tInitial - initial time in range
    *  @param tFinal - final time in range
    */ 
    @Override
    public void getConfigRange(int tInitial, int tFinal){
        //suppressFieldMessages = true;
        resetSimulation(currentSeed);
        DataSaver save = new DataSaver(lattice.getInstId());
        //assert field is right
        for(int t = 0; t<tFinal;t++){
            doOneStep();
            if(t>tInitial){save.saveConfig(lattice,run,t);}
        }
    }
    
    /**
    *   getSimSystem should return the lattice object
    * 
    */
    @Override
    public SimSystem getSimSystem(){return (SimSystem)lattice;}
    
    /**
    *   getRun should return the current run number
    * 
    */
    @Override
    public int getRun(){return run;}
    
    /**
    *         setRun should set the number of the current run to the value given.
    *
    *  @param cr - current run value 
    */
    @Override
    public void setRun(int cr){run = cr;}
    
    /**
    *         setSeed should set the random number seed to the value given
    * 
    *  @param seed - random number seed
    */ 
    public void setSeed(int seed){
        currentSeed = seed;
        Ran.setSeed(seed);
    }
    /**
    *         resetData resets all the accumulators of data and number of mcsteps made.
    */
    public void resetData() {
        mcs = 0;
    }
    /**
    *         resetSimulation should reset all simulation parameters.
    */
    @Override
    public void resetSimulation(){
        currentTime = 0;
        resetData();
        trigger.reset();
        lattice.initialize(param.s);
        percSize = lattice.getMagnetization();
    }
    /**
    *         resetSimulation should reset all simulation parameters but 
    *  set seed to the given value.
    * 
    *  @param seed - random number seed
    */ 
    @Override
    public void resetSimulation(int seed){
        currentTime = 0;
        resetData();
        trigger.reset();
        lattice.initialize(param.s);
        if(param.useHeter){lattice.setFixedLatticeValues();}
        percSize = lattice.getMagnetization();
        setSeed(seed);
    }
    /**
    *         resetSimulation should reset all simulation parameters, set
    *   all lattice values to lattice at the time given, and set the 
    *   seed to the given value
    * 
    *  @param t - time to set lattice to
    *  @param seed - random number seed
    */ 
    @Override
    public void resetSimulation(int t, int seed){
        currentTime = t;
        resetData();
        //change to allow more triggers
        trigger.reset2();//reset but keep std or go into intervention mode
        if(t==0){lattice.initialize(param.s);
        if(param.useHeter){lattice.setFixedLatticeValues();}
        }else{
        lattice.setInitialConfig(t,run,"");}
        percSize = lattice.getMagnetization();
        setSeed(seed);
    }
    
    /**
    *         setMeasuringStartTime sets the time to begin all measurement of all variables
    *   by a system measurer.
    * 
    *  @param tin - initial measurement time of all data
    */ 
    @Override
    public void setMeasuringStartTime(int tin){
        tAccumulate = tin;
    }
    
    /**
    *   getSeed should return the current random number seed
    * 
    */
    @Override
    public int getSeed(){return currentSeed;}

    // test the class
    public static void main(String[] args) {
        DisplayMakerCV dispConf = new DisplayMakerCV( "Configuration" );
        InvasionPercLat3D inv = new InvasionPercLat3D("-perc3D","",0,true);
        //for(int u = 0; u < 15000; u++){
        //    inv.doOneStep();
        //}        
    }
}

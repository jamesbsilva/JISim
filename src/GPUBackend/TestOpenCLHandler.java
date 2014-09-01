package GPUBackend;

import java.io.IOException;
import java.util.Random;
import scikit.opencl.CLHelper;
import scikit.opencl.RandomNumberCL;

/**
 *
 * @author James B Silva
 */
public class TestOpenCLHandler {

    public static void testFloatReduce(){
        CLHelper clhandle = new CLHelper();
        int esize = 2500000;                                  // Length of arrays to process
        int lsize = 100;  // Local work size dimensions
        int gsize = esize;   // rounded up to the nearest multiple of the localWorkSize
        int lsize2 =25;
        int  block= (int)(esize/lsize2);
        
        String kernel = "iterFl";
        String kernel2 = "reduceSumFl";

        clhandle.initializeOpenCL("GPU");

        clhandle.createKernel(kernel);
        clhandle.createKernel(kernel2);

        clhandle.createFloatBuffer(kernel, 0, esize, 0.0f, 0);
        clhandle.setFloatArg(kernel, 0, 2.5f);
        clhandle.setIntArg(kernel, 0,esize);
        clhandle.copyFlBufferAcrossKernel(kernel, 0, kernel2, 0);
        clhandle.setIntArg(kernel2, 0, esize);
        //clhandle.setIntArg(kernel2, 1, block);
        clhandle.createFloatBuffer(kernel2, 1, block, 0.0f, 0);

        // The moving of the buffers into the device doesnt happen until the queue runs.
        clhandle.runKernel(kernel, gsize, lsize);

        clhandle.getFloatBufferAsArray(kernel, 0, 3, true);            

        clhandle.runKernel(kernel, gsize, lsize);

        clhandle.getFloatBufferAsArray(kernel, 0, 3, true);            

        clhandle.runKernel(kernel, gsize, lsize);

        clhandle.getFloatBufferAsArray(kernel, 0, 3, true);            

        System.out.println("Divide array into "+block+" chunks");
        long time =System.nanoTime();
        clhandle.runKernel(kernel2, gsize, lsize2);
        time =System.nanoTime()-time;
        
        float[] subsum=clhandle.getFloatBufferAsArray(kernel2, 1, block, false);
        // Better to use a double if possible to get less rounding errors
        double sum = 0.0f;
        for(int i =0;i<subsum.length;i++){
            sum+= subsum[i];
        }
        System.out.println("GPU Calculation in "+(time/1000)+"  us  with sum : " + sum +"    of "+esize+"  elements");
        
        
        float[] arr =clhandle.getFloatBufferAsArray(kernel2, 0, esize, false);
        
        time =System.nanoTime();
        double sum2 = 0.0f;
        for(int i =0;i<arr.length;i++){
            sum2+= arr[i];
        }
        time =System.nanoTime()-time;
        System.out.println("CPU Calculation in "+(time/1000)+"  us  with sum : " + sum2+"    of "+esize+"  elements");
        
        System.out.println("val theory: "+(esize-1)*7.5);

        // R = 25 is the point of returns for GPU method
        clhandle.closeOpenCL();
    }

    public static void testIntReduce(){
        CLHelper clhandle = new CLHelper();
        clhandle.listAllDevices("");
        int esize = 5*256*256*256;                                  // Length of arrays to process
        int lsize = clhandle.maxLocalSize1D(esize);  // Local work size dimensions
        int gsize = esize;   // rounded up to the nearest multiple of the localWorkSize
        // Local workgroup size dependent on GPU - Too many threads -> bad result
        int lsize2 =25;
        int  block= (int)(esize/lsize2);
        
        String kernel = "iterInt";
        String kernel2 = "reduceSumInt";

        clhandle.initializeOpenCL("GPU");

        clhandle.createKernel(kernel);
        clhandle.createKernel(kernel2);

        clhandle.createIntBuffer(kernel, 0, esize, 0, 0);
        clhandle.setIntArg(kernel, 0, 2);
        clhandle.setIntArg(kernel, 1,esize);
        clhandle.copyIntBufferAcrossKernel(kernel, 0, kernel2, 0);
        clhandle.setIntArg(kernel2, 0, esize);
        //clhandle.setIntArg(kernel2, 1, block);
        clhandle.createIntBuffer(kernel2, 1, block, 0, 0);

        // The moving of the buffers into the device doesnt happen until the queue runs.
        clhandle.runKernel(kernel, gsize, lsize);

        clhandle.getIntBufferAsArray(kernel, 0, 3, true);            

        clhandle.runKernel(kernel, gsize, lsize);

        clhandle.getIntBufferAsArray(kernel, 0, 3, true);            

        clhandle.runKernel(kernel, gsize, lsize);

        clhandle.getIntBufferAsArray(kernel, 0, 3, true);            

        System.out.println("Divide array into "+block+" chunks");
        long time =System.nanoTime();
        int reducedSize = gsize;
        int reducedBy = maxReduce(reducedSize,clhandle);
        clhandle.runKernel(kernel2, reducedSize, reducedBy);
        
        reducedSize = reducedSize/reducedBy;
        System.out.println("Reduced Size | "+reducedSize+"    after reduced by : "+reducedBy);
        time =System.nanoTime()-time;
        
        int[] subsum=clhandle.getIntBufferAsArray(kernel2, 1, reducedSize, false);
        int sum = 0;
        for(int i =0;i<subsum.length;i++){
            sum+= subsum[i];
        }
        System.out.println("GPU Calculation in "+(time/1000)+"  us  with sum : " + sum +"    of "+esize+"  elements");
        
        
        
        time =System.nanoTime();
        int[] arr =clhandle.getIntBufferAsArray(kernel2, 0, esize, false);
        int sum2 = 0;
        for(int i =0;i<arr.length;i++){
            sum2+= arr[i];
        }
        time =System.nanoTime()-time;
        System.out.println("CPU Calculation in "+(time/1000)+"  us  with sum : " + sum2+"    of "+esize+"  elements");
        
        System.out.println("val theory: "+(esize-1)*6);

        // R = 22 is the point of returns for GPU method
        clhandle.closeOpenCL();
    }
    
    
    public static void testLocalSort3D(){
        String RNGKernelName = "mersenne_twister";   
        String kernel = "localSort3D";
        String kernelInit = "makeOrderIntBuffer";
        String kernelInter = "reinterpret1Das3D";
        CLHelper clHandler = new CLHelper();
        clHandler.initializeOpenCL("GPU");
        int mult = 6;
        int subL = 4;
        int L = 200;
        int dim1 = L/subL;
        int lsize = dim1;
        int chunk = (int)(Math.pow(subL, 3.0))*mult;        
        int gsize = (int)(Math.pow(L, 3.0))*mult;
        int nchunk = (gsize/chunk);
        Random ran = new Random();
        
        
        // Initialize order arr
        clHandler.createKernel(kernelInit);
        clHandler.createIntBuffer(kernelInit, 0, gsize, 0, 0);
        clHandler.setIntArg(kernelInit, 0, gsize);
        clHandler.setIntArg(kernelInit, 1, gsize);
        
        // sort setup 
        clHandler.createKernel(kernel);
        clHandler.createFloatBuffer(kernel, 0, gsize, 0.0f, "rw");
        clHandler.copyIntBufferAcrossKernel(kernelInit, 0, kernel, 0);
        clHandler.setIntArg(kernel, 0, subL);
        clHandler.setIntArg(kernel, 1, L);
        clHandler.setIntArg(kernel, 2, chunk);
        clHandler.setIntArg(kernel, 3, dim1);
        
        clHandler.runKernel(kernelInit, gsize, lsize);
        int[] order;// = clHandler.getIntBufferAsArray(kernel, 0, 2*chunk , true);
        
        System.out.println("Chunk Size | "+chunk+ "    nchunk | "+nchunk+"   gsize | "+ gsize);
        
        // throw random 
        RandomNumberCL RNG = new RandomNumberCL(ran.nextInt(),clHandler);
        int RNGKey = RNG.addRandomBuffer(kernel, 0, gsize);
        RNG.fillBufferWithRandom(RNGKey);
        
        // get values init
        System.out.println("Values initially before local sort");
        System.out.println("________________________________________________");
        //clHandler.getIntBufferAsArray(kernel, 0, lsize*lsize*lsize + 15, true);
        
        clHandler.runKernel(kernel, nchunk, lsize);
        order = clHandler.getIntBufferAsArray(kernel, 0, gsize , false);
        float[] values = clHandler.getFloatBufferAsArray(kernel, 0, gsize , false);
        checkChunk(order ,values, 0,1,0, chunk,subL, L, mult);
        
        /*
        int[] orderOld = clHandler.getIntBufferAsArray(kernel, 0, gsize , false);
        float[] valuesOld = clHandler.getFloatBufferAsArray(kernel, 0, gsize , false);
        for(int u = ((nchunk-2)*chunk); u < (nchunk*chunk);u++){
                    System.out.println("Position Ordered | "+(u) 
                            + "   PositionActual | "+ (orderOld[u])+ "    Value | "+valuesOld[orderOld[u]]
                            +"   diff | "+((u-orderOld[u])/chunk));    
        }
        // interpret as 3D 
        clHandler.createKernel(kernelInter);
        clHandler.copyFlBufferAcrossKernel(kernel, 0, kernelInter, 0);
        clHandler.copyIntBufferAcrossKernel(kernel, 0, kernelInter, 0);
        clHandler.createIntBuffer(kernelInter, 1, gsize, 0, "rw");
        clHandler.setIntArg(kernelInter, 0,subL);
        clHandler.setIntArg(kernelInter, 1, L);
        clHandler.createFloatBuffer(kernelInter, 1, 25, 0.0f, "rw");
        
        clHandler.runKernel(kernelInter, gsize, lsize);
    
        System.out.println("Values initially after local sort");
        System.out.println("________________________________________________");
        order = clHandler.getIntBufferAsArray(kernel, 0, gsize , false);
        float[] values = clHandler.getFloatBufferAsArray(kernel, 0, gsize , false);
        float[] err = clHandler.getFloatBufferAsArray(kernelInter, 1, 10 , true);
        checkChunk(order ,values, 0,1,0, chunk,subL, L, mult);
        */
    }
    
    
    public static void testLocalSort(){
        String RNGKernelName = "mersenne_twister";   
        String kernel = "localSort";
        String kernelInit = "makeOrderIntBuffer";
        String kernelInter = "reinterpret1Das3D";
        CLHelper clHandler = new CLHelper();
        clHandler.initializeOpenCL("GPU");
        int mult = 6;
        int lsize = 64;
        int subL = 2;
        int L = 12;
        int chunk = (int)(Math.pow(subL, 3.0))*mult;        
        int gsize = (int)(Math.pow(L, 3.0))*mult;
        int nchunk = (gsize/chunk);
        int dim1 = L/subL;
        Random ran = new Random();
        
        // Initialize order arr
        clHandler.createKernel(kernelInit);
        clHandler.createIntBuffer(kernelInit, 0, gsize, 0, 0);
        clHandler.setIntArg(kernelInit, 0, gsize);
        clHandler.setIntArg(kernelInit, 1, gsize);
        
        // sort setup 
        clHandler.createKernel(kernel);
        clHandler.createFloatBuffer(kernel, 0, gsize, 0.0f, "rw");
        clHandler.copyIntBufferAcrossKernel(kernelInit, 0, kernel, 0);
        clHandler.setIntArg(kernel, 0, chunk);
        clHandler.setIntArg(kernel, 1, nchunk);
        
        clHandler.runKernel(kernelInit, gsize, lsize);
        
        System.out.println("Chunk Size | "+chunk);
        
        // throw random 
        RandomNumberCL RNG = new RandomNumberCL(ran.nextInt(),clHandler);
        int RNGKey = RNG.addRandomBuffer(kernel, 0, gsize);
        RNG.fillBufferWithRandom(RNGKey);
        
        // get values init
        System.out.println("Values initially before local sort");
        System.out.println("________________________________________________");
        //clHandler.getIntBufferAsArray(kernel, 0, lsize*lsize*lsize + 15, true);
        
        clHandler.runKernel(kernel, nchunk, 1);
    }
    
    
    
    public static void checkChunk(int[] order ,float[] vals, int chunkX,int chunkY,int chunkZ,
            int chunk, int subL, int L, int mult){
        int x; int y; int z;int ind; int indSite = 0;
        for(int u = 0; u < (chunk/mult);u++){
            int subX = u%subL;
            int subY = (int)((double)u/(double)subL)%subL;
            int subZ = (int)((double)u/(double)(subL*subL))%subL;
            x = subX + chunkX*subL;
            y = subY + chunkY*subL;
            z = subZ + chunkZ*subL;
            indSite = (x+y*L+ z*(L*L))*mult;
            for(int v = 0; v < mult;v++){
                ind = indSite + v;       
                System.out.println("Position Ordered | "+(ind) 
                        + "   PositionActual | "+ (order[ind])+ "    Value | "+vals[order[ind]]
                        +"   mult | "+v+"   x | "+x+"   y | "+y+"   z | "+z );
            }
        }
    
    }
    
    public static void testInvasion(){
        CLHelper clHandler = new CLHelper();
        clHandler.initializeOpenCL("GPU");
        clHandler.createKernel("invperc3Dreduce");
            
    }

    public static void testOFCMin(){
        CLHelper clHandler = new CLHelper();
        clHandler.listAllDevices("");
        clHandler.initializeOpenCL("GPU");
        String loadCheck = "ofcLoadAndCheck";
        String dist = "ofcDistributeStress";
        String reduce = "reduceMinDiffFl";
        Random ran = new Random();
        clHandler.createKernel(loadCheck);
        clHandler.createKernel(dist);
        clHandler.createKernel(reduce);
        
        clHandler.setPrintMode(false);
        
        int size = 256*256;
        
        //  reduce min diff setup 
        clHandler.createFloatBuffer(reduce, 0, size, 0.0f ,  "rw");
        clHandler.createFloatBuffer(reduce, 1, size, 3.0f ,  "rw");
        clHandler.createFloatBuffer(reduce, 2, (size/2), 0.0f ,  "rw");
        clHandler.createFloatBuffer(reduce, 2, (size/2), 0.0f ,  "rw");
        clHandler.createIntBuffer(reduce, 0, (size/2), 0 ,  "rw");
        clHandler.setIntArg(reduce, 0, 1);
        clHandler.setIntArg(reduce, 1, size);
        
        System.out.println("+++++++++++++++++++++++++++++++++");
        // throw random 
        RandomNumberCL RNG = new RandomNumberCL(985758,clHandler);
        System.out.println("----------------------------------------");
        int RNGKey = RNG.addRandomBuffer(reduce, 0, size);
        System.out.println("****************************************");
        RNG.fillBufferWithRandom(RNGKey);
        System.out.println("+++++++++++++++++++++++++++++++++");

        int reducedSize = size;
        int reducedBy = maxReduce(reducedSize,clHandler);
        
        int[] loc = clHandler.getIntBufferAsArray(reduce, 0, (reducedSize/2), false);
        float[] buff1 = clHandler.getFloatBufferAsArray(reduce, 0, size, false);
        float[] buff2 = clHandler.getFloatBufferAsArray(reduce, 1, size, false);
        float min = 1000000.0f;
        for(int u = 0; u < reducedSize;u++){
            //System.out.println("i: "+u+"   val: "+(buff1[u]-buff2[u]));
            if( (buff1[u]-buff2[u]) < min ){
                min = (buff1[u]-buff2[u]);
            }
        }
        System.out.println("+++++++++++++++++++++++++++++++++");
        System.out.println("Minimum Value : "+min);
        System.out.println("Minimum Value : "+min);
        System.out.println("Minimum Value : "+min);
        System.out.println("Reduced By: "+reducedBy);
        System.out.println("+++++++++++++++++++++++++++++++++");
        
        clHandler.runKernel(reduce, reducedSize, reducedBy);
        clHandler.setIntArg(reduce, 0, 0);
        
        reducedSize = (reducedSize/reducedBy);
        clHandler.getBufferFlAsArrayList(reduce, 2, reducedSize, true);
        System.out.println("+++++++++++++++++++++++++++++++++");
        System.out.println("Minimum Value : "+min);
        System.out.println("+++++++++++++++++++++++++++++++++");
        
        
        reducedBy = maxReduce(reducedSize,clHandler);
        clHandler.runKernel(reduce, reducedSize, reducedBy);
        reducedSize = (reducedSize/reducedBy);
        clHandler.getBufferFlAsArrayList(reduce, 2, reducedSize, true);
        System.out.println("+++++++++++++++++++++++++++++++++");
        System.out.println("Minimum Value : "+min);
        System.out.println("+++++++++++++++++++++++++++++++++");
        
        reducedBy = maxReduce(reducedSize,clHandler);
        clHandler.runKernel(reduce, reducedSize, reducedBy);
        reducedSize = (reducedSize/reducedBy);
        clHandler.getBufferFlAsArrayList(reduce, 2, reducedSize, true);
        System.out.println("+++++++++++++++++++++++++++++++++");
        System.out.println("Minimum Value : "+min);
        System.out.println("+++++++++++++++++++++++++++++++++");
        
        
        
        int dim = 3;
        
        // distribute
        clHandler.copyFlBufferAcrossKernel(reduce, 0, dist, 0);
        clHandler.copyFlBufferAcrossKernel(reduce, 1, dist, 1);
        clHandler.createIntBuffer(dist, 0, size, 0 ,  "rw");
        clHandler.setFloatArg(dist, 0, 0.0f);
        clHandler.setIntArg(dist, 0, size);
        clHandler.setIntArg(dist, 1, size);
        clHandler.setIntArg(dist, 2, (int)(Math.pow(size, 1/dim)));
        clHandler.setIntArg(dist, 3, 1);
        clHandler.setIntArg(dist, 4, dim);
        clHandler.createFloatBuffer(dist, 2, 20, 0.0f ,  "rw");
        
        //load and check
        clHandler.copyFlBufferAcrossKernel(reduce, 0, loadCheck, 0);
        clHandler.copyFlBufferAcrossKernel(reduce, 1, loadCheck, 1);
        clHandler.copyIntBufferAcrossKernel(dist, 0, loadCheck, 0);
        clHandler.setFloatArg(loadCheck, 0, 1.0f);
        clHandler.setIntArg(loadCheck, 0, size);
        clHandler.copyFlBufferAcrossKernel(dist, 2, loadCheck, 2);
        
    }
    
    public static boolean isValThere(float[] values, float val){
        boolean there = false;
        for(int u = 0; u< values.length;u++){
            if(values[u] == val){
                System.out.println("VALUE IS THERE");
                there = true;
            }
        }
        return there;
    }

    private static int maxReduce(int gsize, CLHelper clhandler){      
        int max = clhandler.getCurrentDevice1DMaxWorkItems();
        if(gsize % 2 == 1){return 2;}
        int ind  = 2;
        int lsize = 2;
        while( Math.pow(2, ind) < max  ){
            // scratch space is 2056 for reduction
            if( (gsize % ((int)Math.pow(2, ind)) == 0) ){
                lsize = ((int)Math.pow(2, ind));
            }
            ind++;
        }
        return lsize;
    }
    
    public static void testRNG(){
        String RNGKernelName = "mersenne_twister";   
    
        CLHelper clHandler = new CLHelper();
        clHandler.listAllDevices("");
        Random ran = new Random(37380495);
        
        clHandler.initializeOpenCL("GPU");
        int GlobalWorkSize = 10000;
        int LocalWorkSize = 100;
    
        clHandler.createKernel("", RNGKernelName);
    
        // random numbers
        clHandler.createFloatBuffer(RNGKernelName, 0, GlobalWorkSize, 0.0f, 0);
     
        //seeds
        clHandler.setIntArg(RNGKernelName, 0, Math.abs(ran.nextInt()));
        clHandler.setIntArg(RNGKernelName, 1, Math.abs(ran.nextInt()));
        clHandler.setIntArg(RNGKernelName, 2, 2);
        
        clHandler.setKernelArg(RNGKernelName);

        for(int i=0;i<40;i++){
            clHandler.runKernel(RNGKernelName,GlobalWorkSize,LocalWorkSize);
            clHandler.getFloatBufferAsArray(RNGKernelName, 0, 2,true);
            clHandler.setIntArg(RNGKernelName, 0, Math.abs(ran.nextInt()));
            clHandler.setIntArg(RNGKernelName, 1, Math.abs(ran.nextInt()));
            clHandler.setKernelArg(RNGKernelName, true);
        }
    }
    
    
    public static void testRNGSimple(){
        String RNGKernelName = "simple_rng"; 
        CLHelper clHandler = new CLHelper();
        clHandler.listAllDevices("");
        clHandler.initializeOpenCL("GPU");
        int GlobalWorkSize = 100000000;
        int LocalWorkSize = 100;
        Random ran = new Random();
        
        
        clHandler.createKernel("", RNGKernelName);
    
        
        // random numbers
        clHandler.createFloatBuffer(RNGKernelName, 0, GlobalWorkSize, 0.0f, 0);
        //seeds
        //clHandler.createIntBuffer(RNGKernelName, 0, GlobalWorkSize, getSeedSystemsArray(GlobalWorkSize*2), 0, true);
        clHandler.setLongArg(RNGKernelName, 0,Math.abs(2*ran.nextLong()));
        clHandler.setLongArg(RNGKernelName, 1,Math.abs(2*ran.nextLong()));
        clHandler.setKernelArg(RNGKernelName);
        System.out.println("Using "+clHandler.getDeviceUsedMB(RNGKernelName)+" MB");
        for(int i=0;i<40;i++){
            clHandler.runKernel(RNGKernelName,GlobalWorkSize,LocalWorkSize);
            clHandler.getFloatBufferAsArray(RNGKernelName, 0, 2,true);
            clHandler.setLongArg(RNGKernelName, 0,Math.abs(2*ran.nextLong()));
            clHandler.setLongArg(RNGKernelName, 1,Math.abs(2*ran.nextLong()));
            clHandler.setKernelArg(RNGKernelName, true);
        }
    }
    
    public static long randomLong(long x) {
        x ^= (x << 21);
        x ^= (x >>> 35);
        x ^= (x << 4);
        return x;
    }
    private static int[] getSeedSystemsArray(int sizeMade){
        int[] seeds = new int[sizeMade];
        for(int i=0;i<seeds.length;i++){
        seeds[i] = 202394484;}//(int)(274252423*Math.random());}
        return seeds;
    }
    
    public static void main(String[] args) throws IOException {
            //testRNG();
            testIntReduce();
            //testLocalSort3D();
            //testInvasion();
            //testOFCMin();
            /*long x = 534534394;
            for(int i=0;i<10;i++){
                x= randomLong(x);
                double y = ((double)x)/(Math.pow(2.0,63.0)-1);
            System.out.println("random: "+y);}*/
            
    }
}

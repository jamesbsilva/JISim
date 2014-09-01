/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package GPUBackend;

import java.io.IOException;
import java.util.Random;
import scikit.opencl.CLHelper;
import scikit.opencl.RandomNumberCL;

/**
 *
 * @author James B Silva
 */
public class TestOpenCLHandler2 {
    
    
    public static void testMinFl(){
        CLHelper clhandler = new CLHelper();
        clhandler.initializeOpenCL("GPU");
        String minKernel = "minFl";
        clhandler.createKernel(minKernel);
        clhandler.createKernel("multiPassNgramStatAttack");
        int size = 500;
        int reduced = 100;
        long time;

        clhandler.createFloatBuffer(minKernel, 0, size, 0.0f, "rw");
        clhandler.createFloatBuffer(minKernel, 1, (2*reduced+1),0.0f, "rw");
        clhandler.setIntArg(minKernel, 0, reduced);
        clhandler.setIntArg(minKernel, 1, size);
        
        RandomNumberCL RNG = new RandomNumberCL(985758,clhandler);
        System.out.println("----------------------------------------");
        int RNGKey = RNG.addRandomBuffer(minKernel, 0, size);
        System.out.println("****************************************");
        RNG.fillBufferWithRandom(RNGKey);
        
        System.out.println("Running kernel with | "+reduced);
        time = System.nanoTime();
        clhandler.runKernel(minKernel, reduced, 1);
        printTimeDiff(time,"GPU ");
        
        float[] vals = clhandler.getFloatBufferAsArray(minKernel, 0, size, false);
        
        time = System.nanoTime();
        float min = 1000000;
        int minLoc = 1000000;
        float val;
        for(int u = 0; u < size;u++){
            val = vals[u];
            if(val < min){min = val; minLoc = u;}
        }
        printTimeDiff(time,"CPU ");
        

    }
    
    public static void printTimeDiff(long oldTime,String prefix){
        System.out.println(prefix+"time is "+(System.nanoTime() - oldTime)/1000+"  us");
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
            //testIntReduce();
            //testLocalSort3D();
            //testInvasion();
            testMinFl();
            /*long x = 534534394;
            for(int i=0;i<10;i++){
                x= randomLong(x);
                double y = ((double)x)/(Math.pow(2.0,63.0)-1);
            System.out.println("random: "+y);}*/
            
    }
}

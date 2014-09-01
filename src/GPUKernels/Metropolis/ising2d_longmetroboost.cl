/* 
        Metropolis algorithm implementation for Ising model with dilution and
    static spins. Spins are copied into a local buffer and work is split using
    spatial division and the sum of spins is kept in a buffer to speed up by
    exploiting the probabilistic need to check sum of spins and update. 


    James B. Silva <jbsilva@bu.edu>
*/
__kernel void ising2d_longmetroboost(__global int* spins,__global const int* fixed, __global const float* Params, __global float* random, __global int* dM, __global int* spinSums, __global float* dE, const int numThreads, const int nRuns, __global float* errorBuff) {
	// Get the index of the current element to be processed
	int tId = get_global_id(0);

	// bound check, equivalent to the limit on a 'for' loop
        if (tId >= numThreads){
            return;
        }

        // Get Params
	int numParams = 6;
	float jInteraction = Params[0];    
	float hField = Params[1];
	int L =	(int) Params[2];
	int subL = (int) Params[3];
	int range = (int) Params[4];
	float temp =  Params[5];

        int threads1D = (int)(L/subL);
        int ix = (tId)%(threads1D);
	int iy = (int) (tId/threads1D);
        iy = iy%threads1D;

	int index=0;
	int shift=0;

        int xSub;
	int ySub;
        int x0;
        int y0;
        int s0; 
        int i0; 
        int isCurrPos=0;
	int isSumPos=0;
        int mChange = 0;
        float eChange = 0;
        float delE;

        for(int currRun=0; currRun < nRuns; currRun++ ){
            // synchronize before moving on
            barrier(CLK_GLOBAL_MEM_FENCE);
            barrier(CLK_LOCAL_MEM_FENCE);
            // index for closer to sequential mem access.
            index = (tId+numThreads*2*currRun);
            shift = (int)(random[index]*subL*subL);

            xSub = (shift)%subL;
            ySub = (int) (shift/subL);
            ySub = ySub%subL;

            index = (x0+(y0*L));
            i0 = index;
            s0 = spins[index];

            // No need to do anything if fixed
            if (fixed[index] == 1)  {
                continue;
            }

            delE = 2*s0*(spinSums[i0]*jInteraction+hField);

            // Check to accept flip or not    	
            if((delE <= 0) || (exp(-1.0*delE/temp) > random[tId+(2*currRun+1)*numThreads])){
                    // accepted flip
                    
                    //update spin sums
                    x0 = xSub+subL*ix;
                    y0 = ySub+subL*iy;
                    int leftX = x0+L-range;
                    int leftY = y0+L-range;
                    int index = 0;
                    for(int v=0;v<(2*range+1);v++){
                        for(int u=0;u<(2*range+1);u++){
                            index = ((leftX+u)%L)+((leftY+v)%L)*L;
                            if(i0 != index){
                                spinSums[index] = spinSums[index]+2*(-1)*s0;
                            }
                        }
                    }                    

                    spins[i0] = (-1)*s0;
                    mChange = mChange - (2*s0);
                    eChange = eChange + delE;
            }else{}

        }
        
        dM[tId]=mChange;
        dE[tId]=eChange;
}

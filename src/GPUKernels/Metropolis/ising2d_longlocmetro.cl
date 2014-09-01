/* 
        Metropolis algorithm implementation for Ising model with dilution and
    static spins. Spins are copied into a local buffer and work is split using
    spatial division . 


    James B. Silva <jbsilva@bu.edu>
*/

__kernel void ising2d_longlocmetro(__global int* spins,__global const int* fixed, __global const float* Params, __global float* random, __global int* dM, __global float* dE, const int numThreads, const int nRuns, __global float* errorBuff) {
 
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
        float temp = Params[5];

        if(L > 1024){
            __local int spinsLoc[4194304];
            __local int fixedLoc[4194304];
        }else if(L > 512 && L < 1025){
            __local int spinsLoc[1048576];
            __local int fixedLoc[1048576];
        }else{
            __local int spinsLoc[262144];
            __local int fixedLoc[262144];
        }

        for(int u=0;u<L*L;u++){
            spinsLoc[u]= spins[u];
            fixedLoc[u]= fixed[u];
        }

        int threads1D = (int)(L/subL);
        int ix = (tId)%(threads1D);
	int iy = (int) (tId/threads1D);
        iy = iy%threads1D;

	int index=0;
	int shift=0;

        int sumNeigh = 0;
	int leftX = 0;
	int leftY = 0;
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

            index = (numThreads*nRuns+1+currRun);
            shift = (int)(random[index]*subL*subL);

            xSub = (shift)%subL;
            ySub = (int) (shift/subL);
            ySub = ySub%subL;

            x0 = xSub+subL*ix;
            y0 = ySub+subL*iy;

            index = (x0+(y0*L));
            i0 = index;
            s0 = spinsLoc[index];

            // No need to do anything if fixed
            if (fixedLoc[index] == 1)  {
                continue;
            }


            // Calculate the sum of neighbors	
            sumNeigh = 0;
            leftX = x0+L-range;
            leftY = y0+L-range;

            for(int u=0;u<(2*range+1);u++){
                for(int v=0;v<(2*range+1);v++){
                    index = ((leftX+u)%L)+((leftY+v)%L)*L;
                    if(i0==index){}else{
                    sumNeigh = sumNeigh + spinsLoc[index];}
                }
            }

            // get w index terms
            isCurrPos=0;
            isSumPos=0;
    \
            if(s0 >= 0){isCurrPos=1;}
            if(sumNeigh>=0){isSumPos=1;}

            index = (isCurrPos+(isSumPos*4)+abs(sumNeigh*4*2));

            delE = 2*s0*(sumNeigh*jInteraction+hField);

            // Check to accept flip or not    	
            if((delE<=0)||( (exp(-1.0*delE/temp) >random[tId+numThreads*currRun])){
                    // accepted flip
                    spinsLoc[i0] = (-1)*s0;
                    mChange = mChange - (2*s0);
                    eChange = eChange + delE;
            }else{}

        }
        
        dM[tId]=mChange;
        dE[tId]=eChange;
}

/*
        Metropolis algorithm implementation for nearest neighbor Ising model with dilution and
    static spins. This implementation uses spatial lattice decomposition to split work.

    James B. Silva <jbsilva@bu.edu>
*/
__kernel void ising2d_nnmetro(__global int* spins,__global const int* fixed, __global const float* Params, __global float* random, __global int* dM, __global float* dE, const int numThreads, __global float* errorBuff) {

        // Get the index of the current element to be processed
	int tId = get_global_id(0);

	// bound check, equivalent to the limit on a 'for' loop
        if (tId >= numThreads)  {
            return;
        }
        __local float wTransition[40];

        for(int u=0;u<40;u++){
            wTransition[u] = wMatrix[u];
        }

        // Get Params
	int numParams = 4;
	float jInteraction = Params[0];    
	float hField = Params[1];
	int L =	(int) Params[2];
	int subL = (int) Params[3];
	float temp = Params[4];

        // Calculate number of runs
        int nRuns = (int)(L*L/numThreads);

        int threads1D = (int)(L/subL);
        int ix = (tId)%(threads1D);
	int iy = (int) (tId/threads1D);
        iy=iy%threads1D;

        int index=0;
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
        int shift = 0;
        float delE;
        int mChange = 0;
        float eChange = 0;
        float acceptance;
        float toss;

        for(int currentRun = 0; currentRun<nRuns;currentRun++){
            // synchronize before moving on
            barrier(CLK_GLOBAL_MEM_FENCE);
            barrier(CLK_LOCAL_MEM_FENCE);

            // Get a set of random numbers
            index = (tId+numThreads*2*currRun);
            shift = (int)(random[index]*subL*subL);

            xSub = (shift)%subL;
            ySub = (int) (shift/subL);
            ySub = ySub%subL;

            x0 = xSub+subL*ix;
            y0 = ySub+subL*iy;

            i0 = x0+y0*L;
            s0 = spins[i0];

            // No need to do anything if fixed
            if (fixed[i0] == 1)  {
                continue;
            }

            // Calculate the sum of neighbors	
            sumNeigh = 0;
            leftX = x0-1+L;
            leftY = y0-1+L;

            index = ((leftX)%L)+(y0)*L;
            sumNeigh = sumNeigh +spins[index];

            index = ((leftX+2)%L)+(y0)*L;
            sumNeigh = sumNeigh +spins[index];

            index = (x0)+((leftY)%L)*L;
            sumNeigh = sumNeigh +spins[index];

            index = (x0)+((leftY+2)%L)*L;
            sumNeigh = sumNeigh +spins[index];

            // get w index terms
            isCurrPos=0;
            isSumPos=0;

            if(s0 >= 0){isCurrPos=1;}
            if(sumNeigh>=0){isSumPos=1;}

            index = (isCurrPos+isSumPos*4+abs(sumNeigh*4*2));

            delE = 2*s0*(sumNeigh*jInteraction+hField);

            acceptance = exp(-1.0*delE/temp);

            toss=random[tId+(2*currRun+1)*numThreads];

            // Check to accept flip or not    	
            if((delE<=0) || (acceptance > toss)){
                // accepted flip
                spins[i0] = (-1)*s0;
                mChange = mChange - (2*s0);
                eChange = eChange + delE;
            }else{}

        }

        dM[tId]=mChange;
        dE[tId]=eChange;
}


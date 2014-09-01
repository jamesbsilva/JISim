/* 
        Metropolis algorithm implementation for Ising model with dilution and
    static spins. This implementation performs parallel simulation of whole
    systems.

    James B. Silva <jbsilva@bu.edu>
*/
__kernel void ising2d_longmultmetro(__global int* spins,__global const int* fixed,__global int* spinSums, __global const float* jInteractionS ,__global const float* hFieldS,__global const float* tempS, __global float* random, __global int* M, __global float* E, const int numThreads, const int num, const int len, __global float* errorBuff) { 
	// Get the index of the current element to be processed
	int tId = get_global_id(0);

	// bound check, equivalent to the limit on a 'for' loop
        if (tId >= numThreads){
            return;
        }

        // Get Params
        int L = len;
        int N = num;
	int sysStart = N*tId;
        
        float jInteraction = jInteractionS[tId];    
	float hField = hFieldS[tId];
        float temp = tempS[tId];

	int index=0;

        int sumNeigh = 0;
	int leftX = 0;
	int leftY = 0;
        int x0;
        int y0;
        int s0; 
        int i0; 
        int isCurrPos=0;
	int isSumPos=0;
        int mChange = 0;
        float eChange = 0;
        float delE;

        for(int u = 0; u < (L*L); u++ ){
            index = (2*(u+N*tId));
            index = (int)(random[index]*L*L);
            i0 = index+sysStart;
            s0 = spins[i0];

            // No need to do anything if fixed
            if (fixed[i0] == 1){continue;}

            // Calculate the sum of neighbors	
            sumNeigh = sumSpins[i0];

            // get w index terms
            isCurrPos=0;
            isSumPos=0;
    
            if(s0 >= 0){isCurrPos=1;}
            if(sumNeigh>=0){isSumPos=1;}

            delE = 2*s0*(sumNeigh*jInteraction+hField);

            // Check to accept flip or not    	
            if( (delE <= 0)|| (exp(-1.0*delE/temp) > random[(2*(u+N*tId))+1])){
                // accepted flip
                spins[i0] = (-1)*s0;

                x0 = i0%L;
                y0 = ((int)(i0/L))%L;
                leftX = x0+L-range;
                leftY = y0+L-range;

                for(int u=0;u<(2*range+1);u++){
                    for(int v=0;v<(2*range+1);v++){
                        index = ((leftX+u)%L)+((leftY+v)%L)*L;
                        if(i0 != index){spinSums[index] = spinSums[index]+2*(-1)*s0;}
                    }
                }   
                mChange = mChange - (2*s0);
                eChange = eChange + delE;
            }
        }
        
        M[tId]+=mChange;
        E[tId]+=eChange;
}

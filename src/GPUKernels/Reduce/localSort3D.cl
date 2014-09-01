int getIndexFromChunkSubIndex(int i, int chunkIndex,int subL, int L, int bonds, int chunkL);

__kernel void localSort3D(__global float *bufferVals,__global int *orderLocal,
                                         int subL, int L, int chunk, int chunkL) {
    // Get the index of the current element to be processed
    int iX = get_global_id(0);
    int bonds = 6;
    // bound check, equivalent to the limit on a 'for' loop
    if (iX >= (chunkL*chunkL*chunkL))  {
        return;
    }

    // shell sort
    int increment = chunk / 2;
    while (increment > 0) {
            for (int i = increment; i < chunk; i++) {
                    int j = i;
                    float temp = bufferVals[orderLocal[getIndexFromChunkSubIndex(i,iX,subL,L,bonds,chunkL)]];
                    int temp2 = orderLocal[getIndexFromChunkSubIndex(i,iX,subL,L,bonds,chunkL)];
                    while (j >= increment && bufferVals[orderLocal[getIndexFromChunkSubIndex(j - increment,iX,subL,L,bonds,chunkL)]] > temp) {
                            orderLocal[ getIndexFromChunkSubIndex(j,iX,subL,L,bonds,chunkL) ] = orderLocal[getIndexFromChunkSubIndex(j - increment,iX,subL,L,bonds,chunkL)];
                            j = j - increment;
                    }
                    orderLocal[  getIndexFromChunkSubIndex(j,iX,subL,L,bonds,chunkL) ] = temp2;
            }
            if (increment == 2) {
                    increment = 1;
            } else {
                    increment *= (5.0 / 11);
            }
    }
}

int getIndexFromChunkSubIndex(int i, int chunkIndex,int subL, int L, int bonds, int chunkL){
    int site = (int)((float)i/(float)bonds);
    int bondDir = i % bonds;

    int chunkX = chunkIndex%chunkL;
    int chunkY = ((int)((float)chunkIndex/(float)chunkL) )%chunkL;
    int chunkZ = ((int)((float)chunkIndex/(float)(chunkL*chunkL)) )%chunkL;
    int subX = site % subL;
    int subY = ((int)((float)site/(float)subL) )% subL;
    int subZ = ((int)((float)site/(float)(subL*subL)) )% subL;
    
    return (bonds*((chunkX*subL + subX)+(chunkY*subL + subY)*L+(chunkZ*subL + subZ)*L*L)+bondDir);
}


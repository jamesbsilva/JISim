__kernel void localSort(__global float *bufferVals,__global int *orderLocal, int chunk, int numChunks) {
    // Get the index of the current element to be processed
    int iX = get_global_id(0);

    // bound check, equivalent to the limit on a 'for' loop
    if (iX >= numChunks)  {
        return;
    }
    int indexStart = iX*chunk;

    // shell sort
    int increment = chunk / 2;
    while (increment > 0) {
            for (int i = increment; i < chunk; i++) {
                    int j = i;
                    float temp = bufferVals[orderLocal[indexStart+i]];
                    int temp2 = orderLocal[indexStart+i];
                    while (j >= increment && bufferVals[orderLocal[indexStart + j - increment]] > temp) {
                            orderLocal[ indexStart + j ] = orderLocal[indexStart + j - increment];
                            j = j - increment;
                    }
                    orderLocal[ indexStart + j ] = temp2;
            }
            if (increment == 2) {
                    increment = 1;
            } else {
                    increment *= (5.0 / 11);
            }
    }
}




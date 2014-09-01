/* 
*        invperc3Dminprep .
*
*    James B. Silva <jbsilva@bu.edu>
*/
__kernel void invperc3Dminprep( __global int* localMinLoc, __global float* localMinVal, int numChunks, __global float* errorBuff) { 
    // Get the index of the current element to be processed
    int i = get_global_id(0);
    if( i >= numChunks){
        return;
    }
    localMinLoc[i] = -1;
    localMinVal[i] = 1.0;
}

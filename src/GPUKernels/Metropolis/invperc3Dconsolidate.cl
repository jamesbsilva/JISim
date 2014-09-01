/* 
*        invperc3Dconsolidate consolidates the minimum of bonds given that the  work size
*   divided by the reduction index parameter is bigger than the reduction index parameter.
*
*    James B. Silva <jbsilva@bu.edu>
*/
__kernel void invperc3Dconsolidate( __global int* localMinLoc, __global float* localMinVal, int reduceIndex, int numChunks, __global float* errorBuff) { 
    // Get the index of the current element to be processed
    int i = get_global_id(0);
    int actInd = (1+reduceIndex)*i;
    localMinLoc[i] = localMinLoc[actInd];
    localMinVal[i] = localMinVal[actInd];

}

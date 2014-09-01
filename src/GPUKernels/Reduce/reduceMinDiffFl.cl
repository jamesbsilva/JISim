/*
*              This is a reduce kernel for the minimum difference between 2 float buffers.
*
*
*      Author: James B. Silva (jbsilva@bu.edu) 
*
*/
__kernel void reduceMinDiffFl( __global float* buffer,__global float* buffer2,
             __global float* result,__global float* result2, __global int* resultLoc,
                int initialReduction, const int length) {
    __local float scratch[2056];
    __local int scratchLoc[2056];
    __local int scratch2[2056];
    int global_index = get_global_id(0);
    int local_index = get_local_id(0);
    // Load data into local memory
    if ( global_index < length ){
        if(initialReduction  == 1){
            scratch[local_index] = buffer[global_index]-buffer2[global_index];
            scratch2[local_index] = buffer[global_index];
            scratchLoc[local_index] = global_index;    
        }else{
            scratch[local_index] = buffer[resultLoc[global_index]]-buffer2[resultLoc[global_index]];
            scratch2[local_index] = buffer[resultLoc[global_index]];
            scratchLoc[local_index] = resultLoc[global_index];
        }
    } else {
        // Infinity is the identity element for the min operation
        scratch[local_index] = -1.0*INFINITY;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    for(int offset = 1; offset < get_local_size(0); offset <<= 1) {
        int mask = (offset << 1) - 1;
        if ((local_index & mask) == 0) {
            float other = scratch[local_index + offset];
            float mine = scratch[local_index];
            if(mine < other){
                scratch[local_index] =  mine ;
                scratch2[local_index] =  scratch2[local_index] ;
                scratchLoc[local_index] =  scratchLoc[local_index]  ;
            }else{
                scratch[local_index] =  other;
                scratchLoc[local_index] = scratchLoc[local_index+offset] ;
                scratch2[local_index] = scratch2[local_index+offset] ;
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (local_index == 0) {
        result[get_group_id(0)] = scratch[0];
        resultLoc[get_group_id(0)] = scratchLoc[0];
        result2[get_group_id(0)] = scratch2[0];
    }
}
/* 
        ofc2D_asper


    James B. Silva <jbsilva@bu.edu>
*/
__kernel void ofc2D_asper(__global float* spins, __global float* thresh , __global float* random, __global float* errorBuff) {
    // Get the index of the current element to be processed
    int tId = get_global_id(0);
    // bound check, equivalent to the limit on a 'for' loop
    if (tId >= numThreads){
        return;
    }
    
}

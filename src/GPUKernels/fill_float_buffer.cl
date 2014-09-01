/*
 *              This kernel fills a buffer with value s0.
 *
 *
 *      Author: James B. Silva (jbsilva@bu.edu) 
 *
*/
__kernel void fill_float_buffer(__global float *A, const float s0) {
    // Get the index of the current element to be processed
    int i = get_global_id(0);
    
    // Do the operation
    A[i] = s0;
}


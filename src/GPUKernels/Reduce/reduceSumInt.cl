/*
 *              This is a reduce kernel which partially sums up an array.
 *
 *
 *      Author: James B. Silva (jbsilva@bu.edu) 
 *
*/
__kernel void reduceSumInt( __global int* buffer, const int length, __global int* result) {
    __local int scratch[2056];         
    int global_index = get_global_id(0);
    int local_index = get_local_id(0);
    // Load data into local memory
    if (global_index < length) {
        scratch[local_index] = buffer[global_index];
    } else {
        // ignore element for the sum operation
        scratch[local_index] = 0;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    for(int offset = 1; offset < get_local_size(0);  offset <<= 1) {
        int mask = (offset << 1) - 1;
        if ((local_index & mask) == 0) {
            int other = scratch[local_index + offset];
            int mine = scratch[local_index];
            scratch[local_index] =  mine + other;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if (local_index == 0) {
        result[get_group_id(0)] = scratch[0];
    }
}
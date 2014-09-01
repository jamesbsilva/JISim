__kernel void makeOrderIntBuffer(__global int *A, int chunk, int numElements) { 
    // Get the index of the current element to be processed
    int i = get_global_id(0);
    
    // bound check (equivalent to the limit on a 'for' loop for standard/serial C code
    if ( i == numElements) {
        return;
    } 

    // Do the operation
    A[i] = i%chunk;
}

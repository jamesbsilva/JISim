__kernel void iterFl(__global float *A, float iter, int numElements) {
 
    // Get the index of the current element to be processed
    int i = get_global_id(0);
    
    // bound check, equivalent to the limit on a 'for' loop
    if (i >= numElements)  {
        return;
    }

    // Do the operation
    A[i] = A[i] + iter;
}

__kernel void multFl(__global float *A, float mult, int numElements) {
 
    // Get the index of the current element to be processed
    int i = get_global_id(0);
    
    // bound check, equivalent to the limit on a 'for' loop
    if (i >= numElements)  {
        return;
    }

    // Do the operation
    A[i] = A[i]*mult;
}

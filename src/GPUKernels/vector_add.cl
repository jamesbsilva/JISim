/*
 *              Kernel that adds two buffers 
 *
*/
__kernel void vector_add(global const float* a, global const float* b, global float* c, int numElements) {

        // get index into global data array
        int iGID = get_global_id(0);

        // bound check, equivalent to the limit on a 'for' loop
        if (iGID >= numElements)  {
            return;
        }

        // add the vector elements
        c[iGID] = a[iGID] + b[iGID];
    }
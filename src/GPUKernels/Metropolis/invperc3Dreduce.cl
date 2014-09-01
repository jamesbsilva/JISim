/* 
*        invpercolation3D  performs an invasion percolation step after a local sorting
*   using a localSort3D of the bonds.
*
*    James B. Silva <jbsilva@bu.edu>
*/
__kernel void invperc3Dreduce( __global int* system, __global float* bonds, 
              __global int* localMinLoc, __global float* localMinVal, int reduceIndex, int numChunks, __global float* errorBuff) { 
    // Get the index of the current element to be processed
    int i = get_global_id(0);
    int actInd = (1+reduceIndex)*i;
    int redInd = actInd + reduceIndex;
    
    // bound check, equivalent to the limit on a 'for' loop
    if (actInd >= numChunks){
        return;
    }
    // bound check, equivalent to the limit on a 'for' loop
    if (redInd >= numChunks){
        return;
    }

    // if other is not in fringe do not bother to  compare
    if( localMinLoc[redInd] < 0 ){
        localMinVal[redInd] = 1.0f;
        if( localMinLoc[actInd] < 0 ){
            localMinVal[actInd] = 1.0f;
        }else{
            localMinVal[actInd] = bonds[localMinLoc[actInd]];
        }
        return;
    }

    // if this is not in fringe do not bother to compare
    if( localMinLoc[actInd] < 0 ){
        localMinLoc[actInd] = localMinLoc[redInd];
        localMinVal[actInd] = bonds[localMinLoc[redInd]];
        localMinVal[redInd] = 1.0; 
        localMinLoc[redInd] = -1;
        return;
    }
    
    // compare
    if( bonds[localMinLoc[redInd]] < bonds[localMinLoc[actInd]]  ){
        localMinLoc[actInd] = localMinLoc[redInd]; 
        localMinVal[actInd] = bonds[localMinLoc[redInd]];
        localMinVal[redInd] = 1.0;
        localMinLoc[redInd] = -1; 
        return;
    }else{    
        localMinVal[actInd] = bonds[localMinLoc[actInd]];
        localMinVal[redInd] = 1.0;
        localMinLoc[redInd] = -1;  
    }
}

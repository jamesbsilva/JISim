/*
        This kernel colors a lattice based on the latCol buffer.

    James B. Silva <jbsilva@bu.edu>
*/

__kernel void colorLatticeInt2(__global int *lat,__global float4 *latCol){
    unsigned int i = get_global_id(0);
    // color
    if( lat[i] == 1){
        latCol[i].x = 1.0f;
        latCol[i].y = 1.0f;
        latCol[i].z = 1.0f;
        latCol[i].w = 1.0f;
    }else if( lat[i] == (-1) ){
        latCol[i].x = 1.0f;
        latCol[i].y = 1.0f;
        latCol[i].z = 0.0f;
    }else if(lat[i] == 0){
        latCol[i].x = 0.0f;
        latCol[i].y = 0.0f;
        latCol[i].z = 0.0f;
        latCol[i].w = 0.0f;
    }

    return;
}

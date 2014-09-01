/*
        This kernel draws a 3D BCC lattice

    James B. Silva <jbsilva@bu.edu>
*/

__kernel void createLatticeBCC(__global float4 * latPos,__global float4 *latCol, int L, int geo, float scale){
    unsigned int i = get_global_id(0);

    // calculate uv coordinates
    
    
    return;
}

float4 getPosCubic(int i, int j, int k, float scale, int L){
    return (float4)(i/(L*scale), k/(L*scale), j/(L*scale), 1.0f);
}

float4 getPosCenter(int i, int j, int k, float scale, int L, int face){
    return (float4)((i+0.5f)/(L*scale), (k+0.5f)/(L*scale), (j+0.5f)/(L*scale), 1.0f);    
}

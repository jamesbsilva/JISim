/*
        This kernel draws a 3D FCC lattice

    James B. Silva <jbsilva@bu.edu>
*/
__kernel void createLatticeFCC(__global float4 * latPos,__global float4 *latCol, int L, int geo, float scale){
    unsigned int i = get_global_id(0);

    // calculate uv coordinates
    
    
    return;
}

float4 getPosCubic(int i, int j, int k, float scale, int L){
    return (float4)(i/(L*scale), k/(L*scale), j/(L*scale), 1.0f);
}

float4 getPosFaceBottomCorner(int i, int j, int k, float scale, int L, int face){
    // top face
    if(face ==1){
        return (float4)((i+0.5f)/(L*scale), k/(L*scale), (j+0.5f)/(L*scale), 1.0f);
    }else if(face == 2){
        return (float4)((i+0.5f)/(L*scale), (k+0.5f)/(L*scale), j/(L*scale), 1.0f);
    }else if(face == 3){
        return (float4)(i/(L*scale), (k+0.5f)/(L*scale), (j+0.5f)/(L*scale), 1.0f);
    }    
}

float4 getPosFaceUpperCorner(int i, int j, int k, float scale, int L, int face){
     // top face
    if(face ==1){
        return (float4)((i-0.5f)/(L*scale), k/(L*scale), (j-0.5f)/(L*scale), 1.0f);
    }else if(face == 2){
        return (float4)((i-0.5f)/(L*scale), (k-0.5f)/(L*scale), j/(L*scale), 1.0f);
    }else if(face == 3){
        return (float4)(i/(L*scale), (k-0.5f)/(L*scale), (j-0.5f)/(L*scale), 1.0f);
    }    
}
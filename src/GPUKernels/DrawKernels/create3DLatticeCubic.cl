/*
        This kernel draws a 3D cubic lattice

    James B. Silva <jbsilva@bu.edu>
*/

float4 getPosCubic(int i, int j, int k, float scale, int Lx, int Ly, int Lz);

__kernel void create3DLatticeCubic(__global float4 * latPos,__global float4 *latCol, int Lx, int Ly, int Lz, float scale){
    unsigned int tId = get_global_id(0);
    
    int i = (tId%Lx);
    int j = ((int)((double)tId/(double)Lx))%Ly;
    int k = ((int)((double)tId/(double)(Lx*Ly)))%Lz;

    // calculate uv coordinates
    latPos[tId] = getPosCubic(i,j,k,scale,Lx,Ly,Lz);
    
    return;
}

float4 getPosCubic(int i, int j, int k, float scale, int Lx, int Ly, int Lz){
    return (float4)((i-(Lx/2.0))/(Lx*scale), (k-(Lz/2.0))/(Lz*scale), (j-(Ly/2.0))/(Ly*scale), 1.0f);
}


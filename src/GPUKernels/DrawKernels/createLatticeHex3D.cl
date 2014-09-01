/*
        This kernel draws a 3D hexagonal lattice

    James B. Silva <jbsilva@bu.edu>
*/

float getXrCentered(int u,int v,int L);
float getYrCentered(int u,int v,int L);
float getXcoordTri(int i, int j);
float getYcoordTri(int i, int j);

__kernel void createLatticeHex3D(__global float4 * latPos,__global float4 *latCol, int L, float scale){
    unsigned int i = get_global_id(0);

    // calculate uv coordinates
    int u = (i%L);
    int v = ((i/L)%L);
    int w = ((i/(L*L))%L);
    float x = getXrCentered(u,v,L)/(L*scale);
    float y = getYrCentered(u,v,L)/(L*scale);
    float z = getZrCentered(u,v,L)/(L*scale);

    // write output vertex
    latPos[u+v*L+w*(L*L)] = (float4)(x, z, y, 1.0f);
    latCol[u+v*L+w*(L*L)] = (float4)(0.0f, 1.0f, 1.0f, 255.0f);
    return;
}

float getXrCentered(int u,int v,int L){
    int cenInd = (int)(L/2.0);
    float xcen = getXcoordTri(cenInd,cenInd);
    return (getXr(u,v)-xcen);
}


float getYrCentered(int u,int v, int L){
    int cenInd = (int)(L/2.0);
    float ycen = getYcoordTri(cenInd,cenInd);
    return (getYr(u,v)-ycen);
}

float getZrCentered(int v,int L){
    int cenInd = (int)(L/2.0);
    float zcen = cenInd;
    return (v-zcen);
}

float getXcoordTri(int i, int j){
    return (j%2 == 0) ? (float)i : (float)i + 0.5;
}

float getYcoordTri(int i, int j){
    return ((float)j)*sqrt(3.0f)/2.0;
}

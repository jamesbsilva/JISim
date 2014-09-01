/*
        This kernel draws a 2d lattice of 3 types (tri/square/honey)

    James B. Silva <jbsilva@bu.edu>
*/

float getXrCentered(int u,int v,int Geo,int L);
float getYrCentered(int u,int v,int Geo,int L);
float getXr(int u,int v,int Geo);
float getYr(int u,int v,int Geo);
float getXcoordTri(int i, int j);
float getYcoordTri(int i, int j);
float getXcoordHoney(int i, int j);
float getYcoordHoney(int i, int j);

__kernel void createLattice2D(__global float4 * latPos,__global float4 *latCol, int L, int geo, float scale){
    unsigned int i = get_global_id(0);

    // calculate uv coordinates
    int u = (i%L);
    int v = ((i/L)%L);
    float x = getXrCentered(u,v,geo,L)/(L*scale);
    float y = getYrCentered(u,v,geo,L)/(L*scale);

    // write output vertex
    latPos[u+v*L] = (float4)(x, 0.0f, y, 1.0f);
    latCol[u+v*L] = (float4)(0.0f, 1.0f, 1.0f, 255.0f);
    return;
}

float getXrCentered(int u,int v,int Geo,int L){
    int cenInd = (int)(L/2.0);
    float xcen = getXr(cenInd,cenInd,Geo);
    return (getXr(u,v,Geo)-xcen);
}

float getYrCentered(int u,int v,int Geo, int L){
    int cenInd = (int)(L/2.0);
    float ycen = getYr(cenInd,cenInd,Geo);
    return (getYr(u,v,Geo)-ycen);
}

float getXr(int u,int v,int Geo){
    if(Geo ==  6){
        return getXcoordTri(u,v); 
    }else if(Geo ==  3){
        return getXcoordHoney(u,v); 
    }else{
        return u;
    }
}

float getYr(int u,int v, int Geo){
    if(Geo ==  6){
        return getYcoordTri(u,v); 
    }else if(Geo ==  3){
        return getYcoordHoney(u,v); 
    }else{
        return v;
    }
}

float getXcoordTri(int i, int j){
    return (j%2 == 0) ? (float)i : (float)i + 0.5;
}

float getYcoordTri(int i, int j){
    return ((float)j)*sqrt(3.0f)/2.0;
}

float getXcoordHoney(int i, int j){
    float offset = (j%2 == 0) ? floor((float)(i/2)) : (floor((float)(i/2))+1.5);
    return (((float)i))+offset;
}

float getYcoordHoney(int i, int j){
    return ((float)j)*sqrt(3.0f)/2.0;
}

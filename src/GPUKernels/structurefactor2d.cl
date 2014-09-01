/*
 *              structurefactor2d calculates the structure transform of a LxL lattice of integer spins
 *      interpreted as having a lattice type determined by the geo integer parameter
 *      (6 : triangular, 3 : hexagonal, else square)
 *
 *
 *      Author: James B. Silva (jbsilva@bu.edu) 
 *      Date: 06/29/2013
 *
*/

// forward declarations 
float dotproduct(float a,float b, float c, float d);
float getXrCentered(int u,int v,int Geo,int L);
float getYrCentered(int u,int v,int Geo,int L);
float getXr(int u,int v,int Geo);
float getYr(int u,int v,int Geo);
float getXcoordTri(int i, int j);
float getYcoordTri(int i, int j);
float getXcoordHoney(int i, int j);
float getYcoordHoney(int i, int j);

__kernel void structurefactor2d(__global int* spins, __global float* sFactor, const int Lq, const int L, const int geo, const float PI, __global float* errorBuff) {
    // Get the index of the current element to be processed
    int tId = get_global_id(0);

    // bound check, equivalent to the limit on a 'for' loop
    if (tId >= (int)(Lq*Lq)){
        return;
    }
 
    float q1x = 1;
    float q1y = 0;
    float q2x = 0;
    float q2y = 1;
    float i = (-1.0*Lq/2.0f)+(tId%Lq);
    float j = (-1.0*Lq/2.0f)+(((int)(float)tId/Lq)%Lq);
    float spin;
            
    float qfactor = (2*PI/L);
    float reSum = 0.0;
    float compSum = 0.0;
    float argX = 0.0;
    float argY = 0.0;
    float arg = 0.0;

    for(int v = 0; v < L;v++){for(int u = 0; u < L;u++){
            // calculate components for dot product
            // G1*nx+G2*nx where G1 and G2 are the basis vectors 
            argX = i*q1x+j*q2x;
            argY = i*q1y+j*q2y;
            // calculate distance in functions getXr and getYr to allow for more complicated lattices
            arg = dotproduct(argX,argY,getXrCentered(u,v,geo,L),getYrCentered(u,v,geo,L));
            // qfactor is 2 pi/L
            arg = arg*qfactor;
            //exp(iqr)*s
            spin = spins[u+v*L];
            reSum = reSum + cos(arg)*spin;
            compSum = compSum + sin(arg)*spin;
    }}
    // |e(iqr)|^2/N
    arg = (reSum*reSum+compSum*compSum);
    arg = arg/(L*L);
    sFactor[tId] = arg; 
}

float dotproduct(float a,float b, float c, float d){
    return (a*c+b*d);
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

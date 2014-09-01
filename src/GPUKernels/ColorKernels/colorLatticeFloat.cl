/*
        This kernel colors a lattice based on the latCol buffer.

    James B. Silva <jbsilva@bu.edu>
*/

float getColorX(float val, float max, float min);
float getColorY(float val, float max, float min);
float getColorZ(float val, float max, float min);

__kernel void colorLatticeFloat(__global float *lat,__global float4 *latCol, float max, float min){
    unsigned int i = get_global_id(0);
    bool drawRangeStyle =  true;
    
    // color
    if(drawRangeStyle){
        latCol[i].x = getColorX(lat[i],max,min);
        latCol[i].y = getColorY(lat[i],max,min);
        latCol[i].z = getColorZ(lat[i],max,min);
    }else{
        latCol[i].x = 1.0f;
        latCol[i].y = (lat[i]-min)/(max-min);
        latCol[i].z = (lat[i]-min)/(max-min);
    }
    
    latCol[i].w = (lat[i]-min)/(max-min);

    return;
}

float getColorX(float val, float max, float min){
    float col = (val-min)/(max-min);
    float seg = 7.0f;
    int colScale = (int)(floor(col*seg));
    col  = col*seg-(float)colScale;     
    if(colScale == 2){
        return col;
    }else if(colScale == 1 || colScale == 6){
        return 0.0f;
    }else if(colScale == 0 || colScale == 5){
        return 1.0f-col;
    }else{
        return 1.0f;
    }        
}
float getColorY(float val, float max, float min){
    float col = (val-min)/(max-min);
    float seg = 7.0f;
    int colScale = (int)(floor(col*seg));
    col  = col*seg-(float)colScale;     
    if(colScale == 4){
        return col;
    }else if(colScale == 0 || colScale == 5){
        return 1.0f;
    }else if(colScale == 1 || colScale == 6){
        return 1.0f-col;
    }else{
        return 0.0f;
    }        
}

float getColorZ(float val, float max, float min){
    float col = (val-min)/(max-min);
    float seg = 7.0f;
    int colScale = (int)(floor(col*seg));
    col  = col*seg-(float)colScale;     
    if(colScale == 2 || colScale == 1 || colScale == 0){
        return 1.0f;
    }else if(colScale == 3){
        return 1.0f-col;
    }else{
        return 0.0f;
    }        
}
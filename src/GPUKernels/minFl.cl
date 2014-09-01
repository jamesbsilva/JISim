/*
 *              This is a reduce kernel for the minimum of a float array.
 *
 *      Author: James B. Silva (jbsilva@bu.edu) 
 *
*/
__kernel void minFl( __global float* buffer, __global float* minInfo, const int length, const int chunks) {
    int global_index = get_global_id(0);

    int reduceNum = (int) minInfo[0];
    int i = reduceNum*global_index;
    
    if(i > length){return;}

    if ( (i + reduceNum) > length ){reduceNum = length-i;}
    
    float min = INFINITY;
    int minLoc = 0;
    float val;
    for(int u = 0;u < reduceNum;u++){
        val = buffer[u+i];
        if( val < min ){min = val;minLoc = (u+i);}
    }
    minInfo[1+global_index*2] = min;
    minInfo[2+global_index*2] = (float) minLoc;
    
    barrier(CLK_LOCAL_MEM_FENCE);
    barrier(CLK_GLOBAL_MEM_FENCE);  
   
    if (i == 0) {
        min = INFINITY;
        for( int u = 0; u < chunks;u++){
            val = minInfo[1+u*2];
            if( val < min){
                min = val;
                minLoc = minInfo[2+u*2];
            }
        }
        minInfo[1] = min;
        minInfo[2] = (float)minLoc;
        return;
    }else{
        return;
    }
}

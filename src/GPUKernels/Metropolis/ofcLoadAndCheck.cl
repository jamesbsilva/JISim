/* 
*        ofcLoadAndCheck is a function that loads the system and checks if location is over threshold. 
*
*
*    James B. Silva <jbsilva@bu.edu>
*/
__kernel void ofcLoadAndCheck(__global float* loc, __global float*  thresh,
             __global int*  overThresh, float loadStress, int N, __global float* errorBuff) { 
    // Get the index of the current element to be processed
    int tId = get_global_id(0);

    // bound check, equivalent to the limit on a 'for' loop
    if (tId >= N){
        return;
    }

    loc[tId] = loc[tId]+loadStress;

    if( loc[tId] > thresh[tId] ){
        overThresh[tId] = 1;
    }else{
        overThresh[tId] = 0;
    } 
}
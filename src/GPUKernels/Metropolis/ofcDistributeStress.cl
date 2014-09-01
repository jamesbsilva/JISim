/* 
*        ofcDistributeStress distributes stress across the system with the given interaction range.
*   Ideal for interactions with above 64 points;
*
*    James B. Silva <jbsilva@bu.edu>
*/
__kernel void ofcDistributeStress(__global float* sys, __global float* thresh , __global int* overThresh ,
                    float newVal, float stressDist, int maxLoc,  const int N,
                    const int L, const int R, const int dim, __global float* errorBuff) {
    // Get the index of the current element to be processed
    int tId = get_global_id(0);
    // bound check, equivalent to the limit on a 'for' loop
    if (tId >= N){
        return;
    }
    float diffOut = (stressDist)/((float)pow((float)(2*R+1),dim)-1);
    
    // calculate x y z coordinates of stress max location
    int i = maxLoc%L;
    int j = (dim > 1) ? (maxLoc/L)%L : 0;
    int k = (dim > 2) ? (maxLoc/(L*L))%L : 0;

    // calculate x y z coordinates of stress distribute location
    int a = tId % (2*R+1);
    int b = (dim > 1 ) ? (tId/(2*R+1)) % (2*R+1)  : R;
    int c = (dim > 2 ) ? (tId/((2*R+1)*(2*R+1))) % (2*R+1) : R;
    
    // index of updated stress location
    int ind = ((i+L-R+a)%L) + ((j+L-R+b)%L)*L + ((k+L-R+c)%L)*(L*L);

    if( ind != maxLoc){
        sys[ind] = sys[ind] + diffOut;
        if(sys[ind] > thresh[ind]){overThresh[ind] = 1;}else{overThresh[ind] = 0;}
    }
    
    // update max stress location
    if(tId == 0){
        sys[maxLoc] = newVal;
        if(sys[maxLoc] > thresh[maxLoc]){overThresh[maxLoc] = 1;}else{overThresh[maxLoc] = 0;}
    }
}

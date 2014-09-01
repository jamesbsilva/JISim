/*
 *        This kernel colors a lattice based on the latCol buffer.
 *
 *    James B. Silva <jbsilva@bu.edu>
 */

int getXfromInd(int ind, int L);int getYfromInd(int ind, int L);int getZfromInd(int ind, int L);
int getXfromChunk(int ind, int subL, int L);int getYfromChunk(int ind, int subL, int L);int getZfromChunk(int ind, int subL, int L);

__kernel void boxCounting( __global int *lat, __global int *boxOccupations, int L,
                 int N, int subL,int dim, int divisions){
    unsigned int i = get_global_id(0);
    int x0;int y0;int z0;int x;int y;int z;
    int xs0;int ys0;int zs0;
    
    int boxOccu = 0;int boxInd = 0;
    int currSubL = subL;int currL = L;
    int cL = 0;int j = 0; int n = 0;

    x0 = getXfromChunk(i, subL, L );y0 = getYfromChunk(i, subL, L );
    z0 = getZfromChunk(i, subL, L );
        
    for(int u = divisions; u > 0;u--){
        currSubL = (int)(subL*pow(0.5f,u));cL = (int)((float)subL/(float)currSubL);
        n = (int)(pown((float)cL,dim));
        for(int k = 0; k < n;k++){
            xs0 = getXfromChunk(k, currSubL, subL);ys0 = getYfromChunk(k, currSubL, subL);zs0 = getZfromChunk(k, currSubL, subL);          
            j = -1;boxOccu = 0;
            while(j < ((int)pown((float)currSubL,dim)) && boxOccu == 0 ){
                j++;
                // is box occupied
                x = getXfromInd(j,currSubL);y = getYfromInd(j,currSubL);z = getZfromInd(j,currSubL);
                x = x+xs0+x0;y = y+ys0+y0;z = z+zs0+z0;
                if(lat[x+y*L+z*L*L] == 1){boxOccu = 1;}
            }
            boxOccupations[boxInd] = boxOccu;
            boxInd++;
        }    
    }
    return;
}

int getXfromInd(int ind, int L){
    return (ind%L);
}
int getYfromInd(int ind, int L){
    return ((int)((float)ind/(float)L)%L);
}
int getZfromInd(int ind, int L){
    return ((int)((float)ind/(float)(L*L))%L);
}
int getXfromChunk(int ind, int subL, int L){
    return ((getXfromInd(ind,(int)((float)L/(float)subL)))*subL);
}
int getYfromChunk(int ind, int subL, int L){
    return ((getYfromInd(ind,(int)((float)L/(float)subL)))*subL);
}
int getZfromChunk(int ind, int subL, int L){
    return ((getZfromInd(ind,(int)((float)L/(float)subL)))*subL);
}
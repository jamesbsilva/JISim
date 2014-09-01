/* 
*        invpercolation3D  performs an invasion percolation step after a local sorting
*   using a localSort3D of the bonds.
*
*    James B. Silva <jbsilva@bu.edu>
*/

int getLocationToFromBond3D(int bond, int locX, int locY, int locZ,int L);
int getIndexFromChunkSubIndex(int i, int chunkX, int chunkY, int chunkZ, int subL, int L, int bonds, int chunkL);
int getXfromBondInd(int ind, int L, int nbonds);
int getYfromBondInd(int ind, int L, int nbonds);
int getZfromBondInd(int ind, int L, int nbonds);

__kernel void invperc3Dmin( __global int* system, __global float* bonds, 
        __global int* bondsOrder, __global float* fringe, __global int* localMinLoc, __global int* paramInt,
                 int numThreads, __global float* errorBuff) { 
    // Get the index of the current element to be processed
    int chunkX = get_global_id(0);
    int chunkY = get_global_id(1);
    int chunkZ = get_global_id(2);
    
    int L = paramInt[0];
    int subL = paramInt[1];
    int chunk = paramInt[2];
    int chunkL = paramInt[3];
    int nbonds = 6;
    int tId = chunkX+chunkY*chunkL+chunkZ*(chunkL*chunkL);

    // bound check, equivalent to the limit on a 'for' loop
    if (tId >= numThreads){
        return;
    }

    // looking for a local minimum that is not already occupied
    int searchingGood = 1;
    int ind = 0;
    int sysLoc;
    int actInd = 0;

    while( (searchingGood == 1) && (ind < chunk) ){
        actInd = getIndexFromChunkSubIndex(ind,chunkX,chunkY,chunkZ,subL,L,nbonds,chunkL);
        if(bondsOrder[actInd] > 0){
            sysLoc = getLocationToFromBond3D(bondsOrder[actInd], getXfromBondInd(actInd,L,nbonds),
                     getYfromBondInd(actInd,L,nbonds), getZfromBondInd(actInd,L,nbonds), L);
            if( system[sysLoc] == 1 ){
                bondsOrder[actInd] = -1;
                fringe[actInd] = 0;
            }else{
                if( fringe[actInd] == 1 ){
                    searchingGood = 0;
                    localMinLoc[tId] = bondsOrder[actInd];
                }
            }
        }
        ind++;
    }
    // pick highest val if not in fringe
    if( localMinLoc[tId] == 0 ){
        localMinLoc[tId] = bondsOrder[getIndexFromChunkSubIndex(chunk-1,chunkX,chunkY,chunkZ,subL,L,nbonds,chunkL)];
    }
}

int getXfromBondInd(int ind, int L,int nbonds){
    return ((int)((float)ind/(float)nbonds)%L);
}
int getYfromBondInd(int ind, int L,int nbonds){
    return ((int)((float)ind/(float)(L*nbonds))%L);
}
int getZfromBondInd(int ind, int L,int nbonds){
    return ((int)((float)ind/(float)(L*L*nbonds))%L);
}

int getLocationToFromBond3D(int bond, int locX, int locY, int locZ,int L){
    int bondDir = bond % 6;
    int toLoc = 0;
    switch( bondDir ){
        case 0:
            toLoc = locX + ((locY+1)%L)*L + locZ*(L*L);
        case 1:
            toLoc = ((locX+1)%L) + locY*L + locZ*(L*L);
        case 2:
            toLoc = locX + ((locY+L-1)%L)*L + locZ*(L*L);
        case 3:
            toLoc = ((locY+L-1)%L) + locY*L + locZ*(L*L);
        case 4:
            toLoc = locX + locY*L + ((locZ+1)%L)*(L*L);
        case 5:
            toLoc = locX + locY*L + ((locZ+L-1)%L)*(L*L);
        default :
            toLoc = locX + locY*L + locZ*(L*L);
    }
    return toLoc;
}

int getIndexFromChunkSubIndex(int i, int chunkX, int chunkY, int chunkZ, int subL, int L, int bonds, int chunkL){
    int site = (int)((float)i/(float)bonds);
    int bondDir = i % bonds;

    int subX = site % subL;
    int subY = ((int)((float)site/(float)subL) )% subL;
    int subZ = ((int)((float)site/(float)(subL*subL)) )% subL;
    
    return (bonds*((chunkX*subL + subX)+(chunkY*subL + subY)*L+(chunkZ*subL + subZ)*L*L)+bondDir);
}


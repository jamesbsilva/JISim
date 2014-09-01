/* 
*        invper3Dgrow grows the cluster and adds elements to the fringe if applicable.
*
*    James B. Silva <jbsilva@bu.edu>
*/

int getLocationToFromBond3D(int bond, int locX, int locY, int locZ,int L);
int getXfromBondInd(int ind, int L, int nbonds);
int getYfromBondInd(int ind, int L, int nbonds);
int getZfromBondInd(int ind, int L, int nbonds);

__kernel void invperc3Dgrow( __global int* system, __global int* fringe,
                 __global int* paramInt, int bond, int initLocation,
                 int numThreads, __global float* errorBuff) { 
    // Get the index of the current element to be processed
    int i = get_global_id(0);
    
    // bound check, equivalent to the limit on a 'for' loop
    if ( i >= numThreads){
        return;
    }
    
    int nbonds = 6;
    int L = paramInt[0];

    // remove from fringe
    if(i == 0){fringe[bond] = 0;}
    
    int x = getXfromBondInd(bond,L,nbonds);
    int y = getYfromBondInd(bond,L,nbonds);
    int z = getZfromBondInd(bond,L,nbonds);
    
    int oldLoc = x+y*L+z*(L*L);
    int newLoc = getLocationToFromBond3D(bond, x, y, z, L);
    
    
    // add to fringe
    int newBond = nbonds*newLoc+i;
    x = getXfromBondInd(newBond,L,nbonds);
    y = getYfromBondInd(newBond,L,nbonds);
    z = getZfromBondInd(newBond,L,nbonds);
    int temp;
    
    // add to fringe if good
    temp = getLocationToFromBond3D(newBond, x, y, z , L);
    if( (temp != oldLoc) && (system[temp] == 0) ){
        fringe[ newBond ] = 1; 
    }
    
    if(newLoc == initLocation){
        fringe[ newBond ] = 1; 
    }

    if(i == 0){
        system[newLoc] = 1;
        errorBuff[0] = newLoc;
        errorBuff[1] = oldLoc;
        errorBuff[2] = bond;
        errorBuff[3] = newBond;
        errorBuff[4] = fringe[bond];
        errorBuff[5] = x;
        errorBuff[6] = y;
        errorBuff[7] = z;
        errorBuff[8] = L;
        errorBuff[9] = system[newLoc];
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
    if(bondDir == 0){
            toLoc = locX + ((locY+1)%L)*L + locZ*(L*L);
    }else if(bondDir == 1){
            toLoc = ((locX+1)%L) + locY*L + locZ*(L*L);
    }else if(bondDir == 2){
            toLoc = locX + ((locY+L-1)%L)*L + locZ*(L*L);
    }else if(bondDir == 3){
            toLoc = ((locX+L-1)%L) + locY*L + locZ*(L*L);
    }else if(bondDir == 4){
            toLoc = locX + locY*L + ((locZ+1)%L)*(L*L);
    }else if(bondDir == 5){
            toLoc = locX + locY*L + ((locZ+L-1)%L)*(L*L);
    }
    return toLoc;
}



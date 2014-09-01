/*
        Metropolis algorithm implementation for nearest neighbor Ising model with dilution and
    static spins for a fully connected lattice. This implementation simulates
    multiple systems.

    James B. Silva <jbsilva@bu.edu>
*/
__kernel void ising_fullymetro(__global int* Ms,__global float* Es,__global int* NupS,__global int* NdownS,__global const int* NupFixS,__global const int* NdownFixS, __global const float* hFields,__global const float* temps,__global const float* jInteractions, __global float* random, const int numSys,const int N, __global float* errorBuff) {
        // Get the index of the current element to be processed
	int sysId = get_global_id(0);

	// bound check, equivalent to the limit on a 'for' loop
        if (sysId >= numSys || sysId < 0)  {
            return;
        }

        // Get Params
	float jInteraction = jInteractions[sysId];
	float hField = hFields[sysId];
        float temp = temps[sysId];

        // Get System State values
        int NupFix = NupFixS[sysId];
        int NdownFix = NdownFixS[sysId];
        int M = Ms[sysId];
        int Nup = NupS[sysId];
        int Ndown = NdownS[sysId];
        int Ndynamic = Ndown+Nup;
        int sysN = numSys;

        // Make variables for metropolis
        float Pup; float acceptance;
        bool flipDown;
        int newSpin; int newM;
        float eOld; float eNew; float dE;
        float tossUpDown; float tossAccept;
        float energy = Es[sysId];

        for (int i = 0;i < N;i++){
            //tossUpDown = random[i*2+sysId*N*2+0];
            //tossAccept = random[i*2+sysId*N*2+1];
            tossUpDown = random[sysId+2*sysN*i];
            tossAccept = random[sysId+(2*i+1)*sysN];
            //tossUpDown = random[(int)(tossUpDown*(2*N*(1+numSys)))];
            //tossAccept = random[(int)(tossAccept*(2*N*(1+numSys)))];

            // current state of system
            Pup = ((float)Nup)/((float)Ndynamic);
            // figure out flip direction
            flipDown = (tossUpDown < Pup) ? 1 : 0;
            //if(tossUpDown < Pup){flipDown = 1;}else{flipDown=0;}
            newSpin = (flipDown) ? (-1) : 1;
            //if(flipDown){newSpin = -1;}else{newSpin=1;}
            
            // Metropolis
            newM = M+2*newSpin;
            eOld = -((0.5*jInteraction*M*M)+(hField*M));
            eNew = -((0.5*jInteraction*newM*newM)+(hField*newM));
        
            dE = eNew-eOld;
            acceptance =  exp(-1.0*dE/temp);
            if( (dE<=0) || (acceptance > tossAccept) ){
                    M = newM;
                    energy = eNew;
                    if(newSpin == 1){Nup++;Ndown--;}else{Nup--;Ndown++;}
            }
        }
     
        Ms[sysId] =  M;
        NupS[sysId] =  Nup;
        NdownS[sysId] =  Ndown;
        Es[sysId] = energy;
}


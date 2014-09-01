/*
        This is an interface for MWC64 random number generator.

    James B. Silva <jbsilva@bu.edu>
*/
#include "mwc64x.cl"
__kernel void mwc64x_interface(__global float *rand,unsigned long baseOffset){
    int tid = get_global_id(0);
    baseOffset = baseOffset+tid;
    mwc64xvec4_state_t rng;
    MWC64XVEC4_SeedStreams(&rng, baseOffset, 2);
    ulong4 x=convert_ulong4(MWC64XVEC4_NextUint4(&rng));
    rand[4*tid]   = ((float)x.x)/4294967295.0;
    rand[4*tid+1] = ((float)x.y)/4294967295.0;
    rand[4*tid+2] = ((float)x.z)/4294967295.0;
    rand[4*tid+3] = ((float)x.w)/4294967295.0;
}

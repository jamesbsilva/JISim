
/*
Copyright 2010-2011, D. E. Shaw Research.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of D. E. Shaw Research nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
/*
 * This file is the OpenCL kernel.  It gets preprocessed, munged
 * into a C string declaration and included in pi_opencl, so that
 * running the compiled pi_opencl does not depend on any include
 * files, paths etc.
 */

#include <Random123/threefry.h>

__kernel void random123_interface(__global int* rand, unsigned long seed ) {
    unsigned tid = get_global_id(0);
    seed = seed+tid;
    threefry4x32_key_t k = {{tid, 0xdecafbad, 0xfacebead, 0x12345678}};
    threefry4x32_ctr_t c = {{seed, 0xf00dcafe, 0xdeadbeef, 0xbeeff00d}};
    union {
        threefry4x32_ctr_t c;
        int4 i;
    } u;
    c.v[0]++;
    u.c = threefry4x32(c, k);
    long x1 = u.i.x, x3 = u.i.y;
    long x2 = u.i.z, x4 = u.i.w;

    rand[tid] = x1;
    rand[4*tid+1] = x2;
    rand[4*tid+2] = x3;
    rand[4*tid+3] = x4;
}


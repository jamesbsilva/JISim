__kernel void part2( __global float4* color, float dt, const int size)
{
    unsigned int x = get_global_id(0);
    unsigned int y = get_global_id(1);

    //you can manipulate the color based on properties of the system
    //here we adjust the alpha
    color[x+y*size].w = 1.0f;
    if(color[x+y*size].x < 0.99f){
	color[x+y*size].x = color[x+y*size].x + 0.00075f;
	color[x+y*size].y = color[x+y*size].y + 0.000525f;
	color[x+y*size].w = color[x+y*size].w + 0.000525f;
	color[x+y*size].z = color[x+y*size].z + 0.000625f;
    }else{
	color[x+y*size].x = 0.0f;
    }    
}

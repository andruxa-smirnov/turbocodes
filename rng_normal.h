#ifndef RNG_NORMAL
#define RNG_NORMAL

#include <stdlib.h>
#include <math.h>
#include "mt/mt19937ar.h"
/*
Thanks to my father for this hint.
This method can be found on Knuth Vol. 2 2ND ED., pp. 117 for a N(0,1), pp.127 for a N(u,sigma^2)

*/
float uniform_01(){
	//return((float)rand()/(float)RAND_MAX);
	//Improved with marsenne twister
	return((float)genrand_int32()/(float)ULONG_MAX);
}

float normal_rng (float sigma)
{
	float V1,V2, S;

	do{
		V1 = 2*uniform_01() -1;
		V2 = 2*uniform_01() -1;
		S = V1*V1 + V2*V2;

	}while(S >= 1.0f);
	return sigma * V1 * sqrt(-2.0f * log(S)/S);
}
#endif
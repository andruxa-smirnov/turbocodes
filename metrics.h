#ifndef METRICS
#define METRICS
#include "defs.h"

float euclideanDistance(float * a, float * b, uint size){
	float d=0;
	uint i;
	for(i=0;i<size;i++)
		d+=(a[i]-b[i])*(a[i]-b[i]); //(a_i-b_i) ^2
	return d;
}

float euclideanDistanceCoded(float * a, bit * b, uint size){
	float d=0;
	float temp;
	uint i;
	for(i=0;i<size;i++){
		temp=(float)(2*b[i]-1);
		d+=(a[i]-temp)*(a[i]-temp); //(a_i-b_i) ^2
	}
	return d;
}
float hammingDistance(bit * a, bit * b, uint size){
	uint d=0;
	uint i;
	for(i=0;i<size;i++)
		d+= a[i] ^ b[i]; //bitwise xor
	return ((float)d);
}

#endif
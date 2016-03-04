#include "defs.h"

typedef struct{
	bit * rbit; //register bits
	int clen; //constraint length (number of registers)
}shiftR;

void clearRegister(shiftR* sr){
	memset(sr->rbit,0,sizeof(bit)*sr->clen); //braso a 0 i registri
}

shiftR* allocRegister(int constraint_len){
	shiftR* sr;
	sr = (shiftR*) malloc(sizeof(shiftR));
	sr->rbit = (bit*)malloc(sizeof(bit)*constraint_len);
	memset(sr->rbit,0,sizeof(bit)*constraint_len); //braso a 0 i registri
	sr->clen = constraint_len;
	return sr;
}
bit insertRegister(shiftR * sr, bit newbit){
		int i;
		bit oldbit;
		oldbit=sr->rbit[sr->clen-1];
		for(i=sr->clen-1;i>=0;i--)
			sr->rbit[i]=sr->rbit[i-1];
		sr->rbit[0]=newbit;
		return oldbit;
}
void freeRegister(shiftR * sr){
	free(sr->rbit);
	sr->rbit=NULL;
	free(sr);
}
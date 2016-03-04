/*
Code by Stefano Guerrini. December 2015.
This a part of my work for Coding Theory Course.
University of Ferrara - Master degree in electronics & communication Engineering.

------
This library implements a SINGLE convolutional encoder with or without feedback.

*/
#ifndef CONV_ENCODER
#define CONV_ENCODER


#include "shiftRegister.h"
#include <stdio.h>

typedef struct{
	shiftR* memory;
	gen * generators; //feedforward generator
	gen feedbackGen; //feedback generator
	int n;//n codice
	int L;//vincolo
}convEncoder;

void setEncoderState(convEncoder* ce, gen stateGen){
	int k;
	for(k=0;k<ce->L;k++)
		ce->memory->rbit[ce->L-1-k]=(stateGen >> k) & 1;
}
gen getEncoderState(convEncoder* ce){
	gen g=0;
	int k;
	//for(k=0;k<ce->L-1;k++) //We must not take the first register value!!
	for(k=0;k<ce->L;k++)
		g += ce->memory->rbit[ce->L-1-k] << k;
	return g;
}


/**
* The decoders goes through 2 steps:
* 1) Register shift and (if present) feedback loop [stepEncoder()]
* 2) Encode output by using generator's rule [encodeNowToOut()]
* The steps are done separately and you should call encodeAndGoOut() to get the full steps.
*/

//This not provide buffer control, internal use only!
void encodeNowToOut(convEncoder* ce, bit * out){
	int i=0,k=0;
	gen curgen;
	gen bitgen;
	int cur;
	//out = (bit*)malloc(sizeof(bit) * ce->n);
	for(i=0;i< ce->n ;i++){
		curgen = ce->generators[i];
		cur = 0; //calculate current output bit
		for(k=0;k<ce->memory->clen;k++){
			bitgen = (curgen >> k) & 1; //bit equivalent of the generator
			cur+=ce->memory->rbit[ce->memory->clen-1-k] & bitgen;
		}
		cur = cur%2; //X-OR (Sum in binary domain)
		out[i]=cur;
	}
}
/*
Brief description:
-> gets a bit in input
-> by memory he knows the older history
-> Returns a n-lenght bit string
*/
bit* encodeNow(convEncoder* ce, bit input){
	//I registri sono del tipo rbit[0] = u(i) , rbit[1]=u(i-1), ...
	//I generatori sono del tipo gen[0]= g(i) ,...
	bit * out;
	out = (bit*)malloc(sizeof(bit) * ce->n);
	encodeNowToOut(ce,out);
	return out;
}

//useful for trellis termination in turbo coder
bit stepEncoderFeedbackToForceTerminate(convEncoder * ce){

	bit realIn=0; //realIn=input if no feedback loop exists
	int bitgen;
	int k;
	bit fbk;

	

	if(ce->feedbackGen != 0){ //if feedback != 0, we need to overwrite the first register with the feedback loop one.
		for(k=0;k<ce->memory->clen-1;k++){ //first register never enters loop!
				bitgen = (ce->feedbackGen >> k) & 1; //bit equivalent of the generator
				realIn+=ce->memory->rbit[ce->memory->clen-1-k] & bitgen;
		}
		fbk = realIn % 2 ;
		
		realIn = (realIn + fbk)%2;
		ce->memory->rbit[0]=realIn; //overwrite first register
	}
	insertRegister(ce->memory,0); //shifto i registri!
	return fbk;
}

void stepEncoder(convEncoder * ce, bit input){
	bit realIn=0; //realIn=input if no feedback loop exists
	int bitgen;
	int k;

	insertRegister(ce->memory,input); //shifto i registri!

	if(ce->feedbackGen != 0){ //if feedback != 0, we need to overwrite the first register with the feedback loop one.
		for(k=0;k<ce->memory->clen-1;k++){ //first register never enters loop!
				bitgen = (ce->feedbackGen >> k) & 1; //bit equivalent of the generator
				realIn+=ce->memory->rbit[ce->memory->clen-1-k] & bitgen;
		}
		realIn = ((realIn%2) + input)%2;
		ce->memory->rbit[0]=realIn; //overwrite first register
	}
	
}
//Again, this does not provide buffer overflow control! Internal use only
void encodeAndGoOut(convEncoder* ce, bit input,bit * out){
	stepEncoder(ce,input);
	encodeNowToOut(ce,out);
}
bit* encodeAndGo(convEncoder* ce, bit input){
	bit * out;
	stepEncoder(ce,input);
	out = encodeNow(ce,input);
	return out;
}
void setGenerator(convEncoder* ce, int no, gen gen){
	ce->generators[no]=gen;
}
void setFeedbackGenerator(convEncoder * ce, gen bgen){
	ce->feedbackGen = bgen;
}
//Codice (n,k) L=user def. Unico vincolo: Codici binari: k=1 per def
convEncoder* createEncoder(int code_len, int constraint_len){
	convEncoder * ce;
	ce = (convEncoder*) malloc(sizeof(convEncoder)); //create obj
	ce->generators = (gen*)malloc(sizeof(int)*code_len); //alloc generators, we need n of them
	ce->feedbackGen = 0; //no feedback!
	ce->memory = allocRegister(constraint_len);
	ce->n = code_len;
	ce->L= constraint_len;

	return ce; //All done!
}
void destroyEncoder(convEncoder* ce){
	free(ce->generators); //Destroys generators
	freeRegister(ce->memory); //frees registers

	free(ce);
}

#endif
/*
* This code actually implements a turbo encoder with a nominal code ratio of 1/3 (n=2 in each recursive systematic encoder)
*/

#include "convEncoder.h"
#include <stdlib.h>
#include "graph_turbo.h"

typedef struct{
	convEncoder * c1;
	convEncoder * c2;
	trellis * te;
}turboEncoder;


/*
*
* CCSDS 1/3 Turbo Code implementation.
* See reference: TM SYNCHRONIZATION AND CHANNEL CODING : CCSDS 131.0-B-2 Page 39
* We have two identical RSC encoders in parallel (PCCC)
* Encoder Details:
* Ratio = 1/2 (k=1, n=2) ; Constraint Length = 5 ; No puncturing ; 
* Forward connection: 11011 = 0x1B 
* Backward connection: 10011 = 0x13
*/
#define CCSDS_13_FWD 0x1B
#define CCSDS_13_BWD 0x13
turboEncoder get_CCSDS_13_turboEncoder(){
	turboEncoder te;
	
	te.c1 = createEncoder(1,5);
	te.c2 = createEncoder(1,5);

	setGenerator(te.c1,0,CCSDS_13_FWD);
	setFeedbackGenerator(te.c1,CCSDS_13_BWD);
	setGenerator(te.c2,0,CCSDS_13_FWD);
	setFeedbackGenerator(te.c2,CCSDS_13_BWD);

	te.te = createTrellisFromEncoder(te.c1);
	return te;
}
void freeTurboEncoder (turboEncoder * te){
	destroyEncoder(te->c1);
	destroyEncoder(te->c2);
	freeTrellis(te->te);
	/*free(te->c1);
	free(te->c2);
	free(te->te);*/
}
//p1 = 31; p2 = 37; p3 = 43; p4 = 47; p5 = 53; p6 = 59; p7 = 61; p8 = 67
const uint turbopGen[8]={31,37,46,47,53,59,61,67};

//Blocklen should be chosen from 1784, 3568, 7136, 8920
void CCSDS_turboScrambler_createLookupFile(uint blockLen){
	long m;
	long s;
	long k1=8;
	long k2=blockLen/8;
	long k=blockLen;
	long i,j,t,q,c;
	long ps;

	FILE * fp;
	char outstr[256];
	//char outname[256];
	sprintf(outstr,"CCSDS_Scrambler_Lookup_%d.h",blockLen);
	fp = fopen(outstr,"w+");
	sprintf(outstr,"const unsigned short CCSDS_Scrambler_Lookup_%d [%d] = {\n",blockLen,blockLen);
	fputs(outstr,fp);

	for(s=1;s<=k;s++){
		m = (s-1)%2;
		i = (s-1)/(2*k2);
		j = ((s-1)/2) - i*k2;
		t = (19*i +1) % (k1/2);
		q = (t % 8) +1;
		c = (turbopGen[q-1]*j + 21*m) % k2;
		ps = 2*(t+ c*k1/2 + 1) - m;
		ps--; //perchè nel documento CCSDS è tutto riferito rispetto all'indice partendo da 1, cane cane lech
		sprintf(outstr,"%d",ps);
		
		fputs(outstr,fp);
		if(s<k)
			fputs(",",fp);
		//out[s] = in[ps];
	}
	fputs("};",fp);
	fclose(fp);
}

void fdescramble(float * in, float * out, uint size, const unsigned short * slookup){
	uint i=0;
	for(i=0;i<size;i++)
		out[slookup[i]] = in[i];
}
void fscramble(float * in, float * out, uint size, const unsigned short * slookup){
	uint i=0;
	for(i=0;i<size;i++)
		out[i] = in[slookup[i]];
}


//Valid only for a Rc=1/3 code
//Messi nel più comodo formato
#define trailing_zeros 4
#define trellis_term 4
float stability_factor=0.15;

bit * CCSDS_turboDecodeBlock(turboEncoder * te, float * input, uint blockLen, const unsigned short * slookup, const unsigned short iterations, float noiseVar, int stability, bit * original){

	uint ebl,eil; //effective block length = blockLen + 4 trailing zeros
	float * extrinsic=NULL;
	float * extrinsic2=NULL;
	float * sextr, * dsextr=NULL;
	float * scrambledin;
	float * L;
	float temp;
	//float * input2;
	bit * out;
	uint errcnt;
	uint i;
	uint step;

	float mean,var;
	//FILE * fp = fopen("step_change.csv","w");

	ebl= blockLen + trellis_term;
	eil= blockLen + 2*trellis_term;
	out = (bit*) malloc(sizeof(float)*blockLen);

	L= (float*)malloc(sizeof(float)*ebl);
	scrambledin = (float*) malloc(sizeof(float)*ebl);
	sextr = (float*)malloc(sizeof(float)*ebl);
	dsextr = (float*)malloc(sizeof(float)*ebl);

	fscramble(input,scrambledin,blockLen,slookup); //non devo scramblare i bit di terminazione del trellis!
	memcpy(&scrambledin[blockLen],&input[blockLen+trellis_term],trellis_term*sizeof(float)); //copio i bit di terminazione del trellis (da non scramblare)

	memset(dsextr,0,sizeof(float)*ebl); //set to 0 (initializing)
	memset(sextr,0,sizeof(float)*ebl); //set to 0 (initializing)

	

	
	for(step=0;step<iterations;step++){
		if(extrinsic){
			free(extrinsic);
			free(extrinsic2);
		}
		extrinsic = BCJR_AWGN_TURBO(te->te,input,&input[eil],ebl,noiseVar,dsextr); //compute extrinsic for 1st encoder
		
		fscramble(extrinsic,sextr,blockLen,slookup); //extr ==> (SCRAMBLE) ===> sextr, non devo scramblare i llr di terminazione

		extrinsic2 = BCJR_AWGN_TURBO(te->te,scrambledin,&input[ebl+eil],ebl,noiseVar,sextr); //compute extrinsic
		
		fdescramble(extrinsic2,dsextr,blockLen,slookup); //extr2 ==> (DESCRAMBLE) ===> dsextr, non devo scramblare i llr di terminazione
		
		
	}
	
	//extrinsic = (float*)malloc(sizeof(float)*ebl);
	for(i=0;i<ebl;i++)
		extrinsic[i] = scrambledin[i]*2/noiseVar + sextr[i] + extrinsic2[i] ;

	fdescramble(extrinsic,L,blockLen,slookup); //again, non devo scramblare i llr di terminazione
	
	for(i=0;i<blockLen;i++)
		out[i] = (L[i] >= 0); //i blocchi di terminanzione un me servono

	free(extrinsic);
	free(extrinsic2);
	free(L);
	free(sextr);
	free(dsextr);
	free(scrambledin);
	//fclose(fp);
	return out;
}


//Out sarà di questo tipo: ...,xu[k],yp[k],yq[k],...,
//Alla fine saranno presenti 4*2 *2 = 16 bit di terminazione trellis 
//[SIST][Term1_sist][Term2_sist][Y1_par][Term1_par][Y2_par][Term2_par]
bit * CCSDS_turboEncodeBlock(turboEncoder * te, bit * input, uint blockLen, const unsigned short * slookup){
	bit * out;
	bit curr;
	uint i;
	uint eil = blockLen + 2*trellis_term; //effective information length = blockLen*3
	uint erl = blockLen + trellis_term; //effective redundancy length
	uint ebl; //effective block length = blockLen *3  + trellis term
	//eil = blockLen *3 + 2*trellis_term;

	ebl= blockLen*3 + trellis_term*2*2;
	out = (bit*)malloc(sizeof(bit)*ebl); //Rate 1/3 + 4 padding zeros for trellis terminate

	for(i=0;i<blockLen;i++){
		out[i] = input[i];

		encodeAndGoOut(te->c1,input[i],&curr);
		out[i + eil] = curr;

		encodeAndGoOut(te->c2,input[slookup[i]],&curr);
		out[i + eil + erl] = curr;

	}
	//terminate first trellis
	for(i=blockLen;i<blockLen+trellis_term;i++){
		curr = stepEncoderFeedbackToForceTerminate(te->c1);
		out[i] = curr;
		encodeNowToOut(te->c1,&curr);
		out[i+eil] = curr;
	}

	for(i=blockLen;i<blockLen+trellis_term;i++){
		curr = stepEncoderFeedbackToForceTerminate(te->c2);
		out[i+trellis_term] = curr;
		encodeNowToOut(te->c2,&curr);
		out[i+eil+erl] = curr;
	}
	return out;
}
#ifndef CODER_GRAPHT
#define CODER_GRAPHT
#include "graph.h"
#include <stdlib.h>
#include <string.h>

//rcvd sarà il bit di parità 1 ricevuto (yp[k]), urcvd il corrispondente bit non cod, 
//encout sarà il bit di parità del ramo in esame (p[k]), curtran la transizione corrente 0/1
float gammamet(float urcvd, float rcvd , bit curtran, bit encout, float var, float extrinsic){
	
	return ((float)(2*curtran-1))*extrinsic/2 + ((float)(2*curtran-1))*urcvd/var + ((float)(2*encout-1))*rcvd/var;
//	return extrinsic + ((float)(2*curtran-1))*urcvd/var + ((float)(2*encout-1))*rcvd/var;
}
void BCJR_forward_TURBO(trellis * te, float *rcvd, float * uncodedrcvd, uint steps, float var, float * extrinsic, float ** out_alpha , float ** out_gamma){
	uint i=0;
	uint j=0,k=0;
	float * precalpha;
	float * alpha;
	float * temp;
	float * gamma;
	float * alphaArr;
	float curgamma;

	alpha = (float*)malloc(sizeof(float)*te->no);
	precalpha=(float*)malloc(sizeof(float)*te->no);
	temp= (float*)malloc(sizeof(float)*te->no);
	gamma = (float*)malloc(sizeof(float)*steps*2*te->no); //in generale gamma[passo][stato][0/1bit or equiv: nextstate]
	alphaArr = (float*)malloc(sizeof(float)*(steps+1)*te->no); //in generale gamma[passo][stato][0/1bit or equiv: nextstate]

	memset(alpha,FLOAT_MINF,sizeof(float)*te->no); //fast -INF float assignment
	memset(precalpha,FLOAT_MINF,sizeof(float)*te->no);
	memset(temp,FLOAT_MINF,sizeof(float)*te->no); //-INF assign leads to trascurable terms in max*
	precalpha[0] = 0;
	
	memcpy(alphaArr,precalpha,sizeof(float)*te->no);

	for(i=0;i<steps;i++){ //do subsequent steps
		for(k=0;k<te->no;k++){ //calcolo su tutti i sigma_i
			for(j=0;j<te->no;j++){ //sommo su tutti i sigma_i-1, devo farlo perchè non conosco i sigma_i-1 dato sigma i (conosco solo i rami diretti del grafo!)
				if(te->nodes[j][0] == k){
					curgamma = gammamet(uncodedrcvd[i*te->coden],rcvd[i*te->coden],0,te->seq[(j*2+0)*te->coden],var,extrinsic[i]);
					gamma[((i*te->no+j)*2)+0]=curgamma;
					temp[j] = precalpha[j] + curgamma;
				}
				else if(te->nodes[j][1]==k){

					curgamma = gammamet(uncodedrcvd[i*te->coden],rcvd[i*te->coden],1,te->seq[(j*2+1)*te->coden],var,extrinsic[i]);
					gamma[((i*te->no+j)*2)+1]=curgamma;
					temp[j] = precalpha[j] + curgamma;
				}
			}
			alpha[k]=maxstar(temp,te->no);
			memset(temp,FLOAT_MINF,sizeof(float)*te->no);
		}
		memcpy(&alphaArr[(i+1)*te->no],alpha,sizeof(float)*te->no);
		flip(&precalpha,&alpha); //non ho bisogno di sistemare alpha perchè tanto poi lo ricalcolo e mi preservo il memcpy
		memset(alpha,FLOAT_MINF,sizeof(float)*te->no); //fast -INF float assignment
	}
	free(precalpha);
	free(temp);
	free(alpha);
	*out_gamma = gamma;
	*out_alpha = alphaArr;
//	freeTrellis(supp);
}
void BCJR_backward_TURBO(trellis * te, float * in_gamma, uint steps, float ** out_beta){
	float * beta;
	float * succbeta;
	float * temp;
	float * betaArr;
	int i,j,k;

	beta = (float*)malloc(sizeof(float)*te->no);
	succbeta=(float*)malloc(sizeof(float)*te->no);
	temp= (float*)malloc(sizeof(float)*te->no);
	memset(beta,FLOAT_MINF,sizeof(float)*te->no); //fast -INF float assignment
	memset(succbeta,FLOAT_MINF,sizeof(float)*te->no);
	
	betaArr = (float*)malloc(sizeof(float)*(steps+1)*te->no);

	succbeta[0] = 0; //supposes a 0 terminating trellis
	memcpy(&betaArr[(steps)*te->no],succbeta,sizeof(float)*te->no);
	for(i=steps-1;i>=0;i--){
		for(k=0;k<te->no;k++){ //calcolo su tutti i sigma_i
			for(j=0;j<2;j++){ //qui i rami del grafo li conosco!!!
					//temp[j]=succbeta[te->nodes[k][j]]+in_gamma[(((i-1)*te->no+k)*2)+j];
					temp[j]=succbeta[te->nodes[k][j]]+in_gamma[((i*te->no+k)*2)+j];
			}
			beta[k]=maxstar(temp,2);
		}
		memcpy(&betaArr[(i)*te->no],beta,sizeof(float)*te->no);
		flip(&succbeta,&beta); //non ho bisogno di sistemare alpha perchè tanto poi lo ricalcolo e mi preservo il memcpy
		memset(beta,FLOAT_MINF,sizeof(float)*te->no); 
	}	
	free(succbeta);
	free(beta);
	free(temp);

	*out_beta = betaArr;
}

float bitToFloat(bit b){
	return ((float)(2*b-1));
}

void calcParam(float * in, uint size, float * mean, float* var){
	uint i=0;
	*mean=0;
	*var=0;
	for(i=0;i<size;i++)
		(*mean)+= in[i];
	(*mean)/=(float)size;
	for(i=0;i<size;i++)
		(*var)+= (in[i] - (*mean))*(in[i] - (*mean));
	(*var)/=(float)size;
}
//IT...COULD...WORK! Returns extrinsic only!
float* BCJR_AWGN_TURBO(trellis * te, float * uncodedseq, float * encseq, uint sequenceSize, float noiseVar, float * extrinsic){
	float * gamma;
	float * alpha;
	float * beta;
	float * max1, * max2;
	float * L;
	uint i=1, j=0,k=0,l=0;
	uint steps = sequenceSize/te->coden;
	BCJR_forward_TURBO(te,encseq,uncodedseq,steps,noiseVar,extrinsic,&alpha,&gamma); //calculate gamma & alpha
	BCJR_backward_TURBO(te,gamma,steps, &beta); //calculate beta
	
	max1 = (float *) malloc(sizeof(float)*te->no);
	max2 = (float *) malloc(sizeof(float)*te->no);
	L = (float *)malloc(sizeof(float)*steps);
	//now let's calculate L(ui)
	//ogni u_i
	for(i=1;i<=steps;i++) //passo
	{
		for(k=0;k<te->no;k++){ //tutti i sigma_i-1
			//max1[k]=alpha[(i-1)*te->no +k] + beta[(i)*te->no + te->nodes[k][1]] + encseq[i-1]*(2*te->seq[k*2+1]-1)/noiseVar;//+ gamma[(((i-1)*te->no+k)*2)+1]; //S1
			//max2[k]=alpha[(i-1)*te->no +k] + beta[(i)*te->no + te->nodes[k][0]] +  encseq[i-1]*(2*te->seq[k*2+0]-1)/noiseVar; //S0
			max1[k]=alpha[(i-1)*te->no +k] + beta[(i)*te->no + te->nodes[k][1]] + encseq[i-1]*(bitToFloat(te->seq[k*2+1]))/noiseVar;//+ gamma[(((i-1)*te->no+k)*2)+1]; //S1
			max2[k]=alpha[(i-1)*te->no +k] + beta[(i)*te->no + te->nodes[k][0]] +  encseq[i-1]*(bitToFloat(te->seq[k*2+0]))/noiseVar; //S0
		}
		L[i-1] = maxstar(max1,te->no) - maxstar(max2,te->no);	
	}
	free(max1);
	free(max2);
	free(gamma);
	free(alpha);
	free(beta);
	return L;/**/
}
#endif

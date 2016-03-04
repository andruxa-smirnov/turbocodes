#ifndef CODER_GRAPH
#define CODER_GRAPH
#include "convEncoder.h"
#include "metrics.h"
#include <math.h>
#include <float.h>
#include "max_lookup.h"
//float x[] = {0.000046f};
/*typedef struct{
	int no;
	convEncoder * se; //support Encoder
	float * metrics;
	unsigned int bestNode;
}strellis;
*/
//La posizione di un nodo coincide con il suo stato ovvero nodes[0] => stato 0, nodes[1] => stato 1, nodes [2] => stato 2 = 10
typedef struct{
	uint (*nodes)[2]; //a quale stato vado con ingresso 0/1
	bit *seq; //con che uscita vado??
	int no; //state number
	int coden; //code n value
}trellis; 

trellis * createTrellis(unsigned int statusNo){
	trellis * t;
	t = (trellis*)malloc(sizeof(trellis));
	t->no = statusNo;
	t->nodes = (uint(*)[2])malloc(sizeof(uint)*2*t->no);
	//t->seq = (bit(*)[2])malloc(sizeof(bit)*t->no);
	t->seq=NULL;
	return t;
}
void freeTrellis(trellis * te){
	free(te->nodes);
	if(te->seq)
		free(te->seq);
	free(te);
}
trellis* copyTrellis(trellis * te){
	trellis * nt;
	nt = createTrellis(te->no);
	memcpy(nt->nodes,te->nodes,sizeof(uint)*2*te->no);
	return nt;
}
trellis * createTrellisFromEncoder(convEncoder * ce){
	trellis * t;
	uint cs; //current status
	uint statusNo = (uint)(pow(2.f,ce->L-1)); //numero di stati
	uint future;
	bit * enc;
	int i;
	enc = (bit*)malloc(sizeof(bit)*ce->n);
	t=createTrellis(statusNo);
	t->coden = ce->n;
	t->seq = (bit*)malloc(sizeof(bit)*statusNo*t->coden*2); //per accederci devo fare: t->seq[N.STATO*t->coden + 0/1*t->coden]

	for(cs=0;cs<statusNo;cs++){
		for(i=0;i<2;i++){
			setEncoderState(ce,cs << 1);
			//enc = encodeAndGo(ce,i);
			encodeAndGoOut(ce,i,enc); //Mi risparmio di fare una sbanderva di malloc!
			future = (getEncoderState(ce) >> 1);
			t->nodes[cs][i] = future;
			memcpy(&(t->seq[(cs*2+i)*t->coden]),enc,t->coden);

			//t->seq[cs*t->coden][i] = enc;			
		}
	}
	free(enc);
	setEncoderState(ce,0); //reset status
	return t;
}

int isFloatMin(float x){
	return (*(unsigned int*)(& (x)))==FLOAT_MINF4;
}
/**
MAX* IMPLEMENTATION:
by definition:
1) max*{x,y} = ln(e^x + e^y);
By equation 8.8-31 from Proakis - Digital Communications 5ed. (p.548):
2) max*{x,y} = max{x,y} + ln(1+ exp(-|x-y|))

The second is far more accurate in the float domain because, if for example:
x=-800
y=-850
1) => e^x = 0.f, e^y = 0.f =>  max* = -1.#INF
2) => max{x,y} = x = -800, ln(...) ~= ln(1) = 00 => max* = max = -800 != -INF!!!
**/
float lookup_maxstar(float val){
	if(val >=10.f)
		return 0;
	return max_lookup[(uint)(val*1000)];
	/*uint prova;
	float ret;
	if(val >=10.f)
		return 0;
	prova = val*1000;
	ret = max_lookup[prova];
	return 0.f;*/
}
__inline float maxstar_two(float x, float y){
	float out;
	out = x > y ? x : y ;
	//return out;
//	return out + log(1+exp(-abs(x-y))); //exact max*
#ifdef LOOKUP_MAXSTAR
	return out + lookup_maxstar(fabs(x-y)); //max* using lookup table
#else
	return out; //max* = max (approx.)
#endif
	

	//return out + log(1+exp(-abs(x-y))); //exact max*
	//return log(exp(x)+exp(y)); 
}
//Avoid recursion, not much elegant but much more performant!
float maxstar(float * args, int argc){
	uint i;
	float tmp;
	tmp = maxstar_two(args[0],args[1]);

	if(argc == 2) //perchè mi devi chiamare se arg = 2? MALEDETTO!
		return tmp;
	//max(a,b,c)=max(max(a,b),c)
	for(i=2;i<argc;i++)
		tmp = maxstar_two(tmp,args[i]);
	return tmp;
}
void flip(float ** ptr1, float **ptr2){
	float * tmp;
	tmp = *ptr1;
	*ptr1=*ptr2;
	*ptr2 = tmp;
}

void BCJR_forward(trellis * te, float *rcvd, uint steps, float ** out_gamma, float ** out_alpha, float var){
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
	memset(alpha,FLOAT_MINF,sizeof(float)*te->no); //fast -INF float assignment
	memset(precalpha,FLOAT_MINF,sizeof(float)*te->no);
	memset(temp,FLOAT_MINF,sizeof(float)*te->no); //-INF assign leads to trascurable terms in max*
	precalpha[0] = 0;

	gamma = (float*)malloc(sizeof(float)*steps*2*te->no); //in generale gamma[passo][stato][0/1bit or equiv: nextstate]
	alphaArr = (float*)malloc(sizeof(float)*(steps+1)*te->no); //in generale gamma[passo][stato][0/1bit or equiv: nextstate]
	memcpy(alphaArr,precalpha,sizeof(float)*te->no);

	for(i=0;i<steps;i++){ //do subsequent steps
		for(k=0;k<te->no;k++){ //calcolo su tutti i sigma_i
			for(j=0;j<te->no;j++){ //sommo su tutti i sigma_i-1, devo farlo perchè non conosco i sigma_i-1 dato sigma i (conosco solo i rami diretti del grafo!)
				if(te->nodes[j][0] == k){
					curgamma=-euclideanDistanceCoded(&rcvd[i*te->coden],&te->seq[(j*2+0)*te->coden],te->coden);
					curgamma/=2*var;
					gamma[((i*te->no+j)*2)+0]=curgamma;
					temp[j] = precalpha[j] + curgamma;
				}
				else if(te->nodes[j][1]==k){
					curgamma=-euclideanDistanceCoded(&rcvd[i*te->coden],&te->seq[(j*2+1)*te->coden],te->coden);
					curgamma/=2*var;
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
void BCJR_backward(trellis * te, float * in_gamma, uint steps, float ** out_beta){
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

//IT...COULD...WORK!
float* BCJR_AWGN(trellis * te, float * seq, uint sequenceSize, float noiseVar){
	float * gamma;
	float * alpha;
	float * beta;
	float * max1, * max2;
	float * L;
	uint i=1, j=0,k=0,l=0;
	uint steps = sequenceSize/te->coden;
	BCJR_forward(te,seq,steps,&gamma,&alpha, noiseVar); //calculate gamma & alpha
	BCJR_backward(te,gamma,steps, &beta); //calculate beta

	max1 = (float *) malloc(sizeof(float)*te->no);
	max2 = (float *) malloc(sizeof(float)*te->no);
	L = (float *)malloc(sizeof(float)*steps);
	//now let's calculate L(ui)
	//ogni u_i
	for(i=1;i<=steps;i++)
	{
		for(k=0;k<te->no;k++){ //tutti i sigma_i-1
			max1[k]=alpha[(i-1)*te->no +k] + beta[(i)*te->no + te->nodes[k][1]] + gamma[(((i-1)*te->no+k)*2)+1];
			max2[k]=alpha[(i-1)*te->no +k] + beta[(i)*te->no + te->nodes[k][0]] + gamma[(((i-1)*te->no+k)*2)+0];
		}
		L[i-1] = maxstar(max1,te->no) - maxstar(max2,te->no);	
	}
	free(max1);
	free(max2);
	free(gamma);
	free(alpha);
	free(beta);
	return L;
}
#endif

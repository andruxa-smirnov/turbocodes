//L'utilizzo o meno di questa define stabilisce l'utilizzo o meno della versione approssimata del MAP
#define LOOKUP_MAXSTAR

#include <Windows.h>
#include "include.h"
#include <time.h>
#include "turbo.h"
#include "scramblers\CCSDS_Scrambler_Lookup_1784.h"
#include "scramblers\CCSDS_Scrambler_Lookup_3568.h"
#include "scramblers\CCSDS_Scrambler_Lookup_7136.h"
#include "scramblers\CCSDS_Scrambler_Lookup_8920.h"


/*FILETIME win_time(){
	FILETIME ft;
	GetSystemTimeAsFileTime(&ft);
	return ft;
}
void difftime(FILETIME filetime, FILETIME filetime2){
	ULONGLONG time1,time2;
time1 = (((ULONGLONG) filetime.dwHighDateTime) << 32) + filetime.dwLowDateTime;
time2 = (((ULONGLONG) filetime2.dwHighDateTime) << 32) + filetime2.dwLowDateTime;
}*/
unsigned long win_time(){
	return timeGetTime();
}

//must be time2>time1
unsigned long diff_time(unsigned long time1, unsigned long time2){
	return time2-time1;
}

void test(){
	FILE * fp;
	int i;
	char a[100];
	fp = fopen("prova.csv","w");
	for(i=0;i<10000;i++){
		sprintf(a,"%f\n",normal_rng(1));
		fputs(a,fp);
	}
	fclose(fp);
	
}
bit * generateRandomBits(uint no, uint trailingzeros){
	uint i=0;
	bit * out;
	out = (bit*)malloc(sizeof(float)*no);
	
	for(i=0;i<no;i++)
		//out[i] = rand() % 2;
		out[i] = genrand_int32() % 2; //Marsenne
	for(i=1;i<=trailingzeros;i++)
		out[no - i] = 0; //terminate 0
	
	return out;
}
void corruptBit(float * bs, uint num, float noiseVar){
	uint i=0;
	for(i=0;i<num;i++)
		bs[i] += normal_rng(sqrt(noiseVar));
}

#define simSize 1784
#define frameLen 1784/8
#define DEBUGN (1784*3)+16
#define SCRAMBLER CCSDS_Scrambler_Lookup_1784
/**/
/*

#define simSize 8920
#define frameLen 8920/8
#define DEBUGN (8920*3)+16
#define SCRAMBLER CCSDS_Scrambler_Lookup_8920
/**//*
#define simSize 3568
#define frameLen 3568/8
#define DEBUGN (3568*3)+16
#define SCRAMBLER CCSDS_Scrambler_Lookup_3568
*//*
#define simSize 7136
#define frameLen 7136/8
#define DEBUGN (7136*3)+16
#define SCRAMBLER CCSDS_Scrambler_Lookup_7136
/****
Example: Sending an image
***/
void bytetobit(bit input, bit output[8]){
	uint i=0;
	for(i=0;i<8;i++)
		output[i] = (input >> i)%2;
}
void bittobyte(bit input[8], bit *out){
	uint i=0;

	*out=0;
	for(i=0;i<8;i++)
		*out = *out + (input[i] << i);
}

void createMaxLogAppLookup(){
	//x.xxx
	uint i=0;
	uint size=10000;
	uint len=size/10;
	float curr;
	float ll;
	FILE * fp=fopen("max_lookup.h","w");
	fprintf(fp,"const float max_lookup[]={",size);
	    
	for(i=0;i<size;i++){
		curr=(float)i/(float)len;
		ll=log(1 + exp(-curr));
		fprintf(fp,"%ff",ll);
		if(i<size-1)
			fprintf(fp,",");
	}
	fprintf(fp,"};",size);
	fclose(fp);
}
//Image Send Simulation
int main_image(){
		FILE * fp, *out;
		bit buf[1000000];
	//	bit bufout[1000000];
		uint fsize;
		uint cur;
		uint read;
		uint i;
		bit * turbobits;
		bit * decout;
		bit temp[simSize*8];
		bit byte[8];


		float snr, noisevar;
		float encodedbits[DEBUGN];

		turboEncoder te;

		//Set variables
		snr=0.6f; //Eb/N0 in db
		noisevar = 1.5 / pow(10,snr/10);

		te = get_CCSDS_13_turboEncoder();

		//Read image
		fp = fopen("lena.pbm","rb");
		fsize = fread(buf,sizeof(bit),1000000,fp);
		fclose(fp);
		
		cur=0;
		srand(time(NULL));
		out = fopen("lena_dec3_1784.pbm","wb");
		while(cur < fsize){ //devo splittare in blocchi
			read = fsize-cur;
			if(read > frameLen)	
				read = frameLen;
			else
				read = fsize-cur;
			printf("Simulating...%u of %u (%.2f %%)\n",cur,fsize,100*(float)cur/(float)fsize);
			memset(temp,0,sizeof(bit)*simSize);

			for(i=0;i<read;i++){
				bytetobit(buf[cur+i],byte);
				memcpy(temp+i*8,byte,sizeof(bit)*8);
			}
			

			//Codifico i bit
			turbobits = CCSDS_turboEncodeBlock(&te,temp,simSize,SCRAMBLER);
			
			//Codifico i bit in BPSK
			for(i=0;i<DEBUGN;i++)
				encodedbits[i] = ((float)(turbobits[i]*2))-1;

			//Gli metto il rumore
			corruptBit(encodedbits,DEBUGN,noisevar);

			//Decodifico!
			decout = CCSDS_turboDecodeBlock(&te,encodedbits,simSize,SCRAMBLER,10,noisevar,1,NULL);

			memset(temp,0,sizeof(bit)*simSize);
			for(i=0;i<read;i++)
				bittobyte(decout+i*8,&temp[i]);
			
			fwrite(temp,sizeof(bit),read,out);
			cur += read;

			free(turbobits);
			free(decout);
		}

		
		
		
		fclose(out);
		/*fp = fopen("lena_dec.bmp","wb");
		fwrite(buf,sizeof(char),read,fp);
		fclose(fp);*/
		
return 0;		
}

int main(){
	turboEncoder te;
	bit * inbits;
	bit * turbobits;
	bit * outbits;
	float noisevar;
	float snr;
	uint i,k;
	int it=0,itn;
	int testno = 500; //numero di prove da eseguire per ogni test!
	float err_rate;
	float encodedbits[DEBUGN];

	uint errcnt;
	unsigned long start,end,elapsed,totaltime; //start & end timers
	

	init_genrand(time(NULL));
	te = get_CCSDS_13_turboEncoder();

	printf("*****************************************\n");
	printf("\tCCSDS Turbo Code Simulator\n");
	printf("*****************************************\n");
	printf("***Simulation parameters***\n");
	printf("Simulation size (# of proves): ");
	scanf("%d",&testno);
	printf("SNR (dB): ");
	scanf("%f",&snr);
	printf("Number of decoder iterations: ");
	scanf("%d",&it);	
	printf("Press ENTER to start");
	getch();
			
	printf("\n\n");
	noisevar = 1.5f / pow(10,(snr)/10);

	errcnt=0;
	totaltime=0;
	for(i=0;i<testno;i++){
		inbits = generateRandomBits(simSize,0); //generate random bits

		
				
		turbobits = CCSDS_turboEncodeBlock(&te,inbits,simSize,SCRAMBLER); //let's encode this!
				
		for(k=0;k<DEBUGN;k++) //BPSK
			encodedbits[k] = ((float)(turbobits[k]*2))-1;

		corruptBit(encodedbits,DEBUGN,noisevar); //Send through AWGN
				
		start = win_time();
				
		outbits = CCSDS_turboDecodeBlock(&te,encodedbits,simSize,SCRAMBLER,it,noisevar,0,inbits);
				
		end=win_time();
		elapsed = end-start; //time needed for computation
		totaltime += elapsed;

		for(k=0;k<simSize;k++)
			if(outbits[k] != inbits[k])
				errcnt ++;
		//printf("\tCalls: %u\n",calls);
				
		//Free some memory...
		free(turbobits);
		free(inbits);
		free(outbits);
				
		printf("\#%d of %d) Total Bits: %u, Total Errors: %d, Error rate: %.2E\n",i,testno,(i+1)*simSize,errcnt,((float)errcnt)/((float)(i+1)*simSize));		
	}

	//Ok we now need to calculate some stats!
	err_rate = ((float)errcnt)/((float)testno*simSize);
			
	printf("\n\n*******************\nSimulation Results\n*******************\nIterations No: %d\nSNR: %f dB\nTotal generated bits: %d\nBER: %.2e\nTime needed: %u ms\n",it,snr,i*simSize,err_rate,totaltime/testno);
			



	getch();
	freeTurboEncoder(&te);
	return 0;
}

//Iterations study
int main_sim(){

	
	turboEncoder te;
	bit * inbits;
	bit * turbobits;
	bit * outbits;
	float noisevar;
	float snr;
	uint i,k;
	int it=0,itn;
	int testno = 500; //numero di prove da eseguire per ogni test!
	float err_rate;
	float encodedbits[DEBUGN];
	FILE * fp;
	uint errcnt;
	unsigned long start,end,elapsed,totaltime; //start & end timers
	unsigned int iterations[5] = {1, 2, 4, 10 ,18};
	//createMaxLogAppLookup();
	//exit(0);
	//Initializations
	srand(time(NULL));
	init_genrand(time(NULL));
	te = get_CCSDS_13_turboEncoder();

	fp=fopen("cancellami.csv","w");
	fprintf(fp,"iterations,snr,ber,time_ms\n"); 
	//Let's go!
//	for(stability_factor=0.05f;stability_factor<=1.f;stability_factor+=0.05f){
	//for(it=10;it<=10;it++){ //numero di iterazioni
	for(itn=3;itn<4;itn++){
		it = iterations[itn];
		
		
		snr = -0.1f;
		//if(itn==2)snr=0.8;
		while(snr<=0.4f){ //Set SNR
			/*if(snr <=0.3f)
			testno=20;
			else
				testno=500;
			/*
			if(snr >= 0.5f)
				testno=200;
			/*if(snr > 0.5f)
				testno=2000;
				*/
			noisevar = 1.5f / pow(10,(snr)/10);
		//	noisevar/=2;
			errcnt=0;
			totaltime=0;
			for(i=0;i<testno;i++){
				inbits = generateRandomBits(simSize,0); //generate random bits

				/*******
				****************************
				********UNCODED BPSK********
				****************************
				*******/
				/*for(k=0;k<simSize;k++) //BPSK
					encodedbits[k] = ((float)(inbits[k]*2))-1;
				corruptBit(encodedbits,simSize,noisevar);

				for(k=0;k<simSize;k++)
					if((encodedbits[k]>=0) != inbits[k])
						errcnt ++;*/
				printf("\#%d) #IT: %u, SNR:%.2f, Bit: %u, Errors: %d, Error rate: %.2E\n",i,it,snr,i*simSize,errcnt,((float)errcnt)/((float)i*simSize));
				
				turbobits = CCSDS_turboEncodeBlock(&te,inbits,simSize,SCRAMBLER); //let's encode this!
				
				for(k=0;k<DEBUGN;k++) //BPSK
					encodedbits[k] = ((float)(turbobits[k]*2))-1;

				corruptBit(encodedbits,DEBUGN,noisevar); //Send through AWGN
				
				start = win_time();
				
				outbits = CCSDS_turboDecodeBlock(&te,encodedbits,simSize,SCRAMBLER,it,noisevar,0,inbits);
				
				end=win_time();
				elapsed = end-start; //time needed for computation
				totaltime += elapsed;

				for(k=0;k<simSize;k++)
					if(outbits[k] != inbits[k])
						errcnt ++;
				printf("\tCalls: %u\n",calls);
				
				//Free some memory...
				free(turbobits);/**/
				free(inbits);
				free(outbits);
				//break;
				
			}

			//Ok we now need to calculate some stats!
			err_rate = ((float)errcnt)/((float)testno*simSize);
			fprintf(fp,"%d,%f,%f,%u\n",it,snr,err_rate,totaltime/testno);
			//fprintf(fp,"%f,%f,%f\n",stability_factor,snr,err_rate);
			fflush(fp);
	//		printf("%f,%f,%f\n",stability_factor,snr,err_rate);
			printf("IT: %d,SNR: %f dB,BER: %e, Time needed: %u ms\n",it,snr,err_rate,totaltime/testno);
			//printf("Theta: %f, SNR: %f, Peb: %f\n",stability_factor,snr,err_rate);
			snr +=0.1f;
		}
	}
//	}
	fclose(fp);

	getc(stdin);
	freeTurboEncoder(&te);
	return 0;
}


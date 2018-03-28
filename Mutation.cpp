#define _CRT_SECURE_NO_WARNINGS
#define _CRT_NONSTDC_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>
#include "mkl.h"

#define N 100

void printBits(size_t const size, void const * const ptr)
{
	unsigned char *b = (unsigned char*) ptr;
	unsigned char byte;
	int i, j;

	for (i=size-1;i>=0;i--)
	{
		for (j=7;j>=0;j--)
		{
			byte = (b[i] >> j) & 1;
			printf("%u", byte);
		}
	}
}

void mutate(double ancestor, double **children, unsigned long long index, unsigned long long amount)
{
//mutation chance:   sign exponent mantissa
	double chance[3]={0.00001, 0.0001, 0.05};
	double byte=0, value=0;
	int expo=0;
	double prob=1;
	unsigned int bits[2][2];
	int ra=0;
	
	for(unsigned long long time=0; time<amount; time++)
	{
		//probability of exponenets
		expo=0;
		prob=RAND_MAX*chance[1];
		for(unsigned long long i=0; i<11; i++) //11 is number of bits for the exponent
		{
			ra=rand();
			expo=expo*2; //bit shift left by 1
			if(ra<=(int)prob)expo++;  //change bit
			if(prob<RAND_MAX)prob=prob*1.23284673944;
		}
		byte=DBL_MIN;
		for(int i=0; i<expo; i++)byte=byte*2;
		if(expo==0)byte=0;
		memcpy(bits[0], &byte, sizeof(byte));
		memcpy(bits[1], &ancestor, sizeof(double));
		//printBits(sizeof(byte), &byte);
		//printf("\n");
		//printBits(sizeof(ancestor[origin]), &ancestor[origin]);
		//printf("\n");
		bits[0][0]=bits[0][0]^bits[1][0];
		bits[0][1]=bits[0][1]^bits[1][1];
		memcpy(&value, bits[0], sizeof(double));
		//printBits(sizeof(descend[origin*N+time]), &descend[origin*N+time]);
		//printf("\n");
		//printf("%15e\n", descend[origin*N+time]);

		//probability of mantissa
		byte=0;
		prob=RAND_MAX*chance[2];
		for(int i=0; i<52; i++) //52 is number of bits for the mantissa
		{
			ra=rand();
			byte=byte*(double)2; //bit shift left by 1
			if(ra<=(int)prob)byte++;  //change bit
			if(prob<RAND_MAX)prob=prob*1.04527549532;
		}
		for(int i=0; i<52; i++)byte=byte/2;
		byte=byte*DBL_MIN;
		memcpy(bits[0], &byte, sizeof(byte));
		memcpy(bits[1], &value, sizeof(double));
		/*printBits(sizeof(descend[origin*N+time]), &descend[origin*N+time]);
		printf("\n");
		printBits(sizeof(byte), &byte);
		printf("\n");*/
		bits[0][0]=bits[0][0]^bits[1][0];
		bits[0][1]=bits[0][1]^bits[1][1];
		memcpy(&value, bits, sizeof(double));
		/*printBits(sizeof(descend[origin*N+time]), &descend[origin*N+time]);
		printf("\n");
		printf("%15e\n", descend[origin*N+time]);*/

		//probability of sign
		prob=RAND_MAX*chance[0];
		ra=rand();
		if(ra<=(int)prob)value=-value;
		*(*children+index+time)=value;
	}
}

void main()
{
	double *ancestor=(double *)mkl_malloc(N*sizeof(double), 128);
	double *descend=(double *)mkl_malloc(N*N*sizeof(double), 128);
	FILE *output=fopen("mutate.txt", "w");

	for(unsigned long long i=0; i<N; i++)ancestor[i]=rand();

	for(unsigned long long i=0; i<N; i++)
	{
		mutate(ancestor[i], &descend, i*N, N);
		for(unsigned long long j=0; j<i+1; j++)
		{
			//printf("%15e %15e\n", ancestor[i], descend[i*N+j]);
			fprintf(output, "%15e %15e\n", ancestor[i], descend[i*N+j]);
		}
		/*for(unsigned long long j=0; j<i; j++)
		{
			printf("%15d %15e\n", i, descend[i*N+j]);
			printBits(sizeof(descend[i*N+j]), &descend[i*N+j]);
			printf("\n");
		}*/
	}
}
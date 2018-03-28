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

#define N 100000

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

void crossover(double parent1, double parent2, double **children, unsigned long long index, unsigned long long amount)
{
	unsigned int bits[2][2];
	unsigned int temp[2];
	int ra=0;
	double value;

	for(unsigned long long time=0; time<amount; time++)
	{
		bits[0][0]=0;
		bits[0][1]=0;
		for(int i=0; i<64; i++)
		{
			bits[0][i/32]=bits[0][i/32]*2;
			if(rand()%2==0)
			{
				bits[0][i/32]++;
			}
		}
		bits[1][0]=4294967295-bits[0][0];
		bits[1][1]=4294967295-bits[0][1];

		//parent1's DNA
		memcpy(temp, &parent1, sizeof(temp));
		bits[0][0]=bits[0][0]&temp[0];
		bits[0][1]=bits[0][1]&temp[1];

		//parent2's DNA
		memcpy(temp, &parent2, sizeof(temp));
		bits[1][0]=bits[1][0]&temp[0];
		bits[1][1]=bits[1][1]&temp[1];

		//combined DNA
		bits[0][0]=bits[0][0]+bits[1][0];
		bits[0][1]=bits[0][1]+bits[1][1];
		memcpy(&value, bits[0], sizeof(value));
		*(*children+time)=value;
	}
}

void main()
{
	FILE *output=fopen("crossover.txt", "w");
	double *ancestor=(double *)mkl_malloc(2*sizeof(double), 128);
	double *descend=(double *)mkl_malloc(N*sizeof(double), 128);

	ancestor[0]=79.161896135;
	ancestor[1]=14.123123;
	crossover(ancestor[0], ancestor[1], &descend, 0, N);  //generate N descendants of ancestor[0] and ancestor[1]
	for(int i=0; i<N; i++)fprintf(output, "%15e %15e %15e\n", ancestor[0], ancestor[1], descend[i]);
}
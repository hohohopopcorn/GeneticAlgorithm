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

#define N 30      //population size
#define GEN 100      //number of generations to run
#define BIRTH 0.75  //percent of children
//mutation chance:   sign exponent mantissa
double chance[3]={   0.02,    0.02,    0.04};
VSLStreamStatePtr randstream;

double fitness(double *gene)
{
	double value=exp(-log(sqrt((gene[0]+3)*(gene[0]+3))*sqrt((gene[0]-2)*(gene[0]-2))+1));
	if(isinf(value))value=-DBL_MAX;
	return value;
}

typedef struct GENE
{
	double strand;
	double strength;
}Gene;

int compare(const void * a, const void * b)
{
	if(*(double *)a>*(double *)b)return 1;
	else if(*(double *)a<*(double *)b)return -1;
	else return 0;  
}
int Genecompare(const void * a, const void * b)
{
	Gene *i=(Gene *)a, *j=(Gene *)b;
	if(i->strength>j->strength)return 1;
	else if(i->strength<j->strength)return -1;
	else return 0;  
}

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

void mutate(double **child, int amount)
{
	double *random=(double *)mkl_malloc(amount*64*sizeof(double), 128);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, randstream, amount*64, random, 0, 1);
	double byte=0, value=0;
	int expo=0;
	unsigned int bits[2][2];
	unsigned long long ra=0;
	double original;

	for(int n=0; n<amount; n++)
	{
		original=*(*child+n);

		//probability of exponenets
		expo=0;
		for(unsigned long long i=0; i<11; i++) //11 is number of bits for the exponent
		{
			expo=expo*2; //bit shift left by 1
			if(random[ra]<=chance[1])expo++;  //change bit
			ra++;
		}
		byte=DBL_MIN;
		for(int i=0; i<expo; i++)byte=byte*2;
		if(expo==0)byte=0;
		memcpy(bits[0], &byte, sizeof(byte));
		memcpy(bits[1], &original, sizeof(double));
		bits[0][0]=bits[0][0]^bits[1][0];
		bits[0][1]=bits[0][1]^bits[1][1];
		memcpy(&value, bits[0], sizeof(double));

		//probability of mantissa
		byte=0;
		for(int i=0; i<52; i++) //52 is number of bits for the mantissa
		{
			byte=byte*(double)2; //bit shift left by 1
			if(random[ra]<=chance[2])byte++;  //change bit
			ra++;
		}
		for(int i=0; i<52; i++)byte=byte/2;
		byte=byte*DBL_MIN;
		memcpy(bits[0], &byte, sizeof(byte));
		memcpy(bits[1], &value, sizeof(double));
		bits[0][0]=bits[0][0]^bits[1][0];
		bits[0][1]=bits[0][1]^bits[1][1];
		memcpy(&value, bits, sizeof(double));

		//probability of sign
		if(random[ra]<=chance[0])value=-value;
		ra++;

		*(*child+n)=value;
	}
}

void crossover(double parent1, double parent2, double **children, unsigned long long index, unsigned long long amount)
{
	double *random=(double *)mkl_malloc(amount*64*sizeof(double), 128);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, randstream, amount*64, random, 0, 1);
	int ra=0;
	unsigned int bits[2][2];
	unsigned int temp[2];
	double value;

	for(unsigned long long time=0; time<amount; time++)
	{
		bits[0][0]=0;
		bits[0][1]=0;
		for(int i=0; i<64; i++)
		{
			bits[0][i/32]=bits[0][i/32]*2;
			if(random[ra]<=0.5)bits[0][i/32]++;
			ra++;
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
		*(*children+index+time)=value;
	}
}

void selection(unsigned long long **selected, unsigned long long **notselected, Gene *ancestors)
{
	double *random=(double *)mkl_malloc((N-(int)(N*(double)BIRTH))*sizeof(double), 128);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, randstream, (N-(int)(N*(double)BIRTH)), random, 0, 1);
	double prob, total=0;
	char chosen[N];
	memset(chosen, 0, sizeof(chosen));
	
	unsigned long long j=0;
	for(unsigned long long i=0; i<N; i++)total+=ancestors[i].strength*ancestors[i].strength;
	for(unsigned long long i=0; i<N-(int)(N*(double)BIRTH); i++)
	{
		prob=random[i]*total;
		for(j=0; j<N-1; j++)
		{
			if(chosen[j]==0)prob-=ancestors[j].strength*ancestors[j].strength;
			if(prob<0)break;
		}
		chosen[j]=1;
		*(*notselected+i)=j;
		total-=ancestors[j].strength*ancestors[j].strength;
	}
	j=0;
	for(unsigned long long i=0; i<N&&j<(int)(N*(double)BIRTH); i++)
	{
		if(chosen[i]==0)
		{
			*(*selected+j)=i;
			j++;
		}
	}
}

void reproduction(unsigned long long *notselected, Gene *ancestors, double **descend)
{
	double *random=(double *)mkl_malloc(2*(int)(N*(double)BIRTH)*sizeof(double), 128);
	unsigned long long p1, p2; //parents
	double prob, total=0;

	for(unsigned long long i=0; i<(N-(int)(N*(double)BIRTH)); i++)total+=ancestors[notselected[i]].strength*ancestors[notselected[i]].strength;
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, randstream, 2*(int)(N*(double)BIRTH), random, 0, total);

	for(unsigned long long i=0; i<(int)(N*(double)BIRTH); i++)
	{
		prob=0;
		for(p1=0; prob<random[2*i]&&p1<N-(int)(N*(double)BIRTH)-1; p1++)prob+=ancestors[notselected[p1]].strength*ancestors[notselected[p1]].strength;
		prob=0;
		for(p2=0; prob<random[2*i+1]&&p2<N-(int)(N*(double)BIRTH)-1; p2++)prob+=ancestors[notselected[p2]].strength*ancestors[notselected[p2]].strength;

		crossover(ancestors[notselected[p1]].strand, ancestors[notselected[p2]].strand, descend, i, 1);
	}
	mutate(descend, (int)(N*(double)BIRTH));
}

void main()
{
	vslNewStream(&randstream,VSL_BRNG_MCG31, 1);
	double *random=(double *)mkl_malloc(2*N*sizeof(double), 128);
	Gene *ancestor=(Gene *)mkl_malloc(N*sizeof(Gene), 128);
	double *descend=(double *)mkl_malloc((int)(N*(double)BIRTH)*sizeof(double), 128);
	FILE *output=fopen("generations.txt", "w");  //cost history
	unsigned long long *selected=(unsigned long long *)mkl_malloc((int)(N*(double)BIRTH)*sizeof(unsigned long long), 128);
	unsigned long long *notselected=(unsigned long long *)mkl_malloc((N-(int)(N*(double)BIRTH))*sizeof(unsigned long long), 128);

	for(long long i=0; i<N; i++)ancestor[i].strand=rand();
	for(long long i=0; i<N; i++)ancestor[i].strength=fitness(&ancestor[i].strand);
	qsort(ancestor, N, sizeof(Gene), Genecompare);
	//for(int i=0; i<N; i++)fprintf(output, "%8d %15d %15e %15e\n", 0, i+1, ancestor[i].strand, ancestor[i].strength);
	fprintf(output, "%8d %15d %15e %15e\n", 0, N, ancestor[N-1].strand, ancestor[N-1].strength);

	for(int i=0; i<GEN; i++)
	{
		if((i+1)%10==0)printf("%d\n", i+1);
		selection(&selected, &notselected, ancestor);  //the selected are replaced, the not selected survive
		reproduction(notselected, ancestor, &descend);
		for(int j=0; j<(int)(N*(double)BIRTH); j++)
		{
			ancestor[selected[j]].strand=descend[j];  //BIRTH percent of population is replaced
		}

		for(long long j=0; j<N; j++)ancestor[j].strength=fitness(&ancestor[j].strand);
		qsort(ancestor, N, sizeof(Gene), Genecompare);
		//for(int j=0; j<N; j++)fprintf(output, "%8d %15d %15e %15e\n", i+1, j+1, ancestor[j].strand, ancestor[j].strength);
		fprintf(output, "%8d %15d %15e %15e\n", i+1, N, ancestor[N-1].strand, ancestor[N-1].strength);
	}
	mkl_free(ancestor);
	mkl_free(descend);
	fclose(output);
}
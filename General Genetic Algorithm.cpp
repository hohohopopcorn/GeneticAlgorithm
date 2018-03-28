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

#define N 50          //population size
#define BYTE 1000        //number of bytes needed to represent the chromosome
#define GEN 1000000       //number of generations to run
#define BIRTH 0.75    //percent of children
//mutation chance: not separated by sign exponent mantissa
double change=(double)1/N;
double jump=(double)1/BYTE;
double jumpwidth=50;
double temperature=1;
VSLStreamStatePtr randstream;

typedef struct GENE
{
	double strand[BYTE];
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

void fitness(GENE **species)
{
	double *gene;
	double res, value=0;
	
	for(long long i=0; i<N; i++)
	{
		value=0;
		gene=((*species+i)->strand);

		for(long long j=0; j<BYTE; j++)
		{
			res=gene[j];
			value+=cos(res-j)*cos(res-j)/((res-j)*(res-j)+1);
			if(isinf(value)||isnan(value))
			{
				value=0;
				break;
			}
		}

		(*species+i)->strength=value;
	}

	qsort(*species, N, sizeof(Gene), Genecompare);
}

void mutate(double **child, unsigned long long amount)
{
	double *random=(double *)mkl_malloc(amount*BYTE*sizeof(double), 64);
	double *energy=(double *)mkl_malloc(amount*BYTE*sizeof(double), 64);
	double *mut=(double *)mkl_malloc(amount*BYTE*sizeof(double), 64);
	vdRngWeibull(VSL_RNG_METHOD_WEIBULL_ICDF, randstream, amount*BYTE, energy, 2, 0, temperature*jumpwidth*sqrt(2));
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, randstream, amount*BYTE, mut, 0, 1);
	vdMul(amount*BYTE, energy, mut, mut);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, randstream, BYTE, random, jump, 1.0+jump);
	vdFloor(amount*BYTE, random, random);

	//mutate by gaussian distribution
	vdMul(amount*BYTE, random, mut, mut);
	for(unsigned long long n=0; n<amount*BYTE; n++)mut[n]=((double)random[n])*mut[n];
	vdAdd(amount*BYTE, *child, mut, *child);

	mkl_free(energy);
	mkl_free(mut);
	mkl_free(random);

	//mutate by changing gene position
	unsigned long long a, b, option;
	unsigned long long replace;
	double temp[BYTE];
	random=(double *)mkl_malloc(amount*4*sizeof(double), 64);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, randstream, 4*amount, random, 0, 1);

	for(unsigned long long n=0; n<amount; n++)
	{
		cblas_dcopy(BYTE, *child+n*BYTE, 1, temp, 1);
		option=rand()%4;

		if(random[n*4]<change)       //Swap Position
		{
			a=rand()%(BYTE);
			b=rand()%(BYTE);
			*(*child+n*BYTE+a)=temp[b];
			*(*child+n*BYTE+b)=temp[a];
		}
		if(random[n*4+1]<change)  //Scramble Position
		{
			double store;
			a=rand()%(BYTE-1);
			b=rand()%(BYTE-a-1)+1;
			for(unsigned long long i=0; i<b; i++)
			{
				replace=rand()%(b-i);
				store=temp[replace+a+1+i];
				temp[replace+a+1+i]=temp[a+i];
				temp[a+i]=store;
			}
			cblas_dcopy(BYTE, temp, 1, *child+n*BYTE, 1);
		}
		if(random[n*4+2]<change)  //Displace Position
		{
			a=rand()%(BYTE);
			b=rand()%(BYTE-a);
			replace=rand()%(BYTE-b-1);
			for(unsigned long long i=a; i<=replace; i++)*(*child+n*BYTE+i)=temp[i+b+1];
			for(unsigned long long i=replace; i<replace+b+1; i++)*(*child+n*BYTE+i)=temp[a+i-replace];
			for(unsigned long long i=replace; i<a; i++)*(*child+n*BYTE+b+1+i)=temp[i];
		}
		if(random[n*4+3]<change)  //Inverse Position
		{
			a=rand()%(BYTE-1);
			b=rand()%(BYTE-a-1)+1;
			for(unsigned long long i=0; i<(b+1)/2; i++)*(*child+n*BYTE+a+i)=temp[a+b-i];
		}
	}

	mkl_free(random);
}

void crossover(unsigned long long parent1, unsigned long long parent2, Gene *ancestors, double **child, unsigned long long target)
{
	double *random=(double *)mkl_malloc(BYTE*sizeof(double), 64);
	double *part1=(double *)mkl_malloc(BYTE*sizeof(double), 64);
	double *part2=(double *)mkl_malloc(BYTE*sizeof(double), 64);
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, randstream, BYTE, random, -0.5, 0.5);  //cut two parent DNA randomly
	
	vdFloor(BYTE, random, part1);                            //-1 or 0
	vdCeil(BYTE, random, part2);                             // 1 or 0
	vdMul(BYTE, part1, (ancestors+parent1)->strand, part1);  //negative parts of parent1
	vdMul(BYTE, part2, (ancestors+parent2)->strand, part2);  //positive parts of parent2
	vdSub(BYTE, part2, part1, *child+target*BYTE);           //combine DNA
	
	mkl_free(random);
	mkl_free(part1);
	mkl_free(part2);
}

void selection(unsigned long long **selected, unsigned long long **notselected, Gene *ancestors)
{
	double *random=(double *)mkl_malloc(N*sizeof(double), 64);
	double *probability=(double *)mkl_malloc(N*sizeof(double), 64);
	double total=0;
	char chosen[N];
	memset(chosen, 0, sizeof(chosen));

	for(unsigned long long i=0; i<N; i++)probability[i]=(ancestors[N-1].strength-ancestors[i].strength)/ancestors[N-1].strength;
	vdErfc(N, probability, probability);
	for(unsigned long long i=0; i<N; i++)total+=probability[i];
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, randstream, N, random, 0, total);

	int n=0;
	unsigned long long j=0;
	for(unsigned long long i=0; i<(int)(N*(double)BIRTH); i=i)
	{
		while(chosen[j]==1)j=(j+1)%N;
		if(random[n]>probability[j])
		{
			chosen[j]=1;
			*(*selected+i)=j;
			i++;
		}
		j=(j+1)%N;
		n=(n+1)%N;
	}

	j=0;
	for(unsigned long long i=N-1; i>=0&&j<N-(int)(N*(double)BIRTH); i--)
	{
		if(chosen[i]==0)
		{
			*(*notselected+j)=i;
			j++;
		}
	}
	mkl_free(random);
	mkl_free(probability);
}

void reproduction(unsigned long long *notselected, Gene *ancestors, double **descend)
{
	double *random=(double *)mkl_malloc(2*(int)(N*(double)BIRTH)*sizeof(double), 64);
	double *probability=(double *)mkl_malloc((N-(int)(N*(double)BIRTH))*sizeof(double), 64);
	unsigned long long p1, p2; //parents
	double prob, total=0;

	for(unsigned long long i=0; i<N-(int)(N*(double)BIRTH); i++)
		probability[i]=(ancestors[notselected[N-(int)(N*(double)BIRTH)-1]].strength-ancestors[notselected[i]].strength)
		/ancestors[notselected[N-(int)(N*(double)BIRTH)-1]].strength;
	vdErfc(N-(int)(N*(double)BIRTH), probability, probability);
	for(unsigned long long i=0; i<N-(int)(N*(double)BIRTH); i++)total+=probability[i];
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, randstream, 2*(int)(N*(double)BIRTH), random, 0, total);

	for(unsigned long long i=0; i<(int)(N*(double)BIRTH); i++)
	{
		prob=0;
		for(p1=0; prob<random[2*i]&&p1<N-(int)(N*(double)BIRTH)-1; p1++)prob+=probability[p1];
		prob=0;
		for(p2=0; prob<random[2*i+1]&&p2<N-(int)(N*(double)BIRTH)-1; p2++)prob+=probability[p2];

		crossover(notselected[p1], notselected[p2], ancestors, descend, i);
	}
	mutate(descend, (int)(N*(double)BIRTH));
	mkl_free(random);
	mkl_free(probability);
}

void main()
{
	vslNewStream(&randstream,VSL_BRNG_MCG31, 1);
	double *random=(double*)mkl_malloc(N*BYTE*sizeof(double), 64);
	Gene *ancestor=(Gene *)mkl_malloc(N*sizeof(Gene), 64*(BYTE+1));
	double *descend=(double *)mkl_malloc(N*BYTE*sizeof(double), 64);
	FILE *output=fopen("generations.txt", "w");  //cost history
	unsigned long long *selected=(unsigned long long *)mkl_malloc((int)(N*(double)BIRTH)*sizeof(unsigned long long), 64);
	unsigned long long *notselected=(unsigned long long *)mkl_malloc((N-(int)(N*(double)BIRTH))*sizeof(unsigned long long), 64);

	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, randstream, N*BYTE, random, -1e3, 1e3);
	for(long long i=0; i<N; i++)cblas_dcopy(BYTE, random+i*BYTE, 1, ancestor[i].strand, 1);
	mkl_free(random);
	fitness(&ancestor);
	for(int i=N-1; i<N; i++)
	{
		fprintf(output, "%15e %15d ", (double)0, i+1);
		//for(int j=0; j<BYTE; j++)fprintf(output, "%15e ", ancestor[i].strand[j]);
		fprintf(output, "%15e\n", ancestor[i].strength);
	}
	fclose(output);

	clock_t t1=clock();
	for(unsigned long long i=0; i<GEN; i++)
	{
		temperature=1.0-(double)i/(double)GEN*(1.0-1e-4);
		selection(&selected, &notselected, ancestor);  //the selected are replaced, the not selected survive
		reproduction(notselected, ancestor, &descend);
		for(int j=0; j<(int)(N*(double)BIRTH); j++)
		{
			cblas_dcopy(BYTE, descend+j*BYTE, 1, (ancestor+selected[j])->strand, 1);
		}

		fitness(&ancestor);
		if((i+1)%(int)(GEN/10000)==0)output=fopen("generations.txt", "a+");  //cost history
		for(int k=N-1; k<N&&(i+1)%(int)(GEN/10000)==0; k++)
		{
			fprintf(output, "%15e %15d ", (double)i+1, k+1);
			//for(int j=0; j<BYTE; j++)fprintf(output, "%15e ", ancestor[k].strand[j]);
			fprintf(output, "%15e\n", ancestor[k].strength);
		}
		if((i+1)%(int)(GEN/10000)==0)fclose(output);  //cost history
		if((i+1)%(int)(GEN/1000)==0)printf("%15e %15e %15e %15e\n", (double)i+1, ancestor[N-1].strength, jumpwidth*temperature, (double)(clock()-t1)/(double)(i+1));
	}
	mkl_free(ancestor);
	mkl_free(descend);
	mkl_free(selected);
	mkl_free(notselected);
	fclose(output);
}
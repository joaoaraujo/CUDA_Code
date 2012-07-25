// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <map>
#include <string>
#include <iostream>
#include <fstream>

#include <./include/CkCrypt2.h>
// includes, project
#include <cutil_inline.h>

#include "driver_types.h"
#include "cuda_runtime.h"

#define NUMEVENTS 32
#define NUMATREVENTS 45

__device__ __constant__ int3 events_d[NUMEVENTS*NUMATREVENTS];//Eventos
__device__ __constant__ int2 sizeEvents_d[NUMEVENTS];//Tamanho de cada Evento
__device__ __constant__ float4 valuesEvents_d[NUMEVENTS*NUMATREVENTS*2];//Tamanho de cada Evento

// includes, kernels
#include "mem_info.h"
#include <kernels_1_Event.cu>
#include <kernels_32_Events.cu>
#include <kernels_InsRemSubs.cu>

using namespace std;

#define ACTIVE_BIT(i, j) (i & (1 << j))
#define EDIT_BIT(i, j) (i | (1 << j))

//Host
//#### AGREGADOS ####
unsigned int** clusters;
int numberClusters;

unsigned int* specialCluster;
//#### ATRIBUTOS ####
int numberAtributes;

//map<string, int> atributesRef;
int numberValues;
//map<string, int> words;

//#### SUBSCRIÇÕES ####
int numberSubs;
int numberSubsDisj;
map<int, int3> subscriptions;

int** newSubs;
int* numNewSubs;

map<int, int3> deletedSubs;
int* removedSubs;

//#### PREDICADOS ####
int numberPredicates;
//map<string, int> predicatesRef;
int3* predicates;
float4* valuesPrd;
//#### EVENTOS ####
int3* events;
float4* values;
int2* sizeEvents;
//#### CONFIGURAÇÕES ####
int totalthreads;
int num_blocks;
int num_threadsperblock;
int ntInit;
int ntFinal;
int nbInit;
int nbFinal;
int mulInit;
int mulFinal;
int bit;
int numIter;
int usingshared;
int kernelnum;
double percPerC;
int constant;

//Device
int* d_totalthreads;

__device__ unsigned int* bitmap;

int** d_clusters;
int* d_sizes;

int3* events_dev;//Eventos
int2* sizeEvents_dev;//Tamanho de cada Evento
float4* values_d;//Valores dos Eventos

int3* predicates_d;//Predicados
float4* valuesPrd_d;//Valores dos Predicados
int* numPred_d;//Número de Predicados

//############################ TESTS ######################################
void quickSort(unsigned int arr[], unsigned int left, unsigned int right)
{
	unsigned int i = left, j = right;

	int tmp;

	unsigned int pivot = arr[(left + right) / 2];

	/* partition */

	while (i <= j) {

		while (arr[i] < pivot)

			i++;

		while (arr[j] > pivot)

			j--;

		if (i <= j) {

			tmp = arr[i];

			arr[i] = arr[j];

			arr[j] = tmp;

			i++;

			j--;

		}

	};

	/* recursion */

	if (left < j)

		quickSort(arr, left, j);

	if (i < right)

		quickSort(arr, i, right);

}

void quickSortArrays(unsigned int arr[],unsigned int left,unsigned int right,int poww)
{
	int i = left, j = right;

	int k/* = 1*/;

	int* pivot = (int*)malloc(sizeof(int)*(poww+1));
	memcpy(pivot, &arr[((left + right) / 2)*(poww+1)], sizeof(int)*(poww+1));
	int value = pivot[1];
	/* partition */
	while (i <= j) 
	{	
		while(1)
		{
			k=1;
			while (arr[(i*(poww+1))+1] < (unsigned int)value)
				i++;

			if(arr[(i*(poww+1))+1] > (unsigned int)value)
				break;
			else if(arr[(i*(poww+1))+1] == (unsigned int)value){

				while(arr[(i*(poww+1))+k] == (unsigned int)pivot[k]){
					k++;
					if(k > poww)
						break;
				}
				if(k > poww)
					break;

				if(arr[(i*(poww+1))+k] > (unsigned int)pivot[k])
					break;
				else
					i++;
			}
		}

		while(1)
		{
			k = 1;

			while (arr[(j*(poww+1))+1] > (unsigned int)value)
				j--;

			if(arr[(j*(poww+1))+1] < (unsigned int)value )
				break;
			else if(arr[(j*(poww+1))+1] == (unsigned int)value ){

				while(arr[(j*(poww+1))+k] == pivot[k]){
					k++;

					if(k > poww)
						break;
				}
				if(k > poww)
					break;
				if(arr[(j*(poww+1))+k] < (unsigned int)pivot[k])
					break;
				else
					j--;

			}
		}

		if (i <= j) {
			int* tmp = (int*)malloc(sizeof(int)*(poww+1));
			memcpy(tmp, &arr[(i*(poww+1))],sizeof(int)*(poww+1));
			memcpy(&arr[(i*(poww+1))], &arr[(j*(poww+1))],sizeof(int)*(poww+1));
			memcpy(&arr[(j*(poww+1))], tmp,sizeof(int)*(poww+1));
			free(tmp);

			i++;
			j--;
		}
	};
	free(pivot);
	/* recursion */
	//  free(tmp);
	if (left < (unsigned int)j)

		quickSortArrays(arr, left, j, poww);

	if ((unsigned int)i < right)

		quickSortArrays(arr, i, right,poww);
}

unsigned int randomNumber(int maxV)
{
	return rand() % maxV;
}

unsigned int** randomClusters(int numSubs, int numClus, int numPred, int qs, double perc)
{
	unsigned int** r = (unsigned int**)malloc(sizeof(int)*(numClus+1));
	int i = 1;
	int init = numSubs;
	int tSubs = 0;
	int idSubs = 0;
	int init2;

	srand ( 10000 );

	r[0] = (unsigned int*)malloc(sizeof(int)*(numClus+1));

	while(i != (numClus+1))
	{
		if(i!=numClus)
			init2 = (((int)(init * perc)) >> 5) << 5;
		else
			init2 = (int)(numSubs-tSubs);

		tSubs+=init2;

		int poww = (int)pow(2.0, i);

		cout << init2 << endl;

		r[0][i] = init2*(poww+1);
		int currS = 0;

		r[i] = (unsigned int*)malloc(sizeof(int)*init2*(poww+1));

		while(currS != init2)
		{
			unsigned int* preds = (unsigned int*)malloc(sizeof(int)*(poww+1));
			preds[0] = idSubs;

			for(int j = 1; j <= poww; j++){

				int found = 1;
				unsigned int num;

				while(found == 1)
				{
					num = randomNumber(numPred);
					
					for(int k = 1; k < j; k++)
					{
						if(preds[k] == num)
						{
							found = 2;
							break;
						}
					}

					if(found == 2)
						found = 1;
					else
						found = 0;
				}
				preds[j] = num;
			}

			if(qs)
				quickSort(preds, 1, poww);

			subscriptions[idSubs].x = idSubs; 
			subscriptions[idSubs].y = i;
			subscriptions[idSubs].z = currS*(poww+1);

			memcpy(&r[i][currS*(poww+1)],preds,(poww+1)*sizeof(int));

			currS++;

			idSubs++;

			free(preds);

		}

		init-=init2;

		if(qs)
			quickSortArrays(r[i], 0, r[0][i]/(poww+1)-1, poww);

		i++;
	}
	return r;

}

unsigned int* give32Bit()
{
	unsigned int* ret = (unsigned int*)malloc(sizeof(int)*numberPredicates);
	srand ( 4294967295 );//max 32767

	for(int i = 0; i < numberPredicates; i++)
		ret[i] = (unsigned int)(rand() << 17 | rand() << 2 | rand() >> 13);

	return ret;
}

unsigned int* randomPos(double perc)
{
	srand ( 10000 );

	int size2 = static_cast<int>(numberPredicates * perc);

	//cout << size2 << endl;

	unsigned int* randomPositions = (unsigned int*)malloc(sizeof(int)*size2);

	unsigned int num;

	randomPositions[0] = randomNumber(numberPredicates);
	
		for(int i = 1; i < size2 ; i++)
		{
				int found = 1;

				while(found == 1)
				{
					num = randomNumber(numberPredicates);
					for(int j = 0; j < i; j++)
					{
						if(randomPositions[j] == num)
						{
							found = 2;
							break;
						}
					}
					if(found == 2)
						found = 1;
					else
						found = 0;
				}

				randomPositions[i] = num;
		}

	return randomPositions;

}

unsigned int** randomPercAccepted(double perc)
{
	unsigned int** ret = (unsigned int**)malloc(sizeof(unsigned int)*32);

	for(int i = 0; i < 32; i++)
		ret[i] = randomPos(perc);

	return ret;
}
//############################ PREDICATES ######################################
void initPred()
{
	srand ( 10000 );

	int pos2 = 0;
	for(int j = 0; j < numberPredicates; j++)
	{
		pos2 = j*2;
		predicates[j].x = j;
		predicates[j].y = (int)randomNumber(7);
		predicates[j].z = pos2;
	}
	int val;
	for(int i = 0; i < (numberPredicates*2);i+=2)
	{
		val = randomNumber(numberPredicates);
		valuesPrd[i].x = (float)val;
		valuesPrd[i].y = (float)val;
		valuesPrd[i].z = (float)val;
		valuesPrd[i].w = (float)val;

		valuesPrd[i+1].x = (float)val;
		valuesPrd[i+1].y = (float)val;
		valuesPrd[i+1].z = (float)val;
		valuesPrd[i+1].w = (float)val;
		/*valuesPrd[i].x = (float)i;
		valuesPrd[i].y = (float)i;
		valuesPrd[i].z = (float)i;
		valuesPrd[i].w = (float)i;

		valuesPrd[i+1].x = (float)i;
		valuesPrd[i+1].y = (float)i;
		valuesPrd[i+1].z = (float)i;
		valuesPrd[i+1].w = (float)i;*/
	}
}


//############################ EVENTS ######################################
void initEvents(int numberAtr, int numEvents)
{
	srand ( 10000 );

	int position = 0;

	for(int i = 0; i < numEvents; i++)
	{
		sizeEvents[i].x = position;

		int pos = 0;
		for(int j = 0; j < numberAtr; j++)
		{
			events[i*numberAtr + j].x = (int)randomNumber(numberPredicates); 
			events[i*numberAtr + j].y = i*(numberAtr*2)+pos;
			events[i*numberAtr + j].z = i*(numberAtr*2)+(pos+1);
			pos+=2;
		}
		position += numberAtr;
		sizeEvents[i].y = position;
	}

	for(int j = 0; j < (numEvents*numberAtr)*2; j+=2)
	{
		int ran = randomNumber(numberPredicates*2);
		values[j].x = (float)ran;
		values[j].y = (float)ran;
		values[j].z = (float)ran;
		values[j].w = (float)ran;

		values[j+1].x = (float)ran;
		values[j+1].y = (float)ran;
		values[j+1].z = (float)ran;
		values[j+1].w = (float)ran;
	}

}

//############################ BITMAP ######################################
void generateRandomBitmap(double perc)
{
		unsigned int* r = randomPos(perc);
		
		int s = (int)(numberPredicates*perc);

		int* size;
		cutilSafeCall(cudaMalloc((void**) &size, sizeof(int)));
		cutilSafeCall(cudaMemcpy(size, &s, sizeof(int),cudaMemcpyHostToDevice));

		int* rpos;
		cutilSafeCall(cudaMalloc((void**) &rpos, sizeof(int)*s));
		cutilSafeCall(cudaMemcpy(rpos, r,sizeof(int)*s,cudaMemcpyHostToDevice));
		
		int* d_bitsize;

		if(usingshared)
		{
			int bitsize = (int)ceil((numberPredicates+1)/32.0);
			cutilSafeCall(cudaMalloc((void**) &bitmap, sizeof(unsigned int)*bitsize));
			cutilSafeCall(cudaMalloc((void**) &d_bitsize, sizeof(int)));
			cutilSafeCall(cudaMemcpy(d_bitsize, &bitsize,sizeof(int),cudaMemcpyHostToDevice));	
			kernelFillBit1EventShared<<<1, 1>>>(size, rpos, bitmap);
		}
		else{
			cutilSafeCall(cudaMalloc((void**) &bitmap, sizeof(unsigned int)*(numberPredicates+1)));
			kernelFillBit1Event<<<1, 1>>>(size, rpos, bitmap);
		}
		cudaThreadSynchronize();

		free(r);
		cutilSafeCall(cudaFree(rpos));
		cutilSafeCall(cudaFree(size));
}

void generateRandomBitmap32(double perc)
{
	unsigned int** r = randomPercAccepted(perc);
	
	int s = (int)(numberPredicates*perc);

	int* size;
	cutilSafeCall(cudaMalloc((void**) &size, sizeof(int)));
	cutilSafeCall(cudaMemcpy(size, &s, sizeof(int),cudaMemcpyHostToDevice));

	unsigned int** rpos;
	cutilSafeCall(cudaMalloc((void**) &rpos, sizeof(unsigned int)*s));
	cutilSafeCall(cudaMemcpy(rpos, r,sizeof(unsigned int)*32,cudaMemcpyHostToDevice));

	for(int i = 0; i < 32; i++)
	{
		cutilSafeCall(cudaMalloc((void**) &rpos[i], sizeof(int)*s));
		cutilSafeCall(cudaMemcpy(rpos[i], r[i],sizeof(int)*s,cudaMemcpyHostToDevice));
	}

	cutilSafeCall(cudaMalloc((void**) &bitmap, sizeof(unsigned int)*(numberPredicates+1)));

	int* eNumber;
	cutilSafeCall(cudaMalloc((void**) &eNumber, sizeof(int)));

	for(int i = 0; i < 32; i++)
	{
		cutilSafeCall(cudaMemcpy(eNumber, &i,sizeof(int),cudaMemcpyHostToDevice));
		kernelFillBit32EventShared<<<1, 1>>>( size, eNumber, rpos[i], bitmap);
	}	
	cudaThreadSynchronize();
	
	for(int i = 0; i < 32; i++)
	{
		cutilSafeCall(cudaFree(rpos[i]));
		free(r[i]);
	}

	free(r);
	cutilSafeCall(cudaFree(rpos));
	cutilSafeCall(cudaFree(eNumber));
	cutilSafeCall(cudaFree(size));

}

void generateBitmapVar(int numEvents)
{
	ofstream f;
	f.open ("Time-Events.txt",ios::app|std::ios::out);

	unsigned int timer;
	cutilCheckError( cutCreateTimer(&timer));
	//f << "Number of predicates = " << numberPredicates << endl;
	f << "######################### EVENTS " << numEvents << " #########################" << endl;
	f << "######################### PREDICATES " << numberPredicates << " #########################" << endl;

			int numblocks= 6;
			int numthreads= 128;
		
			int totthreads = numblocks*numthreads;
			int* tt_d;

			cutilSafeCall(cudaMalloc((void**) &tt_d, sizeof(int)));
			cutilSafeCall(cudaMemcpy(tt_d, &totthreads,sizeof(int),cudaMemcpyHostToDevice));	

			for(int i = 1; i <= 45; i++)
			{
				/*if(i%10==0)
					cout << i << endl;*/

				int numberAtr = i;
				double totaltime = 0.0;

				//cutStartTimer(timer);

				events = (int3*)malloc(sizeof(int3)*(numEvents*numberAtr));
				values = (float4*)malloc(sizeof(float4)*(numEvents*numberAtr)*2);
				sizeEvents = (int2*)malloc(sizeof(int2)*numEvents);

				predicates = (int3*)malloc(sizeof(int3)*numberPredicates);
				valuesPrd = (float4*)malloc(sizeof(float4)*(numberPredicates*2));

				initEvents(numberAtr, numEvents);
				initPred();

				cutilSafeCall(cudaMalloc((void**) &predicates_d, sizeof(int3)*numberPredicates));
				cutilSafeCall(cudaMalloc((void**) &valuesPrd_d, sizeof(float4)*(numberPredicates*2)));

				if(constant){
					CUDA_SAFE_CALL(cudaMemcpyToSymbol(events_d, events,sizeof(int3)*(numEvents*numberAtr)));
					CUDA_SAFE_CALL(cudaMemcpyToSymbol(sizeEvents_d, sizeEvents, sizeof(int2)*numEvents));
					CUDA_SAFE_CALL(cudaMemcpyToSymbol(valuesEvents_d, values, sizeof(float4)*((numEvents*numberAtr)*2)));
				}
				else
				{
					cutilSafeCall(cudaMalloc((void**) &events_dev, sizeof(int3)*(numEvents*numberAtr)));
					cutilSafeCall(cudaMalloc((void**) &sizeEvents_dev, sizeof(int2)*numEvents));
					cutilSafeCall(cudaMalloc((void**) &values_d, sizeof(float4)*((numEvents*numberAtr)*2)));

					cutilSafeCall(cudaMemcpy(values_d, values,sizeof(float4)*((numEvents*numberAtr)*2),cudaMemcpyHostToDevice));
					cutilSafeCall(cudaMemcpy(events_dev, events,sizeof(int3)*(numEvents*numberAtr),cudaMemcpyHostToDevice));
					cutilSafeCall(cudaMemcpy(sizeEvents_dev, sizeEvents, sizeof(int2)*numEvents, cudaMemcpyHostToDevice));
				}
				
				cutilSafeCall(cudaMemcpy(predicates_d, predicates,sizeof(int3)*numberPredicates,cudaMemcpyHostToDevice));		
				cutilSafeCall(cudaMemcpy(valuesPrd_d, valuesPrd,sizeof(float4)*(numberPredicates*2),cudaMemcpyHostToDevice));				


				//	cout << numblocks << " " << numthreads << " part = " << usingshared << " pred = " << numberPredicates << " e = " << numEvents << endl;

					for(int k = 0; k < 100; k++)
					{
						/*if(k%10 == 0)
							cout << k << endl;*/

						if(!usingshared || numEvents == 32)
							cutilSafeCall(cudaMalloc((void**) &bitmap, sizeof(unsigned int)*(numberPredicates+1)));
						else
							cutilSafeCall(cudaMalloc((void**) &bitmap, sizeof(unsigned int)*((int)ceil((numberPredicates+1)/32.0))));

						if(numEvents == 32){
							if(!usingshared && !constant)
							{
								cutStartTimer(timer);
								kernelBit32Events<<<numblocks, numthreads>>>(tt_d, predicates_d, numPred_d, valuesPrd_d, events_dev, sizeEvents_dev, values_d, bitmap);
								cudaThreadSynchronize();
								cutStopTimer(timer);
							}else if(usingshared && !constant){
								cutStartTimer(timer);
								kernelBit32EventsShared<<<numblocks, numthreads, sizeof(float4)*((numEvents*numberAtr)*2/* 10 KB = 16*640/1024 */)>>>(tt_d, predicates_d, numPred_d, valuesPrd_d, events_dev, sizeEvents_dev, values_d, bitmap);
								cudaThreadSynchronize();
								cutStopTimer(timer);
							}else
							{
								cutStartTimer(timer);
								kernelBit32EventsConstant<<<numblocks, numthreads>>>(tt_d, predicates_d, numPred_d, valuesPrd_d, bitmap);
								cudaThreadSynchronize();
								cutStopTimer(timer);
							}
						}
						else
						{
							if(!usingshared && !constant){
								cutStartTimer(timer);
								kernelBit1Event<<<numblocks, numthreads>>>(tt_d, predicates_d, numPred_d, valuesPrd_d, events_dev, sizeEvents_dev, values_d, bitmap);
								cudaThreadSynchronize();
								cutStopTimer(timer);
							}else if(usingshared && !constant){
								cutStartTimer(timer);
								kernelBit1EventShared<<<numblocks, numthreads, sizeof(unsigned int)*((int)ceil((numberPredicates+1)/32.0))>>>(tt_d, predicates_d, numPred_d, valuesPrd_d, events_dev, sizeEvents_dev, values_d, bitmap);
								cudaThreadSynchronize();
								cutStopTimer(timer);
							}else
							{
								cutStartTimer(timer);
								kernelBit1EventsConstant<<<numblocks, numthreads>>>(tt_d, predicates_d, numPred_d, valuesPrd_d, bitmap);
								cudaThreadSynchronize();
								cutStopTimer(timer);
							}
						}
						totaltime+=cutGetTimerValue(timer)/(double)CLOCKS_PER_SEC;
						cutResetTimer(timer);
						cutilSafeCall(cudaFree(bitmap));
					}

					
			/*		cudaError_t err = cudaGetLastError();
					printf("Cuda error on bitmap: %s.\n", cudaGetErrorString(err) );	*/
				
				free(events);
				free(values);
				free(sizeEvents);
				free(predicates);
				free(valuesPrd);
				
				if(!constant){
					cutilSafeCall(cudaFree(events_dev));
					cutilSafeCall(cudaFree(sizeEvents_dev));
					cutilSafeCall(cudaFree(values_d));
				}
				cutilSafeCall(cudaFree(predicates_d));
				cutilSafeCall(cudaFree(valuesPrd_d));

				f << totaltime << endl;
		}

	cutilSafeCall(cudaFree(tt_d));
}


void generateBitmap(int numblocks, int numthreads, int numEvents)
{
	//ofstream f;
	//f.open ("Time-Events.txt",ios::app|std::ios::out);

	int filename = numberPredicates+usingshared;
	ofstream f;
	char str[20];
	itoa(filename,str,10);
	f.open (str,ios::app|std::ios::out);

	int numberAtr = 10;
	unsigned int timer;
	cutilCheckError( cutCreateTimer(&timer));
	f << "Number of predicates = " << numberPredicates << endl;
	f << "######################### EVENTS " << numEvents << " #########################" << endl;


				//cutStartTimer(timer);

				events = (int3*)malloc(sizeof(int3)*(numEvents*numberAtr));
				values = (float4*)malloc(sizeof(float4)*(numEvents*numberAtr)*2);
				sizeEvents = (int2*)malloc(sizeof(int2)*numEvents);

				predicates = (int3*)malloc(sizeof(int3)*numberPredicates);
				valuesPrd = (float4*)malloc(sizeof(float4)*(numberPredicates*2));

				initEvents(numberAtr, numEvents);
				initPred();

				cutilSafeCall(cudaMalloc((void**) &predicates_d, sizeof(int3)*numberPredicates));
				cutilSafeCall(cudaMalloc((void**) &valuesPrd_d, sizeof(float4)*(numberPredicates*2)));
	
				if(constant){
					CUDA_SAFE_CALL(cudaMemcpyToSymbol(events_d, events,sizeof(int3)*(numEvents*numberAtr)));
					CUDA_SAFE_CALL(cudaMemcpyToSymbol(sizeEvents_d, sizeEvents, sizeof(int2)*numEvents));
					CUDA_SAFE_CALL(cudaMemcpyToSymbol(valuesEvents_d, values, sizeof(float4)*((numEvents*numberAtr)*2)));
				}
				else
				{
					cutilSafeCall(cudaMalloc((void**) &events_dev, sizeof(int3)*(numEvents*numberAtr)));
					cutilSafeCall(cudaMalloc((void**) &sizeEvents_dev, sizeof(int2)*numEvents));
					cutilSafeCall(cudaMalloc((void**) &values_d, sizeof(float4)*((numEvents*numberAtr)*2)));

					cutilSafeCall(cudaMemcpy(values_d, values,sizeof(float4)*((numEvents*numberAtr)*2),cudaMemcpyHostToDevice));
					cutilSafeCall(cudaMemcpy(events_dev, events,sizeof(int3)*(numEvents*numberAtr),cudaMemcpyHostToDevice));
					cutilSafeCall(cudaMemcpy(sizeEvents_dev, sizeEvents, sizeof(int2)*numEvents, cudaMemcpyHostToDevice));
					
				}
				
				cutilSafeCall(cudaMemcpy(predicates_d, predicates,sizeof(int3)*numberPredicates,cudaMemcpyHostToDevice));		
				cutilSafeCall(cudaMemcpy(valuesPrd_d, valuesPrd,sizeof(float4)*(numberPredicates*2),cudaMemcpyHostToDevice));				

			for(int i = 3; i < 4; i++)
			{
				f << "number blocks = " << i << endl;

				for(int j = 5; j < 6; j++)
				{
					numblocks = i;
					numthreads = pow(2.0, j);
				//	cout << numblocks << " " << numthreads << " part = " << usingshared << " pred = " << numberPredicates << " e = " << numEvents << endl;

					int totthreads = numblocks*numthreads;
					int* tt_d;
					double totaltime = 0.0;

					cutilSafeCall(cudaMalloc((void**) &tt_d, sizeof(int)));
					cutilSafeCall(cudaMemcpy(tt_d, &totthreads,sizeof(int),cudaMemcpyHostToDevice));

					for(int k = 0; k < 1; k++)
					{
						/*if(k%10 == 0)
							cout << k << endl;*/

						if(!usingshared || numEvents == 32)
							cutilSafeCall(cudaMalloc((void**) &bitmap, sizeof(unsigned int)*(numberPredicates+1)));
						else
							cutilSafeCall(cudaMalloc((void**) &bitmap, sizeof(unsigned int)*((int)ceil((numberPredicates+1)/32.0))));

						if(numEvents == 32){
							if(!usingshared && !constant)
							{
								cutStartTimer(timer);
								kernelBit32Events<<<numblocks, numthreads>>>(tt_d, predicates_d, numPred_d, valuesPrd_d, events_dev, sizeEvents_dev, values_d, bitmap);
								cudaThreadSynchronize();
								cutStopTimer(timer);
							}else if(usingshared && !constant){
								cutStartTimer(timer);
								kernelBit32EventsShared<<<numblocks, numthreads, sizeof(float4)*((numEvents*numberAtr)*2/* 10 KB = 16*640/1024 */)>>>(tt_d, predicates_d, numPred_d, valuesPrd_d, events_dev, sizeEvents_dev, values_d, bitmap);
								cudaThreadSynchronize();
								cutStopTimer(timer);
							}else
							{
								cutStartTimer(timer);
								kernelBit32EventsConstant<<<numblocks, numthreads>>>(tt_d, predicates_d, numPred_d, valuesPrd_d, bitmap);
								cudaThreadSynchronize();
								cutStopTimer(timer);
							}
						}
						else
						{
							if(!usingshared && !constant){
								cutStartTimer(timer);
								kernelBit1Event<<<numblocks, numthreads>>>(tt_d, predicates_d, numPred_d, valuesPrd_d, events_dev, sizeEvents_dev, values_d, bitmap);
								cudaThreadSynchronize();
								cutStopTimer(timer);
							}else if(usingshared && !constant){
								cutStartTimer(timer);
								kernelBit1EventShared<<<numblocks, numthreads, sizeof(unsigned int)*((int)ceil((numberPredicates+1)/32.0))>>>(tt_d, predicates_d, numPred_d, valuesPrd_d, events_dev, sizeEvents_dev, values_d, bitmap);
								cudaThreadSynchronize();
								cutStopTimer(timer);
							}else
							{
								cutStartTimer(timer);
								kernelBit1EventsConstant<<<numblocks, numthreads>>>(tt_d, predicates_d, numPred_d, valuesPrd_d, bitmap);
								cudaThreadSynchronize();
								cutStopTimer(timer);
							}
						}
						totaltime+=cutGetTimerValue(timer)/(double)CLOCKS_PER_SEC;
						cutResetTimer(timer);
						//cutilSafeCall(cudaFree(bitmap));
					}

					/*cudaError_t err = cudaGetLastError();
					printf("Cuda error on bitmap: %s.\n", cudaGetErrorString(err) );	*/
					//cout << totaltime << endl;
					f << totaltime << endl;
					cutilSafeCall(cudaFree(tt_d));
				}
			}
				
				free(events);
				free(values);
				free(sizeEvents);
				free(predicates);
				free(valuesPrd);
				
				if(!constant){
					cutilSafeCall(cudaFree(events_dev));
					cutilSafeCall(cudaFree(sizeEvents_dev));
					cutilSafeCall(cudaFree(values_d));
				}
				cutilSafeCall(cudaFree(predicates_d));
				cutilSafeCall(cudaFree(valuesPrd_d));

			//cutilSafeCall(cudaFree(bitmap));
	//cutilSafeCall(cudaFree(tt_d));
}

//############################ SUBSCRIPTIONS ######################################
void initSubs(int ord)
{
	clusters = randomClusters(numberSubs, numberClusters, numberPredicates, ord, percPerC);
	int* sizes;
	sizes = (int*)malloc(sizeof(int)*(numberClusters+1));

	for(int i = 1; i <= numberClusters; i++)
		sizes[i] = (int)pow(2.0, i);//2,4,8,16,32,64,128

	cutilSafeCall(cudaMalloc((void**) &d_sizes, sizeof(int)*(numberClusters+1)));
	cutilSafeCall(cudaMemcpy(d_sizes, sizes,sizeof(int)*(numberClusters+1),cudaMemcpyHostToDevice));

	free(sizes);

	cutilSafeCall(cudaMalloc((void**) &d_clusters, sizeof(int)*(numberClusters+1)));

	cutilSafeCall(cudaMalloc((void**) &d_clusters[0], sizeof(int)*(numberClusters+1)));
	cutilSafeCall(cudaMemcpy(d_clusters[0], clusters[0],sizeof(int)*(numberClusters+1),cudaMemcpyHostToDevice));

	for(int i = 1; i <= numberClusters; i++){
		cutilSafeCall(cudaMalloc((void**) &d_clusters[i], sizeof(int)*clusters[0][i]));
		cutilSafeCall(cudaMemcpy(d_clusters[i], clusters[i],sizeof(int)*clusters[0][i],cudaMemcpyHostToDevice));
	}
}

void removeSubs(int begin, int numSubsToDel)
{
	srand ( 10000 );
	/*ofstrea f;
	f.open ("newSubs.txt",ios::ate);*/

	int2** subsDel;
	int2** subs = (int2**)malloc(sizeof(int2)*(numberClusters+1));

	for(int i = 1; i <= numberClusters; i++)
	{	
		subs[i] = (int2*)malloc(sizeof(int2)*(numSubsToDel+1));
		subs[i][0].x = 1;
	}
		
	int val;

	for(int i = begin; i < numSubsToDel; i++)
	{
		val = i;
		/*while(deletedSubs.count(val) != 0)
			val = rand();

		if(i%2 == 0)
			val = numberSubs - val;*/

		int2 delsub;

		delsub.x = subscriptions[val].y;
		delsub.y = subscriptions[val].z;

		subs[subscriptions[val].y][subs[subscriptions[val].y][0].x] = delsub;
		subs[subscriptions[val].y][0].x++;
		
		deletedSubs[val] = subscriptions[val];

		//numSubsDel[subscriptions[val].y]++;

		subscriptions.erase(val);
	}

	cutilSafeCall(cudaMalloc((void**) &subsDel, sizeof(int2)*(numberClusters+1)));
	cutilSafeCall(cudaMemcpy(subsDel, subs, sizeof(int2)*(numberClusters+1),cudaMemcpyHostToDevice));

	for(int i = 1; i <= numberClusters; i++)
	{
		cutilSafeCall(cudaMalloc((void**) &subsDel[i], sizeof(int2)*subs[i][0].x));
		cutilSafeCall(cudaMemcpy(subsDel[i], subs[i], sizeof(int2)*subs[i][0].x,cudaMemcpyHostToDevice));
	}

	for(int i = 1; i <= numberClusters; i++)
		if(subs[i][0].x > 1)
			kernelRemoveSubs<<<1,1>>>(subsDel[i], d_clusters[i], numPred_d);


	for(int i = 1; i <= numberClusters; i++){
		free(subs[i]);	
		cutilSafeCall(cudaFree(subsDel[i]));
	}

	cutilSafeCall(cudaFree(subsDel));
	free(subs);	

	cudaThreadSynchronize();

//	f.close();

}

void insertNewSubs(int** newSubs, int* sizeSubs, int** positions)
{
	//int newSubs[] = {numberSubs, 3426,9113, numberSubs+1, 3426, 9113};
	int**  d_newSubs;
	cutilSafeCall(cudaMalloc((void**) &d_newSubs, sizeof(int)*(numberClusters+1)));

	int**  d_pos;
	cutilSafeCall(cudaMalloc((void**) &d_pos, sizeof(int)*(numberClusters+1)));
	cutilSafeCall(cudaMemcpy(d_pos, positions, sizeof(int)*(numberClusters+1),cudaMemcpyHostToDevice));
	
	int** d_numsubs;
	cutilSafeCall(cudaMalloc((void**) &d_numsubs, sizeof(int)*(numberClusters+1)));

	for(int i = 1; i <= numberClusters; i++)
	{
		cutilSafeCall(cudaMalloc((void**) &d_pos[i], sizeof(int)*sizeSubs[i]));
		cutilSafeCall(cudaMalloc((void**) &d_newSubs[i], sizeof(int)*(sizeSubs[i]*((int)pow(2.0,i)+1))));
		cutilSafeCall(cudaMemcpy(d_newSubs[i], newSubs[i],sizeof(int)*(sizeSubs[i]*((int)pow(2.0,i)+1)),cudaMemcpyHostToDevice));

		cutilSafeCall(cudaMemcpy(d_pos[i], positions[i], sizeof(int)*sizeSubs[i], cudaMemcpyHostToDevice));
		cutilSafeCall(cudaMalloc((void**) &d_numsubs[i], sizeof(int)*(numberClusters+1)));
		cutilSafeCall(cudaMemcpy(d_numsubs[i], &sizeSubs[i], sizeof(int), cudaMemcpyHostToDevice));
	}

	int* size_d;
	int poww;

	for(int i = 1; i <= numberClusters; i++)
	{
		if(sizeSubs[i] > 0){
			poww = (int)pow(2.0,i);
			cutilSafeCall(cudaMalloc((void**) &size_d, sizeof(int)));
			cutilSafeCall(cudaMemcpy(size_d, &poww, sizeof(int), cudaMemcpyHostToDevice));
			kernelInsertSubs<<<1,1>>>(d_newSubs[i], d_clusters[i], d_pos[i], d_numsubs[i], size_d);	
		}
	}

	/*cudaError_t err = cudaGetLastError();
	printf("Cuda error on bitmap: %s.\n", cudaGetErrorString(err) );*/
	
	cutilSafeCall(cudaFree(size_d));	
	cudaThreadSynchronize();
}



void obtainDisjSubs()
{
	specialCluster = (unsigned int*)malloc(sizeof(int)*(numberSubsDisj*4));

	specialCluster[0] = 0;
	specialCluster[1] = 0;
	specialCluster[2] = 1;
	specialCluster[3] = 2;

	specialCluster[4] = 1;
	specialCluster[5] = 2;
	specialCluster[6] = 5;
	specialCluster[7] = 7;

	specialCluster[8] = 2;
	specialCluster[9] = 2;
	specialCluster[10] = 11;
	specialCluster[11] = 8833;

}


void obtainMatchingDisjSubs(unsigned int* resultSubs)
{
	obtainDisjSubs();
	unsigned int* disjDevice;
	unsigned int* resultSubs_d;
	unsigned int* resultF;
	int* numSubsDisj;
	int* size_d;
	int s = 3;
	int tt = 2*32;

	cutilSafeCall(cudaMalloc((void**)&resultF, sizeof(unsigned int)*numberSubsDisj));
	cutilSafeCall(cudaMalloc((void**)&resultSubs_d, sizeof(unsigned int)*(numberSubs)));
	cutilSafeCall(cudaMalloc((void**)&disjDevice, sizeof(unsigned int)*numberSubsDisj*4));
	cutilSafeCall(cudaMalloc((void**)&d_totalthreads, sizeof(int)));
	cutilSafeCall(cudaMalloc((void**)&numSubsDisj, sizeof(int)));
	cutilSafeCall(cudaMalloc((void**)&size_d, sizeof(int)));

	cutilSafeCall(cudaMemcpy(disjDevice, specialCluster,sizeof(unsigned int)*numberSubsDisj*4,cudaMemcpyHostToDevice));
	cutilSafeCall(cudaMemcpy(resultSubs_d, resultSubs,sizeof(unsigned int)*numberSubs,cudaMemcpyHostToDevice));
	cutilSafeCall(cudaMemcpy(d_totalthreads, &tt,sizeof(int),cudaMemcpyHostToDevice));
	cutilSafeCall(cudaMemcpy(numSubsDisj, &numberSubsDisj,sizeof(int),cudaMemcpyHostToDevice));
	cutilSafeCall(cudaMemcpy(size_d, &s,sizeof(int),cudaMemcpyHostToDevice));

	kernelMatch32EventsDisjunctions<<<2,32>>>(d_totalthreads, resultSubs_d, disjDevice, numSubsDisj, size_d, resultF);

	unsigned int* r = (unsigned int*)malloc(sizeof(int)*numberSubsDisj);
	cutilSafeCall(cudaMemcpy(r, resultF,sizeof(int)*numberSubsDisj,cudaMemcpyDeviceToHost));

	int match = 0;

	for(int j = 0; j < 32; j++)
	{
			//fprintf(f2, "Para o evento %d foram aceites:\n", j);
			for(int i = 0; i < numberSubsDisj; i++)
			{
					unsigned int currSub = r[i];
					//cout << currSub << endl;
					if(ACTIVE_BIT(currSub, j)>0){
						cout << i << " and event = " << j << endl;
						//fprintf(f2,"%d\n",i);
						match++;
					}
			}
	}
				//cout << match << endl;
	cout << "Foram aceites " << match  << "subscricoes\n";
}
void obtainMatchingSubs(int numEvents, int perc)
{
	ofstream f;
	char str[20];
	itoa(numberPredicates,str,10);
	f.open (str,ios::app|std::ios::out);

	char str2[2];
	itoa(numEvents, str2, 10);
	FILE *f2 = fopen(str2, "a");
	
	unsigned int timer;
	cutilCheckError( cutCreateTimer(&timer));
	
	unsigned int* d_Result;

	printf("Phase 2: Start running Kernel...\n");

	cout << "######################### PREDICADOS = " << numberPredicates << " ###############################\n";
	f << "######################### NUMBER EVENTS = " << numEvents << " ###############################\n";

	if(!perc)
		f<< "---------------------------- USING BIT = RANDOM ---------------------------------" << endl;
	else
		f<< "---------------------------- USING BIT = PERC ---------------------------------" << endl;
	
	int tt;
	cudaThreadSynchronize();

	for(int mul = mulInit; mul <= mulFinal; mul +=10)
	{
			double percent = mul/100.0;
			if(perc)
			{
				f << "PERC = " << percent << endl;
				cout << "Starting with " << percent << endl;

				if(numEvents == 32)
					generateRandomBitmap32(percent);
				else
					generateRandomBitmap(percent);
			}
			else{
			//	if(nb*totalthreads <= 224)
					generateBitmap(6, 32, numEvents);
			//	else
			//		generateBitmap(nb, totalthreads, numEvents);
			}
		for(int nb = nbInit; nb < nbFinal; nb++)
		{
			num_blocks = nb;
			f << nb << " blocos" << endl;

			for(int j = ntInit; j < ntFinal; j++)
			{
				double totaltime = 0.0;		

				totalthreads = (int)pow(2.0,j);

				cudaThreadSynchronize();

				tt = totalthreads * num_blocks;
				cutilSafeCall(cudaMemcpy(d_totalthreads, &tt,sizeof(int),cudaMemcpyHostToDevice));

				for(int k = 0; k < numIter; k++)
				{
					cutilSafeCall(cudaMalloc((void**) &d_Result, sizeof(unsigned int)*numberSubs));
					unsigned int* resultFinal = (unsigned int*)malloc(sizeof(unsigned int)*numberSubs);//copiar resultados obtidos no GPU para uma estrutura no Host

					cudaThreadSynchronize();
				
					if(numEvents == 32){
							cutStartTimer(timer);

							for(int i = 1; i <= numberClusters; i++)
								kernelMatch32Events<<<num_blocks, totalthreads>>>(d_totalthreads, bitmap, d_clusters[i], &d_clusters[0][i], &d_sizes[i], d_Result);//assincrona (devolve sem completar a execução)

							cudaThreadSynchronize();
							cutStopTimer(timer);
					}
					else
					{
						if(!usingshared){
							cutStartTimer(timer);

							for(int i = 1; i <= numberClusters; i++)	
								kernelMatch1Event<<<num_blocks, totalthreads>>>(d_totalthreads, bitmap, d_clusters[i], &d_clusters[0][i], &d_sizes[i], d_Result);//assincrona (devolve sem completar a execução)

							cudaThreadSynchronize();
							cutStopTimer(timer);
						}
						else
						{
							int* sizebit;
							cutilSafeCall(cudaMalloc((void**) &sizebit, sizeof(int)));
							int sb = (int)ceil((numberPredicates+1)/32.0);
						//	cout << sb << endl;
							cutilSafeCall(cudaMemcpy(sizebit, &sb,sizeof(int),cudaMemcpyHostToDevice));
							cutStartTimer(timer);

							for(int i = 1; i <= numberClusters; i++)	
								kernelMatch1EventShared<<<num_blocks, totalthreads,	sizeof(int)*(int)ceil((numberPredicates+1)/32.0)>>>(d_totalthreads, bitmap, d_clusters[i], &d_clusters[0][i], &d_sizes[i], sizebit, d_Result);//assincrona (devolve sem completar a execução)

							cudaThreadSynchronize();
							cutStopTimer(timer);
							cutilSafeCall(cudaFree(sizebit));
						}
					}

					/*cudaError_t err = cudaGetLastError();
					printf("Cuda error on bitmap: %s.\n", cudaGetErrorString(err) );*/

					totaltime+=cutGetTimerValue(timer)/(double)CLOCKS_PER_SEC;
					cutResetTimer(timer);

					if(k%10 == 0)
						cout << nb << " " << totalthreads << " " << k << endl;

					cutilSafeCall(cudaMemcpy(resultFinal, d_Result, sizeof(unsigned int)*numberSubs,cudaMemcpyDeviceToHost));

					int match = 0;
					if(numEvents == 32)
					{
						for(int j = 0; j < 32; j++)
						{
							fprintf(f2, "Para o evento %d foram aceites:\n", j);
							for(int i = 0; i < numberSubs; i++)
							{
								unsigned int currSub = resultFinal[i];
								//cout << currSub << endl;
								if(ACTIVE_BIT(currSub, j)>0){
									fprintf(f2,"%d\n",i);
									//cout << currSub << endl;
									match++;
								}
							}
						}
						//cout << match << endl;
						fprintf(f2,"Foram aceites %d subscricoes\n", match);
					}
					else
					{
						for(int i = 0; i < numberSubs; i++)
						{
							int currSub = resultFinal[i];
							if(currSub > 0)
								match++;
						}
						fprintf(f2,"Foram aceites %d subscricoes\n",match);
						//cout << match << endl;
					}

					cutilSafeCall(cudaFree(d_Result));
					//obtainMatchingDisjSubs(resultFinal);
					free(resultFinal);
				}

				cout << totaltime << endl;
				f << totaltime << endl;
			}
		}
		if(bit)
			cutilSafeCall(cudaFree(bitmap));
	}
				
	free(str);
	free(str2);
	f.close();
	fclose(f2);
}




void firstInit(){

	//############ TAMANHO DA EXPERIENCIA ##################
	numberClusters = 7;//indicar o numero de clusters
	//numberPredicates = 1000;//numero de predicados
	numberSubs = 500000;//numero de subscricoes
	numIter = 100;//numero de iterações
	numberSubsDisj = 3;

	//#### PERCENTAGEM DE DIVISAO DOS PREDICADOS PELAS SUBSCRICOES #####
	percPerC = 0.5;
	//############ NUM THREADS DE BASE 2 ##################
	ntInit = 8;//numero de threads iniciais 2^ntInit
	ntFinal = 9;//numero de threads finais 2^ntFinal (exclusive)
	//constant = 1;

	//############ NUM BLOCOS DE BASE ##################
	nbInit = 5;//numero de blocos iniciais
	nbFinal = 6;//numero de blocos finais (exclusive)

	//############ COM OU SEM BITMAP (USANDO OU NÃO REDUÇÃO A NÍVEL DE BITS) ################
	bit = 1;//calcular o bit
	//usingshared = 1;//usar bitmap com bits ou nao

	//##### QUAL KERNEL DE MATCHING #####
	kernelnum = 1;//0 para usar o de divisao de subs pelas threads /1 -  para a versao warp /2 - versao "pior"

	cutilSafeCall(cudaMalloc((void**) &d_totalthreads, sizeof(int)));
	cutilSafeCall(cudaMalloc((void**) &numPred_d, sizeof(int)));
	cutilSafeCall(cudaMemcpy(numPred_d, &numberPredicates,sizeof(int),cudaMemcpyHostToDevice));	
	//numSubsDel = (int*)malloc(sizeof(int)*(numberClusters+1));
	/*for(int i = 0; i < (numberClusters+1); i++)
		numSubsDel[i] = 0;*/
}


void clearMemory()
{
//	free(numSubsDel);
	cutilSafeCall(cudaFree(d_sizes));

	for(int i = 1; i <= numberClusters; i++)
		cutilSafeCall(cudaFree(d_clusters[i]));

	cutilSafeCall(cudaFree(d_clusters));

	for(int i = 0; i <= numberClusters; i++)
		free(clusters[i]);

	free(clusters);

//	cutilSafeCall(cudaFree(bitmap));
	cutilSafeCall(cudaFree(numPred_d));
}

int
main(int argc, char** argv)
{
	printf("Welcome to a matching publish/subscribe system! Made in CUDA!\n");

//	cout << "Running program with " << num_blocks << " blocks and " << num_threadsperblock << " threads per block = " << totalthreads << " total threads." << endl;
	
	unsigned int timer;
	cutilCheckError( cutCreateTimer(&timer));

	cutStartTimer(timer);
	usingshared = 0;
	constant = 0;
	numberPredicates = 10000;
	firstInit();
	initSubs(0);
	generateBitmap(0, 0, 1);
//	obtainMatchingSubs(1,0);
	clearMemory();
	/*cout << "PHASE 1" << endl;
	for(int i = 0; i < 2; i++){

		constant = i;
		cout << "CONSTANT: " << i << endl;
		for(int j = 1; j < 100; j*=10)
		{
			numberPredicates = 1000*j;
			cout << "PREDICATES: " << numberPredicates << endl;
			firstInit();
			initSubs(0);
			cout << "BIT" << endl;
			generateBitmap(0, 0, 1);
			cout << "OTHER BIT" << endl;
			generateBitmapVar(1);
			clearMemory();
		}
			numberPredicates = 30000;
			cout << "PREDICATES: " << numberPredicates << endl;
			firstInit();
			initSubs(0);
			cout << "BIT" << endl;
			generateBitmap(0, 0, 1);
			cout << "OTHER BIT" << endl;
			generateBitmapVar(1);
			clearMemory();

	}

	cout << "PHASE 2" << endl;
	for(int i = 0; i < 2; i++){

		constant = i;
		cout << "CONSTANT: " << i << endl;

		for(int j = 1; j < 100; j*=10)
		{
			numberPredicates = 1000*j;
			cout << "PREDICATES: " << numberPredicates << endl;
			firstInit();
			initSubs(0);
			cout << "BIT" << endl;
			generateBitmap(0, 0, 32);
			cout << "OTHER BIT" << endl;
			generateBitmapVar(32);
			clearMemory();
		}
			numberPredicates = 30000;
			cout << "PREDICATES: " << numberPredicates << endl;
			firstInit();
			initSubs(0);
			cout << "BIT" << endl;
			generateBitmap(0, 0, 32);
			cout << "OTHER BIT" << endl;
			generateBitmapVar(32);
			clearMemory();
	}
*/
	cutStopTimer(timer);
	printf("total time of execution = %f \n",cutGetTimerValue(timer)/(double)CLOCKS_PER_SEC);
	cudaThreadExit();
	//obtainMatchingSubs(32,0);
		//obtainMatchingSubs(32,1);
	//cutilExit(argc, argv);
	exit(EXIT_SUCCESS);

	return 0;
}

/*CUDA_SAFE_CALL( cudaMemcpyToSymbol(c_view_planes_direction, h_view_planes_direction, sizeof(bool)*3) );
-----declaração do array----
_constant_ bool c_view_planes_direction[3];*/
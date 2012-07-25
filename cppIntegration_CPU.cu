//// includes, system
//#include <stdlib.h>
//#include <stdio.h>
//#include <string.h>
//#include <math.h>
//
//#include <map>
//#include <string>
//#include <iostream>
//#include <fstream>
//// includes, project
//#include <cutil_inline.h>
//
//#include "driver_types.h"
//#include "cuda_runtime.h"
//
//using namespace std;
//
//#define ACTIVE_BIT(i, j) (i & (1 << j))
//#define EDIT_BIT(i, j) (i | (1 << j))
//
//using namespace std;
//
////Host
//unsigned int** clusters;
//int numberClusters;
//int numberSubs;
//int numberPredicates;
//int totalthreads;
//int num_blocks;
//int num_threadsperblock;
//int ntInit;
//int ntFinal;
//int nbInit;
//int nbFinal;
//int mulInit;
//int mulFinal;
//int bit;
//int numIter;
//int samePred;
//int same;
//int inverse;
//int usingshared;
//int kernelnum;
//int ran;
//int constant;
//double percPerC;
//int removedSubs;
//unsigned int* bitmap_CPU;
//
//int3* events;
//float4* values;
//int2* sizeEvents;
//int3* predicates;
//float4* valuesPrd;
//
//map<int, int3> subscriptions;
//map<int, int3> deletedSubs;
//
//unsigned int randomNumber(int maxV)
//{
//	return rand() % maxV;
//}
//
//
//unsigned int** randomClusters(int numSubs, int numClus, int numPred, int qs, double perc)
//{
//	unsigned int** r = (unsigned int**)malloc(sizeof(int)*(numClus+1));
//	int i = 1;
//	int init = numSubs;
//	int tSubs = 0;
//	int idSubs = 0;
//	int init2;
//
//	srand ( 10000 );
//
//	r[0] = (unsigned int*)malloc(sizeof(int)*(numClus+1));
//
//	while(i != (numClus+1))
//	{
//		if(i!=numClus)
//			init2 = (((int)(init * perc)) >> 5) << 5;
//		else
//			init2 = (int)(numSubs-tSubs);
//
//		tSubs+=init2;
//
//		int poww = (int)pow(2.0, i);
//
//		cout << init2 << endl;
//
//		r[0][i] = init2*(poww+1);
//		int currS = 0;
//
//		r[i] = (unsigned int*)malloc(sizeof(int)*init2*(poww+1));
//
//		while(currS != init2)
//		{
//			unsigned int* preds = (unsigned int*)malloc(sizeof(int)*(poww+1));
//			preds[0] = idSubs;
//
//			for(int j = 1; j <= poww; j++){
//
//				int found = 1;
//				unsigned int num;
//
//				while(found == 1)
//				{
//					num = randomNumber(numPred);
//					
//					for(int k = 1; k < j; k++)
//					{
//						if(preds[k] == num)
//						{
//							found = 2;
//							break;
//						}
//					}
//
//					if(found == 2)
//						found = 1;
//					else
//						found = 0;
//				}
//				preds[j] = num;
//			}
//
//			subscriptions[idSubs].x = idSubs; 
//			subscriptions[idSubs].y = i;
//			subscriptions[idSubs].z = currS*(poww+1);
//
//			memcpy(&r[i][currS*(poww+1)],preds,(poww+1)*sizeof(int));
//
//			currS++;
//
//			idSubs++;
//
//			free(preds);
//
//		}
//
//		init-=init2;
//
//		i++;
//	}
//	return r;
//
//}
//
//unsigned int* give32Bit()
//{
//	unsigned int* ret = (unsigned int*)malloc(sizeof(int)*numberPredicates);
//	srand ( 4294967295 );//max 32767
//
//	for(int i = 0; i < numberPredicates; i++)
//		ret[i] = (unsigned int)(rand() << 17 | rand() << 2 | rand() >> 13);
//
//	return ret;
//}
//
//unsigned int* randomPos(double perc)
//{
//	srand ( 10000 );
//
//	int size2 = static_cast<int>(numberPredicates * perc);
//
//	cout << size2 << endl;
//
//	unsigned int* randomPositions = (unsigned int*)malloc(sizeof(int)*size2);
//
//	unsigned int num;
//
//	randomPositions[0] = randomNumber(numberPredicates);
//	
//		for(int i = 1; i < size2 ; i++)
//		{
//				int found = 1;
//
//				while(found == 1)
//				{
//					num = randomNumber(numberPredicates);
//					for(int j = 0; j < i; j++)
//					{
//						if(randomPositions[j] == num)
//						{
//							found = 2;
//							break;
//						}
//					}
//					if(found == 2)
//						found = 1;
//					else
//						found = 0;
//				}
//
//				randomPositions[i] = num;
//		}
//
//	return randomPositions;
//
//}
//
//unsigned int** randomPercAccepted(double perc)
//{
//	unsigned int** ret = (unsigned int**)malloc(sizeof(unsigned int)*32);
//
//	for(int i = 0; i < 32; i++)
//		ret[i] = randomPos(perc);
//
//	return ret;
//}
////############################ PREDICATES ######################################
//void initPred()
//{
//	srand ( 10000 );
//
//	int pos2 = 0;
//	for(int j = 0; j < numberPredicates; j++)
//	{
//		pos2 = j*2;
//		predicates[j].x = j;
//		predicates[j].y = (int)randomNumber(7);
//		predicates[j].z = pos2;
//	}
//	for(int i = 0; i < (numberPredicates*2); i+=2)
//	{
//		valuesPrd[i].x = (float)i;
//		valuesPrd[i].y = (float)i;
//		valuesPrd[i].z = (float)i;
//		valuesPrd[i].w = (float)i;
//
//		valuesPrd[i+1].x = (float)i;
//		valuesPrd[i+1].y = (float)i;
//		valuesPrd[i+1].z = (float)i;
//		valuesPrd[i+1].w = (float)i;
//	}
//}
//
//
////############################ EVENTS ######################################
//void initEvents(int numberAtr, int numEvents)
//{
//	srand ( 10000 );
//
//	int position = 0;
//	for(int i = 0; i < numEvents; i++)
//	{
//		sizeEvents[i].x = position;
//
//		int pos = 0;
//		for(int j = 0; j < numberAtr; j++)
//		{
//			events[i*numberAtr + j].x = (int)randomNumber(numberPredicates); 
//			events[i*numberAtr + j].y = i*(numberAtr*2)+pos;
//			events[i*numberAtr + j].z = i*(numberAtr*2)+(pos+1);
//			pos+=2;
//		}
//		position += numberAtr;
//		sizeEvents[i].y = position;
//	}
//
//	for(int j = 0; j < (numEvents*numberAtr)*2; j+=2)
//	{
//		int ran = randomNumber(numberPredicates*2);
//		values[j].x = (float)ran;
//		values[j].y = (float)ran;
//		values[j].z = (float)ran;
//		values[j].w = (float)ran;
//
//		values[j+1].x = (float)ran;
//		values[j+1].y = (float)ran;
//		values[j+1].z = (float)ran;
//		values[j+1].w = (float)ran;
//	}
//
//}
//
////############################ BITMAP ######################################
//void initSubs(int ord)
//{
//	clusters = randomClusters(numberSubs, numberClusters, numberPredicates, ord, percPerC);
//}
//
//void generateRandomBitmapCPU(double perc)
//{
//		unsigned int* r = randomPos(perc);
//		
//		int s = (int)(numberPredicates*perc);
//
//		if(!usingshared){
//			bitmap_CPU = (unsigned int*)malloc(sizeof(unsigned int)*(numberPredicates+1));
//			for(int i = 0; i < numberPredicates; i++)
//				bitmap_CPU[i] = 0;
//		}
//		else
//		{
//			int bitsize = (int)ceil((numberPredicates+1)/32.0);
//			bitmap_CPU = (unsigned int*)malloc(sizeof(unsigned int)*(bitsize));
//			for(int i = 0; i < bitsize; i++)
//				bitmap_CPU[i] = 0;
//		}
//
//		int pos;
//		if(usingshared)
//		{
//			for(int i = 0; i < s; i++)
//			{
//				pos = r[i]/32;
//				bitmap_CPU[pos] = EDIT_BIT(bitmap_CPU[pos], r[i]%32);//bitmap[pos] = 1;
//			}
//		}
//		else
//		{
//			for(int i = 0; i < s; i++)
//				bitmap_CPU[r[i]] = 1;
//		}
//
//		free(r);
//}
//
//void generateRandomBitmap32CPU(double perc)
//{
//		unsigned int** r = randomPercAccepted(perc);
//		
//		int s = (int)(numberPredicates*perc);
//
//		bitmap_CPU = (unsigned int*)malloc(sizeof(unsigned int)*(numberPredicates+1));
//	
//		for(int i = 0; i < numberPredicates; i++)
//			bitmap_CPU[i] = 0;
//
//		int ran;
//
//		for(int i = 0; i < 32; i++)
//		{
//			for(int j = 0; j < s; j++)
//			{
//				ran = r[i][j];
//				bitmap_CPU[ran] = EDIT_BIT(bitmap_CPU[ran], i);//bitmap[pos] = 1;
//			}
//		}
//
//		for(int i = 0; i < 32; i++)
//			free(r[i]);
//
//		free(r);
//}
//
//
//
//
//void cpuKernelgenerateBit1()
//{
//	int oper;
//	int atr;
//	float4 val;
//	float4 val2;
//
//	bitmap_CPU[numberPredicates] = 0;
//
//	for(int i = 0; i < numberPredicates; i++)
//	{
//		int3 currPred = predicates[i];
//		val = valuesPrd[currPred.z];
//		val2 = valuesPrd[currPred.z+1];
//		oper = currPred.y;
//		atr = currPred.x;
//
//		int size = sizeEvents[0].y;
//
//		for(int k = 0; k < size; k++)
//		{
//			int3 elemEvent = events[k];
//
//			if(elemEvent.x == atr)
//			{
//				for(int l = elemEvent.y; l < elemEvent.z; l+=2)
//				{
//					float4 currv = values[l];
//					float4 currv2 = values[l+1];
//
//					switch (oper)
//					{
//					case 0 : 
//						if((val.x == currv.x && val.y == currv.y && val.z == currv.z && val.w == currv.w) && (val2.x == currv2.x && val2.y == currv2.y && val2.z == currv2.z && val2.w == currv2.w))
//							bitmap_CPU[i] = 1;
//						break;
//					case 1 :
//						if(currv.x > val.x)
//							bitmap_CPU[i] = 1;
//						break;
//					case 2 :
//						if(currv.x < val.x)
//							bitmap_CPU[i] = 1;
//						break;
//					case 3 :
//						if((val.x != currv.x || val.y != currv.y || val.z != currv.z || val.w != currv.w) || (val2.x != currv2.x || val2.y != currv2.y || val2.z != currv2.z || val2.w != currv2.w))
//							bitmap_CPU[i] = 1;
//						break;
//					case 4 :
//						if(currv.x <= val.x)
//							bitmap_CPU[i] = 1;
//						break;
//					case 5 :
//						if(currv.x >= val.x)
//							bitmap_CPU[i] = 1;
//						break;
//					case 6 :
//						if((val.x == currv.x && val.y == currv.y && val.z == currv.z && val.w == currv.w) && (val2.x == currv2.x && val2.y == currv2.y && val2.z == currv2.z && val2.w == currv2.w))
//							bitmap_CPU[i] = 1;
//						break;
//					}
//				}
//			}
//		}
//	}
//}
//void cpuKernelgenerateBit32()
//{
//		int oper;
//		int atr;
//		float4 val;
//		float4 val2;
//
//		bitmap_CPU[numberPredicates] = 0;
//
//		for(int i = 0; i < numberPredicates;i++)
//		{
//			int3 currPred = predicates[i];
//			val = valuesPrd[currPred.z];
//			val2 = valuesPrd[currPred.z+1];
//			oper = currPred.y;
//			atr = currPred.x;
//
//			for(unsigned int j = 0; j < 32; j++)
//			{
//				int2 Event = sizeEvents[j];
//
//				for(int k = Event.x; k < Event.y; k++)
//				{
//					int3 elemEvent = events[k];
//
//					if(elemEvent.x == atr)
//					{
//						for(int l = elemEvent.y; l < elemEvent.z; l+=2)
//						{
//							float4 currv = values[l];
//							float4 currv2 = values[l+1];
//
//							switch (oper)
//							{
//								case 0 : 
//									if((val.x == currv.x && val.y == currv.y && val.z == currv.z && val.w == currv.w) && (val2.x == currv2.x && val2.y == currv2.y && val2.z == currv2.z && val2.w == currv2.w))
//											bitmap_CPU[i]=EDIT_BIT(bitmap_CPU[i],j);
//									break;
//								case 1 :
//									if(currv.x > val.x)
//											bitmap_CPU[i]=EDIT_BIT(bitmap_CPU[i],j);
//									break;
//								case 2 :
//									if(currv.x < val.x)
//											bitmap_CPU[i]=EDIT_BIT(bitmap_CPU[i],j);
//									break;
//								case 3 :
//									if((val.x != currv.x || val.y != currv.y || val.z != currv.z || val.w != currv.w) || (val2.x != currv2.x || val2.y != currv2.y || val2.z != currv2.z || val2.w != currv2.w))
//										bitmap_CPU[i]=EDIT_BIT(bitmap_CPU[i],j);
//									break;
//								case 4 :
//									if(currv.x <= val.x)
//										bitmap_CPU[i]=EDIT_BIT(bitmap_CPU[i],j);
//										break;
//								case 5 :
//									if(currv.x >= val.x)
//										bitmap_CPU[i]=EDIT_BIT(bitmap_CPU[i],j);
//										break;
//								case 6 :
//									if((val.x == currv.x && val.y == currv.y && val.z == currv.z && val.w == currv.w) && (val2.x == currv2.x && val2.y == currv2.y && val2.z == currv2.z && val2.w == currv2.w))
//										bitmap_CPU[i]=EDIT_BIT(bitmap_CPU[i],j);
//										break;
//							}
//						}
//					}
//				}
//			}
//		}
//}
//void generateBitmapVar(int numEvents)
//{
//	ofstream f;
//	f.open ("Time-Events.txt",ios::app|std::ios::out);
//
//	unsigned int timer;
//	cutilCheckError( cutCreateTimer(&timer));
//	f << "Number of predicates = " << numberPredicates << endl;
//	f << "######################### EVENTS " << numEvents << " #########################" << endl;
//	//f << "######################### PREDICATES " << numberPredicates << " #########################" << endl;
//
//
//			for(int i = 1; i <= 50; i++)
//			{
//				if(i%10==0)
//					cout << i << endl;
//
//				int numberAtr = i;
//				double totaltime = 0.0;
//
//				//cutStartTimer(timer);
//
//				events = (int3*)malloc(sizeof(int3)*(numEvents*numberAtr));
//				values = (float4*)malloc(sizeof(float4)*(numEvents*numberAtr)*2);
//				sizeEvents = (int2*)malloc(sizeof(int2)*numEvents);
//
//				predicates = (int3*)malloc(sizeof(int3)*numberPredicates);
//				valuesPrd = (float4*)malloc(sizeof(float4)*(numberPredicates*2));
//
//				initEvents(numberAtr, numEvents);
//				initPred();
//
//				if(!usingshared || numEvents == 32){
//					bitmap_CPU = (unsigned int*)malloc(sizeof(unsigned int)*(numberPredicates+1));
//					for(int j = 0; j < numberPredicates; j++)
//						bitmap_CPU[j] = 0;
//				}else{
//					bitmap_CPU = (unsigned int*)malloc(sizeof(unsigned int)*((numberPredicates+1)/32.0));
//					for(int j = 0; j < (numberPredicates+1)/32.0; j++)
//						bitmap_CPU[j] = 0;
//				}
//
//				for(int k = 0; k < 100; k++)
//				{
//					if(numEvents == 32){
//							cutStartTimer(timer);
//							cpuKernelgenerateBit32();				
//							cutStopTimer(timer);
//					}else{
//						if(!usingshared){
//							cutStartTimer(timer);
//							cpuKernelgenerateBit1();
//							cutStopTimer(timer);
//						}else{
//							cutStartTimer(timer);
//							//kernelBit1EventShared<<<numblocks, numthreads, sizeof(unsigned int)*((int)ceil((numberPredicates+1)/32.0))>>>(tt_d, predicates_d, numPred_d, valuesPrd_d, events_d, sizeEvents_d, values_d, bitmap);
//							cutStopTimer(timer);
//						}
//					}
//					totaltime+=cutGetTimerValue(timer)/(double)CLOCKS_PER_SEC;
//					cutResetTimer(timer);
//				}
//
//				free(events);
//				free(values);
//				free(sizeEvents);
//				free(predicates);
//				free(valuesPrd);
//			
//				f << "Numero Atr per Event = " << i << endl;
//				f << totaltime << endl;
//
//		}
//}
//
//
//unsigned int* cpuKernels1Event()
//{
//	unsigned int* resultFinal = (unsigned int*)malloc(sizeof(unsigned int)*numberSubs);
//
//		for(int i = 1; i <= numberClusters; i++)
//		{	
//			int size = (int)pow(2.0, i);
//			int end = (clusters[0][i]/(size+1));
//
//			for(int j = 0; j < end; j++)
//			{
//				int match = 0;
//
//				int ID = clusters[i][j*(size+1)];
//
//				for(int k = 1; k <= size; k++)
//				{
//					if(bitmap_CPU[clusters[i][j*(size+1)+k]]==1)
//						match++;
//					else 
//						break;
//				}
//
//				if(match == size)
//					resultFinal[ID] = 1;	
//			}
//		}	
//
//	return resultFinal;
//}
//
//unsigned int* cpuKernels32Event()
//{
//	unsigned int* resultFinal = (unsigned int*)malloc(sizeof(unsigned int)*numberSubs);
//
//	for(int i = 1; i <= numberClusters; i++)
//	{	
//		int size = (int)pow(2.0, i);
//		int end = (clusters[0][i]/(size+1));
//
//		for(int j = 0; j < end; j++)
//		{
//				int ID = clusters[i][j*(size+1)];
//
//				unsigned int match = bitmap_CPU[clusters[i][j*(size+1)+1]];
//
//				for(int k = 2; k <= size && match != 0; k++)
//					match &= bitmap_CPU[clusters[i][j*(size+1)+k]];
//			
//				resultFinal[ID] = match;
//		}
//	}
//
//	return resultFinal;
//}
//
//unsigned int* cpuKernels1EventShared()
//{
//	unsigned int* resultFinal = (unsigned int*)malloc(sizeof(unsigned int)*numberSubs);
//
//	for(int i = 1; i <= numberClusters; i++)
//	{	
//		int size = (int)pow(2.0, i);
//		int end = (clusters[0][i]/(size+1));
//
//		for(int j = 0; j < end; j++)
//		{
//			int match = 0;
//			int ID = clusters[i][j*(size+1)];
//
//			int pos;
//
//			for(int k = 1; k <= size; k++)
//			{
//				pos = clusters[i][j*(size+1)+k];
//
//				if((unsigned int)ACTIVE_BIT(bitmap_CPU[pos/32], pos%32)>0)
//					match++;
//				else 
//					break;
//			}
//
//			if(match == size)
//				resultFinal[ID] = 1;
//		}
//	}	
//
//	return resultFinal;
//}
//
//void obtainMatchingSubsCPU(int numEvents, int perc)
//{
//	ofstream f;
//	char str[20];
//	itoa(numberPredicates,str,10);
//	f.open (str,ios::app|std::ios::out);
//
//	char str2[2];
//	itoa(numEvents, str2, 10);
//	FILE *f2 = fopen(str2, "a");
//	
//	unsigned int timer;
//	cutilCheckError( cutCreateTimer(&timer));
//
//	cout << "######################### PREDICADOS = " << numberPredicates << " ###############################\n";
//	f << "######################### NUMBER EVENTS = " << numEvents << " ###############################\n";
//
//	if(!perc)
//		f<< "---------------------------- USING BIT = RANDOM ---------------------------------" << endl;
//	else
//		f<< "---------------------------- USING BIT = PERC ---------------------------------" << endl;
//	
//	for(int mul = mulInit; mul <= mulFinal; mul +=10)
//	{
//			double percent = mul/100.0;
//
//			if(perc)
//			{
//				f << "PERC = " << percent << endl;
//				cout << "Starting with " << percent << endl;
//
//				if(numEvents == 32)
//					generateRandomBitmap32CPU(percent);
//				else
//					generateRandomBitmapCPU(percent);
//			}
//			else
//				generateBitmapVar(numEvents);
//
//				double totaltime = 0.0;		
//
//				for(int k = 0; k < numIter; k++)
//				{
//					unsigned int* result;
//
//					if(numEvents == 32){
//							cutStartTimer(timer);
//
//							result = cpuKernels32Event();
//
//							cutStopTimer(timer);
//					}
//					else
//					{
//						if(!usingshared){
//							cutStartTimer(timer);
//
//							result = cpuKernels1Event();
//
//							cutStopTimer(timer);
//						}
//						else
//						{
//							cutStartTimer(timer);
//							
//							result = cpuKernels1EventShared();
//
//							cutStopTimer(timer);
//						}
//					}
//
//					totaltime += cutGetTimerValue(timer)/(double)CLOCKS_PER_SEC;
//					cutResetTimer(timer);
//
//					if(k%10 == 0)
//						cout << k << endl;
//
//					int match = 0;
//
//					if(numEvents == 32)
//					{
//						for(int j = 0; j < 32; j++)
//						{
//							//fprintf(f2, "Para o evento %d foram aceites:\n", j);
//							for(int i = 0; i < numberSubs; i++)
//							{
//								unsigned int currSub = result[i];
//								//cout << currSub << endl;
//								if(ACTIVE_BIT(currSub, j)>0){
//									//fprintf(f2,"%d\n",i);
//									match++;
//								}
//							}
//						}
//						//cout << match << endl;
//						fprintf(f2,"Foram aceites %d subscricoes\n", match);
//					}
//					else
//					{
//						for(int i = 0; i < numberSubs; i++)
//						{
//							int currSub = result[i];
//
//							if(currSub > 0)
//								match++;
//						}
//						fprintf(f2,"Foram aceites %d subscricoes\n",match);
//						//cout << match << endl;
//					}
//					free(result);
//				}
//				cout << totaltime << endl;
//				f << totaltime << endl;
//				free(bitmap_CPU);
//	}
//				
//	/*free(str);
//	free(str2);*/
//	f.close();
//	fclose(f2);
//}
//
//
//
//
//void firstInit(){
//
//	//############ TAMANHO DA EXPERIENCIA ##################
//	numberClusters = 7;//indicar o numero de clusters
//	//numberPredicates = 1000;//numero de predicados
//	numberSubs = 500000;//numero de subscricoes
//	numIter = 100;//numero de iterações
//
//	//#### PERCENTAGEM DE DIVISAO DOS PREDICADOS PELAS SUBSCRICOES #####
//	percPerC = 0.5;
//}
//
//
//void clearMemory()
//{
//	for(int i = 0; i <= numberClusters; i++)
//		free(clusters[i]);
//
//	free(clusters);
//}
//
//
//int
//main(int argc, char** argv)
//{
//	printf("Welcome to a matching publish/subscribe system! Made in CUDA!\n");
//
//	unsigned int timer;
//	cutilCheckError( cutCreateTimer(&timer));
//
//	cutStartTimer(timer);
//
//	usingshared = 0;
//	numberPredicates = 1000;
//	firstInit();
//	initSubs(0);
//	generateBitmapVar(1);
//	generateBitmapVar(32);
//	clearMemory();
//
//	usingshared = 0;
//	numberPredicates = 10000;
//	firstInit();
//	initSubs(0);
//	generateBitmapVar(1);
//	generateBitmapVar(32);
//	clearMemory();
//
//	usingshared = 0;
//	numberPredicates = 30000;
//	firstInit();
//	initSubs(0);
//	generateBitmapVar(1);
//	generateBitmapVar(32);
//	clearMemory();
//
//	cutStopTimer(timer);
//	printf("total time of execution = %f \n",cutGetTimerValue(timer)/(double)CLOCKS_PER_SEC);
//	cudaThreadExit();
//	
//	//cutilExit(argc, argv);
//	exit(EXIT_SUCCESS);
//
//	return 0;
//}
//
///*
////for(int i = 1; i < 100; i*=10)
//	//{
//		usingshared = 0;
//
//		numberPredicates = 10000;
//
//		firstInit();
//
//		initSubs(0);
//
//		mulInit = 100;
//		mulFinal = 100;
//
//		obtainMatchingSubsCPU(32,0);
//
//		//obtainMatchingSubsCPU(32,1);
//
//		clearMemory();
//	//}
//
//
//*/
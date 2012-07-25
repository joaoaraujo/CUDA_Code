#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cutil_inline.h>
#include "driver_types.h"
#include "cuda_runtime.h"
#include "clusters_struct.h"
#include <kernel_clusters.cu>


cluster* initClusters(int numberClusters, int sizePerBlock)
{
	cluster* c;
	cutilSafeCall(cudaMalloc((void**) &c, sizeof(cluster)*numberClusters));//= (cluster*)malloc(sizeof(cluster)*4);
	
	int* d_sizePerBlock;
	cutilSafeCall(cudaMalloc((void**) &d_sizePerBlock, sizeof(int)));
	int* d_numPreds;
	cutilSafeCall(cudaMalloc((void**) &d_numPreds, sizeof(int)));
	
	cutilSafeCall(cudaMemcpy(d_sizePerBlock, &sizePerBlock,sizeof(int),cudaMemcpyHostToDevice));


	for(int i = 1; i <= numberClusters; i++){
		int numPreds = (int)pow((double)2,i);
		cutilSafeCall(cudaMemcpy(d_numPreds, &numPreds,sizeof(int),cudaMemcpyHostToDevice));
		kernelInitCluster<<<1,1>>>(&c[i-1], d_numPreds, d_sizePerBlock);
		cutilSafeCall(cudaMalloc((void**) &c[i-1].subs,sizeof(int)*(numPreds*sizePerBlock+sizePerBlock)));
	}

	cudaThreadSynchronize();

	return c;
}

void expandCluster(cluster* c)
{
	int* r;
	cutilSafeCall(cudaMalloc((void**) &r, sizeof(int)*3));
	int* result = (int*)malloc(sizeof(int)*3);
	kernelClusterValues<<<1,1>>>(c,r);
	cutilSafeCall(cudaMemcpy(result, r, sizeof(int)*3, cudaMemcpyDeviceToHost));
	int* newSubs;
	cutilSafeCall(cudaMalloc((void**) &newSubs, sizeof(int)*((result[1]*result[0] + result[1])+(result[2]*result[0] + result[2]))));
	cutilSafeCall(cudaMemcpy(newSubs, c->subs,((result[1]*result[0] + result[1])*sizeof(int)),cudaMemcpyDeviceToDevice));
	//cutilSafeCall(cudaFree(c->subs));
	cutilSafeCall(cudaMalloc((void**) &c->subs, sizeof(int)*((result[1]*result[0] + result[1])+(result[2]*result[0] + result[2]))));
	cutilSafeCall(cudaMemcpy(c->subs, newSubs,((result[1]*result[0] + result[1])*sizeof(int)),cudaMemcpyDeviceToDevice));
	cutilSafeCall(cudaFree(r));
	cutilSafeCall(cudaFree(newSubs));
	free(result);
	kernelExpandCluster<<<1,1>>>(c);
	cudaThreadSynchronize();
}

void moreClusters(int currClus, int moreClus, cluster** clusters, int sizePerBlock)
{
	int* d_numPreds;
	cutilSafeCall(cudaMalloc((void**) &d_numPreds, sizeof(int)));

	int* d_sizePerBlock;
	cutilSafeCall(cudaMalloc((void**) &d_sizePerBlock, sizeof(int)));
	cutilSafeCall(cudaMemcpy(d_sizePerBlock, &sizePerBlock,sizeof(int),cudaMemcpyHostToDevice));


	for(int i = currClus; i < (moreClus+currClus); i++){
		int numPreds = (int)pow((double)2,i+1);
		cutilSafeCall(cudaMemcpy(d_numPreds, &numPreds,sizeof(int),cudaMemcpyHostToDevice));
		cutilSafeCall(cudaMalloc((void**) &clusters[0][i], sizeof(cluster)));
		kernelInitCluster<<<1,1>>>(&clusters[0][i], d_numPreds, d_sizePerBlock);
		cutilSafeCall(cudaMalloc((void**) &clusters[0][i].subs,sizeof(int)*(numPreds*sizePerBlock+sizePerBlock)));
	}

	cudaThreadSynchronize();
}
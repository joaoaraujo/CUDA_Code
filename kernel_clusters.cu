#ifndef _KERNEL_CLUSTERS_H_
#define _KERNEL_CLUSTERS_H_

__global__ void
kernelInitCluster(cluster* in, int* numPredicates, int* subsPerBlock)
{
	in->size = *numPredicates;
	in->numSubs = 0;
	in->maxSubs = *subsPerBlock;
	in->subsPerBlock = *subsPerBlock;
}
__global__ void
kernelExpandCluster(cluster* in)
{
	in->maxSubs = in->maxSubs + in->subsPerBlock;
}
__global__ void
kernelClusterValues(cluster* c, int* result)
{
	result[0] = c->size;
	result[1] = c->maxSubs;
	result[2] = c->subsPerBlock;
}


#endif
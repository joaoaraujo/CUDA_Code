#ifndef _KERNELS_INSREMSUBS_H_
#define _KERNELS_INSREMSUBS_H_

__global__ void
kernelRemoveSubs(int2* subsToDel, int* cluster, int* n_predicates)
{
	for(int i = 1; i < subsToDel[0].x; i++)
		cluster[subsToDel[i].y+1] = n_predicates[0];
}

__global__ void
kernelInsertSubs(int* subsToIns, int* cluster, int* pos, int* numSubs, int* size)
{
	for(int i = 0; i < numSubs[0]; i++)
		for(int j = 1; j <= size[0]; j++)
			cluster[pos[i]+j] = subsToIns[i*(size[0]+1)+j];
}

#endif // #ifndef _KERNELS_INSREMSUBS_H_
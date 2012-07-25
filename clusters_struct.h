#ifndef _CLUSTERS_STRUCT_H
#define _CLUSTERS_STRUCT_H

typedef struct cluster{
	int* subs;
	int size;
	int maxSubs;
	int numSubs;
	int subsPerBlock;
}cluster;

cluster* initClusters(int numberClusters, int sizePerBlock);
void moreClusters(int currClus, int moreClus, cluster** clusters, int sizePerBlock);
void expandCluster(cluster* c);

#endif


#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include <cutil.h>

#include <cuda.h>

#include "driver_types.h"
#include "cuda_runtime.h"
#include "mem_info.h"

static unsigned long inKB(unsigned long bytes)
{ return bytes/1024; }

static unsigned long inMB(unsigned long bytes)
{ return bytes/(1024*1024); }

static unsigned long inGB(unsigned long bytes)
{ return bytes/(1024*1024*1024); }


static unsigned long convert_units(unsigned long bytes, char** units)
{
	unsigned long res;
	*units = "bytes";
	if(bytes>1024)
		if(bytes>(1024*1024))
		{
			res = inMB(bytes);
			*units = "MB";
			if(bytes>(1024*1024*1024))
			{
				res = inGB(bytes);	
				*units = "GB";
			}
		}
		else
		{
			res = inKB(bytes);
			*units = "KB";
		}
	return res;
}

static void printStats(unsigned long free, unsigned long total)
{
	char* units = (char*)malloc(sizeof(char)*4);

	unsigned long converted_free = convert_units(free, &units);

	printf(" Free : %lu %s\n",converted_free, units);

	unsigned long converted_total = convert_units(total, &units);

	printf(" Total : %lu %s\n" , converted_total, units);

	printf("%g%% free, %g%% used\n\n", 100.0*free/(double)total, 100.0*(total - free)/(double)total);
}

void printGPUMemoryInfo()
{
	unsigned int free, total;
	cuMemGetInfo(&free, &total);
	printStats(free,total);
}

void printCUDAErrors(const char* id)
{
	cudaError_t status = cudaGetLastError();

	if( status!= CUDA_SUCCESS )
	{
		printGPUMemoryInfo();
		printf("%s - %s\n",id,cudaGetErrorString(status));
	//	CUT_EXIT(d_argc,d_argv);
	}
}
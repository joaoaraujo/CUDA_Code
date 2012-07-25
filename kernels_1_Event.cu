#ifndef _KERNELS_1_EVENTS_H_
#define _KERNELS_1_EVENTS_H_

#define ACTIVE_BIT(i, j) (i & (1 << j))
#define EDIT_BIT(i, j) (i | (1 << j))

///////////////////////////////////////////////////////////////////////////////
//! Simple test kernel for device functionality
//! @param g_odata  memory to process (in and out)
///////////////////////////////////////////////////////////////////////////////


__global__ void
kernelBit1Event(int* tthreads, int3* predicates, int* n_predicates, float4* values, int3* anEvent, int2* sizeEvent, float4* valuesEvent, unsigned int* bitmap)
{
	const unsigned int tid = threadIdx.x;
	const unsigned int bid = blockIdx.x;
	const unsigned int pos = bid*blockDim.x + tid;
	const unsigned int nWarps = *tthreads/32;
	int oper;
	int atr;
	float4 val;
	float4 val2;

	if(pos == 0)
		bitmap[n_predicates[0]] = 0;

	for(int i = pos; i < *n_predicates; i+=(nWarps*32))
	{
		int3 currPred = predicates[i];
		val = values[currPred.z];
		val2 = values[currPred.z+1];
		oper = currPred.y;
		atr = currPred.x;

		int size = sizeEvent[0].y;

		for(int k = 0; k < size; k++)
		{
			int3 elemEvent = anEvent[k];

			if(elemEvent.x == atr)
			{
				for(int l = elemEvent.y; l < elemEvent.z; l+=2)
				{
					float4 currv = valuesEvent[l];
					float4 currv2 = valuesEvent[l+1];

					switch (oper)
					{
					case 0 : 
						if(val.x == currv.x)
							bitmap[i] = 1;
						break;
					case 1 :
						if(currv.x > val.x)
							bitmap[i] = 1;
						break;
					case 2 :
						if(currv.x < val.x)
							bitmap[i] = 1;
						break;
					case 3 :
						if(val.x != currv.x)
							bitmap[i] = 1;
						break;
					case 4 :
						if(currv.x <= val.x)
							bitmap[i] = 1;
						break;
					case 5 :
						if(currv.x >= val.x)
							bitmap[i] = 1;
						break;
					case 6 :
						if((val.x == currv.x && val.y == currv.y && val.z == currv.z && val.w == currv.w) && (val2.x == currv2.x && val2.y == currv2.y && val2.z == currv2.z && val2.w == currv2.w))
							bitmap[i] = 1;
						break;
					}
				}
			}
		}
	}
}


__global__ void
kernelFillBit1Event(int* size, int* randomPos, unsigned int* bitmap)
{
	for(int i = 0; i < size[0]; i++)
	{
		bitmap[randomPos[i]] = 1;
	}
}

__global__ void
kernelMatch1Event(int* tthreads, unsigned int* d_bitmap, int* subs, int* clusterSize, int* size, unsigned int* result)
{
	const unsigned int tid = threadIdx.x;
	const unsigned int bid = blockIdx.x;
	const unsigned int pos = bid*blockDim.x + tid;
	const unsigned int s = (clusterSize[0]/(size[0]+1));
	const unsigned int nWarps = *tthreads/32;

	for(unsigned int i = pos; i < s;)
	{
		int match = 0;
		int ID = subs[i*(size[0]+1)];

		for(unsigned int j = 1; j <= size[0]; j++)
		{
			if(d_bitmap[subs[i*(size[0]+1)+j]]==1)
				match++;
			else 
				break;
		}

		if(match == size[0]){
			result[ID] = 1;
		}

		i += (nWarps*32);
	}
}

__global__ void
kernelBit1EventShared(int* tthreads, int3* predicates, int* n_predicates, float4* values, int3* anEvent, int2* sizeEvent, float4* valuesEvent, unsigned int* bitmap)
{
	const unsigned int tid = threadIdx.x;
	const unsigned int bid = blockIdx.x;
	const unsigned int pos = bid*blockDim.x + tid;
	const unsigned int nWarps = *tthreads/32;

	int oper;
	int atr;
	float4 val;
	float4 val2;

	if(pos == 0)
		bitmap[n_predicates[0]] = 0;

	for(int i = pos; i < *n_predicates;i+=(nWarps*32))
	{
		int3 currPred = predicates[i];
		val = values[currPred.z];
		val2 = values[currPred.z+1];
		oper = currPred.y;
		atr = currPred.x;

		int size = sizeEvent[0].y;
		for(int k = 0; k < size; k++)
		{
			int3 elemEvent = anEvent[k];

			if(elemEvent.x == atr)
			{
				for(int l = elemEvent.y; l < elemEvent.z; l+=2)
				{
					float4 currv = valuesEvent[l];
					float4 currv2 = valuesEvent[l+1];

					switch (oper)
					{
					case 0 : 
						if((val.x == currv.x && val.y == currv.y && val.z == currv.z && val.w == currv.w) && (val2.x == currv2.x && val2.y == currv2.y && val2.z == currv2.z && val2.w == currv2.w))
							bitmap[i]= EDIT_BIT(bitmap[i],i);
						break;
					case 1 :
						if(currv.x > val.x)
							bitmap[i]= EDIT_BIT(bitmap[i],i);
						break;
					case 2 :
						if(currv.x < val.x)
							bitmap[i]= EDIT_BIT(bitmap[i],i);
						break;
					case 3 :
						if((val.x != currv.x || val.y != currv.y || val.z != currv.z || val.w != currv.w) || (val2.x != currv2.x || val2.y != currv2.y || val2.z != currv2.z || val2.w != currv2.w))
							bitmap[i]= EDIT_BIT(bitmap[i],i);
						break;
					case 4 :
						if(currv.x <= val.x)
							bitmap[i]= EDIT_BIT(bitmap[i],i);
						break;
					case 5 :
						if(currv.x >= val.x)
							bitmap[i]= EDIT_BIT(bitmap[i],i);
						break;
					case 6 :
						if((val.x == currv.x && val.y == currv.y && val.z == currv.z && val.w == currv.w) && (val2.x == currv2.x && val2.y == currv2.y && val2.z == currv2.z && val2.w == currv2.w))
							bitmap[i]= EDIT_BIT(bitmap[i],i);
						break;
					}
				}
			}
		}
	}
}


__global__ void
kernelFillBit1EventShared(int* size, int* randomPos, unsigned int* bitmap)
{
	int pos;

	for(int i = 0; i < size[0]; i++)
	{
		pos = randomPos[i]/32;
		bitmap[pos] = EDIT_BIT(bitmap[pos], randomPos[i]%32);//bitmap[pos] = 1;
	}
}


__global__ void
kernelMatch1EventShared(int* tthreads, unsigned int* d_bitmap, int* subs, int* c, int* size, int* bitsize, unsigned int* result)
{
	const unsigned int tid = threadIdx.x;
	const unsigned int bid = blockIdx.x;
	const unsigned int pos = bid*blockDim.x + tid;
	const unsigned int s = (c[0]/(size[0]+1));
	const unsigned int siz = size[0];
	const unsigned int nWarps = *tthreads/32;
	const unsigned int bits = *bitsize;

	extern __shared__ unsigned int b[];

	if(tid == 0)
	{
		for(int i = 0; i < bits; i++)
			b[i] = d_bitmap[i];
	}

	__syncthreads();

	for(unsigned int i = pos; i < s;)
	{
		int match = 0;
		int ID = subs[i*(siz+1)];
		int pos;

		for(unsigned int j = 1; j <= siz; j++)
		{
			pos = subs[i*(siz+1)+j];
			if((unsigned int)ACTIVE_BIT(b[pos/32], pos%32)>0)
				match++;
			else 
				break;
		}

		if(match == size[0]){
			result[ID] = 1;
		}

		i += (nWarps*32);
	}
}


__global__ void
kernelBit1EventsConstant(int* tthreads, int3* predicates, int* n_predicates, float4* values, unsigned int* bitmap)
{
		const unsigned int tid = threadIdx.x;
		const unsigned int bid = blockIdx.x;
		const unsigned int pos = bid*blockDim.x + tid;
		const unsigned int nWarps = *tthreads/32;

		int oper;
		int atr;
		float4 val;
		float4 val2;

		for(int i = pos; i < *n_predicates;)
		{
			int3 currPred = predicates[i];
			val = values[currPred.z];
			val2 = values[currPred.z+1];
			oper = currPred.y;
			atr = currPred.x;

				int size = sizeEvents_d[0].y;

				for(int k = 0; k < size; k++)
				{
					int3 elemEvent = events_d[k];

					if(elemEvent.x == atr)
					{
						for(int l = elemEvent.y; l < elemEvent.z; l+=2)
						{
							float4 currv = valuesEvents_d[l];
							float4 currv2 = valuesEvents_d[l+1];

							switch (oper)
							{
								case 0 : 
									if(val.x == currv.x)
										bitmap[i] = 1;
									break;
								case 1 :
									if(currv.x > val.x)
										bitmap[i] = 1;
									break;
								case 2 :
									if(currv.x < val.x)
										bitmap[i] = 1;
									break;
								case 3 :
									if(val.x != currv.x)
										bitmap[i] = 1;
									break;
								case 4 :
									if(currv.x <= val.x)
										bitmap[i] = 1;
									break;
								case 5 :
									if(currv.x >= val.x)
										bitmap[i] = 1;
									break;
								case 6 :
									if((val.x == currv.x && val.y == currv.y && val.z == currv.z && val.w == currv.w) && (val2.x == currv2.x && val2.y == currv2.y && val2.z == currv2.z && val2.w == currv2.w))
										bitmap[i] = 1;
									break;
							}
						}
					}
				}

			i += (nWarps*32);
		}

}

#endif // #ifndef _CPP_INTEGRATION_KERNEL_H_
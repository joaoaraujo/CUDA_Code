/* Example of integrating CUDA functions into an existing 
 * application / framework.
 * Device code.
 */

#ifndef _KERNELS_32_EVENTS_H_
#define _KERNELS_32_EVENTS_H_

#define ACTIVE_BIT(i, j) (i & (1 << j))
#define EDIT_BIT(i, j) (i | (1 << j))

/*extern __constant__ int3 events_d[320];//Eventos
extern __constant__ int2 sizeEvents_d[32];//Tamanho de cada Evento*/

///////////////////////////////////////////////////////////////////////////////
//! Simple test kernel for device functionality
//! @param g_odata  memory to process (in and out)
///////////////////////////////////////////////////////////////////////////////
__global__ void
kernelBit32EventsShared(int* tthreads, int3* predicates, int* n_predicates, float4* values, int3* anEvent, int2* sizeEvent, float4* vEvent, unsigned int* bitmap)
{
		const unsigned int tid = threadIdx.x;
		const unsigned int bid = blockIdx.x;
		const unsigned int pos = bid*blockDim.x + tid;
		const unsigned int nWarps = *tthreads/32;

		extern __shared__ float4 valuesEvent[];

		int oper;
		int atr;
		float4 val;
		float4 val2;

		if(pos == 0)
			bitmap[n_predicates[0]] = 0;


		if(tid == 0)
		{
			for(int i = 0; i < 32; i++)
			{
				int2 Event = sizeEvent[i];

				for(int k = Event.x; k < Event.y; k++)
				{
					int3 elemEvent = anEvent[k];

					for(int l = elemEvent.y; l < elemEvent.z; l+=2){
						valuesEvent[l]= vEvent[l];
						valuesEvent[l+1]= vEvent[l+1];
					}
				}
			}
		}

		__syncthreads();

		for(int i = pos; i < *n_predicates;)
		{
			int3 currPred = predicates[i];
			val = values[currPred.z];
			val2 = values[currPred.z+1];
			oper = currPred.y;
			atr = currPred.x;

			for(unsigned int j = 0; j < 32; j++)
			{
				int2 Event = sizeEvent[j];
				for(int k = Event.x; k < Event.y; k++)
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
											bitmap[i]=EDIT_BIT(bitmap[i],j);
									break;
								case 1 :
									if(currv.x > val.x)
											bitmap[i]=EDIT_BIT(bitmap[i],j);
									break;
								case 2 :
									if(currv.x < val.x)
											bitmap[i]=EDIT_BIT(bitmap[i],j);
									break;
								case 3 :
									if((val.x != currv.x || val.y != currv.y || val.z != currv.z || val.w != currv.w) || (val2.x != currv2.x || val2.y != currv2.y || val2.z != currv2.z || val2.w != currv2.w))
											bitmap[i]=EDIT_BIT(bitmap[i],j);
									break;
								case 4 :
									if(currv.x <= val.x)
											bitmap[i]=EDIT_BIT(bitmap[i],j);
										break;
								case 5 :
									if(currv.x >= val.x)
											bitmap[i]=EDIT_BIT(bitmap[i],j);
										break;
								case 6 :
									if((val.x == currv.x && val.y == currv.y && val.z == currv.z && val.w == currv.w) && (val2.x == currv2.x && val2.y == currv2.y && val2.z == currv2.z && val2.w == currv2.w))
											bitmap[i]=EDIT_BIT(bitmap[i],j);
										break;
							}
						}
					}
				}
			}

			i += (nWarps*32);
		}
}

__global__ void
kernelBit32EventsConstant(int* tthreads, int3* predicates, int* n_predicates, float4* values, unsigned int* bitmap)
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

			for(unsigned int j = 0; j < 32; j++)
			{
				int2 Event = sizeEvents_d[j];
				for(int k = Event.x; k < Event.y; k++)
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
									if((val.x == currv.x && val.y == currv.y && val.z == currv.z && val.w == currv.w) && (val2.x == currv2.x && val2.y == currv2.y && val2.z == currv2.z && val2.w == currv2.w))
											bitmap[i]=EDIT_BIT(bitmap[i],j);
									break;
								case 1 :
									if(currv.x > val.x)
											bitmap[i]=EDIT_BIT(bitmap[i],j);
									break;
								case 2 :
									if(currv.x < val.x)
											bitmap[i]=EDIT_BIT(bitmap[i],j);
									break;
								case 3 :
									if((val.x != currv.x || val.y != currv.y || val.z != currv.z || val.w != currv.w) || (val2.x != currv2.x || val2.y != currv2.y || val2.z != currv2.z || val2.w != currv2.w))
											bitmap[i]=EDIT_BIT(bitmap[i],j);
									break;
								case 4 :
									if(currv.x <= val.x)
											bitmap[i]=EDIT_BIT(bitmap[i],j);
										break;
								case 5 :
									if(currv.x >= val.x)
											bitmap[i]=EDIT_BIT(bitmap[i],j);
										break;
								case 6 :
									if((val.x == currv.x && val.y == currv.y && val.z == currv.z && val.w == currv.w) && (val2.x == currv2.x && val2.y == currv2.y && val2.z == currv2.z && val2.w == currv2.w))
											bitmap[i]=EDIT_BIT(bitmap[i],j);
										break;
							}
						}
					}
				}
			}

			i += (nWarps*32);
		}

}

__global__ void
kernelBit32Events(int* tthreads, int3* predicates, int* n_predicates, float4* values, int3* events, int2* sizeEvent, float4* valuesEvent, unsigned int* bitmap)
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

			for(unsigned int j = 0; j < 32; j++)
			{
				int2 Event = sizeEvent[j];
				for(int k = Event.x; k < Event.y; k++)
				{
					int3 elemEvent = events[k];
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
											bitmap[i]=EDIT_BIT(bitmap[i],j);
									break;
								case 1 :
									if(currv.x > val.x)
											bitmap[i]=EDIT_BIT(bitmap[i],j);
									break;
								case 2 :
									if(currv.x < val.x)
											bitmap[i]=EDIT_BIT(bitmap[i],j);
									break;
								case 3 :
									if((val.x != currv.x || val.y != currv.y || val.z != currv.z || val.w != currv.w) || (val2.x != currv2.x || val2.y != currv2.y || val2.z != currv2.z || val2.w != currv2.w))
											bitmap[i]=EDIT_BIT(bitmap[i],j);
									break;
								case 4 :
									if(currv.x <= val.x)
											bitmap[i]=EDIT_BIT(bitmap[i],j);
										break;
								case 5 :
									if(currv.x >= val.x)
											bitmap[i]=EDIT_BIT(bitmap[i],j);
										break;
								case 6 :
									if((val.x == currv.x && val.y == currv.y && val.z == currv.z && val.w == currv.w) && (val2.x == currv2.x && val2.y == currv2.y && val2.z == currv2.z && val2.w == currv2.w))
											bitmap[i]=EDIT_BIT(bitmap[i],j);
										break;
							}
						}
					}
				}
			}

			i += (nWarps*32);
		}

}

__global__ void
kernelFillBit32EventShared(int* n_values, int* eventNumber, unsigned int* randomPos, unsigned int* bitmap)
{
	const unsigned int tid = threadIdx.x;
	const unsigned int bid = blockIdx.x;
	const unsigned int pos = bid*blockDim.x + tid;

		int ran;

		for(int i = pos; i < n_values[0];i++)
		{
			ran = randomPos[i];
			bitmap[ran] = EDIT_BIT(bitmap[ran], eventNumber[0]);//bitmap[pos] = 1;
		}

}

__global__ void
kernelMatch32Events(int* tthreads, unsigned int* d_bitmap, int* subs, int* c, int* size, unsigned int* result)
{
		const unsigned int tid = threadIdx.x;
		const unsigned int bid = blockIdx.x;
		const unsigned int pos = bid*blockDim.x + tid;
		const unsigned int s = (c[0]/(size[0]+1));
		const unsigned int nWarps = *tthreads/32;
		
		for(unsigned int i = pos; i < s;)
		{
			int ID = subs[i*(size[0]+1)];
			unsigned int match = d_bitmap[subs[i*(size[0]+1)+1]];

			for(unsigned int j = 2; j <= size[0] && match != 0; j++)
				match &= d_bitmap[subs[i*(size[0]+1)+j]];
		
			result[ID] = match;

			i += (nWarps*32);
		}
}

__global__ void
kernelMatch32EventsDisjunctions(int* tthreads, unsigned int* subsAccepted,unsigned int* subs, int* numSubs, int* size, unsigned int* result)
{
		const unsigned int tid = threadIdx.x;
		const unsigned int bid = blockIdx.x;
		const unsigned int pos = bid*blockDim.x + tid;
		const unsigned int nWarps = *tthreads/32;
		
		for(unsigned int i = pos; i < *numSubs;)
		{
			int ID = subs[i*(size[0]+1)];
			unsigned int match = subsAccepted[subs[i*(size[0]+1)+1]];

			for(unsigned int j = 2; j <= size[0]; j++)
				match |= subsAccepted[subs[i*(size[0]+1)+j]];
		
			result[ID] = match;

			i += (nWarps*32);
		}
}


#endif // #ifndef _CPP_INTEGRATION_KERNEL_H_
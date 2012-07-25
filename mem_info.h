#ifndef _MEM_INFO_H_
#define _MEM_INFO_H_

static unsigned long inKB(unsigned long bytes);

static unsigned long inMB(unsigned long bytes);

static unsigned long inGB(unsigned long bytes);

static unsigned long convert_units(unsigned long bytes, char** units);

static void printStats(unsigned long free, unsigned long total);

void printGPUMemoryInfo();

void printCUDAErrors(const char* id);

#endif

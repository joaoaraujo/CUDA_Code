// includes, system
#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <string>
#include <fstream>
#include <map>
#include <math.h>

#include <./include/CkCrypt2.h>

// Required to include CUDA vector types
#include <vector_types.h>
#include "cutil_inline.h"
#include <./clusters_struct.h>
//#include <Functions.cu>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// declaration, forward

extern "C" 
void readFile(char* fileName, char*** strings, int* numberStrings);
extern "C"
int hexToInt(char* str);
extern "C"
int4* convertString(char* s);
extern "C"
int operatorNumber(char* operat);
extern "C"
int logBase2(int x);

int hexToInt(char* str)
{
	int res[8];
	char values[] = {"0123456789ABCDEF"};

	for(int i = 0; i < 8; i++)
	{
		for(int j = 0; j < 16; j++)
		{
			if(values[j] == str[i]){
					res[i] = j;
					break;
			}
		}
	}

	int result = res[0] << 28 | res[1] << 24 | res[2] << 20 | res[3] << 16 | res[4] << 12 | res[5] << 8 | res[6] << 4 | res[7];

	return result;
}

int4* convertString(char* s)
{

		CkCrypt2 crypt;
		//  Any string argument automatically begins the 30-day trial.
		bool success;
		success = crypt.UnlockComponent("30-day trial");
		if (success != true) {
			printf("Crypt component unlock failed\n");
			return 0;
		}

		crypt.put_HashAlgorithm("sha256");
		crypt.put_EncodingMode("hex");
		char* hash;
		hash = (char*)crypt.hashStringENC(s);

		int4* finalV = (int4*)malloc(sizeof(int4)*2);

		for(int i = 0; i < 64; i+=8){

			char* integer = (char*)malloc(sizeof(char)*4);

			for(int j = 0; j < 8; j++)
				integer[j] = hash[j+i];

			integer[4] = '\0';
			int finalInt = hexToInt(integer);

			if(i < 32){
				if(i == 0)
					finalV[0].x = finalInt;
				else if(i == 8)
					finalV[0].y = finalInt;
				else if(i == 16)
					finalV[0].z = finalInt;
				else if(i == 24)
					finalV[0].w = finalInt;
			}
			else{
				if(i == 32)
					finalV[1].x = finalInt;
				else if(i == 40)
					finalV[1].y = finalInt;
				else if(i == 48)
					finalV[1].z = finalInt;
				else if(i == 56)
					finalV[1].w = finalInt;
			}
		}
			//crypt.ClearEncryptCerts();
			return finalV;
}

int operatorNumber(char* operat)
{
	int o;

	if(strcmp(operat, "=") == 0)
		o = 0;
	else if(strcmp(operat, ">") == 0)
		o = 1;
	else if(strcmp(operat, "<") == 0)
		o = 2;
	else if(strcmp(operat, "!=") == 0)
		o = 3;
	else if(strcmp(operat, "<=") == 0)
		o = 4;
	else if(strcmp(operat, ">=") == 0)
		o = 5;
	else if(strcmp(operat, "C|") == 0)
		o = 6;

	return o;
}

//log x para base 2 = ln(x) / ln(2)
int logBase2(int x)
{
	double val = log((double)x) / log(2.0);
	double res;
	double dec = modf(val,&res);

	if(dec == 0.0 && res != 0)
		res--;

	return (int)res;
}

void readFile(char* fileName, char*** strings, int* numberStrings){

	FILE* fp = fopen(fileName,"r+");
	numberStrings[0] = 0;
	strings[0] = (char**)malloc(500000*sizeof(int));

	if(fp == NULL){
		printf("Please insert a correct file name!\n");
		return;
	}

	char s[500];

	while(fgets(s,sizeof(s),fp))
	{
		strings[0][numberStrings[0]] = (char*)malloc(sizeof(char)*(strlen(s)+1));
		
		strcpy(strings[0][numberStrings[0]],s);
		numberStrings[0]++;
	}

	fclose(fp);
}


//TODO

//Por a biblioteca de clusters_struct no lado do device, enviando apenas o int* criado no host, como ja esta....

//testar a estrutura em memoria enviando novas subscricoes e afins...

//criar leitor de ficheiros com configuracoes...

//Fazer de modo a que dê para correr correr vários eventos de cada vez, reaproveitando as estruturas.
//Dinamico? entrada e saida de subscrições, por consequencia predicados.

//Tratar dos erros do numero de threads

//testes

/* 

Jornal = Record
Noticia |= O Benfica ganhou o campeonato mais uma vez, desta vez com 30 pontos de avanco!
*/
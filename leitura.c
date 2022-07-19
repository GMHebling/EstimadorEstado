#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "matriz.h"
#include "leitura.h"


char *config()
{
    FILE *arquivo=NULL;
    arquivo = fopen("config.txt","r");

    if (arquivo == NULL) return "-1";
    char p;
    char *path;
    path=(char*) malloc(100*sizeof(char));
    while ((p=fgetc(arquivo))!=EOF)
    {
        if (p!='\n')strncat(path,&p,1);
        
    }
    return path;
    
}


int r_dmatfcsv(double ***mat,char path[100],char sep,int *nlin,int *ncol)//Lê uma matriz de double em csv
{
    
    FILE *arquivo=NULL;
    char p;
    char num[100];
    num[0]='\0';
    *nlin=0;
    *ncol=1;
    int i =0,j=0;
    arquivo = fopen(path,"r");
   
    if (arquivo==NULL) return 1;
    while ((p=fgetc(arquivo))!=EOF)
    {
        if (p==sep && (*nlin)==0) (*ncol)++;
        if (p=='\n') (*nlin)++;
    }

    (void) matrizDinamica(mat,*nlin,*ncol);

    rewind(arquivo);


    while ((p=fgetc(arquivo))!=EOF)
    {
         
        if (p==sep){
            
            (*mat)[i][j]=(double) atof(num);
            j++;
            num[0]='\0';
        }
        else if (p=='\n') 
        {
            (*mat)[i][j]=(double) atof(num);
            i++;
            j=0;
            num[0]='\0';
        }
        else 
        {
            strncat(num,&p,1);
        }
    }

    fclose(arquivo);
    return 0;
}

void imprimirmat(double** mat, int linha, int col)
{
	int i, j;
	for (i = 0; i < linha; i++)
	{

		for (j = 0; j < col; j++)
		{
			printf("%.2e\t ", mat[i][j]);
		}
		printf("\n");
	}
}

void fimprimirmat(FILE* arquivo, double** mat, int linha, int col)
{
	int i, j;
	for (i = 0; i < linha; i++)
	{

		for (j = 0; j < col; j++)
		{
            if (j+1==col) fprintf(arquivo, "%.2e\n ", mat[i][j]);
            else fprintf(arquivo, "%.2e\t ", mat[i][j]);

		}
		
	}
}
void imprimirvet_double(double* vet, int linha)
{
	int i;
	for (i = 0; i < linha; i++)
	{
		printf("%e ", vet[i]);
		printf("\n");
	}
}

void fimprimirvet_double(FILE* arquivo, double* vet, int linha)
{
	int i;
	for (i = 0; i < linha; i++)
	{
		fprintf(arquivo, "%e", vet[i]);
		fprintf(arquivo, "\n");
	}
}

void imprimirvet_int(int* vet, int linha)
{
	int i;
	for (i = 0; i < linha; i++)
	{
		printf("%d ", vet[i]);
		printf("\n");
	}
}

void fimprimirvet_int(FILE* arquivo, int* vet, int linha)
{
	int i;
	for (i = 0; i < linha; i++)
	{
		fprintf(arquivo, "%d", vet[i]);
		fprintf(arquivo, "\n");
	}
}


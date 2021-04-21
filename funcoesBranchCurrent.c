/* 
 * File:   funcoesBranchCurrent.c
 * Author: Gustavo Hebling
 *
 * Created on 15 de Abril de 2021, 18:20
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <complex.h>


#include "data_structures.h"
#include "funcoesWLS.h"
#include "funcoesTopologia.h"
#include "funcoesCalculoEletrico.h"
#include "funcoesOtimizacao.h"
#include "funcoesMatematicas.h"
#include "funcoesFluxoVarredura.h"

//#include "mmio.h"
#include <cholmod.h>
#include "SuiteSparseQR_C.h"


#include "data_structures.h"




int *montaRNP(ALIMENTADOR alimentadores){
    int *RNP;
    int k = 0;
    FILABARRAS *barraAtual;

    RNP = aloca_vetor_int(alimentadores.numeroNos+1);
    barraAtual = &alimentadores.rnp[0];
    while(barraAtual != NULL){
        RNP[k] = barraAtual->idNo;
        k++;
        barraAtual = barraAtual->prox;
    }

    return RNP;
}

void estimadorBC_RECT(GRAFO *grafo, long int numeroRamos, long int numeroBarras, DMED *medidas, long int **numeroMedidas, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos, double Sbase)
{
    long int nmed, nvar;
    int i, j, k, r;
    double *z = NULL, **h = NULL, ***H = NULL, **W = NULL, *x = NULL, *regua = NULL, aux = 0;

    printf("Estimador de Estado Branch Current em Coordenadas retangulares...\n");
    //--------------------------------------------------------------------------
    //Alocação de memória das variáveis do estimador de estado
    nmed = 0;
    for (i = 0; i < 9; i++)
    {
        for (j = 0; j < 8; j++)
        {
            nmed = nmed + numeroMedidas[i][j];
        }
    }
    nvar = 0;
    //printf("numero barras: %d\n", numeroBarras);
    for (i = 0; i < numeroBarras; i++)
    {
        switch (grafo[i].fases)
        {
        case 1:
            nvar += 2;
            break;
        case 2:
            nvar += 2;
            break;
        case 3:
            nvar += 2;
            break;
        case 4:
            nvar += 4;
            break;
        case 5:
            nvar += 4;
            break;
        case 6:
            nvar += 4;
            break;
        case 7:
            nvar += 6;
            break;
        }
    }
    //printf("nmed: %d\n", nmed);
    //printf("nvar: %d\n", nvar);
    if ((z = (double *)malloc((nmed) * sizeof(double))) == NULL)
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor z!!!!");
        exit(1);
    }
    if ((h = malloc((nmed) * sizeof(double *))) == NULL)
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor h!!!!");
        exit(1);
    }
    if ((x = (double *)malloc((nvar) * sizeof(double))) == NULL)
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor x!!!!");
        exit(1);
    }
    if ((regua = (double *)malloc((nvar) * sizeof(double))) == NULL)
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor regua!!!!");
        exit(1);
    }

    H = (double ***)malloc(nmed * sizeof(double **));
    for (i = 0; i < nmed; i++)
    {
        H[i] = (double **)malloc(nvar * sizeof(double *));
        for (j = 0; j < nvar; j++)
        {
            H[i][j] = &aux;
        }
    }

    //Inicializa vetor x (correntes)
    //utilizar variavel numeroRamos

    
    //RNP - a partir do alimentador
    int *RNP;
    RNP = montaRNP(*alimentadores);

    //inicializar vetor de variaveis de estado

    
    

    //TODO: Montar a regua conforme a estrutura DRAM - Ramos 
    //--------------------------------------------------------------------------
    // Direcionamento dos ponteiros que compõem a régua, vetor h e matriz H
    //monta a regua (ordenaçao do vetor x) - V e teta conforme o grafo
    j = 0;
    for (i = 0; i < numeroBarras; i++)
    {
        aux = (double)grafo[i].idNo;
        aux += 0.01;
        switch (grafo[i].fases)
        {
        case 1:
            regua[j] = aux;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            break;
        case 2:
            regua[j] = aux + 0.1;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            break;
        case 3:
            regua[j] = aux + 0.2;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            break;
        case 4:
            regua[j] = aux;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            regua[j] = aux + 0.1;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            break;
        case 5:
            regua[j] = aux;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            regua[j] = aux + 0.2;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            break;
        case 6:
            regua[j] = aux + 0.1;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            regua[j] = aux + 0.2;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            break;
        case 7:
            regua[j] = aux;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            regua[j] = aux + 0.1;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            regua[j] = aux + 0.2;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            break;
        }
    }
    aux = 0;

    //Tratamento da referência
    long int ref_1, ref_2;
    tratamento_referencia(&ref_1, &ref_2, &alimentadores[0], regua, nvar);
    tira_refs_regua(nvar, ref_1, ref_2, regua);


    nvar = nvar - (ref_2 - ref_1 + 1);
    //printf("tira refs\n");
    //vetor h aponta para a estrutura de dados das medidas
    for (i = 0; i < nmed; i++)
    {
        h[i] = &medidas[i].h;
    }
    //Matriz H aponta para a estrutura de dados das medidas
    for (i = 0; i < nmed; i++)
    {
        for (j = 0; j < medidas[i].nvar; j++)
        {
            for (r = 0; r < nvar; r++)
            {
                if (cabs(medidas[i].reguaH[j] - regua[r]) < EPS)
                {
                    H[i][r] = &medidas[i].H[j];
                    break;
                }
            }
        }
    }
    //--------------------------------------------------------------------------
    //Estimação de Estado
    monta_z(z, nmed, medidas);
    monta_W(NULL, nmed, medidas);
    //monta_W_cte(W,nmed,medidas);
    //monta_W_Ident(NULL,nmed,medidas);

    incializa_vetor_x(grafo, numeroBarras, alimentadores, numeroAlimentadores, x, regua, nvar);

    double tol = 0.000001;
    clock_t tIni = clock();
    int conv = otimiza_Gauss_NewtonQR(z, h, H, grafo, numeroBarras, ramos, medidas, nvar, nmed, regua, x, tol, ref_1, ref_2);

    clock_t t1 = clock();
    double tempoWLS = (double)(t1 - tIni) / CLOCKS_PER_SEC;
    printf("\nEstimação WLS: %lf", tempoWLS);

    exportaCasoReferencia(grafo, numeroBarras, Sbase);
    exportaEstado(grafo, regua, nvar);
    imprimeEstado(grafo, numeroBarras);

    FILE *residuo;
    residuo = fopen("residuo.txt", "wt");
    for (i = 0; i < nmed; i++)
    {
        fprintf(residuo, "%.10lf\n", z[i] - medidas[i].h);
    }
    fclose(residuo);

    free(z);
    free(h);
    free(x);
    free(regua);
    for (i = 0; i < nmed; i++)
        free(H[i]);
    free(H);
}

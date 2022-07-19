#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <complex.h>

#include "data_structures.h"
#include "funcoesWLS.h"
#include "fluxoNRQR.h"
#include "leitura.h"
#include "matriz.h"
#include "funcoesTopologia.h"
#include "funcoesCalculoEletrico.h"
#include "funcoesOtimizacao.h"
#include "funcoesMatematicas.h"



void fluxoNRQR(GRAFO *grafo, long int numeroBarras, DMED *medidas, long int **numeroMedidas, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos, double Sbase)
{
    long int nmed, nvar,tratref;
    int i, j, k, r,nref,pseu,nang=0,nmodv=0;
    long int ref_1, ref_2;
    double *z = NULL, **h = NULL, ***H = NULL, **W = NULL, *x = NULL, *regua = NULL,*regua_o=NULL, aux = 0,*regua_rem=NULL;
    FILE *reguaout;
    reguaout=fopen("regua.csv","w");


    
    printf("Fluxo de Potência Trifásico...\n");
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
        if( grafo[i].tipo==1)
        {
            if(grafo[i].fases==1||grafo[i].fases==2||grafo[i].fases==3)
            {
                nmodv++;
            }
            else if(grafo[i].fases==4||grafo[i].fases==5||grafo[i].fases==6)
            {
                nmodv=nmodv+2;
            }
            else if(grafo[i].fases==7)
            {
                nmodv=nmodv+3;
            }
        }
        else if (grafo[i].tipo==2)
        {
            if(grafo[i].fases==1||grafo[i].fases==2||grafo[i].fases==3)
            {
                nang++;
                nmodv++;
            }
            else if(grafo[i].fases==4||grafo[i].fases==5||grafo[i].fases==6)
            {
                nang=nang+2;
                nmodv=nmodv+2;
            }
            else if(grafo[i].fases==7)
            {
                nang=nang+3;
                nmodv=nmodv+3;
            }   
        }
        

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

    if ((z = (double *)malloc((nmed+2) * sizeof(double))) == NULL)
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor z!!!!");
        exit(1);
    }
    if ((h = malloc((nmed+2) * sizeof(double *))) == NULL)
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
    if ((regua_o = (double *)malloc((nvar) * sizeof(double))) == NULL)
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor regua!!!!");
        exit(1);
    }
    if ((regua_rem = (double *)malloc((nmodv+nang) * sizeof(double))) == NULL)
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor regua de remocao!!!!");
        exit(1);
    }


    H = (double ***)malloc((nmed+2) * sizeof(double **));
    for (i = 0; i < nmed+2; i++)
    {
        H[i] = (double **)malloc(nvar * sizeof(double *));
        for (j = 0; j < nvar; j++)
        {
            H[i][j] = &aux;
        }
    }
    //--------------------------------------------------------------------------
    // Direcionamento dos ponteiros que compõem a régua, vetor h e matriz H
    //monta a regua (ordenaçao do vetor x) - V e teta conforme o grafo
    j = 0;
    k=0;
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
        // regua com as barras PV
        if( grafo[i].tipo==1)
        {
            switch (grafo[i].fases)
            {
            case 1:
                regua_rem[k] = aux;
                k++;
                break;
            case 2:
                regua_rem[k] = aux + 0.1;
                k++;
                break;
            case 3:
                regua_rem[k] = aux + 0.2;
                k++;
                break;
            case 4:
                regua_rem[k] = aux;
                k++;
                regua_rem[k] = aux + 0.1;
                k++;
                break;
            case 5:
                regua_rem[k] = aux;
                k++;
                regua_rem[k] = aux + 0.2;
                k++;
                break;
            case 6:
                regua_rem[k] = aux + 0.1;
                k++;
                regua_rem[k] = aux + 0.2;
                k++;
                break;
            case 7:
                regua_rem[k] = aux;
                k++;
                regua_rem[k] = aux + 0.1;
                k++;
                regua_rem[k] = aux + 0.2;
                k++;
                break;
            }   
        }
        //Regua com as variáveis da slack
        else if (grafo[i].tipo==2)
        {
        switch (grafo[i].fases)
        {
        case 1:
            regua_rem[k] = aux;
            k++;
            regua_rem[k] = -regua_rem[k-1];
            k++;
            break;
        case 2:
            regua_rem[k] = aux + 0.1;
            k++;
            regua_rem[k] = -regua_rem[k-1];
            k++;
            break;
        case 3:
            regua_rem[k] = aux + 0.2;
            k++;
            regua_rem[k] = -regua_rem[k-1];
            k++;
            break;
        case 4:
            regua_rem[k] = aux;
            k++;
            regua_rem[k] = -regua_rem[k-1];
            k++;
            regua_rem[k] = aux + 0.1;
            k++;
            regua_rem[k] = -regua_rem[k-1];
            k++;
            break;
        case 5:
            regua_rem[k] = aux;
            k++;
            regua_rem[k] = -regua_rem[k-1];
            k++;
            regua_rem[k] = aux + 0.2;
            k++;
            regua_rem[k] = -regua_rem[k-1];
            k++;
            break;
        case 6:
            regua_rem[k] = aux + 0.1;
            k++;
            regua_rem[k] = -regua_rem[k-1];
            k++;
            regua_rem[k] = aux + 0.2;
            k++;
            regua_rem[k] = -regua_rem[k-1];
            k++;
            break;
        case 7:
            regua_rem[k] = aux;
            k++;
            regua_rem[k] = -regua_rem[k-1];
            k++;
            regua_rem[k] = aux + 0.1;
            k++;
            regua_rem[k] = -regua_rem[k-1];
            k++;
            regua_rem[k] = aux + 0.2;
            k++;
            regua_rem[k] = -regua_rem[k-1];
            k++;
            break;
        }
        }
    }
    for(i=0;i<(int)nvar / 2;i++)
    {
        regua_o[i] = regua[i];
        regua_o[i + (int)nvar / 2] = -regua_o[i];
    }
    aux = 0;

    //for (i=0;i<nvar;i++)
    //{
    //    fprintf(reguaout,"%f\n",regua_o[i]);
    //}
    //for (i=0;i<nmodv+nang;i++)
    //{
    //    printf("%.3f\n",regua_rem[i]);
    //}

    

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

   
    // To Do, inicia barras V das PV e coloca 1 pra todo mundo

    //conferir tensoes iniciais & fazer inicialização FP
    incializa_vetor_x_FP(grafo, numeroBarras, alimentadores, numeroAlimentadores, x, regua, nvar);

    //Retira da regua os angulos de referencia e os modulos de tensão


    
    

    atualiza_H(grafo, numeroBarras, ramos, medidas, nmed);

    
    //--------------------------------------------------------------------------
 

    

    FILE *TWLS;
    TWLS=fopen("twls.csv","a");
    double tol = 1e-6;
    clock_t tic  = clock();
    
    //alterar o sistema de equacoes
    //(void) otimiza_Gauss_NewtonQR(z, h, H, grafo, numeroBarras, ramos, medidas, nvar, nmed, regua, x, tol, ref_1, ref_2,nref);

    (void) otimiza_Gauss_NewtonQR_FP(z,h,H,grafo,numeroBarras,ramos, medidas, nvar,nmed, regua, x,tol,regua_rem,nmodv+nang);

    clock_t toc = clock();

    double tempoWLS = (double)(toc - tic) / CLOCKS_PER_SEC;
    printf("\nFluxo de Potencia: %lf", tempoWLS);
    
    ///fprintf(TWLS,"%.2e\n",tempoWLS);

    exportaCasoReferencia(grafo, numeroBarras, Sbase);
    exportaRamoscorrentes(grafo, numeroBarras, Sbase);
    exportaEstado(grafo, regua, nvar);
    imprimeEstado(grafo, numeroBarras);
    //imprimeInjecoes(grafo, numeroBarras);

    FILE *residuo;
    residuo = fopen("residuo.txt", "wt");
    for (i = 0; i < nmed; i++)
    {
        fprintf(residuo, "%.10lf\n", z[i] - medidas[i].h);
    }
    fclose(residuo);
    fclose(TWLS);
    free(z);
    free(h);
    free(x);
    free(regua);
    for (i = 0; i < nmed; i++)
        free(H[i]);
    free(H);
}

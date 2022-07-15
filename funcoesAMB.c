/*
 * File:   funcoesBranchCurrent.c
 * Author: Gustavo Hebling
 *
 * Created on 15 de Julho de 2022, 08:45
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
#include "funcoesBranchCurrent.h"

//#include "mmio.h"
#include <cholmod.h>
#include "SuiteSparseQR_C.h"

#include "data_structures.h"


void estimadorAMB(GRAFO *grafo, long int numeroRamos, long int numeroBarras, DMED *medidas, long int **numeroMedidas, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos, double Sbase, DBAR *barra)
{
    long int nmed, nvar;
    int i, j, k, r;
    double *z = NULL, **h = NULL, ***H = NULL, **W = NULL, *x = NULL, *regua = NULL, aux = 0;
    double tol = 0.000001;
    //__complex__ double *x_bc = NULL;

    long int nmed_BC;
    long int nmed_T;

    int *caminho = NULL;
    int *nadj_proxbarra = NULL;
    caminho = malloc(numeroBarras * sizeof(int));
    nadj_proxbarra = malloc(numeroBarras * sizeof(int));
    busca_loop_grafo(grafo, numeroRamos, numeroBarras, caminho, nadj_proxbarra);

    printf("Estimador de Estado Branch Current em Coordenadas retangulares com medidas de tensão...\n");
    //--------------------------------------------------------------------------
    // Alocação de memória das variáveis do estimador de estado

    nmed = 0;
    nvar = 0;
    for (i = 0; i < 9; i++)
    {
        for (j = 0; j < 8; j++)
        {
            nmed = nmed + numeroMedidas[i][j];
        }
    }
    printf("numero barras: %d\n", numeroBarras);
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

    printf("nmed: %d\n", nmed);
    printf("nvar: %d\n", nvar);
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

    // Inicializa vetor x (correntes)
    // utilizar variavel numeroRamos

    // RNP - a partir do alimentador
    int *RNP;
    // RNP = aloca_vetor(numeroBarras);
    k = 0;
    FILABARRAS *barraAtual;
    RNP = aloca_vetor_int(numeroBarras);
    barraAtual = &(alimentadores[0].rnp[0]);
    while (barraAtual != NULL)
    {
        RNP[k] = barraAtual->idNo;
        k++;
        barraAtual = barraAtual->prox;
    }

    nmed_BC = conta_medidas_BC(medidas, nmed);
    // cria_B_Z_ramos(grafo, numeroRamos, ramos, Sbase);

    // x_bc = aloca_vetor(numeroRamos);
    // inicializa_vetor_estados_BC(x_bc, 3*numeroRamos);
    // inicializar vetor de variaveis de estado

    DMED_COMPLEX *medidas_complexas = NULL;
    medidas_complexas = (DMED_COMPLEX *)malloc((nmed_BC) * sizeof(DMED_COMPLEX));

    double *x_bc = NULL;
    x_bc = aloca_vetor(6 * numeroRamos);

    double *delta_x_bc = NULL;
    delta_x_bc = aloca_vetor(6 * numeroRamos);
    // x_bc = (double *)malloc((6*numeroRamos) * sizeof(double));

    DMED_COMPLEX *medidas_equivalentes = NULL;
    medidas_equivalentes = (DMED_COMPLEX *)malloc((nmed_BC) * sizeof(DMED_COMPLEX));

    __complex__ double *z_eq = NULL;
    z_eq = c_vetAloca(3 * nmed_BC);

    double *z_eq_tensao = NULL;
    z_eq_tensao = aloca_vetor(6 * nmed_BC + nmed_T);

    DMED_COMPLEX *medidas_tensao = NULL;
    nmed_T = conta_medidas_Tensao(medidas, nmed);
    medidas_tensao = (DMED_COMPLEX *)malloc((nmed_T) * sizeof(DMED_COMPLEX));

    // converte medidas tensao (magnitude) em complexas

    // z_eq = (__complex__ double *)malloc(3*nmed_BC * sizeof(__complex__ double));

    double *regua_x = NULL;
    double *regua_med = NULL;
    double *regua_med_inv = NULL;
    double **H_BC = NULL;
    double **H_T = NULL;

    double *x_anterior = NULL;
    x_anterior = aloca_vetor(6 * numeroRamos);
    double *dif_x = NULL;
    dif_x = aloca_vetor(6 * numeroRamos);
    // vetor de estados: 1 para cada ramo e fase;
    regua_med = aloca_vetor(3 * nmed_BC);
    regua_med_inv = aloca_vetor(3 * nmed_BC);
    regua_x = aloca_vetor(3 * numeroRamos);
    H_BC = aloca_matriz(3 * nmed_BC, 3 * numeroRamos);
    H_T = aloca_matriz(nmed_T, 6 * numeroRamos);

    incializa_tensoes_grafo(grafo, numeroBarras, alimentadores, numeroAlimentadores);
    // printf("1\n");
    medidas_complexas = converte_medidas_para_complexo(medidas, nmed);
    //
    medidas_tensao = calcula_medida_tensao_complexa(medidas, nmed, grafo, numeroBarras);

    double *regua_V = NULL;
    regua_V = (double *)malloc(nmed_T * sizeof(double));
    // regua das medidas de tensao
    monta_regua_medidas_tensao(medidas_tensao, nmed_T, regua_V);

    // printf("2\n");
    medidas_equivalentes = divide_medidas_por_tensao(medidas_complexas, nmed_BC, numeroBarras, grafo);

    printf("3\n");

    // regua vetor x
    monta_regua_x(numeroRamos, regua_x, ramos);
    // printf("4\n");
    monta_regua_medidas(nmed_BC, regua_med, regua_med_inv, medidas_equivalentes);

    double *regua_caminho = NULL;
    regua_caminho = aloca_vetor(numeroBarras);
    monta_regua_caminho(numeroRamos, numeroBarras, regua_caminho, caminho, grafo);
    // printf("5\n");
    H_BC = monta_matriz_H(numeroRamos, nmed_BC, regua_x, regua_med, regua_med_inv);
    //vetor h(x) das medidas de corrente
    double *hx_I = NULL;
    hx_I = (double*)malloc(6*nmed_BC*sizeof(double));
    //vetor h(x) das medidas de tensão
    double *hx_V = NULL;
    hx_V = (double*)malloc(nmed_T*sizeof(double));
    // montar H das medidas de tensão
    H_T = monta_matriz_H_tensao(numeroBarras, numeroRamos, nmed_T, caminho, medidas_tensao, regua_x, regua_caminho, ramos, grafo, hx_V);


    // printf("6\n");
    int it = 0;
    int conv = 0;
    while (conv < 1)
    {

        // monta_z_complexa(medidas_equivalentes, z_eq, nmed_BC);

        monta_z_real_e_imag(medidas_equivalentes, z_eq_tensao, nmed_BC, medidas_tensao, nmed_T);

        int st = 0;

        //x_anterior = x_bc;
        
        // x_bc = resolve_linear_QR(H_BC, z_eq, numeroRamos, nmed_BC);
        //
        delta_x_bc = resolve_linear_QR_Tensao(H_BC, H_T, z_eq_tensao, numeroRamos, nmed_BC, nmed_T, hx_I, hx_V);

        atualiza_vetor_x(x_bc, delta_x_bc, numeroRamos);


        calcula_hx_corrente(H_BC, x_bc, hx_I, nmed_BC, numeroRamos);

        // atualiza matriz jacobiana de derivadas das medidas de tensao
        H_T = monta_matriz_H_tensao(numeroBarras, numeroRamos, nmed_T, caminho, medidas_tensao, regua_x, regua_caminho, ramos, grafo, hx_V);

        // for (int cx = 0; cx < 6 * numeroRamos; cx++)
        // {
        //     dif_x[cx] = x_anterior[cx] - x_bc[cx];
        // }
        double nfx;

        nfx = norma_inf(delta_x_bc, 6 * numeroRamos);
        printf("\n\nIteracao:  %d \t|Dx|_inf =  %.17lf \t  \n", it, nfx);

        // mudar atualiza rede para receber complexo
        // atualiza_Rede_BC(grafo, numeroBarras, barra, regua_x, numeroRamos, x_bc);
        atualiza_Rede_BC_Tensao(grafo, numeroBarras, barra, regua_x, numeroRamos, x_bc);


        for (int k = 0; k < numeroBarras; k++)
        {
            int ct = forward_sweep(&grafo[RNP[k]], grafo);
        }
        if (nfx<tol | it> 70)
        {
            conv = 10;
        }

        medidas_equivalentes = divide_medidas_por_tensao(medidas_complexas, nmed_BC, numeroBarras, grafo);

        it++;
    }
}

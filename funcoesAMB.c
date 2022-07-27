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

//#include "mmio.h"
#include <cholmod.h>
#include "SuiteSparseQR_C.h"

#include "data_structures.h"

double *monta_regua_x_AMB(GRAFO *grafo, long int numeroBarras){
    double *regua_x;
    regua_x = (double *)malloc(numeroBarras * sizeof(double));

    for (int i = 0; i < numeroBarras; i++){
        regua_x[i] = grafo[i].barra->ID;
    }
    return regua_x;
}

double *monta_vetor_x_inicial_AMB(long int numeroBarras){
    double *x_init;
    x_init = (double *)malloc(6*numeroBarras * sizeof(double));
    for (int i = 0; i < numeroBarras; i++){
        x_init[3*i] = 1.0;
        x_init[3*i + (3*numeroBarras)] = 0.0;
    }
    return x_init;
}

DMED_COMPLEX *calcula_medida_tensao_complexa_AMB(DMED *medidas, long int numeroMedidas, GRAFO *grafo, int numeroBarras)
{
    DMED_COMPLEX *medida_tensao = NULL;

    int cont_med_tensao;
    cont_med_tensao = conta_medidas_Tensao(medidas, numeroMedidas);

    printf("nmed_tensao = %d\n", cont_med_tensao);
    if (((medida_tensao) = (DMED_COMPLEX *)malloc((cont_med_tensao) * sizeof(DMED_COMPLEX))) == NULL)
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para as medidas complexas !!!!");
        exit(1);
    }

    int aux_contador = 0;
    int med_found = 0;
    for (int contador = 0; contador < numeroMedidas; contador++)
    {
        if (medidas[contador].tipo == 4)
        {
            medida_tensao[aux_contador].ligado = medidas[contador].ligado;
            medida_tensao[aux_contador].tipo = medidas[contador].tipo;
            medida_tensao[aux_contador].DE = medidas[contador].DE;
            medida_tensao[aux_contador].PARA = medidas[contador].PARA;
            medida_tensao[aux_contador].fases = medidas[contador].fases;
            medida_tensao[aux_contador].id = medidas[contador].id;
            medida_tensao[aux_contador].par = medidas[contador].par;

            medida_tensao[aux_contador].sigma = medidas[contador].sigma;
            medida_tensao[aux_contador].prec = medidas[contador].prec;

            double magnitude_tensao = medidas[contador].zmed;

            int idx_barra_DE = calcula_idx_ID(grafo, medida_tensao[aux_contador].DE, numeroBarras);

            int f_idx = 0;

            switch (medida_tensao[aux_contador].fases)
            {
            case 1:
                f_idx = 0;
                break;
            case 2:
                f_idx = 1;
                break;
            case 3:
                f_idx = 2;
                break;
            }
            __complex__ double tensao_de = grafo[idx_barra_DE].V[f_idx];
            double abs_tensao_de = abs(tensao_de);

            double parte_real = creal(tensao_de)/abs_tensao_de;
            double parte_imag = cimag(tensao_de)/abs_tensao_de;

            medida_tensao[aux_contador].zmed = (magnitude_tensao*parte_real) + (magnitude_tensao*parte_imag) * I;

            aux_contador += 1;
        }
    }
    return medida_tensao;
}

double calcula_G_AMB(int fase, int i_ramo, DRAM *ramos){
    double R, X;
    R = creal(ramos[i_ramo].Z[fase][fase]);
    X = cimag(ramos[i_ramo].Z[fase][fase]);
    double G;
    G = R/((R*R) + (X*X));
    return G;
}
double calcula_B_AMB(int fase, int i_ramo, DRAM *ramos){
    double R, X;
    R = creal(ramos[i_ramo].Z[fase][fase]);
    X = cimag(ramos[i_ramo].Z[fase][fase]);
    double B;
    B = -X/((R*R) + (X*X));
    return B;
}

double **monta_matriz_H_AMB(long int numeroBarras, long int numeroRamos, int nmed_AMB, int *caminho, DMED_COMPLEX *medidas_equivalentes, double *regua_x, double *regua_caminho, DRAM *ramos, GRAFO *grafo, double *hx_V)
{
    //__complex__ double **H_BC = NULL;
    double **H_T = NULL;

    H_T = aloca_matriz(6 * nmed_AMB, 6 * numeroBarras);

    __complex__ double ramo_Z;
    double R = 0;
    double X = 0;
    double G = 0;
    double B = 0;
    int i, j;

    

    for (i = 0; i < nmed_AMB; i++)
    {
        long int DE = medidas_equivalentes[i].DE;
        long int PARA = medidas_equivalentes[i].PARA;
        long int i_ramo = medidas_equivalentes[i].ramo;
        for (j = 0; j < numeroBarras; j++)
        {
            int x_atual = regua_x[j];

            //x_atual = ID da barra atual na regua
            //se a medida de potencia conter a barra (DE ou PARA) - calcular o Y para a fase

            //é medida de fluxo 
            if (medidas_equivalentes[i].tipo == 0 || medidas_equivalentes[i].tipo == 1){
                if ((DE == x_atual) || (PARA == x_atual))
                {
                    
                    switch (medidas_equivalentes[i].fases)
                    {
                    case 1:
                        //calcula Y
                        //FASE A
                        //IR/VR
                        G = calcula_G_AMB(0,i_ramo,ramos);
                        B = calcula_B_AMB(0,i_ramo,ramos);
                        H_T[3*i][3*j] = G;
                        //IR/VX
                        H_T[3*i][3*j + (3*numeroBarras)] = -B;
                        //IX/VR
                        H_T[3*i +(3*nmed_AMB)][3*j] = B;
                        //IX/VX
                        H_T[3*i +(3*nmed_AMB)][3*j + (3*numeroBarras)] = G;

                        break;
                    case 2:
                        //calcula Y
                        //FASE B
                        G = calcula_G_AMB(1,i_ramo,ramos);
                        B = calcula_B_AMB(1,i_ramo,ramos);
                        //IR/VR
                        H_T[3*i+1][3*j+1] = G;
                        //IR/VX
                        H_T[3*i+1][3*j+1 + (3*numeroBarras)] = -B;
                        //IX/VR
                        H_T[3*i+1 +(3*nmed_AMB)][3*j+1] = B;
                        //IX/VX
                        H_T[3*i+1 +(3*nmed_AMB)][3*j+1 + (3*numeroBarras)] = G;
                        break;
                    case 3:
                        //calcula Y
                        //FASE C
                        G = calcula_G_AMB(2,i_ramo,ramos);
                        B = calcula_B_AMB(2,i_ramo,ramos);
                        //IR/VR
                        H_T[3*i+2][3*j+2] = G;
                        //IR/VX
                        H_T[3*i+2][3*j+2 + (3*numeroBarras)] = -B;
                        //IX/VR
                        H_T[3*i+2 +(3*nmed_AMB)][3*j+2] = B;
                        //IX/VX
                        H_T[3*i+2 +(3*nmed_AMB)][3*j+2 + (3*numeroBarras)] = G;
                        break;
                    case 4:
                        //FASE AB
                        //IR/VR
                        G = calcula_G_AMB(0,i_ramo,ramos);
                        B = calcula_B_AMB(0,i_ramo,ramos);
                        H_T[3*i][3*j] = G;
                        //IR/VX
                        H_T[3*i][3*j + (3*numeroBarras)] = -B;
                        //IX/VR
                        H_T[3*i +(3*nmed_AMB)][3*j] = B;
                        //IX/VX
                        H_T[3*i +(3*nmed_AMB)][3*j + (3*numeroBarras)] = G;

                        //IR/VR
                        G = calcula_G_AMB(1,i_ramo,ramos);
                        B = calcula_B_AMB(1,i_ramo,ramos);
                        H_T[3*i+1][3*j+1] = G;
                        //IR/VX
                        H_T[3*i+1][3*j+1 + (3*numeroBarras)] = -B;
                        //IX/VR
                        H_T[3*i+1 +(3*nmed_AMB)][3*j+1] = B;
                        //IX/VX
                        H_T[3*i+1 +(3*nmed_AMB)][3*j+1 + (3*numeroBarras)] = G;
                        break;
                    case 5:
                        //FASE CA
                        //IR/VR
                        G = calcula_G_AMB(0,i_ramo,ramos);
                        B = calcula_B_AMB(0,i_ramo,ramos);
                        H_T[3*i][3*j] = G;
                        //IR/VX
                        H_T[3*i][3*j + (3*numeroBarras)] = -B;
                        //IX/VR
                        H_T[3*i +(3*nmed_AMB)][3*j] = B;
                        //IX/VX
                        H_T[3*i +(3*nmed_AMB)][3*j + (3*numeroBarras)] = G;

                        //IR/VR
                        G = calcula_G_AMB(2,i_ramo,ramos);
                        B = calcula_B_AMB(2,i_ramo,ramos);
                        H_T[3*i+2][3*j+2] = G;
                        //IR/VX
                        H_T[3*i+2][3*j+2 + (3*numeroBarras)] = -B;
                        //IX/VR
                        H_T[3*i+2 +(3*nmed_AMB)][3*j+2] = B;
                        //IX/VX
                        H_T[3*i+2 +(3*nmed_AMB)][3*j+2 + (3*numeroBarras)] = G;
                        break;
                    case 6:
                        //FASE BC
                        //IR/VR
                        G = calcula_G_AMB(1,i_ramo,ramos);
                        B = calcula_B_AMB(1,i_ramo,ramos);
                        H_T[3*i+1][3*j+1] = G;
                        //IR/VX
                        H_T[3*i+1][3*j+1 + (3*numeroBarras)] = -B;
                        //IX/VR
                        H_T[3*i+1 +(3*nmed_AMB)][3*j+1] = B;
                        //IX/VX
                        H_T[3*i+1 +(3*nmed_AMB)][3*j+1 + (3*numeroBarras)] = G;


                        //IR/VR
                        G = calcula_G_AMB(2,i_ramo,ramos);
                        B = calcula_B_AMB(2,i_ramo,ramos);
                        H_T[3*i+2][3*j+2] = G;
                        //IR/VX
                        H_T[3*i+2][3*j+2 + (3*numeroBarras)] = -B;
                        //IX/VR
                        H_T[3*i+2 +(3*nmed_AMB)][3*j+2] = B;
                        //IX/VX
                        H_T[3*i+2 +(3*nmed_AMB)][3*j+2 + (3*numeroBarras)] = G;
                        break;
                    case 7:
                        //FASE ABC
                        //IR/VR
                        G = calcula_G_AMB(0,i_ramo,ramos);
                        B = calcula_B_AMB(0,i_ramo,ramos);
                        H_T[3*i][3*j] = G;
                        //IR/VX
                        H_T[3*i][3*j + (3*numeroBarras)] = -B;
                        //IX/VR
                        H_T[3*i +(3*nmed_AMB)][3*j] = B;
                        //IX/VX
                        H_T[3*i +(3*nmed_AMB)][3*j + (3*numeroBarras)] = G;

                        //IR/VR
                        G = calcula_G_AMB(1,i_ramo,ramos);
                        B = calcula_B_AMB(1,i_ramo,ramos);
                        H_T[3*i+1][3*j+1] = G;
                        //IR/VX
                        H_T[3*i+1][3*j+1 + (3*numeroBarras)] = -B;
                        //IX/VR
                        H_T[3*i+1 +(3*nmed_AMB)][3*j+1] = B;
                        //IX/VX
                        H_T[3*i+1 +(3*nmed_AMB)][3*j+1 + (3*numeroBarras)] = G;

                        
                        //IR/VR
                        G = calcula_G_AMB(2,i_ramo,ramos);
                        B = calcula_B_AMB(2,i_ramo,ramos);
                        H_T[3*i+2][3*j+2] = G;
                        //IR/VX
                        H_T[3*i+2][3*j+2 + (3*numeroBarras)] = -B;
                        //IX/VR
                        H_T[3*i+2 +(3*nmed_AMB)][3*j+2] = B;
                        //IX/VX
                        H_T[3*i+2 +(3*nmed_AMB)][3*j+2 + (3*numeroBarras)] = G;
                        break;
                    }
                } 
            }
            //é medida de injecao
            if (medidas_equivalentes[i].tipo == 2 || medidas_equivalentes[i].tipo == 3){
                //TODO
            }
        }
    }

    return H_T;
}

double **monta_matriz_H_AMB_Tensao(DMED_COMPLEX *medidas_tensao, GRAFO *grafo, long int numeroBarras, int nmed_T){
    double **H_Tensao;
    H_Tensao = aloca_matriz(nmed_T, 6*numeroBarras);

    int i, j;
    for (i = 0; i< nmed_T; i++){
        int DE_med = medidas_tensao[i].DE;
        for (j = 0; j<numeroBarras; j++){
            //TODO
            int DE_barra = grafo[j].barra->ID;
            if (DE_med == DE_barra){
                switch (medidas_tensao[i].fases)
                {
                case 1:
                    /* code */
                    H_Tensao[i][3*j] = 1;
                    H_Tensao[i][3*j + 3*numeroBarras] = 1;
                    break;
                case 2:
                    /* code */
                    H_Tensao[i][3*j+1] = 1;
                    H_Tensao[i][3*j+1 + 3*numeroBarras] = 1;
                    break;
                case 3:
                    /* code */
                    H_Tensao[i][3*j+2] = 1;
                    H_Tensao[i][3*j+2 + 3*numeroBarras] = 1;
                    break;
                }
            }

        }
    }

    return H_Tensao;
}


void estimadorAMB(GRAFO *grafo, long int numeroRamos, long int numeroBarras, DMED *medidas, long int **numeroMedidas, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos, double Sbase, DBAR *barra)
{
    long int nmed, nvar;
    int i, j, k, r;
    double *z = NULL, **h = NULL, ***H = NULL, **W = NULL, *x = NULL, *regua = NULL, aux = 0;
    double tol = 0.000001;
    //__complex__ double *x_bc = NULL;

    long int nmed_AMB;
    long int nmed_T;

    int *caminho = NULL;
    int *nadj_proxbarra = NULL;
    caminho = malloc(numeroBarras * sizeof(int));
    nadj_proxbarra = malloc(numeroBarras * sizeof(int));
    busca_loop_grafo(grafo, numeroRamos, numeroBarras, caminho, nadj_proxbarra);

    printf("Estimador de Estado AMB...\n");
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

    nmed_AMB = conta_medidas_BC(medidas, nmed);
    // cria_B_Z_ramos(grafo, numeroRamos, ramos, Sbase);

    // x_bc = aloca_vetor(numeroRamos);
    // inicializa_vetor_estados_BC(x_bc, 3*numeroRamos);
    // inicializar vetor de variaveis de estado

    DMED_COMPLEX *medidas_complexas = NULL;
    medidas_complexas = (DMED_COMPLEX *)malloc((nmed_AMB) * sizeof(DMED_COMPLEX));

    double *x_bc = NULL;
    x_bc = monta_vetor_x_inicial_AMB(numeroBarras);

    double *delta_x_bc = NULL;
    delta_x_bc = aloca_vetor(6 * numeroBarras);
    // x_bc = (double *)malloc((6*numeroRamos) * sizeof(double));

    DMED_COMPLEX *medidas_equivalentes = NULL;
    medidas_equivalentes = (DMED_COMPLEX *)malloc((nmed_AMB) * sizeof(DMED_COMPLEX));


    double *z_AMB = NULL;
    z_AMB = aloca_vetor(6 * nmed_AMB + nmed_T);

    DMED_COMPLEX *medidas_tensao = NULL;
    nmed_T = conta_medidas_Tensao(medidas, nmed);
    medidas_tensao = (DMED_COMPLEX *)malloc((nmed_T) * sizeof(DMED_COMPLEX));

    // converte medidas tensao (magnitude) em complexas

    // z_eq = (__complex__ double *)malloc(3*nmed_BC * sizeof(__complex__ double));

    double *regua_x = NULL;
    double *regua_med = NULL;
    double *regua_med_inv = NULL;
    double **H_AMB = NULL;
    double **H_T = NULL;

    double *x_anterior = NULL;
    x_anterior = aloca_vetor(6 * numeroBarras);
    double *dif_x = NULL;
    dif_x = aloca_vetor(6 * numeroBarras);
    // vetor de estados: 1 para cada ramo e fase;
    regua_med = aloca_vetor(3 * nmed_AMB);
    regua_med_inv = aloca_vetor(3 * nmed_AMB);
    regua_x = aloca_vetor(3 * numeroBarras);
    H_AMB = aloca_matriz(6 * nmed_AMB, 6 * numeroBarras);
    H_T = aloca_matriz(nmed_T, 6 * numeroRamos);

    incializa_tensoes_grafo(grafo, numeroBarras, alimentadores, numeroAlimentadores);
    // printf("1\n");
    medidas_complexas = converte_medidas_para_complexo(medidas, nmed);


    //
    
    medidas_tensao = calcula_medida_tensao_complexa_AMB(medidas, nmed, grafo, numeroBarras);

    double *regua_V = NULL;
    regua_V = (double *)malloc(nmed_T * sizeof(double));
    // regua das medidas de tensao
    monta_regua_medidas_tensao(medidas_tensao, nmed_T, regua_V);

    // printf("2\n");
    medidas_equivalentes = divide_medidas_por_tensao(medidas_complexas, nmed_AMB, numeroBarras, grafo);

    printf("3\n");

    // regua vetor x
    //monta_regua_x(numeroRamos, regua_x, ramos);
    regua_x = monta_regua_x_AMB(grafo, numeroBarras);
    // printf("4\n");
    monta_regua_medidas(nmed_AMB, regua_med, regua_med_inv, medidas_equivalentes);

    double *regua_caminho = NULL;
    regua_caminho = aloca_vetor(numeroBarras);
    monta_regua_caminho(numeroRamos, numeroBarras, regua_caminho, caminho, grafo);
    // printf("5\n");

    //TODO: montar matriz H para o AMB;
    //Medidas de potência
    //H_BC = monta_matriz_H(numeroRamos, nmed_BC, regua_x, regua_med, regua_med_inv);
    //vetor h(x) das medidas de corrente
    double *hx_I = NULL;
    hx_I = (double*)malloc(6*nmed_AMB*sizeof(double));
    //vetor h(x) das medidas de tensão
    double *hx_V = NULL;
    hx_V = (double*)malloc(nmed_T*sizeof(double));
    // montar H das medidas de tensão
    //TODO: montar matriz H para o AMB;
    //Medidas de tensão
    H_T = monta_matriz_H_tensao(numeroBarras, numeroRamos, nmed_T, caminho, medidas_tensao, regua_x, regua_caminho, ramos, grafo, hx_V);


    // printf("6\n");
    int it = 0;
    int conv = 0;
    while (conv < 1)
    {

        // monta_z_complexa(medidas_equivalentes, z_eq, nmed_BC);

        monta_z_real_e_imag(medidas_equivalentes, z_AMB, nmed_AMB, medidas_tensao, nmed_T);

        int st = 0;

        //x_anterior = x_bc;
        
        // x_bc = resolve_linear_QR(H_BC, z_eq, numeroRamos, nmed_BC);
        //
        delta_x_bc = resolve_linear_QR_Tensao(H_AMB, H_T, z_AMB, numeroRamos, nmed_AMB, nmed_T, hx_I, hx_V);

        atualiza_vetor_x(x_bc, delta_x_bc, numeroRamos);

        //calcula_hx_corrente(H_BC, x_bc, hx_I, nmed_BC, numeroRamos);

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

        medidas_equivalentes = divide_medidas_por_tensao(medidas_complexas, nmed_AMB, numeroBarras, grafo);

        it++;
    }
}

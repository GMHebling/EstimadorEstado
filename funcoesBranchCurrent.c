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

int *montaRNP(ALIMENTADOR alimentadores)
{
    int *RNP;
    int k = 0;
    FILABARRAS *barraAtual;

    RNP = aloca_vetor_int(alimentadores.numeroNos + 1);
    barraAtual = &alimentadores.rnp[0];
    while (barraAtual != NULL)
    {
        RNP[k] = barraAtual->idNo;
        k++;
        barraAtual = barraAtual->prox;
    }

    return RNP;
}

void inicializa_vetor_estados_BC(double *x_bc, long int numeroRamos)
{
    // x_bc[0] = 1;
    for (int i = 0; i < numeroRamos; i++)
    {
        x_bc[i] = 0;
    }
}

DMED_COMPLEX *converte_medidas_para_complexo(DMED *medidas, long int numeroMedidas)
{
    DMED_COMPLEX *medidas_complexas = NULL;

    int cont_med_complex;
    cont_med_complex = conta_medidas_BC(medidas, numeroMedidas);

    printf("nmed_complex = %d\n", cont_med_complex);
    if (((medidas_complexas) = (DMED_COMPLEX *)malloc((cont_med_complex) * sizeof(DMED_COMPLEX))) == NULL)
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para as medidas complexas !!!!");
        exit(1);
    }

    int aux_contador = 0;
    int med_found = 0;
    for (int contador = 0; contador < numeroMedidas; contador++)
    {
        if (medidas[contador].tipo == 0 || medidas[contador].tipo == 2)
        {
            medidas_complexas[aux_contador].ligado = medidas[contador].ligado;
            medidas_complexas[aux_contador].tipo = medidas[contador].tipo;
            medidas_complexas[aux_contador].DE = medidas[contador].DE;
            medidas_complexas[aux_contador].PARA = medidas[contador].PARA;
            medidas_complexas[aux_contador].fases = medidas[contador].fases;
            medidas_complexas[aux_contador].id = medidas[contador].id;
            medidas_complexas[aux_contador].par = medidas[contador].par;

            medidas_complexas[aux_contador].sigma = medidas[contador].sigma;
            medidas_complexas[aux_contador].prec = medidas[contador].prec;

            double parte_real = medidas[contador].zmed;

            medidas_complexas[aux_contador].zmed = parte_real + 0 * I;

            for (int j = 0; j < numeroMedidas; j++)
            {
                if (medidas[j].tipo == 1 || medidas[j].tipo == 3)
                {
                    if (medidas[j].DE == medidas[contador].DE && medidas[j].PARA == medidas[contador].PARA && medidas[j].fases == medidas[contador].fases)
                    {
                        double parte_imag = medidas[j].zmed;

                        medidas_complexas[aux_contador].zmed = parte_real + parte_imag * I;
                        med_found += 1;
                        break;
                    }
                }
            }

            // printf("medida %d: %f + i*%f\n", aux_contador, creal(medidas_complexas[aux_contador].zmed)), cimag(medidas_complexas[aux_contador].zmed);
            aux_contador += 1;
        }
    }
    // printf("med_found: %d\n", med_found);

    return medidas_complexas;
}

int conta_medidas_BC(DMED *medidas, long int numeroMedidas)
{
    int cont_med_complex = 0;
    for (int contador = 0; contador < numeroMedidas; contador++)
    {
        if (medidas[contador].tipo == 0 || medidas[contador].tipo == 2)
        {
            cont_med_complex += 1;
        }
    }
    return cont_med_complex;
}

DMED_COMPLEX *divide_medidas_por_tensao(DMED_COMPLEX *medidas_complexas, long int numeroMedidas, long int numeroBarras, GRAFO *grafo)
{

    DMED_COMPLEX *medidas_div = NULL;
    if (((medidas_div) = (DMED_COMPLEX *)malloc((numeroMedidas) * sizeof(DMED_COMPLEX))) == NULL)
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para as medidas complexas !!!!");
        exit(1);
    }

    for (int cont = 0; cont < numeroMedidas; cont++)
    {
        medidas_div[cont] = medidas_complexas[cont];
        __complex__ double aux = medidas_complexas[cont].zmed;

        // printf("%f + i*%f\n", creal(medidas_div[cont].zmed), cimag(medidas_div[cont].zmed));

        long int id_de = medidas_div[cont].DE;

        for (int i = 0; i < numeroBarras; i++)
        {
            long int id_grafo_barra = grafo[i].barra->ID;

            if (id_grafo_barra == id_de)
            {
                switch (medidas_div[cont].fases)
                {
                case 1:
                    medidas_div[cont].zmed = aux / grafo[i].V[0];
                    break;
                case 2:
                    medidas_div[cont].zmed = aux / grafo[i].V[1];
                    break;
                case 3:
                    medidas_div[cont].zmed = aux / grafo[i].V[2];
                    break;
                case 4:
                    medidas_div[cont].zmed = aux / grafo[i].V[0];
                    break;
                case 5:
                    medidas_div[cont].zmed = aux / grafo[i].V[0];
                    break;
                case 6:
                    medidas_div[cont].zmed = aux / grafo[i].V[0];
                    break;
                case 7:
                    medidas_div[cont].zmed = aux / grafo[i].V[0];
                    break;
                }
            }
        }
        // printf("%f + i*%f\n", creal(medidas_div[cont].zmed), cimag(medidas_div[cont].zmed));
        // printf("\n");
        // TODO: Qual fase utilizar pra dividir a medida?
    }

    return medidas_div;
}

void monta_regua_x(long int numeroRamos, double *regua_x, DRAM *ramos)
{
    for (int nr = 0; nr < numeroRamos; nr++)
    {
        switch (ramos[nr].fases)
        {
        case 1:
            regua_x[3 * nr] = ramos[nr].DE + ramos[nr].PARA / 10000.0;
            break;
        case 2:
            regua_x[3 * nr + 1] = 2 * (ramos[nr].DE + ramos[nr].PARA / 10000.0);
            break;
        case 3:
            regua_x[3 * nr + 2] = 3 * (ramos[nr].DE + ramos[nr].PARA / 10000.0);
            break;
        case 4:
            regua_x[3 * nr] = ramos[nr].DE + ramos[nr].PARA / 10000.0;
            regua_x[3 * nr + 1] = 2 * (ramos[nr].DE + ((ramos[nr].PARA) / 10000.0));
            break;
        case 5:
            regua_x[3 * nr] = ramos[nr].DE + ramos[nr].PARA / 10000.0;
            regua_x[3 * nr + 2] = 3 * (ramos[nr].DE + ((ramos[nr].PARA) / 10000.0));
            break;
        case 6:
            regua_x[3 * nr + 1] = 2 * (ramos[nr].DE + ramos[nr].PARA / 10000.0);
            regua_x[3 * nr + 2] = 3 * (ramos[nr].DE + ((ramos[nr].PARA) / 10000.0));
            break;
        case 7:
            regua_x[3 * nr] = ramos[nr].DE + ramos[nr].PARA / 10000.0;
            regua_x[3 * nr + 1] = 2 * (ramos[nr].DE + ((ramos[nr].PARA) / 10000.0));
            regua_x[3 * nr + 2] = 3 * (ramos[nr].DE + ((ramos[nr].PARA) / 10000.0));
            break;
        }
        // printf("de: %ld, para: %ld, regua: %.4f\n", ramos[nr].DE, ramos[nr].PARA, regua_x[3*nr]);
    }
}

void monta_regua_medidas(long int nmed_BC, double *regua_med, double *regua_med_inv, DMED_COMPLEX *medidas_equivalentes)
{
    for (int nm = 0; nm < nmed_BC; nm++)
    {
        if (medidas_equivalentes[nm].tipo == 0 | medidas_equivalentes[nm].tipo == 1)
        {
            switch (medidas_equivalentes[nm].fases)
            {
            case 1:
                regua_med[3 * nm] = medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA / 10000.0;
                regua_med_inv[3 * nm] = -1 * (medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE / 10000.0);
                break;
            case 2:
                regua_med[3 * nm + 1] = 2.0 * (medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA / 10000.0);
                regua_med_inv[3 * nm + 1] = 2.0 * (-1 * (medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE / 10000.0));
                break;
            case 3:
                regua_med[3 * nm + 2] = 3.0 * (medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA / 10000.0);
                regua_med_inv[3 * nm + 2] = 3.0 * (-1 * (medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE / 10000.0));
                break;
            case 4:
                regua_med[3 * nm] = medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA / 10000.0;
                regua_med_inv[3 * nm] = -1 * (medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE / 10000.0);

                regua_med[3 * nm + 1] = 2 * (medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA / 10000.0);
                regua_med_inv[3 * nm + 1] = 2 * (-1 * (medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE / 10000.0));

                break;
            case 5:
                regua_med[3 * nm] = medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA / 10000.0;
                regua_med_inv[3 * nm] = -1 * (medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE / 10000.0);
                regua_med[3 * nm + 2] = 3 * (medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA / 10000.0);
                regua_med_inv[3 * nm + 2] = 3 * (-1 * (medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE / 10000.0));
                break;
            case 6:
                regua_med[3 * nm + 1] = 2 * (medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA / 10000.0);
                regua_med_inv[3 * nm + 1] = 2 * (-1 * (medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE / 10000.0));
                regua_med[3 * nm + 2] = 3 * (medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA / 10000.0);
                regua_med_inv[3 * nm + 2] = 3 * (-1 * (medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE / 10000.0));
                break;
            case 7:
                regua_med[3 * nm] = medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA / 10000.0;
                regua_med_inv[3 * nm] = -1 * (medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE / 10000.0);
                regua_med[3 * nm + 1] = 2 * (medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA / 10000.0);
                regua_med_inv[3 * nm + 1] = 2 * (-1 * (medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE / 10000.0));
                regua_med[3 * nm + 2] = 3 * (medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA / 10000.0);
                regua_med_inv[3 * nm + 2] = 3 * (-1 * (medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE / 10000.0));
                break;
            }
        }
        else
        {

            switch (medidas_equivalentes[nm].fases)
            {
            case 1:
                regua_med[3 * nm] = medidas_equivalentes[nm].DE;
                regua_med_inv[3 * nm] = medidas_equivalentes[nm].DE;
                break;
            case 2:
                regua_med[3 * nm + 1] = 2 * (medidas_equivalentes[nm].DE);
                regua_med_inv[3 * nm + 1] = 2 * (medidas_equivalentes[nm].DE);
                break;
            case 3:
                regua_med[3 * nm + 2] = 3 * (medidas_equivalentes[nm].DE);
                regua_med_inv[3 * nm + 2] = 3 * (medidas_equivalentes[nm].DE);
                break;
            case 4:
                regua_med[3 * nm] = medidas_equivalentes[nm].DE;
                regua_med_inv[3 * nm] = medidas_equivalentes[nm].DE;
                regua_med[3 * nm + 1] = 2 * (medidas_equivalentes[nm].DE);
                regua_med_inv[3 * nm + 1] = 2 * (medidas_equivalentes[nm].DE);
                break;
            case 5:
                regua_med[3 * nm] = medidas_equivalentes[nm].DE;
                regua_med_inv[3 * nm] = medidas_equivalentes[nm].DE;
                regua_med[3 * nm + 2] = 3 * (medidas_equivalentes[nm].DE);
                regua_med_inv[3 * nm + 2] = 3 * (medidas_equivalentes[nm].DE);
                break;
            case 6:
                regua_med[3 * nm + 1] = 2 * (medidas_equivalentes[nm].DE);
                regua_med_inv[3 * nm + 1] = 2 * (medidas_equivalentes[nm].DE);
                regua_med[3 * nm + 2] = 3 * (medidas_equivalentes[nm].DE);
                regua_med_inv[3 * nm + 2] = 3 * (medidas_equivalentes[nm].DE);
                break;
            case 7:
                regua_med[3 * nm] = medidas_equivalentes[nm].DE;
                regua_med_inv[3 * nm] = medidas_equivalentes[nm].DE;
                regua_med[3 * nm + 1] = 2 * (medidas_equivalentes[nm].DE);
                regua_med_inv[3 * nm + 1] = 2 * (medidas_equivalentes[nm].DE);
                regua_med[3 * nm + 2] = 3 * (medidas_equivalentes[nm].DE);
                regua_med_inv[3 * nm + 2] = 3 * (medidas_equivalentes[nm].DE);
                break;
            }
        }
        // printf("%.4f\n", regua_med_inv[nm]);
    }
}

double **monta_matriz_H(long int numeroRamos, long int nmed_BC, double *regua_x, double *regua_med, double *regua_med_inv)
{
    double **H_BC = NULL;
    int ctnz = 0;
    H_BC = aloca_matriz(3 * nmed_BC, 3 * numeroRamos);

    for (int nm = 0; nm < 3 * nmed_BC; nm++)
    {
        for (int nv = 0; nv < 3 * numeroRamos; nv++)
        {

            if (regua_med[nm] != 0.0 && regua_x[nv] != 0.0)
            {
                // ctnz += 1;

                // printf("x: %f ---- med: %f\n", regua_x[nv], regua_med[nm]);
                if (fabs(regua_med[nm] - regua_x[nv]) < EPS)
                {
                    // printf("1");
                    H_BC[nm][nv] = 1.0;
                    // printf("x: %f ---- med: %f ---- H: %f\n", regua_x[nv], regua_med[nm], H_BC[nm][nv]);
                }
                if (fabs(regua_med_inv[nm] + regua_x[nv]) < EPS)
                {

                    H_BC[nm][nv] = -1.0;
                    // printf("x: %f ---- med: %f ---- H: %f\n", regua_x[nv], regua_med[nm], H_BC[nm][nv]);
                }
                if (regua_x[nv] - regua_med[nm] > 0 && regua_x[nv] - regua_med[nm] < 1.0)
                {

                    H_BC[nm][nv] = 1.0;
                    // printf("x: %f ---- med: %f ---- H: %f\n", regua_x[nv], regua_med[nm], H_BC[nm][nv]);
                }
                if (fabs(((regua_x[nv] - (int)regua_x[nv]) * 10000.0) - regua_med[nm]) < EPS)
                {

                    H_BC[nm][nv] = -1.0;
                    // printf("x: %f ---- med: %f ---- H: %f\n", regua_x[nv], regua_med[nm], H_BC[nm][nv]);
                }
            }
        }
    }
    // printf("ctzn: %d\n", ctnz);

    return H_BC;
}
int verifica_ramo_caminho(double *regua_caminho, double *regua_x, int idx_ramo, long int numeroRamos, int limite_caminho)
{

    int i, j;
    double regua_vetor_x = 0.0;
    double regua_c = 0.0;
    regua_vetor_x = regua_x[3 * idx_ramo];
    for (j = 0; j < limite_caminho; j++)
    {
        regua_c = regua_caminho[j];
        if (regua_vetor_x != 0.0 && regua_c != 0.0)
        {
            if (fabs(regua_vetor_x - regua_c) < EPS)
            {
                return 1;
            }
        }
    }
    return 0;
}

void monta_regua_caminho(long int numeroRamos, long int numeroBarras, double *regua_caminho, int *caminho, GRAFO *grafo)
{
    // regua_caminho = aloca_vetor(numeroRamos);

    int i, j, cont = 0;
    int DE, PARA;
    for (i = 0; i < numeroBarras; i++)
    {
        DE = grafo[caminho[i]].barra->ID;
        for (j = i + 1; j < numeroBarras; j++)
        {
            PARA = grafo[caminho[j]].barra->ID;
            regua_caminho[cont] = DE + PARA / 10000.0;
            cont += 1;
        }
    }
}

int busca_id_caminho(int *caminho, int numeroBarraMedida, long int numeroBarras, GRAFO *grafo)
{
    int i;
    for (i = 0; i < numeroBarras; i++)
    {
        if (grafo[caminho[i]].barra->ID == numeroBarraMedida)
        {
            return i;
        }
    }
}

double calculo_theta_dv(GRAFO *grafo, int i_grafo, int i_fase, double *valor_hx_v)
{
    __complex__ double V_barra = grafo[i_grafo].V[i_fase];

    double V_real = creal(V_barra);
    double V_imag = cimag(V_barra);

    // double theta = atan(V_real / V_imag);
    double theta = atan(V_imag / V_real);

    valor_hx_v[0] = sqrt(V_real * V_real + V_imag * V_imag);
    return theta;
}

double _dVdIr(double R, double X, double theta)
{
    return (X * sin(theta) - R * cos(theta));
}
double _dVdIx(double R, double X, double theta)
{
    return (-1 * R * sin(theta) - X * cos(theta));
}

// TODO: trocar o tipo da funcao para void e adicionar as matrizes H_real e H_imag nos parametros
double **monta_matriz_H_tensao(long int numeroBarras, long int numeroRamos, int nmed_T, int *caminho, DMED_COMPLEX *medidas_tensao, double *regua_x, double *regua_caminho, DRAM *ramos, GRAFO *grafo, double *hx_V)
{
    //__complex__ double **H_BC = NULL;
    double **H_T = NULL;

    H_T = aloca_matriz(3 * nmed_T, 6 * numeroRamos);

    __complex__ double ramo_Z;
    double R = 0;
    double X = 0;

    double dVdIr = 0;
    double dVdIx = 0;
    int i, j;
    int i_grafo, i_fase;
    double theta;
    double *valor_hx_v = NULL;
    valor_hx_v = (double *)malloc(sizeof(double));

    for (i = 0; i < nmed_T; i++)
    {
        for (j = 0; j < numeroRamos; j++)
        {
            long int DE = medidas_tensao[i].DE;

            // barra 800 ou barra 0?
            int limite_caminho = busca_id_caminho(caminho, DE, numeroBarras, grafo);
            // todas as barras antes do limite devem entrar na jacobiana

            if (verifica_ramo_caminho(regua_caminho, regua_x, j, numeroRamos, limite_caminho) == 1)
            {
                // se o DE-PARA do ramo está antes do limite (barra da medida de tensao)
                switch (medidas_tensao[i].fases)
                {
                case 1:
                    i_fase = 0;
                    i_grafo = calcula_idx_ID(grafo, DE, numeroBarras);
                    theta = calculo_theta_dv(grafo, i_grafo, i_fase, valor_hx_v);
                    hx_V[i] = valor_hx_v[0];
                    if (ramos[j].Z == NULL)
                    {
                        R = 0.0;
                        X = 0.0;
                    }
                    else
                    {
                        ramo_Z = ramos[j].Z[0][0];
                        R = creal(ramo_Z);
                        X = cimag(ramo_Z);
                    }

                    dVdIr = _dVdIr(R, X, theta);
                    dVdIx = _dVdIx(R, X, theta);
                    H_T[i][3 * j] = dVdIr;
                    H_T[i][(3 * j) + 3 * numeroRamos] = dVdIx;
                    break;
                case 2:
                    i_fase = 1;
                    i_grafo = calcula_idx_ID(grafo, DE, numeroBarras);
                    theta = calculo_theta_dv(grafo, i_grafo, i_fase, hx_V);
                    hx_V[i] = valor_hx_v[0];
                    if (ramos[j].Z == NULL)
                    {
                        R = 0.0;
                        X = 0.0;
                    }
                    else
                    {
                        ramo_Z = ramos[j].Z[1][1];
                        R = creal(ramo_Z);
                        X = cimag(ramo_Z);
                    }

                    dVdIr = _dVdIr(R, X, theta);
                    dVdIx = _dVdIx(R, X, theta);
                    H_T[i][3 * j + 1] = dVdIr;
                    H_T[i][(3 * j + 1) + 3 * numeroRamos] = dVdIx;
                    break;
                case 3:
                    i_fase = 2;
                    i_grafo = calcula_idx_ID(grafo, DE, numeroBarras);
                    theta = calculo_theta_dv(grafo, i_grafo, i_fase, hx_V);
                    hx_V[i] = valor_hx_v[0];
                    if (ramos[j].Z == NULL)
                    {
                        R = 0.0;
                        X = 0.0;
                    }
                    else
                    {
                        ramo_Z = ramos[j].Z[2][2];
                        R = creal(ramo_Z);
                        X = cimag(ramo_Z);
                    }

                    dVdIr = _dVdIr(R, X, theta);
                    dVdIx = _dVdIx(R, X, theta);
                    H_T[i][3 * j + 2] = dVdIr;
                    H_T[i][(3 * j + 2) + 3 * numeroRamos] = dVdIx;
                    break;
                }
            }
        }
    }

    return H_T;
}

double *resolve_linear_QR(double **H_BC, double *z, long int numeroRamos, long int nmed_BC)
{
    cholmod_sparse *A = NULL;
    cholmod_triplet *T = NULL;
    cholmod_dense *b = NULL;
    cholmod_dense *X = NULL;
    cholmod_common Common, *c;
    c = &Common;
    cholmod_l_start(c);
    T = cholmod_l_allocate_triplet(3 * nmed_BC, 3 * numeroRamos, 3 * nmed_BC * 3 * numeroRamos, 0, CHOLMOD_COMPLEX, c);
    A = cholmod_l_allocate_sparse(3 * nmed_BC, 3 * numeroRamos, 3 * nmed_BC * 3 * numeroRamos, 0, 0, 0, CHOLMOD_COMPLEX, c);
    b = cholmod_l_allocate_dense(3 * nmed_BC, 1, 3 * nmed_BC, CHOLMOD_COMPLEX, c);
    X = cholmod_l_allocate_dense(3 * numeroRamos, 1, 3 * numeroRamos, CHOLMOD_COMPLEX, c);
    int index = 0;
    for (int i = 0; i < 3 * nmed_BC; i++)
    {
        for (int r = 0; r < 3 * numeroRamos; r++)
        {
            if (H_BC[i][r] != 0)
            {
                ((long int *)T->i)[index] = i;
                ((long int *)T->j)[index] = r;
                ((double *)T->x)[index] = H_BC[i][r];
                T->nnz += 1;
                index += 1;
            }
        }
    }
    for (int i = 0; i < 3 * nmed_BC; i++)
    {
        ((double *)b->x)[(2 * i)] = creal(z[i]);
        ((double *)b->x)[(2 * i) + 1] = cimag(z[i]);
    }
    A = cholmod_l_triplet_to_sparse(T, 3 * nmed_BC * 3 * numeroRamos, c);
    X = SuiteSparseQR_C_backslash(SPQR_ORDERING_BEST, SPQR_DEFAULT_TOL, A, b, c);
    double *ponto;
    ponto = aloca_vetor(6 * numeroRamos);
    ponto = (double *)X->x;
    // for (int ctz = 0; ctz < 10; ctz ++){
    //         printf("x[%d] = %f + i*%f\n", ctz, ponto[2*ctz], ponto[2*ctz+1]);
    // }
    return ponto;
}

double *resolve_linear_QR_Tensao(double **H_BC, double **H_T, double *z, long int numeroRamos, long int nmed_BC, int nmed_T, double *hx_I, double *hx_V)
{
    cholmod_sparse *A = NULL;
    cholmod_sparse *AT = NULL;
    cholmod_sparse *G = NULL;
    cholmod_triplet *T = NULL;
    cholmod_dense *b = NULL;
    cholmod_dense *bH = NULL;
    cholmod_dense *X = NULL;
    cholmod_common Common, *c;
    c = &Common;
    cholmod_l_start(c);
    T = cholmod_l_allocate_triplet(6 * nmed_BC + nmed_T, 6 * numeroRamos, (6 * nmed_BC + nmed_T) * 3 * numeroRamos, 0, CHOLMOD_REAL, c);

    // G = cholmod_l_allocate_sparse(6 * numeroRamos,6 * numeroRamos, (6 * numeroRamos* 6*numeroRamos), 0, 0, 1, CHOLMOD_REAL, c);
    A = cholmod_l_allocate_sparse(6 * nmed_BC + nmed_T, 6 * numeroRamos, (6 * nmed_BC + nmed_T) * 3 * numeroRamos, 0, 0, 0, CHOLMOD_REAL, c);
    // AT = cholmod_l_allocate_sparse( 6 * numeroRamos,6 * nmed_BC + nmed_T, (6 * nmed_BC + nmed_T) * 3 * numeroRamos, 0, 0, 0, CHOLMOD_REAL, c);
    b = cholmod_l_allocate_dense(6 * nmed_BC + nmed_T, 1, 6 * nmed_BC + nmed_T, CHOLMOD_REAL, c);
    bH = cholmod_l_allocate_dense(6 * numeroRamos, 1, 6 * numeroRamos, CHOLMOD_REAL, c);
    X = cholmod_l_allocate_dense(6 * numeroRamos, 1, 6 * numeroRamos, CHOLMOD_REAL, c);
    int index = 0;
    for (int i = 0; i < 3 * nmed_BC; i++)
    {
        for (int r = 0; r < 3 * numeroRamos; r++)
        {
            if (H_BC[i][r] != 0)
            {
                ((long int *)T->i)[index] = i;
                ((long int *)T->j)[index] = r;
                ((double *)T->x)[index] = H_BC[i][r];
                T->nnz += 1;
                index += 1;

                ((long int *)T->i)[index] = i + 3 * nmed_BC;
                ((long int *)T->j)[index] = r + 3 * numeroRamos;
                ((double *)T->x)[index] = H_BC[i][r];
                T->nnz += 1;
                index += 1;
            }
        }
    }

    for (int t = 0; t < nmed_T; t++)
    {
        for (int cv = 0; cv < 6 * numeroRamos; cv++)
        {
            if (H_T[t][cv] != 0)
            {
                ((long int *)T->i)[index] = t + 6 * nmed_BC;
                ((long int *)T->j)[index] = cv;
                ((double *)T->x)[index] = H_T[t][cv];
                T->nnz += 1;
                index += 1;
            }
        }
    }

    for (int i = 0; i < 6 * nmed_BC; i++)
    {
        ((double *)b->x)[(i)] = (z[i] - hx_I[i]);
    }

    for (int i = 6 * nmed_BC; i < 6 * nmed_BC + nmed_T; i++)
    {
        ((double *)b->x)[(i)] = (z[i] - hx_V[i]);
    }
    A = cholmod_l_triplet_to_sparse(T, (6 * nmed_BC + nmed_T) * 3 * numeroRamos, c);
    AT = cholmod_l_transpose(A, 1, c);
    G = cholmod_l_ssmult(AT, A, 0, 1, 0, c);
    double one[2] = {1, 0};
    double m1[2] = {0, 0};
    cholmod_l_sdmult(A, 1, one, m1, b, bH, c);
    X = SuiteSparseQR_C_backslash(SPQR_ORDERING_AMD, SPQR_DEFAULT_TOL, A, b, c);
    // X = SuiteSparseQR_C_backslash(SPQR_ORDERING_BEST, SPQR_DEFAULT_TOL, G, bH, c);
    double *ponto;
    ponto = aloca_vetor(6 * numeroRamos);
    ponto = (double *)X->x;
    // for (int ctz = 0; ctz < 10; ctz ++){
    //         printf("x[%d] = %f + i*%f\n", ctz, ponto[2*ctz], ponto[2*ctz+1]);
    // }
    return ponto;
}

void monta_z_complexa(DMED_COMPLEX *medidas_eq, __complex__ double *z, long int nmed_BC)
{
    for (int i = 0; i < nmed_BC; i++)
    {
        switch (medidas_eq[i].fases)
        {
        case 1:
            z[3 * i] = medidas_eq[i].zmed;
            break;
        case 2:
            z[3 * i + 1] = medidas_eq[i].zmed;
            break;
        case 3:
            z[3 * i + 2] = medidas_eq[i].zmed;
            break;
        case 4:
            z[3 * i] = medidas_eq[i].zmed;
            z[3 * i + 1] = medidas_eq[i].zmed;
            break;
        case 5:
            z[3 * i] = medidas_eq[i].zmed;
            z[3 * i + 2] = medidas_eq[i].zmed;
            break;
        case 6:
            z[3 * i + 1] = medidas_eq[i].zmed;
            z[3 * i + 2] = medidas_eq[i].zmed;
            break;
        case 7:
            z[3 * i] = medidas_eq[i].zmed;
            z[3 * i + 1] = medidas_eq[i].zmed;
            z[3 * i + 2] = medidas_eq[i].zmed;
            break;
        }
    }
}

void monta_z_real_e_imag(DMED_COMPLEX *medidas_eq, double *z, long int nmed_BC, DMED_COMPLEX *medidas_tensao, int nmed_T)
{
    for (int i = 0; i < nmed_BC; i++)
    {
        switch (medidas_eq[i].fases)
        {
        case 1:
            z[3 * i] = creal(medidas_eq[i].zmed);
            z[3 * i + 3 * nmed_BC] = cimag(medidas_eq[i].zmed);
            break;
        case 2:
            z[3 * i + 1] = creal(medidas_eq[i].zmed);
            z[3 * i + 1 + 3 * nmed_BC] = cimag(medidas_eq[i].zmed);
            break;
        case 3:
            z[3 * i + 2] = creal(medidas_eq[i].zmed);
            z[3 * i + 2 + 3 * nmed_BC] = cimag(medidas_eq[i].zmed);
            break;
        case 4:
            z[3 * i] = creal(medidas_eq[i].zmed);
            z[3 * i + 1] = creal(medidas_eq[i].zmed);

            z[3 * i + 3 * nmed_BC] = cimag(medidas_eq[i].zmed);
            z[3 * i + 1 + 3 * nmed_BC] = cimag(medidas_eq[i].zmed);
            break;
        case 5:
            z[3 * i] = creal(medidas_eq[i].zmed);
            z[3 * i + 2] = creal(medidas_eq[i].zmed);

            z[3 * i + 3 * nmed_BC] = cimag(medidas_eq[i].zmed);
            z[3 * i + 2 + 3 * nmed_BC] = cimag(medidas_eq[i].zmed);
            break;
        case 6:
            z[3 * i + 1] = creal(medidas_eq[i].zmed);
            z[3 * i + 2] = creal(medidas_eq[i].zmed);

            z[3 * i + 1 + 3 * nmed_BC] = cimag(medidas_eq[i].zmed);
            z[3 * i + 2 + 3 * nmed_BC] = cimag(medidas_eq[i].zmed);
            break;
        case 7:
            z[3 * i] = creal(medidas_eq[i].zmed);
            z[3 * i + 1] = creal(medidas_eq[i].zmed);
            z[3 * i + 2] = creal(medidas_eq[i].zmed);

            z[3 * i + 3 * nmed_BC] = cimag(medidas_eq[i].zmed);
            z[3 * i + 1 + 3 * nmed_BC] = cimag(medidas_eq[i].zmed);
            z[3 * i + 2 + 3 * nmed_BC] = cimag(medidas_eq[i].zmed);
            break;
        }
    }

    for (int j = 0; j < nmed_T; j++)
    {
        z[j + (6 * nmed_BC)] = creal(medidas_tensao[j].zmed);
    }
}

void incializa_tensoes_grafo(GRAFO *grafo, long int numeroBarras, ALIMENTADOR *alimentadores, long int numeroAlimentadores)
{
    int i, k, fase;
    BOOL visitado[numeroBarras];
    __complex__ double V0[3], **Yaux;

    Yaux = c_matAloca(3);

    // Flat start trifásico (Va = Vb = Vc = 1p.u.  Ta = 0  Tb = -120  Tc = 120) - com busca em profundidade para atualizar taps iniciais
    for (i = 0; i < numeroBarras; i++)
    {
        visitado[i] = false;
    }
    for (i = 0; i < numeroAlimentadores; i++)
    {
        // Tensão Inicial da subestação
        V0[0] = grafo[alimentadores[i].noRaiz].barra->Vinicial[0];
        V0[1] = grafo[alimentadores[i].noRaiz].barra->Vinicial[1];
        V0[2] = grafo[alimentadores[i].noRaiz].barra->Vinicial[2];

        FILABARRAS *barraAtual = &alimentadores[i].rnp[0];

        int de = barraAtual->idNo;
        grafo[de].V[0] = V0[0];
        grafo[de].V[1] = V0[1];
        grafo[de].V[2] = V0[2];

        while (barraAtual != NULL)
        {
            de = barraAtual->idNo;
            int n_adj = grafo[de].numeroAdjacentes;
            for (k = 0; k < n_adj; k++)
            {
                int para = grafo[de].adjacentes[k].idNo;
                if (visitado[para] == false)
                {
                    if ((grafo[de].adjacentes[k].tipo == 1))
                    { // Atualiza o V0 para trafo visto a ligação e tap
                        grafo[para].V[0] = grafo[de].V[0];
                        grafo[para].V[1] = grafo[de].V[1];
                        grafo[para].V[2] = grafo[de].V[2];

                        if ((grafo[de].adjacentes[k].ramo->trafo.lig_pri == 1) && (grafo[de].adjacentes[k].ramo->trafo.lig_sec == 2))
                        {
                            grafo[para].V[0] = cabs(grafo[de].V[0]) * (cos(-30 * PI / 180) + I * sin(-30 * PI / 180));
                            grafo[para].V[1] = cabs(grafo[de].V[1]) * (cos(-150 * PI / 180) + I * sin(-150 * PI / 180));
                            grafo[para].V[2] = cabs(grafo[de].V[2]) * (cos(90 * PI / 180) + I * sin(90 * PI / 180));
                        }
                        else if ((grafo[de].adjacentes[k].ramo->trafo.lig_pri == 3) && (grafo[de].adjacentes[k].ramo->trafo.lig_sec == 2))
                        {
                            grafo[para].V[0] = cabs(grafo[de].V[0]) * (cos(-30 * PI / 180) + I * sin(-30 * PI / 180));
                            grafo[para].V[1] = cabs(grafo[de].V[1]) * (cos(-150 * PI / 180) + I * sin(-150 * PI / 180));
                            grafo[para].V[2] = cabs(grafo[de].V[2]) * (cos(90 * PI / 180) + I * sin(90 * PI / 180));
                        }
                        else if ((grafo[de].adjacentes[k].ramo->trafo.lig_pri == 2) && (grafo[de].adjacentes[k].ramo->trafo.lig_sec == 1))
                        {
                            if (grafo[de].adjacentes[k].ramo->k == de)
                            {
                                grafo[para].V[0] = cabs(grafo[de].V[0]) * (cos(-30 * PI / 180) + I * sin(-30 * PI / 180));
                                grafo[para].V[1] = cabs(grafo[de].V[1]) * (cos(-150 * PI / 180) + I * sin(-150 * PI / 180));
                                grafo[para].V[2] = cabs(grafo[de].V[2]) * (cos(90 * PI / 180) + I * sin(90 * PI / 180));
                            }
                            else
                            {
                                grafo[para].V[0] = cabs(grafo[de].V[0]) * (cos(0 * PI / 180) + I * sin(0 * PI / 180));
                                grafo[para].V[1] = cabs(grafo[de].V[1]) * (cos(-120 * PI / 180) + I * sin(-120 * PI / 180));
                                grafo[para].V[2] = cabs(grafo[de].V[2]) * (cos(120 * PI / 180) + I * sin(120 * PI / 180));
                            }
                        }
                    }
                    else if (grafo[de].adjacentes[k].tipo == 2)
                    { // Para o caso de regulador de tensão
                        grafo[para].V[0] = grafo[de].V[0] * grafo[de].adjacentes[k].ramo->tap_pri[0] * grafo[de].adjacentes[k].ramo->tap_sec[0];
                        grafo[para].V[1] = grafo[de].V[1] * grafo[de].adjacentes[k].ramo->tap_pri[1] * grafo[de].adjacentes[k].ramo->tap_sec[1];
                        grafo[para].V[2] = grafo[de].V[2] * grafo[de].adjacentes[k].ramo->tap_pri[2] * grafo[de].adjacentes[k].ramo->tap_sec[2];
                    }
                    else
                    {
                        grafo[para].V[0] = grafo[de].V[0];
                        grafo[para].V[1] = grafo[de].V[1];
                        grafo[para].V[2] = grafo[de].V[2];
                    }
                }
            }
            visitado[de] = true;
            barraAtual = barraAtual->prox;
        }
    }
}

void atualiza_Rede_BC(GRAFO *grafo, long int numeroBarras, DBAR *barra, double *regua_x, long int numeroRamos, double *x_bc)
{
    int i, j, k, idMed, de, para, ramo, fase;
    __complex__ double *Saux, *Iaux;
    double aux_regua;
    double aux_regua_inv;
    BOOL visitado[numeroBarras];

    // Percorre o grafo atualizando o cálculo de h(x))
    for (i = 0; i < numeroBarras; i++)
    {
        // Percorre os ramos adjacentes
        for (k = 0; k < (grafo)[i].numeroAdjacentes; k++)
        {
            aux_regua = barra[(grafo)[i].idNo].ID + barra[(grafo)[i].adjacentes[k].idNo].ID / 10000.0;
            aux_regua_inv = -(barra[(grafo)[i].adjacentes[k].idNo].ID + barra[(grafo)[i].idNo].ID / 10000.0);

            for (j = 0; j < numeroRamos; j++)
            {
                // printf("\naux: %f, regua: %f", aux_regua, regua_x[3*j]);
                if (fabs(aux_regua - regua_x[3 * j]) < EPS)
                {
                    switch (grafo[i].adjacentes[k].ramo->fases)
                    {
                    case 1:
                        (grafo)[i].adjacentes[k].Cur[0] = x_bc[6 * j] + I * x_bc[6 * j + 1];
                        break;
                    case 2:
                        (grafo)[i].adjacentes[k].Cur[1] = x_bc[6 * j + 2] + I * x_bc[6 * j + 3];
                        break;
                    case 3:
                        (grafo)[i].adjacentes[k].Cur[2] = x_bc[6 * j + 4] + I * x_bc[6 * j + 5];
                        break;
                    case 4:
                        (grafo)[i].adjacentes[k].Cur[0] = x_bc[6 * j] + I * x_bc[6 * j + 1];
                        (grafo)[i].adjacentes[k].Cur[1] = x_bc[6 * j + 2] + I * x_bc[6 * j + 3];
                        break;
                    case 5:
                        (grafo)[i].adjacentes[k].Cur[0] = x_bc[6 * j] + I * x_bc[6 * j + 1];
                        (grafo)[i].adjacentes[k].Cur[2] = x_bc[6 * j + 4] + I * x_bc[6 * j + 5];
                        break;
                    case 6:
                        (grafo)[i].adjacentes[k].Cur[1] = x_bc[6 * j + 2] + I * x_bc[6 * j + 3];
                        (grafo)[i].adjacentes[k].Cur[2] = x_bc[6 * j + 4] + I * x_bc[6 * j + 5];
                        break;
                    case 7:
                        (grafo)[i].adjacentes[k].Cur[0] = x_bc[6 * j] + I * x_bc[6 * j + 1];
                        (grafo)[i].adjacentes[k].Cur[1] = x_bc[6 * j + 2] + I * x_bc[6 * j + 3];
                        (grafo)[i].adjacentes[k].Cur[2] = x_bc[6 * j + 4] + I * x_bc[6 * j + 5];
                        break;
                    }
                }
                if (fabs(regua_x[3 * j] + aux_regua_inv) < EPS)
                {

                    switch (grafo[i].adjacentes[k].ramo->fases)
                    {
                    case 1:
                        (grafo)[i].adjacentes[k].Cur[0] = -(x_bc[6 * j] + I * x_bc[6 * j + 1]);
                        break;
                    case 2:
                        (grafo)[i].adjacentes[k].Cur[1] = -(x_bc[6 * j + 2] + I * x_bc[6 * j + 3]);
                        break;
                    case 3:
                        (grafo)[i].adjacentes[k].Cur[2] = -(x_bc[6 * j + 4] + I * x_bc[6 * j + 5]);
                        break;
                    case 4:
                        (grafo)[i].adjacentes[k].Cur[0] = -(x_bc[6 * j] + I * x_bc[6 * j + 1]);
                        (grafo)[i].adjacentes[k].Cur[1] = -(x_bc[6 * j + 2] + I * x_bc[6 * j + 3]);
                        break;
                    case 5:
                        (grafo)[i].adjacentes[k].Cur[0] = -(x_bc[6 * j] + I * x_bc[6 * j + 1]);
                        (grafo)[i].adjacentes[k].Cur[2] = -(x_bc[6 * j + 4] + I * x_bc[6 * j + 5]);
                        break;
                    case 6:
                        (grafo)[i].adjacentes[k].Cur[1] = -(x_bc[6 * j + 2] + I * x_bc[6 * j + 3]);
                        (grafo)[i].adjacentes[k].Cur[2] = -(x_bc[6 * j + 4] + I * x_bc[6 * j + 5]);
                        break;
                    case 7:
                        (grafo)[i].adjacentes[k].Cur[0] = -(x_bc[6 * j] + I * x_bc[6 * j + 1]);
                        (grafo)[i].adjacentes[k].Cur[1] = -(x_bc[6 * j + 2] + I * x_bc[6 * j + 3]);
                        (grafo)[i].adjacentes[k].Cur[2] = -(x_bc[6 * j + 4] + I * x_bc[6 * j + 5]);
                        break;
                    }
                }
            }
        }
    }
    // printf("\ni: %d -> k: %d = %f + j*%f", 0, 0, creal(grafo[0].adjacentes[0].Cur[0]), cimag(grafo[0].adjacentes[0].Cur[0]));
}

void atualiza_Rede_BC_Tensao(GRAFO *grafo, long int numeroBarras, DBAR *barra, double *regua_x, long int numeroRamos, double *x_bc)
{
    int i, j, k, idMed, de, para, ramo, fase;
    __complex__ double *Saux, *Iaux;
    double aux_regua;
    double aux_regua_inv;
    BOOL visitado[numeroBarras];

    // Percorre o grafo atualizando o cálculo de h(x))
    for (i = 0; i < numeroBarras; i++)
    {
        // Percorre os ramos adjacentes
        for (k = 0; k < (grafo)[i].numeroAdjacentes; k++)
        {
            aux_regua = barra[(grafo)[i].idNo].ID + barra[(grafo)[i].adjacentes[k].idNo].ID / 10000.0;
            aux_regua_inv = -(barra[(grafo)[i].adjacentes[k].idNo].ID + barra[(grafo)[i].idNo].ID / 10000.0);

            for (j = 0; j < numeroRamos; j++)
            {
                // printf("\naux: %f, regua: %f", aux_regua, regua_x[3*j]);
                if (fabs(aux_regua - regua_x[3 * j]) < EPS)
                {
                    switch (grafo[i].adjacentes[k].ramo->fases)
                    {
                    case 1:
                        (grafo)[i].adjacentes[k].Cur[0] = x_bc[3 * j] + I * x_bc[3 * j + numeroRamos];
                        break;
                    case 2:
                        (grafo)[i].adjacentes[k].Cur[1] = x_bc[3 * j + 1] + I * x_bc[3 * j + 1 + numeroRamos];
                        break;
                    case 3:
                        (grafo)[i].adjacentes[k].Cur[2] = x_bc[3 * j + 2] + I * x_bc[3 * j + 2 + numeroRamos];
                        break;
                    case 4:
                        (grafo)[i].adjacentes[k].Cur[0] = x_bc[3 * j] + I * x_bc[3 * j + numeroRamos];
                        (grafo)[i].adjacentes[k].Cur[1] = x_bc[3 * j + 1] + I * x_bc[3 * j + 1 + numeroRamos];
                        break;
                    case 5:
                        (grafo)[i].adjacentes[k].Cur[0] = x_bc[3 * j] + I * x_bc[3 * j + 1 + numeroRamos];
                        (grafo)[i].adjacentes[k].Cur[2] = x_bc[3 * j + 2] + I * x_bc[3 * j + 2 + numeroRamos];
                        break;
                    case 6:
                        (grafo)[i].adjacentes[k].Cur[1] = x_bc[3 * j + 1] + I * x_bc[3 * j + 1 + numeroRamos];
                        (grafo)[i].adjacentes[k].Cur[2] = x_bc[3 * j + 2] + I * x_bc[3 * j + 2 + numeroRamos];
                        break;
                    case 7:
                        (grafo)[i].adjacentes[k].Cur[0] = x_bc[3 * j] + I * x_bc[3 * j + numeroRamos];
                        (grafo)[i].adjacentes[k].Cur[1] = x_bc[3 * j + 1] + I * x_bc[3 * j + 1 + numeroRamos];
                        (grafo)[i].adjacentes[k].Cur[2] = x_bc[3 * j + 2] + I * x_bc[3 * j + 2 + numeroRamos];
                        break;
                    }
                }
                if (fabs(regua_x[3 * j] + aux_regua_inv) < EPS)
                {

                    switch (grafo[i].adjacentes[k].ramo->fases)
                    {
                    case 1:
                        (grafo)[i].adjacentes[k].Cur[0] = -(x_bc[3 * j] + I * x_bc[3 * j + numeroRamos]);
                        break;
                    case 2:
                        (grafo)[i].adjacentes[k].Cur[1] = -(x_bc[3 * j + 1] + I * x_bc[3 * j + 1 + numeroRamos]);
                        break;
                    case 3:
                        (grafo)[i].adjacentes[k].Cur[2] = -(x_bc[3 * j + 2] + I * x_bc[3 * j + 2 + numeroRamos]);
                        break;
                    case 4:
                        (grafo)[i].adjacentes[k].Cur[0] = -(x_bc[3 * j] + I * x_bc[3 * j + numeroRamos]);
                        (grafo)[i].adjacentes[k].Cur[1] = -(x_bc[3 * j + 1] + I * x_bc[3 * j + 1 + numeroRamos]);
                        break;
                    case 5:
                        (grafo)[i].adjacentes[k].Cur[0] = -(x_bc[3 * j] + I * x_bc[3 * j + 1 + numeroRamos]);
                        (grafo)[i].adjacentes[k].Cur[2] = -(x_bc[3 * j + 2] + I * x_bc[3 * j + 2 + numeroRamos]);
                        break;
                    case 6:
                        (grafo)[i].adjacentes[k].Cur[1] = -(x_bc[3 * j + 1] + I * x_bc[3 * j + 1 + numeroRamos]);
                        (grafo)[i].adjacentes[k].Cur[2] = -(x_bc[3 * j + 2] + I * x_bc[3 * j + 2 + numeroRamos]);
                        break;
                    case 7:
                        (grafo)[i].adjacentes[k].Cur[0] = -(x_bc[3 * j] + I * x_bc[3 * j + numeroRamos]);
                        (grafo)[i].adjacentes[k].Cur[1] = -(x_bc[3 * j + 1] + I * x_bc[3 * j + 1 + numeroRamos]);
                        (grafo)[i].adjacentes[k].Cur[2] = -(x_bc[3 * j + 2] + I * x_bc[3 * j + 2 + numeroRamos]);
                        break;
                    }
                }
            }
        }
    }
    // printf("\ni: %d -> k: %d = %f + j*%f", 0, 0, creal(grafo[0].adjacentes[0].Cur[0]), cimag(grafo[0].adjacentes[0].Cur[0]));
}

void exportaEstado_BC(GRAFO *grafo, long int numeroBarras)
{
    int i, k, fase;
    FILE *arqout;

    arqout = fopen("state.txt", "w+");

    for (i = 0; i < numeroBarras; i++)
    {
        switch (grafo[i].fases)
        {
        case 1:
            fprintf(arqout, "%d\t%.15f\n", i, cabs(grafo[i].V[0]));
            break;
        case 2:
            fprintf(arqout, "%d\t%.15f\n", i, cabs(grafo[i].V[1]));
            break;
        case 3:
            fprintf(arqout, "%d\t%.15f\n", i, cabs(grafo[i].V[2]));
            break;
        case 4:
            fprintf(arqout, "%d\t%.15f\n", i, cabs(grafo[i].V[0]));
            fprintf(arqout, "%d\t%.15f\n", i, cabs(grafo[i].V[1]));
            break;
        case 5:
            fprintf(arqout, "%d\t%.15f\n", i, cabs(grafo[i].V[0]));
            fprintf(arqout, "%d\t%.15f\n", i, cabs(grafo[i].V[2]));
            break;
        case 6:
            fprintf(arqout, "%d\t%.15f\n", i, cabs(grafo[i].V[1]));
            fprintf(arqout, "%d\t%.15f\n", i, cabs(grafo[i].V[2]));
            break;
        case 7:
            fprintf(arqout, "%d\t%.15f\n", i, cabs(grafo[i].V[0]));
            fprintf(arqout, "%d\t%.15f\n", i, cabs(grafo[i].V[1]));
            fprintf(arqout, "%d\t%.15f\n", i, cabs(grafo[i].V[2]));
            break;
        }
    }

    for (i = 0; i < numeroBarras; i++)
    {
        switch (grafo[i].fases)
        {
        case 1:
            fprintf(arqout, "%d\t%.15f\n", i, carg(grafo[i].V[0]));
            break;
        case 2:
            fprintf(arqout, "%d\t%.15f\n", i, carg(grafo[i].V[1]));
            break;
        case 3:
            fprintf(arqout, "%d\t%.15f\n", i, carg(grafo[i].V[2]));
            break;
        case 4:
            fprintf(arqout, "%d\t%.15f\n", i, carg(grafo[i].V[0]));
            fprintf(arqout, "%d\t%.15f\n", i, carg(grafo[i].V[1]));
            break;
        case 5:
            fprintf(arqout, "%d\t%.15f\n", i, carg(grafo[i].V[0]));
            fprintf(arqout, "%d\t%.15f\n", i, carg(grafo[i].V[2]));
            break;
        case 6:
            fprintf(arqout, "%d\t%.15f\n", i, carg(grafo[i].V[1]));
            fprintf(arqout, "%d\t%.15f\n", i, carg(grafo[i].V[2]));
            break;
        case 7:
            fprintf(arqout, "%d\t%.15f\n", i, carg(grafo[i].V[0]));
            fprintf(arqout, "%d\t%.15f\n", i, carg(grafo[i].V[1]));
            fprintf(arqout, "%d\t%.15f\n", i, carg(grafo[i].V[2]));
            break;
        }
    }
    fclose(arqout);
}

void estimadorBC_RECT(GRAFO *grafo, long int numeroRamos, long int numeroBarras, DMED *medidas, long int **numeroMedidas, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos, double Sbase, DBAR *barra)
{
    long int nmed, nvar;
    int i, j, k, r;
    double *z = NULL, **h = NULL, ***H = NULL, **W = NULL, *x = NULL, *regua = NULL, aux = 0;
    double tol = 0.000001;
    //__complex__ double *x_bc = NULL;

    long int nmed_BC;

    printf("Estimador de Estado Branch Current em Coordenadas retangulares...\n");
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
    // printf("numero barras: %d\n", numeroBarras);
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

    // printf("nmed: %d\n", nmed);
    // printf("nvar: %d\n", nvar);
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
    // x_bc = (double *)malloc((6*numeroRamos) * sizeof(double));

    DMED_COMPLEX *medidas_equivalentes = NULL;
    medidas_equivalentes = (DMED_COMPLEX *)malloc((nmed_BC) * sizeof(DMED_COMPLEX));

    __complex__ double *z_eq = NULL;
    z_eq = c_vetAloca(3 * nmed_BC);
    // z_eq = (__complex__ double *)malloc(3*nmed_BC * sizeof(__complex__ double));

    double *regua_x = NULL;
    double *regua_med = NULL;
    double *regua_med_inv = NULL;
    double **H_BC = NULL;

    double *x_anterior = NULL;
    x_anterior = aloca_vetor(6 * numeroRamos);
    double *dif_x = NULL;
    dif_x = aloca_vetor(6 * numeroRamos);
    // vetor de estados: 1 para cada ramo e fase;
    regua_med = aloca_vetor(3 * nmed_BC);
    regua_med_inv = aloca_vetor(3 * nmed_BC);
    regua_x = aloca_vetor(3 * numeroRamos);
    H_BC = aloca_matriz(3 * nmed_BC, 3 * numeroRamos);

    incializa_tensoes_grafo(grafo, numeroBarras, alimentadores, numeroAlimentadores);
    printf("1\n");
    medidas_complexas = converte_medidas_para_complexo(medidas, nmed);
    printf("2\n");
    medidas_equivalentes = divide_medidas_por_tensao(medidas_complexas, nmed_BC, numeroBarras, grafo);

    printf("3\n");

    monta_regua_x(numeroRamos, regua_x, ramos);
    printf("4\n");
    monta_regua_medidas(nmed_BC, regua_med, regua_med_inv, medidas_equivalentes);
    printf("5\n");
    H_BC = monta_matriz_H(numeroRamos, nmed_BC, regua_x, regua_med, regua_med_inv);
    printf("6\n");
    int it = 0;
    int conv = 0;
    while (conv < 1)
    {

        monta_z_complexa(medidas_equivalentes, z_eq, nmed_BC);
        printf("7\n");

        // printf("\n");

        // for (int ctz = 0; ctz < 20; ctz++)
        // {
        //     //printf("z[%d] = %f + i*%f\n", ctz, creal(z_eq[ctz]), cimag(z_eq[ctz]));
        //     //printf("reguax : %f\n", regua_x[ctz]);
        // }
        // //printf("\n");

        // monta matriz Jacobiana
        // H_BC = monta_matriz_H(numeroRamos, nmed_BC, regua_x, regua_med, regua_med_inv);
        int st = 0;

        x_anterior = x_bc;
        x_bc = resolve_linear_QR(H_BC, z_eq, numeroRamos, nmed_BC);
        // for (int cx = 0; cx < 2; cx++){
        //     printf("xbc: %f\n", x_bc[cx]);
        //     printf("\n");
        // }

        for (int cx = 0; cx < 6 * numeroRamos; cx++)
        {
            dif_x[cx] = x_anterior[cx] - x_bc[cx];
        }
        double nfx;

        nfx = norma_inf(dif_x, 6 * numeroRamos);
        printf("\n\nIteracao:  %d \t|Dx|_inf =  %.17lf \t  \n", it, nfx);

        // mudar atualiza rede para receber complexo
        atualiza_Rede_BC(grafo, numeroBarras, barra, regua_x, numeroRamos, x_bc);

        // for (int nb = 0; nb < numeroBarras; nb++){
        //     for (int nj = 0; nj< grafo[nb].numeroAdjacentes; nj++){
        //         printf("\ni: %d -> k: %d = %f + j*%f", nb, nj, creal(grafo[nb].adjacentes[nj].Cur[0]), cimag(grafo[nb].adjacentes[nj].Cur[0]));
        //
        //     }
        // }

        for (int k = 0; k < numeroBarras; k++)
        {
            int ct = forward_sweep(&grafo[RNP[k]], grafo);
        }
        if (nfx<tol | it> 30)
        {
            conv = 10;

            // for (int nm = 0; nm < 3 * nmed_BC; nm++)
            // {
            //     for (int nv = 0; nv < 3 * numeroRamos; nv++)
            //     {

            //         if (regua_med[nm] != 0.0 && regua_x[nv] != 0.0)
            //         {
            //             if (fabs(regua_med[nm] - regua_x[nv]) < EPS)
            //             {
            //                 __complex__ double residuo = x_bc[2 * nv] + I * x_bc[2 * nv + 1] - z_eq[nm];
            //                 //printf("residuo: %f + i*%f\n", creal(residuo), cimag(residuo));
            //             }
            //         }
            //     }
            // }
        }

        medidas_equivalentes = divide_medidas_por_tensao(medidas_complexas, nmed_BC, numeroBarras, grafo);

        it++;
    }
}
void buscaProfundidadeLoop(GRAFO *grafo, int idNo, int barraAnterior, BOOL *visitado, int *caminho, int *contBarras, int *barraEntrada, int *nadj_proxbarra)
{
    int i;
    int proxBarra;

    int cont_loop = 0;
    visitado[idNo] = true;
    barraEntrada[idNo] = barraAnterior;
    for (i = 0; i < grafo[idNo].numeroAdjacentes; i++)
    {
        proxBarra = grafo[idNo].adjacentes[i].idNo;
        if ((proxBarra != barraEntrada[idNo]))
        {
            if (visitado[proxBarra] == true)
            {
                cont_loop += 1;
                // desc_loop[(2*n_loops[0])] = proxBarra;
                // desc_loop[(2*n_loops[0] + 1)] = idNo;
                // n_loops[0] += 1;
                // printf("idNo: %d | proxBarra? %d\n", grafo[idNo].barra->ID, grafo[proxBarra].barra->ID);
                // primeiro adjacente pode ser diferente da barra imediatamente antes
                // contar os loops e salvar na estrutura
                continue;
            }
            else
            {
                barraAnterior = idNo;
                contBarras[0] += 1;
                caminho[contBarras[0]] = proxBarra;
                nadj_proxbarra[contBarras[0]] = grafo[proxBarra].numeroAdjacentes;
                // printf("%d\n", caminho[contBarras[0]]);
                buscaProfundidadeLoop(grafo, proxBarra, barraAnterior, visitado, caminho, contBarras, barraEntrada, nadj_proxbarra);
            }
        }
        else
        {
            continue;
        }
    }
}

void busca_loop_grafo(GRAFO *grafo, long int numeroRamos, long int numeroBarras, int *caminho, int *nadj_proxbarra)
{
    // funcao recursiva
    // usar buscaprofundidade
    BOOL *visitado = NULL;
    int *barraEntrada = NULL;
    // int *caminho = NULL;
    // int *nadj_proxbarra = NULL;
    int i;
    int noRaiz = 0;
    int barraAnterior = 0;
    int *contBarras;
    // int *n_loops = NULL;
    // int *desc_loop = NULL;

    contBarras = malloc(1 * sizeof(int));
    visitado = malloc(numeroBarras * sizeof(BOOL));
    // caminho = malloc(numeroBarras*sizeof(int));
    // nadj_proxbarra = malloc(numeroBarras*sizeof(int));
    barraEntrada = malloc(numeroBarras * sizeof(int));

    contBarras[0] = 0;
    caminho[0] = noRaiz;
    nadj_proxbarra[0] = grafo[noRaiz].numeroAdjacentes;
    barraEntrada[0] = noRaiz;

    for (i = 0; i < numeroBarras; i++)
    {
        visitado[i] = false;
    }

    buscaProfundidadeLoop(grafo, noRaiz, barraAnterior, visitado, caminho, contBarras, barraEntrada, nadj_proxbarra);
    // printf("N_loops: %d\n", n_loops[0]);
    printf("\n");
    // for (i = 0; i < numeroBarras; i++){
    //     printf(" id: %d, n_adj: %d ", caminho[i], nadj_proxbarra[i]);
    // }
    printf("\n");
    printf("busca Completa\n");

    free(contBarras);
    // free(caminho);
    free(visitado);
    free(barraEntrada);

    return;
}

int conta_medidas_Tensao(DMED *medidas, long int numeroMedidas)
{
    int nmed_tensao = 0;
    for (int i = 0; i < numeroMedidas; i++)
    {
        if (medidas[i].tipo == 4)
        {
            nmed_tensao += 1;
        }
    }
    return nmed_tensao;
}

// retorna i do grafo a partir do ID da barra
int calcula_idx_ID(GRAFO *grafo, int ID, int numeroBarras)
{
    for (int i = 0; i < numeroBarras; i++)
    {
        if (grafo[i].barra->ID == ID)
        {
            return i;
        }
    }
}

DMED_COMPLEX *calcula_medida_tensao_complexa(DMED *medidas, long int numeroMedidas, GRAFO *grafo, int numeroBarras)
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

            // __complex__ double v_fase;
            // v_fase = grafo[idx_barra_DE].V[f_idx];

            // double v_real = creal(v_fase);
            // double v_imag = cimag(v_fase);

            // double v_ang = atan(v_imag/v_real);

            // // fase da medida
            // // tensao da barra na fase

            // double parte_real_z = magnitude_tensao * cos(v_ang);
            // double parte_imag_z = magnitude_tensao * sin(v_ang);

            medida_tensao[aux_contador].zmed = magnitude_tensao + 0 * I;

            aux_contador += 1;
        }
    }
    return medida_tensao;
}

void monta_regua_medidas_tensao(DMED_COMPLEX *medidas_tensao, int nmed_T, double *regua_V)
{
    // regua = DE + 0.FASE
    for (int i = 0; i < nmed_T; i++)
    {

        switch (medidas_tensao[i].fases)
        {
        case 1:
            regua_V[i] = (double)medidas_tensao[i].DE + 0.1;
            break;
        case 2:
            regua_V[i] = (double)medidas_tensao[i].DE + 0.2;
            break;
        case 3:
            regua_V[i] = (double)medidas_tensao[i].DE + 0.3;
            break;
        }
    }
}

void calcula_hx_corrente(double **H_BC, double *x, double *hx, int nmed_BC, int numeroRamos)
{
    int i, j, k;
    double aux_r = 0;
    double aux_x = 0;
    for (i = 0; i < 3 * nmed_BC; i++)
    {
        aux_r = 0;
        aux_x = 0;
        for (j = 0; j < 3 * numeroRamos; j++)
        {
            aux_r += H_BC[i][j] * x[j];
            aux_x += H_BC[i][j] * x[j + 3 * numeroRamos];
        }
        hx[i] = aux_r;
        hx[i + 3 * nmed_BC] = aux_x;
    }
}

void atualiza_vetor_x(double *x, double *dx, long int numeroRamos)
{
    for (int i = 0; i < 6 * numeroRamos; i++)
    {
        x[i] = x[i] + dx[i];
    }
}

void estimadorBC_RECT_Malhado(GRAFO *grafo, long int numeroRamos, long int numeroBarras, DMED *medidas, long int **numeroMedidas, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos, double Sbase, DBAR *barra)
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
    // vetor h(x) das medidas de corrente
    double *hx_I = NULL;
    hx_I = (double *)malloc(6 * nmed_BC * sizeof(double));
    // vetor h(x) das medidas de tensão
    double *hx_V = NULL;
    hx_V = (double *)malloc(nmed_T * sizeof(double));
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

        // x_anterior = x_bc;

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

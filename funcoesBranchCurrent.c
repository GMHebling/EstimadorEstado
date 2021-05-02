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

void inicializa_vetor_estados_BC(double *x_bc, long int numeroRamos){
    //x_bc[0] = 1;
    for (int i = 0; i<numeroRamos; i++){
        x_bc[i] = 0;
    }
}


DMED_COMPLEX *converte_medidas_para_complexo(DMED *medidas, long int numeroMedidas){
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
    for (int contador = 0; contador < numeroMedidas; contador++){
        if (medidas[contador].tipo == 0 || medidas[contador].tipo == 2) {
            medidas_complexas[aux_contador].ligado = medidas[contador].ligado;
            medidas_complexas[aux_contador].tipo = medidas[contador].tipo; 
            medidas_complexas[aux_contador].DE = medidas[contador].DE;
            medidas_complexas[aux_contador].PARA = medidas[contador].PARA;
            medidas_complexas[aux_contador].fases = medidas[contador].fases;
            medidas_complexas[aux_contador].id = medidas[contador].id;
            medidas_complexas[aux_contador].par = medidas[contador].par;
            
            medidas_complexas[aux_contador].sigma = medidas[contador].sigma;
            medidas_complexas[aux_contador].prec = medidas[contador].prec;

            double parte_real =  medidas[contador].zmed;

            medidas_complexas[aux_contador].zmed = parte_real + 0*I;

            
            for (int j = 0; j < numeroMedidas; j++){
                if (medidas[j].tipo == 1 || medidas[j].tipo == 3){
                    if (medidas[j].DE == medidas[contador].DE && medidas[j].PARA == medidas[contador].PARA && medidas[j].fases == medidas[contador].fases){
                        double parte_imag = medidas[j].zmed;

                        medidas_complexas[aux_contador].zmed = parte_real + parte_imag*I;
                        med_found += 1;
                        break;
                    }
                    
                }
            }
            
            //printf("medida %d: %f + i*%f\n", aux_contador, creal(medidas_complexas[aux_contador].zmed)), cimag(medidas_complexas[aux_contador].zmed);            
            aux_contador += 1;
        }
    }
    //printf("med_found: %d\n", med_found);

    return medidas_complexas;
}


int conta_medidas_BC(DMED *medidas, long int numeroMedidas){
    int cont_med_complex = 0;
    for (int contador = 0; contador < numeroMedidas; contador++){
        
        if (medidas[contador].tipo == 0 || medidas[contador].tipo == 2) {
             cont_med_complex += 1;
        }
    }
    return cont_med_complex;

}


DMED_COMPLEX *divide_medidas_por_tensao(DMED_COMPLEX *medidas_complexas, long int numeroMedidas, long int numeroBarras, GRAFO *grafo){

    DMED_COMPLEX *medidas_div = NULL;
    if (((medidas_div) = (DMED_COMPLEX *)malloc((numeroMedidas) * sizeof(DMED_COMPLEX))) == NULL)
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para as medidas complexas !!!!");
        exit(1);
    }

    for (int cont = 0; cont < numeroMedidas; cont++){
        medidas_div[cont] = medidas_complexas[cont];
        __complex__ double aux = medidas_complexas[cont].zmed;
        
        long int id_de = medidas_div[cont].DE;
        

        for (int i = 0; i < numeroBarras; i++){
            long int id_grafo_barra = grafo[i].barra->ID;
            
            if (id_grafo_barra == id_de){
                medidas_div[cont].zmed = aux/grafo[i].V[0];
                
                break;
            }
        }
        //TODO: Qual fase utilizar pra dividir a medida?
    }

    return medidas_div;
}


void monta_regua_x(long int numeroRamos, double *regua_x, DRAM *ramos){
    for (int nr = 0; nr < numeroRamos; nr++){
        regua_x[3*nr] = ramos[nr].DE + ramos[nr].PARA/10000.0;
        regua_x[3*nr+1] = ramos[nr].DE + ((ramos[nr].PARA)/10000.0);
        regua_x[3*nr+2] = ramos[nr].DE + ((ramos[nr].PARA)/10000.0);
        //printf("de: %ld, para: %ld, regua: %.4f\n", ramos[nr].DE, ramos[nr].PARA, regua_x[3*nr]);
        
    }
}

void monta_regua_medidas(long int nmed_BC, double *regua_med, double *regua_med_inv, DMED_COMPLEX *medidas_equivalentes){
    for (int nm = 0; nm < nmed_BC; nm++){
        if (medidas_equivalentes[nm].tipo == 0){
            regua_med[nm] = medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA/10000.0;
            regua_med_inv[nm] = -1*(medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE/10000.0);
        }
        else {
            regua_med[nm] = medidas_equivalentes[nm].DE;
            regua_med_inv[nm] = medidas_equivalentes[nm].DE;
        }
        //printf("%.4f\n", regua_med_inv[nm]);
    }
}

void estimadorBC_RECT(GRAFO *grafo, long int numeroRamos, long int numeroBarras, DMED *medidas, long int **numeroMedidas, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos, double Sbase)
{
    long int nmed, nvar;
    int i, j, k, r;
    double *z = NULL, **h = NULL, ***H = NULL, **W = NULL, *x = NULL, *regua = NULL, aux = 0;

    double *x_bc = NULL;

    DMED_COMPLEX *medidas_complexas = NULL;
    long int nmed_BC;
    

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

    
    x_bc = aloca_vetor(numeroRamos);
    inicializa_vetor_estados_BC(x_bc, 3*numeroRamos);
    //inicializar vetor de variaveis de estado

    nmed_BC = conta_medidas_BC(medidas, nmed);
    medidas_complexas = converte_medidas_para_complexo(medidas, nmed);

    DMED_COMPLEX *medidas_equivalentes = NULL;
    medidas_equivalentes = divide_medidas_por_tensao(medidas_complexas, nmed_BC, numeroBarras, grafo);


    double *regua_x = NULL;
    //vetor de estados: 1 para cada ramo e fase;
    regua_x = aloca_vetor(3*numeroRamos);
    monta_regua_x(numeroRamos, regua_x, ramos);
    
    double *regua_med = NULL;
    double *regua_med_inv = NULL;
    regua_med = aloca_vetor(nmed_BC);
    regua_med_inv = aloca_vetor(nmed_BC);

    monta_regua_medidas(nmed_BC, regua_med, regua_med_inv, medidas_equivalentes);

    double **H_BC = NULL;
    H_BC = aloca_matriz(nmed_BC, 3*numeroRamos);

    FILE *matrizH = NULL;
    matrizH = fopen("matrizH.txt", "w+");

    //monta Jacobiana para medidas de fluxo
    for (int nm = 0; nm < nmed_BC; nm++){
        for (int nv = 0; nv < 3*numeroRamos; nv++){
            if (regua_med[nm]/regua_x[nv] == 1.0){
                H_BC[nm][nv] = 1.0;                
            }
            if (regua_med_inv[nm]/regua_x[nv] == -1.0){
                H_BC[nm][nv] = -1.0;                
            }
            if (regua_x[nv] - regua_med[nm] > 0 && regua_x[nv] - regua_med[nm] < 1.0){
                H_BC[nm][nv] = 1.0; 
            }
            if (((regua_x[nv] - (int)regua_x[nv]) * 10000.0) - regua_med[nm] == 0.0){
                H_BC[nm][nv] = -1.0;
            }
            if (H_BC[nm][nv] != 0){
                fprintf(matrizH, "%d,%d,%f\n", nm, nv, H_BC[nm][nv]);
            }
            
        }
    }
    fclose(matrizH);

    
    //Todo: estimador WLS para as correntes

    //Todo: atualizar as correntes e fazer a varredura forward
    //atualizaRede
    //grafo[i].adjacentes[k].Cur[0] = corrente_estimada[i,k][0]
    //forward_sweep()
    

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
    int conv = 1;//otimiza_Gauss_NewtonQR(z, h, H, grafo, numeroBarras, ramos, medidas, nvar, nmed, regua, x, tol, ref_1, ref_2);

    clock_t t1 = clock();
    double tempoWLS = (double)(t1 - tIni) / CLOCKS_PER_SEC;
    printf("\nEstimação WLS: %lf", tempoWLS);

    exportaCasoReferencia(grafo, numeroBarras, Sbase);
    //exportaEstado(grafo, regua, nvar);
    //imprimeEstado(grafo, numeroBarras);

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

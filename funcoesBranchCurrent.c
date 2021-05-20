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
        switch (ramos[nr].fases){
            case 1:
                regua_x[3*nr] = ramos[nr].DE + ramos[nr].PARA/10000.0;
                break;
            case 2:
                regua_x[3*nr+1] = 2*(ramos[nr].DE + ramos[nr].PARA/10000.0);
                break;
            case 3:
                regua_x[3*nr+2] = 3*(ramos[nr].DE + ramos[nr].PARA/10000.0);
                break;
            case 4:
                regua_x[3*nr] = ramos[nr].DE + ramos[nr].PARA/10000.0;
                regua_x[3*nr+1] = 2*(ramos[nr].DE + ((ramos[nr].PARA)/10000.0));
                break;    
            case 5:
                regua_x[3*nr] = ramos[nr].DE + ramos[nr].PARA/10000.0;
                regua_x[3*nr+2] = 3*(ramos[nr].DE + ((ramos[nr].PARA)/10000.0));
                break;
            case 6:
                regua_x[3*nr+1] = 2*(ramos[nr].DE + ramos[nr].PARA/10000.0);
                regua_x[3*nr+2] = 3*(ramos[nr].DE + ((ramos[nr].PARA)/10000.0));
                break;
            case 7:
                regua_x[3*nr] = ramos[nr].DE + ramos[nr].PARA/10000.0;
                regua_x[3*nr+1] = 2*(ramos[nr].DE + ((ramos[nr].PARA)/10000.0));
                regua_x[3*nr+2] = 3*(ramos[nr].DE + ((ramos[nr].PARA)/10000.0));
                break;
        }
        //printf("de: %ld, para: %ld, regua: %.4f\n", ramos[nr].DE, ramos[nr].PARA, regua_x[3*nr]);
        
    }
}

void monta_regua_medidas(long int nmed_BC, double *regua_med, double *regua_med_inv, DMED_COMPLEX *medidas_equivalentes){
    for (int nm = 0; nm < nmed_BC; nm++){
        if (medidas_equivalentes[nm].tipo == 0){
            switch (medidas_equivalentes[nm].fases)
            {
            case 1:
                regua_med[3*nm] = medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA/10000.0;
                regua_med_inv[3*nm] = -1*(medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE/10000.0);
                break;
            case 2:
                regua_med[3*nm + 1] = 2.0*(medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA/10000.0);
                regua_med_inv[3*nm + 1] = 2.0*(-1*(medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE/10000.0));
                break;
            case 3:
                regua_med[3*nm + 2] = 3.0*(medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA/10000.0);
                regua_med_inv[3*nm + 2] = 3.0*(-1*(medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE/10000.0));
                break;
            case 4:
                regua_med[3*nm] = medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA/10000.0;
                regua_med_inv[3*nm] = -1*(medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE/10000.0);
                regua_med[3*nm + 1] = 2*(medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA/10000.0);
                regua_med_inv[3*nm + 1] = 2*(-1*(medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE/10000.0));
                break;
            case 5:
                regua_med[3*nm] = medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA/10000.0;
                regua_med_inv[3*nm] = -1*(medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE/10000.0);
                regua_med[3*nm + 2] = 3*(medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA/10000.0);
                regua_med_inv[3*nm + 2] = 3*(-1*(medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE/10000.0));
                break;
            case 6:
                regua_med[3*nm + 1] = 2*(medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA/10000.0);
                regua_med_inv[3*nm + 1] = 2*(-1*(medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE/10000.0));
                regua_med[3*nm + 2] = 3*(medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA/10000.0);
                regua_med_inv[3*nm + 2] = 3*(-1*(medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE/10000.0));
                break;
            case 7:
                regua_med[3*nm] = medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA/10000.0;
                regua_med_inv[3*nm] = -1*(medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE/10000.0);
                regua_med[3*nm + 1] = 2*(medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA/10000.0);
                regua_med_inv[3*nm + 1] = 2*(-1*(medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE/10000.0));
                regua_med[3*nm + 2] = 3*(medidas_equivalentes[nm].DE + medidas_equivalentes[nm].PARA/10000.0);
                regua_med_inv[3*nm + 2] = 3*(-1*(medidas_equivalentes[nm].PARA + medidas_equivalentes[nm].DE/10000.0));
                break;
            }
            
        }
        else {
            
            switch (medidas_equivalentes[nm].fases)
            {
            case 1:
                regua_med[3*nm] = medidas_equivalentes[nm].DE;
                regua_med_inv[3*nm] = medidas_equivalentes[nm].DE;
                break;
            case 2:
                regua_med[3*nm + 1] = 2*(medidas_equivalentes[nm].DE);
                regua_med_inv[3*nm + 1] = 2*(medidas_equivalentes[nm].DE);
                break;
            case 3:
                regua_med[3*nm + 2] = 3*(medidas_equivalentes[nm].DE);
                regua_med_inv[3*nm + 2] = 3*(medidas_equivalentes[nm].DE);
                break;
            case 4:
                regua_med[3*nm] = medidas_equivalentes[nm].DE;
                regua_med_inv[3*nm] = medidas_equivalentes[nm].DE;
                regua_med[3*nm + 1] = 2*(medidas_equivalentes[nm].DE);
                regua_med_inv[3*nm + 1] = 2*(medidas_equivalentes[nm].DE);
                break;
            case 5:
                regua_med[3*nm] = medidas_equivalentes[nm].DE;
                regua_med_inv[3*nm] = medidas_equivalentes[nm].DE;
                regua_med[3*nm + 2] = 3*(medidas_equivalentes[nm].DE);
                regua_med_inv[3*nm + 2] = 3*(medidas_equivalentes[nm].DE);
                break;
            case 6:
                regua_med[3*nm + 1] = 2*(medidas_equivalentes[nm].DE);
                regua_med_inv[3*nm + 1] = 2*(medidas_equivalentes[nm].DE);
                regua_med[3*nm + 2] = 3*(medidas_equivalentes[nm].DE);
                regua_med_inv[3*nm + 2] = 3*(medidas_equivalentes[nm].DE);
                break;
            case 7:
                regua_med[3*nm] = medidas_equivalentes[nm].DE;
                regua_med_inv[3*nm] = medidas_equivalentes[nm].DE;
                regua_med[3*nm + 1] = 2*(medidas_equivalentes[nm].DE);
                regua_med_inv[3*nm + 1] = 2*(medidas_equivalentes[nm].DE);
                regua_med[3*nm + 2] = 3*(medidas_equivalentes[nm].DE);
                regua_med_inv[3*nm + 2] = 3*(medidas_equivalentes[nm].DE);
                break;
            }
        }
        //printf("%.4f\n", regua_med_inv[nm]);
    }
}

double  **monta_matriz_H(long int numeroRamos, long int nmed_BC, double *regua_x, double *regua_med, double *regua_med_inv){
    double **H_BC = NULL;
    int ctnz = 0;
    H_BC = aloca_matriz(3*nmed_BC, 3*numeroRamos);
    
    for (int nm = 0; nm < 3*nmed_BC; nm++){
        for (int nv = 0; nv < 3*numeroRamos; nv++){
            
            if (regua_med[nm] != 0.0 && regua_x[nv] != 0.0){
                ctnz += 1;
                if (fabs(regua_med[nm] - regua_x[nv]) < EPS){
                    
                    H_BC[nm][nv] = 1.0;                
                }
                if (fabs(regua_med_inv[nm] + regua_x[nv]) < EPS){
                    
                    H_BC[nm][nv] = -1.0;                
                }
                if (regua_x[nv] - regua_med[nm] > 0 && regua_x[nv] - regua_med[nm] < 1.0){
                    
                    H_BC[nm][nv] = 1.0; 
                }
                if (fabs(((regua_x[nv] - (int)regua_x[nv]) * 10000.0) - regua_med[nm]) < EPS){
                    
                    H_BC[nm][nv] = -1.0;
                }
            }
        }
    }
    //printf("ctzn: %d\n", ctnz);

    return H_BC;
}

double *resolve_linear_QR(double **H_BC, double *z, long int numeroRamos, long int nmed_BC){
    cholmod_sparse *A = NULL;
    cholmod_triplet *T = NULL;
    cholmod_dense *b = NULL;
    cholmod_dense *X = NULL;
    
    cholmod_common Common, *c;

    c = &Common;
    cholmod_l_start(c);

    T = cholmod_l_allocate_triplet(3*nmed_BC, 3*numeroRamos, 3*nmed_BC*3*numeroRamos, 0, CHOLMOD_COMPLEX, c);
    A = cholmod_l_allocate_sparse(3*nmed_BC, 3*numeroRamos, 3*nmed_BC*3*numeroRamos, 0, 0, 0, CHOLMOD_COMPLEX, c);
    b = cholmod_l_allocate_dense(3*nmed_BC, 1, 3*nmed_BC, CHOLMOD_COMPLEX, c);
    X = cholmod_l_allocate_dense(3*numeroRamos, 1, 3*numeroRamos, CHOLMOD_COMPLEX, c);
    
    int index = 0;
    for (int i = 0; i < 3*nmed_BC; i++)
    {
        for (int r = 0; r < 3*numeroRamos; r++)
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

    for (int i = 0; i < 3*nmed_BC; i++)
    {
        ((double *)b->x)[(2*i)] = creal(z[i]);
        ((double *)b->x)[(2*i)+1] = cimag(z[i]);
    }
    

    A = cholmod_l_triplet_to_sparse(T, 3*nmed_BC*3*numeroRamos, c);
      
    X = SuiteSparseQR_C_backslash(SPQR_ORDERING_BEST, SPQR_DEFAULT_TOL, A, b, c);

    double *ponto;
    ponto = aloca_vetor(6*numeroRamos);
    ponto = (double *)X->x;
    //for (int ctz = 0; ctz < 10; ctz ++){
    //        printf("x[%d] = %f + i*%f\n", ctz, ponto[2*ctz], ponto[2*ctz+1]);
    //}

    
    return ponto;

}

void monta_z_complexa(DMED_COMPLEX *medidas_eq, __complex__ double *z, long int nmed_BC){
    for (int i = 0; i < nmed_BC; i++){
        switch (medidas_eq[i].fases){
            case 1:
                z[3*i] = medidas_eq[i].zmed;
                break;
            case 2:
                z[3*i+1] = medidas_eq[i].zmed;
                break;
            case 3:
                z[3*i+2] = medidas_eq[i].zmed;
                break;
            case 4:
                z[3*i] = medidas_eq[i].zmed;
                z[3*i+1] = medidas_eq[i].zmed;
                break;
            case 5:
                z[3*i] = medidas_eq[i].zmed;
                z[3*i+2] = medidas_eq[i].zmed;
                break;
            case 6:
                z[3*i+1] = medidas_eq[i].zmed;
                z[3*i+2] = medidas_eq[i].zmed;
                break;
            case 7:
                z[3*i] = medidas_eq[i].zmed;
                z[3*i+1] = medidas_eq[i].zmed;
                z[3*i+2] = medidas_eq[i].zmed;
                break;
        }
        
    } 
}



void incializa_tensoes_grafo(GRAFO *grafo, long int numeroBarras, ALIMENTADOR *alimentadores, long int numeroAlimentadores)
{
    int i, k, fase;
    BOOL visitado[numeroBarras];
    __complex__ double V0[3], **Yaux;

    Yaux = c_matAloca(3);

    //Flat start trifásico (Va = Vb = Vc = 1p.u.  Ta = 0  Tb = -120  Tc = 120) - com busca em profundidade para atualizar taps iniciais
    for (i = 0; i < numeroBarras; i++)
    {
        visitado[i] = false;
    }
    for (i = 0; i < numeroAlimentadores; i++)
    {
        //Tensão Inicial da subestação
        V0[0] =  grafo[alimentadores[i].noRaiz].barra->Vinicial[0];
        V0[1] =  grafo[alimentadores[i].noRaiz].barra->Vinicial[1];
        V0[2] =  grafo[alimentadores[i].noRaiz].barra->Vinicial[2];

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
                    { //Atualiza o V0 para trafo visto a ligação e tap
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
                    { //Para o caso de regulador de tensão
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


    //Percorre o grafo atualizando o cálculo de h(x))
    for (i = 0; i < numeroBarras; i++)
    {
        //Percorre os ramos adjacentes
        for (k = 0; k < (grafo)[i].numeroAdjacentes; k++)
        {
            aux_regua = barra[(grafo)[i].idNo].ID + barra[(grafo)[i].adjacentes[k].idNo].ID/10000.0;
            aux_regua_inv = -(barra[(grafo)[i].adjacentes[k].idNo].ID + barra[(grafo)[i].idNo].ID/10000.0);
            

            for (j = 0; j < numeroRamos; j++){
                //printf("\naux: %f, regua: %f", aux_regua, regua_x[3*j]);
                if (fabs(aux_regua - regua_x[3*j]) < EPS){
                    
                    (grafo)[i].adjacentes[k].Cur[0] = x_bc[6*j] + I * x_bc[6*j+1];
                    //printf("\ni: %d -> k: %d = %f + j*%f", i, j, x_bc[6*j], x_bc[6*j+1]);
                    //printf("\ni: %d -> k: %d = %f + j*%f", i, k, creal(grafo[i].adjacentes[k].Cur[0]), cimag(grafo[i].adjacentes[k].Cur[0]));
                    //printf("\ni: %d -> k: %d = %f + j*%f", 0, 0, creal(grafo[0].adjacentes[0].Cur[0]), cimag(grafo[0].adjacentes[0].Cur[0]));
    
                    (grafo)[i].adjacentes[k].Cur[1] = x_bc[6*j+2] + I * x_bc[6*j+3];
                    (grafo)[i].adjacentes[k].Cur[2] = x_bc[6*j+4] + I * x_bc[6*j+5];
                }
                if (fabs(regua_x[3*j] + aux_regua_inv) < EPS){
                    (grafo)[i].adjacentes[k].Cur[0] = -(x_bc[6*j] + I * x_bc[6*j+1]);
                    //printf("\ni: %d -> k: %d = %f + j*%f", i, j, x_bc[6*j], x_bc[6*j+1]);
                    //printf("\ninv i: %d -> k: %d = %f + j*%f", i, k, creal(grafo[i].adjacentes[k].Cur[0]), cimag(grafo[i].adjacentes[k].Cur[0]));
                    //printf("\ni: %d -> k: %d = %f + j*%f", 0, 0, creal(grafo[0].adjacentes[0].Cur[0]), cimag(grafo[0].adjacentes[0].Cur[0]));
    
                    (grafo)[i].adjacentes[k].Cur[1] = -(x_bc[6*j+2] + I * x_bc[6*j+3]);
                    (grafo)[i].adjacentes[k].Cur[2] = -(x_bc[6*j+4] + I * x_bc[6*j+5]);

                }

                
            }
        }
    }
    //printf("\ni: %d -> k: %d = %f + j*%f", 0, 0, creal(grafo[0].adjacentes[0].Cur[0]), cimag(grafo[0].adjacentes[0].Cur[0]));
                    
}

void exportaEstado_BC(GRAFO *grafo,  long int numeroBarras){
    int i, k, fase;
    FILE *arqout;
    
    arqout = fopen("state.txt","w+");
    
    for(i=0;i<numeroBarras;i++){
        switch (grafo[i].fases){
            case 1:
                fprintf(arqout,"%d\t%.15f\n",i,cabs(grafo[i].V[0]));
                break;
            case 2:
                fprintf(arqout,"%d\t%.15f\n",i,cabs(grafo[i].V[1]));
                break;
            case 3:
                fprintf(arqout,"%d\t%.15f\n",i,cabs(grafo[i].V[2]));
                break;
            case 4:
                fprintf(arqout,"%d\t%.15f\n",i,cabs(grafo[i].V[0]));
                fprintf(arqout,"%d\t%.15f\n",i,cabs(grafo[i].V[1]));
                break;
            case 5:
                fprintf(arqout,"%d\t%.15f\n",i,cabs(grafo[i].V[0]));
                fprintf(arqout,"%d\t%.15f\n",i,cabs(grafo[i].V[2]));
                break;
            case 6:
                fprintf(arqout,"%d\t%.15f\n",i,cabs(grafo[i].V[1]));
                fprintf(arqout,"%d\t%.15f\n",i,cabs(grafo[i].V[2]));
                break;
            case 7:
                fprintf(arqout,"%d\t%.15f\n",i,cabs(grafo[i].V[0]));
                fprintf(arqout,"%d\t%.15f\n",i,cabs(grafo[i].V[1]));
                fprintf(arqout,"%d\t%.15f\n",i,cabs(grafo[i].V[2]));
                break;
        }
    }

    for(i=0;i<numeroBarras;i++){
        switch (grafo[i].fases){
            case 1:
                fprintf(arqout,"%d\t%.15f\n",i,carg(grafo[i].V[0]));
                break;
            case 2:
                fprintf(arqout,"%d\t%.15f\n",i,carg(grafo[i].V[1]));
                break;
            case 3:
                fprintf(arqout,"%d\t%.15f\n",i,carg(grafo[i].V[2]));
                break;
            case 4:
                fprintf(arqout,"%d\t%.15f\n",i,carg(grafo[i].V[0]));
                fprintf(arqout,"%d\t%.15f\n",i,carg(grafo[i].V[1]));
                break;
            case 5:
                fprintf(arqout,"%d\t%.15f\n",i,carg(grafo[i].V[0]));
                fprintf(arqout,"%d\t%.15f\n",i,carg(grafo[i].V[2]));
                break;
            case 6:
                fprintf(arqout,"%d\t%.15f\n",i,carg(grafo[i].V[1]));
                fprintf(arqout,"%d\t%.15f\n",i,carg(grafo[i].V[2]));
                break;
            case 7:
                fprintf(arqout,"%d\t%.15f\n",i,carg(grafo[i].V[0]));
                fprintf(arqout,"%d\t%.15f\n",i,carg(grafo[i].V[1]));
                fprintf(arqout,"%d\t%.15f\n",i,carg(grafo[i].V[2]));
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
    //Alocação de memória das variáveis do estimador de estado

    nmed = 0;
    for (i = 0; i < 9; i++)
    {
        for (j = 0; j < 8; j++)
        {
            nmed = nmed + numeroMedidas[i][j];
        }
    }
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

    
    //Inicializa vetor x (correntes)
    //utilizar variavel numeroRamos

    //RNP - a partir do alimentador
    int *RNP;
    //RNP = aloca_vetor(numeroBarras);
    k = 0;
    FILABARRAS *barraAtual;
    RNP = aloca_vetor_int(numeroBarras);
    barraAtual = &(alimentadores[0].rnp[0]);
    while(barraAtual != NULL){
        RNP[k] = barraAtual->idNo;
        k++;
        barraAtual = barraAtual->prox;
    }

    nmed_BC = conta_medidas_BC(medidas, nmed);
    //cria_B_Z_ramos(grafo, numeroRamos, ramos, Sbase);
    
    //x_bc = aloca_vetor(numeroRamos);
    //inicializa_vetor_estados_BC(x_bc, 3*numeroRamos);
    //inicializar vetor de variaveis de estado
    
    DMED_COMPLEX *medidas_complexas = NULL;
    medidas_complexas = (DMED_COMPLEX *)malloc((nmed_BC) * sizeof(DMED_COMPLEX));
    
    
    double *x_bc = NULL;
    x_bc = aloca_vetor(6*numeroRamos);
    //x_bc = (double *)malloc((6*numeroRamos) * sizeof(double));
    

    DMED_COMPLEX *medidas_equivalentes = NULL;
    medidas_equivalentes = (DMED_COMPLEX *)malloc((nmed_BC) * sizeof(DMED_COMPLEX));

    __complex__ double *z_eq = NULL;
    z_eq = c_vetAloca(3*nmed_BC);
    //z_eq = (__complex__ double *)malloc(3*nmed_BC * sizeof(__complex__ double));


    double *regua_x = NULL;
    double *regua_med = NULL;
    double *regua_med_inv = NULL;
    double **H_BC = NULL;
    
    double *x_anterior = NULL;
    x_anterior = aloca_vetor(6*numeroRamos);
    double *dif_x = NULL;
    dif_x = aloca_vetor(6*numeroRamos);
    //vetor de estados: 1 para cada ramo e fase;
    regua_med = aloca_vetor(3*nmed_BC);
    regua_med_inv = aloca_vetor(3*nmed_BC);
    regua_x = aloca_vetor(3*numeroRamos);
    H_BC = aloca_matriz(3*nmed_BC, 3*numeroRamos);

    

    

    
    incializa_tensoes_grafo(grafo, numeroBarras, alimentadores, numeroAlimentadores);
    medidas_complexas = converte_medidas_para_complexo(medidas, nmed);
    medidas_equivalentes = divide_medidas_por_tensao(medidas_complexas, nmed_BC, numeroBarras, grafo);

    monta_regua_x(numeroRamos, regua_x, ramos);
    monta_regua_medidas(nmed_BC, regua_med, regua_med_inv, medidas_equivalentes);
    
    int it = 0;
    int conv = 0;
    while (conv < 1){

        monta_z_complexa(medidas_equivalentes, z_eq, nmed_BC);
        //printf("\n");

        for (int ctz = 0; ctz < 20; ctz ++){
            //printf("z[%d] = %f + i*%f\n", ctz, creal(z_eq[ctz]), cimag(z_eq[ctz]));
            //printf("reguax : %f\n", regua_x[ctz]);
        }
        //printf("\n");
        
        //monta matriz Jacobiana
        H_BC = monta_matriz_H(numeroRamos, nmed_BC, regua_x, regua_med, regua_med_inv);
        int st = 0;

        
        x_anterior = x_bc;
        x_bc = resolve_linear_QR(H_BC, z_eq, numeroRamos, nmed_BC);
        //for (int cx = 0; cx < 2; cx++){
        //    printf("xbc: %f\n", x_bc[cx]);
        //    printf("\n");
        //}

        for (int cx = 0; cx < 6*numeroRamos; cx ++){
            dif_x[cx] = x_anterior[cx] - x_bc[cx];
        }
        double nfx;

        nfx = norma_inf(dif_x, 6*numeroRamos);
        printf("\n\nIteracao:  %d \t|Dx|_inf =  %.17lf \t  \n", it, nfx);


        
        //mudar atualiza rede para receber complexo
        atualiza_Rede_BC(grafo, numeroBarras, barra, regua_x, numeroRamos, x_bc);
        
        

        //for (int nb = 0; nb < numeroBarras; nb++){
        //    for (int nj = 0; nj< grafo[nb].numeroAdjacentes; nj++){
        //        printf("\ni: %d -> k: %d = %f + j*%f", nb, nj, creal(grafo[nb].adjacentes[nj].Cur[0]), cimag(grafo[nb].adjacentes[nj].Cur[0]));
        //
        //    }
        //}
        
        for (int k = 0; k<numeroBarras; k++){
            int ct = forward_sweep(&grafo[RNP[k]], grafo);
        }
        if (nfx < tol | it > 30){
            conv = 10;
        }

        medidas_equivalentes = divide_medidas_por_tensao(medidas_complexas, nmed_BC, numeroBarras, grafo);
        
        it++;
    }


}


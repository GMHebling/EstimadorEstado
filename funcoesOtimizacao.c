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

//#include "mmio.h"
#include <cholmod.h>
#include "SuiteSparseQR_C.h"

int otimiza_Gauss_NewtonQR(double *z, double **h, double ***H, GRAFO *grafo, long int numeroBarras, DRAM *ramos, DMED *medidas, long int nvar, long int nmed, double *regua_comp, double *ponto, double tol, long int ref1, long int ref2)
{
    long int it, r;
    long int NAV = 0;
    double alpha0 = 1;

    double *Dz = NULL;
    double *b = NULL;
    double *Dx = NULL;
    double **H_rf = NULL;
    double **H_T = NULL;
    double **Gain = NULL;
    double *regua = NULL;
    int i, j, conv = 0;
    double nGx, nFx;

    double *rN = NULL, *bHat = NULL;

    FILE *arquivo;
    arquivo = fopen("iteracoesWLS.rtf", "wr+");
    //printf("nmed: %d\n", nmed);
    //printf("nvar: %d\n", nvar);

    //----------------------------------------------------------------------------
    //ALOCAÇÃO
    //
    //----------------------------------------------------------------------------
    H_rf = aloca_matriz(nmed, nvar);
    H_T = aloca_matriz(nvar, nmed);
    Gain = aloca_matriz(nvar, nvar);
    Dz = aloca_vetor(nmed);
    b = aloca_vetor(nvar);
    rN = aloca_vetor(nmed);
    bHat = aloca_vetor(nmed);

    //Inicializa o modelo no ponto inicial
    atualiza_Rede(grafo, numeroBarras);                  //atualiza a condição da rede conforme o estado atual
    atualiza_Modelo(grafo, numeroBarras, nmed, medidas); // atualiza modelo de medição conforme a condição atual da rede

    for (i = 0; i < nmed; i++)
    {
        fprintf(arquivo, "%.7lf\t%.7lf\t%.7f\n", z[i], medidas[i].h, z[i] - medidas[i].h);
    }

    //-----> LOOP DO ALGORITIMO WLS
    it = 0;
    conv = 0;
    clock_t tIni = clock();

    //inicializa vetores esparsos
    long int *i_sparse = NULL;
    long int *j_sparse = NULL;
    double *x_sparse = NULL;

    int nzeros = 0;
    i_sparse = aloca_vetor(nmed);
    j_sparse = aloca_vetor(nvar);
    x_sparse = aloca_vetor(nmed * nvar);
    int escrito = 0;
    //atualiza_H(grafo, numeroBarras, ramos, medidas, nmed);
    while ((it < 30))
    {
        clock_t t0 = clock();
        //************************************************************************
        //MONTAGEM DE H(x)
        //
        //************************************************************************
        //TODO: montar nova atualiza_H_SS

        atualiza_H(grafo, numeroBarras, ramos, medidas, nmed); //atualiza Jacobiana do modelo de medição de acordo com o estado atual

        //printf("H update\n");
        //Multiplica R_1/2*H e o modelo de medição - formulação do estimador via método QR
        for (i = 0; i < nmed; i++)
        {
            for (j = 0; j < medidas[i].nvar; j++)
            {
                medidas[i].H[j] = medidas[i].H[j] / medidas[i].sigma;
            }
            Dz[i] = (medidas[i].zmed - medidas[i].h) / medidas[i].sigma;
        }

        //Tira a coluna de angulos da referencia na matriz H(x)

        clock_t tMontaH = clock();
        free(Dx);

        //=========================================================================
        //Convergência do algoritimo - exporta resultados finais e processa erros grosseiros
        if (conv == 1)
        {
            clock_t tFim = clock();
            double tempoIt = (double)(tFim - tIni) / CLOCKS_PER_SEC;
            printf("\n\nIteracao:  %d \t|Dx|_inf =  %.17lf \t |Grad|_inf =  %.17lf \n", it, nFx, nGx);
            printf("\n\n Convergência em %d iteracoes e tempo: %.4lf", it, tempoIt);

            fprintf(arquivo, "\n\n Convergência em %d iteracoes e tempo: %.4lf", it, tempoIt);
            fprintf(arquivo, "\n\nVetor z \t\tVetor h(x) \t\t\tVetor Dz\n");
            for (i = 0; i < nmed; i++)
            {
                fprintf(arquivo, "%.7lf\t%.7lf\t%.7f\n", z[i], *h[i], z[i] - *h[i]);
            }
            saidaEstado(grafo, numeroBarras, it, tempoIt, nFx, nGx);
            free(Dz);
            free(b); //free(regua);
            //for (i=0;i<nmed;i++) free(H_rf[i]);
            //free(H_rf);
            fclose(arquivo);
            return conv;
        }
        //=========================================================================

        int nz = 0;
        //inicializacao de variaveis
        //printf("INIT SPARSE\n");
        cholmod_sparse *A_SS = NULL;
        cholmod_sparse *A_T = NULL;
        cholmod_sparse *G_S = NULL;
        cholmod_dense *b_SS = NULL;
        cholmod_dense *b_WLS = NULL;
        cholmod_dense *X_SS = NULL;
        cholmod_triplet *T_SS = NULL;
        cholmod_factor *L = NULL;

        cholmod_common Common, *c;

        c = &Common;
        cholmod_l_start(c);

        //alocacao de memoria das struturas do suitesparse
        T_SS = cholmod_l_allocate_triplet(nmed, nvar, nmed * nvar, 0, CHOLMOD_REAL, c);
        A_SS = cholmod_l_allocate_sparse(nmed, nvar, nmed * nvar, 0, 0, 0, CHOLMOD_REAL, c);
        A_T = cholmod_l_allocate_sparse(nmed, nvar, nmed * nvar, 0, 0, 0, CHOLMOD_REAL, c);
        G_S = cholmod_l_allocate_sparse(nvar, nvar, nvar * nvar, 0, 0, 0, CHOLMOD_REAL, c);
        b_SS = cholmod_l_allocate_dense(nmed, 1, nmed, CHOLMOD_REAL, c);
        b_WLS = cholmod_l_allocate_dense(nvar, 1, nmed, CHOLMOD_REAL, c);
        X_SS = cholmod_l_allocate_dense(nvar, 1, nvar, CHOLMOD_REAL, c);
        L = cholmod_l_allocate_factor(nmed, c);

        int index = 0;
        for (int i = 0; i < nmed; i++)
        {
            for (int r = 0; r < nvar; r++)
            {
                if (*H[i][r] != 0)
                {
                    ((long int *)T_SS->i)[index] = i;
                    ((long int *)T_SS->j)[index] = r;
                    ((double *)T_SS->x)[index] = *H[i][r];
                    //fprintf(mat,"%ld,%ld,%f\n",i,r,*H[i][r]);

                    T_SS->nnz += 1;
                    index += 1;
                }
            }
        }

        //escreve o vetor Dz no formato Dense
        for (int i = 0; i < nmed; i++)
        {
            ((double *)b_SS->x)[i] = Dz[i];
        }

        //converte a matrix triplet para sparse

        A_SS = cholmod_l_triplet_to_sparse(T_SS, nvar * nvar, c);
        //A_T = cholmod_l_transpose(A_SS, 2, c);

        int m1[2] = {0, 1};
        int m2[2] = {0, 1};

        cholmod_l_free_triplet(&T_SS, c);
        clock_t WriteMatrix = clock();
        // //Solucao via SuiteSparse
        int mtype = 0;

        clock_t tHouse = clock();

        X_SS = SuiteSparseQR_C_backslash(SPQR_ORDERING_BEST, SPQR_NO_TOL, A_SS, b_SS, c);

        clock_t tSolve = clock();
        Dx = (double *)X_SS->x;

        cholmod_l_finish(c);

        float passo = 1;

        for (i = 0; i < nvar; i++)
        {
            b[i] = 0;
            ponto[i] = ponto[i] + Dx[i];
            //Dx[i] += passo * Dx[i];
            for (j = 0; j < nmed; j++)
            {
                b[i] = b[i] + (*H[j][i]) * Dz[j];
            }
        }

        nGx = norma_inf(b, nvar);
        nFx = norma_inf(Dx, nvar);
        clock_t t1 = clock();
        double tempo_it = (double)(t1 - t0) / CLOCKS_PER_SEC;

        fprintf(arquivo, "\n\nIteracao:  %ld \t|Dx|_inf =  %.7lf \t |Grad|_inf =  %.7lf \n", it, nFx, nGx);
        printf("\n\nIteracao:  %ld \t|Dx|_inf =  %.17lf \t |Grad|_inf =  %.17lf \n", it, nFx, nGx);
        double tempoHouseholder = (double)(t1 - tHouse) / CLOCKS_PER_SEC;
        //printf("\nSolve QR: %lf",tempoHouseholder);
        double tempoTrataMat = (double)(tHouse - tMontaH) / CLOCKS_PER_SEC;
        //printf("\nOperacoes Matrizes: %lf",tempoTrataMat);
        double tempoMontaH = (double)(tMontaH - t0) / CLOCKS_PER_SEC;
        //printf("\nMonta H: %lf",tempoMontaH);

        //************************************************************************
        //TESTE DE CONVERGÊNCIA
        //
        //************************************************************************
        //printf("conv\n");
        for (i = 0; i < nvar; i++)
        {
            if (cabs(Dx[i]) >= tol)
            {
                conv = 0;
                break;
            }
            else
                conv = 1;
        }
        //printf("FIM conv \n");

        //************************************************************************
        //ATUALIZAÇÃO DO ESTADO DA REDE
        //
        //************************************************************************
        atualiza_estado(grafo, ponto, regua_comp, nvar); //atualiza o estado atual do grafo conforme o vetor x calculado
        //printf("estado ok\n");
        atualiza_Rede(grafo, numeroBarras);
        //printf("rede ok\n");
        atualiza_Modelo(grafo, numeroBarras, nmed, medidas);
        //printf("modelo ok\n");

        it++;
        printf(".");
    }
    //******************************************************************
    //FIM DO LOOP DO WLS
    //
    //******************************************************************
    fprintf(arquivo, "\n\nNumero maximo de iteracoes atingido %d \t|Dx|_inf =  %.7lf \t |Grad|_inf =  %.7lf \n", it, nFx, nGx);
    printf("\n\nNumero maximo de iteracoes atingido %ld \t|Dx|_inf =  %.17lf \t |Grad|_inf =  %.17lf \n", it, nFx, nGx);

    atualiza_estado(grafo, ponto, regua_comp, nvar);
    atualiza_Rede(grafo, numeroBarras);
    atualiza_Modelo(grafo, numeroBarras, nmed, medidas);

    free(Dz);
    free(b); //free(regua);
    //for (i=0;i<nmed;i++) free(H_rf[i]);
    //free(H_rf);
    fclose(arquivo);
    return (conv);
}

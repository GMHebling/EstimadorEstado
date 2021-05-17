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
#include "funcoesBadData.h"

//Calcula o resíduo normalizao para todas as medidas
void residuosNormalizados(double *rN, double *bHat, long int m, long int n, double *Dz, double **H, double **W, double *z){
    int i,j;
    double *col, *UI, **G_inv, **CovR;
    double **R, **Rt, **Ht, **Mtmp, **Mtmp2;
    
    FILE *arquivo;
    arquivo = fopen("residuoNormalizado.txt","w+");
    
    //Alocação
    col = aloca_vetor(m);
    UI = aloca_vetor(m);
    G_inv = aloca_matriz(n,n);
    CovR = aloca_matriz(m,m);
    R = aloca_matriz(m,n);
    Rt = aloca_matriz(n,m);
    Ht = aloca_matriz(n,m);
    Mtmp = aloca_matriz(m,n);
    Mtmp2 = aloca_matriz(m,m);
    
    //Calcula a matriz R via fatoração QR da matriz W^1/2.H salva em H
    QRfactorization(H, m, n, R);
    for(i=0;i<m;i++){
        for(j=i;j<n;j++){
            Rt[j][i] = R[i][j];
        }
    }
    //Retorna a matriz H original sem ponderação
    for (i=0;i<m;i++){
        for(j=0;j<n;j++){
            H[i][j] = H[i][j]/sqrt(W[i][i]);
        }
        Dz[i] = Dz[i]/sqrt(W[i][i]);
    }
    
    //Calcula o inverso da Ganho
    for (i=0;i<n;i++){
        col[i] = 1;
        
        forwardSubs(Rt,n,n,col);
        backwardSubs(R,n,n,col);
        
        for (j=0;j<n;j++){
            G_inv[j][i] = col[j];
            col[j] = 0;
        }
    }
    
    //Matriz de covariancia do resisduo estimado
    matTransp(H, m, n, Ht);
    mult_mat(H,m,n,G_inv,n,n,Mtmp);
    mult_mat(Mtmp,m,n,Ht,n,m,Mtmp2);
    
    for( i = 0; i < m; i++ ){
        for( j = 0; j<m; j++ ){
            CovR[i][j] =  - Mtmp2[i][j];
        }
        CovR[i][i] = 1/W[i][i] + CovR[i][i];  //R - H(Ginv)H'
    }
    
    for (i=0;i<m;i++){
        rN[i] = Dz[i]/pow(CovR[i][i],0.5);
        bHat[i] = (rN[i]/pow(W[i][i],0.5))/pow(CovR[i][i],0.5);
        UI[i] = pow((1-CovR[i][i]*W[i][i])/(CovR[i][i]*W[i][i]),0.5);
//        if (z[i] == 0){ //medida virtual
//            rN[i] = 0;
//            UI[i] = 0;
//            bHat[i] = 0;
//        }
    }
    
//    fprintf(arquivo,"index\tDz\trN\tbHat\tUI\n");
    for (i=0;i<m;i++){
        fprintf(arquivo,"%d\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n",i,Dz[i],rN[i],bHat[i],UI[i]);
    }
        
    fclose(arquivo);
    free(R);free(Rt);free(Ht);free(Mtmp);free(Mtmp2);free(col);
    
}

//Calcula o resíduo normalizao para todas as medidas
void residuosNormalizadosQR(double *rN, double *bHat, long int m, long int n, double *Dz, double **H, double **W, double *z){
    int i,j;
    double *col, *UI, *CovR;
    double **R, **Rt;
    
    FILE *arquivo;
    arquivo = fopen("residuoNormalizado.txt","w+");
    
    //Alocação
    col = aloca_vetor(m);
    UI = aloca_vetor(m);
    CovR = aloca_vetor(m);
    R = aloca_matriz(m,n);
    Rt = aloca_matriz(n,m);
    
    //Calcula a matriz R via fatoração QR da matriz W^1/2.H salva em H
    QRfactorization(H, m, n, R);
    for(i=0;i<m;i++){
        for(j=i;j<n;j++){
            Rt[j][i] = R[i][j];
        }
    }
    
    //Retorna a matriz H original sem ponderação
    for (i=0;i<m;i++){
        for(j=0;j<n;j++){
            if (H[i][j] != 0) H[i][j] = H[i][j]/sqrt(W[i][i]);
        }
        Dz[i] = Dz[i]/sqrt(W[i][i]);
    }
    
    double aux_K;
    for (i=0;i<m;i++){
        for (j=0;j<n;j++) col[j] = H[i][j];
        forwardSubs(Rt,n,n,col);
        
        aux_K = 0;
        for (j=0;j<n;j++) aux_K += pow(col[j],2);
        
        CovR[i] = 1/W[i][i] - aux_K;        
    }
    
    
    for (i=0;i<m;i++){
        rN[i] = Dz[i]/pow(CovR[i],0.5);
        bHat[i] = (rN[i]/pow(W[i][i],0.5))/pow(CovR[i],0.5);
        UI[i] = pow((1-CovR[i]*W[i][i])/(CovR[i]*W[i][i]),0.5);
//        if (z[i] == 0){ //medida virtual
//            rN[i] = 0;
//            UI[i] = 0;
//            bHat[i] = 0;
//        }
    }
    
    fprintf(arquivo,"index\tDz\trN\tbHat\tUI\n");
    for (i=0;i<m;i++){
        fprintf(arquivo,"%d\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n",i,Dz[i],rN[i],bHat[i],UI[i]);
    }
        
    fclose(arquivo);
    free(R);free(col);
    
}



/*

//---------------------------------------------------------------------------------------
//
//Detecção de Erros Grosseiros
//
//---------------------------------------------------------------------------------------
void erros_grosseiros(){
    int i,j,maior,maior2;
    double soma1;

    //******************************************************************
    //DETEC��O DE ERROS GROSSEIROS
    //
    //******************************************************************
    printf("\nIdentificacao de Erros Grosseiros...\n");

    //Calculo da inversa da matriz ganho G-1
    G_inv = (double **)calloc(nvar,sizeof(double*));
    for (i=0;i<nvar;i++)
        G_inv[i] = (double *)calloc(nvar,sizeof(double));

    for (i=0;i<nvar;i++)
    {
        free(aux_inv);
        free(col_inv);
        aux_inv = (double *)calloc(nvar,sizeof(double));
        aux_inv[i] = 1;
        col_inv = solve_Crout(G,nvar,nvar,aux_inv);
        for (j=0;j<nvar;j++)
        {
            G_inv[j][i] = col_inv[j];
        }
    }
    

    printf("\n\nMax Resid.Norm. = %.6f  ===> Med[%d] %s",med[maior].resid_norm,maior,med[maior].ident);
    
    //--------------------------------------------------------------------
    //
    //TESTE DOS RES�DUOS NORMALIZADOS
    //
    //--------------------------------------------------------------------
    if (fabs(med[maior].resid_norm)>3)
    {
        printf("\nErro Grosseiro Detectado pelo Teste do Residuo Normalizado");
        printf("\nRetirar medida com maior residuo normalizado: ");
        printf("\n>>>> Med[%d] %s Valor: %.3f  Residuo: %.3f",maior,med[maior].ident,med[maior].valor,med[maior].resid_norm);
        gross_err1 = maior;
    }
    else{
        printf("\nSem erros grosseiros detectados pelo Teste do Residuo Normalizado");
        gross_err1 = -1;
    }

    //--------------------------------------------------------------------
    //
    //TESTE CHI QUADRADO
    //Redundancia = smed - nvar (Graus de Liberdade teste Chi-Quadrado)
    //--------------------------------------------------------------------
    lamb = lambda((smed-nvar));
    ind_J = indice_J_est(smed);
    printf("\n\nInd_J = %.4f    lamb = %.4f",ind_J,lamb);

    if (ind_J>lamb){
        printf("\nErro Grosseiro Detectado pelo Teste Chi-2");
        printf("\nRetirar medida com maior residuo normalizado: ");
        printf("\n>>>> Med[%d] %s Valor: %.4f  Residuo: %.4f",maior,med[maior].tipo,med[maior].valor,med[maior].resid_norm);
        gross_err2 = maior;
    }
    else{
        printf("\nSem erros grosseiros detectados pelo Teste Chi-2");
        gross_err2 = -1;
    }



}*/







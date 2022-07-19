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
#include "funcoesMatematicas.h"
#include "Observabilidade.h"

#define zero 1e-15


int observabilidade_H(double ***H,int nvar, int nmed,double* regua,GRAFO *grafo,long int numeroBarras, DMED *medidas, long int **numeroMedidas,DRAM *ramos,int nref){
    int i,j,k,soma;
    int *reguam, **adj_obs;
    int obs = 0; 
    double **H_rf_o, **Ht_o;
    FILE *matHdel,*matHT;
    //arqout=fopen("Observabilidade.txt","w");
    matHdel=fopen("matHdel.txt","w");
    matHT=fopen("matHT.txt","w");
   // matHaum=fopen("matHaum.txt","w");
    printf("\n\nAnalise de Observabilidade pela matriz H(x)...\n");

    //==============================================================================
    //An�lise de Observabilidade
    //
    //
    //==============================================================================
    //fprintf(arqout,"\n\n\n\n========================================================");
    //fprintf(arqout,"\nAn�lise de Observabilidade - M�todo de Fatora��o da Matriz H(x)");
    //fprintf(arqout,"\n========================================================\n");

    //==============================================================================
    //Flat Start
    //
    //==============================================================================
  



    //Observabilidade SCADA

    //Matriz Jacobiana
    H_rf_o = aloca_matriz(nmed, nvar);
    Ht_o = aloca_matriz(nvar,nmed);

    //Matriz Jacobiana Aumentada
    //double **Haum = aloca_matriz(nmed,nvar);
    //Matriz dos Fatores Triangulares
    double **F = aloca_matriz(nvar,nmed);
    
    //************************************************************************
    //MONTAGEM DE H(x)
    //
    //************************************************************************
    atualiza_Rede(grafo, numeroBarras); //atualiza a condição da rede conforme o estado atual
    atualiza_Modelo(grafo, numeroBarras, nmed, medidas);
    atualiza_H(grafo, numeroBarras, ramos, medidas, nmed);
    
    mat_ig(H,nmed,nvar,H_rf_o);
    
    //monta_H(2);
    //monta_Hpseudo(2);

    //cat_vert(H,smed,2*NB,Hpseudo,pNIP+pNIQ+pNFP+pNFQ+pNMV,2*NB, Haum);

//    fimp_mat(matHaum,H_rf_o,nmed,nvar);
   // fclose(matHaum);
    //Tira a coluna de angulo da referencia na matriz H(x)
   /* if (NAV == 0){
        tira_col(Haum,2*NB,nmed,ref,H_rf_o);
    }
    else{
        mat_ig(Haum,nmed,2*NB, H_rf_o);
    }
    */
    fimp_mat(matHT,H_rf_o,nmed,nvar);
    //transp(H_rf_o,nmed,nvar,Ht_o);

 

    //fimp_mat(matHT,Ht_o,nvar,nmed);
    
    fclose(matHT);
    //Fatora��o da Matriz Jacobiana Transposta concatenada com a Matriz de Pseudo-Medidas
    //reguam = fatora_observabilidade(Ht_o, F,nvar, nmed,0,nref,1);

    //fprintf(arqout,"\n\nMatriz Jacobiana Fatorada Ht_delta(x)\n");
   /* for(i = 0; i < smed; i++ ){
        if (regua[i] < smed - psmed) fprintf(arqout,"%s\t\t",med[regua[i]].ident);
        if (regua[i] >= smed - psmed) fprintf(arqout,"*%s\t\t",pmed[regua[i] - (smed - psmed)].ident);
    }
    fprintf(arqout,"\n");
    for(i = 0; i < smed; i++ )
        fprintf(arqout,"%d\t\t",regua[i]);
    fprintf(arqout,"\n");
    */
   // fimp_mat(arqout,Ht_o,nvar,nmed);
    fimp_mat(matHdel,Ht_o,nvar,nmed);
    fclose(matHdel);
   // fprintf(arqout,"\n\nMatriz de Fatores Triangulares\n");
    //fimp_mat(arqout,F,nvar,nmed);
    //************************************************************************
    //AN�LISE DA MATRIZ FATORADA
    //
    //************************************************************************
    //Identifica se o sistema � observ�vel como um todo
    //obs = 1;
    //for(i = 0; i < nvar; i++ )
    //{
    //    if (fabs(Ht_o[i][i])<zero) 
    //    {
    //        if(i<nvar-3)
    //        {
    //        obs = 0;
    //        break;
    //        }
    //        else if (i>nvar-3)
    //        {
    //            printf("sistema observavel considerando referencia equilibrada\n\n");
    //            obs=1;
    //            break;
    //        }
    //    }
    //}
    //if (obs==1){
    //    //-------------------------------------------------------------------------
    //    //Identificar Pseudo-Medidas Necess�rias
    // //   fprintf(arqout,"\nPseudo-Medidas Neces�rias para restaurarar a observabilidade:\n");
    //    //printf("\tPseudo-Medidas Necesarias para restaurarar a observabilidade:\n");
    //   /* for(i = 0; i < smed - psmed; i++ )
    //    {
    //        if (regua[i] >= smed - psmed) {
    //            obs = 2;
    //            fprintf(arqout,"*%s,",pmed[regua[i] - (smed - psmed)].ident);
    //            printf("*%s,",pmed[regua[i] - (smed - psmed)].ident);
    //        }
    //    }*/
    //    //-------------------------------------------------------------------------
    //    //Identificar Medidas Cr�ticas e Medidas Redundantes
    ////    fprintf(arqout,"\nMedidas Cr�ticas Identificadas:\n");
    //    //printf("\tMedidas Criticas Identificadas:\n");
    //  /*  for(i = 0; i < nvar; i++ )
    //    {
    //        soma = 0;
    //        for(j = nvar+1; j < nmed; j++ ){
    //            if ((fabs(Ht_o[i][j])>zero) && (regua[j] < nmed))
    //                soma++;
    //        }
    //        if (soma == 0){
    //            if (reguam[j] < nmed ){
    //                fprintf(arqout,"%s,",medidas[regua[i]].id);
    //                printf("%s,",medidas[regua[i]].id);
    //            }
//
    //        }
    //    }*/
    //    //-------------------------------------------------------------------------
    //    //Identificar Conjuntos Cr�ticos de Medidas
//
//
    //}
    ////-------------------------------------------------------------------------
    ////Identificar Ilhas Observ�veis
//


    //Resultado Final
    //switch (obs){
    //case 0:
    ////    fprintf(arqout,"\n\tSistema n�o observavel como um todo.\n");
    //    printf("\n\tSistema nao observavel como um todo.\n");
    //    break;
    //case 1:
    ////    fprintf(arqout,"\n\tSistema observavel.\n");
    //    printf("\n\tSistema observavel.\n");
    //    break;
    //case 2:
    // //   fprintf(arqout,"\n\tSistema com observabilidade recuperada atrav�s de pseudomedidas.\n");
    //    printf("\n\tSistema com observabilidade recuperada atraves de pseudomedidas.\n");
    //    break;
// //}
    //
    for(i=0;i<nmed;i++) free(H_rf_o[i]);
    free(H_rf_o);
    
    for(i=0;i<nvar;i++) free(Ht_o[i]);
    free(Ht_o);
    
//    for(i=0;i<nmed;i++) free(Haum[i]);
   // free(Haum);
    
    for(i=0;i<nvar;i++) free(F[i]);
    free(F);
  //  fclose(arqout);
 
  
    return obs;
}


void printH(double ***H,int nvar, int nmed,double* regua,int nref,int it,char *buf){
    //ImprimeAmatrizJacobiana em uma iteracao especifica, no arquivo de entrada com o nome indicado
    
    int i,j,k,soma;
    int *reguam, **adj_obs;
    int obs = 0; 
    double **H_rf_o, **Ht_o;
    float **Ht_f,**H_rf_f;
    FILE* matHT;
    matHT=fopen(buf,"w");
 
    nvar=nvar+nref;



    //Matriz Jacobiana
    H_rf_o = aloca_matriz(nmed, nvar);
    
    mat_ig(H,nmed,nvar,H_rf_o);


    fimp_mat(matHT,H_rf_o,nmed,nvar);

    
    fclose(matHT);


    for(i=0;i<nmed;i++) free(H_rf_o[i]);
    free(H_rf_o);
    
}


int observabilidade_Hit(double ***H,int nvar, int nmed,double* regua,int nref,int it){
    int i,j,k,soma;
    int *reguam, **adj_obs;
    int obs = 0; 
    double **H_rf_o, **Ht_o;
    float **Ht_f,**H_rf_f;
    FILE *matHdel,*matHT,*matHdel_f;
    char buf[50];
    char buf2[50];
       //arqout=fopen("Observabilidade.txt","w");
    snprintf(buf,50,"matHdel%d.txt",it);
    snprintf(buf2,50,"matHdel_f%d.txt",it);
    matHdel=fopen(buf,"w");
    matHdel_f=fopen(buf2,"w");
    snprintf(buf,50,"matHT%d.txt",it);
    matHT=fopen(buf,"w");
   // matHaum=fopen("matHaum.txt","w");
    printf("\n\nAnalise de Observabilidade pela matriz H(x)...\n");

    //==============================================================================
    //An�lise de Observabilidade
    //
    //
    //==============================================================================
    //fprintf(arqout,"\n\n\n\n========================================================");
    //fprintf(arqout,"\nAn�lise de Observabilidade - M�todo de Fatora��o da Matriz H(x)");
    //fprintf(arqout,"\n========================================================\n");

    //==============================================================================
    //Flat Start
    //
    //==============================================================================
  
    nvar=nvar+nref;


    //Observabilidade SCADA

    //Matriz Jacobiana
    H_rf_o = aloca_matriz(nmed, nvar);
    Ht_o = aloca_matriz(nvar,nmed);
    H_rf_f = aloca_matrizf(nmed, nvar);
    Ht_f = aloca_matrizf(nvar,nmed);


    //Matriz Jacobiana Aumentada
    //double **Haum = aloca_matriz(nmed,nvar);
    //Matriz dos Fatores Triangulares
    double **F = aloca_matriz(nvar,nmed);
    
    //************************************************************************
    //MONTAGEM DE H(x)
    //
    //************************************************************************
    
    mat_ig(H,nmed,nvar,H_rf_o);
    mat_igf(H,nmed,nvar,H_rf_f);
    
    //monta_H(2);
    //monta_Hpseudo(2);

    //cat_vert(H,smed,2*NB,Hpseudo,pNIP+pNIQ+pNFP+pNFQ+pNMV,2*NB, Haum);

//    fimp_mat(matHaum,H_rf_o,nmed,nvar);
   // fclose(matHaum);
    //Tira a coluna de angulo da referencia na matriz H(x)
   /* if (NAV == 0){
        tira_col(Haum,2*NB,nmed,ref,H_rf_o);
    }
    else{
        mat_ig(Haum,nmed,2*NB, H_rf_o);
    }
    */
    fimp_mat(matHT,H_rf_o,nmed,nvar);
    transp(H_rf_o,nmed,nvar,Ht_o);
    transpf(H_rf_f,nmed,nvar,Ht_f);

 

    //fimp_mat(matHT,Ht_o,nvar,nmed);
    
    fclose(matHT);
    //Fatora��o da Matriz Jacobiana Transposta concatenada com a Matriz de Pseudo-Medidas
    //reguam = fatora_observabilidade(Ht_o, F,nvar, nmed,0,nref,it);
    //reguam = fatora_observabilidadef(Ht_f,F,nvar, nmed,0,nref,it);

    //fprintf(arqout,"\n\nMatriz Jacobiana Fatorada Ht_delta(x)\n");
   /* for(i = 0; i < smed; i++ ){
        if (regua[i] < smed - psmed) fprintf(arqout,"%s\t\t",med[regua[i]].ident);
        if (regua[i] >= smed - psmed) fprintf(arqout,"*%s\t\t",pmed[regua[i] - (smed - psmed)].ident);
    }
    fprintf(arqout,"\n");
    for(i = 0; i < smed; i++ )
        fprintf(arqout,"%d\t\t",regua[i]);
    fprintf(arqout,"\n");
    */
   // fimp_mat(arqout,Ht_o,nvar,nmed);
    fimp_mate(matHdel,Ht_o,nvar,nmed);
    fimp_matef(matHdel_f,Ht_f,nvar,nmed);
    fclose(matHdel);
    fclose(matHdel_f);
   // fprintf(arqout,"\n\nMatriz de Fatores Triangulares\n");
    //fimp_mat(arqout,F,nvar,nmed);
    //************************************************************************
    //AN�LISE DA MATRIZ FATORADA
    //
    //************************************************************************
    //Identifica se o sistema � observ�vel como um todo
    obs = 1;
    for(i = 0; i < nvar; i++ )
    {
        if (fabs(Ht_o[i][i])<zero) 
        {
            if(i<nvar-3)
            {
            obs = 0;
            break;
            }
            else if (i>nvar-3)
            {
                printf("sistema observavel considerando referencia equilibrada\n\n");
                obs=1;
                break;
            }
        }
    }
    if (obs==1){
        //-------------------------------------------------------------------------
        //Identificar Pseudo-Medidas Necess�rias
     //   fprintf(arqout,"\nPseudo-Medidas Neces�rias para restaurarar a observabilidade:\n");
        //printf("\tPseudo-Medidas Necesarias para restaurarar a observabilidade:\n");
       /* for(i = 0; i < smed - psmed; i++ )
        {
            if (regua[i] >= smed - psmed) {
                obs = 2;
                fprintf(arqout,"*%s,",pmed[regua[i] - (smed - psmed)].ident);
                printf("*%s,",pmed[regua[i] - (smed - psmed)].ident);
            }
        }*/
        //-------------------------------------------------------------------------
        //Identificar Medidas Cr�ticas e Medidas Redundantes
    //    fprintf(arqout,"\nMedidas Cr�ticas Identificadas:\n");
        //printf("\tMedidas Criticas Identificadas:\n");
      /*  for(i = 0; i < nvar; i++ )
        {
            soma = 0;
            for(j = nvar+1; j < nmed; j++ ){
                if ((fabs(Ht_o[i][j])>zero) && (regua[j] < nmed))
                    soma++;
            }
            if (soma == 0){
                if (reguam[j] < nmed ){
                    fprintf(arqout,"%s,",medidas[regua[i]].id);
                    printf("%s,",medidas[regua[i]].id);
                }

            }
        }*/
        //-------------------------------------------------------------------------
        //Identificar Conjuntos Cr�ticos de Medidas


    }
    //-------------------------------------------------------------------------
    //Identificar Ilhas Observ�veis



    //Resultado Final
    switch (obs){
    case 0:
    //    fprintf(arqout,"\n\tSistema n�o observavel como um todo.\n");
        printf("\n\tSistema nao observavel como um todo.\n");
        break;
    case 1:
    //    fprintf(arqout,"\n\tSistema observavel.\n");
        printf("\n\tSistema observavel.\n");
        break;
    case 2:
     //   fprintf(arqout,"\n\tSistema com observabilidade recuperada atrav�s de pseudomedidas.\n");
        printf("\n\tSistema com observabilidade recuperada atraves de pseudomedidas.\n");
        break;

    }
    
    for(i=0;i<nmed;i++) free(H_rf_o[i]);
    free(H_rf_o);
    
    for(i=0;i<nvar;i++) free(Ht_o[i]);
    free(Ht_o);
    
//    for(i=0;i<nmed;i++) free(Haum[i]);
   // free(Haum);
    
    for(i=0;i<nvar;i++) free(F[i]);
    free(F);
  //  fclose(arqout);
 
  
    return obs;
}
int *fatora_observabilidade(double **Ht, double **F,int nvar, int nmed, int npseud,int nref,int it)
{
 
    int i,j,k, temp;
    int *regua, *pre_ord;
    char buf[50];
    FILE *texto,*teste,*fatoracao,*fatoracao2;
    texto = fopen("fatoracao.txt","w+");//a+
    teste = fopen("MatHTteste.txt","w");
    
    //R�gua de medidas
    regua=(int*)calloc(nmed,sizeof(int));
    for(i=0;i<nmed;i++)
        regua[i] = i;

/*
    for(i=0;i<smed;i++){
        if (regua[i] < smed - npseud) fprintf(texto,"%s\t\t",med[regua[i]].ident);
        if (regua[i] >= smed - npseud) fprintf(texto,"*%s\t\t",pmed[regua[i] - (smed - npseud)].ident);
    }
    fprintf(texto,"\n");
    for(i=0;i<smed;i++){
        fprintf(texto,"%d\t\t",regua[i]);
    }
    fprintf(texto,"\n");
    fimp_mat(texto,Ht,nvar,smed);
*/

    //Pr� Ordena��o da Matriz para melhorar a fatora��o
  /* pre_ord=(int*)calloc(nmed,sizeof(int));
    for(i=0;i<nmed-npseud;i++)
    {
        for(j = 0; j < nvar; j++ ){
            if (fabs(Ht[j][i])>zero) pre_ord[i]++;//conta o número de elementos não nulos na linha
        }
    }
    for(i=0;i<nvar;i++)
    {
       if (fabs(Ht[i][i])>zero)
            temp = pre_ord[i];//se pivo é não nulo temp=numero de não nulos na linha
       else
            temp = 1000; //se pivo é nulo temp =1000
       k= i;
       for(j=i+1;j < nmed-npseud;j++){
            if ((fabs(Ht[i][j])>zero) && (pre_ord[j]<temp)){
                //se o elemento da coluna seguinte é não nulo, troca de lugar as colunas
                temp = pre_ord[j];
                k = j;
            }
       }
       if (k != i){
            perm_mat(Ht,nvar,nmed,i,k);
            temp = regua[i];
            regua[i] = regua[k];
            regua[k] = temp;
            temp = pre_ord[i];
            pre_ord[i] = pre_ord[k];
            pre_ord[k] = temp;
            //temp = regua[30 - (smed - npseud)];
       }
    }*/
 
  //  fprintf(texto,"\n\nAplica��o da Pr�-Ordena��o\n");
   /* for(i=0;i<nmed;i++){
        if (regua[i] < nmed - npseud) fprintf(texto,"%s\t\t",medidas[regua[i]].id);
        if (regua[i] >= nmed - npseud) fprintf(texto,"*%s\t\t",medidas[regua[i] - (nmed - npseud)].id);
    }
    fprintf(texto,"\n");
    for(i=0;i<nmed;i++){
        fprintf(texto,"%d\t\t",regua[i]);
    }
    fprintf(texto,"\n");*/
//    fimp_mat(texto,Ht,nvar,nmed);

    //---------------------------------------------------
    // Etapa Forward
    nvar=nvar-1;
    for(j = 0; j < nvar; j++ )
    {
      
        //Busca outra medida quando encontra pivo nulo
        if (fabs(Ht[j][j])<zero){
            for(k=j+1;k<nmed; k++ ){
                if (fabs(Ht[j][k])> zero){
                    perm_mat(Ht,nvar,nmed,j,k);
                    temp = regua[j];
                    regua[j] = regua[k];
                    regua[k] = temp;
                    break;
                }
            }
        }
        /*if(it==0)
        {
        snprintf(buf,50,"analise/fatoracao_etapas%d.txt",j);    
        fatoracao=fopen(buf,"w");
        fimp_mate(fatoracao,Ht,nvar,nmed);
        fclose(fatoracao);
        }
        */
        if (fabs(Ht[j][j])>zero){
            for( i = (j+1); i<nvar; i++ )
            {
                for( k = (j+1); k<nmed; k++ )
                {
                    if(fabs(Ht[i][j])>zero && fabs(Ht[j][k])>zero){
                        Ht[i][k] = -(Ht[i][j]/Ht[j][j])*Ht[j][k] + Ht[i][k];
                        if(fabs(Ht[i][k])<zero)
                        {
                            Ht[i][k]=0;
                        }
                    }
                
                }
                F[i][j] = -Ht[i][j]/Ht[j][j];
                Ht[i][j] = 0;
            }
        }
        /*if(it==0)
        {
        snprintf(buf,50,"analise/fatoracao_etapaspos%d.txt",j);    
        fatoracao2=fopen(buf,"w");
        fimp_mate(fatoracao2,Ht,nvar,nmed);
        fclose(fatoracao2);
        }*/
        
    }
    

    // Etapa Diagonal
  /*  for(i = 0; i < nvar; i++ )
    {
        if (fabs(Ht[i][i])>zero){
            F[i][i] = 1/Ht[i][i];
            for (j = (i);j<nmed;j++)
            {
                if (fabs(Ht[i][i])>zero)
                    Ht[i][j] = (1/Ht[i][i])*Ht[i][j];
            }
           // Ht[i][i] = 1;
        }
    }*/

    //Etapa Backward
    /*for( j = (nvar-1); j >= 1; j-- )
	{
        if (fabs(Ht[j][j])>zero){
            for( i = (j-1); i>=0; i-- )
            {
                F[i][j] = -Ht[i][j]/Ht[j][j];
                for( k = j+1; k<nmed; k++)
                {
                    Ht[i][k] = -(Ht[i][j]/Ht[j][j])*Ht[j][k] + Ht[i][k];
                     if(fabs(Ht[i][k])<zero)
                    {
                        Ht[i][k]=0;
                    }
                }
                Ht[i][j] = 0;
            }
        }
	}*/
    //fimp_mat(teste,Ht,nvar,nmed);
 
//	fprintf(texto,"\n\nMatriz Fatorada\n");
   /*for(i=0;i<nmed;i++){
        if (regua[i] < nmed - npseud) fprintf(texto,"%s\t\t",med[regua[i]].ident);
        if (regua[i] >= nmed - npseud) fprintf(texto,"*%s\t\t",pmed[regua[i] - (nmed - npseud)].ident);
    }
    fprintf(texto,"\n");
    for(i=0;i<smed;i++){
        fprintf(texto,"%d\t\t",regua[i]);
    }
    fprintf(texto,"\n");
    fimp_mat(texto,Ht,nvar,smed);

    fprintf(texto,"\n\nMatriz de Fatores Trianguladres\n");
    for(i=0;i<smed;i++){
        if (regua[i] < smed - npseud) fprintf(texto,"%s\t\t",med[regua[i]].ident);
        if (regua[i] >= smed - npseud) fprintf(texto,"*%s\t\t",pmed[regua[i] - (smed - npseud)].ident);
    }
    fprintf(texto,"\n");
    for(i=0;i<smed;i++){
        fprintf(texto,"%d\t\t",regua[i]);
    }
    fprintf(texto,"\n");*/
    //fimp_mat(texto,F,nvar,nmed);

    fclose(texto);
    fclose(teste);
	return regua;
}

int *fatora_observabilidadef(float **Ht, double **F,int nvar, int nmed, int npseud,int nref,int it)
{
 
    int i,j,k, temp;
    int *regua, *pre_ord;
    char buf[50];
    FILE *texto,*teste,*fatoracao,*fatoracao2;
    texto = fopen("fatoracao.txt","w+");//a+
    teste = fopen("MatHTteste.txt","w");
    
    //R�gua de medidas
    regua=(int*)calloc(nmed,sizeof(int));
    for(i=0;i<nmed;i++)
        regua[i] = i;

/*
    for(i=0;i<smed;i++){
        if (regua[i] < smed - npseud) fprintf(texto,"%s\t\t",med[regua[i]].ident);
        if (regua[i] >= smed - npseud) fprintf(texto,"*%s\t\t",pmed[regua[i] - (smed - npseud)].ident);
    }
    fprintf(texto,"\n");
    for(i=0;i<smed;i++){
        fprintf(texto,"%d\t\t",regua[i]);
    }
    fprintf(texto,"\n");
    fimp_mat(texto,Ht,nvar,smed);
*/

    //Pr� Ordena��o da Matriz para melhorar a fatora��o
  /* pre_ord=(int*)calloc(nmed,sizeof(int));
    for(i=0;i<nmed-npseud;i++)
    {
        for(j = 0; j < nvar; j++ ){
            if (fabs(Ht[j][i])>zero) pre_ord[i]++;//conta o número de elementos não nulos na linha
        }
    }
    for(i=0;i<nvar;i++)
    {
       if (fabs(Ht[i][i])>zero)
            temp = pre_ord[i];//se pivo é não nulo temp=numero de não nulos na linha
       else
            temp = 1000; //se pivo é nulo temp =1000
       k= i;
       for(j=i+1;j < nmed-npseud;j++){
            if ((fabs(Ht[i][j])>zero) && (pre_ord[j]<temp)){
                //se o elemento da coluna seguinte é não nulo, troca de lugar as colunas
                temp = pre_ord[j];
                k = j;
            }
       }
       if (k != i){
            perm_mat(Ht,nvar,nmed,i,k);
            temp = regua[i];
            regua[i] = regua[k];
            regua[k] = temp;
            temp = pre_ord[i];
            pre_ord[i] = pre_ord[k];
            pre_ord[k] = temp;
            //temp = regua[30 - (smed - npseud)];
       }
    }*/
 
  //  fprintf(texto,"\n\nAplica��o da Pr�-Ordena��o\n");
   /* for(i=0;i<nmed;i++){
        if (regua[i] < nmed - npseud) fprintf(texto,"%s\t\t",medidas[regua[i]].id);
        if (regua[i] >= nmed - npseud) fprintf(texto,"*%s\t\t",medidas[regua[i] - (nmed - npseud)].id);
    }
    fprintf(texto,"\n");
    for(i=0;i<nmed;i++){
        fprintf(texto,"%d\t\t",regua[i]);
    }
    fprintf(texto,"\n");*/
//    fimp_mat(texto,Ht,nvar,nmed);

    //---------------------------------------------------
    // Etapa Forward
    for(j = 0; j < nvar; j++ )
    {
      
        //Busca outra medida quando encontra pivo nulo
        if (fabs(Ht[j][j])<zero){
            for(k=j+1;k<nmed; k++ ){
                if (fabs(Ht[j][k])> zero){
                    perm_matf(Ht,nvar,nmed,j,k);
                    temp = regua[j];
                    regua[j] = regua[k];
                    regua[k] = temp;
                    break;
                }
            }
        }
        if(it==0)
        {
        snprintf(buf,50,"analise/fatoracao_etapas%d.txt",j);    
        fatoracao=fopen(buf,"w");
        fimp_matef(fatoracao,Ht,nvar,nmed);
        fclose(fatoracao);
        }
        
        if (fabs(Ht[j][j])>zero){
            for( i = (j+1); i<nvar; i++ )
            {
                for( k = (j+1); k<nmed; k++ )
                {
                    if(fabs(Ht[i][j])>zero && fabs(Ht[j][k])>zero){
                        Ht[i][k] = -(Ht[i][j]/Ht[j][j])*Ht[j][k] + Ht[i][k];
                        if(fabs(Ht[i][k])<zero)
                        {
                            Ht[i][k]=0;
                        }
                    }
                
                }
               // F[i][j] = -Ht[i][j]/Ht[j][j];
                Ht[i][j] = 0;
            }
        }
        if(it==0)
        {
        snprintf(buf,50,"analise/fatoracao_etapaspos%d.txt",j);    
        fatoracao2=fopen(buf,"w");
        fimp_matef(fatoracao2,Ht,nvar,nmed);
        fclose(fatoracao2);
        }
        
    }
    

    // Etapa Diagonal
    /*for(i = 0; i < nvar; i++ )
    {
        if (fabs(Ht[i][i])>zero){
           // F[i][i] = 1/Ht[i][i];
            for (j = (i);j<nmed;j++)
            {
                if (fabs(Ht[i][i])>zero)
                    Ht[i][j] = (1/Ht[i][i])*Ht[i][j];
            }
            //Ht[i][i] = 1;
        }
    }*/

    //Etapa Backward
    /*for( j = (nvar-1); j >= 1; j-- )
	{
        if (fabs(Ht[j][j])>zero){
            for( i = (j-1); i>=0; i-- )
            {
              //  F[i][j] = -Ht[i][j]/Ht[j][j];
                for( k = j+1; k<nmed; k++)
                {
                    Ht[i][k] = -(Ht[i][j]/Ht[j][j])*Ht[j][k] + Ht[i][k];
                     if(fabs(Ht[i][k])<zero)
                    {
                        Ht[i][k]=0;
                    }
                }
                Ht[i][j] = 0;
            }
        }
	}*/
    fimp_matef(teste,Ht,nvar,nmed);
 
//	fprintf(texto,"\n\nMatriz Fatorada\n");
   /*for(i=0;i<nmed;i++){
        if (regua[i] < nmed - npseud) fprintf(texto,"%s\t\t",med[regua[i]].ident);
        if (regua[i] >= nmed - npseud) fprintf(texto,"*%s\t\t",pmed[regua[i] - (nmed - npseud)].ident);
    }
    fprintf(texto,"\n");
    for(i=0;i<smed;i++){
        fprintf(texto,"%d\t\t",regua[i]);
    }
    fprintf(texto,"\n");
    fimp_mat(texto,Ht,nvar,smed);

    fprintf(texto,"\n\nMatriz de Fatores Trianguladres\n");
    for(i=0;i<smed;i++){
        if (regua[i] < smed - npseud) fprintf(texto,"%s\t\t",med[regua[i]].ident);
        if (regua[i] >= smed - npseud) fprintf(texto,"*%s\t\t",pmed[regua[i] - (smed - npseud)].ident);
    }
    fprintf(texto,"\n");
    for(i=0;i<smed;i++){
        fprintf(texto,"%d\t\t",regua[i]);
    }
    fprintf(texto,"\n");*/
    //fimp_mat(texto,F,nvar,nmed);

    fclose(texto);
    fclose(teste);
	return regua;
}
void fimp_mat(FILE *arquivo, double **A,int m, int n)
{
    int i,j;

    for (i=0;i<m;i++)
    {
        for (j=0;j<n-1;j++)
        {
            fprintf(arquivo,"%.15e,",A[i][j]);
        }
        fprintf(arquivo,"%.15e\n",A[i][j]);
    }

}

void fimp_mate(FILE *arquivo, double **A,int m, int n)
{
    int i,j;

    for (i=0;i<m;i++)
    {
        for (j=0;j<n;j++)
        {
            fprintf(arquivo,"%.12e\t",A[i][j]);
        }
        fprintf(arquivo,"\n");
    }

}
void fimp_matef(FILE *arquivo, float **A,int m, int n)
{
    int i,j;

    for (i=0;i<m;i++)
    {
        for (j=0;j<n;j++)
        {
            fprintf(arquivo,"%.7e\t",A[i][j]);
        }
        fprintf(arquivo,"\n");
    }

}
void fimp_vet(FILE *arquivo, double *A,int m)
{
    int i;

    for (i=0;i<m;i++)
    {
        fprintf(arquivo,"%.6f\n",A[i]);
    }
}
void perm_mat(double **A,int m, int n, int col1, int col2)
{
    int i;
    double temp;

    for (i=0;i<m;i++)
    {
        temp = A[i][col1];
        A[i][col1] = A[i][col2];
        A[i][col2] = temp;
    }
}
void perm_matf(float **A,int m, int n, int col1, int col2)
{
    int i;
    float temp;

    for (i=0;i<m;i++)
    {
        temp = A[i][col1];
        A[i][col1] = A[i][col2];
        A[i][col2] = temp;
    }
}
void transp(double **A,int m,int n, double **B){
     int i,j;

     for (i=0;i<n;i++)
     {
         for (j=0;j<m;j++)
         {
             B[i][j]=A[j][i];
         }//if (fabs(Ht[i][k])>zero)

     }

}
void transpf(float **A,int m,int n, float **B){
     int i,j;

     for (i=0;i<n;i++)
     {
         for (j=0;j<m;j++)
         {
             B[i][j]=A[j][i];
         }//if (fabs(Ht[i][k])>zero)

     }

}

void mat_igf(double ***A, int m, int n, float **B)
{
    int i, j;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            B[i][j] = (float)*A[i][j];
        }
    }
}

float **aloca_matrizf(int m, int n)
{
    int i, j;
    float **A;

    A = (float**)malloc(m * sizeof(float *));
    for (i = 0; i < m; i++)
        A[i] = (float *)malloc(n * sizeof(float));

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            A[i][j] = 0;
        }
    }
    return (A);
}
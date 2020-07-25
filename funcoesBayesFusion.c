#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>

#include "data_structures.h"
#include "funcoesWLS.h"
#include "funcoesLeitura.h"
#include "funcoesTopologia.h"
#include "funcoesCalculoEletrico.h"
#include "funcoesOtimizacao.h"
#include "funcoesMatematicas.h"
#include "funcoesBadData.h"

#include "/home/julio/SuiteSparse-5.6.0/include/cholmod.h"
#include "/home/julio/SuiteSparse-5.6.0/include/SuiteSparseQR_C.h"

#define Long SuiteSparse_long

#define WITH_PRIOR 1

#define maxIT 100

// Importa arquivos de Texto com as matrizes de covariância a priori (formato SuiteSparse) e vetor x a priori
void importaPrioriQR_ss(double *x0, cholmod_sparse **R_SS,cholmod_common *c){
    int i,j,k,r, fase;
    long int nvar, nmed;
    char blocoLeitura[2000]; 
    char *dados;
    
    cholmod_triplet *T_SS = NULL;
    
    //Importa arquivo com a matriz de covariância do estado (R da fatoração, onde Cov_X = Rt * R)
    FILE *file_R;
    file_R = fopen("Priori_MatrizR.txt","r");
    // *R_SS = cholmod_l_read_sparse(file_R, c); // Aparentemente fica mais rápido se ler com o fgets e converter

    fgets(blocoLeitura, 2000, file_R);
    fgets(blocoLeitura, 2000, file_R);
    dados = blocoLeitura;
    nvar = (getfield_int(dados,1));
    int nz = (getfield_int(dados,3));
    T_SS = cholmod_l_allocate_triplet(nvar, nvar, nz, 0, CHOLMOD_REAL, c);

    int index = 0;
    while( (fgets(blocoLeitura, 2000, file_R))!= NULL ){
        dados = blocoLeitura;
        ((long int*)T_SS->i)[index] = (getfield_int(dados,1)) -1 ;
        ((long int*)T_SS->j)[index] = (getfield_int(dados,2)) -1;
        ((double*)T_SS->x)[index] = (getfield_double(dados,3));
        T_SS->nnz += 1;
        index += 1;
        *dados = NULL;
    }
    // cholmod_l_print_triplet(T_SS,"R",c);
    *R_SS = cholmod_l_triplet_to_sparse(T_SS, T_SS->nnz, c);   
    // cholmod_l_print_sparse(R_SS,"R",c);
    fclose(file_R);

    //Importa informações para estimação a posteriori
    FILE *priori;
    priori = fopen("Priori_VetorX.txt","r");
    fgets(blocoLeitura, 2000, priori);
    dados = blocoLeitura;
    nvar = (getfield_int(dados,1));
    nmed = (getfield_int(dados,2));  
    i=0;
    while( (fgets(blocoLeitura, 2000, priori))!= NULL ){
        dados = blocoLeitura;
        x0[i] = (getfield_double(dados,1));
        i++;
        if (i==nvar) i = 0;     
    }
    fclose(priori);
}


// Exporta arquivos de Texto com as matrizes de covariância a priori (formato SuiteSparse) e vetor x a priori
void exportaPrioriQR_ss(GRAFO *grafo,double *regua, long int nvar, long int nmed, long int polar, cholmod_sparse *A_SS, cholmod_common *c){
    int i,j,k,r, fase;
    
    // Alocação suitesparse
    cholmod_sparse **Q;     // m-by-e sparse matrix
    cholmod_sparse **R;     // e-by-n sparse matrix
    cholmod_dense *R_dense; 
    int **E;

    //Fatora a Matriz Jacobiana (recebida na entrada no formato esparso) para obter Q e R
    int rank = SuiteSparseQR_C_QR(SPQR_ORDERING_BEST, SPQR_NO_TOL, nvar, A_SS, &Q, &R, &E, c);
    
    // clock_t tIni = clock();
    
    //Exporta arquivo com a matriz de covariância do estado (R da fatoração, onde Cov_X = Rt * R)
    FILE *file_R;
    file_R = fopen("Priori_MatrizR.txt","w");
    // cholmod_l_write_sparse(file_R, R, NULL, NULL, c); // Aparentemente fica mais rápido se imprimir direto no fprintf
    R_dense = cholmod_l_sparse_to_dense(R,c);
    fprintf(file_R,"%%MatrixMarket matrix coordinate real general\n");
    fprintf(file_R,"%d %d %d\n",nvar,nvar,cholmod_l_nnz(R,c));    
    for(i=0;i<nvar*nvar;i++){
        if (((double*)R_dense->x)[i] !=0) fprintf(file_R,"%d %d %lf\n", (int) i % nvar + 1 , (int) i / nvar + 1, ((double*)R_dense->x)[i]);
    }
    fclose(file_R);

    // clock_t tPrint = clock();
    // double tempo = (double)(tPrint-tIni)/CLOCKS_PER_SEC;
    // printf("\nExporta: %lf",tempo);

    //Exporta informações para estimação a posteriori
    FILE *priori;
    priori = fopen("Priori_VetorX.txt","w");
    
    fprintf(priori,"%d\t%d",nvar,nmed);
    if (polar == 1){
        for(i=0; i<nvar; i++){ 
            if (regua[i]>0){
                k = (int) regua[i];
                fase = (int) cabs((regua[i] - k)*10);                    
                fprintf(priori,"\n%.14lf",cabs(grafo[k].V[fase]));
            }
            else{
                k = (int) regua[i];
                fase = (int)cabs((regua[i] - k)*10);                    
                k = -k;
                fprintf(priori,"\n%.14lf",carg(grafo[k].V[fase]));
            }
        }
    }
    else{
        for(i=0; i<nvar; i++){ 
            if (regua[i]>0){
                k = (int) regua[i];
                fase = (int) cabs((regua[i] - k)*10);                    
                fprintf(priori,"\n%.14lf",__real__ grafo[k].V[fase]);
            }
            else{
                k = (int) regua[i];
                fase = (int) cabs((regua[i] - k)*10);                    
                k = -k;
                fprintf(priori,"\n%.14lf",__imag__ grafo[k].V[fase]);
            }
        }
    }
    fclose(priori);
}

//Função que atualiza os valores da matriz Jacobiana do modelo de medição (H esparsa)
void atualiza_H_ss(GRAFO *grafo, long int numeroBarras, DRAM *ramos, DMED *medidas, long int numeroMedidas, cholmod_triplet *T_SS){
    int i,j,k,idMed, de, para,adj,ramo,opt,fase;
    __complex__ double *dSaux, *dIaux;
    int it = 0, r;
    int index = 0;

    dSaux = (__complex__ double*)malloc(3 * sizeof(__complex__ double));
    dSaux[0] = 0;
    dSaux[1] = 0;
    dSaux[2] = 0;
    dIaux = (__complex__ double*)malloc(3 * sizeof(__complex__ double));
    dIaux[0] = 0;
    dIaux[1] = 0;
    dIaux[2] = 0;
    
    //Inicia a contagem de nonzero elements na matriz esparsa
    T_SS->nnz = 0;
    //Atualiza a matriz H de acordo com a régua
    for(i=0;i<numeroMedidas;i++){
        switch (medidas[i].tipo){
            case 0: //0: Fluxo de Potência Ativa - kW
                for(j=0;j<medidas[i].nvar;j++){
                    k = (int) medidas[i].reguaH[j];
                    fase = (int) cabs((medidas[i].reguaH[j] - k)*10);
                    
                    de = medidas[i].k;
                    para = medidas[i].m;
                    ramo = medidas[i].ramo;
                    
                    if (ramos[ramo].k == de){
                        if(medidas[i].reguaH[j]>=0){
                            if (de == k) opt = 0;
                            else opt = 2;
                        }
                        else{
                            if (de == -k) opt = 1;
                            else opt = 3;
                        }
                        dSkm(&grafo[de], &grafo[para], &ramos[ramo], dSaux, opt, fase);
                        medidas[i].H[j] = __real__ dSaux[medidas[i].fases - 1];
                    }
                    else{
                        if(medidas[i].reguaH[j]>=0){
                            if (de == k) opt = 2;
                            else opt = 0;
                        }
                        else{
                            if (de == -k) opt = 3;
                            else opt = 1;
                        }
                        dSmk(&grafo[para], &grafo[de], &ramos[ramo], dSaux, opt, fase);
                        medidas[i].H[j] = __real__ dSaux[medidas[i].fases - 1];
                    }
                    
                    // Preenche a matrix esparsa
                    r = medidas[i].reguaH_loc[j];
                    if (r != -1){
                        ((long int*)T_SS->i)[index] = i;
                        ((long int*)T_SS->j)[index] = r;
                        
                        ((double*)T_SS->x)[index] = medidas[i].H[j]/medidas[i].sigma;

                        T_SS->nnz += 1;
                        index += 1;
                    }
                }
                break;
            case 1: //1: Fluxo de Potência Reativa - kVAr
                for(j=0;j<medidas[i].nvar;j++){
                    k = (int) medidas[i].reguaH[j];
                    fase = (int) cabs((medidas[i].reguaH[j] - k)*10);
                    
                    de = medidas[i].k;
                    para = medidas[i].m;
                    ramo = medidas[i].ramo;
                    
                    if (ramos[ramo].k == de){
                        if(medidas[i].reguaH[j]>=0){
                            if (de == k) opt = 0;
                            else opt = 2;
                        }
                        else{
                            if (de == -k) opt = 1;
                            else opt = 3;
                        }
                        dSkm(&grafo[de], &grafo[para], &ramos[ramo], dSaux, opt, fase);
                        medidas[i].H[j] = __imag__ dSaux[medidas[i].fases - 1];
                    }
                    else{
                        if(medidas[i].reguaH[j]>=0){
                            if (de == k) opt = 2;
                            else opt = 0;
                        }
                        else{
                            if (de == -k) opt = 3;
                            else opt = 1;
                        }
                        dSmk(&grafo[para], &grafo[de], &ramos[ramo], dSaux, opt, fase);
                        medidas[i].H[j] = __imag__ dSaux[medidas[i].fases - 1];
                    }
                    
                    // Preenche a matrix esparsa
                    r = medidas[i].reguaH_loc[j];
                    if (r != -1){
                        ((long int*)T_SS->i)[index] = i;
                        ((long int*)T_SS->j)[index] = r;
                        
                        ((double*)T_SS->x)[index] = medidas[i].H[j]/medidas[i].sigma;
                        T_SS->nnz += 1;
                        index += 1;
                    }
                }
                break;    
            case 2: //2: Injeção de Potência Ativa - kW
                for(j=0;j<medidas[i].nvar;j++){
                    k = (int) medidas[i].reguaH[j];
                    fase = (int) cabs((medidas[i].reguaH[j] - k)*10);
                    
                    de = medidas[i].k;
                    
                    if(medidas[i].reguaH[j]>=0){
                        opt = 0;
                    }
                    else{
                        opt = 1;
                    }
                    dSk(grafo, de, dSaux, opt, cabs(k), fase);
                   
                    medidas[i].H[j] = __real__ dSaux[medidas[i].fases - 1];
                    
                    // Preenche a matrix esparsa
                    r = medidas[i].reguaH_loc[j];
                    if (r != -1){
                        ((long int*)T_SS->i)[index] = i;
                        ((long int*)T_SS->j)[index] = r;
                        
                        ((double*)T_SS->x)[index] = medidas[i].H[j]/medidas[i].sigma;
                        T_SS->nnz += 1;
                        index += 1;
                    }
                }
                break;     
            case 3: //3: Injeção de Potência Reativa - kVAr
                for(j=0;j<medidas[i].nvar;j++){
                    k = (int) medidas[i].reguaH[j];
                    fase = (int) cabs((medidas[i].reguaH[j] - k)*10);
                    
                    de = medidas[i].k;
                    
                    if(medidas[i].reguaH[j]>=0){
                        opt = 0;
                    }
                    else{
                        opt = 1;
                    }
                    dSk(grafo, de, dSaux, opt, cabs(k), fase);
                   
                    medidas[i].H[j] = __imag__ dSaux[medidas[i].fases - 1];
                    
                    // Preenche a matrix esparsa
                    r = medidas[i].reguaH_loc[j];
                    if (r != -1){
                        ((long int*)T_SS->i)[index] = i;
                        ((long int*)T_SS->j)[index] = r;
                        
                        ((double*)T_SS->x)[index] = medidas[i].H[j]/medidas[i].sigma;
                        T_SS->nnz += 1;
                        index += 1;
                    }
                }break;
            case 4: //4: Magnitude de tensão
                for(j=0;j<medidas[i].nvar;j++){
                    k = (int) medidas[i].reguaH[j];
                    fase = (int) cabs((medidas[i].reguaH[j] - k)*10);
                    
                    if ((medidas[i].fases -1 == fase)&&(medidas[i].reguaH[j]>0))
                        medidas[i].H[j] = 1;
                    
                    // Preenche a matrix esparsa
                    r = medidas[i].reguaH_loc[j];
                    if (r != -1){
                        ((long int*)T_SS->i)[index] = i;
                        ((long int*)T_SS->j)[index] = r;
                        
                        ((double*)T_SS->x)[index] = medidas[i].H[j]/medidas[i].sigma;
                        T_SS->nnz += 1;
                        index += 1;
                    }
                }
                break;
            case 5: //5: Ângulo de tensão
                for(j=0;j<medidas[i].nvar;j++){
                    k = (int) medidas[i].reguaH[j];
                    fase = (int) cabs((medidas[i].reguaH[j] - k)*10);
                    
                    if ((medidas[i].fases -1 == fase)&&(medidas[i].reguaH[j]<0))
                        medidas[i].H[j] = 1;
                    
                    // Preenche a matrix esparsa
                    r = medidas[i].reguaH_loc[j];
                    if (r != -1){
                        ((long int*)T_SS->i)[index] = i;
                        ((long int*)T_SS->j)[index] = r;
                        
                        ((double*)T_SS->x)[index] = medidas[i].H[j]/medidas[i].sigma;
                        T_SS->nnz += 1;
                        index += 1;
                    }
                }
                break;
            case 7: //7: Magnitude de Corrente
                for(j=0;j<medidas[i].nvar;j++){
                    k = (int) medidas[i].reguaH[j];
                    fase = (int) cabs((medidas[i].reguaH[j] - k)*10);
                    
                    de = medidas[i].k;
                    para = medidas[i].m;
                    ramo = medidas[i].ramo;
                    
                    if (ramos[ramo].k == de){
                        if(medidas[i].reguaH[j]>=0){
                            if (de == k) opt = 0;
                            else opt = 2;
                        }
                        else{
                            if (de == -k) opt = 1;
                            else opt = 3;
                        }
                        if (it == 0) dSkm(&grafo[de], &grafo[para], &ramos[ramo], dSaux, opt, fase);
                        else dIkm(&grafo[de], &grafo[para], &ramos[ramo], dSaux, opt, fase);
                        medidas[i].H[j] = dSaux[medidas[i].fases - 1];
                    }
                    else{
                        if(medidas[i].reguaH[j]>=0){
                            if (de == k) opt = 2;
                            else opt = 0;
                        }
                        else{
                            if (de == -k) opt = 3;
                            else opt = 1;
                        }
                        dImk(&grafo[para], &grafo[de], &ramos[ramo], dSaux, opt, fase);
                        medidas[i].H[j] = __real__ dSaux[medidas[i].fases - 1];
                    }
                    
                    // Preenche a matrix esparsa
                    r = medidas[i].reguaH_loc[j];
                    if (r != -1){
                        ((long int*)T_SS->i)[index] = i;
                        ((long int*)T_SS->j)[index] = r;
                        
                        ((double*)T_SS->x)[index] = medidas[i].H[j]/medidas[i].sigma;
                        T_SS->nnz += 1;
                        index += 1;
                    }
                }
                break;    
                
        }

        
             
    }
    free(dSaux);free(dIaux);
}

//Análise residual do Estimador de Estado Bayesiano - Rotina com o suitesparse
// Calcula resíduos, resíduos normalizados, b-chapéu (distância do Cook), inidce UI, Intervalo de Credibilidade, Modelo de Variância a Posteriori
//
// recebe:
// 
// A_SS: Jacobiana calculada no vetor x estimado
void annaliseResidual_MAP_QR_ss(double *rN, double *bHat, double *IC, long int m, long int n, double *x, double *x0, DMED *medidas, cholmod_sparse *A_SS){
    int i,j, r;
    double *UI, *CovR, *Dz;
    
    //Alocação
    Dz = aloca_vetor(m);
    UI = aloca_vetor(m);
    CovR = aloca_vetor(m);

    
    
    // Alocação suitesparse
    cholmod_dense *Ht_SS = NULL;
    cholmod_dense *X_SS= NULL;
    SuiteSparseQR_C_factorization *QR;
    cholmod_common Common, *c;

    c = &Common;
    cholmod_l_start(c);    

    // Ht_SS = cholmod_l_allocate_dense(n, m, n, CHOLMOD_REAL, c);
    Ht_SS = cholmod_l_zeros(n,m + n,CHOLMOD_REAL, c);
    X_SS = cholmod_l_allocate_dense(n, m + n, n, CHOLMOD_REAL, c);

    //Retorna a matriz H original sem ponderação
    for (i=0;i<n;i++){
        ((double*)Ht_SS->x)[i + i*n] = 1;
        Dz[i + n] = x[i] - x0[i];
    }
    for (i=0;i<m;i++){
        for(j=0;j<medidas[i].nvar;j++){
            r = medidas[i].reguaH_loc[j];
            if (r != -1) ((double*)Ht_SS->x)[(r + n) + i*n] =  medidas[i].H[j];
        }
        Dz[i + n] = medidas[i].zmed - medidas[i].h;
    }

    //Fatora a Matriz Jacobiana (recebida na entrada no formato esparso) para obter Q e R
    QR = SuiteSparseQR_C_factorize(SPQR_ORDERING_BEST, SPQR_NO_TOL, A_SS,c); 
    
    //Soluciona X = Rt\(Et*B) (para cálculo somente da diagonal da matriz de covariância dos resíduos e matriz de permutação)
    X_SS = SuiteSparseQR_C_solve(SPQR_RTX_EQUALS_ETB, QR, Ht_SS, c);
    
    //Cálculo da diagonal de matriz de covariância do resíduos
    double aux_K;
    for (i=0;i<m;i++){
        aux_K = 0;
        for (j=0;j<n;j++){
            aux_K += pow(((double*)X_SS->x)[j + i*n],2);
        }
        // printf("\n %d\t%.15lf\t%.15lf ", i , medidas[i].sigma ,aux_K);
        CovR[i] = pow(medidas[i].sigma,2) - aux_K;  
    }
    SuiteSparseQR_C_free(&QR,c);
    
    cholmod_l_finish(c);

    //Cálculo dos resíduos normalizados, b-chapeu e indice UI
    for (i=0;i<m;i++){
        if (fabs(Dz[i]) > 0.00000001 ){
            rN[i] = fabs(Dz[i])/pow(CovR[i],0.5);
            bHat[i] = (rN[i]*medidas[i].sigma)/pow(CovR[i],0.5);
            UI[i] = pow((1-CovR[i]/pow(medidas[i].sigma,2))/(CovR[i]/pow(medidas[i].sigma,2)),0.5);
        }
        else if (medidas[i].zmed == 0){ //medida virtual
           rN[i] = 0;
           UI[i] = 0;
           bHat[i] = 0;
       }
       else{ //medida crítica
           rN[i] = 0;
           UI[i] = 0;
           bHat[i] = 0;

       }
    }
    
    //Imprime arquivo de saída Arquivo de saída de resultados
    FILE *arquivo;
    arquivo = fopen("residuoNormalizado.txt","w+");
    fprintf(arquivo,"index\tDz\trN\tbHat\tUI\n");
    for (i=0;i<m;i++){
        fprintf(arquivo,"%d\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n",i,Dz[i],rN[i],bHat[i],UI[i]);
    }
        
    fclose(arquivo);
    free(Dz);free(UI);free(CovR);
    
}


//------------------------------------------------------------------------------
//
//Método Ortogonal para a solução do Estimador - Gauss Newton via QR factorization Esparso
//
int otimiza_Gauss_NewtonQR_ss(double *z, double **h, GRAFO *grafo, long int numeroBarras, DRAM *ramos, DMED *medidas, long int nvar, long int nmed, double *regua_comp, double *ponto, double tol){
    long int it,r;
    long int NAV = 0;
    double alpha0 = 1;
    double *Dz = NULL;
    double *Dx = NULL;
    double *regua = NULL;
    int i,j, conv = 0;
    double nGx,nFx;
    double *rN = NULL, *bHat = NULL;
    
    //Create basic scalars
    double one [2], zero [2];
    one [0] = 1;
    one [1] = 0 ;
    zero [0] = 0 ;
    zero [1] = 0 ;

    //----------------------------------------------------------------------------
    //ALOCAÇÃO
    //
    //----------------------------------------------------------------------------
    Dz = aloca_vetor(nmed);
    rN = aloca_vetor(nmed);
    bHat = aloca_vetor(nmed);

    //=========================================================================
    int nz = 0;
    //inicializacao de variaveis
    cholmod_sparse *A_SS = NULL;
    // cholmod_sparse *A_T = NULL;
    // cholmod_sparse *G_S = NULL;
    cholmod_dense *b_SS = NULL;
    cholmod_dense *b_WLS = NULL;
    cholmod_dense *X_SS = NULL;
    cholmod_triplet *T_SS = NULL;
    // cholmod_factor *L = NULL;

    cholmod_common Common, *c;
    
    c = &Common;
    cholmod_l_start(c);

    //alocacao de memoria das struturas do suitesparse
    T_SS = cholmod_l_allocate_triplet(nmed, nvar, nmed*nvar, 0, CHOLMOD_REAL, c);
    A_SS = cholmod_l_allocate_sparse(nmed, nvar, nmed*nvar, 0, 0, 0, CHOLMOD_REAL, c);
    // A_T = cholmod_l_allocate_sparse(nmed, nvar, nmed*nvar, 0, 0, 0, CHOLMOD_REAL, c);
    // G_S = cholmod_l_allocate_sparse(nvar, nvar, nvar*nvar, 0, 0, 0, CHOLMOD_REAL, c);
    
    b_SS = cholmod_l_allocate_dense(nmed, 1, nmed, CHOLMOD_REAL, c);
    b_WLS = cholmod_l_allocate_dense(nvar, 1, nmed, CHOLMOD_REAL, c);
    X_SS = cholmod_l_allocate_dense(nvar, 1, nvar, CHOLMOD_REAL, c);
    // L = cholmod_l_allocate_factor(nmed, c);

    
    //Inicializa o modelo no ponto inicial
    atualiza_Rede(grafo, numeroBarras); //atualiza a condição da rede conforme o estado atual
    atualiza_Modelo(grafo, numeroBarras, nmed, medidas); // atualiza modelo de medição conforme a condição atual da rede
    
    
    //-----> LOOP DO ALGORITIMO WLS
    it=0;
    conv = 0;
    clock_t tIni = clock();

    while((it < maxIT)){
        clock_t t0 = clock();
        
        free(Dx);
                
        //Atualiza a Jacobiana esparsa em formato triplet
        atualiza_H_ss(grafo,numeroBarras,ramos,medidas,nmed,T_SS);
        clock_t WriteMatrix = clock();
        
        //converte a matrix triplet para sparse
        A_SS = cholmod_l_triplet_to_sparse(T_SS, T_SS->nnz, c);
        // A_T = cholmod_l_transpose(A_SS, 2, c);
        // G_S = cholmod_l_ssmult(A_T, A_SS, 1, true, false, c);
        // L = cholmod_l_analyze(G_S, c);
        // cholmod_l_factorize(G_S, L, c);  

        //escreve o vetor Dz no formato Dense
        for(int i=0;i<nmed;i++){
            Dz[i] = (medidas[i].zmed - medidas[i].h)/medidas[i].sigma;
            ((double*)b_SS->x)[i] = Dz[i];
        }
        
        int m1 [2] = {0,1};
        int m2 [2] = {0,1};
        

        // cholmod_l_free_triplet(&T_SS, c);
        // clock_t WriteMatrix = clock();
        // //Solucao via SuiteSparse
        int mtype = 0;

        //=========================================================================
        //Convergência do algoritimo - exporta resultados finais e processa erros grosseiros
        if (conv == 1){
            clock_t tFim = clock();
            double tempoIt = (double)(tFim-tIni)/CLOCKS_PER_SEC;
            
            printf("\n\nTempo: %.4lf  Iteracao:  %d \t|Dx|_inf =  %.17lf \t |Grad|_inf =  %.17lf \n",tempoIt, it,nFx,nGx);
            
            exportaPrioriQR_ss(grafo,regua_comp, nvar, nmed, 1, A_SS,c); //exporta arquivo de texto com distribuição a posteriori

            cholmod_l_finish(c);
            //Processamento de Erros Grosseiros via Análise dos Resíduos
            // residuosNormalizadosQR_ss(rN, bHat, nmed, nvar, Dz, H, medidas, A_SS);
            // residuosNormalizadosQR(rN, bHat, nmed, nvar, Dz, H_rf, NULL,z);
            // clock_t tFim3 = clock();
            // double tempoResiduos2 = (double)(tFim3-tFim)/CLOCKS_PER_SEC;
            // printf("\n\n Tempo calculo residuos normalizados: %.4lf",tempoResiduos2);
            
            
            // saidaEstado(grafo, numeroBarras, it, tempoIt, nFx, nGx);
            free(Dz);
            return conv;
        }

        //************************************************************************
        //SOLUCAO DO SISTEMA LINEAR
        //Rt*Dx = Qt*(z-h(x))
        //************************************************************************
        clock_t tHouse = clock();
        //X_SS = cholmod_solve(CHOLMOD_A, L, b_SS, c);
        X_SS = SuiteSparseQR_C_backslash(SPQR_ORDERING_BEST, SPQR_NO_TOL, A_SS, b_SS, c);
        //X_SS = cholmod_l_solve(CHOLMOD_A, L, b_WLS, c);
        
        
        clock_t tSolve = clock();
        Dx = (double*)X_SS->x;
    
        //Atualiza o vetor x
        // double passo = NewtonMod_calcPasso(alpha0, b, Dx, ponto, nvar, Armijo_c1, grafo,numeroBarras, medidas, nmed, ponto, nvar, z,h, W, regua);
        double passo = 1;
        for(i=0;i<nvar;i++){
            ponto[i] = ponto[i] + passo*Dx[i];
        }
        //Calcula vetor gradiente: Multiplica Ht * W * Dz - salva em b_WLS - Y = AX + bY
        cholmod_l_sdmult(A_SS, 1, one, zero,b_SS,b_WLS,c);

        // nGx = cholmod_norm_dense(b_WLS,0,c);
        nGx = norma_inf((double*)b_WLS->x,nmed);
        nFx = norma_inf(Dx,nvar);
        
        clock_t t1 = clock();
        double tempo_it = (double)(t1-t0)/CLOCKS_PER_SEC;
        
        // printf("\n\nIteracao:  %d \t|Dx|_inf =  %.17lf \t |Grad|_inf =  %.17lf \n",it,nFx,nGx);
        // double tempoHouseholder = (double)(t1-tHouse)/CLOCKS_PER_SEC;
        // printf("\nSolve QR: %lf",tempoHouseholder);
        // double tempoTrataMat = (double)(tHouse-WriteMatrix)/CLOCKS_PER_SEC;
        // printf("\nOperacoes Matrizes: %lf",tempoTrataMat);
        // double tempoMontaH = (double)(WriteMatrix-t0)/CLOCKS_PER_SEC;
        // printf("\nMonta H: %lf",tempoMontaH);
        
        
        //************************************************************************
        //TESTE DE CONVERGÊNCIA
        //
        //************************************************************************
        for (i=0;i<nvar;i++){
            if (cabs(Dx[i]) >= tol)
            {
                conv = 0;
                break;
            }
            else
                conv = 1;
        }
        
        //************************************************************************
        //ATUALIZAÇÃO DO ESTADO DA REDE
        //
        //************************************************************************
        atualiza_estado(grafo, ponto, regua_comp, nvar); //atualiza o estado atual do grafo conforme o vetor x calculado
        atualiza_Rede(grafo, numeroBarras);
        atualiza_Modelo(grafo, numeroBarras, nmed, medidas);
        
        it++;
        printf(".");
    }
    //******************************************************************
    //FIM DO LOOP DO WLS
    //
    //******************************************************************
    printf("\n\nNumero maximo de iteracoes atingido %d  %d \t|Dx|_inf =  %.17lf \t |Grad|_inf =  %.17lf \n",it,nFx,nGx);
        
    atualiza_estado(grafo, ponto, regua_comp, nvar);
    atualiza_Rede(grafo, numeroBarras);
    atualiza_Modelo(grafo, numeroBarras, nmed, medidas);

    cholmod_l_finish(c);

    free(Dz);//free(regua);
    //for (i=0;i<nmed;i++) free(H_rf[i]);
    //free(H_rf);
    return (conv);
}

//------------------------------------------------------------------------------
//
//Método Ortogonal para a solução do Estimador - Gauss Newton via QR factorization Esparso
// Utiliza a informação a priori na estimação
//
int otimiza_Gauss_NewtonQR_MAP_ss(double *z, double **h, GRAFO *grafo, long int numeroBarras, DRAM *ramos, DMED *medidas, long int nvar, long int nmed, double *regua_comp, double *ponto, double tol){
    long int it,r;
    long int NAV = 0;
    double alpha0 = 1;
    double *Dz = NULL;
    double *Dx = NULL;
    double *regua = NULL;
    int i,j, conv = 0;
    double nGx,nFx;
    double *rN = NULL, *bHat = NULL;
    double *x0;
    
    //Create basic scalars
    double minusone[2], one [2], zero [2];
    one [0] = 1;
    one [1] = 0 ;
    minusone [0] = -1;
    minusone [1] = 0 ;
    zero [0] = 0 ;
    zero [1] = 0 ;

    //----------------------------------------------------------------------------
    //ALOCAÇÃO
    //
    //----------------------------------------------------------------------------
    Dz = aloca_vetor(nmed);
    rN = aloca_vetor(nmed);
    bHat = aloca_vetor(nmed);
    x0 = aloca_vetor(nmed);

    //=========================================================================
    int nz = 0;
    //inicializacao de variaveis
    cholmod_sparse *A_SS = NULL;
    cholmod_sparse *H_SS = NULL;
    cholmod_sparse *Q_SS = NULL;

    cholmod_dense *x0_SS = NULL;
    cholmod_dense *ponto_SS = NULL;
    cholmod_dense *b_Pr = NULL;
    cholmod_dense *b0_Pr = NULL;
    cholmod_dense *b_SS = NULL;
    cholmod_dense *b_WLS = NULL;
    cholmod_dense *X_SS = NULL;

    cholmod_triplet *T_SS = NULL;

    cholmod_common Common, *c;
    
    c = &Common;
    cholmod_l_start(c);

    //alocacao de memoria das struturas do suitesparse
    T_SS = cholmod_l_allocate_triplet(nmed, nvar, nmed*nvar, 0, CHOLMOD_REAL, c);
    // H_SS = cholmod_l_allocate_sparse(nmed, nvar, nmed*nvar, 0, 0, 0, CHOLMOD_REAL, c);
    // Q_SS = cholmod_l_allocate_sparse(nmed, nvar, nmed*nvar, 0, 0, 0, CHOLMOD_REAL, c);
    // A_T = cholmod_l_allocate_sparse(nmed, nvar, nmed*nvar, 0, 0, 0, CHOLMOD_REAL, c);
    // G_S = cholmod_l_allocate_sparse(nvar, nvar, nvar*nvar, 0, 0, 0, CHOLMOD_REAL, c);    
    x0_SS = cholmod_l_allocate_dense(nvar, 1, nvar, CHOLMOD_REAL, c);
    ponto_SS = cholmod_l_allocate_dense(nvar, 1, nvar, CHOLMOD_REAL, c);
    b_Pr = cholmod_l_allocate_dense(nvar, 1, nvar, CHOLMOD_REAL, c);
    b0_Pr = cholmod_l_allocate_dense(nvar, 1, nvar, CHOLMOD_REAL, c);
    b_SS = cholmod_l_allocate_dense(nmed + nvar, 1, nmed + nvar, CHOLMOD_REAL, c);
    b_WLS = cholmod_l_allocate_dense(nvar, 1, nvar, CHOLMOD_REAL, c);
    X_SS = cholmod_l_allocate_dense(nvar, 1, nvar, CHOLMOD_REAL, c);
    // L = cholmod_l_allocate_factor(nmed, c);
    
    //Importa Modelo a Priori para as Variáveis de Estado
    importaPrioriQR_ss(x0, &Q_SS, c);
    for(i=0;i<nvar;i++) {
        ((double*)x0_SS->x)[i] = x0[i];
        ((double*)ponto_SS->x)[i] = ponto[i];
    }
    free(x0);
    // cholmod_l_print_sparse(Q_SS,"Q",c);
    
    //Multiplica Q * x0 para montar o vetor b
    cholmod_l_sdmult(Q_SS, 0, one, zero,x0_SS,b0_Pr,c);
        
    //Inicializa o modelo no ponto inicial
    atualiza_Rede(grafo, numeroBarras); //atualiza a condição da rede conforme o estado atual
    atualiza_Modelo(grafo, numeroBarras, nmed, medidas); // atualiza modelo de medição conforme a condição atual da rede
    
    //-----> LOOP DO ALGORITIMO WLS
    it=0;
    conv = 0;
    clock_t tIni = clock();

    while((it < maxIT)){
        clock_t t0 = clock();
        
        free(Dx);
        // cholmod_l_free_dense(&X_SS,c);
        // cholmod_l_free_sparse(&H_SS,c);
        // cholmod_l_free_sparse(&A_SS,c);
        
        //Atualiza a Jacobiana esparsa em formato triplet
        atualiza_H_ss(grafo,numeroBarras,ramos,medidas,nmed,T_SS);
        clock_t WriteMatrix = clock();
        
        //converte a matrix triplet para sparse
        H_SS = cholmod_l_triplet_to_sparse(T_SS, T_SS->nnz, c);
        
        
        // A_T = cholmod_l_transpose(A_SS, 2, c);
        // G_S = cholmod_l_ssmult(A_T, A_SS, 1, true, false, c);
        // L = cholmod_l_analyze(G_S, c);
        // cholmod_l_factorize(G_S, L, c);  

        // cholmod_l_print_sparse(Q_SS,"Q", c);
        // cholmod_l_print_sparse(H_SS,"H", c);

        //Concatena a matriz Jacobiana com a matriz R a priori
        A_SS = cholmod_l_vertcat(Q_SS,H_SS,1,c);
        // cholmod_l_print_sparse(A_SS,"A", c);
        
        
        //escreve o vetor b no formato Dense
        cholmod_l_sdmult(Q_SS, 0, one, zero, ponto_SS,b_Pr,c);
        for(int i=0;i<nvar;i++){
            ((double*)b_SS->x)[i] = - (((double*)b_Pr->x)[i] - ((double*)b0_Pr->x)[i]);
        }
        for(int i=0;i<nmed;i++){
            Dz[i] = (medidas[i].zmed - medidas[i].h)/medidas[i].sigma;
            ((double*)b_SS->x)[i + nvar] = Dz[i];
        }
        int m1 [2] = {0,1};
        int m2 [2] = {0,1};
        
        
        // clock_t WriteMatrix = clock();
        // //Solucao via SuiteSparse
        int mtype = 0;

        //=========================================================================
        //Convergência do algoritimo - exporta resultados finais e processa erros grosseiros
        if (conv == 1){
            clock_t tFim = clock();
            double tempoIt = (double)(tFim-tIni)/CLOCKS_PER_SEC;
            
            printf("\n\nTempo: %.4lf  Iteracao:  %d \t|Dx|_inf =  %.17lf \t |Grad|_inf =  %.17lf \n",tempoIt, it,nFx,nGx);
            
            // exportaPrioriQR_ss(grafo,regua_comp, nvar, nmed, 1, A_SS,c); //exporta arquivo de texto com distribuição a posteriori

            cholmod_l_finish(c);
            //Processamento de Erros Grosseiros via Análise dos Resíduos
            // residuosNormalizadosQR_ss(rN, bHat, nmed, nvar, Dz, H, medidas, A_SS);
            // residuosNormalizadosQR(rN, bHat, nmed, nvar, Dz, H_rf, NULL,z);
            // clock_t tFim3 = clock();
            // double tempoResiduos2 = (double)(tFim3-tFim)/CLOCKS_PER_SEC;
            // printf("\n\n Tempo calculo residuos normalizados: %.4lf",tempoResiduos2);
            
            
            // saidaEstado(grafo, numeroBarras, it, tempoIt, nFx, nGx);
            free(Dz);
            return conv;
        }
        
        //************************************************************************
        //SOLUCAO DO SISTEMA LINEAR
        //Rt*Dx = Qt*(z-h(x))
        //************************************************************************
        clock_t tHouse = clock();
        //X_SS = cholmod_solve(CHOLMOD_A, L, b_SS, c);
        X_SS = SuiteSparseQR_C_backslash(SPQR_ORDERING_BEST, SPQR_NO_TOL, A_SS, b_SS, c);

        // SuiteSparseQR_C_factorization *QR;
        // cholmod_dense *Y ;
        // if (it==0) QR = SuiteSparseQR_C_factorize(SPQR_ORDERING_BEST, SPQR_NO_TOL, A_SS,c); 
        // // cholmod_l_print_dense(X_SS,"X",c);
        // Y = SuiteSparseQR_C_qmult(SPQR_QTX,QR,b_SS,c);
        // X_SS = SuiteSparseQR_C_solve(SPQR_RETX_EQUALS_B, QR, Y, c);
        
        //X_SS = cholmod_l_solve(CHOLMOD_A, L, b_WLS, c);
        
        clock_t tSolve = clock();
        Dx = (double*)X_SS->x;

        //Atualiza o vetor x
        // double passo = NewtonMod_calcPasso(alpha0, b, Dx, ponto, nvar, Armijo_c1, grafo,numeroBarras, medidas, nmed, ponto, nvar, z,h, W, regua);
        double passo = 1;
        for(i=0;i<nvar;i++){
            ponto[i] = ponto[i] + passo*Dx[i];
            ((double*)ponto_SS->x)[i] = ponto[i];
        }
        //Calcula vetor gradiente: Multiplica Ht * W * Dz - salva em b_WLS - Y = AX + bY
        cholmod_l_sdmult(A_SS, 1, one, zero,b_SS,b_WLS,c);

        // nGx = cholmod_norm_dense(b_WLS,0,c);
        nGx = norma_inf((double*)b_WLS->x,nvar);
        nFx = norma_inf(Dx,nvar);
        
        clock_t t1 = clock();
        double tempo_it = (double)(t1-t0)/CLOCKS_PER_SEC;
        
        
        printf("\n\nIteracao:  %d \t|Dx|_inf =  %.17lf \t |Grad|_inf =  %.17lf \n",it,nFx,nGx);
        double tempoHouseholder = (double)(t1-tHouse)/CLOCKS_PER_SEC;
        printf("\nSolve QR: %lf",tempoHouseholder);
        double tempoTrataMat = (double)(tHouse-WriteMatrix)/CLOCKS_PER_SEC;
        printf("\nOperacoes Matrizes: %lf",tempoTrataMat);
        double tempoMontaH = (double)(WriteMatrix-t0)/CLOCKS_PER_SEC;
        printf("\nMonta H: %lf",tempoMontaH);      
          
        
        //************************************************************************
        //TESTE DE CONVERGÊNCIA
        //
        //************************************************************************
        for (i=0;i<nvar;i++){
            if (cabs(Dx[i]) >= tol)
            {
                conv = 0;
                break;
            }
            else
                conv = 1;
        }
        // conv = 1;
        
        //************************************************************************
        //ATUALIZAÇÃO DO ESTADO DA REDE
        //
        //************************************************************************
        atualiza_estado(grafo, ponto, regua_comp, nvar); //atualiza o estado atual do grafo conforme o vetor x calculado
        atualiza_Rede(grafo, numeroBarras);
        atualiza_Modelo(grafo, numeroBarras, nmed, medidas);
        
        it++;
        printf(".");
    }
    //******************************************************************
    //FIM DO LOOP DO WLS
    //
    //******************************************************************
    printf("\n\nNumero maximo de iteracoes atingido %d  %d \t|Dx|_inf =  %.17lf \t |Grad|_inf =  %.17lf \n",it,nFx,nGx);
        
    atualiza_estado(grafo, ponto, regua_comp, nvar);
    atualiza_Rede(grafo, numeroBarras);
    atualiza_Modelo(grafo, numeroBarras, nmed, medidas);

    cholmod_l_finish(c);

    free(Dz);//free(regua);
    //for (i=0;i<nmed;i++) free(H_rf[i]);
    //free(H_rf);
    return (conv);
}


//------------------------------------------------------------------------------
//
// ESTIMADOR COM FUSÃO BAYESIANA
// Versão Bayesian Information Fusion - via Maximum a Posteriori Estimation
//
//------------------------------------------------------------------------------
void estimadorBayesFusion_MAP(GRAFO *grafo, long int numeroBarras, DMED *medidas, long int **numeroMedidas, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos,double Sbase){
    long int nmed,nvar;
    int i,j,k, r;
    int conv,choice;
    double *z = NULL,**h = NULL,***H = NULL,**W = NULL, *x = NULL, *regua = NULL, aux = 0;
    
    printf("Estimador de Estado MAP Trifásico via Bayesian Information Fusion...\n");
    //--------------------------------------------------------------------------
    //Alocação de memória das variáveis do estimador de estado
    nmed = 0;
    for (i = 0; i < 14; i++){ 
        for (j = 0; j < 8; j++){
            nmed = nmed + numeroMedidas[i][j];
        }
    }
    nvar = 0;
    //printf("numero barras: %d\n", numeroBarras);
    for (i = 0; i < numeroBarras; i++){
        switch (grafo[i].fases){
            case 1:
                nvar +=2;
                break;
            case 2:
                nvar +=2;
                break;
            case 3:
                nvar +=2;
                break;
            case 4:
                nvar +=4;
                break;    
            case 5:
                nvar +=4;
                break;    
            case 6:
                nvar +=4;
                break;    
            case 7:
                nvar +=6;
                break;    
        }
    }    
    if ((z = (double *)malloc( (nmed) * sizeof(double)))==NULL){
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor z!!!!");
        exit(1); 
    }
    if ((h = malloc( (nmed) * sizeof(double*)))==NULL){
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor h!!!!");
        exit(1); 
    }
    if ((x = (double *)malloc( (nvar) * sizeof(double)))==NULL){
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor x!!!!");
        exit(1); 
    }
    if ((regua = (double *)malloc( (nvar) * sizeof(double)))==NULL){
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor regua!!!!");
        exit(1); 
    }
    
    
    //--------------------------------------------------------------------------
    // Direcionamento dos ponteiros que compõem a régua, vetor h e matriz H
    //monta a regua (ordenaçao do vetor x) - V e teta conforme o grafo
    j=0;
    for(i=0;i<numeroBarras;i++){
        aux = (double) grafo[i].idNo;
        aux += 0.01;
        switch (grafo[i].fases){
            case 1:
                regua[j] = aux;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 2:
                regua[j] = aux + 0.1;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 3:
                regua[j] = aux + 0.2;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 4:
                regua[j] = aux;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.1;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 5:
                regua[j] = aux;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.2;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 6:
                regua[j] = aux+0.1;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.2;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 7:
                regua[j] = aux;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.1;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.2;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;    
        }
    }
    aux = 0;

    //Tratamento da referência
    long int ref_1, ref_2;
    tratamento_referencia(&ref_1, &ref_2, &alimentadores[0], regua, nvar);
    tira_refs_regua(nvar, ref_1, ref_2, regua); 
    nvar = nvar - (ref_2 - ref_1 +1);  

    //vetor h aponta para a estrutura de dados das medidas
    for(i=0;i<nmed;i++){
        h[i] = &medidas[i].h;
    }  
    //Matriz H aponta para a estrutura de dados da régua
    int max_nz_H = 0;
    for(i=0;i<nmed;i++){
        for(j=0;j<medidas[i].nvar;j++){
            for(r = 0;r<nvar;r++){
                medidas[i].reguaH_loc[j] = -1;
                if (cabs(medidas[i].reguaH[j]-regua[r]) < EPS){
                    //Atualiza a Matriz H
                    medidas[i].reguaH_loc[j] = r;
                    break;
                }
            }
        }
    }
    
    //--------------------------------------------------------------------------
    //Estimação de Estado    
    monta_z(z,nmed,medidas);
    monta_W(NULL,nmed,medidas);
    //monta_W_cte(W,nmed,medidas);
    // monta_W_Ident(NULL,nmed,medidas);
    
    incializa_vetor_x(grafo, numeroBarras, alimentadores, numeroAlimentadores,x,regua,nvar);
//    incializa_vetor_x_leitura(grafo, numeroBarras, alimentadores, numeroAlimentadores,x,regua,nvar); //Função que le de arquivo externo a condição inicial


    double tol = 0.000001;
    clock_t tIni = clock();
    choice =  WITH_PRIOR; //0
    if (choice) 
        conv = otimiza_Gauss_NewtonQR_MAP_ss(z, h, grafo, numeroBarras, ramos, medidas, nvar, nmed, regua, x, tol);
    else
        conv = otimiza_Gauss_NewtonQR_ss(z, h, grafo, numeroBarras, ramos, medidas, nvar, nmed, regua, x, tol);
        
    
    exportaCasoReferencia(grafo, numeroBarras, Sbase);
    exportaEstado(grafo,regua,nvar);
    imprimeEstado(grafo,numeroBarras);
    
//    atualiza_H_ret(grafo, numeroBarras, ramos, medidas, nmed);
//    exportaPrioriSCADA(grafo,regua, H, W, nvar, nmed, 0);
//    exportaPrioriQR(grafo,regua,medidas, nvar, nmed, 0);

    clock_t t1 = clock();
    double tempoWLS = (double)(t1-tIni)/CLOCKS_PER_SEC;
    printf("\nEstimação Bayesian MAP: %lf\n\n",tempoWLS);

    free(z);free(h);free(x);free(regua);
    free(H);
}
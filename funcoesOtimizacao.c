#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "data_structures.h"
#include "funcoesBadData.h"
#include "funcoesCalculoEletrico.h"
#include "funcoesMatematicas.h"
#include "funcoesOtimizacao.h"
#include "funcoesTopologia.h"
#include "funcoesWLS.h"

// #include "mmio.h"
#include "SuiteSparseQR_C.h"
#include <cholmod.h>

int chamaFuncao;
int chamaGrad;
int chamaHessiana;
long int optionHachtel; // Fala se é Hachtel ou Normal

void atualiza_x(double *ponto, double *mu, double *x, long int nmed,
                long int nvar) {
  long int i;

  if (optionHachtel == 0) {
    for (i = 0; i < nvar; i++) {
      ponto[i] = ponto[i] + x[i];
    }
  } else if (optionHachtel == 1) {
    for (i = 0; i < nmed; i++) {
      mu[i] = mu[i] + x[i];
    }
    for (i = 0; i < nvar; i++) {
      ponto[i] = ponto[i] + x[i + nmed];
    }
  }
}

double valFun(double *z, double **h, double **W, long int nmed) {
  chamaFuncao++;
  double soma;
  double parcela1;
  double *Dz = aloca_vetor(nmed);
  int i;

  if (optionHachtel == 1) {
    //        for(i=0;i<nmed;i++) r[i] = 0;
    //        sparseMultMatVet(Cov, mu, r);
    //
    //        for(i=0;i<nmed;i++){
    //            aux[i] = r[i] - z[i] + *h[i];
    //        }
    //
    //        parcela1 = prod_escalar(r,mu,nmed);
    //        parcela2 = 0;
    //        parcela3 = prod_escalar(mu,aux,nmed);
    //
    //
    //        soma = 0.5 * parcela1 - parcela2 - parcela3;
    //        free(r);
    //        free(aux);
  } else if (optionHachtel == 0) {
    for (i = 0; i < nmed; i++)
      Dz[i] = (z[i] - *h[i]) / sqrt(W[i][i]);

    parcela1 = prod_escalar(Dz, Dz, nmed);

    soma = 0.5 * parcela1;
  }
  free(Dz);
  return (soma);
}

// Calcula o precodicionador via Cholesky Incompleto para a matriz de Hachtel
void IncCholeskyHachtelPCG(CCSPARSE *Hcc, CCSPARSE *Hcct, double **W,
                           CCSPARSE *L, int nmed, int nvar) {
  int i, j, k;
  double valor;
  SPARSE *tmp;

  // Para o cálculo da diagonal da matriz ganho
  double diagG[nvar];
  for (i = 0; i < nvar; i++)
    diagG[i] = 0;
  for (j = 0; j < nvar; j++) {
    tmp = Hcc[j].val;
    while (tmp != NULL) {
      diagG[j] += tmp->ij[0] * 1 / W[tmp->j][tmp->j] * tmp->ij[0];
      tmp = tmp->prox;
    }
  }

  // Montagem da matriz Lower da fatoração Cholesky
  for (j = 1; j < nmed; j++) {
    L[j].n++;
    valor = 1 / pow(W[j][j], 0.5);
    sparseAdd(&L[j].val, j, j, &valor);

    tmp = Hcct[j].val;
    while (tmp != NULL) {
      L[j].n++;
      valor = 1 / pow(W[j][j], 0.5) * tmp->x;
      sparseAdd(&L[j].val, tmp->i + nmed, j, &valor);
      tmp = tmp->prox;
    }
  }
  for (j = 1; j < nvar; j++) {
    L[j + nmed].n++;
    valor = pow(diagG[j], 0.5);
    sparseAdd(&L[j].val, j + nmed, j + nmed, &valor);
  }
}

// backtracking (interpolacao quadratica) com Armijo. c1 -> constante de Armijo
double NewtonMod_calcPasso(double alpha0, double *grad, double *dir, double *x,
                           int n, double c1, GRAFO *grafo,
                           long int numeroBarras, DMED *medidas, long int nmed,
                           double *ponto, long int nvar, double *z, double **h,
                           double **W, double *regua) {
  int i = 0;
  int loop = 0;
  double passo = alpha0;
  double aux0 = 0;
  double aux1 = 0;
  double fx = 0;
  double *point;
  double trainspotting;

  atualiza_estado(grafo, ponto, regua, nvar);
  atualiza_Rede(grafo, numeroBarras);
  atualiza_Modelo(grafo, numeroBarras, nmed, medidas);

  fx = valFun(z, h, W, nmed);

  point = (double *)malloc(n * sizeof(double));

  for (i = 0; i < n; i++) {
    aux0 += -grad[i] * dir[i];
  }

  do {
    // Atualiza o vetor x
    for (i = 0; i < n; i++) {
      point[i] = x[i] + passo * dir[i];
    }
    atualiza_estado(grafo, point, regua, nvar);
    atualiza_Rede(grafo, numeroBarras);
    atualiza_Modelo(grafo, numeroBarras, nmed, medidas);
    aux1 = valFun(z, h, W, nmed);
    printf("f(x) = %lf   f(x+1)=%lf   Armijo=%lf \t passo = %1f", fx, aux1,
           (fx - c1 * passo * aux0), passo);
    printf("\n");
    // getchar();
    // getchar();

    if (aux1 <= fx - c1 * passo * aux0) {
      loop = 1;
    }

    else {
      trainspotting =
          ((-aux0 * passo * passo) / (2 * (aux1 - fx - passo * aux0)));

      if (trainspotting <= 0.01) {
        passo = passo / 2;
      } else {
        passo = trainspotting;
      }

      //               passo = (-aux0*passo*passo/(2*(aux1-fx-passo*aux0)));
      // passo = passo/2; //da pra tentar usar uma redução não interpolada tb.
    }
  } while (loop == 0);
  atualiza_estado(grafo, ponto, regua, nvar);
  atualiza_Rede(grafo, numeroBarras);
  atualiza_Modelo(grafo, numeroBarras, nmed, medidas);
  free(point);

  return passo;
}

SPARSE *sparseCatHor(SPARSE *A, long int n, SPARSE *B) {
  SPARSE *atual = A;
  SPARSE *C = NULL;

  while (atual != NULL) {
    sparseAdd(&C, atual->i, atual->j, atual->ij);
    atual = atual->prox;
  }
  atual = B;
  while (atual != NULL) {
    sparseAdd(&C, atual->i, atual->j + n, atual->ij);
    atual = atual->prox;
  }
  return (C);
}
SPARSE *sparseCatVert(SPARSE *A, long int n, SPARSE *B) {
  SPARSE *atual = A;
  SPARSE *C = NULL;
  while (atual != NULL) {
    sparseAdd(&C, atual->i, atual->j, atual->ij);
    atual = atual->prox;
  }
  atual = B;
  while (atual != NULL) {
    sparseAdd(&C, atual->i + n, atual->j, atual->ij);
    atual = atual->prox;
  }
  return (C);
}
/*SPARSE* atualizaHsparse(DMED *medidas,long int nmed, double *regua, long int
nvar){ int i,j,r; SPARSE *Hs = NULL;

    for(i=0;i<nmed;i++){
        for(j=0;j<medidas[i].nvar;j++){
            for(r = 0;r<nvar;r++){
                if (cabs(medidas[i].reguaH[j]-regua[r]) < EPS){
                    //Atualiza a Matriz H
                    if (medidas[i].H[j] != 0) sparseAdd(&Hs, i, r,
medidas[i].H[j]); break;
                }
            }
        }
    }
    return(Hs);
}*/

// Monta a matriz Ganho e o vetor b (serve somente para W diagonal)
void montaGb(double **H, double **W, double **h, double *z, int nvar, int smed,
             double **G, double *b, int **Gsimb) {
  int i, j, k;
  double aux_sum = 0;
  int iG, jG, iHt, jHt, iH, jH;

  // Monta Matriz G
  // Somente para W diagonal
  for (i = 0; i < nvar; i++) {
    for (j = 0; j < nvar; j++) {
      G[i][j] = 0;
      for (k = 0; k < smed; k++) {
        if ((H[k][i] != 0) && (H[k][j] != 0))
          G[i][j] += H[k][i] * W[k][k] * H[k][j];
      }
    }
  }
  //    for(i=0;i<nvar;i++){
  //        for(j=0;j<nvar;j++){
  //            G[i][j] = 0;
  //        }
  //    }
  //    for(i=0;i<nvar*nvar*smed;i++){
  //        if (Gsimb[i][0] == NAN) break;
  //        else{
  //            iG = Gsimb[i][0];
  //            jG = Gsimb[i][1];
  //            iHt = Gsimb[i][2];
  //            jHt = Gsimb[i][3];
  //            iH = Gsimb[i][4];
  //            jH = Gsimb[i][5];
  //            G[iG][jG] += H[iHt][jHt]*W[iHt][iHt]*H[iH][jH];
  //        }
  //    }

  // Monta vetor b = Ht.W.(z-h)
  for (i = 0; i < nvar; i++) {
    b[i] = 0;
    for (k = 0; k < smed; k++) {
      b[i] = b[i] + H[k][i] * W[k][k] * (z[k] - (*h[k]));
    }
  }
}

// Salva em arquivo cada iteração
void salva_sol(FILE *arqout, GRAFO *grafo, long int numeroBarras, double **h,
               double ***H, double **G, double *z, double *x, double *Dx,
               double *regua, int it, long int smed, long int nvar, double fx,
               double nGx, double tempo) {
  // extern double **H;
  // extern double **G;

  int i, j;

  // if(it == 0){
  fprintf(arqout, "\n========================================================");
  //    fprintf(arqout,"\nIteracao %d",it);
  //    fprintf(arqout,"\nMatriz Jacobiana H(x)\n");
  //    for(j=0;j<nvar+3;j++)
  //    {
  //        fprintf(arqout,"%.1f\t",regua[j]);
  //    }
  //
  //    fprintf(arqout,"\n\n");
  //    for(i=0;i<smed;i++)
  //    {
  //        for(j=0;j<nvar+3;j++)
  //        {
  //            fprintf(arqout,"%.6f\t",*H[i][j]);
  //        }
  //        fprintf(arqout,"\n");
  //    }
  //    fprintf(arqout,"\n\n");
  //    for(i=0;i<nvar;i++)
  //    {
  //        for(j=0;j<nvar;j++)
  //        {
  //            fprintf(arqout,"%.15f\t",G[i][j]);
  //        }
  //        fprintf(arqout,"\n");
  //    }
  //
  //    fprintf(arqout,"\n\nVetor z \t\tVetor h(x) \t\t\tVetor Dz\n");
  //    for(i=0;i<smed;i++)
  //    {
  //        fprintf(arqout,"%.7lf\t%.7lf\t%.7f\n",z[i],*h[i],z[i] - *h[i]);
  //    }

  fprintf(arqout, "\n\nVetor x e Dx \n");
  for (i = 0; i < nvar; i++) {
    fprintf(arqout, "%.1f\t%.6f\t%.6f\n", regua[i], x[i], Dx[i]);
  }
  fprintf(arqout, "\nTensoes Nodais: Fase-Terra\n");
  for (i = 0; i < numeroBarras; i++) {
    // Retangulares
    // printf("%d\tVa: %.5lf + j%.5lf\tVb: %.5lf + j%.5lf\tVc: %.5lf +
    // j%.5lf\n",grafo[i].barra->ID,__real__ grafo[i].V[0],__imag__
    // grafo[i].V[0],__real__ grafo[i].V[1],__imag__ grafo[i].V[1],__real__
    // grafo[i].V[2],__imag__ grafo[i].V[2]); Polares
    fprintf(arqout, "%d\t%.5lf\t%.5lf\t%.5lf\t %.3lf\t%.3lf\t%.3lf\n",
            grafo[i].barra->ID, cabs(grafo[i].V[0]), cabs(grafo[i].V[1]),
            cabs(grafo[i].V[2]), carg(grafo[i].V[0]) * 180 / PI,
            carg(grafo[i].V[1]) * 180 / PI, carg(grafo[i].V[2]) * 180 / PI);
  }
  fprintf(arqout, "\nTensoes Nodais: Fase-Fase\n");
  for (i = 0; i < numeroBarras; i++) {
    // Retangulares
    // printf("%d\tVa: %.5lf + j%.5lf\tVb: %.5lf + j%.5lf\tVc: %.5lf +
    // j%.5lf\n",grafo[i].barra->ID,__real__ grafo[i].V[0],__imag__
    // grafo[i].V[0],__real__ grafo[i].V[1],__imag__ grafo[i].V[1],__real__
    // grafo[i].V[2],__imag__ grafo[i].V[2]); Polares
    fprintf(arqout, "%d\t%.5lf\t%.5lf\t%.5lf\t %.3lf\t%.3lf\t%.3lf\n",
            grafo[i].barra->ID, cabs(grafo[i].V[0] - grafo[i].V[1]),
            cabs(grafo[i].V[1] - grafo[i].V[2]),
            cabs(grafo[i].V[2] - grafo[i].V[0]),
            carg(grafo[i].V[0] - grafo[i].V[1]) * 180 / PI,
            carg(grafo[i].V[1] - grafo[i].V[2]) * 180 / PI,
            carg(grafo[i].V[2] - grafo[i].V[0]) * 180 / PI);
  }
}

//------------------------------------------------------------------------------
//
// Método Clássico para a solução do Estimador WLS - Gauss Newton
//
int otimiza_Gauss_Newton(double *z, double **h, double ***H, double **W,
                         GRAFO *grafo, long int numeroBarras, DRAM *ramos,
                         DMED *medidas, long int nvar, long int nmed,
                         double *regua_comp, double *ponto, double tol,
                         long int ref1, long int ref2, int **Gsimb) {
  long int it;
  long int NAV = 0, ref;
  double alpha0 = 1;

  double *Dz;
  double *b;
  double *Dx = NULL;
  double **H_rf;
  double *regua;
  double **Ht;
  double **G;
  double **Mtmp;
  int i, j, conv;
  double nGx, nFx;

  optionHachtel = 0;

  FILE *arquivo;
  arquivo = fopen("iteracoesWLS.rtf", "w");

  //    FILE *arqout;
  //    arqout = fopen("ganho.txt","w");

  //----------------------------------------------------------------------------
  //
  // ALOCAÇÃO
  //
  //----------------------------------------------------------------------------
  // free(H_rf);free(Ht);free(Mtmp);free(G);free(Dz);free(b);
  H_rf = aloca_matriz(nmed, nvar);
  Ht = aloca_matriz(nvar, nmed);
  Mtmp = aloca_matriz(nvar, nmed);
  G = aloca_matriz(nvar, nvar);
  Dz = aloca_vetor(nmed);
  b = aloca_vetor(nvar);

  // Preserva a régua
  regua = aloca_vetor(nvar + 3);
  for (i = 0; i < nvar + 3; i++)
    regua[i] = regua_comp[i];

  // atualiza_h(grafo, numeroBarras, nmed, medidas);
  atualiza_Rede(grafo, numeroBarras);
  atualiza_Modelo(grafo, numeroBarras, nmed, medidas);

  //    fprintf(arquivo,"\n\nVetor z \t\tVetor h(x) \t\t\tVetor Dz\n");
  //    for(i=0;i<nmed;i++)
  //    {
  //        fprintf(arquivo,"%.12lf\t
  //        %.12lf\t%.12lf\t%.12f\n",z[i],medidas[i].sigma,*h[i],z[i] - *h[i]);
  //    }

  //-----> LOOP DO ALGORITIMO WLS
  it = 0;
  conv = 0;
  clock_t tIni = clock();

  while ((it < 30)) {
    clock_t t0 = clock();
    //************************************************************************
    // MONTAGEM DE H(x)
    //
    //************************************************************************
    atualiza_H(grafo, numeroBarras, ramos, medidas,
               nmed); // monta_H(2);
                      //         if (it <= 1) monta_W_Ident(W,nmed,medidas);
                      //         else monta_W(W,nmed,medidas);

    clock_t tMontaH = clock();
    // Tira a coluna de angulo da referencia na matriz H(x)
    if (NAV == 0) {
      // tira_refs(H,nvar+3,nmed,12,14,H_rf,regua,ponto,it); // IEEE4
      // tira_refs(H,nvar+3,nmed,35,37,H_rf,regua,ponto,it); //IEEE13
      // tira_refs(H,nvar+3,nmed,92,94,H_rf,regua,ponto,it); //IEEE34

      // tira_refs(H,nvar+3,nmed,6,8,H_rf,regua,ponto,it); //Subestação

      tira_refs(H, nvar + 3, nmed, ref1, ref2, H_rf, regua, ponto, it); // IEEE4
    } else {
      //            if (it>=1){
      //                tira_refs(H,nvar+3,nmed,6,6,H_rf,regua,ponto,it);
      //                //Subestação nvar = nvar + 2;
      //            }
      //            else{
      //                tira_refs(H,nvar+3,nmed,6,8,H_rf,regua,ponto,it);
      //                //Subestação
      //            }
      mat_ig(H, nmed, nvar, H_rf);
    }

    //************************************************************************
    // MONTAGEM DE G(x)
    //
    //************************************************************************
    free(Dx);

    montaGb(H_rf, W, h, z, nvar, nmed, G, b, Gsimb);

    //        fprintf(arqout,"Ganho it: %d nvar = %d\n",it,nvar);
    //        for(i=0;i<nvar;i++){
    //            for(j=0;j<nvar;j++) {
    //                fprintf(arqout,"%.15f\t",G[i][j]);
    //            }
    //            fprintf(arqout,"\n");
    //        }

    //        double eig= eigenvalue_largest(G, nvar, nvar);
    //
    //        printf("\n\n eigen = %.12lf",eig);
    //        getchar();
    //        getchar();

    // fimp_vet(arqmat, b,nvar);
    // fclose(arqmat);
    // chamaFuncao++;
    // chamaGrad++;
    // chamaHessiana++;

    //=========================================================================
    // Algoritimo convergiu
    if (conv == 1) {
      // relatorio_conv(1,tol,tIni);

      clock_t tFim = clock();
      double tempoIt = (double)(tFim - tIni) / CLOCKS_PER_SEC;

      fprintf(arquivo,
              "\n\nIteracao:  %d \t|Dx|_inf =  %.7lf \t |Grad|_inf =  %.7lf \n",
              it, nFx, nGx);

      printf("\n\n Convergência em %d iteracoes e tempo: %.4lf", it, tempoIt);
      fprintf(arquivo, "\n\n Convergência em %d iteracoes e tempo: %.4lf", it,
              tempoIt);
      salva_sol(arquivo, grafo, numeroBarras, h, H, G, z, ponto, Dx, regua, it,
                nmed, nvar, nFx, nGx, tempoIt);

      //            fprintf(arquivo,"\n\nCONVERGENCIA em %d iteracoes\n\n",it);
      //            fprintf(arquivo,"\nTensoes Nodais: Fase-Terra\n");
      //            for(i=0; i<numeroBarras; i++){
      //                //Retangulares
      //                //printf("%d\tVa: %.5lf + j%.5lf\tVb: %.5lf +
      //                j%.5lf\tVc: %.5lf +
      //                j%.5lf\n",grafo[i].barra->ID,__real__
      //                grafo[i].V[0],__imag__ grafo[i].V[0],__real__
      //                grafo[i].V[1],__imag__ grafo[i].V[1],__real__
      //                grafo[i].V[2],__imag__ grafo[i].V[2]);
      //                //Polares
      //                fprintf(arquivo,"%d\t%.5lf\t%.5lf\t%.5lf\t
      //                %.3lf\t%.3lf\t%.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[0]),cabs(grafo[i].V[1]),cabs(grafo[i].V[2]),carg(grafo[i].V[0])*180/PI,carg(grafo[i].V[1])*180/PI,carg(grafo[i].V[2])*180/PI);
      //            }
      /*fprintf(arquivo,"\nTensoes Nodais Convergido (V):\n");
      for(i=0; i<numeroBarras; i++){
          //Polares
          fprintf(arquivo,"%d\tVa: %.5lf | %.3lf \tVb: %.5lf | %.3lf\tVc: %.5lf
      |
      %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[0])*grafo[i].Vbase,carg(grafo[i].V[0])*180/PI,cabs(grafo[i].V[1])*grafo[i].Vbase,carg(grafo[i].V[1])*180/PI,cabs(grafo[i].V[2])*grafo[i].Vbase,carg(grafo[i].V[2])*180/PI);
      }*/

      fprintf(arquivo, "\n\nVetor z \t\tVetor h(x) \t\t\tVetor Dz\n");
      for (i = 0; i < nmed; i++) {
        fprintf(arquivo, "%.7lf\t%.7lf\t%.7f\n", z[i], *h[i], z[i] - *h[i]);
      }
      saidaEstado(grafo, numeroBarras, it, tempoIt, nFx, nGx);
      free(Dz);
      free(b);
      free(H_rf);
      free(Ht);
      free(G);
      free(Mtmp);
      free(regua);
      fclose(arquivo);
      //            fclose(arqout);
      return conv;
    }
    //=========================================================================

    //************************************************************************
    // SOLUCAO DO SISTEMA LINEAR: Fatores Trinagulares
    // G(x)*Dx = b
    //************************************************************************
    clock_t tHouse = clock();
    //        Dx = solve_Crout(G,nvar,nvar,b);
    Dx = solve_Householder(G, nvar, nvar, b);

    nGx = norma_inf(b, nvar);
    nFx = norma_inf(Dx, nvar);
    clock_t t1 = clock();
    double tempo_it = (double)(t1 - t0) / CLOCKS_PER_SEC;

    //        double Armijo_c1 = 0.01;
    //        double passo = NewtonMod_calcPasso(alpha0, b, Dx, ponto, nvar,
    //        Armijo_c1, grafo,numeroBarras, medidas, nmed, ponto, nvar, z,h, W,
    //        regua); double passo =1; if (it>=4) passo = 0.8; if (it>=8) passo
    //        = 0.0000010; for(i=0;i<nvar;i++) Dx[i] = passo*Dx[i];

    //        salva_sol(arquivo,grafo, numeroBarras, h, H, G, z, ponto, Dx,
    //        regua, it, nmed,nvar, nFx, nGx, 0); fprintf(arquivo,"\n\nVetor b
    //        \n"); for(i=0;i<nvar;i++){
    //            fprintf(arquivo,"%.15lf\n",b[i]);
    //        }
    //
    fprintf(arquivo,
            "\n\nIteracao:  %d \t|Dx|_inf =  %.7lf \t |Grad|_inf =  %.7lf \n",
            it, nFx, nGx);
    printf("\n\nIteracao:  %d \t|Dx|_inf =  %.17lf \t |Grad|_inf =  %.17lf \n",
           it, nFx, nGx);

    double tempoHouseholder = (double)(t1 - tHouse) / CLOCKS_PER_SEC;
    printf("\nSolve Householder: %lf", tempoHouseholder);
    double tempoTrataMat = (double)(tHouse - tMontaH) / CLOCKS_PER_SEC;
    printf("\nOperacoes Matrizes: %lf", tempoTrataMat);
    double tempoMontaH = (double)(tMontaH - t0) / CLOCKS_PER_SEC;
    printf("\nMonta H: %lf", tempoMontaH);

    for (i = 0; i < nvar; i++) {
      ponto[i] = ponto[i] + Dx[i];
    }

    //************************************************************************
    // TESTE DE CONVERGÊNCIA
    //
    //************************************************************************
    for (i = 0; i < nvar; i++) {
      if (cabs(Dx[i]) >= tol) {
        conv = 0;
        break;
      } else
        conv = 1;
    }
    // Convergência no gradiente
    //         if (nGx >= 0.001) conv = 0;
    //         if (it <= 3) conv = 0;

    //************************************************************************
    // ATUALIZAÇÃO DO VETOR X
    //
    //************************************************************************
    // atualiza_x(ponto);
    atualiza_estado(grafo, ponto, regua, nvar);
    //        atualiza_h(grafo, numeroBarras, nmed, medidas);
    atualiza_Rede(grafo, numeroBarras);
    atualiza_Modelo(grafo, numeroBarras, nmed, medidas);

    it++;
    printf(".");
  }
  //******************************************************************
  // FIM DO LOOP DO WLS
  //
  //******************************************************************
  // atualiza_x(ponto);
  atualiza_estado(grafo, ponto, regua, nvar);
  // atualiza_h(grafo, numeroBarras, nmed, medidas);
  atualiza_Rede(grafo, numeroBarras);
  atualiza_Modelo(grafo, numeroBarras, nmed, medidas);

  // relatorio_conv(0,tol,tIni);
  //     fclose(arquivo);
  free(Dz);
  free(b);
  free(H_rf);
  free(Ht);
  free(G);
  free(Mtmp);
  free(regua);
  return (conv);
}

//------------------------------------------------------------------------------
//
// Método Clássico para a solução do Estimador WLS - Gauss Newton via QR
// factorization
//
int otimiza_Gauss_NewtonQR(double *z, double **h, double ***H, GRAFO *grafo,
                           long int numeroBarras, DRAM *ramos, DMED *medidas,
                           long int nvar, long int nmed, double *regua_comp,
                           double *ponto, double tol, long int ref1,
                           long int ref2) {
  long int it, r;
  long int NAV = 0;
  double alpha0 = 1;

  double *Dz = NULL;
  double *b = NULL;
  double *Dx = NULL;
  double **H_rf = NULL;
  double *regua = NULL;
  int i, j, conv = 0;
  double nGx, nFx;

  double *rN = NULL, *bHat = NULL;

  optionHachtel = 0;

  FILE *arquivo;
  arquivo = fopen("iteracoesWLS.rtf", "wr+");
  // printf("nmed: %d\n", nmed);
  // printf("nvar: %d\n", nvar);

  //----------------------------------------------------------------------------
  // ALOCAÇÃO
  //
  //----------------------------------------------------------------------------
  // H_rf = aloca_matriz(nmed,nvar);
  Dz = aloca_vetor(nmed);
  b = aloca_vetor(nvar);
  rN = aloca_vetor(nmed);
  bHat = aloca_vetor(nmed);

  // Preserva a régua (ordenação do vetor de estado)
  // regua = aloca_vetor(nvar+3);
  // for (i=0;i<nvar+3;i++) regua[i] = regua_comp[i];

  // Inicializa o modelo no ponto inicial
  atualiza_Rede(
      grafo,
      numeroBarras); // atualiza a condição da rede conforme o estado atual
  atualiza_Modelo(
      grafo, numeroBarras, nmed,
      medidas); // atualiza modelo de medição conforme a condição atual da rede

  for (i = 0; i < nmed; i++) {
    fprintf(arquivo, "%.7lf\t%.7lf\t%.7f\n", z[i], medidas[i].h,
            z[i] - medidas[i].h);
  }

  //-----> LOOP DO ALGORITIMO WLS
  it = 0;
  conv = 0;
  clock_t tIni = clock();

  // inicializa vetores esparsos
  long int *i_sparse = NULL;
  long int *j_sparse = NULL;
  double *x_sparse = NULL;

  int nzeros = 0;
  i_sparse = aloca_vetor(nmed);
  j_sparse = aloca_vetor(nvar);
  x_sparse = aloca_vetor(nmed * nvar);
  int escrito = 0;
  // atualiza_H(grafo, numeroBarras, ramos, medidas, nmed);
  while ((it < 30)) {
    clock_t t0 = clock();
    //************************************************************************
    // MONTAGEM DE H(x)
    //
    //************************************************************************
    // TODO: montar nova atualiza_H_SS

    atualiza_H(grafo, numeroBarras, ramos, medidas,
               nmed); // atualiza Jacobiana do modelo de medição de acordo com o
                      // estado atual

    // printf("H update\n");
    // Multiplica R_1/2*H e o modelo de medição - formulação do estimador via
    // método QR
    for (i = 0; i < nmed; i++) {
      for (j = 0; j < medidas[i].nvar; j++) {
        medidas[i].H[j] = medidas[i].H[j] / medidas[i].sigma;
      }
      Dz[i] = (medidas[i].zmed - medidas[i].h) / medidas[i].sigma;
    }

    // Tira a coluna de angulos da referencia na matriz H(x)
    if (NAV == 0) {
      // tira_refs(H,nvar+3,nmed,ref1,ref2,H_rf,regua,ponto,it);

      // implementa todos os loops de montagem da matrix H (e T_esparsa) em uma
      // mesma função nzeros =
      // tira_refs_sparse(medidas,H,nvar+3,nmed,ref1,ref2,i_sparse, j_sparse,
      // x_sparse,regua,ponto,it);
    } else {
      // mat_ig(H,nmed,nvar,H_rf);
      // nzeros = mat_ig_sparse(H, nmed, nvar, i_sparse, j_sparse, x_sparse);
    }
    clock_t tMontaH = clock();
    free(Dx);

    //=========================================================================
    // Convergência do algoritimo - exporta resultados finais e processa erros
    // grosseiros
    if (conv == 1) {
      clock_t tFim = clock();
      double tempoIt = (double)(tFim - tIni) / CLOCKS_PER_SEC;
      printf(
          "\n\nIteracao:  %d \t|Dx|_inf =  %.17lf \t |Grad|_inf =  %.17lf \n",
          it, nFx, nGx);
      printf("\n\n Convergência em %d iteracoes e tempo: %.4lf", it, tempoIt);

      printf("\n\n Nmed: %d e Nvar: %d \n\n", nmed, nvar);
      // Processamento de Erros Grosseiros via Análise dos Resíduos
      //             residuosNormalizadosQR(rN, bHat, nmed, nvar, Dz, H_rf,
      //             NULL,z); clock_t tFim3 = clock(); double tempoResiduos2 =
      //             (double)(tFim3-tFim)/CLOCKS_PER_SEC; printf("\n\n Tempo
      //             calculo residuos normalizados: %.4lf",tempoResiduos2);

      fprintf(arquivo, "\n\n Convergência em %d iteracoes e tempo: %.4lf", it,
              tempoIt);
      fprintf(arquivo, "\n\nVetor z \t\tVetor h(x) \t\t\tVetor Dz\n");
      for (i = 0; i < nmed; i++) {
        fprintf(arquivo, "%.7lf\t%.7lf\t%.7f\n", z[i], *h[i], z[i] - *h[i]);
      }
      saidaEstado(grafo, numeroBarras, it, tempoIt, nFx, nGx);
      free(Dz);
      free(b); // free(regua);
      // for (i=0;i<nmed;i++) free(H_rf[i]);
      // free(H_rf);
      fclose(arquivo);
      return conv;
    }
    //=========================================================================
    // Impressão das Matrizes esparsas para o Suite SParse - GUSTAVO MIRANDA

    // Entender esta parte do código que é onde eu exporto a cada iteração as
    // matrizes e vetores a serem utilizados na solução iterativa do estimador
    // Modificar aqui para utilizar a SPQR

    int nz = 0;

    // inicializacao de variaveis
    // printf("INIT SPARSE\n");
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

    // alocacao de memoria das struturas do suitesparse
    // printf("nmed: %d\n", nmed);
    // printf("nvar: %d\n", nvar);
    // alocacao de memoria das struturas do suitesparse
    T_SS =
        cholmod_l_allocate_triplet(nmed, nvar, nmed * nvar, 0, CHOLMOD_REAL, c);
    A_SS = cholmod_l_allocate_sparse(nmed, nvar, nmed * nvar, 0, 0, 0,
                                     CHOLMOD_REAL, c);
    A_T = cholmod_l_allocate_sparse(nmed, nvar, nmed * nvar, 0, 0, 0,
                                    CHOLMOD_REAL, c);
    G_S = cholmod_l_allocate_sparse(nvar, nvar, nvar * nvar, 0, 0, 0,
                                    CHOLMOD_REAL, c);
    b_SS = cholmod_l_allocate_dense(nmed, 1, nmed, CHOLMOD_REAL, c);
    b_WLS = cholmod_l_allocate_dense(nvar, 1, nmed, CHOLMOD_REAL, c);
    X_SS = cholmod_l_allocate_dense(nvar, 1, nvar, CHOLMOD_REAL, c);
    L = cholmod_l_allocate_factor(nmed, c);

    int index = 0;
    for (int i = 0; i < nmed; i++) {
      for (int r = 0; r < nvar; r++) {
        if (*H[i][r] != 0) {
          ((long int *)T_SS->i)[index] = i;
          ((long int *)T_SS->j)[index] = r;
          ((double *)T_SS->x)[index] = *H[i][r];
          // fprintf(mat,"%ld,%ld,%f\n",i,r,*H[i][r]);

          T_SS->nnz += 1;
          index += 1;
        }
      }
    }

    // escreve o vetor Dz no formato Dense
    for (int i = 0; i < nmed; i++) {
      ((double *)b_SS->x)[i] = Dz[i];
    }

    // converte a matrix triplet para sparse

    A_SS = cholmod_l_triplet_to_sparse(T_SS, nvar * nvar, c);
    // A_T = cholmod_l_transpose(A_SS, 2, c);

    int m1[2] = {0, 1};
    int m2[2] = {0, 1};

    cholmod_l_free_triplet(&T_SS, c);
    clock_t WriteMatrix = clock();
    // //Solucao via SuiteSparse
    int mtype = 0;

    clock_t tHouse = clock();

    X_SS = SuiteSparseQR_C_backslash(SPQR_ORDERING_BEST, SPQR_DEFAULT_TOL, A_SS,
                                     b_SS, c);

    clock_t tSolve = clock();
    Dx = (double *)X_SS->x;

    cholmod_l_finish(c);

    float passo = 1;

    for (i = 0; i < nvar; i++) {
      b[i] = 0;
      ponto[i] = ponto[i] + Dx[i];
      // Dx[i] += passo * Dx[i];
      for (j = 0; j < nmed; j++) {
        b[i] = b[i] + (*H[j][i]) * Dz[j];
      }
    }

    nGx = norma_inf(b, nvar);
    nFx = norma_inf(Dx, nvar);
    clock_t t1 = clock();
    double tempo_it = (double)(t1 - t0) / CLOCKS_PER_SEC;

    fprintf(arquivo,
            "\n\nIteracao:  %d \t|Dx|_inf =  %.7lf \t |Grad|_inf =  %.7lf \n",
            it, nFx, nGx);
    printf("\n\nIteracao:  %d \t|Dx|_inf =  %.17lf \t |Grad|_inf =  %.17lf \n",
           it, nFx, nGx);
    double tempoHouseholder = (double)(t1 - tHouse) / CLOCKS_PER_SEC;
    // printf("\nSolve QR: %lf",tempoHouseholder);
    double tempoTrataMat = (double)(tHouse - tMontaH) / CLOCKS_PER_SEC;
    // printf("\nOperacoes Matrizes: %lf",tempoTrataMat);
    double tempoMontaH = (double)(tMontaH - t0) / CLOCKS_PER_SEC;
    // printf("\nMonta H: %lf",tempoMontaH);

    //************************************************************************
    // TESTE DE CONVERGÊNCIA
    //
    //************************************************************************
    // printf("conv\n");
    for (i = 0; i < nvar; i++) {
      if (cabs(Dx[i]) >= tol) {
        conv = 0;
        break;
      } else
        conv = 1;
    }
    // printf("FIM conv \n");

    //************************************************************************
    // ATUALIZAÇÃO DO ESTADO DA REDE
    //
    //************************************************************************
    atualiza_estado(
        grafo, ponto, regua_comp,
        nvar); // atualiza o estado atual do grafo conforme o vetor x calculado
    // printf("estado ok\n");
    atualiza_Rede(grafo, numeroBarras);
    // printf("rede ok\n");
    atualiza_Modelo(grafo, numeroBarras, nmed, medidas);
    // printf("modelo ok\n");
    atualiza_h(grafo, numeroBarras, nmed, medidas);

    it++;
    printf(".");
  }
  //******************************************************************
  // FIM DO LOOP DO WLS
  //
  //******************************************************************
  fprintf(arquivo,
          "\n\nNumero maximo de iteracoes atingido %d \t|Dx|_inf =  %.7lf \t "
          "|Grad|_inf =  %.7lf \n",
          it, nFx, nGx);
  printf("\n\nNumero maximo de iteracoes atingido %d  %d \t|Dx|_inf =  %.17lf "
         "\t |Grad|_inf =  %.17lf \n",
         it, nFx, nGx);

  atualiza_estado(grafo, ponto, regua_comp, nvar);
  atualiza_Rede(grafo, numeroBarras);
  atualiza_Modelo(grafo, numeroBarras, nmed, medidas);

  free(Dz);
  free(b); // free(regua);
  // for (i=0;i<nmed;i++) free(H_rf[i]);
  // free(H_rf);
  fclose(arquivo);
  return (conv);
}

// Matriz de Hachtel Aumentada - sem virtuais
int otimiza_Gauss_Newton_Hachtel(double *z, double **h, double ***H, double **W,
                                 GRAFO *grafo, long int numeroBarras,
                                 DRAM *ramos, DMED *medidas, long int nvar,
                                 long int nmed, double *regua, double *ponto,
                                 double tol, long int ref1, long int ref2) {
  long int it;
  long int NAV = 0, ref;

  double *Dz;
  double *b;
  double *Dx = NULL;
  double **H_rf;
  double **Ht;
  double **G;
  double **Mtmp;
  double *mu, *b_hach, *zeros_vet;
  double **Hachtel, **Zeros, **Tmp1;
  int i, j, conv;
  double nGx, nFx;

  FILE *arquivo;
  arquivo = fopen("iteracoes.txt", "w");

  //----------------------------------------------------------------------------
  //
  // ALOCAÇÃO
  //
  //----------------------------------------------------------------------------
  // free(H_rf);free(Ht);free(Mtmp);free(G);free(Dz);free(b);
  H_rf = aloca_matriz(nmed, nvar);
  Ht = aloca_matriz(nvar, nmed);
  Mtmp = aloca_matriz(nvar, nmed);
  G = aloca_matriz(nvar, nvar);
  Dz = aloca_vetor(nmed);
  b = aloca_vetor(nvar);

  Zeros = aloca_matriz(nvar, nvar);
  Tmp1 = aloca_matriz(nvar, nvar + nmed);
  Hachtel = aloca_matriz(nvar + nmed, nvar + nmed);
  b_hach = aloca_vetor(nmed + nvar);
  mu = aloca_vetor(nmed);
  zeros_vet = aloca_vetor(nvar);
  // clock_t tIni = clock();

  // atualiza_x(ponto); //ponto inicial já inicializado no corpo do main
  atualiza_h(grafo, numeroBarras, nmed, medidas);

  //    atualiza_H(grafo, numeroBarras, ramos, medidas, nmed);
  //    "tira_refs"(H,nvar+3,nmed,92,94,H_rf,regua,ponto,it); //IEEE34
  //    montaGb(H_rf, W, h, z,nvar, nmed, G,b);
  // tira_refs(H,nvar+3,nmed,12,14,H_rf,regua,ponto,it);

  //-----> LOOP DO ALGORITIMO WLS
  it = 0;
  conv = 0;

  while ((it < 30)) {
    clock_t t0 = clock();
    //************************************************************************
    // MONTAGEM DE H(x)
    //
    //************************************************************************
    atualiza_H(grafo, numeroBarras, ramos, medidas, nmed); // monta_H(2);

    clock_t tMontaH = clock();
    // Tira a coluna de angulo da referencia na matriz H(x)
    if (NAV == 0) {
      // tira_refs(H,nvar+3,nmed,12,14,H_rf,regua,ponto,it); // IEEE4
      // tira_refs(H,nvar+3,nmed,35,37,H_rf,regua,ponto,it); //IEEE13
      // tira_refs(H,nvar+3,nmed,92,94,H_rf,regua,ponto,it); //IEEE34

      // tira_refs(H,nvar+3,nmed,6,8,H_rf,regua,ponto,it); //Subestação

      tira_refs(H, nvar + 3, nmed, ref1, ref2, H_rf, regua, ponto, it); // IEEE4
    } else {
      //            if (it>=1){
      //                tira_refs(H,nvar+3,nmed,6,6,H_rf,regua,ponto,it);
      //                //Subestação nvar = nvar + 2;
      //            }
      //            else{
      //                tira_refs(H,nvar+3,nmed,6,8,H_rf,regua,ponto,it);
      //                //Subestação
      //            }
      mat_ig(H, nmed, nvar, H_rf);
    }

    //************************************************************************
    // MONTAGEM DE G(x)
    //
    //************************************************************************
    free(Dx);

    // montaGb(H_rf, W, h, z,nvar, nmed, G,b);
    for (i = 0; i < nmed; i++)
      Dz[i] = z[i] - *h[i];
    matTransp(H_rf, nmed, nvar, Ht);
    cat_hor(W, nmed, nmed, H_rf, nmed, nvar, Hachtel);
    cat_hor(Ht, nvar, nmed, Zeros, nvar, nvar, Tmp1);
    cat_vert(Hachtel, nmed, nmed + nvar, Tmp1, nvar, nvar + nmed, Hachtel);

    // fimp_vet(arqmat, b,nvar);
    // fclose(arqmat);
    // chamaFuncao++;
    // chamaGrad++;
    // chamaHessiana++;

    //=========================================================================
    // Algoritimo convergiu
    if (conv == 1) {
      // relatorio_conv(1,tol,tIni);
      printf("\n\nCONVERGENCIA em %d iteracoes\n\n", it);
      fprintf(arquivo, "\n\nCONVERGENCIA em %d iteracoes\n\n", it);
      fprintf(arquivo, "\nTensoes Nodais: Fase-Terra\n");
      for (i = 0; i < numeroBarras; i++) {
        // Retangulares
        // printf("%d\tVa: %.5lf + j%.5lf\tVb: %.5lf + j%.5lf\tVc: %.5lf +
        // j%.5lf\n",grafo[i].barra->ID,__real__ grafo[i].V[0],__imag__
        // grafo[i].V[0],__real__ grafo[i].V[1],__imag__ grafo[i].V[1],__real__
        // grafo[i].V[2],__imag__ grafo[i].V[2]); Polares
        fprintf(arquivo, "%d\t%.5lf\t%.5lf\t%.5lf\t %.3lf\t%.3lf\t%.3lf\n",
                grafo[i].barra->ID, cabs(grafo[i].V[0]), cabs(grafo[i].V[1]),
                cabs(grafo[i].V[2]), carg(grafo[i].V[0]) * 180 / PI,
                carg(grafo[i].V[1]) * 180 / PI, carg(grafo[i].V[2]) * 180 / PI);
      }
      /*fprintf(arquivo,"\nTensoes Nodais Convergido (V):\n");
      for(i=0; i<numeroBarras; i++){
          //Polares
          fprintf(arquivo,"%d\tVa: %.5lf | %.3lf \tVb: %.5lf | %.3lf\tVc: %.5lf
      |
      %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[0])*grafo[i].Vbase,carg(grafo[i].V[0])*180/PI,cabs(grafo[i].V[1])*grafo[i].Vbase,carg(grafo[i].V[1])*180/PI,cabs(grafo[i].V[2])*grafo[i].Vbase,carg(grafo[i].V[2])*180/PI);
      }*/

      fprintf(arquivo, "\n\nVetor z \t\tVetor h(x) \t\t\tVetor Dz\n");
      for (i = 0; i < nmed; i++) {
        fprintf(arquivo, "%.7lf\t%.7lf\t%.7f\n", z[i], *h[i], z[i] - *h[i]);
      }
      return conv;
    }
    //=========================================================================

    //************************************************************************
    // SOLUCAO DO SISTEMA LINEAR: Fatores Trinagulares
    // G(x)*Dx = b
    //************************************************************************
    cat_vert_vet(Dz, nmed, zeros_vet, nvar, b_hach);
    clock_t tHouse = clock();
    Dx = solve_Householder(Hachtel, nvar + nmed, nvar + nmed, b_hach);

    nGx = norma_euc(b, nvar);
    nFx = norma_inf(Dx, nvar);
    clock_t t1 = clock();
    double tempo_it = (double)(t1 - t0) / CLOCKS_PER_SEC;

    double tempoHouseholder = (double)(t1 - tHouse) / CLOCKS_PER_SEC;
    printf("\nSolve Householder: %lf", tempoHouseholder);
    double tempoTrataMat = (double)(tHouse - tMontaH) / CLOCKS_PER_SEC;
    printf("\nOperacoes Matrizes: %lf", tempoTrataMat);
    double tempoMontaH = (double)(tMontaH - t0) / CLOCKS_PER_SEC;
    printf("\nMonta H: %lf", tempoMontaH);

    // Atualiza os mutilicadores de lagrange mu
    for (i = 0; i < nmed; i++) {
      mu[i] = mu[i] + Dx[i];
    }
    // Atualiza o vetor x
    for (i = 0; i < nvar; i++) {
      ponto[i] = ponto[i] + Dx[i + nmed];
    }

    //************************************************************************
    // TESTE DE CONVERGÊNCIA
    //
    //************************************************************************
    for (i = 0; i < nvar; i++) {
      if (cabs(Dx[i]) >= tol) {
        conv = 0;
        break;
      } else
        conv = 1;
    }
    //************************************************************************
    // ATUALIZAÇÃO DO VETOR X
    //
    //************************************************************************
    // atualiza_x(ponto);
    atualiza_estado(grafo, ponto, regua, nvar);
    atualiza_h(grafo, numeroBarras, nmed, medidas);

    it++;
    printf(".");
  }
  //******************************************************************
  // FIM DO LOOP DO WLS
  //
  //******************************************************************
  // atualiza_x(ponto);
  atualiza_estado(grafo, ponto, regua, nvar);
  atualiza_h(grafo, numeroBarras, nmed, medidas);

  // relatorio_conv(0,tol,tIni);
  fclose(arquivo);
  return (conv);
}

// Solução do problema de estimação de estado via gradiente conjugado para Least
// Squares lineares https://math.aalto.fi/opetus/inv/CGalgorithm.pdf int
// otimiza_CGLS_sparseNewton(double *z, double **h, double **W, GRAFO *grafo,
// long int numeroBarras, DRAM *ramos, DMED *medidas, long int nvar, long int
// nmed, double *regua, double *ponto, double tol, long int ref1, long int
// ref2){
//     long int it = 0, conv = 0, i, j,k, NAV = 0;
//     double alpha0 = 1;
//     double Armijo_c1 = 0.01;
//
//     optionHachtel = 0;
//
//     FILE *arquivo;
//     arquivo = fopen("iteracoes.txt","w");
//
//     SPARSE *Ht = NULL, *Hs = NULL,*Cov = NULL,*Tmp = NULL, *M = NULL, *Ws =
//     NULL;
//
//
//     //----------------------------------------------------------------------------
//     //Matriz de covariância esparsa
//     double teste[1];
//     teste[0] = 1;
//     for(i=0;i<nmed;i++){
//         sparseAdd(&Cov, i, i, &W[i][i]);   //Matriz de Covariância ao invés
//         da peso para o Hachtel - atualizar no valFun
//     }
//     // Matriz W esparsa
//     for(i=0;i<nmed;i++){
//         sparseAdd(&Ws, i, i, &W[i][i]);
//     }
//
//     ///----------------------------------------------------------------------------
//     //Alocação de vetores
//     double *Dz = aloca_vetor(nmed);
//     double *Dx = aloca_vetor(nvar);
//     double *dir = aloca_vetor(nvar);
//     double *b = aloca_vetor(nvar);
//     double *regua_comp = aloca_vetor(nvar+3);
//     for(i=0;i<nvar+3;i++) regua_comp[i] = regua[i];
//
//     double *mu  = aloca_vetor(nmed);
//     double *zeros_vet = aloca_vetor(nvar);
//
//     atualiza_h(grafo, numeroBarras, nmed, medidas);
//     atualiza_H(grafo, numeroBarras, ramos, medidas, nmed);
//
//     ///----------------------------------------------------------------------------
//     //Alocação de vetores
//     if (NAV==0){
//         //tira_refs(H,nvar+3,nmed,ref1,ref2,H_rf,regua,ponto,it); // IEEE4
//         for (j=0;j<nvar+3;j++)
//             {
//                 if (j<ref1){
//                     regua[j] = regua[j];
//                     ponto[j] = ponto[j];
//                 }
//                 else if (j>ref2){
//                     regua[(j-(ref2-ref1+1))] = regua[j];
//                     ponto[(j-(ref2-ref1+1))] = ponto[j];
//                 }
//             }
//     }
//     ///----------------------------------------------------------------------------
//     //Monta matriz esparsa de Hachtel
//     int r;
//     for(i=0;i<nmed;i++){
//         for(j=0;j<medidas[i].nvar;j++){
//             for(r = 0;r<nvar;r++){
//                 if (cabs(medidas[i].reguaH[j]-regua[r]) < EPS){
//                     //Atualiza a Matriz H
//                     if (medidas[i].H[j] != 0){
//                         sparseAdd(&Hs, i, r, &medidas[i].H[j]);
//                         sparseAdd(&Ht, r, i, &medidas[i].H[j]);
//                     }
//                     break;
//                 }
//             }
//         }
//     }
//     //Multiplica pela ponderação
//     for(i=0;i<nmed;i++){
//         for(j=0;j<medidas[i].nvar;j++){
//             medidas[i].H[j] = pow(W[i][i],0.5) * medidas[i].H[j];
//         }
//     }
//
//     //Encontra a diagonal da matriz ganho para usar no precondicionamento
//     //Jacobi Preconditioner
//     double diagG[nmed];
//     for(i=0;i<nvar;i++) diagG[i] = 0;
//     SPARSE *tmp = Hs;
//     while (tmp != NULL)
//     {
//         diagG[tmp->j] += tmp->ij[0] * tmp->ij[0];
//         tmp = tmp->prox;
//     }
//     for(i=0;i<nvar;i++) diagG[i] = 1/diagG[i];
//
//     FILE *arquivo2;
//     arquivo2 = fopen("CGLS.txt","w");
//
//     SPARSE *atual = Hs;
//     fprintf(arquivo2,"\n\nCGLS its: %d",j);
//     fprintf(arquivo2,"\nMatriz H:");
//     while(atual != NULL)
//     {
//         fprintf(arquivo2,"\n%d \t %d \t
//         %.18lf",atual->i+1,atual->j+1,*atual->ij); atual = atual->prox;
//     }
//     fprintf(arquivo2,"\nMatriz W:");
//     for(i=0;i<nmed;i++) fprintf(arquivo2,"\n%d \t %d \t
//     %.18lf",i+1,i+1,W[i][i]); fclose(arquivo2);
//
//     ///----------------------------------------------------------------------------
//     //Monta matriz H esparsa no formato compressed column
//     CCSPARSE Hcc[nvar];
//     for(i=0;i<nmed;i++){
//         for(j=0;j<medidas[i].nvar;j++){
//             for(r = 0;r<nvar;r++){
//                 if (cabs(medidas[i].reguaH[j]-regua[r]) < EPS){
//                     //Atualiza a Matriz H
//                     if (medidas[i].H[j] != 0){
//                         Hcc[r].n++;
//                         sparseAdd(&Hcc[r].val, i, r, &medidas[i].H[j]);
//                     }
//                     break;
//                 }
//             }
//         }
//     }
//     CCSPARSE L[nvar];
//     sparseCCcholesky(Hcc,L,nvar);
//
//     //Leitura da Lower
//     FILE *arquivo3;
//     arquivo3 = fopen("Lower.txt","r");
//
//     for (i=0;i<nvar;i++){
//         L[i].n = 0;
//         L[i].val = NULL;
//     }
//     sparseCCread(arquivo3, L);
//
//     fclose(arquivo3);
////    //Matrizes esparsas no formato compressed column
////    CCSPARSE Jac[nvar], Ganho[nvar];
////    for (i=0;i<nvar;i++){
////        Jac[i].val = NULL;
////        Ganho[i].val = NULL;
////    }
////    tmp = Hs;
////    while (tmp != NULL){
////        sparseAdd(&Jac[tmp->j].val, tmp->i, tmp->j, tmp->ij);
////        tmp = tmp->prox;
////    }
////    for (i=0;i<nvar;i++){
////        sparsePrint(Jac[i].val);
////    }
//
//    ///----------------------------------------------------------------------------
//    //Inicio do loop do algoritmo
//    double ajuste = 1;
//    double* gradiente_n = aloca_vetor(nvar);
//    double* gradiente = aloca_vetor(nvar);
//
//    clock_t tIni = clock();
//    double norma = 1000;
//    do{
//        clock_t t0 = clock();
//        for (i=0;i<nmed;i++) Dz[i] = pow(W[i][i],0.5)*(z[i] - *h[i]);
//        double valorFuncao = valFun(z,h,mu, Cov, nmed);
//
//        for (i=0;i<nvar;i++) gradiente_n[i] = 0;
//        sparseMultMatVet(Ht, Dz, gradiente_n); //Tem que mudar para colocar a
//        covariância diferente da identidade
//
//        for(i=0;i<nvar;i++) gradiente[i] = -gradiente_n[i];
//
//        double norma = norma_inf(gradiente,nvar);
////        printf("it: %d \t norma = %2f \t f(x) = %2f \n", it,
///norma,valorFuncao);
//
//        if ((norma <= tol)&&(it>=3)){
//            //Algoritmo convergiu
//            conv = 1;
//            clock_t t1 = clock();
//            double tempoIt = (double)(t1-tIni)/CLOCKS_PER_SEC;
//
//            printf("it: %d \t norma = %2f \t f(x) = %2f \n", it,
//            norma,valorFuncao); printf("\n\n Convergência em %d iteracoes e
//            tempo: %.4lf",it,tempoIt);
//
//            //relatorio_conv(conv,tol,tIni);
//            break;
//        }else{
//            //Calcula o proximo ponto
//            //double** hess = hessianaAprox(); //Hachtel já vai estar
//            atualizada
//
//            //double* dir = solve_Householder(hess,nvar,nvar,gradiente_n);
//            for(i=0;i<nvar;i++) Dx[i] = 0;
////            sparseJacobiPrecon(Hs, &M);
////            sparseCGLS_preJacobi(diagG,Hs,Ht,Dz,Dx,nmed,nvar);
//            sparseCGLS_preCholesky(L,Hs,Ht,Dz,Dx,nmed,nvar);
////            sparseCGLS(Hs,Ht,Dz,Dx,nmed,nvar);
//
//            //double passo =
//            NewtonMod_calcPasso(alpha0,gradiente,dir,ponto,nvar,Armijo_c1);
//            //double passo = NewtonMod_calcPasso(alpha0,gradiente,dir,Dx,
//            nvar, Armijo_c1, grafo, numeroBarras, medidas, nmed, mu, ponto,
//            nvar, z,h, Cov, regua);
//            //ajusta o tamanho da pertubacao com base no passo
//            /*if(passo<=0.2){
//                ajuste = ajuste * 5;
//            }
//            if(passo>=0.9){
//                ajuste = ajuste / 5;
//            }*/
//
//            clock_t t1 = clock();
//            double tempoIt = (double)(t1-t0)/CLOCKS_PER_SEC;
//
//            //Atualiza a variável do problema de otimização
//            //for(i=0;i<nvar;i++) Dx[i] = Dx[i] - passo*gradiente[i];
//
//            //norma = norma_inf(Dx,nvar);
//
//            //Atualiza os mutilicadores de lagrange mu
//            atualiza_x(ponto, mu, Dx, nmed, nvar);
//
//            //atualiza o modelo
//            atualiza_estado(grafo, ponto, regua, nvar);
//            atualiza_h(grafo, numeroBarras, nmed, medidas);
//            atualiza_H(grafo, numeroBarras, ramos, medidas, nmed);
//            //Multiplica pela ponderação
//            for(i=0;i<nmed;i++){
//                for(j=0;j<medidas[i].nvar;j++){
//                    medidas[i].H[j] = W[i][i] * medidas[i].H[j];
//                }
//            }
//
////            for(i=0;i<nvar;i++) diagG[i] = 0;
////            SPARSE *tmp = Hs;
////            while (tmp != NULL)
////            {
////                diagG[tmp->j] += tmp->ij[0] * tmp->ij[0];
////                tmp = tmp->prox;
////            }
////            for(i=0;i<nvar;i++) L[i].val->x = pow(diagG[i],0.5);
//
//            it++;
//        }
//    } while( it < 20);
//
//    fclose(arquivo);
//
//    free(Dz);free(Dx);free(dir);free(b);free(regua_comp);
//    free(mu);free(zeros_vet);
//
//    free(gradiente_n);
//    free(gradiente);
//    return conv;
//}

// Matriz de Hachtel Aumentada - sem virtuais
int otimiza_Gauss_Newton_sparseHachtel(double *z, double **h, double **W,
                                       GRAFO *grafo, long int numeroBarras,
                                       DRAM *ramos, DMED *medidas,
                                       long int nvar, long int nmed,
                                       double *regua, double *ponto, double tol,
                                       long int ref1, long int ref2) {
  long int it;
  long int NAV = 0, ref;
  optionHachtel = 1;

  double *Dz;
  double *b;
  double *Dx = NULL;
  double **G;
  double *mu, *b_hach, *zeros_vet, *regua_comp;
  int i, j, conv;
  double nGx, nFx;

  SPARSE *Ht = NULL, *Hs = NULL, *Hachtel = NULL, *Cov = NULL, *Tmp = NULL;

  FILE *arquivo;
  arquivo = fopen("iteracoes.txt", "w");

  //----------------------------------------------------------------------------
  //
  // ALOCAÇÃO
  //
  //----------------------------------------------------------------------------
  // free(H_rf);free(Ht);free(Mtmp);free(G);free(Dz);free(b);
  Dz = aloca_vetor(nmed);
  Dx = aloca_vetor(nmed + nvar);
  b = aloca_vetor(nvar);
  regua_comp = aloca_vetor(nvar + 3);
  for (i = 0; i < nvar + 3; i++)
    regua_comp[i] = regua[i];

  b_hach = aloca_vetor(nmed + nvar);
  mu = aloca_vetor(nmed);
  zeros_vet = aloca_vetor(nvar);

  atualiza_h(grafo, numeroBarras, nmed, medidas);
  atualiza_H(grafo, numeroBarras, ramos, medidas, nmed);

  // Tira a coluna de angulo da referencia na matriz H(x)
  if (NAV == 0) {
    // tira_refs(H,nvar+3,nmed,ref1,ref2,H_rf,regua,ponto,it); // IEEE4
    for (j = 0; j < nvar + 3; j++) {
      if (j < ref1) {
        regua[j] = regua[j];
        ponto[j] = ponto[j];
      } else if (j > ref2) {
        regua[(j - (ref2 - ref1 + 1))] = regua[j];
        ponto[(j - (ref2 - ref1 + 1))] = ponto[j];
      }
    }
  }
  // Para Fluxo de carga manter constante  tensão slack
  for (j = 0; j < nvar; j++) {
    if (j > 2) {
      regua[(j - 3)] = regua[j];
      ponto[(j - 3)] = ponto[j];
    }
  }
  nvar = nvar - 3; // Se comentar tirar do final desta função também

  //----------------------------------------------------------------------------
  // Insere a Matriz de Covariancia na Matriz de hachtel
  for (i = 0; i < nmed; i++) {
    W[i][i] = 1 / W[i][i]; // Na Hachtel usa a matriz de covariância
    sparseAdd(&Hachtel, i, i, &W[i][i]);
  }

  // Monta matriz H e Ht esparsa no formato compressed column
  // Matriz H esparsa
  CCSPARSE Hcc[nvar], Hcct[nmed];
  for (i = 0; i < nmed; i++) {
    Hcct[i].n = 0;
    Hcct[i].val = NULL;
  }
  for (i = 0; i < nvar; i++) {
    Hcc[i].n = 0;
    Hcc[i].val = NULL;
  }
  int r;
  for (i = 0; i < nmed; i++) {
    for (j = 0; j < medidas[i].nvar; j++) {
      for (r = 0; r < nvar; r++) {
        if (cabs(medidas[i].reguaH[j] - regua[r]) < EPS) {
          // Atualiza a Matriz H
          if (medidas[i].H[j] != 0) {
            sparseAdd(&Hs, i, r, &medidas[i].H[j]);
            sparseAdd(&Ht, r, i, &medidas[i].H[j]);

            sparseAdd(&Hachtel, i, r + nmed, &medidas[i].H[j]);
            sparseAdd(&Hachtel, r + nmed, i, &medidas[i].H[j]);

            Hcc[r].n++;
            sparseAdd(&Hcc[r].val, i, r, &medidas[i].H[j]);
            Hcct[i].n++;
            sparseAdd(&Hcct[r].val, r, i, &medidas[i].H[j]);
          }
          break;
        }
      }
    }
  }

  //    //Para o cálculo da diagonal da matriz ganho
  SPARSE *tmp;
  double diagG[nvar], valor;
  for (i = 0; i < nvar; i++)
    diagG[i] = 0;
  for (j = 0; j < nvar; j++) {
    tmp = Hcc[j].val;
    while (tmp != NULL) {
      diagG[j] += tmp->ij[0] * 1 / W[tmp->j][tmp->j] * tmp->ij[0];
      tmp = tmp->prox;
    }
  }
  CCSPARSE L[nvar + nmed];
  double teste[1];
  teste[0] = 1;
  for (i = 0; i < nvar + nmed; i++) {
    L[i].n = 0;
    L[i].val = NULL;
    sparseAdd(&L[i].val, i, i, &teste[0]);
  }
  for (i = 0; i < nmed; i++) {
    valor = pow(W[i][i], 0.5);
    L[i].val->x = valor;
  }
  for (i = 0; i < nvar; i++) {
    valor = -pow(diagG[i], 0.5);
    L[i].val->x = valor;
  }

  // Montagem da matriz Lower da fatoração Cholesky
  //     for (j=0;j<nmed;j++){
  //         L[j].n++;
  //         valor = 1/pow(W[j][j],0.5);
  //         sparseAdd(&L[j].val, j, j, &valor);
  //
  //         tmp=Hcct[j].val;
  //         while(tmp!=NULL){
  //             L[j].n++;
  //             valor = 1/pow(W[j][j],0.5)*tmp->x;
  //             sparseAdd(&L[j].val, tmp->i + nmed, j, &valor);
  //             tmp = tmp->prox;
  //         }
  //     }
  //     for (j=0;j<nvar;j++){
  //         L[j+nmed].n++;
  //         valor = pow(diagG[j],0.5);
  //         sparseAdd(&L[j+nmed].val, j+nmed, j+nmed, &valor);
  //     }

  //-----> LOOP DO ALGORITIMO WLS
  it = 0;
  conv = 0;

  while ((it < 100)) {
    clock_t t0 = clock();

    // free(Dx);

    for (i = 0; i < nmed; i++)
      Dz[i] = z[i] - *h[i];

    //=========================================================================
    // Algoritimo convergiu
    if (conv == 1) {
      // relatorio_conv(1,tol,tIni);
      printf("\n\nCONVERGENCIA em %d iteracoes\n\n", it);
      fprintf(arquivo, "\n\nCONVERGENCIA em %d iteracoes\n\n", it);
      fprintf(arquivo, "\nTensoes Nodais: Fase-Terra\n");
      for (i = 0; i < numeroBarras; i++) {
        // Retangulares
        // printf("%d\tVa: %.5lf + j%.5lf\tVb: %.5lf + j%.5lf\tVc: %.5lf +
        // j%.5lf\n",grafo[i].barra->ID,__real__ grafo[i].V[0],__imag__
        // grafo[i].V[0],__real__ grafo[i].V[1],__imag__ grafo[i].V[1],__real__
        // grafo[i].V[2],__imag__ grafo[i].V[2]); Polares
        fprintf(arquivo, "%d\t%.5lf\t%.5lf\t%.5lf\t %.3lf\t%.3lf\t%.3lf\n",
                grafo[i].barra->ID, cabs(grafo[i].V[0]), cabs(grafo[i].V[1]),
                cabs(grafo[i].V[2]), carg(grafo[i].V[0]) * 180 / PI,
                carg(grafo[i].V[1]) * 180 / PI, carg(grafo[i].V[2]) * 180 / PI);
      }
      /*fprintf(arquivo,"\nTensoes Nodais Convergido (V):\n");
      for(i=0; i<numeroBarras; i++){
          //Polares
          fprintf(arquivo,"%d\tVa: %.5lf | %.3lf \tVb: %.5lf | %.3lf\tVc: %.5lf
      |
      %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[0])*grafo[i].Vbase,carg(grafo[i].V[0])*180/PI,cabs(grafo[i].V[1])*grafo[i].Vbase,carg(grafo[i].V[1])*180/PI,cabs(grafo[i].V[2])*grafo[i].Vbase,carg(grafo[i].V[2])*180/PI);
      }*/
      fprintf(arquivo, "\n\nVetor z \t\tVetor h(x) \t\t\tVetor Dz\n");
      for (i = 0; i < nmed; i++) {
        fprintf(arquivo, "%.7lf\t%.7lf\t%.7f\n", z[i], *h[i], z[i] - *h[i]);
      }
      nvar = nvar + 3;
      fclose(arquivo);
      return conv;
    }
    //=========================================================================

    //************************************************************************
    // SOLUCAO DO SISTEMA LINEAR: Fatores Trinagulares
    // G(x)*Dx = b
    //************************************************************************
    cat_vert_vet(Dz, nmed, zeros_vet, nvar, b_hach);
    clock_t tHouse = clock();
    // Dx = solve_Householder(Hachtel,nvar + nmed,nvar + nmed,b_hach);
    for (i = 0; i < nmed + nvar; i++)
      Dx[i] = 0;
    sparseCG_preCholesky(L, Hachtel, b_hach, Dx, nvar + nmed);

    nGx = norma_inf(b_hach, nvar);
    nFx = norma_inf(Dx, nvar);
    clock_t t1 = clock();
    double tempo_it = (double)(t1 - t0) / CLOCKS_PER_SEC;

    double tempoHouseholder = (double)(t1 - tHouse) / CLOCKS_PER_SEC;
    printf("\nSolve CG: %lf", tempoHouseholder);
    double tempoTrataMat = (double)(tHouse - t0) / CLOCKS_PER_SEC;
    printf("\nOperacoes Matrizes: %lf", tempoTrataMat);
    printf("\nIteracao: %d \t |Dx| = %f \t |G(x)| = %f \n", it, nFx, nGx);

    // Atualiza os mutilicadores de lagrange mu
    for (i = 0; i < nmed; i++) {
      mu[i] = mu[i] + Dx[i];
    }
    // Atualiza o vetor x
    for (i = 0; i < nvar; i++) {
      ponto[i] = ponto[i] + Dx[i + nmed];
    }

    //************************************************************************
    // TESTE DE CONVERGÊNCIA
    //
    //************************************************************************
    for (i = 0; i < nvar; i++) {
      if (cabs(Dx[i]) >= tol) {
        conv = 0;
        break;
      } else
        conv = 1;
    }
    if (nGx >= 0.00001)
      conv = 0;

    //************************************************************************
    // ATUALIZAÇÃO DO VETOR X
    //
    //************************************************************************
    // atualiza_x(ponto);
    clock_t tM1 = clock();
    atualiza_estado(grafo, ponto, regua, nvar);
    atualiza_h(grafo, numeroBarras, nmed, medidas);
    atualiza_H(grafo, numeroBarras, ramos, medidas, nmed);
    clock_t tM2 = clock();

    double tempoMontaH = (double)(tM2 - tM1) / CLOCKS_PER_SEC;
    printf("\nMonta H: %lf", tempoMontaH);

    it++;
    printf(".");
  }
  //******************************************************************
  // FIM DO LOOP DO WLS
  //
  //******************************************************************
  // atualiza_x(ponto);
  atualiza_estado(grafo, ponto, regua, nvar);
  atualiza_h(grafo, numeroBarras, nmed, medidas);

  // Para o fluxo de carga manter a slack constante
  nvar = nvar + 3;
  // relatorio_conv(0,tol,tIni);

  fprintf(arquivo, "\n\nVetor z \t\tVetor h(x) \t\t\tVetor Dz\n");
  for (i = 0; i < nmed; i++) {
    fprintf(arquivo, "%.7lf\t%.7lf\t%.7f\n", z[i], *h[i], z[i] - *h[i]);
  }

  fclose(arquivo);
  free(Dz);
  free(Dx);
  free(b);
  free(regua_comp);
  free(mu);
  free(zeros_vet);

  return (conv);
}

// Matriz de Hachtel Aumentada - sem virtuais
int otimiza_Gauss_Newton_sparsePCGLS(double *z, double **h, double ***H,
                                     double **W, GRAFO *grafo,
                                     long int numeroBarras, DRAM *ramos,
                                     DMED *medidas, long int nvar,
                                     long int nmed, double *regua,
                                     double *ponto, double tol, long int ref1,
                                     long int ref2) {
  long int it;
  long int NAV = 0, ref;
  optionHachtel = 0;

  double *Dz;
  double *b, *c2;
  double *Dx = NULL;
  double **G;
  double *zeros_vet, *regua_comp;
  int i, j, conv;
  double nGx, nFx;

  SPARSE *Ht = NULL, *Hs = NULL, *Tmp = NULL;

  FILE *arquivo;
  arquivo = fopen("iteracoes.txt", "w");

  FILE *arqout;
  arqout = fopen("matrizHW.txt", "w");

  //----------------------------------------------------------------------------
  //
  // ALOCAÇÃO
  //
  //----------------------------------------------------------------------------
  // free(H_rf);free(Ht);free(Mtmp);free(G);free(Dz);free(b);
  Dz = aloca_vetor(nmed);
  Dx = aloca_vetor(nvar);
  b = aloca_vetor(nvar);
  regua_comp = aloca_vetor(nvar + 3);
  for (i = 0; i < nvar + 3; i++)
    regua_comp[i] = regua[i];

  atualiza_Rede(grafo, numeroBarras);
  // Atualiza modelo das medidas
  atualiza_Modelo(grafo, numeroBarras, nmed, medidas);
  atualiza_H(grafo, numeroBarras, ramos, medidas, nmed);

  //--------------------------------------------------------------------------
  //
  // Precondicionador QR
  CCSPARSE *L;
  L = (CCSPARSE *)malloc(nvar * sizeof(CCSPARSE));

  for (i = 0; i < nvar; i++) {
    L[i].n = 0;
    L[i].val = NULL;
    L[i].diag = 1;
  }
  //    double **R, **H_rf;
  //
  //    H_rf = aloca_matriz(nmed,nvar);
  //    R = aloca_matriz(nmed,nvar);
  //
  //    //Cálculo do precondicionador
  //    for (i=0;i<nmed;i++){
  //        for(j=0;j<medidas[i].nvar;j++){
  //            medidas[i].H[j] = medidas[i].H[j]/medidas[i].sigma;
  //        }
  //    }
  //
  //    it =1;
  //    tira_refs(H,nvar+3,nmed,ref1,ref2,H_rf,regua,c2,it);
  //
  //    fprintf(arqout,"Matriz W1/2*H\n");
  //    for(i=0;i<nmed;i++){
  //        for(j=0;j<nvar;j++) {
  //            fprintf(arqout,"%.15f\t",H_rf[i][j]);
  //        }
  //        fprintf(arqout,"\n");
  //    }
  //
  //    //Fatoração QR para obter a matriz R
  //    QRfactorization(H_rf, nmed, nvar, R);
  //
  //    fprintf(arqout,"Precondicionador\n");
  //    for(i=0;i<nmed;i++){
  //        for(j=0;j<nvar;j++) {
  //            fprintf(arqout,"%.15f\t",R[i][j]);
  //        }
  //        fprintf(arqout,"\n");
  //    }
  //
  //    //Precondicionador QR esparso
  //    double teste[1];
  //    teste[0] = 1;
  //    for(i=0;i<nvar;i++){
  //        L[i].n = 0;
  //        L[i].val = NULL;
  //        L[i].diag = 1;
  ////        teste[0] = R[i][i];
  ////        sparseAdd(&L[i].val, j, i, &teste[0]);
  //        for (j=i;j<nvar;j++){
  //            if (cabs(R[i][j])>0.000000001){
  //                L[i].n++;
  //                teste[0] = R[i][j];
  //                sparseAdd(&L[i].val, j, i, &teste[0]);
  //            }
  //        }
  //    }

  // Opção 2 leitura de arquivo externo
  FILE *precon;
  precon = fopen("precon.txt", "r+");
  int *permut;
  int aux1_precon, aux_i, aux_j;

  // Leitura do vetor de ordenação das variáveis de estado
  fscanf(precon, "%d", &aux1_precon);
  permut = (int *)malloc(aux1_precon * sizeof(int));
  for (i = 0; i < aux1_precon; i++)
    fscanf(precon, "%d", &permut[i]);

  // Aplica permutação na regua

  fscanf(precon, "%d", &aux1_precon);
  //    aux1_precon = 7802;
  // Leitura da matriz R
  double teste[1];
  teste[0] = 1;
  for (i = 0; i < aux1_precon; i++) {
    fscanf(precon, "%d\t%d\t%lf\n", &aux_i, &aux_j, &teste[0]);
    L[aux_i].n++;
    sparseAdd(&L[aux_i - 1].val, aux_j - 1, aux_i - 1, &teste[0]);
  }

  fclose(precon);

  //--------------------------------------------------------------------------

  // Tira a coluna de angulo da referencia na matriz H(x)
  if (NAV == 0) {
    // tira_refs(H,nvar+3,nmed,ref1,ref2,H_rf,regua,ponto,it); // IEEE4

    for (j = 0; j < nvar + 3; j++) {
      if (j < ref1) {
        regua[j] = regua[j];
        ponto[j] = ponto[j];
      } else if (j > ref2) {
        regua[(j - (ref2 - ref1 + 1))] = regua[j];
        ponto[(j - (ref2 - ref1 + 1))] = ponto[j];
      }
    }
  }

  int r;
  for (i = 0; i < nmed; i++) {
    for (j = 0; j < medidas[i].nvar; j++) {
      for (r = 0; r < nvar; r++) {
        if (cabs(medidas[i].reguaH[j] - regua[r]) < EPS) {
          // Atualiza a Matriz H
          if (medidas[i].H[j] != 0) {
            sparseAdd(&Hs, i, r, &medidas[i].H[j]);
          }
          break;
        }
      }
    }
  }

  // Precondicionador QR
  //     CCSPARSE L[nvar];
  //     double teste[1];
  //     teste[0] = 1;
  //     for(i=0;i<nvar;i++){
  //         L[i].n = 0;
  //         L[i].val = NULL;
  //         L[i].diag = 1;
  //         sparseAdd(&L[i].val, i, i, &teste[0]);
  //     }

  //-----> LOOP DO ALGORITIMO WLS
  it = 0;
  conv = 0;
  clock_t tIni = clock();

  while ((it < 50)) {
    clock_t t0 = clock();
    for (i = 0; i < nmed; i++) {
      Dz[i] = (z[i] - *h[i]) / medidas[i].sigma;
    }

    fprintf(arqout, "vetor\n");
    for (i = 0; i < nmed; i++) {
      fprintf(arqout, "%.15f\n", Dz[i]);
    }
    fprintf(arqout, "\n");
    fclose(arqout);
    //=========================================================================
    // Algoritimo convergiu
    if (conv == 1) {
      clock_t tFim = clock();
      double tempoIt = (double)(tFim - tIni) / CLOCKS_PER_SEC;
      printf("\n\n Convergência em %d iteracoes e tempo: %.4lf", it, tempoIt);

      fprintf(arquivo, "\n\nCONVERGENCIA em %d iteracoes\n\n", it);
      fprintf(arquivo, "\nTensoes Nodais: Fase-Terra\n");
      for (i = 0; i < numeroBarras; i++) {
        // Retangulares
        // printf("%d\tVa: %.5lf + j%.5lf\tVb: %.5lf + j%.5lf\tVc: %.5lf +
        // j%.5lf\n",grafo[i].barra->ID,__real__ grafo[i].V[0],__imag__
        // grafo[i].V[0],__real__ grafo[i].V[1],__imag__ grafo[i].V[1],__real__
        // grafo[i].V[2],__imag__ grafo[i].V[2]); Polares
        fprintf(arquivo, "%d\t%.5lf\t%.5lf\t%.5lf\t %.3lf\t%.3lf\t%.3lf\n",
                grafo[i].barra->ID, cabs(grafo[i].V[0]), cabs(grafo[i].V[1]),
                cabs(grafo[i].V[2]), carg(grafo[i].V[0]) * 180 / PI,
                carg(grafo[i].V[1]) * 180 / PI, carg(grafo[i].V[2]) * 180 / PI);
      }
      /*fprintf(arquivo,"\nTensoes Nodais Convergido (V):\n");
      for(i=0; i<numeroBarras; i++){
          //Polares
          fprintf(arquivo,"%d\tVa: %.5lf | %.3lf \tVb: %.5lf | %.3lf\tVc: %.5lf
      |
      %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[0])*grafo[i].Vbase,carg(grafo[i].V[0])*180/PI,cabs(grafo[i].V[1])*grafo[i].Vbase,carg(grafo[i].V[1])*180/PI,cabs(grafo[i].V[2])*grafo[i].Vbase,carg(grafo[i].V[2])*180/PI);
      }*/
      fprintf(arquivo, "\n\nVetor z \t\tVetor h(x) \t\t\tVetor Dz\n");
      for (i = 0; i < nmed; i++) {
        fprintf(arquivo, "%.7lf\t%.7lf\t%.7f\n", z[i], *h[i], z[i] - *h[i]);
      }

      fclose(arquivo);
      return conv;
    }
    //=========================================================================

    //************************************************************************
    // SOLUCAO DO SISTEMA LINEAR: Fatores Trinagulares
    // G(x)*Dx = b
    //************************************************************************
    for (i = 0; i < nvar; i++)
      Dx[i] = 0;
    clock_t tHouse = clock();
    //        sparseCG_preCholesky(L,Hachtel,b_hach,Dx,nvar+nmed+nvir);
    sparsePCGLS(L, Hs, Dz, Dx, nmed, nvar);

    for (i = 0; i < nvar; i++)
      b[i] = 0;
    sparseMultMatTransVet(Hs, Dz, b); // Grad = Ht*W*Dz

    nGx = norma_inf(b, nvar);
    nFx = norma_inf(Dx, nvar);
    clock_t t1 = clock();
    double tempo_it = (double)(t1 - t0) / CLOCKS_PER_SEC;

    double tempoHouseholder = (double)(t1 - tHouse) / CLOCKS_PER_SEC;
    printf("\nSolve PCGLS: %lf", tempoHouseholder);
    double tempoTrataMat = (double)(tHouse - t0) / CLOCKS_PER_SEC;
    printf("\nOperacoes Matrizes: %lf", tempoTrataMat);
    printf("\nIteracao: %d \t |Dx| = %.17lf \t |G(x)| = %.17lf \n", it, nFx,
           nGx);

    double passo = 1.0;
    //        double Armijo_c1 = 0.01;
    //        double alpha0 = 2;
    ////        double passo = NewtonMod_calcPasso(alpha0, b, Dx, ponto, nvar,
    ///Armijo_c1, grafo,numeroBarras, medidas, nmed, ponto, nvar, z,h, W,
    ///regua);
    //        if (it>=4) passo = 0.8;
    //        if (it>=10) passo = 0.0010;

    // Atualiza o vetor x
    for (i = 0; i < nvar; i++) {
      ponto[i] = ponto[i] + passo * Dx[i];
    }

    //************************************************************************
    // TESTE DE CONVERGÊNCIA
    //
    //************************************************************************
    for (i = 0; i < nvar; i++) {
      if (cabs(Dx[i]) >= tol) {
        conv = 0;
        break;
      } else
        conv = 1;
    }
    //        if (nGx >= 0.01)
    //            conv = 0;

    //************************************************************************
    // ATUALIZAÇÃO DO VETOR X
    //
    //************************************************************************
    // atualiza_x(ponto);
    clock_t tM1 = clock();
    atualiza_estado(grafo, ponto, regua, nvar);
    atualiza_Rede(grafo, numeroBarras);
    // Atualiza modelo das medidas
    atualiza_Modelo(grafo, numeroBarras, nmed, medidas);
    atualiza_H(grafo, numeroBarras, ramos, medidas, nmed);
    for (i = 0; i < nmed; i++) {
      for (j = 0; j < medidas[i].nvar; j++) {
        medidas[i].H[j] = medidas[i].H[j] / medidas[i].sigma;
      }
    }

    clock_t tM2 = clock();

    double tempoMontaH = (double)(tM2 - tM1) / CLOCKS_PER_SEC;
    printf("\nMonta H: %lf", tempoMontaH);

    it++;
    printf(".");
  }
  //******************************************************************
  // FIM DO LOOP DO WLS
  //
  //******************************************************************
  // atualiza_x(ponto);
  atualiza_estado(grafo, ponto, regua, nvar);
  atualiza_Rede(grafo, numeroBarras);
  atualiza_Modelo(grafo, numeroBarras, nmed, medidas);

  fprintf(arquivo, "\n\nVetor z \t\tVetor h(x) \t\t\tVetor Dz\n");
  for (i = 0; i < nmed; i++) {
    fprintf(arquivo, "%.7lf\t%.7lf\t%.7f\n", z[i], *h[i], z[i] - *h[i]);
  }

  fclose(arquivo);
  //    free(Dz);free(Dx);free(b);free(regua_comp);

  return (conv);
}

// Matriz de Hachtel Aumentada - sem virtuais
int otimiza_Gauss_Newton_sparseHachtel_Virtuais(
    double *z, double **h, double **c, double **W, GRAFO *grafo,
    long int numeroBarras, DRAM *ramos, DMED *medidas, DMED *virtuais,
    long int nvar, long int nmed, long int nvir, double *regua, double *ponto,
    double tol, long int ref1, long int ref2) {
  long int it;
  long int NAV = 0, ref;
  optionHachtel = 1;

  double *Dz;
  double *b, *c2;
  double *Dx = NULL;
  double **G;
  double *mu, *lamb, *b_hach, *zeros_vet, *regua_comp;
  int i, j, conv;
  double nGx, nFx;

  SPARSE *Ht = NULL, *Hs = NULL, *Hachtel = NULL, *Cov = NULL, *Tmp = NULL;

  FILE *arquivo;
  arquivo = fopen("iteracoes.txt", "w");

  //----------------------------------------------------------------------------
  //
  // ALOCAÇÃO
  //
  //----------------------------------------------------------------------------
  // free(H_rf);free(Ht);free(Mtmp);free(G);free(Dz);free(b);
  Dz = aloca_vetor(nmed);
  Dx = aloca_vetor(nmed + nvar + nvir);
  // b = aloca_vetor(nvar);
  regua_comp = aloca_vetor(nvar + 3);
  for (i = 0; i < nvar + 3; i++)
    regua_comp[i] = regua[i];

  b_hach = aloca_vetor(nmed + nvar + nvir);
  mu = aloca_vetor(nmed);
  lamb = aloca_vetor(nvir);
  c2 = aloca_vetor(nvir);
  zeros_vet = aloca_vetor(nvar);

  atualiza_Rede(grafo, numeroBarras);
  // Atualiza modelo das medidas
  atualiza_Modelo(grafo, numeroBarras, nmed, medidas);
  atualiza_H(grafo, numeroBarras, ramos, medidas, nmed);
  // Atualiza modelo das virtuais
  atualiza_Modelo(grafo, numeroBarras, nvir, virtuais);
  atualiza_H(grafo, numeroBarras, ramos, virtuais, nvir);

  // Tira a coluna de angulo da referencia na matriz H(x)
  if (NAV == 0) {
    // tira_refs(H,nvar+3,nmed,ref1,ref2,H_rf,regua,ponto,it); // IEEE4
    for (j = 0; j < nvar + 3; j++) {
      if (j < ref1) {
        regua[j] = regua[j];
        ponto[j] = ponto[j];
      } else if (j > ref2) {
        regua[(j - (ref2 - ref1 + 1))] = regua[j];
        ponto[(j - (ref2 - ref1 + 1))] = ponto[j];
      }
    }
  }

  //----------------------------------------------------------------------------
  // Insere a Matriz de Covariancia na Matriz de hachtel
  for (i = 0; i < nmed; i++) {
    W[i][i] = 1 / W[i][i]; // Na Hachtel usa a matriz de covariância
    sparseAdd(&Hachtel, i, i, &W[i][i]);
  }

  // Monta matriz H e Ht esparsa no formato compressed column
  // Matriz H esparsa
  CCSPARSE Hcc[nvar], Hcct[nmed];
  for (i = 0; i < nmed; i++) {
    Hcct[i].n = 0;
    Hcct[i].val = NULL;
  }
  for (i = 0; i < nvar; i++) {
    Hcc[i].n = 0;
    Hcc[i].val = NULL;
  }
  int r;
  for (i = 0; i < nmed; i++) {
    for (j = 0; j < medidas[i].nvar; j++) {
      for (r = 0; r < nvar; r++) {
        if (cabs(medidas[i].reguaH[j] - regua[r]) < EPS) {
          // Atualiza a Matriz H
          if (medidas[i].H[j] != 0) {
            sparseAdd(&Hs, i, r, &medidas[i].H[j]);
            sparseAdd(&Ht, r, i, &medidas[i].H[j]);

            sparseAdd(&Hachtel, i, r + nmed, &medidas[i].H[j]);
            sparseAdd(&Hachtel, r + nmed, i, &medidas[i].H[j]);

            //                        Hcc[r].n++;
            //                        sparseAdd(&Hcc[r].val, i, r,
            //                        &medidas[i].H[j]); Hcct[i].n++;
            //                        sparseAdd(&Hcct[r].val, r, i,
            //                        &medidas[i].H[j]);
          }
          break;
        }
      }
    }
  }

  // Inseri as medidas virtuais na matriz de Hachtel
  for (i = 0; i < nvir; i++) {
    for (j = 0; j < virtuais[i].nvar; j++) {
      for (r = 0; r < nvar; r++) {
        if (cabs(virtuais[i].reguaH[j] - regua[r]) < EPS) {
          // Atualiza a Matriz H
          if (virtuais[i].H[j] != 0) {
            sparseAdd(&Hachtel, i + nmed + nvar, r + nmed, &virtuais[i].H[j]);
            sparseAdd(&Hachtel, r + nmed, i + nmed + nvar, &virtuais[i].H[j]);
          }
          break;
        }
      }
    }
  }

  //    //Para o cálculo da diagonal da matriz ganho
  SPARSE *tmp;
  double diagG[nvar], valor;
  //    for(i=0;i<nvar;i++) diagG[i] = 0;
  //    for(j=0;j<nvar;j++){
  //        tmp = Hcc[j].val;
  //        while (tmp != NULL){
  //            diagG[j] += tmp->ij[0] * 1/W[tmp->j][tmp->j] * tmp->ij[0];
  //            tmp = tmp->prox;
  //        }
  //    }
  CCSPARSE L[nvar + nmed + nvir];
  double teste[1];
  teste[0] = 1;
  for (i = 0; i < nvar + nmed + nvir; i++) {
    L[i].n = 0;
    L[i].val = NULL;
    L[i].diag = 1;
    sparseAdd(&L[i].val, i, i, &teste[0]);
  }

  //    SPARSE *atual = &Hachtel;
  //    while(atual != NULL){
  //        i = atual->i;
  //        j = atual->j;
  //
  //        if (i == j){
  //            valor = atual->x;
  //            L[i].val->x = pow(valor,0.5);
  //            L[i].diag = pow(valor,0.5);
  //        }
  //        else if (j > i){
  //            sparseAdd(&L[i].val, i, j, &teste[0]);
  //            valor = atual->x;
  //            L[i].val->x = ;
  //        }
  //        printf("\n[%d , %d]: %.5lf",atual->i,atual->j,atual->x);
  //        atual = atual->prox;
  //    }

  //    for(i=0;i<nmed;i++) {
  //        valor = pow(W[i][i],0.5);
  //        L[i].val->x = valor;
  //    }
  //    for(i=0;i<nvar;i++) {
  //        valor = -pow(diagG[i],0.5);
  //        L[i].val->x = valor;
  //    }

  // Montagem da matriz Lower da fatoração Cholesky
  //     for (j=0;j<nmed;j++){
  //         L[j].n++;
  //         valor = 1/pow(W[j][j],0.5);
  //         sparseAdd(&L[j].val, j, j, &valor);
  //
  //         tmp=Hcct[j].val;
  //         while(tmp!=NULL){
  //             L[j].n++;
  //             valor = 1/pow(W[j][j],0.5)*tmp->x;
  //             sparseAdd(&L[j].val, tmp->i + nmed, j, &valor);
  //             tmp = tmp->prox;
  //         }
  //     }
  //     for (j=0;j<nvar;j++){
  //         L[j+nmed].n++;
  //         valor = pow(diagG[j],0.5);
  //         sparseAdd(&L[j+nmed].val, j+nmed, j+nmed, &valor);
  //     }

  //-----> LOOP DO ALGORITIMO WLS
  it = 0;
  conv = 0;

  while ((it < 50)) {
    clock_t t0 = clock();

    // free(Dx);

    for (i = 0; i < nmed; i++)
      Dz[i] = z[i] - *h[i];
    for (i = 0; i < nvir; i++)
      c2[i] = -*c[i];

    //=========================================================================
    // Algoritimo convergiu
    if (conv == 1) {
      // relatorio_conv(1,tol,tIni);
      printf("\n\nCONVERGENCIA em %d iteracoes\n\n", it);
      fprintf(arquivo, "\n\nCONVERGENCIA em %d iteracoes\n\n", it);
      fprintf(arquivo, "\nTensoes Nodais: Fase-Terra\n");
      for (i = 0; i < numeroBarras; i++) {
        // Retangulares
        // printf("%d\tVa: %.5lf + j%.5lf\tVb: %.5lf + j%.5lf\tVc: %.5lf +
        // j%.5lf\n",grafo[i].barra->ID,__real__ grafo[i].V[0],__imag__
        // grafo[i].V[0],__real__ grafo[i].V[1],__imag__ grafo[i].V[1],__real__
        // grafo[i].V[2],__imag__ grafo[i].V[2]); Polares
        fprintf(arquivo, "%d\t%.5lf\t%.5lf\t%.5lf\t %.3lf\t%.3lf\t%.3lf\n",
                grafo[i].barra->ID, cabs(grafo[i].V[0]), cabs(grafo[i].V[1]),
                cabs(grafo[i].V[2]), carg(grafo[i].V[0]) * 180 / PI,
                carg(grafo[i].V[1]) * 180 / PI, carg(grafo[i].V[2]) * 180 / PI);
      }
      /*fprintf(arquivo,"\nTensoes Nodais Convergido (V):\n");
      for(i=0; i<numeroBarras; i++){
          //Polares
          fprintf(arquivo,"%d\tVa: %.5lf | %.3lf \tVb: %.5lf | %.3lf\tVc: %.5lf
      |
      %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[0])*grafo[i].Vbase,carg(grafo[i].V[0])*180/PI,cabs(grafo[i].V[1])*grafo[i].Vbase,carg(grafo[i].V[1])*180/PI,cabs(grafo[i].V[2])*grafo[i].Vbase,carg(grafo[i].V[2])*180/PI);
      }*/
      fprintf(arquivo, "\n\nVetor z \t\tVetor h(x) \t\t\tVetor Dz\n");
      for (i = 0; i < nmed; i++) {
        fprintf(arquivo, "%.7lf\t%.7lf\t%.7f\n", z[i], *h[i], z[i] - *h[i]);
      }
      fprintf(arquivo, "\n\nVetor c\n");
      for (i = 0; i < nvir; i++) {
        fprintf(arquivo, "%.7lf\t\n", *c[i]);
      }
      //            nvar = nvar +3;
      fclose(arquivo);
      return conv;
    }
    //=========================================================================

    fprintf(arquivo, "\n\nVetor z \t\tVetor h(x) \t\t\tVetor Dz\n");
    for (i = 0; i < nmed; i++) {
      fprintf(arquivo, "%.7lf\t%.7lf\t%.7f\n", z[i], *h[i], z[i] - *h[i]);
    }
    fprintf(arquivo, "\n\nVetor c\n");
    for (i = 0; i < nvir; i++) {
      fprintf(arquivo, "%.7lf\t\n", c2[i]);
    }

    //************************************************************************
    // SOLUCAO DO SISTEMA LINEAR: Fatores Trinagulares
    // G(x)*Dx = b
    //************************************************************************
    cat_vert_vet(Dz, nmed, zeros_vet, nvar, b_hach);
    cat_vert_vet(b_hach, nmed + nvar, c2, nvir, b_hach);
    clock_t tHouse = clock();
    // Dx = solve_Householder(Hachtel,nvar + nmed,nvar + nmed,b_hach);
    for (i = 0; i < nmed + nvar + nvir; i++)
      Dx[i] = 0;
    sparseCG_preCholesky(L, Hachtel, b_hach, Dx, nvar + nmed + nvir);

    nGx = norma_inf(b_hach, nvar + nmed + nvir);
    nFx = norma_inf(Dx, nvar + nmed + nvir);
    clock_t t1 = clock();
    double tempo_it = (double)(t1 - t0) / CLOCKS_PER_SEC;

    double tempoHouseholder = (double)(t1 - tHouse) / CLOCKS_PER_SEC;
    printf("\nSolve CG Hachtel: %lf", tempoHouseholder);
    double tempoTrataMat = (double)(tHouse - t0) / CLOCKS_PER_SEC;
    printf("\nOperacoes Matrizes: %lf", tempoTrataMat);
    printf("\nIteracao: %d \t |Dx| = %f \t |G(x)| = %f \n", it, nFx, nGx);

    // Atualiza os mutilicadores de lagrange mu
    for (i = 0; i < nmed; i++) {
      mu[i] = mu[i] + Dx[i];
    }
    // Atualiza o vetor x
    for (i = 0; i < nvar; i++) {
      ponto[i] = ponto[i] + Dx[i + nmed];
    }
    // Atualiza os mutilicadores de lagrange lambda
    for (i = 0; i < nvir; i++) {
      lamb[i] = lamb[i] + Dx[i + nmed + nvar];
    }
    //************************************************************************
    // TESTE DE CONVERGÊNCIA
    //
    //************************************************************************
    for (i = 0; i < nvar + nmed + nvir; i++) {
      if (cabs(Dx[i]) >= tol) {
        conv = 0;
        break;
      } else
        conv = 1;
    }
    //        if (nGx >= tol)
    //            conv = 0;

    //************************************************************************
    // ATUALIZAÇÃO DO VETOR X
    //
    //************************************************************************
    // atualiza_x(ponto);
    clock_t tM1 = clock();
    atualiza_estado(grafo, ponto, regua, nvar);

    atualiza_Rede(grafo, numeroBarras);
    // Atualiza modelo das medidas
    atualiza_Modelo(grafo, numeroBarras, nmed, medidas);
    atualiza_H(grafo, numeroBarras, ramos, medidas, nmed);
    // Atualiza modelo das virtuais
    atualiza_Modelo(grafo, numeroBarras, nvir, virtuais);
    atualiza_H(grafo, numeroBarras, ramos, virtuais, nvir);
    clock_t tM2 = clock();

    double tempoMontaH = (double)(tM2 - tM1) / CLOCKS_PER_SEC;
    printf("\nMonta H: %lf", tempoMontaH);

    it++;
    printf(".");
  }
  //******************************************************************
  // FIM DO LOOP DO WLS
  //
  //******************************************************************
  // atualiza_x(ponto);
  atualiza_estado(grafo, ponto, regua, nvar);
  atualiza_Rede(grafo, numeroBarras);
  atualiza_Modelo(grafo, numeroBarras, nmed, medidas);

  fprintf(arquivo, "\n\nVetor z \t\tVetor h(x) \t\t\tVetor Dz\n");
  for (i = 0; i < nmed; i++) {
    fprintf(arquivo, "%.7lf\t%.7lf\t%.7f\n", z[i], *h[i], z[i] - *h[i]);
  }

  fclose(arquivo);
  //    free(Dz);free(Dx);free(b);free(regua_comp);
  free(mu);
  free(zeros_vet);

  return (conv);
}

int otimizaNEC(double *z, double **h, double ***H, double ***C, GRAFO *grafo,
               long int numeroBarras, DRAM *ramos, DMED *medidas,
               DMED *virtuais, long int nvir, long int nvar, long int nmed,
               double *regua_comp, double *ponto, double tol, long int ref1,
               long int ref2) {
  double **C_rf;
  double **Haum;
  double *Dz;
  double *Dv;
  double *Dx;
  double nGx, nFx;
  double **H_rf;
  double **H_T;
  double *Dz_aux;
  int i, j, k;
  double *b = NULL;
  H_rf = aloca_matriz(nmed, nvar);
  H_T = aloca_matriz(nvar, nmed);
  C_rf = aloca_matriz(nvir, nvar);
  Dz = aloca_vetor(nmed);
  Dz_aux = aloca_vetor(nvar);
  Dv = aloca_vetor(nvir);
  Dx = aloca_vetor(nvar + nvir);
  Haum = aloca_matriz(nvar, nvar);
  b = aloca_vetor(nvar);

  atualiza_Rede(
      grafo,
      numeroBarras); // atualiza a condição da rede conforme o estado atual
  atualiza_Modelo(
      grafo, numeroBarras, nmed,
      medidas); // atualiza modelo de medição conforme a condição atual da rede
  atualiza_Modelo(grafo, numeroBarras, nvir, virtuais);

  int it = 0;
  int conv = 0;

  while (it < 30) {
    if (conv == 1) {
      FILE *res = NULL;
      res = fopen("residuo.txt", "w+");
      int ct = 0;
      for (int i = 0; i < nmed + nvir; i++) {
        if (ct < nmed) {
          float resid = medidas[ct].zmed - medidas[ct].h;
          fprintf(res, "%f\n", resid);
        } else {
          float resid = virtuais[ct].zmed - virtuais[ct].h;
          fprintf(res, "%f\n", resid);
        }
        ct += 1;
      }
      fclose(res);
      return conv;
    }
    atualiza_H(grafo, numeroBarras, ramos, medidas, nmed);
    atualiza_H(grafo, numeroBarras, ramos, virtuais, nvir);
    double max_sigma = 0;
    for (i = 0; i < nmed; i++) {
      for (j = 0; j < medidas[i].nvar; j++) {
        // H.T * H / (sigma^2)
        medidas[i].H[j] = (medidas[i].H[j]) / medidas[i].sigma;
      }
      // H.T*deltaZ/(sigma^2)
      if (medidas[i].sigma > max_sigma) {
        max_sigma = medidas[i].sigma;
      }
      Dz[i] = (medidas[i].zmed - medidas[i].h) / medidas[i].sigma;
    }
    for (i = 0; i < nvir; i++) {
      Dv[i] = virtuais[i].h;
    }
    double soma;
    float somaDz;
    int _i, _j;
    // printf("max sigma: %f\n", (max_sigma));

    float somaVl = 0;
    for (int i = 0; i < nmed; i++) {
      somaVl += medidas[i].sigma;
    }

    double vl = (max_sigma);
    // double vl = sqrt(nmed)*somaVl;
    // double vl = 1;
    // printf("vl: %f\n", vl);
    for (i = 0; i < nmed; i++) {
      for (j = 0; j < nvar; j++) {
        H_rf[i][j] = *H[i][j];
      }
    }
    // printf("vl: %f\n", vl);

    // transposta da matriz w1/2H
    for (i = 0; i < nvar; i++) {
      for (j = 0; j < nmed; j++) {
        H_T[i][j] = H_rf[j][i];
      }
    }
    printf("vl: %f\n", vl);

    ////multiplica H.T por H
    for (int k = 0; k < nvar; k++) {

      somaDz = 0;
      for (i = 0; i < nmed; i++) {
        if (H_rf[i][k] != 0 && Dz[i] != 0) {
          somaDz += H_rf[i][k] * Dz[i];
        }
      }
      Dz_aux[k] = pow(vl, 2) * somaDz;
    }
    //
    cholmod_triplet *T_nec = NULL;

    cholmod_triplet *T_aux = NULL;
    cholmod_sparse *Gaux = NULL;

    cholmod_triplet *C_T = NULL;
    cholmod_triplet *Cn = NULL;
    cholmod_sparse *sC_T = NULL;
    cholmod_sparse *sCn = NULL;

    cholmod_sparse *horzC = NULL;

    cholmod_triplet *T_tr = NULL;

    cholmod_sparse *A_aux = NULL;
    cholmod_sparse *A_T = NULL;

    cholmod_dense *b_nec = NULL;
    cholmod_dense *X_nec = NULL;

    cholmod_sparse *A_nec = NULL;
    cholmod_factor *L = NULL;
    cholmod_common Common, *c;

    c = &Common;
    cholmod_l_start(c);

    T_nec = cholmod_l_allocate_triplet(nvar + nvir, nvar + nvir,
                                       (nvar + nvir) * (nvar + nvir), 0,
                                       CHOLMOD_REAL, c);

    T_aux =
        cholmod_l_allocate_triplet(nmed, nvar, nmed * nvar, 0, CHOLMOD_REAL, c);
    Gaux = cholmod_l_allocate_sparse(nvar, nvar, nvar * nvar, 0, 0, 0,
                                     CHOLMOD_REAL, c);

    C_T =
        cholmod_l_allocate_triplet(nvar, nvir, nvar * nvir, 0, CHOLMOD_REAL, c);
    sC_T = cholmod_l_allocate_sparse(nvar, nvir, nvar * nvir, 0, 0, 0,
                                     CHOLMOD_REAL, c);

    Cn = cholmod_l_allocate_triplet(nvir, nvar + nvir, nvar * nvir, 0,
                                    CHOLMOD_REAL, c);
    sCn = cholmod_l_allocate_sparse(nvir, nvar + nvir, nvar * nvir, 0, 0, 0,
                                    CHOLMOD_REAL, c);

    T_tr =
        cholmod_l_allocate_triplet(nvar, nvar, nvar * nvar, 0, CHOLMOD_REAL, c);
    A_aux = cholmod_l_allocate_sparse(nmed, nvar, nmed * nvar, 0, 0, 0,
                                      CHOLMOD_REAL, c);
    A_T = cholmod_l_allocate_sparse(nvar, nmed, nmed * nvar, 0, 0, 0,
                                    CHOLMOD_REAL, c);

    A_nec = cholmod_l_allocate_sparse(nvar + nvir, nvar + nvir,
                                      (nvar + nvir) * (nvar + nvir), 0, 0, 0,
                                      CHOLMOD_REAL, c);
    b_nec = cholmod_l_allocate_dense(nvar + nvir, 1, (nvar + nvir),
                                     CHOLMOD_REAL, c);
    X_nec = cholmod_l_allocate_dense(nvar + nvir, 1, (nvar + nvir),
                                     CHOLMOD_REAL, c);

    // calculo de transposta via rotinas esparsas
    int index = 0;
    for (int i = 0; i < nmed; i++) {
      for (int r = 0; r < nvar; r++) {
        if ((H_rf[i][r]) != 0) {
          ((long int *)T_aux->i)[index] = i;
          ((long int *)T_aux->j)[index] = r;
          ((double *)T_aux->x)[index] = vl * H_rf[i][r];
          T_aux->nnz += 1;
          index += 1;
        }
      }
    }

    A_aux = cholmod_l_triplet_to_sparse(T_aux, nmed * nvar, c);
    A_T = cholmod_l_transpose(A_aux, 1, c);
    Gaux = cholmod_l_ssmult(A_T, A_aux, CHOLMOD_REAL, true, false, c);

    // T_tr = cholmod_l_sparse_to_triplet(Gaux, c);

    // int nnz_gain = T_tr->nnz;
    // for (int k=0; k<nnz_gain;k++){
    //     long int _i = ((long int *)T_tr->i)[k];
    //     long int _j = ((long int *)T_tr->j)[k];
    //     double _x = ((double *)T_tr->x)[k];
    //
    //    Haum[_i][_j] = _x;
    //
    //}
    //

    // for (int i = 0; i < nvar; i++)
    //{
    //     for (int r = 0; r < nvar; r++)
    //     {
    //         if ((Haum[i][r]) != 0)
    //         {
    //             ((long int *)T_nec->i)[index] = i;
    //             ((long int *)T_nec->j)[index] = r;
    //             ((double *)T_nec->x)[index] = Haum[i][r];
    //             T_nec->nnz += 1;
    //             index += 1;
    //             ((long int *)T_nec->i)[index] = r;
    //             ((long int *)T_nec->j)[index] = i;
    //             ((double *)T_nec->x)[index] = Haum[i][r];
    //             T_nec->nnz += 1;
    //             index += 1;
    //         }
    //     }
    // }
    index = 0;
    for (i = 0; i < nvir; i++) {
      for (int r = 0; r < nvar; r++) {
        if ((*C[i][r]) != 0) {
          ((long int *)Cn->i)[index] = i;
          ((long int *)Cn->j)[index] = r;
          ((double *)Cn->x)[index] = -1 * (*C[i][r]);
          Cn->nnz += 1;
          index += 1;
        }
      }
    }
    index = 0;
    for (i = 0; i < nvar; i++) {
      for (int r = 0; r < nvir; r++) {
        if ((*C[r][i]) != 0) {
          ((long int *)C_T->i)[index] = i;
          ((long int *)C_T->j)[index] = r;
          ((double *)C_T->x)[index] = -1 * (*C[r][i]);
          C_T->nnz += 1;
          index += 1;
        }
      }
    }
    sCn = cholmod_l_triplet_to_sparse(Cn, nvar * nvir, c);
    sC_T = cholmod_l_triplet_to_sparse(C_T, nvar * nvir, c);

    horzC = cholmod_l_horzcat(Gaux, sC_T, 1, c);
    A_nec = cholmod_l_vertcat(horzC, sCn, 1, c);
    T_nec = cholmod_l_sparse_to_triplet(A_nec, c);

    printf("nmed: %d, nvir: %d\n", nmed, nvir);

    BOOL write = false;

    if (write && it == 0) {
      FILE *matNEC;
      matNEC = fopen("matnec_123_alfa0.txt", "w+");
      for (i = 0; i < T_nec->nnz; i++) {
        long int _i, _j;
        double v;
        _i = ((long int *)T_nec->i)[i];
        _j = ((long int *)T_nec->j)[i];
        v = ((double *)T_nec->x)[i];
        fprintf(matNEC, "%d,%d,%f\n", _i, _j, v);
      }
    }

    for (int i = 0; i < nvar; i++) {
      ((double *)b_nec->x)[i] = Dz_aux[i];
    }
    for (int i = 0; i < nvir; i++) {
      ((double *)b_nec->x)[i + nvar] = Dv[i];
    }

    // A_nec = cholmod_l_triplet_to_sparse(T_nec, (nvar + nvir) * (nvar + nvir),
    // c); L = cholmod_l_analyze(A_nec, c); cholmod_l_factorize(A_nec, L, c);
    // X_nec = cholmod_l_solve(CHOLMOD_A, L, b_nec, c);
    // printf("nmed: %d, nvir: %d, nvar: %d\n", nmed, nvir, nvar);
    // printf("(A) nrow: %ld, ncol: %ld.\n", (long int)A_nec->nrow, (long
    // int)A_nec->ncol); printf("(B) nrow: %ld\n", (long int)b_nec->nrow);
    X_nec = SuiteSparseQR_C_backslash(SPQR_ORDERING_AMD, SPQR_DEFAULT_TOL,
                                      A_nec, b_nec, c);

    Dx = (double *)X_nec->x;

    for (i = 0; i < nvar; i++) {
      ponto[i] = ponto[i] + Dx[i];

      for (j = 0; j < nmed; j++) {
        b[i] = b[i] + (*H[j][i]) * Dz_aux[j];
      }
    }
    for (i = 0; i < nvar; i++) {
      if (fabs(Dx[i]) >= tol) {
        conv = 0;
        break;
      } else
        conv = 1;
    }

    atualiza_estado(
        grafo, ponto, regua_comp,
        nvar); // atualiza o estado atual do grafo conforme o vetor x calculado
    atualiza_Rede(grafo, numeroBarras);
    atualiza_Modelo(grafo, numeroBarras, nmed, medidas);
    atualiza_Modelo(grafo, numeroBarras, nvir, virtuais);

    // atualiza_h(grafo, numeroBarras, nmed, medidas);
    // atualiza_h(grafo, numeroBarras, nvir, virtuais);

    nFx = norma_inf(Dx, nvar);
    nGx = norma_inf(b, nvar);
    printf("\n\nIteracao:  %d \t|Dx|_inf =  %.17lf \t |Grad|_inf =  %.17lf \n",
           it, nFx, nGx);

    it++;
    printf(".");
  }
  atualiza_estado(
      grafo, ponto, regua_comp,
      nvar); // atualiza o estado atual do grafo conforme o vetor x calculado
  atualiza_Rede(grafo, numeroBarras);
  atualiza_Modelo(grafo, numeroBarras, nmed, medidas);
  atualiza_Modelo(grafo, numeroBarras, nvir, virtuais);
  // atualiza_h(grafo, numeroBarras, nmed, medidas);
  // atualiza_h(grafo, numeroBarras, nvir, virtuais);

  return conv;
}
int otimizaHatchel(double *z, double **h, double ***H, double ***C,
                   GRAFO *grafo, long int numeroBarras, DRAM *ramos,
                   DMED *medidas, DMED *virtuais, long int nvir, long int nvar,
                   long int nmed, double *regua_comp, double *ponto, double tol,
                   long int ref1, long int ref2) {
  double **C_rf;
  double **Haum;
  double *Dz;
  double *Dv;
  double *Dx;
  double nGx, nFx;
  double **H_rf;
  double **H_T;
  double *Dz_aux;
  int i, j, k;
  double *b = NULL;
  H_rf = aloca_matriz(nmed, nvar);
  H_T = aloca_matriz(nvar, nmed);
  C_rf = aloca_matriz(nvir, nvar);
  Dz = aloca_vetor(nmed);
  Dz_aux = aloca_vetor(nvar);
  Dv = aloca_vetor(nvir);
  Dx = aloca_vetor(nvar + nvir);
  Haum = aloca_matriz(nvar, nvar);
  b = aloca_vetor(nvar);

  atualiza_Rede(
      grafo,
      numeroBarras); // atualiza a condição da rede conforme o estado atual
  atualiza_Modelo(
      grafo, numeroBarras, nmed,
      medidas); // atualiza modelo de medição conforme a condição atual da rede
  atualiza_Modelo(grafo, numeroBarras, nvir, virtuais);

  int it = 0;
  int conv = 0;

  while (it < 30) {
    if (conv == 1) {
      FILE *res = NULL;
      res = fopen("residuo.txt", "w+");
      int ct = 0;
      for (int i = 0; i < nmed + nvir; i++) {
        if (ct < nmed) {
          float resid = medidas[ct].zmed - medidas[ct].h;
          fprintf(res, "%f\n", resid);
        } else {
          float resid = virtuais[ct].zmed - virtuais[ct].h;
          fprintf(res, "%f\n", resid);
        }
        ct += 1;
      }
      fclose(res);
      return conv;
    }
    atualiza_H(grafo, numeroBarras, ramos, medidas, nmed);
    atualiza_H(grafo, numeroBarras, ramos, virtuais, nvir);
    double max_sigma = 0;
    for (i = 0; i < nmed; i++) {
      for (j = 0; j < medidas[i].nvar; j++) {
        // H.T * H / (sigma^2)
        medidas[i].H[j] = (medidas[i].H[j]); // /  medidas[i].sigma;
      }
      // H.T*deltaZ/(sigma^2)
      if (medidas[i].sigma > max_sigma) {
        max_sigma = medidas[i].sigma;
      }
      Dz[i] = (medidas[i].zmed - medidas[i].h); /// medidas[i].sigma;
    }
    for (i = 0; i < nvir; i++) {
      Dv[i] = virtuais[i].h;
    }
    double soma;
    float somaDz;
    int _i, _j;
    // printf("max sigma: %f\n", (max_sigma));

    float somaVl = 0;
    for (int i = 0; i < nmed; i++) {
      somaVl += medidas[i].sigma;
    }

    // double vl = (max_sigma);
    // double vl = sqrt(nmed)*somaVl;
    double vl = 1;
    // printf("vl: %f\n", vl);
    for (i = 0; i < nmed; i++) {
      for (j = 0; j < nvar; j++) {
        H_rf[i][j] = *H[i][j];
      }
    }

    cholmod_common Common, *c;

    cholmod_triplet *A_col1_T = NULL;
    cholmod_triplet *A_col2_T = NULL;
    cholmod_triplet *A_col3_T = NULL;

    cholmod_sparse *A_col1 = NULL;
    cholmod_sparse *A_col2 = NULL;
    cholmod_sparse *A_col3 = NULL;

    cholmod_sparse *A_aux_col_12 = NULL;

    cholmod_sparse *A_hatchel = NULL;
    cholmod_dense *b_hatchel = NULL;
    cholmod_dense *X_hatchel = NULL;

    cholmod_triplet *T_hatchel = NULL;

    c = &Common;
    cholmod_l_start(c);
    clock_t tIniMatriz = clock();
    A_hatchel = cholmod_l_allocate_sparse(
        (nmed + nvar + nvir), (nmed + nvar + nvir),
        (nmed + nvar + nvir) * (nmed + nvar + nvir), 0, 0, 0, CHOLMOD_REAL, c);
    b_hatchel = cholmod_l_allocate_dense((nmed + nvar + nvir), 1,
                                         (nmed + nvar + nvir), CHOLMOD_REAL, c);
    X_hatchel = cholmod_l_allocate_dense((nmed + nvar + nvir), 1,
                                         (nmed + nvar + nvir), CHOLMOD_REAL, c);

    A_col1 = cholmod_l_allocate_sparse((nmed + nvar + nvir), nmed,
                                       ((nmed + nvar + nvir) * nmed), 0, 0, 0,
                                       CHOLMOD_REAL, c);
    A_col2 = cholmod_l_allocate_sparse((nmed + nvar + nvir), nvar,
                                       ((nmed + nvar + nvir) * nvar), 0, 0, 0,
                                       CHOLMOD_REAL, c);
    A_col3 = cholmod_l_allocate_sparse((nmed + nvar + nvir), nvir,
                                       ((nmed + nvar + nvir) * nvir), 0, 0, 0,
                                       CHOLMOD_REAL, c);

    A_aux_col_12 = cholmod_l_allocate_sparse(
        (nmed + nvar + nvir), (nmed + nvar),
        ((nmed + nvar + nvir) * (nmed + nvar)), 0, 0, 0, CHOLMOD_REAL, c);

    A_col1_T = cholmod_l_allocate_triplet((nmed + nvar + nvir), nmed,
                                          (nmed + nvar + nvir) * nmed, 0,
                                          CHOLMOD_REAL, c);
    A_col2_T = cholmod_l_allocate_triplet((nmed + nvar + nvir), nvar,
                                          (nmed + nvar + nvir) * nvar, 0,
                                          CHOLMOD_REAL, c);
    A_col3_T = cholmod_l_allocate_triplet((nmed + nvar + nvir), nvir,
                                          (nmed + nvar + nvir) * nvir, 0,
                                          CHOLMOD_REAL, c);

    T_hatchel = cholmod_l_allocate_triplet(
        (nmed + nvar + nvir), (nmed + nvar + nvir),
        (nmed + nvar + nvir) * (nmed + nvar + nvir), 0, CHOLMOD_REAL, c);

    // Coluna 1 da matriz aumentada de Hatchel
    int index1 = 0;
    for (int i = 0; i < nmed; i++) {
      ((long int *)A_col1_T->i)[index1] = i;
      ((long int *)A_col1_T->j)[index1] = i;
      ((double *)A_col1_T->x)[index1] = medidas[i].sigma;
      A_col1_T->nnz += 1;
      index1 += 1;
    }
    for (int i = 0; i < nvar; i++) {
      for (int j = 0; j < nmed; j++) {
        if (*H[j][i] != 0) {
          ((long int *)A_col1_T->i)[index1] = i + nmed;
          ((long int *)A_col1_T->j)[index1] = j;
          ((double *)A_col1_T->x)[index1] = *H[j][i];
          A_col1_T->nnz += 1;
          index1 += 1;
        }
      }
    }

    // Coluna 2 da matriz aumentada de Hatchel
    int index2 = 0;
    for (int i = 0; i < nmed; i++) {
      for (int j = 0; j < nvar; j++) {
        if (*H[i][j] != 0) {
          ((long int *)A_col2_T->i)[index2] = i;
          ((long int *)A_col2_T->j)[index2] = j;
          ((double *)A_col2_T->x)[index2] = *H[i][j];
          A_col2_T->nnz += 1;
          index2 += 1;
        }
      }
    }
    for (int i = 0; i < nvir; i++) {
      for (int j = 0; j < nvar; j++) {
        if (*C[i][j] != 0) {
          ((long int *)A_col2_T->i)[index2] = i + nmed + nvar;
          ((long int *)A_col2_T->j)[index2] = j;
          ((double *)A_col2_T->x)[index2] = *C[i][j];
          A_col2_T->nnz += 1;
          index2 += 1;
        }
      }
    }

    // Coluna 3 da matriz aumentada de Hatchel
    int index3 = 0;
    for (int i = 0; i < nvar; i++) {
      for (int j = 0; j < nvir; j++) {
        if (*C[j][i] != 0) {
          ((long int *)A_col3_T->i)[index3] = i + nmed;
          ((long int *)A_col3_T->j)[index3] = j;
          ((double *)A_col3_T->x)[index3] = *C[j][i];
          A_col3_T->nnz += 1;
          index3 += 1;
        }
      }
    }

    printf("nmed: %d, nvar: %d, nvir: %d\n", nmed, nvar, nvir);
    A_col1 = cholmod_l_triplet_to_sparse(A_col1_T, A_col1_T->nnz, c);

    A_col2 = cholmod_l_triplet_to_sparse(A_col2_T, A_col2_T->nnz, c);

    A_col3 = cholmod_l_triplet_to_sparse(A_col3_T, A_col3_T->nnz, c);

    A_aux_col_12 = cholmod_l_horzcat(A_col1, A_col2, 1, c);
    A_hatchel = cholmod_l_horzcat(A_aux_col_12, A_col3, 1, c);
    T_hatchel = cholmod_l_sparse_to_triplet(A_hatchel, c);

    BOOL write = false;

    if (write && it == 0) {
      FILE *matNEC;
      matNEC = fopen("mathatchel_342_alfa0.txt", "w+");
      for (i = 0; i < T_hatchel->nnz; i++) {
        long int _i, _j;
        double v;
        _i = ((long int *)T_hatchel->i)[i];
        _j = ((long int *)T_hatchel->j)[i];
        v = ((double *)T_hatchel->x)[i];
        fprintf(matNEC, "%d,%d,%f\n", _i, _j, v);
      }
    }

    // Lado direito - formulação com matriz aumentada de Hatchel
    for (int i = 0; i < nmed; i++) {
      ((double *)b_hatchel->x)[i] = (medidas[i].zmed - medidas[i].h);
    }
    for (int i = 0; i < nvar; i++) {
      ((double *)b_hatchel->x)[i + nmed] = 0;
    }
    for (int i = 0; i < nvir; i++) {
      ((double *)b_hatchel->x)[i + nmed + nvar] = -virtuais[i].h;
    }

    clock_t t1Matriz = clock();
    double tempoMatriz = (double)(t1Matriz - tIniMatriz) / CLOCKS_PER_SEC;
    printf("\nAlocacao e construcao das matrizes: %lf", tempoMatriz);

    // TODO: exportar matriz para plot spy
    // if (write && it == 0)
    //{
    //     FILE *matHat;
    //     matHat = fopen("matHatchel_123_alfa0.txt", "w+");
    //     for (i = 0; i < T_nec->nnz; i++)
    //     {
    //         long int _i, _j;
    //         double v;
    //         _i = ((long int *)T_nec->i)[i];
    //         _j = ((long int *)T_nec->j)[i];
    //         v = ((double *)T_nec->x)[i];
    //         fprintf(matHat, "%d,%d,%f\n", _i, _j, v);
    //     }
    // }
    clock_t tIni = clock();
    X_hatchel = SuiteSparseQR_C_backslash(SPQR_ORDERING_BEST, SPQR_DEFAULT_TOL,
                                          A_hatchel, b_hatchel, c);
    clock_t t1 = clock();
    double tempoWLS = (double)(t1 - tIni) / CLOCKS_PER_SEC;
    printf("\nResolucao sistema linear: %lf", tempoWLS);
    Dx = (double *)X_hatchel->x;

    for (i = 0; i < nvar; i++) {
      ponto[i] = ponto[i] + Dx[i + nmed];
    }
    for (i = 0; i < nvar; i++) {
      if (fabs(Dx[i + nmed]) >= tol) {
        conv = 0;
        break;
      } else
        conv = 1;
    }

    atualiza_estado(
        grafo, ponto, regua_comp,
        nvar); // atualiza o estado atual do grafo conforme o vetor x calculado
    atualiza_Rede(grafo, numeroBarras);
    atualiza_Modelo(grafo, numeroBarras, nmed, medidas);
    atualiza_Modelo(grafo, numeroBarras, nvir, virtuais);

    nFx = norma_inf(Dx, nvar);
    // nGx = norma_inf(b, nvar);
    printf("\n\nIteracao:  %d \t|Dx|_inf =  %.17lf \t", it, nFx);

    it++;
    printf(".");
  }
  atualiza_estado(
      grafo, ponto, regua_comp,
      nvar); // atualiza o estado atual do grafo conforme o vetor x calculado
  atualiza_Rede(grafo, numeroBarras);
  atualiza_Modelo(grafo, numeroBarras, nmed, medidas);
  atualiza_Modelo(grafo, numeroBarras, nvir, virtuais);
  // atualiza_h(grafo, numeroBarras, nmed, medidas);
  // atualiza_h(grafo, numeroBarras, nvir, virtuais);

  return conv;
}
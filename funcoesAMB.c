/*
 * File:   funcoesBranchCurrent.c
 * Author: Gustavo Hebling
 *
 * Created on 15 de Julho de 2022, 08:45
 */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "data_structures.h"
#include "funcoesBranchCurrent.h"
#include "funcoesCalculoEletrico.h"
#include "funcoesFluxoVarredura.h"
#include "funcoesMatematicas.h"
#include "funcoesOtimizacao.h"
#include "funcoesTopologia.h"
#include "funcoesWLS.h"

// #include "mmio.h"
#include "SuiteSparseQR_C.h"
#include <cholmod.h>

#include "data_structures.h"

double *monta_regua_x_AMB(GRAFO *grafo, long int numeroBarras) {
  double *regua_x;
  regua_x = (double *)malloc(numeroBarras * sizeof(double));

  for (int i = 0; i < numeroBarras; i++) {
    regua_x[i] = grafo[i].barra->ID;
  }
  return regua_x;
}

double *monta_vetor_x_inicial_AMB(long int numeroBarras, GRAFO *grafo) {
  double *x_init;
  x_init = (double *)malloc(6 * numeroBarras * sizeof(double));
  for (int i = 0; i < numeroBarras; i++) {
    x_init[3 * i] = creal(grafo[i].V[0]);
    x_init[3 * i + (3 * numeroBarras)] = cimag(grafo[i].V[0]);

    x_init[3 * i + 1] = creal(grafo[i].V[1]);
    x_init[3 * i + 1 + (3 * numeroBarras)] = cimag(grafo[i].V[1]);

    x_init[3 * i + 2] = creal(grafo[i].V[2]);
    x_init[3 * i + 2 + (3 * numeroBarras)] = cimag(grafo[i].V[2]);
  }

  return x_init;
}

DMED_COMPLEX *calcula_medida_tensao_complexa_AMB(DMED *medidas,
                                                 long int numeroMedidas,
                                                 GRAFO *grafo,
                                                 int numeroBarras) {
  DMED_COMPLEX *medida_tensao = NULL;

  int cont_med_tensao;
  cont_med_tensao = conta_medidas_Tensao(medidas, numeroMedidas);

  printf("nmed_tensao = %d\n", cont_med_tensao);
  if (((medida_tensao) = (DMED_COMPLEX *)malloc(
           (cont_med_tensao) * sizeof(DMED_COMPLEX))) == NULL) {
    printf("Erro -- Nao foi possivel alocar espaco de memoria para as medidas "
           "complexas !!!!");
    exit(1);
  }

  int aux_contador = 0;
  int med_found = 0;
  for (int contador = 0; contador < numeroMedidas; contador++) {
    if (medidas[contador].tipo == 4) {
      medida_tensao[aux_contador].ligado = medidas[contador].ligado;
      medida_tensao[aux_contador].tipo = medidas[contador].tipo;
      medida_tensao[aux_contador].DE = medidas[contador].DE;
      medida_tensao[aux_contador].PARA = medidas[contador].PARA;
      medida_tensao[aux_contador].fases = medidas[contador].fases;
      medida_tensao[aux_contador].id = medidas[contador].id;
      medida_tensao[aux_contador].par = medidas[contador].par;
      medida_tensao[aux_contador].k = medidas[contador].k;
      medida_tensao[aux_contador].m = medidas[contador].m;
      medida_tensao[aux_contador].ramo = medidas[contador].ramo;
      medida_tensao[aux_contador].idAlim = medidas[contador].idAlim;
      medida_tensao[aux_contador].idArea = medidas[contador].idArea;
      medida_tensao[aux_contador].idFront = medidas[contador].idFront;

      medida_tensao[aux_contador].sigma = medidas[contador].sigma;
      medida_tensao[aux_contador].prec = medidas[contador].prec;

      double magnitude_tensao = medidas[contador].zmed;

      int idx_barra_DE = medida_tensao[aux_contador].k;

      int f_idx = 0;

      double phi_k;
      switch (medida_tensao[aux_contador].fases) {
      case 1:
        f_idx = 0;
        phi_k = 0.0;
        break;
      case 2:
        f_idx = 1;
        phi_k = -2 * PI / 3;
        break;
      case 3:
        f_idx = 2;
        phi_k = 2 * PI / 3;
        break;
      }
      __complex__ double tensao_de = grafo[idx_barra_DE].V[f_idx];
      double abs_tensao_de = sqrt(creal(tensao_de) * creal(tensao_de) +
                                  cimag(tensao_de) * cimag(tensao_de));

      double parte_real = creal(tensao_de) / abs_tensao_de;
      double parte_imag = cimag(tensao_de) / abs_tensao_de;

      double ang = atan(parte_imag / parte_real);
      //       medida_tensao[aux_contador].zmed =
      //           (magnitude_tensao * parte_real) + (magnitude_tensao *
      //           parte_imag) * I;

      // medida_tensao[aux_contador].zmed = abs_tensao_de + I*0;
      medida_tensao[aux_contador].zmed =
          ((magnitude_tensao * parte_real) * cos(phi_k) +
           (magnitude_tensao * parte_imag) * -sin(phi_k)) +
          I * 0;

      aux_contador += 1;
    }
  }
  return medida_tensao;
}

DMED_COMPLEX *converte_medidas_para_complexo_AMB(DMED *medidas,
                                                 long int numeroMedidas) {
  DMED_COMPLEX *medidas_complexas = NULL;

  int cont_med_complex;
  cont_med_complex = conta_medidas_BC(medidas, numeroMedidas);

  printf("nmed_complex = %d\n", cont_med_complex);
  if (((medidas_complexas) = (DMED_COMPLEX *)malloc(
           (cont_med_complex) * sizeof(DMED_COMPLEX))) == NULL) {
    printf("Erro -- Nao foi possivel alocar espaco de memoria para as medidas "
           "complexas !!!!");
    exit(1);
  }

  int aux_contador = 0;
  int med_found = 0;
  for (int contador = 0; contador < numeroMedidas; contador++) {
    if (medidas[contador].tipo == 0 || medidas[contador].tipo == 2) {
      medidas_complexas[aux_contador].DE = medidas[contador].DE;
      medidas_complexas[aux_contador].PARA = medidas[contador].PARA;
      medidas_complexas[aux_contador].id = medidas[contador].id;
      medidas_complexas[aux_contador].k = medidas[contador].k;
      medidas_complexas[aux_contador].m = medidas[contador].m;
      medidas_complexas[aux_contador].ramo = medidas[contador].ramo;
      medidas_complexas[aux_contador].fases = medidas[contador].fases;
      medidas_complexas[aux_contador].tipo = medidas[contador].tipo;
      medidas_complexas[aux_contador].idAlim = medidas[contador].idAlim;
      medidas_complexas[aux_contador].idArea = medidas[contador].idArea;
      medidas_complexas[aux_contador].idFront = medidas[contador].idFront;
      medidas_complexas[aux_contador].ligado = medidas[contador].ligado;
      medidas_complexas[aux_contador].par = medidas[contador].par;
      // zmed
      medidas_complexas[aux_contador].sigma = medidas[contador].sigma;
      medidas_complexas[aux_contador].prec = medidas[contador].prec;

      double parte_real = medidas[contador].zmed;

      medidas_complexas[aux_contador].zmed = parte_real + 0 * I;

      for (int j = 0; j < numeroMedidas; j++) {
        if (medidas[j].tipo == 1 || medidas[j].tipo == 3) {
          if (medidas[j].DE == medidas[contador].DE &&
              medidas[j].PARA == medidas[contador].PARA &&
              medidas[j].fases == medidas[contador].fases) {
            double parte_imag = medidas[j].zmed;

            medidas_complexas[aux_contador].zmed = parte_real + parte_imag * I;
            med_found += 1;
            break;
          }
        }
      }

      // printf("medida %d: %f + i*%f\n", aux_contador,
      // creal(medidas_complexas[aux_contador].zmed)),
      // cimag(medidas_complexas[aux_contador].zmed);
      aux_contador += 1;
    }
  }
  // printf("med_found: %d\n", med_found);

  return medidas_complexas;
}

double calcula_G_AMB(int fase, int i_ramo, DRAM *ramos) {
  double R, X, G;
  if (ramos[i_ramo].Z != NULL) {
    R = creal(ramos[i_ramo].Z[fase][fase]);
    X = cimag(ramos[i_ramo].Z[fase][fase]);
    G = R / ((R * R) + (X * X));
    return G;
  } else {
    G = 0.1;
    return G;
  }
}
double calcula_B_AMB(int fase, int i_ramo, DRAM *ramos) {
  double R, X;
  double B;
  if (ramos[i_ramo].Z != NULL) {
    R = creal(ramos[i_ramo].Z[fase][fase]);
    X = cimag(ramos[i_ramo].Z[fase][fase]);
    B = -X / ((R * R) + (X * X));
    return B;
  } else {
    B = 0.1;
    return B;
  }
}

double **monta_matriz_H_AMB(long int numeroBarras, long int numeroRamos,
                            int nmed_AMB, DMED_COMPLEX *medidas_equivalentes,
                            double *regua_x, DRAM *ramos, GRAFO *grafo) {
  //__complex__ double **H_BC = NULL;
  double **H_T = NULL;
  double *hx_I_2 = NULL;

  H_T = aloca_matriz(2 * nmed_AMB, 6 * numeroBarras);
  hx_I_2 = aloca_vetor(2 * nmed_AMB);

  __complex__ double ramo_Z;
  double R = 0;
  double X = 0;
  double G = 0;
  double B = 0;
  int i, j;

  for (i = 0; i < nmed_AMB; i++) {
    long int DE = medidas_equivalentes[i].DE;
    long int PARA = medidas_equivalentes[i].PARA;
    long int i_ramo = medidas_equivalentes[i].ramo;
    for (j = 0; j < numeroBarras; j++) {
      int x_atual = regua_x[j];

      // x_atual = ID da barra atual na regua
      // se a medida de potencia conter a barra (DE ou PARA) - calcular o Y para
      // a fase

      // é medida de fluxo
      if (medidas_equivalentes[i].tipo == 0 ||
          medidas_equivalentes[i].tipo == 1) {
        if ((medidas_equivalentes[i].k == j) ||
            (medidas_equivalentes[i].m == j)) {
          double sg = 1.0;
          if ((medidas_equivalentes[i].k == j)) {
            sg = 1;
          } else {
            sg = -1;
          }

          switch (medidas_equivalentes[i].fases) {
          case 1:
            // calcula Y
            // FASE A
            // IR/VR
            // G = calcula_G_AMB(0,i_ramo,ramos);
            // B = calcula_B_AMB(0,i_ramo,ramos);

            G = creal(ramos[i_ramo].Ypp[0][0]);
            B = cimag(ramos[i_ramo].Ypp[0][0]);

            H_T[i][3 * j] = sg * G;
            // IR/VX
            H_T[i][3 * j + (3 * numeroBarras)] = -sg * B;
            // IX/VR
            H_T[i + nmed_AMB][3 * j] = sg * B;
            // IX/VX
            H_T[i + nmed_AMB][(3 * j) + (3 * numeroBarras)] = sg * G;

            break;
          case 2:
            // calcula Y
            // FASE B
            // G = calcula_G_AMB(1,i_ramo,ramos);
            // B = calcula_B_AMB(1,i_ramo,ramos);

            G = creal(ramos[i_ramo].Ypp[1][1]);
            B = cimag(ramos[i_ramo].Ypp[1][1]);

            H_T[i][(3 * j) + 1] = sg * G;
            // IR/VX
            H_T[i][(3 * j) + 1 + (3 * numeroBarras)] = -sg * B;
            // IX/VR
            H_T[i + nmed_AMB][(3 * j) + 1] = sg * B;
            // IX/VX
            H_T[i + nmed_AMB][(3 * j) + 1 + (3 * numeroBarras)] = sg * G;

            break;
          case 3:
            // calcula Y
            // FASE C
            // G = calcula_G_AMB(2,i_ramo,ramos);
            // B = calcula_B_AMB(2,i_ramo,ramos);

            G = creal(ramos[i_ramo].Ypp[2][2]);
            B = cimag(ramos[i_ramo].Ypp[2][2]);

            H_T[i][(3 * j) + 2] = sg * G;
            // IR/VX
            H_T[i][(3 * j) + 2 + (3 * numeroBarras)] = -sg * B;
            // IX/VR
            H_T[i + nmed_AMB][(3 * j) + 2] = sg * B;
            // IX/VX
            H_T[i + nmed_AMB][(3 * j) + 2 + (3 * numeroBarras)] = sg * G;

            break;
          case 4:
            break;
          case 5:
            break;
          case 6:
            break;
          case 7:
            break;
          }
        }
      }
      // é medida de injecao
      if (medidas_equivalentes[i].tipo == 2 ||
          medidas_equivalentes[i].tipo == 3) {
        // TODO
        // medida i e barra j
        // percorre os adjacentes da barra i e caso encontre a barra da medida
        // j, associa a impedancia na matriz
        int barra_atual = medidas_equivalentes[i].k;
        for (int k = 0; k < grafo[barra_atual].numeroAdjacentes; k++) {
          if (j == grafo[barra_atual].adjacentes[k].idNo) {
            switch (medidas_equivalentes[i].fases) {
            case 1:
              if ((grafo[barra_atual].adjacentes[k].ramo->Ypp) == NULL) {
                G = 1;
                B = 1;
                printf("NULL\n");
              } else {
                // G = creal(1/(grafo[j].adjacentes[k].ramo->Z[0][0]));
                // B = -cimag(1/(grafo[j].adjacentes[k].ramo->Z[0][0]));

                G = -creal(grafo[barra_atual].adjacentes[k].ramo->Ypp[0][0]);
                B = -cimag(grafo[barra_atual].adjacentes[k].ramo->Ypp[0][0]);
              }
              H_T[i][3 * j] = G;
              // IR/VX
              H_T[i][3 * j + (3 * numeroBarras)] = -B;
              // IX/VR
              H_T[i + (nmed_AMB)][3 * j] = B;
              // IX/VX
              H_T[i + (nmed_AMB)][3 * j + (3 * numeroBarras)] = G;
              break;

            case 2:
              if ((grafo[barra_atual].adjacentes[k].ramo->Ypp) == NULL) {
                G = 1;
                B = 1;
                printf("NULL\n");
              } else {
                // G = creal(1/(grafo[j].adjacentes[k].ramo->Z[1][1]));
                // B = -cimag(1/(grafo[j].adjacentes[k].ramo->Z[1][1]));

                G = -creal(grafo[barra_atual].adjacentes[k].ramo->Ypp[1][1]);
                B = -cimag(grafo[barra_atual].adjacentes[k].ramo->Ypp[1][1]);
              }
              H_T[i][3 * j + 1] = G;
              // IR/VX
              H_T[i][3 * j + 1 + (3 * numeroBarras)] = -B;
              // IX/VR
              H_T[i + (nmed_AMB)][3 * j + 1] = B;
              // IX/VX
              H_T[i + (nmed_AMB)][3 * j + 1 + (3 * numeroBarras)] = G;
              break;
            case 3:

              if ((grafo[barra_atual].adjacentes[k].ramo->Ypp) == NULL) {
                G = 1;
                B = 1;
              } else {
                // G = creal(1/(grafo[j].adjacentes[k].ramo->Z[2][2]));
                // B = cimag(1/(grafo[j].adjacentes[k].ramo->Z[2][2]));

                G = -creal(grafo[barra_atual].adjacentes[k].ramo->Ypp[2][2]);
                B = -cimag(grafo[barra_atual].adjacentes[k].ramo->Ypp[2][2]);

                // IR/VR
                H_T[i][3 * j + 2] = G;
                // IR/VX
                H_T[i][3 * j + 2 + (3 * numeroBarras)] = -B;
                // IX/VR
                H_T[i + (nmed_AMB)][3 * j + 2] = B;
                // IX/VX
                H_T[i + (nmed_AMB)][3 * j + 2 + (3 * numeroBarras)] = G;
              }
              break;
            }
          }
        }
      }
    }
  }

  double soma1, soma2, soma3, soma4;
  for (i = 0; i < nmed_AMB; i++) {
    if (medidas_equivalentes[i].tipo == 2 ||
        medidas_equivalentes[i].tipo == 3) {
      soma1 = 0.0;
      soma2 = 0.0;
      soma3 = 0.0;
      soma4 = 0.0;

      for (j = 0; j < numeroBarras; j++) {
        switch (medidas_equivalentes[i].fases) {
        case 1:
          soma1 += H_T[i][3 * j];
          // IR/VX
          soma2 += H_T[i][3 * j + (3 * numeroBarras)];
          // IX/VR
          soma3 += H_T[i + nmed_AMB][3 * j];
          // IX/VX
          soma4 += H_T[i + nmed_AMB][3 * j + (3 * numeroBarras)];
          break;

        case 2:
          soma1 += H_T[i][3 * j + 1];
          // IR/VX
          soma2 += H_T[i][3 * j + 1 + (3 * numeroBarras)];
          // IX/VR
          soma3 += H_T[i + nmed_AMB][3 * j + 1];
          // IX/VX
          soma4 += H_T[i + nmed_AMB][3 * j + 1 + (3 * numeroBarras)];
          break;

        case 3:
          soma1 += H_T[i][3 * j + 2];
          // IR/VX
          soma2 += H_T[i][3 * j + 2 + (3 * numeroBarras)];
          // IX/VR
          soma3 += H_T[i + nmed_AMB][3 * j + 2];
          // IX/VX
          soma4 += H_T[i + nmed_AMB][3 * j + 2 + (3 * numeroBarras)];
          break;
        }
      }
      int k;
      k = medidas_equivalentes[i].k;

      switch (medidas_equivalentes[i].fases) {

      case 1:
        H_T[i][3 * k] = -soma1;
        // IR/
        H_T[i][3 * k + (3 * numeroBarras)] = -soma2;
        //        //IX/
        H_T[i + nmed_AMB][3 * k] = -soma3;
        // IX/
        H_T[i + nmed_AMB][3 * k + (3 * numeroBarras)] = -soma4;
        break;
      case 2:
        H_T[i][3 * k + 1] = -soma1;
        // IR/
        H_T[i][3 * k + 1 + (3 * numeroBarras)] = -soma2;
        // IX/
        H_T[i + nmed_AMB][3 * k + 1] = -soma3;
        // IX/
        H_T[i + nmed_AMB][3 * k + 1 + (3 * numeroBarras)] = -soma4;
        break;
      case 3:
        H_T[i][3 * k + 2] = -soma1;
        // IR/
        H_T[i][3 * k + 2 + (3 * numeroBarras)] = -soma2;
        // IX/
        H_T[i + nmed_AMB][3 * k + 2] = -soma3;
        // IX/
        H_T[i + nmed_AMB][3 * k + 2 + (3 * numeroBarras)] = -soma4;
        break;
      }
    }
  }

  return H_T;
}

double **monta_matriz_H_AMB_Tensao(DMED_COMPLEX *medidas_tensao, GRAFO *grafo,
                                   long int numeroBarras, int nmed_T) {
  double **H_Tensao;
  H_Tensao = aloca_matriz(2 * nmed_T, 6 * numeroBarras);

  int i, j;
  double phi_k, ang_tensao, abs_tensao;
  for (i = 0; i < nmed_T; i++) {
    int DE_med = medidas_tensao[i].k;
    for (j = 1; j < numeroBarras; j++) {
      // TODO
      // int DE_barra = grafo[j].barra->ID;
      if (DE_med == j && j != 0) {
        switch (medidas_tensao[i].fases) {
        case 1:
          /* code */
          abs_tensao =
              sqrt(creal(grafo[DE_med].V[0]) * creal(grafo[DE_med].V[0]) +
                   cimag(grafo[DE_med].V[0]) * cimag(grafo[DE_med].V[0]));
          phi_k = 0.0;
          H_Tensao[i][3 * j] =
              cos(phi_k); // creal(grafo[DE_med].V[0])/abs_tensao;
          H_Tensao[i][3 * j + 3 * numeroBarras] =
              -sin(phi_k); // creal(grafo[DE_med].V[0])/abs_tensao;
          break;
        case 2:
          /* code */
          abs_tensao =
              sqrt(creal(grafo[DE_med].V[1]) * creal(grafo[DE_med].V[1]) +
                   cimag(grafo[DE_med].V[1]) * cimag(grafo[DE_med].V[1]));
          phi_k = -2 * PI / 3;
          H_Tensao[i][3 * j + 1] =
              cos(phi_k); // creal(grafo[DE_med].V[1])/abs_tensao;
          H_Tensao[i][3 * j + 1 + 3 * numeroBarras] =
              -sin(phi_k); // creal(grafo[DE_med].V[1])/abs_tensao;
          break;
        case 3:
          /* code */
          abs_tensao =
              sqrt(creal(grafo[DE_med].V[2]) * creal(grafo[DE_med].V[2]) +
                   cimag(grafo[DE_med].V[2]) * cimag(grafo[DE_med].V[2]));
          phi_k = 2 * PI / 3;
          H_Tensao[i][3 * j + 2] =
              cos(phi_k); // creal(grafo[DE_med].V[2])/abs_tensao;
          H_Tensao[i][3 * j + 2 + 3 * numeroBarras] =
              -sin(phi_k); // creal(grafo[DE_med].V[2])/abs_tensao;
          break;
        }
      }
      //      if (DE_med == j && j == 0){
      //        switch (medidas_tensao[i].fases){
      //            case 1:
      //              H_Tensao[i][3 * j] = 1;
      //              H_Tensao[i][(3 * j) + (3 * numeroBarras)] = 1;
      //              break;
      //            case 2:
      //              H_Tensao[i][3 * j] = 1;
      //              H_Tensao[i][(3 * j) + 1 + (3 * numeroBarras)] = 1;
      //              break;
      //            case 3:
      //              H_Tensao[i][3 * j] = 1;
      //              H_Tensao[i][(3 * j) + 2 + (3 * numeroBarras)] = 1;
      //              break;
      //        }
      //      }
    }
  }

  return H_Tensao;
}

double *monta_hx_V(long int nmed_T, DMED_COMPLEX *medidas_tensao, GRAFO *grafo,
                   long int numeroBarras) {
  double *hx_v;
  hx_v = aloca_vetor(2 * nmed_T);

  int i;
  double tensao_barra;
  double tensao_imag;
  double phi_k;
  double abs_tensao;
  for (i = 0; i < nmed_T; i++) {
    int barra_de = medidas_tensao[i].DE;
    int i_grafo = medidas_tensao[i].k;
    // int i_grafo = calcula_idx_ID(grafo, barra_de, numeroBarras);

    switch (medidas_tensao[i].fases) {
    case 1:
      // tensao_barra =
      // sqrt(creal(grafo[i_grafo].V[0])*creal(grafo[i_grafo].V[0]) +
      // cimag(grafo[i_grafo].V[0])*cimag(grafo[i_grafo].V[0]));
      phi_k = 0;
      tensao_barra = creal(grafo[i_grafo].V[0]);
      tensao_imag = cimag(grafo[i_grafo].V[0]);
      abs_tensao =
          sqrt(tensao_barra * tensao_barra + tensao_imag * tensao_imag);
      hx_v[i] = tensao_barra * cos(phi_k) + tensao_imag * -sin(phi_k);
      break;
    case 2:
      // tensao_barra =
      // sqrt(creal(grafo[i_grafo].V[1])*creal(grafo[i_grafo].V[1]) +
      // cimag(grafo[i_grafo].V[1])*cimag(grafo[i_grafo].V[1]));
      phi_k = -2 * PI / 3;
      tensao_barra = creal(grafo[i_grafo].V[1]);
      tensao_imag = cimag(grafo[i_grafo].V[1]);
      abs_tensao =
          sqrt(tensao_barra * tensao_barra + tensao_imag * tensao_imag);
      hx_v[i] = tensao_barra * cos(phi_k) + tensao_imag * -sin(phi_k);
      break;
    case 3:
      // tensao_barra =
      // sqrt(creal(grafo[i_grafo].V[2])*creal(grafo[i_grafo].V[2]) +
      // cimag(grafo[i_grafo].V[2])*cimag(grafo[i_grafo].V[2]));
      phi_k = 2 * PI / 3;
      tensao_barra = creal(grafo[i_grafo].V[2]);
      tensao_imag = cimag(grafo[i_grafo].V[2]);
      abs_tensao =
          sqrt(tensao_barra * tensao_barra + tensao_imag * tensao_imag);
      hx_v[i] = tensao_barra * cos(phi_k) + tensao_imag * -sin(phi_k);
      break;
    }
    //    hx_v[i] = abs_tensao;
    //    hx_v[i + nmed_T] = tensao_imag;
  }
  return hx_v;
}

double *monta_hx_I(int nmed_AMB, DMED_COMPLEX *medidas_equivalentes,
                   GRAFO *grafo) {
  double *hx_i;
  hx_i = aloca_vetor(2 * nmed_AMB);
  int de, para;
  int n_adj;
  int cont_adj;
  int i, j;

  double soma_cur_real = 0.0;
  double soma_cur_imag = 0.0;
  // para cada medida é necessário buscar a corrente no ramo DE-PARA
  // ou a soma das correntes no caso de injecao
  // assumindo que as correntes foram atualizadas via processo backwards
  for (i = 0; i < nmed_AMB; i++) {
    soma_cur_real = 0.0;
    soma_cur_imag = 0.0;
    if (medidas_equivalentes[i].tipo == 0 ||
        medidas_equivalentes[i].tipo == 1) {
      // para medidas de fluxo
      de = medidas_equivalentes[i].k;
      para = medidas_equivalentes[i].m;
      n_adj = grafo[de].numeroAdjacentes;
      for (cont_adj = 0; cont_adj < n_adj; cont_adj++) {
        if (grafo[de].adjacentes[cont_adj].idNo == para) {
          switch (medidas_equivalentes[i].fases) {
          case 1:
            hx_i[i] = creal(grafo[de].adjacentes[cont_adj].Cur[0]);
            hx_i[i + nmed_AMB] = cimag(grafo[de].adjacentes[cont_adj].Cur[0]);
            break;

          case 2:
            hx_i[i] = creal(grafo[de].adjacentes[cont_adj].Cur[1]);
            hx_i[i + nmed_AMB] = cimag(grafo[de].adjacentes[cont_adj].Cur[1]);
            break;

          case 3:
            hx_i[i] = creal(grafo[de].adjacentes[cont_adj].Cur[2]);
            hx_i[i + nmed_AMB] = cimag(grafo[de].adjacentes[cont_adj].Cur[2]);
            break;
          }
        }
      }

    } else {
      // para medidas de injecao
      de = medidas_equivalentes[i].k;
      n_adj = grafo[de].numeroAdjacentes;
      for (cont_adj = 0; cont_adj < n_adj; cont_adj++) {
        switch (medidas_equivalentes[i].fases) {
        case 1:
          soma_cur_real += creal(grafo[de].adjacentes[cont_adj].Cur[0]);
          soma_cur_imag += cimag(grafo[de].adjacentes[cont_adj].Cur[0]);
          break;

        case 2:
          soma_cur_real += creal(grafo[de].adjacentes[cont_adj].Cur[1]);
          soma_cur_imag += cimag(grafo[de].adjacentes[cont_adj].Cur[1]);
          break;

        case 3:
          soma_cur_real += creal(grafo[de].adjacentes[cont_adj].Cur[2]);
          soma_cur_imag += cimag(grafo[de].adjacentes[cont_adj].Cur[2]);
          break;
        }
      }
      switch (medidas_equivalentes[i].fases) {
      case 1:
        hx_i[i] = soma_cur_real;
        hx_i[i + nmed_AMB] = soma_cur_imag;
        break;

      case 2:
        hx_i[i] = soma_cur_real;
        hx_i[i + nmed_AMB] = soma_cur_imag;
        break;

      case 3:
        hx_i[i] = soma_cur_real;
        hx_i[i + nmed_AMB] = soma_cur_imag;
        break;
      }
    }
  }
  return hx_i;
}

double *resolve_linear_QR_AMB(double **H_AMB, double **H_T, double *z,
                              long int numeroBarras, long int nmed_AMB,
                              int nmed_T, double *hx_I, double *hx_V,
                              double **W) {
  cholmod_sparse *A = NULL;
  cholmod_sparse *AT = NULL;
  cholmod_sparse *G = NULL;
  cholmod_sparse *WH = NULL;

  cholmod_sparse *W_sparse = NULL;
  cholmod_triplet *W_triplet = NULL;

  cholmod_triplet *T = NULL;
  cholmod_dense *b = NULL;
  cholmod_dense *bH = NULL;
  cholmod_dense *X = NULL;
  cholmod_common Common, *c;
  c = &Common;
  cholmod_l_start(c);
  T = cholmod_l_allocate_triplet(2 * nmed_AMB + 2 * nmed_T, 6 * numeroBarras,
                                 (2 * nmed_AMB + 2 * nmed_T) * 6 * numeroBarras,
                                 0, CHOLMOD_REAL, c);
  W_triplet = cholmod_l_allocate_triplet(
      2 * nmed_AMB + 2 * nmed_T, 2 * nmed_AMB + 2 * nmed_T,
      (2 * nmed_AMB + 2 * nmed_T) * 6 * numeroBarras, 0, CHOLMOD_REAL, c);

  // G = cholmod_l_allocate_sparse(6 * numeroRamos,6 * numeroRamos, (6 *
  // numeroRamos* 6*numeroRamos), 0, 0, 1, CHOLMOD_REAL, c);
  A = cholmod_l_allocate_sparse(2 * nmed_AMB + 2 * nmed_T, 6 * numeroBarras,
                                (2 * nmed_AMB + 2 * nmed_T) * 6 * numeroBarras,
                                0, 0, 0, CHOLMOD_REAL, c);
  W_sparse = cholmod_l_allocate_sparse(
      2 * nmed_AMB + 2 * nmed_T, 2 * nmed_AMB + 2 * nmed_T,
      (2 * nmed_AMB + 2 * nmed_T) * 6 * numeroBarras, 0, 0, 0, CHOLMOD_REAL, c);
  // AT = cholmod_l_allocate_sparse( 6 * numeroRamos,6 * nmed_BC + nmed_T, (6 *
  // nmed_BC + nmed_T) * 3 * numeroRamos, 0, 0, 0, CHOLMOD_REAL, c);
  b = cholmod_l_allocate_dense(2 * nmed_AMB + 2 * nmed_T, 1,
                               2 * nmed_AMB + 2 * nmed_T, CHOLMOD_REAL, c);
  bH = cholmod_l_allocate_dense(6 * numeroBarras, 1, 6 * numeroBarras,
                                CHOLMOD_REAL, c);
  X = cholmod_l_allocate_dense(6 * numeroBarras, 1, 6 * numeroBarras,
                               CHOLMOD_REAL, c);
  int index = 0;
  for (int i = 0; i < 2 * nmed_AMB; i++) {
    for (int r = 3; r < 3 * numeroBarras; r++) {
      if ((H_AMB[i][r] != 0)) {
        ((long int *)T->i)[index] = i;
        ((long int *)T->j)[index] = r;
        ((double *)T->x)[index] = H_AMB[i][r];
        T->nnz += 1;
        index += 1;
      }
    }

    for (int r = 3 * numeroBarras + 3; r < 6 * numeroBarras; r++) {
      if ((H_AMB[i][r] != 0)) {
        ((long int *)T->i)[index] = i;
        ((long int *)T->j)[index] = r;
        ((double *)T->x)[index] = H_AMB[i][r];
        T->nnz += 1;
        index += 1;
      }
    }
  }

  for (int t = 0; t < nmed_T; t++) {
    for (int cv = 0; cv < 6 * numeroBarras; cv++) {
      if ((H_T[t][cv] != 0)) {
        ((long int *)T->i)[index] = t + 2 * nmed_AMB;
        ((long int *)T->j)[index] = cv;
        ((double *)T->x)[index] = H_T[t][cv];
        T->nnz += 1;
        index += 1;
      }
    }
  }

  int i_w = 0;
  for (int q = 0; q < (2 * nmed_AMB) + nmed_T; q++) {
    for (int t = 0; t < (2 * nmed_AMB) + nmed_T; t++) {
      if (W[q][t] != 0) {
        ((long int *)W_triplet->i)[i_w] = q;
        ((long int *)W_triplet->j)[i_w] = t;
        ((long int *)W_triplet->x)[i_w] = W[q][t];
        W_triplet->nnz += 1;
        i_w += 1;
      }
    }
  }

  for (int i = 0; i < 2 * nmed_AMB; i++) {
    ((double *)b->x)[(i)] = (z[i] - hx_I[i]);
  }

  for (int i = 0; i < nmed_T; i++) {
    ((double *)b->x)[(i + 2 * nmed_AMB)] = (z[i + 2 * nmed_AMB] - hx_V[i]);
  }

  A = cholmod_l_triplet_to_sparse(
      T, (2 * nmed_AMB + 2 * nmed_T) * 3 * numeroBarras, c);
  W_sparse = cholmod_l_triplet_to_sparse(
      W_triplet, (2 * nmed_AMB + 2 * nmed_T) * (2 * nmed_AMB + 2 * nmed_T), c);

  AT = cholmod_l_transpose(A, 1, c);
  WH = cholmod_l_ssmult(W_sparse, A, 0, 1, 0, c);
  G = cholmod_l_ssmult(AT, A, 0, 1, 0, c);
  double one[2] = {1, 0};
  double m1[2] = {0, 0};
  cholmod_l_sdmult(AT, 0, one, m1, b, bH, c);
//  X = SuiteSparseQR_C_backslash(SPQR_ORDERING_AMD, SPQR_DEFAULT_TOL, A, b, c);
    X = SuiteSparseQR_C_backslash(SPQR_ORDERING_AMD, SPQR_DEFAULT_TOL, G, bH, c);
  // c);
  double *ponto;
  ponto = aloca_vetor(6 * numeroBarras);
  ponto = (double *)X->x;
  // for (int ctz = 0; ctz < 10; ctz ++){
  //         printf("x[%d] = %f + i*%f\n", ctz, ponto[2*ctz], ponto[2*ctz+1]);
  // }
  return ponto;
}

void atualiza_estado_AMB(GRAFO *grafo, long int numeroBarras, double *x_amb) {
  int i;
  double real, imag;

  for (i = 0; i < numeroBarras; i++) {
    real = x_amb[3 * i];
    imag = x_amb[3 * i + 3 * numeroBarras];
    grafo[i].V[0] = real + I * imag;

    real = x_amb[3 * i + 1];
    imag = x_amb[(3 * i) + 1 + (3 * numeroBarras)];
    grafo[i].V[1] = real + I * imag;

    real = x_amb[3 * i + 2];
    imag = x_amb[(3 * i) + 2 + (3 * numeroBarras)];
    grafo[i].V[2] = real + I * imag;
  }
}

double *monta_hx_I_v2(GRAFO *grafo, DMED_COMPLEX *medidas_equivalentes,
                      long int numerobarras, long int nmed_AMB, DRAM *ramos) {
  double *hx_I = NULL;
  hx_I = aloca_vetor(6 * nmed_AMB);

  int i, j, i_ramo, de, para;

  for (i = 0; i < nmed_AMB; i++) {
    if ((medidas_equivalentes[i].tipo == 0) ||
        (medidas_equivalentes[i].tipo == 1)) {
      de = medidas_equivalentes[i].k;
      para = medidas_equivalentes[i].m;
      i_ramo = medidas_equivalentes[i].ramo;

      switch (medidas_equivalentes[i].fases) {
      case 1:
        hx_I[i] = creal(ramos[i_ramo].Ypp[0][0] *
                        (grafo[de].V[0] - grafo[para].V[0]));
        hx_I[i + nmed_AMB] = cimag(ramos[i_ramo].Ypp[0][0] *
                                   (grafo[de].V[0] - grafo[para].V[0]));
        break;
      case 2:
        hx_I[i] = creal(ramos[i_ramo].Ypp[1][1] *
                        (grafo[de].V[1] - grafo[para].V[1]));
        hx_I[i + nmed_AMB] = cimag(ramos[i_ramo].Ypp[1][1] *
                                   (grafo[de].V[1] - grafo[para].V[1]));
        break;
      case 3:
        hx_I[i] = creal(ramos[i_ramo].Ypp[2][2] *
                        (grafo[de].V[2] - grafo[para].V[2]));
        hx_I[i + nmed_AMB] = cimag(ramos[i_ramo].Ypp[2][2] *
                                   (grafo[de].V[2] - grafo[para].V[2]));
        break;
      }

    } else if ((medidas_equivalentes[i].tipo == 2) ||
               (medidas_equivalentes[i].tipo == 3)) {
      int n_adj;
      de = medidas_equivalentes[i].k;
      int ramo_adj;
      int para;
      double soma_real = 0.0;
      double soma_imag = 0.0;
      n_adj = grafo[de].numeroAdjacentes;
      for (j = 0; j < n_adj; j++) {
        ramo_adj = grafo[de].adjacentes[j].idram;
        para = grafo[de].adjacentes[j].idNo;
        switch (medidas_equivalentes[i].fases) {
        case 1:
          soma_real += creal(ramos[ramo_adj].Ypp[0][0] *
                             (grafo[de].V[0] - grafo[para].V[0]));
          soma_imag += cimag(ramos[ramo_adj].Ypp[0][0] *
                             (grafo[de].V[0] - grafo[para].V[0]));
          break;
        case 2:
          soma_real += creal(ramos[ramo_adj].Ypp[1][1] *
                             (grafo[de].V[1] - grafo[para].V[1]));
          soma_imag += cimag(ramos[ramo_adj].Ypp[1][1] *
                             (grafo[de].V[1] - grafo[para].V[1]));
          break;
        case 3:
          soma_real += creal(ramos[ramo_adj].Ypp[2][2] *
                             (grafo[de].V[2] - grafo[para].V[2]));
          soma_imag += cimag(ramos[ramo_adj].Ypp[2][2] *
                             (grafo[de].V[2] - grafo[para].V[2]));
          break;
        }
      }
      switch (medidas_equivalentes[i].fases) {
      case 1:
        hx_I[i] = soma_real;
        hx_I[i + nmed_AMB] = soma_imag;
        break;
      case 2:
        hx_I[i] = soma_real;
        hx_I[i + nmed_AMB] = soma_imag;
        break;
      case 3:
        hx_I[i] = soma_real;
        hx_I[i + nmed_AMB] = soma_imag;
        break;
      }
      /* code */
    }
  }
  return hx_I;
}

double calcula_c_covariancia(int i_medida, int fase, GRAFO *grafo) {
  double parte_real, parte_imag;
  parte_real = creal(grafo[i_medida].V[fase]);
  parte_imag = cimag(grafo[i_medida].V[fase]);

  return parte_real / ((parte_real * parte_real) + (parte_imag * parte_imag));
}

double calcula_d_covariancia(int i_medida, int fase, GRAFO *grafo) {
  double parte_real, parte_imag;
  parte_real = creal(grafo[i_medida].V[fase]);
  parte_imag = cimag(grafo[i_medida].V[fase]);

  return parte_imag / ((parte_real * parte_real) + (parte_imag * parte_imag));
}

double calcula_A(double sigmaP, double sigmaQ, double c, double d) {
  double result;
  result = (c * c) * (sigmaP * sigmaP) + (d * d) * (sigmaQ * sigmaQ);
  return result;
}

double calcula_B(double sigmaP, double sigmaQ, double c, double d) {
  return (c * d) * ((sigmaP * sigmaP) - (sigmaQ * sigmaQ));
}

double calcula_C(double sigmaP, double sigmaQ, double c, double d) {
  return (c * d) * ((sigmaP * sigmaP) - (sigmaQ * sigmaQ));
}

double calcula_D(double sigmaP, double sigmaQ, double c, double d) {
  return (d * d) * (sigmaP * sigmaP) + (c * c) * (sigmaQ * sigmaQ);
}

double **monta_W_AMB(GRAFO *grafo, DMED_COMPLEX *medidas_equivalentes,
                     long int nmed_AMB, DMED_COMPLEX *medidas_tensao,
                     long int nmed_T) {
  double **W = NULL;
  W = aloca_matriz(2 * nmed_AMB + 2 * nmed_T, 2 * nmed_AMB + 2 * nmed_T);
  double val = 0.01;
  int fase;
  double c, d;
  double A, B, C, D;
  double sigmaP = 0.01;
  double sigmaQ = 0.01;
  double det;
  int i;
  int id_barra;
  for (i = 0; i < nmed_AMB; i++) {
    id_barra = medidas_equivalentes[i].k;
    switch (medidas_equivalentes[i].fases) {
    case 1:
      fase = 0;
      c = calcula_c_covariancia(id_barra, fase, grafo);
      d = calcula_d_covariancia(id_barra, fase, grafo);

      sigmaP = medidas_equivalentes[i].sigma;
      sigmaQ = medidas_equivalentes[i].sigma;

      A = calcula_A(sigmaP, sigmaQ, c, d);
      B = calcula_B(sigmaP, sigmaQ, c, d);
      C = calcula_C(sigmaP, sigmaQ, c, d);
      D = calcula_D(sigmaP, sigmaQ, c, d);

      det = 1 / (A * D - B * C);

      W[i][i] = D * det;
      W[i][i + nmed_AMB] = -B * det;
      W[i + nmed_AMB][i] = -C * det;
      W[i + nmed_AMB][i + nmed_AMB] = A * det;
      break;
    case 2:
      fase = 1;
      c = calcula_c_covariancia(id_barra, fase, grafo);
      d = calcula_d_covariancia(id_barra, fase, grafo);

      sigmaP = medidas_equivalentes[i].sigma;
      sigmaQ = medidas_equivalentes[i].sigma;

      A = calcula_A(sigmaP, sigmaQ, c, d);
      B = calcula_B(sigmaP, sigmaQ, c, d);
      C = calcula_C(sigmaP, sigmaQ, c, d);
      D = calcula_D(sigmaP, sigmaQ, c, d);

      det = 1 / (A * D - B * C);
      W[i][i] = D * det;
      W[i][i + nmed_AMB] = -B * det;
      W[i + nmed_AMB][i] = -C * det;
      W[i + nmed_AMB][i + nmed_AMB] = A * det;

      break;
    case 3:
      fase = 2;
      c = calcula_c_covariancia(id_barra, fase, grafo);
      d = calcula_d_covariancia(id_barra, fase, grafo);

      sigmaP = medidas_equivalentes[i].sigma;
      sigmaQ = medidas_equivalentes[i].sigma;

      A = calcula_A(sigmaP, sigmaQ, c, d);
      B = calcula_B(sigmaP, sigmaQ, c, d);
      C = calcula_C(sigmaP, sigmaQ, c, d);
      D = calcula_D(sigmaP, sigmaQ, c, d);

      det = 1 / (A * D - B * C);

      W[i][i] = D * det;
      W[i][i + nmed_AMB] = -B * det;
      W[i + nmed_AMB][i] = -C * det;
      W[i + nmed_AMB][i + nmed_AMB] = A * det;

      break;
    }
  }
  int j;
  for (j = 0; j < nmed_T; j++) {
    switch (medidas_tensao[j].fases) {
    case 1:
      W[2 * nmed_AMB + j][2 * nmed_AMB + j] = val;
      W[2 * nmed_AMB + j + nmed_T][2 * nmed_AMB + j + nmed_T] = val;
      break;
    case 2:
      W[2 * nmed_AMB + j][2 * nmed_AMB + j] = val;
      W[2 * nmed_AMB + j + nmed_T][2 * nmed_AMB + j + nmed_T] = val;

      break;
    case 3:
      W[2 * nmed_AMB + j][2 * nmed_AMB + j] = val;
      W[2 * nmed_AMB + j + nmed_T][2 * nmed_AMB + j + nmed_T] = val;

      break;
    }
    // W[6*nmed_AMB + j][6*nmed_AMB + j] = 1.0;
  }

  return W;
}

void monta_z_real_imag_AMB(DMED_COMPLEX *medidas_eq, double *z_AMB,
                           long int nmed_AMB, DMED_COMPLEX *medidas_tensao,
                           int nmed_T) {
  for (int i = 0; i < nmed_AMB; i++) {
    z_AMB[i] = creal(medidas_eq[i].zmed);
    z_AMB[i + nmed_AMB] = cimag(medidas_eq[i].zmed);
  }
  for (int j = 0; j < nmed_T; j++) {
    z_AMB[j + 2 * nmed_AMB] = creal(medidas_tensao[j].zmed);
    // z_AMB[j + 2 * nmed_AMB + nmed_T] = cimag(medidas_tensao[j].zmed);
  }
};

void estimadorAMB(GRAFO *grafo, long int numeroRamos, long int numeroBarras,
                  DMED *medidas, long int **numeroMedidas,
                  ALIMENTADOR *alimentadores, long int numeroAlimentadores,
                  DRAM *ramos, double Sbase, DBAR *barra) {
  long int nmed, nvar;
  int i, j, k, r;
  double *z = NULL, **h = NULL, ***H = NULL, **W = NULL, *x = NULL,
         *regua = NULL, aux = 0;

  double tol = 0.000001;
  //__complex__ double *x_bc = NULL;

  long int nmed_AMB;
  long int nmed_T;

  int *caminho = NULL;
  int *nadj_proxbarra = NULL;
  // caminho = malloc(numeroBarras * sizeof(int));
  // nadj_proxbarra = malloc(numeroBarras * sizeof(int));
  // busca_loop_grafo(grafo, numeroRamos, numeroBarras, caminho,
  // nadj_proxbarra);

  printf("Estimador de Estado AMB...\n");
  //--------------------------------------------------------------------------
  // Alocação de memória das variáveis do estimador de estado

  nmed = 0;
  nvar = 0;
  for (i = 0; i < 9; i++) {
    for (j = 0; j < 8; j++) {
      nmed = nmed + numeroMedidas[i][j];
    }
  }
  printf("numero barras: %d\n", numeroBarras);
  for (i = 0; i < numeroBarras; i++) {
    switch (grafo[i].fases) {
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
  if ((z = (double *)malloc((nmed) * sizeof(double))) == NULL) {
    printf(
        "Erro -- Nao foi possivel alocar espaco de memoria para o vetor z!!!!");
    exit(1);
  }
  if ((h = malloc((nmed) * sizeof(double *))) == NULL) {
    printf(
        "Erro -- Nao foi possivel alocar espaco de memoria para o vetor h!!!!");
    exit(1);
  }
  if ((x = (double *)malloc((nvar) * sizeof(double))) == NULL) {
    printf(
        "Erro -- Nao foi possivel alocar espaco de memoria para o vetor x!!!!");
    exit(1);
  }
  if ((regua = (double *)malloc((nvar) * sizeof(double))) == NULL) {
    printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor "
           "regua!!!!");
    exit(1);
  }

  aux = 0;

  int *RNP;
  // RNP = aloca_vetor(numeroBarras);
  k = 0;
  FILABARRAS *barraAtual;
  RNP = aloca_vetor_int(numeroBarras);
  barraAtual = &(alimentadores[0].rnp[0]);
  while (barraAtual != NULL) {
    RNP[k] = barraAtual->idNo;
    k++;
    barraAtual = barraAtual->prox;
  }

  nmed_AMB = conta_medidas_BC(medidas, nmed);

  DMED_COMPLEX *medidas_complexas = NULL;
  medidas_complexas = (DMED_COMPLEX *)malloc((nmed_AMB) * sizeof(DMED_COMPLEX));

  double *x_bc = NULL;
  x_bc = aloca_vetor(6 * numeroBarras);

  double *delta_x_bc = NULL;
  delta_x_bc = aloca_vetor(6 * numeroBarras);
  // x_bc = (double *)malloc((6*numeroRamos) * sizeof(double));

  DMED_COMPLEX *medidas_equivalentes = NULL;
  medidas_equivalentes =
      (DMED_COMPLEX *)malloc((nmed_AMB) * sizeof(DMED_COMPLEX));

  double *z_AMB = NULL;
  z_AMB = aloca_vetor(2 * nmed_AMB + 2 * nmed_T);

  DMED_COMPLEX *medidas_tensao = NULL;
  nmed_T = conta_medidas_Tensao(medidas, nmed);
  medidas_tensao = (DMED_COMPLEX *)malloc((nmed_T) * sizeof(DMED_COMPLEX));

  double *regua_x = NULL;
  double *regua_med = NULL;
  double *regua_med_inv = NULL;
  double **H_AMB = NULL;
  double **H_T = NULL;

  regua_med = aloca_vetor(3 * nmed_AMB);
  regua_med_inv = aloca_vetor(3 * nmed_AMB);
  regua_x = aloca_vetor(3 * numeroBarras);
  H_AMB = aloca_matriz(2 * nmed_AMB, 6 * numeroBarras);
  H_T = aloca_matriz(2 * nmed_T, 6 * numeroBarras);

  incializa_tensoes_grafo(grafo, numeroBarras, alimentadores,
                          numeroAlimentadores);
  x_bc = monta_vetor_x_inicial_AMB(numeroBarras, grafo);
  // printf("1\n");
  medidas_complexas = converte_medidas_para_complexo_AMB(medidas, nmed);
  //
  medidas_tensao =
      calcula_medida_tensao_complexa_AMB(medidas, nmed, grafo, numeroBarras);

  medidas_equivalentes = divide_medidas_por_tensao(medidas_complexas, nmed_AMB,
                                                   numeroBarras, grafo);
  // regua vetor x
  // monta_regua_x(numeroRamos, regua_x, ramos);
  regua_x = monta_regua_x_AMB(grafo, numeroBarras);
  // printf("4\n");

  double *hx_I = NULL;
  hx_I = (double *)malloc(2 * nmed_AMB * sizeof(double));
  // vetor h(x) das medidas de tensão
  double *hx_V = NULL;
  hx_V = (double *)malloc(2 * nmed_T * sizeof(double));

  for (int k = numeroBarras - 1; k >= 0; k--) {
    backward(&grafo[RNP[k]], grafo);
  }

  hx_V = monta_hx_V(nmed_T, medidas_tensao, grafo, numeroBarras);
  // montar H das medidas de tensão

  H_AMB = monta_matriz_H_AMB(numeroBarras, numeroRamos, nmed_AMB,
                             medidas_equivalentes, regua_x, ramos, grafo);

  // Medidas de tensão
  H_T = monta_matriz_H_AMB_Tensao(medidas_tensao, grafo, numeroBarras, nmed_T);
  // int i,j;

  // print elementos nao nulos de H_AMB
  // for(i = 0;i<2*nmed_AMB;i++){
  //     for (j = 0; j < 6*numeroBarras; j++){
  //         if (H_AMB[i][j]!=0.0){
  //             printf("%d, %d, %f\n", i, j, H_AMB[i][j]);
  //         }
  //
  //     }
  //  }

  // vetor h(x) para medidas de potencia
  // hx_I = monta_hx_I(nmed_AMB, medidas_equivalentes, grafo);

  hx_I =
      monta_hx_I_v2(grafo, medidas_equivalentes, numeroBarras, nmed_AMB, ramos);

  // monta matriz de covariancia
  W = monta_W_AMB(grafo, medidas_equivalentes, nmed_AMB, medidas_tensao,
                  nmed_T);

  // int i;

  int it = 0;
  int conv = 0;
  while (conv < 1) {
    // monta_z_complexa(medidas_equivalentes, z_eq, nmed_BC);
    // monta_z_real_e_imag(medidas_equivalentes, z_AMB, nmed_AMB,
    // medidas_tensao, nmed_T);
    monta_z_real_imag_AMB(medidas_equivalentes, z_AMB, nmed_AMB, medidas_tensao,
                          nmed_T);
    int st = 0;
    delta_x_bc = resolve_linear_QR_AMB(H_AMB, H_T, z_AMB, numeroBarras,
                                       nmed_AMB, nmed_T, hx_I, hx_V, W);
    atualiza_vetor_x(x_bc, delta_x_bc, numeroBarras);

    double nfx;
    nfx = norma_inf(delta_x_bc, 6 * numeroBarras);
    printf("\n\nIteracao:  %d \t|Dx|_inf =  %.17lf \t  \n", it, nfx);
    atualiza_estado_AMB(grafo, numeroBarras, x_bc);

    for (k = numeroBarras - 1; k >= 0; k--) {
      backward(&grafo[RNP[k]], grafo);
    }
    if (nfx<tol | it> 20) {
      conv = 10;
    }
    medidas_tensao =
        calcula_medida_tensao_complexa_AMB(medidas, nmed, grafo, numeroBarras);

//    H_T =
//        monta_matriz_H_AMB_Tensao(medidas_tensao, grafo, numeroBarras, nmed_T);

    medidas_equivalentes = divide_medidas_por_tensao(
        medidas_complexas, nmed_AMB, numeroBarras, grafo);
    hx_V = monta_hx_V(nmed_T, medidas_tensao, grafo, numeroBarras);
    // hx_I = monta_hx_I(nmed_AMB, medidas_equivalentes, grafo);
    hx_I = monta_hx_I_v2(grafo, medidas_equivalentes, numeroBarras, nmed_AMB,
                         ramos);
    it++;
  }
}

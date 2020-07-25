/* 
 * File:   funcoesBayesFusion.h
 * Author: Julio Massignan
 *
 * Created on 25 de março de 2020, 19:25
 */

#ifndef FUNCOESBAYESFUSION_H
#define	FUNCOESBAYESFUSION_H

#include "/home/julio/SuiteSparse-5.6.0/include/cholmod.h"
#include "/home/julio/SuiteSparse-5.6.0/include/SuiteSparseQR_C.h"
#include "data_structures.h"

// construtor para a matriz Jacobiana esparsa
void atualiza_H_ss(GRAFO *grafo, long int numeroBarras, DRAM *ramos, DMED *medidas, long int numeroMedidas, cholmod_triplet *T_SS);

// Exporta arquivos de Texto com as matrizes de covariância a priori (formato SuiteSparse) e vetor x a priori
void exportaPrioriQR_ss(GRAFO *grafo,double *regua, long int nvar, long int nmed, long int polar, cholmod_sparse *A_SS, cholmod_common *c);

// Importa arquivos de Texto com as matrizes de covariância a priori (formato SuiteSparse) e vetor x a priori
void importaPrioriQR_ss(double *x0, cholmod_sparse **R_SS,cholmod_common *c);

// Otimizador via QR esparso
int otimiza_Gauss_NewtonQR_ss(double *z, double **h, GRAFO *grafo, long int numeroBarras, DRAM *ramos, DMED *medidas, long int nvar, long int nmed, double *regua_comp, double *ponto, double tol);

// Função principal do estimador via fusão baeysiana
void estimadorBayesFusion_MAP(GRAFO *grafo, long int numeroBarras, DMED *medidas, long int **numeroMedidas, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos,double Sbase);


#endif	/* FUNCOESBAYESFUSION_H */
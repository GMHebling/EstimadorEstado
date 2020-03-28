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

// construtor para a matriz Jacobiana
void construtor_H(GRAFO *grafo, long int numeroBarras, DRAM *ramos, DMED medida);

// Otimizador via QR esparso
int otimiza_Gauss_NewtonQR_ss(double *z, double **h, GRAFO *grafo, long int numeroBarras, DRAM *ramos, DMED *medidas, long int nvar, long int nmed, double *regua_comp, double *ponto, double tol);

// Função principal do estimador via fusão baeysiana
void estimadorBayesFusion_MAP(GRAFO *grafo, long int numeroBarras, DMED *medidas, long int **numeroMedidas, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos,double Sbase);


#endif	/* FUNCOESBAYESFUSION_H */
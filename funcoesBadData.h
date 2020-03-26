/* 
 * File:   funcoesBadData.h
 * Author: Julio
 *
 * Created on 29 de Janeiro de 2019, 19:25
 */

#ifndef FUNCOESBADDATA_H
#define	FUNCOESBADDATA_H

#include "/home/julio/SuiteSparse-5.6.0/include/cholmod.h"
#include "/home/julio/SuiteSparse-5.6.0/include/SuiteSparseQR_C.h"

void residuosNormalizados(double *rN, double *bHat, long int m, long int n, double *Dz, double **H, double **W, double *z);
void residuosNormalizadosQR(double *rN, double *bHat, long int m, long int n, double *Dz, double **H, double **W, double *z);

void residuosNormalizadosQR_ss(double *rN, double *bHat, long int m, long int n, double *Dz, double ***H, DMED *medidas, cholmod_sparse *A_SS);

#endif	/* FUNCOESBADDATA_H */


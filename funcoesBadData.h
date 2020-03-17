/* 
 * File:   funcoesBadData.h
 * Author: Julio
 *
 * Created on 29 de Janeiro de 2019, 19:25
 */

#ifndef FUNCOESBADDATA_H
#define	FUNCOESBADDATA_H

void residuosNormalizados(double *rN, double *bHat, long int m, long int n, double *Dz, double **H, double **W, double *z);
void residuosNormalizadosQR(double *rN, double *bHat, long int m, long int n, double *Dz, double **H, double **W, double *z);

#endif	/* FUNCOESBADDATA_H */


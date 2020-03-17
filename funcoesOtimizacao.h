/* 
 * File:   funcoesOtimizacao.h
 * Author: Julio Massignan
 *
 * Created on 14 de Março de 2018, 12:44
 */

#ifndef FUNCOESOTIMIZACAO_H
#define	FUNCOESOTIMIZACAO_H

void montaGb(double **H,double **W, double **h, double *z,int nvar, int smed, double **G,double *b, int **Gsimb);
void salva_sol(FILE *arqout, GRAFO *grafo, long int numeroBarras, double **h, double ***H, double **G, double *z, double *x, double *Dx,  double *regua, int it, long int smed, long int nvar, double fx, double nGx, double tempo);


int otimiza_Gauss_Newton(double *z, double **h, double ***H, double **W, GRAFO *grafo, long int numeroBarras, DRAM *ramos, DMED *medidas, long int nvar, long int nmed, double *regua_comp, double *ponto, double tol, long int ref1, long int ref2, int **Gsimb);
int otimiza_Gauss_Newton_Hachtel(double *z, double **h, double ***H, double **W, GRAFO *grafo, long int numeroBarras, DRAM *ramos, DMED *medidas, long int nvar, long int nmed, double *regua, double *ponto, double tol, long int ref1, long int ref2);
int otimiza_Gauss_Newton_sparseHachtel_Virtuais(double *z, double **h, double **c, double **W, GRAFO *grafo, long int numeroBarras, DRAM *ramos, DMED *medidas,DMED *virtuais, long int nvar, long int nmed,long int nvir, double *regua, double *ponto, double tol, long int ref1, long int ref2);

#endif	/* FUNCOESOTIMIZACAO_H */


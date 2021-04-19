/* 
 * File:   funcoesMatematicas.h
 * Author: Julio Massignan
 *
 * Created on 19 de Julho de 2017, 17:20
 */

#ifndef FUNCOESMATEMATICAS_H
#define	FUNCOESMATEMATICAS_H


//Funções com matrizes complexas
__complex__ double *c_vetAloca(int n);
__complex__ double **c_matAloca(int n);
void c_matInversaZ(__complex__ double **A, int n);
void c_matConj(__complex__ double **A, int n);
void c_matIgual(__complex__ double **A, __complex__ double **B, int n);
void c_matTransp(__complex__ double **A, int n);
void c_matMultEsc(__complex__ double **A, __complex__ double b, int n);
void c_multMatMat(__complex__ double **A, __complex__ double **B, int n);
void c_multMatVet(__complex__ double **A, __complex__ double *B, int n);

void c_matImprime(__complex__ double **A, int n);


//Funções com matrizes reais
double **aloca_matriz(int m, int n);
double *aloca_vetor(int m);
double prod_escalar(double *a,double *b,int n);
double norma_inf(double *a,int n);
double norma_euc(double *a,int n);
void matTransp( double **A, int m, int n, double **At);

void mat_ig(double ***A,int m,int n, double **B);
void tira_refs(double ***A,int m,int n,int col1, int col2, double **temp, double *regua, double *x, long int it);

long int tira_refs_sparse(DMED *medidas, double ***A,int m,int n,int col1, int col2,long int *sparse_i, long int *sparse_j, double *sparse_x, double *regua, double *x, long int it);
long int mat_ig_sparse(double ***A,int m,int n, long int *i_sparse, long int *j_sparse, double  *x_sparse);

void tira_refs_regua(int m, int col1, int col2, double *regua);

void cat_hor(double **A,int m1,int n1,double **B,int m2,int n2, double **temp);
void cat_vert(double **A,int m1,int n1,double **B,int m2,int n2, double **temp);
void cat_vert_vet(double *A,int m1,double *B,int m2, double *temp);


//Funções de Algebra Linear
double eigenvalue_largest(double **A, long int m, long int n);

void QRfactorization(double **A,int m,int n, double **R);
void backwardSubs (double **A, long int m, long int n, double *b);
void forwardSubs (double **A, long int m, long int n, double *b);

double *solve_Householder(double **A,int m,int n, double *b);
double *solve_Householder_LS(double **A,int m,int n, double *b);
double *solve_PreconHouseholder(double **A,int m,int n, double *b);
double *solve_Crout(double **A,int m,int n, double *b);

long int *aloca_vetor_int(int m);

#endif	/* FUNCOESMATEMATICAS_H */


/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Observabilidade.h
 * Author: vitor
 *
 * Created on 20 de junho de 2020, 09:48
 */

#ifndef OBSERVABILIDADE_H
#define OBSERVABILIDADE_H
int observabilidade_H(double ***H,int nvar, int nmed,double* regua,GRAFO *grafo,long int numeroBarras, DMED *medidas, long int **numeroMedidas,DRAM *ramos,int nref);
int observabilidade_Hit(double ***H,int nvar, int nmed,double* regua,int nref,int it);
int *fatora_observabilidade(double **Ht, double **F,int nvar, int nmed, int npseud,int nref,int it);
void fimp_mat(FILE *arquivo, double **A,int m, int n);
void fimp_mate(FILE *arquivo, double **A,int m, int n);
void fimp_vet(FILE *arquivo, double *A,int m);
void perm_mat(double **A,int m, int n, int col1, int col2);
void transp(double **A,int m,int n, double **B);
void mat_igf(double ***A, int m, int n, float **B);
void printH(double ***H,int nvar, int nmed,double* regua,int nref,int it,char *buf);
int *fatora_observabilidadef(float **Ht, double **F,int nvar, int nmed, int npseud,int nref,int it);
float **aloca_matrizf(int m, int n);
void fimp_matef(FILE *arquivo, float **A,int m, int n);
void perm_matf(float **A,int m, int n, int col1, int col2);
void transpf(float **A,int m,int n, float **B);


#endif /* OBSERVABILIDADE_H */


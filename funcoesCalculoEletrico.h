/* 
 * File:   Leitura_Dados.h
 * Author: Julio Massignan
 *
 * Created on 24 de Fevereiro de 2017, 14:03
 */

#ifndef funcoesCalculoEletrico_H
#define	funcoesCalculoEletrico_H

//------------------------------------------------------------------------------
//Funções das grandezas elétricas
void Skm(GRAFO *noP, GRAFO *noS, DRAM *ramo, __complex__ double *S);
void Smk(GRAFO *noP, GRAFO *noS, DRAM *ramo, __complex__ double *S);
void Sk(GRAFO *grafo, long int k, __complex__ double *S);
void Ikm(GRAFO *noP, GRAFO *noS, DRAM *ramo, __complex__ double *Ikm);
void Imk(GRAFO *noP, GRAFO *noS, DRAM *ramo, __complex__ double *Imk);

//------------------------------------------------------------------------------
//Funções das derivadas
void dSkm(GRAFO *noP, GRAFO *noS, DRAM *ramo, __complex__ double *dS, long int opt, long int i);
void dSmk(GRAFO *noP, GRAFO *noS, DRAM *ramo, __complex__ double *dS, long int opt, long int i);
void dSk(GRAFO *grafo, long int k, __complex__ double *dS, long int opt, long int barra, long int fase);

void dSkm_ret(GRAFO *noP, GRAFO *noS, DRAM *ramo, __complex__ double *dS, long int opt, long int i);
void dSmk_ret(GRAFO *noP, GRAFO *noS, DRAM *ramo, __complex__ double *dS, long int opt, long int i);
void dSk_ret(GRAFO *grafo, long int k, __complex__ double *dS, long int opt, long int barra, long int fase);


#endif	/* funcoesCalculoEletrico_H */
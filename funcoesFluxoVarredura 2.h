/* 
 * File:   funcoesFluxoVarredura.h
 * Author: Julio Massignan
 *
 * Created on 11 de Dezembro de 2019, 18:29
 */

#ifndef FUNCOESFLUXOVARREDURA_H
#define	FUNCOESFLUXOVARREDURA_H


void incializaTensoesVarredura(GRAFO *grafo, long int numeroBarras, ALIMENTADOR *alimentadores, long int numeroAlimentadores);
    
int fluxoPotencia_BFS_Alimentador(GRAFO *grafo, long int numeroBarras, ALIMENTADOR alimentador, DRAM *ramos,double Sbase);
void fluxoPotencia_BFS_Multiplos(GRAFO *grafo, long int numeroBarras, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos,double Sbase);

void atualizaTapRegulador(DRAM *ramo);

void imprimeTensoesNodais(GRAFO *grafo);
void imprimeTaps(GRAFO *grafo);
void imprimeiInjecoesCorrentes(GRAFO *grafo);
void imprimeCorrentes(NOADJACENTE *noAdj);

void incializaTensoesRaiz(GRAFO *grafo, long int numeroBarras, ALIMENTADOR *alimentadores, long int numeroAlimentadores);
void incializaTensoesVarredura(GRAFO *grafo, long int numeroBarras, ALIMENTADOR *alimentadores, long int numeroAlimentadores);

void atualizaTapRegulador(DRAM *ramo);
void atualizaInjecoes(GRAFO *no);

void calculaCorrenteMontante();
void calculaTensaoJusante(complex double *Vp, complex double *Vs, complex double *Ips, DRAM *ramo);

void backward(GRAFO *noP, GRAFO *grafo);
BOOL controleReguladorTensao_LDC(double Vbase, double Ibase, complex double *Vp, complex double *Vs, complex double *Ips, complex double *Isp, DRAM *ramo);

BOOL forward_sweep(GRAFO *noP, GRAFO *grafo);

#endif	/* FUNCOESFLUXOVARREDURA_H */


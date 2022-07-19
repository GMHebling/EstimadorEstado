/* 
 * File:   Leitura_Dados.h
 * Author: Julio Massignan
 *
 * Created on 24 de Fevereiro de 2017, 14:03
 */

#ifndef funcoesWLS_H
#define	funcoesWLS_H

void exportaPrioriQR(GRAFO *grafo,double *regua, DMED *medidas, long int nvar, long int nmed, long int polar);

void monta_z(double *z, long int nmed, DMED *medida);
void monta_W(double **W, long int nmed, DMED *medidas);
void monta_W_Ident(double **W, long int nmed, DMED *medidas);
void incializa_vetor_x(GRAFO *grafo, long int numeroBarras, ALIMENTADOR *alimentadores, long int numeroAlimentadores, double *x, double *regua, long int nVariaveis);
void incializa_vetor_x_leitura(GRAFO *grafo, long int numeroBarras, ALIMENTADOR *alimentadores, long int numeroAlimentadores, double *x, double *regua, long int nVariaveis);
void incializa_vetor_x_perturbado(GRAFO *grafo, long int numeroBarras, ALIMENTADOR *alimentadores, long int numeroAlimentadores, double *x, double *regua, long int nVariaveis);
void incializa_vetor_x_FP(GRAFO *grafo,long int numeroBarras, ALIMENTADOR *alimentadores, long int numeroAlimentadores, double *x, double *regua, long int nVariaveis);
void saidaEstado(GRAFO *grafo, long int numeroBarras, int it, double tempoIt, double nFx, double nGx);
void atualiza_Rede(GRAFO *grafo, long int numeroBarras);
void atualiza_Modelo(GRAFO *grafo, long int numeroBarras, long int nmed, DMED *medidas);
void atualiza_h(GRAFO *grafo, long int numeroBarras, long int nmed, DMED *medidas);
void atualiza_H(GRAFO *grafo, long int numeroBarras, DRAM *ramos, DMED *medidas, long int numeroMedidas);
void atualiza_H_ret(GRAFO *grafo, long int numeroBarras, DRAM *ramos, DMED *medidas, long int numeroMedidas);
void atualiza_estado(GRAFO *grafo, double *x, double *regua, long int nVariaveis);
void exportaEstadoIT(GRAFO *grafo,int it, double *regua, long int nvar,int bus);
void tratamento_referencia(long int *ref_1, long int *ref_2, ALIMENTADOR *alimentador, double *regua, long int nVariaveis);
void tratamento_referencia_improv(long int *ref_1, long int *ref_2, ALIMENTADOR *alimentador, double *regua, long int nvar,double ***H,GRAFO *grafo, long int numeroBarras,DRAM *ramos, DMED *medidas, long int *nmed,int* pseu);
void add_pseudo(ALIMENTADOR *alimentador,GRAFO *grafo,DRAM *ramos, DMED *medidas, long int nmed);
void exportaRamoscorrentes(GRAFO *grafo, long int numeroBarras, double Sbase);
void estimadorWLS(GRAFO *grafo, long int numeroBarras, DMED *medidas, long int **numeroMedidas, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos, double Sbase,char **argv);
//void estimadorWLShachtel(GRAFO *grafo, long int numeroBarras, DMED *medidas, long int **numeroMedidas,DMED *virtuais, long int **numeroVirtuais, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos);

#endif	/* funcoesWLS_H */


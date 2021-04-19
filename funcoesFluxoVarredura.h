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
#endif	/* FUNCOESFLUXOVARREDURA_H */


/* 
 * File:   funcoesTopologia.h
 * Author: Julio Massignan
 *
 * Created on 3 de Março de 2017, 10:51
 */

#ifndef FUNCOESTOPOLOGIA_H
#define	FUNCOESTOPOLOGIA_H

void adicionaNo(FILABARRAS **setor, long int idNo);
void adicionaNoNaFila(FILABARRAS ** fila, long int idNo);
void apontaProxNoNaFila(FILABARRAS ** fila);
int retiraNoDaFila(FILABARRAS ** fila);
BOOL filaNaoVazia(FILABARRAS * fila);
BOOL ramoLigado(NOADJACENTE adjacente);
BOOL estaLista(FILABARRAS *setor, int idNo);


void buscaLargura(GRAFO * grafo, ALIMENTADOR *alimentador, long int idAlim, long int idNoRaiz, BOOL * visitado);
void buscaProfundidade(FILABARRAS *barraAtual, long int idNo, int profundidade,  BOOL *visitado, GRAFO * grafo, long int idAlim);
void buscaProfundidadeAlimentadores(GRAFO *grafo, long int numeroBarras, ALIMENTADOR **alimentadores, long int numeroAlimentadores);

void calculaPU(GRAFO *grafo, long int numeroBarras, DRAM *ramos, long int numeroRamos, double Sbase);
void tranformaPU(GRAFO *grafo, long int numeroBarras, double Sbase, ALIMENTADOR *alimentadores, long int numeroAlimentadores);
void atualizaTaps(DRAM *ramos, long int numeroRamos);

#endif	/* FUNCOESTOPOLOGIA_H */


/* 
 * File:   Leitura_Dados.h
 * Author: Julio Massignan
 *
 * Created on 24 de Fevereiro de 2017, 14:03
 */

#ifndef funcoesLeitura_H
#define	funcoesLeitura_H

//------------------------------------------------------------------------------
// Auxiliares e Leitura de csv
char* replace(char *st);
char* getfield(char* line, int num);

char* charLigacao(LIGACAO num);


//------------------------------------------------------------------------------
// Leitura de dados da rede elétrica
char *leituraDados(DBAR **barra, DRAM **ramo, long int *numeroBarras, long int *numeroRamos, long int *numeroAlimentadores);
void leituraDBAR(FILE *arquivo, DBAR **barras, long int *numeroBarras, long int *numeroAlimentadores);
void leituraDSHNT(FILE *arquivo, DBAR **barras, long int *numeroBarras);
void leituraDGD(FILE *arquivo, DBAR **barras, long int *numeroBarras);
void leituraDLIN(FILE *arquivo, DRAM **ramos, long int *numeroRamos, DBAR **barras, long int *numeroBarras);
void leituraDTRF(FILE *arquivo, DRAM **ramos, long int *numeroRamos, DBAR **barras, long int *numeroBarras);
void leituraDREG(FILE *arquivo, DRAM **ramos, long int *numeroRamos, DBAR **barras, long int *numeroBarras);
void leituraDSWTC(FILE *arquivo, DRAM **ramos, long int *numeroRamos, DBAR **barras, long int *numeroBarras);
long int **leituraMedidas(char *folder,char *file, DMED **medidas, DRAM *ramos, long int numeroRamos, DBAR *barras, long int numeroBarras, GRAFO *grafo, double Sbase);

void leituraVinicial(FILE *arquivo, DBAR **barras, long int *numeroBarras);

void geraGrafo(GRAFO ** grafo, DBAR *barras, long int numeroBarras,DRAM *ramos,long int numeroRamos);

//------------------------------------------------------------------------------
// Impressão de dados da rede elétrica
void salvaDadosRedeEletrica(DBAR *barras, long int numeroBarras, DRAM *ramos, long int numeroRamos, DMED *medidas, long int **numeroMedidas);
void salvaMedidasRedeEletrica(DMED *medidas, long int **numeroMedidas);

#endif	/* funcoesLeitura_H */


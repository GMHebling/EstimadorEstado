/* 
 * File:   funcoesBranchCurrent.h
 * Author: Gustavo Hebling
 *
 * Created on 15 de Abril de 2021, 18:20
 */


void estimadorBC_RECT(GRAFO *grafo, long int numeroRamos, long int numeroBarras, DMED *medidas, long int **numeroMedidas, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos,double Sbase, DBAR *barra);
void estimadorBC_RECT_Malhado(GRAFO *grafo, long int numeroRamos, long int numeroBarras, DMED *medidas, long int **numeroMedidas, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos,double Sbase, DBAR *barra);

int *montaRNP(ALIMENTADOR alimentadores);
void inicializa_vetor_estados_BC(double *x_bc, long int numeroRamos);

DMED_COMPLEX *converte_medidas_para_complexo(DMED *medidas, long int numeroMedidas);

void buscaProfundidadeLoop(GRAFO *grafo, int idNo, int barraAnterior, BOOL *visitado, int *caminho, int *contBarras, int *barraEntrada);
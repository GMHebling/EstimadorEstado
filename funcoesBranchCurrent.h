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

void buscaProfundidadeLoop(GRAFO *grafo, int idNo, int barraAnterior, BOOL *visitado, int *caminho, int *contBarras, int *barraEntrada, int *nadj_proxbarra);
void busca_loop_grafo(GRAFO *grafo, long int numeroRamos, long int numeroBarras, int *caminho, int *nadj_proxbarra);

int busca_id_caminho(int *caminho, int *numeroBarraMedida, long int numeroBarras);
DMED_COMPLEX *calcula_medida_tensao_complexa(DMED *medidas, long int numeroMedidas);
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

int conta_medidas_BC(DMED *medidas, long int numeroMedidas);
DMED_COMPLEX *divide_medidas_por_tensao(DMED_COMPLEX *medidas_complexas, long int numeroMedidas, long int numeroBarras, GRAFO *grafo);
void monta_regua_x(long int numeroRamos, double *regua_x, DRAM *ramos);
void monta_regua_medidas(long int nmed_BC, double *regua_med, double *regua_med_inv, DMED_COMPLEX *medidas_equivalentes);
double **monta_matriz_H(long int numeroRamos, long int nmed_BC, double *regua_x, double *regua_med, double *regua_med_inv);
int verifica_ramo_caminho(double *regua_caminho, double *regua_x, int idx_ramo, long int numeroRamos, int limite_caminho);
void monta_regua_caminho(long int numeroRamos, long int numeroBarras, double *regua_caminho, int *caminho, GRAFO *grafo);
int busca_id_caminho(int *caminho, int numeroBarraMedida, long int numeroBarras, GRAFO *grafo);
double calculo_theta_dv(GRAFO *grafo, int i_grafo, int i_fase, double *valor_hx_v);

double _dVdIr(double R, double X, double theta);
double _dVdIx(double R, double X, double theta);

double **monta_matriz_H_tensao(long int numeroBarras, long int numeroRamos, int nmed_T, int *caminho, DMED_COMPLEX *medidas_tensao, double *regua_x, double *regua_caminho, DRAM *ramos, GRAFO *grafo, double *hx_V);

double *resolve_linear_QR(double **H_BC, double *z, long int numeroRamos, long int nmed_BC);
double *resolve_linear_QR_Tensao(double **H_BC, double **H_T, double *z, long int numeroRamos, long int nmed_BC, int nmed_T, double *hx_I, double *hx_V);
void monta_z_complexa(DMED_COMPLEX *medidas_eq, __complex__ double *z, long int nmed_BC);
void monta_z_real_e_imag(DMED_COMPLEX *medidas_eq, double *z, long int nmed_BC, DMED_COMPLEX *medidas_tensao, int nmed_T);

void incializa_tensoes_grafo(GRAFO *grafo, long int numeroBarras, ALIMENTADOR *alimentadores, long int numeroAlimentadores);

void atualiza_Rede_BC(GRAFO *grafo, long int numeroBarras, DBAR *barra, double *regua_x, long int numeroRamos, double *x_bc);
void atualiza_Rede_BC_Tensao(GRAFO *grafo, long int numeroBarras, DBAR *barra, double *regua_x, long int numeroRamos, double *x_bc);

void exportaEstado_BC(GRAFO *grafo, long int numeroBarras);

int conta_medidas_Tensao(DMED *medidas, long int numeroMedidas);
int calcula_idx_ID(GRAFO *grafo, int ID, int numeroBarras);

void monta_regua_medidas_tensao(DMED_COMPLEX *medidas_tensao, int nmed_T, double *regua_V);
void calcula_hx_corrente(double **H_BC, double *x, double *hx, int nmed_BC, int numeroRamos);
void atualiza_vetor_x(double *x, double *dx, long int numeroRamos);

void buscaProfundidadeLoop(GRAFO *grafo, int idNo, int barraAnterior, BOOL *visitado, int *caminho, int *contBarras, int *barraEntrada, int *nadj_proxbarra);
void busca_loop_grafo(GRAFO *grafo, long int numeroRamos, long int numeroBarras, int *caminho, int *nadj_proxbarra);


DMED_COMPLEX *calcula_medida_tensao_complexa(DMED *medidas, long int numeroMedidas);
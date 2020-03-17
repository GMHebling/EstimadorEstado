/* 
 * File:   data_structures.h
 * Author: Julio Massignan
 *
 * Created on 24 de Fevereiro de 2017, 12:58
 */

#ifndef DATA_STRUCTURES_H
#define	DATA_STRUCTURES_H

#define PI 3.14159265359
#define EPS 1E-9

//------------------------------------------------------------------------------
//
//Tipos de dados especiais
//
//------------------------------------------------------------------------------
//Define variável BOOL
enum BOOLEAN {
	true = 1, /**< Valor verdadeiro para os testes booleanos. */
	false = 0 /**< Valor falso para os testes booleanos. */
};
typedef enum BOOLEAN BOOL;

// Tipo de ramo
typedef enum  {
	ramal, // trecho de circuito
	trafo, // transformador de potência
        regulador, //regulador de tensão
        chave //chave
} RAMO;

// Estado do ramo
typedef enum  {
	aberto, 
        fechado
} ESTADO;

// Tipo de ligação trifásica
typedef enum  {
        O,
        YN, // Estrela Aterrado
        D, //Delta
        Y, // Estrela Aberta
        OY, //Estrela aberto
        OD //Delta aberto
} LIGACAO;

// Tipo conexão das fases de ligação trifásica
typedef enum  {
        N,
        A, // Fase a
	B, // Fase b
        C, //Fase c
        AB, //Fase ab
        CA, //Fase ca
        BC, //Fase bc
        ABC //Fase abc
                
} FASES;

//------------------------------------------------------------------------------
//
//Estruturas de Dados para dados de equipamentos e da rede elétrica
//
//------------------------------------------------------------------------------
// Define dados de cargas
typedef struct {
    long int ID;
    FASES fases;
    LIGACAO lig;
    double Vbase;
    
    double Pnom[3];
    double Qnom[3];
    double ZIP;
    
} DLOAD;

// Define dados de shunts (bancos de capacitores fixos e automáticos)
typedef struct {
    long int ID;
    FASES fases;
    LIGACAO lig;
    double Vbase;
    
    double Qnom[3];
    long int controle;
    long int num;
    double DV;
    double Vset[3];   
} DSHNT;

// Define dados de geradores distribuidos
typedef struct {
    long int ID;  
    FASES fases;
    LIGACAO lig;
    double Vbase;
    double Snominal;
    
    double Pnom[3];
    double Qnom[3];
    long int controle;
    
    double Qmin;
    double Qmax;
    double Vset[3];
    long int controlePV;
} DGD;

// Define dados de barras
typedef struct {
    long int ID;
    long int i;
    FASES fases;
    LIGACAO ligacao;
    double Vbase;
    long int tipo; // Barra  0 - PQ; 1 - PV; 2 - Ref
    
    __complex__ double Vref[3];
    __complex__ double Vinicial[3]; //recebe de arquivo de entrada a tensão inicial
    
    DLOAD loads[10]; //cargas conectadas à barra - limite de 10 cargas por barra para não ter que trabalhar na alocação dinâmica dentro do array de struct
    long int nloads;
    DSHNT *shunts; //shunts conectadas à barra
    long int nshunts;
    DGD *gds; //geradores distribuídos conectadas à barra
    long int ngds;
    
} DBAR;

// Define os dados de Linhas de Transmissão (circuitos)
typedef struct {
    FASES fases;
    double comprimento;
    __complex__ double Zaa, Zab, Zac, Zbb, Zbc, Zcc;
    double Baa, Bab, Bac, Bbb, Bbc, Bcc;
  
} DLIN;

// Define os dados de transformadores de potência
typedef struct {
  FASES fases;
  double Vpri, Vsec, Snominal, R, X;
  LIGACAO lig_pri, lig_sec;
  long int defasamento;
  double tap_pri, tap_sec;
  
} DTRF;

// Define os dados de reguladores de tensão
typedef struct {
  FASES fases;
  double Vnom, regulacao, Snominal, R, X;
  long int ntaps;
  LIGACAO lig;
  
  //Parametros do controlador
  double TP, TC;
  double deltaV, R1,X1,R2,X2,R3,X3, V1,V2,V3;
  long int controle; // (0 = Somente Forward sem Restrição / 1 = Locked Forward / 2 = Locked Reverse / 3 = Bidirectional / 4 = Idle / 5 = Neutral Reverse / 6 = Cogenartion / Ver manual Siemens MJ4A)
  double tap[3];
  double deltaVr, R1r,X1r,R2r,X2r,R3r,X3r, V1r,V2r,V3r; //Parametros de controle reverso
  
} DREG;

// Define os dados dos ramos da rede elétrica
typedef struct {
    long int DE;
    long int PARA;
    long int k;
    long int m;
    FASES fases;
    RAMO tipo; // 0= linha  1= trafo  2= regulador  3= chave
    ESTADO estado;
    
    //Dados
    DLIN linha;
    DTRF trafo;
    DREG regulador;
    
    //Matrizes de Impedância
    __complex__ double **Ypp;
    __complex__ double **Yps;
    __complex__ double **Ysp;
    __complex__ double **Yss;
    double tap_pri[3];
    double tap_sec[3];
    
} DRAM;

// Define os dados dos medidores da rede elétrica
typedef struct {
    long int DE;
    long int PARA;
    long int id;
    long int k;
    long int m;
    long int ramo;
    FASES fases;
    long int tipo; //
    long int idAlim; //
    long int idArea; //
    long int idFront;
    ESTADO ligado; //aberto = on, 
    
    double zmed;
    double sigma;
    double prec;
    
    long int par; //para indicar o par PQ, par VTeta, e par IdeltaI
    
    double h;
    long int nvar;
    double *reguaH;
    double *H; //salva o gradiente da respectiva medida
    
} DMED;

//------------------------------------------------------------------------------
//
//Estruturas de Dados para dados topológicos e RNP
//
//------------------------------------------------------------------------------
typedef struct Fila {
    long int idNo;
    long int profundidade;
    struct Fila * prox;
} FILABARRAS;

typedef struct {
    long int noRaiz;
    long int idRaiz;
    long int idAlim;
    long int numeroNos;
    FILABARRAS rnp[1];
} ALIMENTADOR;

typedef struct {
  long int idNo;/**< Inteiro que identifica cada nó do SDR. */
  ESTADO estado; // 0 = desligado e 1 = ligado
  RAMO tipo;
  double relacao;
  DRAM *ramo;
  long int idram;
  
  //Medidores
  long int nmed;
  DMED **medidores; /**< Lista dos medidores. Limite de 20 para evitar alocação dinamica*/
  
  //Grandezas elétricas que caracterizam o ramo
  __complex__ double S[3];
  __complex__ double Cur[3];
  
} NOADJACENTE;

typedef struct {
  long int idNo;
  long int profundidade;
  long int idAlim;
  long int idArea;
  long int idFront;
  FASES fases;
  double Vbase;
  long int tipo; // Barra  0 - PQ; 1 - PV; 2 - Ref
  DBAR *barra; //ponteiro para dados completos de barras
  long int numeroAdjacentes; /**< Indica o tamanho da lista de adjacentes do nó. */
  NOADJACENTE adjacentes[15]; /**< Lista dos nós adjacentes. Limite de 15 para evitar alocação dinamica*/
  
  //Medidores
  long int nmed;
  DMED **medidores; /**< Lista dos medidores. Limite de 20 para evitar alocação dinamica*/
  
  //Impedância Shunt Totoal
  __complex__ Ysh[3][3];
  
  //Grandezas elétricas que caracterizam a barra
  __complex__ double S[3];
  __complex__ double V[3];
  __complex__ double Cur[3];
  __complex__ double V_aux[3];
  double V_aux_Sig[3];
  __complex__ double V_est[3];
  __complex__ double V_coord[3];
  int estimado;
} GRAFO;

typedef struct sparseMatrix {
    long int i, j;
    double *ij; //Ponteiro para o valor do elemento i j
    double x; //Valor do elemento i j
    struct sparseMatrix * prox;
} SPARSE;

typedef struct compressedColumnSPARSE {
    int n;
    double diag;
    SPARSE *val;
} CCSPARSE;

#endif	/* DATA_STRUCTURES_H */


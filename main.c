#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#include "data_structures.h"
#include "funcoesLeitura.h"
#include "funcoesTopologia.h"
#include "funcoesMatematicas.h"
#include "funcoesWLS.h"
#include "funcoesAMB.h"
#include "funcoesBranchCurrent.h"

int main(int argc, char **argv)
{
    long int numeroBarras = 0;
    long int numeroAlimentadores = 0, numeroAreas = 0;
    long int numeroRamos = 0;
    long int **numeroMedidas = NULL, **numeroMedidasPMU = NULL, **numeroVirtuais = NULL;
    //long int i = 0;
    double Sbase = 1000000; //VA
    char *folder = NULL;
    DBAR *barra = NULL;
    DRAM *ramo = NULL;
    GRAFO *grafo = NULL;
    DMED *medida = NULL, *medidaPMU = NULL, *virtuais = NULL;
    ALIMENTADOR *alimentador = NULL, *areas = NULL;


    // Leitura dos dados da rede elétrica
    folder = leituraDados(&barra, &ramo, &numeroBarras, &numeroRamos, &numeroAlimentadores);
    printf("leituraDados ok\n");
    //Melhorar o tratamento dos taps de reguladores e trafos - leitura separa de estados lógicos
    // Cria estrutura de dados da rede elétrica
    geraGrafo(&grafo, barra, numeroBarras, ramo, numeroRamos);
    printf("geraGrafo ok\n");
    // Cria as listas encadeadas radiais dos alimentadores
    buscaProfundidadeAlimentadores(grafo, numeroBarras, &alimentador, numeroAlimentadores);
    printf("buscaProfundidadeAlimentadores ok\n");
    // Transforma em pu e cria matrizes de admitância iniciais
    calculaPU(grafo, numeroBarras, ramo, numeroRamos, Sbase);
    printf("calculaPU ok\n");
    //Atualiza Valores de taps
    atualizaTaps(ramo, numeroRamos); //terminar o tap de trafos
    printf("atualizaTaps ok \n");
    // Leitura das Medidas e associa os medidores ao grafo da rede - Numero de medidas retorna matriz tipo de media/por fase
    // Leitura das Medidas e associa os medidores ao grafo da rede - Numero de medidas retorna matriz tipo de media/por fase
    numeroMedidas = leituraMedidas(folder, "DMED.csv", &medida, ramo, numeroRamos, barra, numeroBarras, grafo, Sbase); //Melhorar o tratamento de chaves

    
    //estimadorWLS(grafo, numeroBarras, medida, numeroMedidas, alimentador, numeroAlimentadores, ramo, Sbase / 1000);
    clock_t tIni = clock();

    //constroi o caminho do no raiz ate cada barra
    
    
    
    
    //estimadorBC_RECT(grafo, numeroRamos, numeroBarras, medida, numeroMedidas, alimentador, numeroAlimentadores, ramo, Sbase / 1000, barra);
    
    

    estimadorAMB(grafo, numeroRamos, numeroBarras, medida, numeroMedidas, alimentador, numeroAlimentadores, ramo, Sbase / 1000, barra);
    clock_t t1 = clock();
    double tempoWLS = (double)(t1 - tIni) / CLOCKS_PER_SEC;
    printf("\nEstimação BC: %lf", tempoWLS);

    exportaCasoReferencia(grafo, numeroBarras, Sbase);
    exportaEstado_BC(grafo, numeroBarras);
    imprimeEstado(grafo, numeroBarras);
    
    //    salvaDadosRedeEletrica(barra, numeroBarras, ramo, numeroRamos, medida, numeroMedidas);
    free(barra);
    free(ramo);
    free(grafo);
    free(alimentador);
    return (EXIT_SUCCESS);
}

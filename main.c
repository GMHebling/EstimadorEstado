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

/*
 * 
 */
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
    //printf("leituraDados ok\n");
    //Melhorar o tratamento dos taps de reguladores e trafos - leitura separa de estados lógicos
    // Cria estrutura de dados da rede elétrica
    geraGrafo(&grafo, barra, numeroBarras, ramo, numeroRamos);
    //printf("geraGrafo ok\n");
    // Cria as listas encadeadas radiais dos alimentadores
    buscaProfundidadeAlimentadores(grafo, numeroBarras, &alimentador, numeroAlimentadores);
    //printf("buscaProfundidadeAlimentadores ok\n");
    // Transforma em pu e cria matrizes de admitância iniciais
    calculaPU(grafo, numeroBarras, ramo, numeroRamos, Sbase);
    //printf("calculaPU ok\n");
    //Atualiza Valores de taps
    atualizaTaps(ramo, numeroRamos); //terminar o tap de trafos
    //printf("atualizaTaps ok \n");
    // Leitura das Medidas e associa os medidores ao grafo da rede - Numero de medidas retorna matriz tipo de media/por fase
    numeroMedidas = leituraMedidas(folder, "DMED.csv", &medida, ramo, numeroRamos, barra, numeroBarras, grafo, Sbase); //Melhorar o tratamento de chaves
    numeroVirtuais = leituraMedidas(folder,"DVMED.csv", &virtuais, ramo, numeroRamos, barra, numeroBarras,grafo,Sbase);  //Melhorar o tratamento de chaves
    
    //Função Atualiza Status Lógicos - chaves e taps no arquivo DSTATUS.csv
    //printf("leituraMedidas ok\n");
    // Observabilidade e Seleção de Pseudo-Medidas - Fatoração da matriz H -  Tratamento das referências
    // Observabilidade Hdelta

    // Fluxo de Potência via Newton Raphson Trifásico QR
    //fluxoPotencia_NRQR(grafo, numeroBarras, medida, numeroMedidas, alimentador, numeroAlimentadores,ramo,Sbase/1000);

    // Estimador WLS Convencional
     
    estimadorWLS(grafo, numeroBarras, medida, numeroMedidas, alimentador, numeroAlimentadores, ramo, Sbase / 1000);

    //    salvaDadosRedeEletrica(barra, numeroBarras, ramo, numeroRamos, medida, numeroMedidas);
    free(barra);
    free(ramo);
    free(grafo);
    free(alimentador);
    return (EXIT_SUCCESS);
}

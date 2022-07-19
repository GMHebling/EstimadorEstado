#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <complex.h>

#include "data_structures.h"
#include "funcoesWLS.h"
#include "leitura.h"
#include "matriz.h"
#include "funcoesTopologia.h"
#include "funcoesCalculoEletrico.h"
#include "funcoesOtimizacao.h"
#include "funcoesMatematicas.h"

//Imprime na tela o estado atual da rede (tensões complexas nodais em cada fase)
void imprimeEstado(GRAFO *grafo, long int numeroBarras)
{
    int i;

    printf("\nTensoes Nodais (p.u.):\n");
    for (i = 0; i < numeroBarras; i++)
    {
        //Retangulares
        //printf("%d\tVa: %.5lf + j%.5lf\tVb: %.5lf + j%.5lf\tVc: %.5lf + j%.5lf\n",grafo[i].barra->ID,__real__ grafo[i].V[0],__imag__ grafo[i].V[0],__real__ grafo[i].V[1],__imag__ grafo[i].V[1],__real__ grafo[i].V[2],__imag__ grafo[i].V[2]);
        //Polares
        switch (grafo[i].fases)
        {
        case 1:
            printf("%ld\tVa: %.5lf | %.3lf \tVb:    -    |    -   \tVc:    -    |    -   \n", grafo[i].barra->ID, cabs(grafo[i].V[0]), carg(grafo[i].V[0]) * 180 / PI);
            break;
        case 2:
            printf("%ld\tVa:    -    |    -    \tVb: %.5lf | %.3lf\tVc:    -    |    -   \n", grafo[i].barra->ID, cabs(grafo[i].V[1]), carg(grafo[i].V[1]) * 180 / PI);
            break;
        case 3:
            printf("%ld\tVa:    -    |    -    \tVb:    -    |    -   \tVc: %.5lf | %.3lf\n", grafo[i].barra->ID, cabs(grafo[i].V[0]), carg(grafo[i].V[2]) * 180 / PI);
            break;
        case 4:
            printf("%ld\tVa: %.5lf | %.3lf \tVb: %.5lf | %.3lf\tVc:    -    |    -   \n", grafo[i].barra->ID, cabs(grafo[i].V[0]), carg(grafo[i].V[0]) * 180 / PI, cabs(grafo[i].V[1]), carg(grafo[i].V[1]) * 180 / PI);
            break;
        case 5:
            printf("%ld\tVa: %.5lf | %.3lf \tVb:    -    |    -   \tVc: %.5lf | %.3lf\n", grafo[i].barra->ID, cabs(grafo[i].V[0]), carg(grafo[i].V[0]) * 180 / PI, cabs(grafo[i].V[2]), carg(grafo[i].V[2]) * 180 / PI);
            break;
        case 6:
            printf("%ld\tVa:    -    |    -    \tVb: %.5lf | %.3lf\tVc: %.5lf | %.3lf\n", grafo[i].barra->ID, cabs(grafo[i].V[1]), carg(grafo[i].V[1]) * 180 / PI, cabs(grafo[i].V[2]), carg(grafo[i].V[2]) * 180 / PI);
            break;
        case 7:
            printf("%ld\tVa: %.5lf | %.3lf \tVb: %.5lf | %.3lf\tVc: %.5lf | %.3lf\n", grafo[i].barra->ID, cabs(grafo[i].V[0]), carg(grafo[i].V[0]) * 180 / PI, cabs(grafo[i].V[1]), carg(grafo[i].V[1]) * 180 / PI, cabs(grafo[i].V[2]), carg(grafo[i].V[2]) * 180 / PI);
            break;
        }
    }
}

void imprimeInjecoes(GRAFO *grafo, long int numeroBarras)
{
    int i;

    printf("\nInjecoes Nodais (p.u.):\n");
    for (i = 0; i < numeroBarras; i++)
    {
        //Retangulares
        //printf("%d\tVa: %.5lf + j%.5lf\tVb: %.5lf + j%.5lf\tVc: %.5lf + j%.5lf\n",grafo[i].barra->ID,__real__ grafo[i].V[0],__imag__ grafo[i].V[0],__real__ grafo[i].V[1],__imag__ grafo[i].V[1],__real__ grafo[i].V[2],__imag__ grafo[i].V[2]);
        //Polares
        switch (grafo[i].fases)
        {
        case 1:
            printf("%ld\tSa: %.5lf | %.3lf \tSb:    -    |    -   \tSc:    -    |    -   \n", grafo[i].barra->ID, creal(grafo[i].S[0]), cimag(grafo[i].S[0])) ;
            break;
        case 2:
            printf("%ld\tSa:    -    |    -    \tSb: %.5lf | %.3lf\tSc:    -    |    -   \n", grafo[i].barra->ID, creal(grafo[i].S[1]), cimag(grafo[i].S[1]) );
            break;
        case 3:
            printf("%ld\tSa:    -    |    -    \tSb:    -    |    -   \tSc: %.5lf | %.3lf\n", grafo[i].barra->ID, creal(grafo[i].S[0]), cimag(grafo[i].S[2]) );
            break;
        case 4:
            printf("%ld\tSa: %.5lf | %.3lf \tSb: %.5lf | %.3lf\tSc:    -    |    -   \n", grafo[i].barra->ID, creal(grafo[i].S[0]), cimag(grafo[i].S[0]) , creal(grafo[i].S[1]), cimag(grafo[i].S[1]) );
            break;
        case 5:
            printf("%ld\tSa: %.5lf | %.3lf \tSb:    -    |    -   \tSc: %.5lf | %.3lf\n", grafo[i].barra->ID, creal(grafo[i].S[0]), cimag(grafo[i].S[0]) , creal(grafo[i].S[2]), cimag(grafo[i].S[2]) );
            break;
        case 6:
            printf("%ld\tSa:    -    |    -    \tSb: %.5lf | %.3lf\tSc: %.5lf | %.3lf\n", grafo[i].barra->ID, creal(grafo[i].S[1]), cimag(grafo[i].S[1]) , creal(grafo[i].S[2]), cimag(grafo[i].S[2]) );
            break;
        case 7:
            printf("%ld\tSa: %.5lf | %.3lf \tSb: %.5lf | %.3lf\tSc: %.5lf | %.3lf\n", grafo[i].barra->ID, creal(grafo[i].S[0]), cimag(grafo[i].S[0]) , creal(grafo[i].S[1]), cimag(grafo[i].S[1]) , creal(grafo[i].S[2]), cimag(grafo[i].S[2]) );
            break;
        }
    }
}

//SalSa arquivo com o estado atual da rede (vetor x) conforme régua
void exportaEstado(GRAFO *grafo, double *regua, long int nvar)
{
    int i, k, fase;
    FILE *arqout;

    arqout = fopen("state.txt", "w+");
    for (i = 0; i < nvar; i++)
    {
        if (regua[i] > 0)
        {
            k = (int)regua[i];
            fase = (int)cabs((regua[i] - k) * 10);
            fprintf(arqout, "%.2f\t%.15f\n", regua[i], cabs(grafo[k].V[fase]));
        }
        else
        {
            k = (int)regua[i];
            fase = (int)cabs((regua[i] - k) * 10);
            k = -k;
            fprintf(arqout, "%.2f\t%.15f\n", regua[i], carg(grafo[k].V[fase]));
        }
    }
    fclose(arqout);
}

void exportaRamoscorrentes(GRAFO *grafo, long int numeroBarras, double Sbase)
{
    int i, k, fase;
    FILE *arquivo;
    complex A;
    arquivo = fopen("correntes.txt", "w+");
    for (i = 0; i < numeroBarras; i++)
    {
        //Percorre os ramos adjacentes
        for (k = 0; k < grafo[i].numeroAdjacentes; k++)
        {

            switch (grafo[i].adjacentes[k].ramo->fases)
            {
            case 1:
                fprintf(arquivo, "1,0,%ld,%ld,1,%.10lf,%.10lf\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, cabs((grafo[i].adjacentes[k].S[0]*1000000)/(grafo[i].Vbase*grafo[i].V[0])) ,carg((grafo[i].adjacentes[k].S[0]*1000000/(grafo[i].Vbase*grafo[i].V[0]))));
                break;
            case 2:
                fprintf(arquivo, "1,0,%ld,%ld,2,%.10lf,%.10lf\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, cabs((grafo[i].adjacentes[k].S[1]*1000000)/(grafo[i].Vbase*grafo[i].V[1])) ,carg((grafo[i].adjacentes[k].S[1]*1000000/(grafo[i].Vbase*grafo[i].V[1]))));
                break;
            case 3:
                fprintf(arquivo, "1,0,%ld,%ld,3,%.10lf,%.10lf\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, cabs((grafo[i].adjacentes[k].S[2]*1000000)/(grafo[i].Vbase*grafo[i].V[2])) ,carg((grafo[i].adjacentes[k].S[2]*1000000/(grafo[i].Vbase*grafo[i].V[2]))));
                break;
            case 4:
                fprintf(arquivo, "1,0,%ld,%ld,1,%.10lf,%.10lf\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, cabs((grafo[i].adjacentes[k].S[0]*1000000)/(grafo[i].Vbase*grafo[i].V[0])) ,carg((grafo[i].adjacentes[k].S[0]*1000000/(grafo[i].Vbase*grafo[i].V[0]))));
                fprintf(arquivo, "1,0,%ld,%ld,2,%.10lf,%.10lf\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, cabs((grafo[i].adjacentes[k].S[1]*1000000)/(grafo[i].Vbase*grafo[i].V[1])) ,carg((grafo[i].adjacentes[k].S[1]*1000000/(grafo[i].Vbase*grafo[i].V[1]))));
                break;
            case 5:
                fprintf(arquivo, "1,0,%ld,%ld,1,%.10lf,%.10lf\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, cabs((grafo[i].adjacentes[k].S[0]*1000000)/(grafo[i].Vbase*grafo[i].V[0])) ,carg((grafo[i].adjacentes[k].S[0]*1000000/(grafo[i].Vbase*grafo[i].V[0]))));
                fprintf(arquivo, "1,0,%ld,%ld,3,%.10lf,%.10lf\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, cabs((grafo[i].adjacentes[k].S[2]*1000000)/(grafo[i].Vbase*grafo[i].V[2])) ,carg((grafo[i].adjacentes[k].S[2]*1000000/(grafo[i].Vbase*grafo[i].V[2]))));
                break;
                fprintf(arquivo, "1,0,%ld,%ld,2,%.10lf,%.10lf\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, cabs((grafo[i].adjacentes[k].S[1]*1000000)/(grafo[i].Vbase*grafo[i].V[1])) ,carg((grafo[i].adjacentes[k].S[1]*1000000/(grafo[i].Vbase*grafo[i].V[1]))));
                fprintf(arquivo, "1,0,%ld,%ld,3,%.10lf,%.10lf\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, cabs((grafo[i].adjacentes[k].S[2]*1000000)/(grafo[i].Vbase*grafo[i].V[2])) ,carg((grafo[i].adjacentes[k].S[2]*1000000/(grafo[i].Vbase*grafo[i].V[2]))));
                break;
            case 7:
                fprintf(arquivo, "1,0,%ld,%ld,1,%.10lf,%.10lf\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, cabs((grafo[i].adjacentes[k].S[0]*1000000)/(grafo[i].Vbase*grafo[i].V[0])) ,carg((grafo[i].adjacentes[k].S[0]*1000000/(grafo[i].Vbase*grafo[i].V[0]))));
                fprintf(arquivo, "1,0,%ld,%ld,2,%.10lf,%.10lf\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, cabs((grafo[i].adjacentes[k].S[1]*1000000)/(grafo[i].Vbase*grafo[i].V[1])) ,carg((grafo[i].adjacentes[k].S[1]*1000000/(grafo[i].Vbase*grafo[i].V[1]))));
                fprintf(arquivo, "1,0,%ld,%ld,3,%.10lf,%.10lf\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, cabs((grafo[i].adjacentes[k].S[2]*1000000)/(grafo[i].Vbase*grafo[i].V[2])) ,carg((grafo[i].adjacentes[k].S[2]*1000000/(grafo[i].Vbase*grafo[i].V[2]))));
                break;
            }
        }
    }
    fclose(arquivo);
}

void exportaEstadoIT(GRAFO *grafo,int it, double *regua, long int nvar,int bus)
{
    int i, k, fase;
    FILE *arqout;
    char buff[50];
    sprintf(buff,"stateiTbus%d.txt",bus);
    if (it==0)
    {
        arqout = fopen(buff, "w+");
    }
    else
    {
        arqout = fopen(buff, "a");
    }
    
    fprintf(arqout, "%d,%.15f,%.15f,%.15f,%.15f,%.15f,%.15f\n", it, cabs(grafo[bus].V[0]),cabs(grafo[bus].V[1]),cabs(grafo[bus].V[2]),carg(grafo[bus].V[0]),carg(grafo[bus].V[1]),carg(grafo[bus].V[2]),"a");

    fclose(arqout);
}
void saidaEstado(GRAFO *grafo, long int numeroBarras, int it, double tempoIt, double nFx, double nGx)
{
    int i, j, k;

    FILE *arquivo;
    arquivo = fopen("resultadoEstado.txt", "wt");

    fprintf(arquivo, "|Dx|_inf =  %.7lf \t |Grad|_inf =  %.7lf \n", nFx, nGx);
    fprintf(arquivo, "\n\n Convergência em %d iteracoes e tempo: %.4lf", it, tempoIt);

    fprintf(arquivo, "\nTensoes Nodais: Fase-Terra\n");
    for (i = 0; i < numeroBarras; i++)
    {
        fprintf(arquivo, "Vbase:%.3lf\t", grafo[i].Vbase);
        switch (grafo[i].fases)
        {
        case 1:
            fprintf(arquivo, "A\t%ld\tVa: %.5lf | %.3lf \tSb:    -    |    -   \tVc:    -    |    -   \n", grafo[i].barra->ID, cabs(grafo[i].V[0]), carg(grafo[i].V[0]) );
            break;
        case 2:
            fprintf(arquivo, "B\t%ld\tVa:    -    |    -    \tSb: %.5lf | %.3lf\tVc:    -    |    -   \n", grafo[i].barra->ID, cabs(grafo[i].V[1]), carg(grafo[i].V[1]) );
            break;
        case 3:
            fprintf(arquivo, "C\t%ld\tVa:    -    |    -    \tSb:    -    |    -   \tVc: %.5lf | %.3lf\n", grafo[i].barra->ID, cabs(grafo[i].V[2]), carg(grafo[i].V[2]) );
            break;
        case 4:
            fprintf(arquivo, "AB\t%ld\tVa: %.5lf | %.3lf \tSb: %.5lf | %.3lf\tVc:    -    |    -   \n", grafo[i].barra->ID, cabs(grafo[i].V[0]), carg(grafo[i].V[0]) , cabs(grafo[i].V[1]), carg(grafo[i].V[1]) );
            break;
        case 5:
            fprintf(arquivo, "CA\t%ld\tVa: %.5lf | %.3lf \tSb:    -    |    -   \tVc: %.5lf | %.3lf\n", grafo[i].barra->ID, cabs(grafo[i].V[0]), carg(grafo[i].V[0]) , cabs(grafo[i].V[2]), carg(grafo[i].V[2]) );
            break;
        case 6:
            fprintf(arquivo, "BC\t%ld\tVa:    -    |    -    \tSb: %.5lf | %.3lf\tVc: %.5lf | %.3lf\n", grafo[i].barra->ID, cabs(grafo[i].V[1]), carg(grafo[i].V[1]) , cabs(grafo[i].V[2]), carg(grafo[i].V[2]) );
            break;
        case 7:
            fprintf(arquivo, "ABC\t%ld\tVa: %.5lf | %.3lf \tSb: %.5lf | %.3lf\tVc: %.5lf | %.3lf\n", grafo[i].barra->ID, cabs(grafo[i].V[0]), carg(grafo[i].V[0]) , cabs(grafo[i].V[1]), carg(grafo[i].V[1]) , cabs(grafo[i].V[2]), carg(grafo[i].V[2]) );
            break;
        }
    }
    fprintf(arquivo, "\nTensoes Nodais: Fase-Fase\n");
    for (i = 0; i < numeroBarras; i++)
    {
        fprintf(arquivo, "Vbase:%.3lf\t", grafo[i].Vbase);
        switch (grafo[i].fases)
        {
        case 1:
            fprintf(arquivo, "A\t%ld\tVan: %.5lf | %.3lf \tVbn:    -    |    -   \tVcn:    -    |    -   \n", grafo[i].barra->ID, cabs(grafo[i].V[0]), carg(grafo[i].V[0]) );
            break;
        case 2:
            fprintf(arquivo, "B\t%ld\tVan:    -    |    -    \tVbn: %.5lf | %.3lf\tVcn:    -    |    -   \n", grafo[i].barra->ID, cabs(grafo[i].V[1]), carg(grafo[i].V[1]) );
            break;
        case 3:
            fprintf(arquivo, "C\t%ld\tVan:    -    |    -    \tVbn:    -    |    -   \tVcn: %.5lf | %.3lf\n", grafo[i].barra->ID, cabs(grafo[i].V[2]), carg(grafo[i].V[2]) );
            break;
        case 4:
            fprintf(arquivo, "AB\t%ld\tVab: %.5lf | %.3lf \tVbc:   -    |   -  \tVca:    -    |    -   \n", grafo[i].barra->ID, cabs(grafo[i].V[0] - grafo[i].V[1]), carg(grafo[i].V[0] - grafo[i].V[1]) );
            break;
        case 5:
            fprintf(arquivo, "CA\t%ld\tVab:   -   |   -  \tVbc:    -    |    -   \tVca: %.5lf | %.3lf\n", grafo[i].barra->ID, cabs(grafo[i].V[2] - grafo[i].V[0]), carg(grafo[i].V[2] - grafo[i].V[0]) );
            break;
        case 6:
            fprintf(arquivo, "BC\t%ld\tVab:    -    |    -    \tVbc: %.5lf | %.3lf\tVca:   -    |   -  \n", grafo[i].barra->ID, cabs(grafo[i].V[1] - grafo[i].V[2]), carg(grafo[i].V[1] - grafo[i].V[2]) );
            break;
        case 7:
            fprintf(arquivo, "ABC\t%ld\tVab: %.5lf | %.3lf \tVbc: %.5lf | %.3lf\tVca: %.5lf | %.3lf\n", grafo[i].barra->ID, cabs(grafo[i].V[0] - grafo[i].V[1]), carg(grafo[i].V[0] - grafo[i].V[1]) , cabs(grafo[i].V[1] - grafo[i].V[2]), carg(grafo[i].V[1] - grafo[i].V[2]) , cabs(grafo[i].V[2] - grafo[i].V[0]), carg(grafo[i].V[2] - grafo[i].V[0]) );
            break;
        }
    }
    fclose(arquivo);
}

//Exporta arquivo de texto com as informações da distribuição a Priori SCADA
void exportaCasoReferencia(GRAFO *grafo, long int numeroBarras, double Sbase)
{
    int i, j, k;

    FILE *arquivo;
    arquivo = fopen("referencia.txt", "wt");

    //----------------------------------------------------------------------
    //Fluxo de potência ativa em kW
    for (i = 0; i < numeroBarras; i++)
    {
        //Percorre os ramos adjacentes
        for (k = 0; k < grafo[i].numeroAdjacentes; k++)
        {

            switch (grafo[i].adjacentes[k].ramo->fases)
            {
            case 1:
                fprintf(arquivo, "1,0,%ld,%ld,1,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__real__ grafo[i].adjacentes[k].S[0]) * Sbase);
                break;
            case 2:
                fprintf(arquivo, "1,0,%ld,%ld,2,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__real__ grafo[i].adjacentes[k].S[1]) * Sbase);
                break;
            case 3:
                fprintf(arquivo, "1,0,%ld,%ld,3,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__real__ grafo[i].adjacentes[k].S[2]) * Sbase);
                break;
            case 4:
                fprintf(arquivo, "1,0,%ld,%ld,1,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__real__ grafo[i].adjacentes[k].S[0]) * Sbase);
                fprintf(arquivo, "1,0,%ld,%ld,2,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__real__ grafo[i].adjacentes[k].S[1]) * Sbase);
                break;
            case 5:
                fprintf(arquivo, "1,0,%ld,%ld,1,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__real__ grafo[i].adjacentes[k].S[0]) * Sbase);
                fprintf(arquivo, "1,0,%ld,%ld,3,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__real__ grafo[i].adjacentes[k].S[2]) * Sbase);
                break;
            case 6:
                fprintf(arquivo, "1,0,%ld,%ld,2,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__real__ grafo[i].adjacentes[k].S[1]) * Sbase);
                fprintf(arquivo, "1,0,%ld,%ld,3,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__real__ grafo[i].adjacentes[k].S[2]) * Sbase);
                break;
            case 7:
                fprintf(arquivo, "1,0,%ld,%ld,1,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__real__ grafo[i].adjacentes[k].S[0]) * Sbase);
                fprintf(arquivo, "1,0,%ld,%ld,2,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__real__ grafo[i].adjacentes[k].S[1]) * Sbase);
                fprintf(arquivo, "1,0,%ld,%ld,3,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__real__ grafo[i].adjacentes[k].S[2]) * Sbase);
                break;
            }
        }
    }
    //----------------------------------------------------------------------
    //Fluxo de potência reativa em kVAr
    for (i = 0; i < numeroBarras; i++)
    {
        //Percorre os ramos adjacentes
        for (k = 0; k < grafo[i].numeroAdjacentes; k++)
        {

            switch (grafo[i].adjacentes[k].ramo->fases)
            {
            case 1:
                fprintf(arquivo, "1,1,%ld,%ld,1,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__imag__ grafo[i].adjacentes[k].S[0]) * Sbase);
                break;
            case 2:
                fprintf(arquivo, "1,1,%ld,%ld,2,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__imag__ grafo[i].adjacentes[k].S[1]) * Sbase);
                break;
            case 3:
                fprintf(arquivo, "1,1,%ld,%ld,3,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__imag__ grafo[i].adjacentes[k].S[2]) * Sbase);
                break;
            case 4:
                fprintf(arquivo, "1,1,%ld,%ld,1,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__imag__ grafo[i].adjacentes[k].S[0]) * Sbase);
                fprintf(arquivo, "1,1,%ld,%ld,2,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__imag__ grafo[i].adjacentes[k].S[1]) * Sbase);
                break;
            case 5:
                fprintf(arquivo, "1,1,%ld,%ld,1,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__imag__ grafo[i].adjacentes[k].S[0]) * Sbase);
                fprintf(arquivo, "1,1,%ld,%ld,3,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__imag__ grafo[i].adjacentes[k].S[2]) * Sbase);
                break;
            case 6:
                fprintf(arquivo, "1,1,%ld,%ld,2,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__imag__ grafo[i].adjacentes[k].S[1]) * Sbase);
                fprintf(arquivo, "1,1,%ld,%ld,3,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__imag__ grafo[i].adjacentes[k].S[2]) * Sbase);
                break;
            case 7:
                fprintf(arquivo, "1,1,%ld,%ld,1,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__imag__ grafo[i].adjacentes[k].S[0]) * Sbase);
                fprintf(arquivo, "1,1,%ld,%ld,2,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__imag__ grafo[i].adjacentes[k].S[1]) * Sbase);
                fprintf(arquivo, "1,1,%ld,%ld,3,%.10lf,0.020000\n", grafo[i].barra->ID, grafo[grafo[i].adjacentes[k].idNo].barra->ID, (__imag__ grafo[i].adjacentes[k].S[2]) * Sbase);
                break;
            }
        }
    }

    //----------------------------------------------------------------------
    //Injeção de potência ativa em kW
    for (i = 0; i < numeroBarras; i++)
    {
        switch (grafo[i].fases)
        {
        case 1:
            fprintf(arquivo, "1,2,%ld,-1,1,%.10lf,0.020000\n", grafo[i].barra->ID, (__real__ grafo[i].S[0]) * Sbase);
            break;
        case 2:
            fprintf(arquivo, "1,2,%ld,-1,2,%.10lf,0.020000\n", grafo[i].barra->ID, (__real__ grafo[i].S[1]) * Sbase);
            break;
        case 3:
            fprintf(arquivo, "1,2,%ld,-1,3,%.10lf,0.020000\n", grafo[i].barra->ID, (__real__ grafo[i].S[2]) * Sbase);
            break;
        case 4:
            fprintf(arquivo, "1,2,%ld,-1,1,%.10lf,0.020000\n", grafo[i].barra->ID, (__real__ grafo[i].S[0]) * Sbase);
            fprintf(arquivo, "1,2,%ld,-1,2,%.10lf,0.020000\n", grafo[i].barra->ID, (__real__ grafo[i].S[1]) * Sbase);
            break;
        case 5:
            fprintf(arquivo, "1,2,%ld,-1,1,%.10lf,0.020000\n", grafo[i].barra->ID, (__real__ grafo[i].S[0]) * Sbase);
            fprintf(arquivo, "1,2,%ld,-1,3,%.10lf,0.020000\n", grafo[i].barra->ID, (__real__ grafo[i].S[2]) * Sbase);
            break;
        case 6:
            fprintf(arquivo, "1,2,%ld,-1,2,%.10lf,0.020000\n", grafo[i].barra->ID, (__real__ grafo[i].S[1]) * Sbase);
            fprintf(arquivo, "1,2,%ld,-1,3,%.10lf,0.020000\n", grafo[i].barra->ID, (__real__ grafo[i].S[2]) * Sbase);
            break;
        case 7:
            fprintf(arquivo, "1,2,%ld,-1,1,%.10lf,0.020000\n", grafo[i].barra->ID, (__real__ grafo[i].S[0]) * Sbase);
            fprintf(arquivo, "1,2,%ld,-1,2,%.10lf,0.020000\n", grafo[i].barra->ID, (__real__ grafo[i].S[1]) * Sbase);
            fprintf(arquivo, "1,2,%ld,-1,3,%.10lf,0.020000\n", grafo[i].barra->ID, (__real__ grafo[i].S[2]) * Sbase);
            break;
        }
    }
    //----------------------------------------------------------------------
    //Injeção de potência ativa em kVAr
    for (i = 0; i < numeroBarras; i++)
    {
        switch (grafo[i].fases)
        {
        case 1:
            fprintf(arquivo, "1,3,%ld,-1,1,%.10lf,0.020000\n", grafo[i].barra->ID, (__imag__ grafo[i].S[0]) * Sbase);
            break;
        case 2:
            fprintf(arquivo, "1,3,%ld,-1,2,%.10lf,0.020000\n", grafo[i].barra->ID, (__imag__ grafo[i].S[1]) * Sbase);
            break;
        case 3:
            fprintf(arquivo, "1,3,%ld,-1,3,%.10lf,0.020000\n", grafo[i].barra->ID, (__imag__ grafo[i].S[2]) * Sbase);
            break;
        case 4:
            fprintf(arquivo, "1,3,%ld,-1,1,%.10lf,0.020000\n", grafo[i].barra->ID, (__imag__ grafo[i].S[0]) * Sbase);
            fprintf(arquivo, "1,3,%ld,-1,2,%.10lf,0.020000\n", grafo[i].barra->ID, (__imag__ grafo[i].S[1]) * Sbase);
            break;
        case 5:
            fprintf(arquivo, "1,3,%ld,-1,1,%.10lf,0.020000\n", grafo[i].barra->ID, (__imag__ grafo[i].S[0]) * Sbase);
            fprintf(arquivo, "1,3,%ld,-1,3,%.10lf,0.020000\n", grafo[i].barra->ID, (__imag__ grafo[i].S[2]) * Sbase);
            break;
        case 6:
            fprintf(arquivo, "1,3,%ld,-1,2,%.10lf,0.020000\n", grafo[i].barra->ID, (__imag__ grafo[i].S[1]) * Sbase);
            fprintf(arquivo, "1,3,%ld,-1,3,%.10lf,0.020000\n", grafo[i].barra->ID, (__imag__ grafo[i].S[2]) * Sbase);
            break;
        case 7:
            fprintf(arquivo, "1,3,%ld,-1,1,%.10lf,0.020000\n", grafo[i].barra->ID, (__imag__ grafo[i].S[0]) * Sbase);
            fprintf(arquivo, "1,3,%ld,-1,2,%.10lf,0.020000\n", grafo[i].barra->ID, (__imag__ grafo[i].S[1]) * Sbase);
            fprintf(arquivo, "1,3,%ld,-1,3,%.10lf,0.020000\n", grafo[i].barra->ID, (__imag__ grafo[i].S[2]) * Sbase);
            break;
        }
    }
    //----------------------------------------------------------------------
    //Magnitudes de tensão em kV
    for (i = 0; i < numeroBarras; i++)
    {
        switch (grafo[i].fases)
        {
        case 1:
            fprintf(arquivo, "1,4,%ld,-1,1,%.10lf,0.010000\n", grafo[i].barra->ID, grafo[i].Vbase / 1000 * cabs(grafo[i].V[0]));
            break;
        case 2:
            fprintf(arquivo, "1,4,%ld,-1,2,%.10lf,0.010000\n", grafo[i].barra->ID, grafo[i].Vbase / 1000 * cabs(grafo[i].V[1]));
            break;
        case 3:
            fprintf(arquivo, "1,4,%ld,-1,3,%.10lf,0.010000\n", grafo[i].barra->ID, grafo[i].Vbase / 1000 * cabs(grafo[i].V[2]));
            break;
        case 4:
            fprintf(arquivo, "1,4,%ld,-1,1,%.10lf,0.010000\n", grafo[i].barra->ID, grafo[i].Vbase / 1000 * cabs(grafo[i].V[0]));
            fprintf(arquivo, "1,4,%ld,-1,2,%.10lf,0.010000\n", grafo[i].barra->ID, grafo[i].Vbase / 1000 * cabs(grafo[i].V[1]));
            break;
        case 5:
            fprintf(arquivo, "1,4,%ld,-1,1,%.10lf,0.010000\n", grafo[i].barra->ID, grafo[i].Vbase / 1000 * cabs(grafo[i].V[0]));
            fprintf(arquivo, "1,4,%ld,-1,3,%.10lf,0.010000\n", grafo[i].barra->ID, grafo[i].Vbase / 1000 * cabs(grafo[i].V[2]));
            break;
        case 6:
            fprintf(arquivo, "1,4,%ld,-1,2,%.10lf,0.010000\n", grafo[i].barra->ID, grafo[i].Vbase / 1000 * cabs(grafo[i].V[1]));
            fprintf(arquivo, "1,4,%ld,-1,3,%.10lf,0.010000\n", grafo[i].barra->ID, grafo[i].Vbase / 1000 * cabs(grafo[i].V[2]));
            break;
        case 7:
            fprintf(arquivo, "1,4,%ld,-1,1,%.10lf,0.010000\n", grafo[i].barra->ID, grafo[i].Vbase / 1000 * cabs(grafo[i].V[0]));
            fprintf(arquivo, "1,4,%ld,-1,2,%.10lf,0.010000\n", grafo[i].barra->ID, grafo[i].Vbase / 1000 * cabs(grafo[i].V[1]));
            fprintf(arquivo, "1,4,%ld,-1,3,%.10lf,0.010000\n", grafo[i].barra->ID, grafo[i].Vbase / 1000 * cabs(grafo[i].V[2]));
            break;
        }
    }
    //    //----------------------------------------------------------------------
    //    //Ângulo de tensão em graus
    for (i = 0; i < numeroBarras; i++)
    {
        switch (grafo[i].fases)
        {
        case 1:
            fprintf(arquivo, "1,5,%ld,-1,1,%.10lf,0.001000\n", grafo[i].barra->ID, 180 / PI * carg(grafo[i].V[0]));
            break;
        case 2:
            fprintf(arquivo, "1,5,%ld,-1,2,%.10lf,0.001000\n", grafo[i].barra->ID, 180 / PI * carg(grafo[i].V[1]));
            break;
        case 3:
            fprintf(arquivo, "1,5,%ld,-1,3,%.10lf,0.001000\n", grafo[i].barra->ID, 180 / PI * carg(grafo[i].V[2]));
            break;
        case 4:
            fprintf(arquivo, "1,5,%ld,-1,1,%.10lf,0.001000\n", grafo[i].barra->ID, 180 / PI * carg(grafo[i].V[0]));
            fprintf(arquivo, "1,5,%ld,-1,2,%.10lf,0.001000\n", grafo[i].barra->ID, 180 / PI * carg(grafo[i].V[1]));
            break;
        case 5:
            fprintf(arquivo, "1,5,%ld,-1,1,%.10lf,0.001000\n", grafo[i].barra->ID, 180 / PI * carg(grafo[i].V[0]));
            fprintf(arquivo, "1,5,%ld,-1,3,%.10lf,0.001000\n", grafo[i].barra->ID, 180 / PI * carg(grafo[i].V[2]));
            break;
        case 6:
            fprintf(arquivo, "1,5,%ld,-1,2,%.10lf,0.001000\n", grafo[i].barra->ID, 180 / PI * carg(grafo[i].V[1]));
            fprintf(arquivo, "1,5,%ld,-1,3,%.10lf,0.001000\n", grafo[i].barra->ID, 180 / PI * carg(grafo[i].V[2]));
            break;
        case 7:
            fprintf(arquivo, "1,5,%ld,-1,1,%.10lf,0.001000\n", grafo[i].barra->ID, 180 / PI * carg(grafo[i].V[0]));
            fprintf(arquivo, "1,5,%ld,-1,2,%.10lf,0.001000\n", grafo[i].barra->ID, 180 / PI * carg(grafo[i].V[1]));
            fprintf(arquivo, "1,5,%ld,-1,3,%.10lf,0.001000\n", grafo[i].barra->ID, 180 / PI * carg(grafo[i].V[2]));
            break;
        }
    }
    fclose(arquivo);
}

//------------------------------------------------------------------------------
//
// FUNÇÕES DO ESTIMADOR WLS
//
//------------------------------------------------------------------------------
//Função monta vetor z
void monta_z(double *z, long int nmed, DMED *medidas)
{
    int i;

    for (i = 0; i < nmed; i++)
    {
        z[i] = medidas[i].zmed;
    }
}

//Função monta matriz W
void monta_W(double **W, long int nmed, DMED *medidas)
{
    int i, j;
    double prec, menorSigma = 1000000;
    double fundoEscala, auxFE;

    //Matriz W diagonal - inverso da variância
    for (i = 0; i < nmed; i++)
    {
        auxFE = round(medidas[i].zmed * 1.25 * 1000);
        fundoEscala = auxFE / 1000;
        switch (medidas[i].tipo)
        {
        case 0:
        case 1:
            //                medidas[i].sigma = 0.002;
            //                break;
        case 2:
        case 3:
            //W[i][i] = 40000;
            //                medidas[i].sigma = 0.010;
            prec = medidas[i].prec; //0.02; //5% para SCADA de potência
            medidas[i].sigma = 0.33333 * prec * cabs(medidas[i].zmed);
            //                medidas[i].sigma = 0.33333*prec*cabs(fundoEscala);
            break;
        case 4:
        case 6:
            //W[i][i] = 40000;
            //                medidas[i].sigma = 0.001;
            prec = medidas[i].prec; //0.01; //1% para magnitude de tensão ou corrente
            medidas[i].sigma = 0.33333 * prec * cabs(medidas[i].zmed);
            //                medidas[i].sigma = 0.33333*prec*cabs(fundoEscala);
            break;
        case 5:
        case 7:
            //W[i][i] = 1000000;
            //medidas[i].sigma = 0.0001;
            prec = medidas[i].prec; //0.001; //0.1% para PMU de ângulo de tensão ou corrente
            medidas[i].sigma = 0.33333 * prec * cabs(medidas[i].zmed);
            break;
        case 8:
        case 9:
        case 10:
        case 11:
        case 12:
        case 13:
            //W[i][i] = 1000000;
            //                medidas[i].sigma = 0.0001;
            prec = medidas[i].prec; //0.001; //0.1% para PMUs retangulares
            medidas[i].sigma = 0.33333 * prec * cabs(medidas[i].zmed);
            break;
        }
        //medidas[i].sigma = 0.33333*prec*cabs(medidas[i].zmed); //Ponderação de acordo com o valor medido (Fórmula B. Pal)
        if (cabs(medidas[i].zmed) > 0.00001)
        {
            if (medidas[i].sigma < menorSigma)
                menorSigma = medidas[i].sigma;
        }
        medidas[i].sigma = medidas[i].sigma;
        //        W[i][i] = 1/(pow(medidas[i].sigma,2));
    }
    for (i = 0; i < nmed; i++)
    { //Tratamento da medida virtual e medidas proximo de zero
        if (cabs(medidas[i].zmed) < 0.00001)
        {
            medidas[i].sigma = 0.01 * menorSigma;
            //            W[i][i] = 1/(pow(medidas[i].sigma,2));
        }
    }
}

//Função monta matriz W
void monta_W_cte(double **W, long int nmed, DMED *medidas)
{
    int i, j;
    double prec, menorSigma = 1000000;
    double fundoEscala, auxFE;

    //Matriz W diagonal - inverso da variância
    for (i = 0; i < nmed; i++)
    {
        switch (medidas[i].tipo)
        {
        case 0:
        case 1:
        case 2:
        case 3:
            if (medidas[i].prec == 0.3)
                medidas[i].sigma = 0.02; //pseudo medida (20 kVA)
            if (medidas[i].prec == 0.02)
                medidas[i].sigma = 0.002; //SCADA (2 kVA)
            if (medidas[i].prec == 0.001)
                medidas[i].sigma = 0.0001; //virtual (0.1 kVA)

            break;
        case 4:
        case 6:
            medidas[i].sigma = 0.01; //SCADA de tensão (1% da tensão base)
            break;
        case 5:
        case 7:
            //W[i][i] = 1000000;
            medidas[i].sigma = 0.0001;
            prec = medidas[i].prec; //0.001; //0.1% para PMU de ângulo de tensão ou corrente
            break;
        case 8:
        case 9:
        case 10:
        case 11:
        case 12:
        case 13:
            //W[i][i] = 1000000;
            medidas[i].sigma = 0.0001;
            prec = medidas[i].prec; //0.001; //0.1% para PMUs retangulares
            break;
        }
        //medidas[i].sigma = 0.33333*prec*cabs(medidas[i].zmed); //Ponderação de acordo com o valor medido (Fórmula B. Pal)
        if (cabs(medidas[i].zmed) > 0.00001)
        {
            if (medidas[i].sigma < menorSigma)
                menorSigma = medidas[i].sigma;
        }
        W[i][i] = 1 / (pow(medidas[i].sigma, 2));
    }
    //    for(i=0;i<nmed;i++){ //Tratamento da medida virtual e medidas proximo de zero
    //        if (cabs(medidas[i].zmed) < 0.00001){
    //            medidas[i].sigma = menorSigma;
    //            W[i][i] = 1/(pow(medidas[i].sigma,2));
    //        }
    //    }
}

void monta_W_Ident(double **W, long int nmed, DMED *medidas)
{
    int i, j;

    //Matriz W diagonal - inverso da variância
    for (i = 0; i < nmed; i++)
    {
        medidas[i].sigma = 1;
        //        W[i][i] = 1;
    }
}

//Função inicialização do vetor x
void incializa_vetor_x(GRAFO *grafo, long int numeroBarras, ALIMENTADOR *alimentadores, long int numeroAlimentadores, double *x, double *regua, long int nVariaveis)
{
    int i, k, fase;
    BOOL visitado[numeroBarras];
    __complex__ double V0[3], **Yaux;

    Yaux = c_matAloca(3);

    //Flat start trifásico (Va = Vb = Vc = 1p.u.  Ta = 0  Tb = -120  Tc = 120) - com busca em profundidade para atualizar taps iniciais
    for (i = 0; i < numeroBarras; i++)
    {
        visitado[i] = false;
    }
    for (i = 0; i < numeroAlimentadores; i++)
    {
        //Tensão Inicial da subestação
        V0[0] = grafo[alimentadores[i].noRaiz].barra->Vinicial[0]; //1.0*(cos(0) + I*sin(0));
        V0[1] = grafo[alimentadores[i].noRaiz].barra->Vinicial[1]; //1.0*(cos(-120*PI/180) + I*sin(-120*PI/180));
        V0[2] = grafo[alimentadores[i].noRaiz].barra->Vinicial[2]; //1.0*(cos(120*PI/180) + I*sin(120*PI/180));

        FILABARRAS *barraAtual = &alimentadores[i].rnp[0];

        int de = barraAtual->idNo;
        grafo[de].V[0] = V0[0];
        grafo[de].V[1] = V0[1];
        grafo[de].V[2] = V0[2];

        while (barraAtual != NULL)
        {
            de = barraAtual->idNo;
            int n_adj = grafo[de].numeroAdjacentes;
            for (k = 0; k < n_adj; k++)
            {
                int para = grafo[de].adjacentes[k].idNo;
                if (visitado[para] == false)
                {
                    if ((grafo[de].adjacentes[k].tipo == 1))
                    { //Atualiza o V0 para trafo visto a ligação e tap
                        grafo[para].V[0] = grafo[de].V[0];
                        grafo[para].V[1] = grafo[de].V[1];
                        grafo[para].V[2] = grafo[de].V[2];

                        if ((grafo[de].adjacentes[k].ramo->trafo.lig_pri == 1) && (grafo[de].adjacentes[k].ramo->trafo.lig_sec == 2))
                        {
                            grafo[para].V[0] = cabs(grafo[de].V[0]) * (cos(-30 * PI / 180) + I * sin(-30 * PI / 180));
                            grafo[para].V[1] = cabs(grafo[de].V[1]) * (cos(-150 * PI / 180) + I * sin(-150 * PI / 180));
                            grafo[para].V[2] = cabs(grafo[de].V[2]) * (cos(90 * PI / 180) + I * sin(90 * PI / 180));
                        }
                        else if ((grafo[de].adjacentes[k].ramo->trafo.lig_pri == 3) && (grafo[de].adjacentes[k].ramo->trafo.lig_sec == 2))
                        {
                            grafo[para].V[0] = cabs(grafo[de].V[0]) * (cos(-30 * PI / 180) + I * sin(-30 * PI / 180));
                            grafo[para].V[1] = cabs(grafo[de].V[1]) * (cos(-150 * PI / 180) + I * sin(-150 * PI / 180));
                            grafo[para].V[2] = cabs(grafo[de].V[2]) * (cos(90 * PI / 180) + I * sin(90 * PI / 180));
                        }
                        else if ((grafo[de].adjacentes[k].ramo->trafo.lig_pri == 2) && (grafo[de].adjacentes[k].ramo->trafo.lig_sec == 1))
                        {
                            if (grafo[de].adjacentes[k].ramo->k == de)
                            {
                                grafo[para].V[0] = cabs(grafo[de].V[0]) * (cos(-30 * PI / 180) + I * sin(-30 * PI / 180));
                                grafo[para].V[1] = cabs(grafo[de].V[1]) * (cos(-150 * PI / 180) + I * sin(-150 * PI / 180));
                                grafo[para].V[2] = cabs(grafo[de].V[2]) * (cos(90 * PI / 180) + I * sin(90 * PI / 180));
                            }
                            else
                            {
                                grafo[para].V[0] = cabs(grafo[de].V[0]) * (cos(0 * PI / 180) + I * sin(0 * PI / 180));
                                grafo[para].V[1] = cabs(grafo[de].V[1]) * (cos(-120 * PI / 180) + I * sin(-120 * PI / 180));
                                grafo[para].V[2] = cabs(grafo[de].V[2]) * (cos(120 * PI / 180) + I * sin(120 * PI / 180));
                            }
                        }
                    }
                    else if (grafo[de].adjacentes[k].tipo == 2)
                    { //Para o caso de regulador de tensão
                        grafo[para].V[0] = grafo[de].V[0] * grafo[de].adjacentes[k].ramo->tap_pri[0] * grafo[de].adjacentes[k].ramo->tap_sec[0];
                        grafo[para].V[1] = grafo[de].V[1] * grafo[de].adjacentes[k].ramo->tap_pri[1] * grafo[de].adjacentes[k].ramo->tap_sec[1];
                        grafo[para].V[2] = grafo[de].V[2] * grafo[de].adjacentes[k].ramo->tap_pri[2] * grafo[de].adjacentes[k].ramo->tap_sec[2];
                    }
                    else
                    {
                        grafo[para].V[0] = grafo[de].V[0];
                        grafo[para].V[1] = grafo[de].V[1];
                        grafo[para].V[2] = grafo[de].V[2];
                    }
                }
            }
            visitado[de] = true;
            barraAtual = barraAtual->prox;
        }
    }

    //Montagem do vetor x e da régua
    for (i = 0; i < nVariaveis; i++)
    {
        k = (int)regua[i];
        fase = (int)cabs((regua[i] - k) * 10);
        if (regua[i] > 0)
        {                                  //magnitude de tensão
            x[i] = cabs(grafo[k].V[fase]); //mantém o ângulo anterior e altera a magnitude
        }
        else
        { //ângulo de tensão
            fase = fase;
            k = -k;
            x[i] = carg(grafo[k].V[fase]);
        }
    }
}


void incializa_vetor_x_FP(GRAFO *grafo,long int numeroBarras, ALIMENTADOR *alimentadores, long int numeroAlimentadores, double *x, double *regua, long int nVariaveis)
{

    //Flat start trifásico (Va = Vb = Vc = 1p.u.  Ta = 0  Tb = -120  Tc = 120) - com busca em profundidade para atualizar taps iniciais
    int i, k, fase;
    BOOL visitado[numeroBarras];
    __complex__ double V0[3], **Yaux;

    Yaux = c_matAloca(3);

    //Flat start trifásico (Va = Vb = Vc = 1p.u.  Ta = 0  Tb = -120  Tc = 120) - com busca em profundidade para atualizar taps iniciais
    for (i = 0; i < numeroBarras; i++)
    {
        visitado[i] = false;
    }
    for (i = 0; i < numeroAlimentadores; i++)
    {
        //Tensão Inicial da subestação
        V0[0] =grafo[alimentadores[i].noRaiz].barra->Vinicial[0];
        V0[1] =grafo[alimentadores[i].noRaiz].barra->Vinicial[1]; //1.0*(cos(-120*PI/180) + I*sin(-120*PI/180));
        V0[2] =grafo[alimentadores[i].noRaiz].barra->Vinicial[2];//1.0*(cos(120*PI/180) + I*sin(120*PI/180));

        FILABARRAS *barraAtual = &alimentadores[i].rnp[0];

        int de = barraAtual->idNo;
        grafo[de].V[0] = V0[0];
        grafo[de].V[1] = V0[1];
        grafo[de].V[2] = V0[2];

        while (barraAtual != NULL)
        {
            de = barraAtual->idNo;
            int n_adj = grafo[de].numeroAdjacentes;
            for (k = 0; k < n_adj; k++)
            {
                int para = grafo[de].adjacentes[k].idNo;
                if (visitado[para] == false)
                {
                    if ((grafo[de].adjacentes[k].tipo == 1))
                    { //Atualiza o V0 para trafo visto a ligação e tap
                        grafo[para].V[0] = grafo[de].V[0];
                        grafo[para].V[1] = grafo[de].V[1];
                        grafo[para].V[2] = grafo[de].V[2];

                        if ((grafo[de].adjacentes[k].ramo->trafo.lig_pri == 1) && (grafo[de].adjacentes[k].ramo->trafo.lig_sec == 2))
                        {
                            grafo[para].V[0] = cabs(grafo[de].V[0]) * (cos(-30 * PI / 180) + I * sin(-30 * PI / 180));
                            grafo[para].V[1] = cabs(grafo[de].V[1]) * (cos(-150 * PI / 180) + I * sin(-150 * PI / 180));
                            grafo[para].V[2] = cabs(grafo[de].V[2]) * (cos(90 * PI / 180) + I * sin(90 * PI / 180));
                        }
                        else if ((grafo[de].adjacentes[k].ramo->trafo.lig_pri == 3) && (grafo[de].adjacentes[k].ramo->trafo.lig_sec == 2))
                        {
                            grafo[para].V[0] = cabs(grafo[de].V[0]) * (cos(-30 * PI / 180) + I * sin(-30 * PI / 180));
                            grafo[para].V[1] = cabs(grafo[de].V[1]) * (cos(-150 * PI / 180) + I * sin(-150 * PI / 180));
                            grafo[para].V[2] = cabs(grafo[de].V[2]) * (cos(90 * PI / 180) + I * sin(90 * PI / 180));
                        }
                        else if ((grafo[de].adjacentes[k].ramo->trafo.lig_pri == 2) && (grafo[de].adjacentes[k].ramo->trafo.lig_sec == 1))
                        {
                            if (grafo[de].adjacentes[k].ramo->k == de)
                            {
                                grafo[para].V[0] = cabs(grafo[de].V[0]) * (cos(-30 * PI / 180) + I * sin(-30 * PI / 180));
                                grafo[para].V[1] = cabs(grafo[de].V[1]) * (cos(-150 * PI / 180) + I * sin(-150 * PI / 180));
                                grafo[para].V[2] = cabs(grafo[de].V[2]) * (cos(90 * PI / 180) + I * sin(90 * PI / 180));
                            }
                            else
                            {
                                grafo[para].V[0] = cabs(grafo[de].V[0]) * (cos(0 * PI / 180) + I * sin(0 * PI / 180));
                                grafo[para].V[1] = cabs(grafo[de].V[1]) * (cos(-120 * PI / 180) + I * sin(-120 * PI / 180));
                                grafo[para].V[2] = cabs(grafo[de].V[2]) * (cos(120 * PI / 180) + I * sin(120 * PI / 180));
                            }
                        }
                    }
                    else if (grafo[de].adjacentes[k].tipo == 2)
                    { //Para o caso de regulador de tensão
                        grafo[para].V[0] = grafo[de].V[0] * grafo[de].adjacentes[k].ramo->tap_pri[0] * grafo[de].adjacentes[k].ramo->tap_sec[0];
                        grafo[para].V[1] = grafo[de].V[1] * grafo[de].adjacentes[k].ramo->tap_pri[1] * grafo[de].adjacentes[k].ramo->tap_sec[1];
                        grafo[para].V[2] = grafo[de].V[2] * grafo[de].adjacentes[k].ramo->tap_pri[2] * grafo[de].adjacentes[k].ramo->tap_sec[2];
                    }
                    else
                    {
                        grafo[para].V[0] = grafo[de].V[0];
                        grafo[para].V[1] = grafo[de].V[1];
                        grafo[para].V[2] = grafo[de].V[2];
                    }
                }
            }
            visitado[de] = true;
            barraAtual = barraAtual->prox;
        }
    }

    for (i = 0; i < numeroBarras; i++)
    {
        //Tensão Inicial da subestação
        if (grafo[i].tipo==1)
        {
            grafo[i].V[0]=(grafo[i].V[0]/cabs(grafo[i].V[0]))*cabs(grafo[i].barra->Vinicial[0]);
            grafo[i].V[1]=(grafo[i].V[1]/cabs(grafo[i].V[1]))*cabs(grafo[i].barra->Vinicial[1]);
            grafo[i].V[2]=(grafo[i].V[2]/cabs(grafo[i].V[2]))*cabs(grafo[i].barra->Vinicial[2]);
        }

    }

    //Montagem do vetor x e da régua
    for (i = 0; i < nVariaveis; i++)
    {
        k = (int)regua[i];
        fase = (int)cabs((regua[i] - k) * 10);
        if (regua[i] > 0)
        {                                  //magnitude de tensão
            x[i] = cabs(grafo[k].V[fase]); //mantém o ângulo anterior e altera a magnitude
        }
        else
        { //ângulo de tensão
            fase = fase;
            k = -k;
            x[i] = carg(grafo[k].V[fase]);
        }
    }
}


void incializa_vetor_x_perturbado(GRAFO *grafo, long int numeroBarras, ALIMENTADOR *alimentadores, long int numeroAlimentadores, double *x, double *regua, long int nVariaveis)
{
    int i, k, fase;
    BOOL visitado[numeroBarras];
    __complex__ double V0[3], **Yaux;

    Yaux = c_matAloca(3);

    //Flat start trifásico (Va = Vb = Vc = 1p.u.  Ta = 0  Tb = -120  Tc = 120) - com busca em profundidade para atualizar taps iniciais
    for (i = 0; i < numeroBarras; i++)
    {
        visitado[i] = false;
    }
    for (i = 0; i < numeroAlimentadores; i++)
    {
        //Tensão Inicial da subestação
        V0[0] = grafo[alimentadores[i].noRaiz].barra->Vinicial[0]; //1.0*(cos(0) + I*sin(0));
        V0[1] = grafo[alimentadores[i].noRaiz].barra->Vinicial[1]; //1.0*(cos(-120*PI/180) + I*sin(-120*PI/180));
        V0[2] = grafo[alimentadores[i].noRaiz].barra->Vinicial[2]; //1.0*(cos(120*PI/180) + I*sin(120*PI/180));

        FILABARRAS *barraAtual = &alimentadores[i].rnp[0];

        int de = barraAtual->idNo;
        grafo[de].V[0] = V0[0];
        grafo[de].V[1] = V0[1];
        grafo[de].V[2] = V0[2];

        while (barraAtual != NULL)
        {
            de = barraAtual->idNo;
            int n_adj = grafo[de].numeroAdjacentes;
            for (k = 0; k < n_adj; k++)
            {
                int para = grafo[de].adjacentes[k].idNo;
                if (visitado[para] == false)
                {
                    if ((grafo[de].adjacentes[k].tipo == 1))
                    { //Atualiza o V0 para trafo visto a ligação e tap
                        grafo[para].V[0] = grafo[de].V[0];
                        grafo[para].V[1] = grafo[de].V[1];
                        grafo[para].V[2] = grafo[de].V[2];


                        if ((grafo[de].adjacentes[k].ramo->trafo.lig_pri == 1) && (grafo[de].adjacentes[k].ramo->trafo.lig_sec == 2))
                        {
                            grafo[para].V[0] = cabs(grafo[de].V[0]) * (cos(-30 * PI / 180) + I * sin(-30 * PI / 180));
                            grafo[para].V[1] = cabs(grafo[de].V[1]) * (cos(-150 * PI / 180) + I * sin(-150 * PI / 180));
                            grafo[para].V[2] = cabs(grafo[de].V[2]) * (cos(90 * PI / 180) + I * sin(90 * PI / 180));
                        }
                        else if ((grafo[de].adjacentes[k].ramo->trafo.lig_pri == 3) && (grafo[de].adjacentes[k].ramo->trafo.lig_sec == 2))
                        {
                            grafo[para].V[0] = cabs(grafo[de].V[0]) * (cos(-30 * PI / 180) + I * sin(-30 * PI / 180));
                            grafo[para].V[1] = cabs(grafo[de].V[1]) * (cos(-150 * PI / 180) + I * sin(-150 * PI / 180));
                            grafo[para].V[2] = cabs(grafo[de].V[2]) * (cos(90 * PI / 180) + I * sin(90 * PI / 180));
                        }
                        else if ((grafo[de].adjacentes[k].ramo->trafo.lig_pri == 2) && (grafo[de].adjacentes[k].ramo->trafo.lig_sec == 1))
                        {
                            if (grafo[de].adjacentes[k].ramo->k == de)
                            {
                                grafo[para].V[0] = cabs(grafo[de].V[0]) * (cos(-30 * PI / 180) + I * sin(-30 * PI / 180));
                                grafo[para].V[1] = cabs(grafo[de].V[1]) * (cos(-150 * PI / 180) + I * sin(-150 * PI / 180));
                                grafo[para].V[2] = cabs(grafo[de].V[2]) * (cos(90 * PI / 180) + I * sin(90 * PI / 180));
                            }
                            else
                            {
                                grafo[para].V[0] = cabs(grafo[de].V[0]) * (cos(0 * PI / 180) + I * sin(0 * PI / 180));
                                grafo[para].V[1] = cabs(grafo[de].V[1]) * (cos(-120 * PI / 180) + I * sin(-120 * PI / 180));
                                grafo[para].V[2] = cabs(grafo[de].V[2]) * (cos(120 * PI / 180) + I * sin(120 * PI / 180));
                            }
                        }
                    }
                    else if (grafo[de].adjacentes[k].tipo == 2)
                    { //Para o caso de regulador de tensão
                        grafo[para].V[0] = grafo[de].V[0] * grafo[de].adjacentes[k].ramo->tap_pri[0] * grafo[de].adjacentes[k].ramo->tap_sec[0];
                        grafo[para].V[1] = grafo[de].V[1] * grafo[de].adjacentes[k].ramo->tap_pri[1] * grafo[de].adjacentes[k].ramo->tap_sec[1];
                        grafo[para].V[2] = grafo[de].V[2] * grafo[de].adjacentes[k].ramo->tap_pri[2] * grafo[de].adjacentes[k].ramo->tap_sec[2];
                    }
                    else
                    {
                        grafo[para].V[0] = grafo[de].V[0];
                        grafo[para].V[1] = grafo[de].V[1];
                        grafo[para].V[2] = grafo[de].V[2];
                    }
                    grafo[para].V[0] = grafo[para].V[0]+((double)rand()/RAND_MAX*2.0-1.0)/1000+I*((double)rand()/RAND_MAX*2.0-1.0)/1000;
                    grafo[para].V[1] = grafo[para].V[1]+((double)rand()/RAND_MAX*2.0-1.0)/1000+I*((double)rand()/RAND_MAX*2.0-1.0)/1000;
                    grafo[para].V[2] = grafo[para].V[2]+((double)rand()/RAND_MAX*2.0-1.0)/1000+I*((double)rand()/RAND_MAX*2.0-1.0)/1000;
                }
            }
            visitado[de] = true;
            barraAtual = barraAtual->prox;
        }
    }

    //Montagem do vetor x e da régua
    for (i = 0; i < nVariaveis; i++)
    {
        k = (int)regua[i];
        fase = (int)cabs((regua[i] - k) * 10);
        if (regua[i] > 0)
        {                                  //magnitude de tensão
            x[i] = cabs(grafo[k].V[fase]); //mantém o ângulo anterior e altera a magnitude
        }
        else
        { //ângulo de tensão
            fase = fase;
            k = -k;
            x[i] = carg(grafo[k].V[fase]);
        }
    }
}

//Função inicialização do vetor x a partir de arquivo de entrada
void incializa_vetor_x_leitura(GRAFO *grafo, long int numeroBarras, ALIMENTADOR *alimentadores, long int numeroAlimentadores, double *x, double *regua, long int nVariaveis)
{
    int i, k, fase;
    BOOL visitado[numeroBarras];

    //Inicialização com base nas tensões lidas no arquivo Vinicial.csv
    for (i = 0; i < numeroBarras; i++)
    {
        visitado[i] = false;
    }
    for (i = 0; i < numeroAlimentadores; i++)
    {
        //Tensão Inicial da subestação
        FILABARRAS *barraAtual = &alimentadores[i].rnp[0];
        while (barraAtual != NULL)
        {
            int de = barraAtual->idNo;
            if (visitado[de] == false)
            {
                grafo[de].V[0] = grafo[de].barra->Vinicial[0];
                grafo[de].V[1] = grafo[de].barra->Vinicial[1];
                grafo[de].V[2] = grafo[de].barra->Vinicial[2];

                visitado[de] = true;
                barraAtual = barraAtual->prox;
            }
        }
    }

    //Montagem do vetor x e da régua
    for (i = 0; i < nVariaveis; i++)
    {
        k = (int)regua[i];
        fase = (int)cabs((regua[i] - k) * 10);
        if (regua[i] > 0)
        {                                  //magnitude de tensão
            x[i] = cabs(grafo[k].V[fase]); //mantém o ângulo anterior e altera a magnitude
        }
        else
        { //ângulo de tensão
            fase = fase;
            k = -k;
            x[i] = carg(grafo[k].V[fase]);
        }
    }
}

//Função atualiza as grandezas elétricas da rede (fluxos, injeções, correntes)
void atualiza_Rede(GRAFO *grafo, long int numeroBarras)
{
    int i, j, k, idMed, de, para, ramo, fase;
    __complex__ double *Saux, *Iaux;
    BOOL visitado[numeroBarras];

    Saux = (__complex__ double *)malloc(3 * sizeof(__complex__ double));
    Saux[0] = 0;
    Saux[1] = 0;
    Saux[2] = 0;
    Iaux = (__complex__ double *)malloc(3 * sizeof(__complex__ double));
    Iaux[0] = 0;
    Iaux[1] = 0;
    Iaux[2] = 0;

    //Percorre o grafo atualizando o cálculo de h(x))
    for (i = 0; i < numeroBarras; i++)
    {
        Sk(grafo, i, Saux);
        grafo[i].S[0] = Saux[0];
        grafo[i].S[1] = Saux[1];
        grafo[i].S[2] = Saux[2];

        grafo[i].Cur[0] = conj(Saux[0] / grafo[i].V[0]);
        grafo[i].Cur[1] = conj(Saux[1] / grafo[i].V[1]);
        grafo[i].Cur[2] = conj(Saux[2] / grafo[i].V[2]);

        //Percorre os ramos adjacentes
        for (k = 0; k < grafo[i].numeroAdjacentes; k++)
        {
            if (i == grafo[i].adjacentes[k].ramo->k)
            {
                Skm(&grafo[i], &grafo[grafo[i].adjacentes[k].idNo], grafo[i].adjacentes[k].ramo, Saux);
                Ikm(&grafo[i], &grafo[grafo[i].adjacentes[k].idNo], grafo[i].adjacentes[k].ramo, Iaux);
            }
            else
            {
                Smk(&grafo[grafo[i].adjacentes[k].idNo], &grafo[i], grafo[i].adjacentes[k].ramo, Saux);
                Imk(&grafo[grafo[i].adjacentes[k].idNo], &grafo[i], grafo[i].adjacentes[k].ramo, Iaux);
            }
            grafo[i].adjacentes[k].S[0] = Saux[0];
            grafo[i].adjacentes[k].S[1] = Saux[1];
            grafo[i].adjacentes[k].S[2] = Saux[2];

            grafo[i].adjacentes[k].Cur[0] = Iaux[0];
            grafo[i].adjacentes[k].Cur[1] = Iaux[1];
            grafo[i].adjacentes[k].Cur[2] = Iaux[2];
        }
    }
    free(Saux);
    free(Iaux);
}

//Função atualiza modelo de medição conforme as grandezas elétricas calculadas
void atualiza_Modelo(GRAFO *grafo, long int numeroBarras, long int nmed, DMED *medidas)
{
    int i, j, k, idMed, de, para, ramo, fase;
    __complex__ double *Saux, *Iaux;
    BOOL visitado[numeroBarras];

    for (idMed = 0; idMed < nmed; idMed++)
    {
        if (medidas[idMed].PARA == -1)
        { //Medidor instalado em uma barra
            i = medidas[idMed].k;
            fase = medidas[idMed].fases - 1; //Revisar para trifásico genérico - atual somente A ou B ou C - medida no Delta por exemplo

            switch (medidas[idMed].tipo)
            {
            case 2: //2: Injeção de Potência Ativa - kW
                medidas[idMed].h = __real__ grafo[i].S[fase];
                break;
            case 3: //3: Injeção de Potência Reativa - kVAr
                medidas[idMed].h = __imag__ grafo[i].S[fase];
                break;
            case 4: //4: Magnitude de Tensão - kV
                medidas[idMed].h = cabs(grafo[i].V[fase]);
                break;
            case 5: //5: Ângulo de Tensão - graus
                medidas[idMed].h = carg(grafo[i].V[fase]);
                break;
            case 8: //8: Injeção Magnitude de Corrente - A  //REVISAR A MEDIDA DE INJEÇÂO DE PMU NO ESTIMADOR HIBRIDO
                medidas[idMed].h = cabs(grafo[i].Cur[fase]);
                break;
            case 9: //9: Injeção Ângulo de Corrente) - graus
                medidas[idMed].h = carg(grafo[i].Cur[fase]);
                break;
            }
        }
        else
        { //Medidor instalado em um ramo
            i = medidas[idMed].k;
            fase = medidas[idMed].fases - 1; //Revisar para trifásico genérico - atual somente A ou B ou C - medida no Delta por exemplo
            for (k = 0; k < grafo[i].numeroAdjacentes; k++)
            {
                if (grafo[i].adjacentes[k].idNo == medidas[idMed].m)
                {
                    switch (medidas[idMed].tipo)
                    {
                    case 0: //0: Fluxo de Potência Ativa - kW
                        medidas[idMed].h = __real__ grafo[i].adjacentes[k].S[fase];
                        break;
                    case 1: //1: Fluxo de Potência Reativa - kVAr
                        medidas[idMed].h = __imag__ grafo[i].adjacentes[k].S[fase];
                        break;
                    case 6: //6: Fluxo Magnitude de Corrente - A
                        medidas[idMed].h = cabs(grafo[i].adjacentes[k].Cur[fase]);
                        break;
                    case 7: //7: Fluxo Ângulo de Corrente) - graus
                        medidas[idMed].h = carg(grafo[i].adjacentes[k].Cur[fase]);
                        break;
                    }
                }
            }
        }
    }
}

//Função atualiza vetor h
void atualiza_h(GRAFO *grafo, long int numeroBarras, long int nmed, DMED *medidas)
{
    int i, j, k, idMed, de, para, ramo, fase;
    __complex__ double *Saux, *Iaux;
    BOOL visitado[numeroBarras];

    Saux = (__complex__ double *)malloc(3 * sizeof(__complex__ double));
    Saux[0] = 0;
    Saux[1] = 0;
    Saux[2] = 0;
    Iaux = (__complex__ double *)malloc(3 * sizeof(__complex__ double));
    Iaux[0] = 0;
    Iaux[1] = 0;
    Iaux[2] = 0;

    //Percorre o grafo atualizando o cálculo de h(x))
    for (i = 0; i < numeroBarras; i++)
    {
        //Percorre os medidores da barra
        if (grafo[i].nmed > 0)
        {
            Sk(grafo, i, Saux);
            grafo[i].S[0] = Saux[0];
            grafo[i].S[1] = Saux[1];
            grafo[i].S[2] = Saux[2];

            grafo[i].Cur[0] = conj(Saux[0] / grafo[i].V[0]);
            grafo[i].Cur[1] = conj(Saux[1] / grafo[i].V[1]);
            grafo[i].Cur[2] = conj(Saux[2] / grafo[i].V[2]);
        }

        for (j = 0; j < grafo[i].nmed; j++)
        {
            idMed = grafo[i].medidores[j]->id;
            fase = grafo[i].medidores[j]->fases - 1; //Revisar para trifásico genérico - atual somente A ou B ou C

            switch (medidas[idMed].tipo)
            {
            case 2: //2: Injeção de Potência Ativa - kW
                medidas[idMed].h = __real__ grafo[i].S[fase];
                break;
            case 3: //3: Injeção de Potência Reativa - kVAr
                medidas[idMed].h = __imag__ grafo[i].S[fase];
                break;
            case 4: //4: Magnitude de Tensão - kV
                medidas[idMed].h = cabs(grafo[i].V[fase]);
                break;
            case 5: //5: Ângulo de Tensão - graus
                medidas[idMed].h = carg(grafo[i].V[fase]);
                break;
            case 8: //8: Tensão real PMU
                medidas[idMed].h = __real__ grafo[i].V[fase];
                break;
            case 9: //8: Tensão imag PMU
                medidas[idMed].h = __imag__ grafo[i].V[fase];
                break;
            case 10: //10: Injeção Corrente real PMU
                medidas[idMed].h = __real__ grafo[i].Cur[fase];
                break;
            case 11: //11: Injeção Corrente imag PMU
                medidas[idMed].h = __imag__ grafo[i].Cur[fase];
                break;
            }
        }

        //Percorre os medidores dos ramos adjacentes
        for (k = 0; k < grafo[i].numeroAdjacentes; k++)
        {
            if (grafo[i].adjacentes[k].nmed > 0)
            {
                if (i == grafo[i].adjacentes[k].ramo->k)
                {
                    Skm(&grafo[i], &grafo[grafo[i].adjacentes[k].idNo], grafo[i].adjacentes[k].ramo, Saux);
                    Ikm(&grafo[i], &grafo[grafo[i].adjacentes[k].idNo], grafo[i].adjacentes[k].ramo, Iaux);
                }
                else
                {
                    Smk(&grafo[grafo[i].adjacentes[k].idNo], &grafo[i], grafo[i].adjacentes[k].ramo, Saux);
                    Imk(&grafo[grafo[i].adjacentes[k].idNo], &grafo[i], grafo[i].adjacentes[k].ramo, Iaux);
                }
                grafo[i].adjacentes[k].S[0] = Saux[0];
                grafo[i].adjacentes[k].S[1] = Saux[1];
                grafo[i].adjacentes[k].S[2] = Saux[2];

                grafo[i].adjacentes[k].Cur[0] = Iaux[0];
                grafo[i].adjacentes[k].Cur[1] = Iaux[1];
                grafo[i].adjacentes[k].Cur[2] = Iaux[2];
            }

            for (j = 0; j < grafo[i].adjacentes[k].nmed; j++)
            {
                idMed = grafo[i].adjacentes[k].medidores[j]->id;
                fase = grafo[i].adjacentes[k].medidores[j]->fases - 1; //Revisar para trifásico genérico - atual somente A ou B ou C

                switch (medidas[idMed].tipo)
                {
                case 0: //0: Fluxo de Potência Ativa - kW
                    medidas[idMed].h = __real__ grafo[i].adjacentes[k].S[fase];
                    break;
                case 1: //1: Fluxo de Potência Reativa - kVAr
                    medidas[idMed].h = __imag__ grafo[i].adjacentes[k].S[fase];
                    break;
                case 6: //6: Fluxo Magnitude de Corrente - A
                    medidas[idMed].h = cabs(grafo[i].adjacentes[k].Cur[fase]);
                    break;
                case 7: //7: Fluxo Ângulo de Corrente) - graus
                    medidas[idMed].h = carg(grafo[i].adjacentes[k].Cur[fase]);
                    break;
                case 12: //12: Corrente Real PMU
                    medidas[idMed].h = __real__ grafo[i].adjacentes[k].Cur[fase];
                    break;
                case 13: //12: Corrente Imag PMU
                    medidas[idMed].h = __imag__ grafo[i].adjacentes[k].Cur[fase];
                    break;
                }
            }
        }
    }
    free(Saux);
    free(Iaux);
}

//Função atualiza matriz H
void atualiza_H(GRAFO *grafo, long int numeroBarras, DRAM *ramos, DMED *medidas, long int numeroMedidas)
{
    int i, j, k, idMed, de, para, adj, ramo, opt, fase;
    __complex__ double *dSaux, *dIaux;
    int it = 0;

    dSaux = (__complex__ double *)malloc(3 * sizeof(__complex__ double));
    dSaux[0] = 0;
    dSaux[1] = 0;
    dSaux[2] = 0;
    dIaux = (__complex__ double *)malloc(3 * sizeof(__complex__ double));
    dIaux[0] = 0;
    dIaux[1] = 0;
    dIaux[2] = 0;

    //Atualiza a matriz H de acordo com a régua
    for (i = 0; i < numeroMedidas; i++)
    {
        switch (medidas[i].tipo)
        {
        case 0: //0: Fluxo de Potência Ativa - kW
            for (j = 0; j < medidas[i].nvar; j++)
            {
                k = (int)medidas[i].reguaH[j];
                fase = (int)cabs((medidas[i].reguaH[j] - k) * 10);

                de = medidas[i].k;
                para = medidas[i].m;
                ramo = medidas[i].ramo;

                if (ramos[ramo].k == de)
                {
                    if (medidas[i].reguaH[j] >= 0)
                    {
                        if (de == k)
                            opt = 0;
                        else
                            opt = 2;
                    }
                    else
                    {
                        if (de == -k)
                            opt = 1;
                        else
                            opt = 3;
                    }
                    dSkm(&grafo[de], &grafo[para], &ramos[ramo], dSaux, opt, fase);
                    medidas[i].H[j] = __real__ dSaux[medidas[i].fases - 1];
                }
                else
                {
                    if (medidas[i].reguaH[j] >= 0)
                    {
                        if (de == k)
                            opt = 2;
                        else
                            opt = 0;
                    }
                    else
                    {
                        if (de == -k)
                            opt = 3;
                        else
                            opt = 1;
                    }
                    dSmk(&grafo[para], &grafo[de], &ramos[ramo], dSaux, opt, fase);
                    medidas[i].H[j] = __real__ dSaux[medidas[i].fases - 1];
                }
            }
            break;
        case 1: //1: Fluxo de Potência Reativa - kVAr
            for (j = 0; j < medidas[i].nvar; j++)
            {
                k = (int)medidas[i].reguaH[j];
                fase = (int)cabs((medidas[i].reguaH[j] - k) * 10);

                de = medidas[i].k;
                para = medidas[i].m;
                ramo = medidas[i].ramo;

                if (ramos[ramo].k == de)
                {
                    if (medidas[i].reguaH[j] >= 0)
                    {
                        if (de == k)
                            opt = 0;
                        else
                            opt = 2;
                    }
                    else
                    {
                        if (de == -k)
                            opt = 1;
                        else
                            opt = 3;
                    }
                    dSkm(&grafo[de], &grafo[para], &ramos[ramo], dSaux, opt, fase);
                    medidas[i].H[j] = __imag__ dSaux[medidas[i].fases - 1];
                }
                else
                {
                    if (medidas[i].reguaH[j] >= 0)
                    {
                        if (de == k)
                            opt = 2;
                        else
                            opt = 0;
                    }
                    else
                    {
                        if (de == -k)
                            opt = 3;
                        else
                            opt = 1;
                    }
                    dSmk(&grafo[para], &grafo[de], &ramos[ramo], dSaux, opt, fase);
                    medidas[i].H[j] = __imag__ dSaux[medidas[i].fases - 1];
                }
            }
            break;
        case 2: //2: Injeção de Potência Ativa - kW
            for (j = 0; j < medidas[i].nvar; j++)
            {
                k = (int)medidas[i].reguaH[j];
                fase = (int)cabs((medidas[i].reguaH[j] - k) * 10);

                de = medidas[i].k;

                if (medidas[i].reguaH[j] >= 0)
                {
                    opt = 0;
                }
                else
                {
                    opt = 1;
                }
                dSk(grafo, de, dSaux, opt, cabs(k), fase);

                medidas[i].H[j] = __real__ dSaux[medidas[i].fases - 1];
            }
            break;
        case 3: //3: Injeção de Potência Reativa - kVAr
            for (j = 0; j < medidas[i].nvar; j++)
            {
                k = (int)medidas[i].reguaH[j];
                fase = (int)cabs((medidas[i].reguaH[j] - k) * 10);

                de = medidas[i].k;

                if (medidas[i].reguaH[j] >= 0)
                {
                    opt = 0;
                }
                else
                {
                    opt = 1;
                }
                dSk(grafo, de, dSaux, opt, cabs(k), fase);

                medidas[i].H[j] = __imag__ dSaux[medidas[i].fases - 1];
            }
            break;
        case 4: //4: Magnitude de tensão
            for (j = 0; j < medidas[i].nvar; j++)
            {
                k = (int)medidas[i].reguaH[j];
                fase = (int)cabs((medidas[i].reguaH[j] - k) * 10);

                if ((medidas[i].fases - 1 == fase) && (medidas[i].reguaH[j] > 0))
                    medidas[i].H[j] = 1;
            }
            break;
        case 5: //5: Ângulo de tensão
            for (j = 0; j < medidas[i].nvar; j++)
            {
                k = (int)medidas[i].reguaH[j];
                fase = (int)cabs((medidas[i].reguaH[j] - k) * 10);

                if ((medidas[i].fases - 1 == fase) && (medidas[i].reguaH[j] < 0))
                    medidas[i].H[j] = 1;
            }
            break;
        case 7: //7: Magnitude de Corrente
            for (j = 0; j < medidas[i].nvar; j++)
            {
                k = (int)medidas[i].reguaH[j];
                fase = (int)cabs((medidas[i].reguaH[j] - k) * 10);

                de = medidas[i].k;
                para = medidas[i].m;
                ramo = medidas[i].ramo;

                if (ramos[ramo].k == de)
                {
                    if (medidas[i].reguaH[j] >= 0)
                    {
                        if (de == k)
                            opt = 0;
                        else
                            opt = 2;
                    }
                    else
                    {
                        if (de == -k)
                            opt = 1;
                        else
                            opt = 3;
                    }
                    if (it == 0)
                        dSkm(&grafo[de], &grafo[para], &ramos[ramo], dSaux, opt, fase);
                    else
                        dIkm(&grafo[de], &grafo[para], &ramos[ramo], dSaux, opt, fase);
                    medidas[i].H[j] = dSaux[medidas[i].fases - 1];
                }
                else
                {
                    if (medidas[i].reguaH[j] >= 0)
                    {
                        if (de == k)
                            opt = 2;
                        else
                            opt = 0;
                    }
                    else
                    {
                        if (de == -k)
                            opt = 3;
                        else
                            opt = 1;
                    }
                    dSmk(&grafo[para], &grafo[de], &ramos[ramo], dSaux, opt, fase);
                    medidas[i].H[j] = __real__ dSaux[medidas[i].fases - 1];
                }
            }
            break;
        }
    }
    free(dSaux);
    free(dIaux);
}

//Função atualiza matriz H
void atualiza_H_ret(GRAFO *grafo, long int numeroBarras, DRAM *ramos, DMED *medidas, long int numeroMedidas)
{
    int i, j, k, idMed, de, para, adj, ramo, opt, fase;
    __complex__ double *dSaux, *dIaux;
    BOOL visitado[numeroBarras];

    dSaux = (__complex__ double *)malloc(3 * sizeof(__complex__ double));
    dSaux[0] = 0;
    dSaux[1] = 0;
    dSaux[2] = 0;
    dIaux = (__complex__ double *)malloc(3 * sizeof(__complex__ double));
    dIaux[0] = 0;
    dIaux[1] = 0;
    dIaux[2] = 0;

    //Atualiza a matriz H de acordo com a régua
    for (i = 0; i < numeroMedidas; i++)
    {
        switch (medidas[i].tipo)
        {
        case 0: //0: Fluxo de Potência Ativa - kW
            for (j = 0; j < medidas[i].nvar; j++)
            {
                k = (int)medidas[i].reguaH[j];
                fase = (int)cabs((medidas[i].reguaH[j] - k) * 10);

                de = medidas[i].k;
                para = medidas[i].m;
                ramo = medidas[i].ramo;

                if (ramos[ramo].k == de)
                {
                    if (medidas[i].reguaH[j] >= 0)
                    {
                        if (de == k)
                            opt = 0;
                        else
                            opt = 2;
                    }
                    else
                    {
                        if (de == -k)
                            opt = 1;
                        else
                            opt = 3;
                    }
                    dSkm_ret(&grafo[de], &grafo[para], &ramos[ramo], dSaux, opt, fase);
                    medidas[i].H[j] = __real__ dSaux[medidas[i].fases - 1];
                }
                else
                {
                    if (medidas[i].reguaH[j] >= 0)
                    {
                        if (de == k)
                            opt = 2;
                        else
                            opt = 0;
                    }
                    else
                    {
                        if (de == -k)
                            opt = 3;
                        else
                            opt = 1;
                    }
                    dSmk_ret(&grafo[para], &grafo[de], &ramos[ramo], dSaux, opt, fase);
                    medidas[i].H[j] = __real__ dSaux[medidas[i].fases - 1];
                }
            }
            break;
        case 1: //1: Fluxo de Potência Reativa - kVAr
            for (j = 0; j < medidas[i].nvar; j++)
            {
                k = (int)medidas[i].reguaH[j];
                fase = (int)cabs((medidas[i].reguaH[j] - k) * 10);

                de = medidas[i].k;
                para = medidas[i].m;
                ramo = medidas[i].ramo;

                if (ramos[ramo].k == de)
                {
                    if (medidas[i].reguaH[j] >= 0)
                    {
                        if (de == k)
                            opt = 0;
                        else
                            opt = 2;
                    }
                    else
                    {
                        if (de == -k)
                            opt = 1;
                        else
                            opt = 3;
                    }
                    dSkm_ret(&grafo[de], &grafo[para], &ramos[ramo], dSaux, opt, fase);
                    medidas[i].H[j] = __imag__ dSaux[medidas[i].fases - 1];
                }
                else
                {
                    if (medidas[i].reguaH[j] >= 0)
                    {
                        if (de == k)
                            opt = 2;
                        else
                            opt = 0;
                    }
                    else
                    {
                        if (de == -k)
                            opt = 3;
                        else
                            opt = 1;
                    }
                    dSmk_ret(&grafo[para], &grafo[de], &ramos[ramo], dSaux, opt, fase);
                    medidas[i].H[j] = __imag__ dSaux[medidas[i].fases - 1];
                }
            }
            break;
        case 2: //2: Injeção de Potência Ativa - kW
            for (j = 0; j < medidas[i].nvar; j++)
            {
                k = (int)medidas[i].reguaH[j];
                fase = (int)cabs((medidas[i].reguaH[j] - k) * 10);

                de = medidas[i].k;

                if (medidas[i].reguaH[j] >= 0)
                {
                    opt = 0;
                }
                else
                {
                    opt = 1;
                }
                dSk_ret(grafo, de, dSaux, opt, cabs(k), fase);

                medidas[i].H[j] = __real__ dSaux[medidas[i].fases - 1];
            }
            break;
        case 3: //3: Injeção de Potência Reativa - kVAr
            for (j = 0; j < medidas[i].nvar; j++)
            {
                k = (int)medidas[i].reguaH[j];
                fase = (int)cabs((medidas[i].reguaH[j] - k) * 10);

                de = medidas[i].k;

                if (medidas[i].reguaH[j] >= 0)
                {
                    opt = 0;
                }
                else
                {
                    opt = 1;
                }
                dSk_ret(grafo, de, dSaux, opt, cabs(k), fase);

                medidas[i].H[j] = __imag__ dSaux[medidas[i].fases - 1];
            }
            break;
        case 4: //4: Magnitude de tensão
            for (j = 0; j < medidas[i].nvar; j++)
            {
                k = (int)medidas[i].reguaH[j];
                fase = (int)cabs((medidas[i].reguaH[j] - k) * 10);

                if (medidas[i].reguaH[j] >= 0)
                {
                    if (medidas[i].fases - 1 == fase)
                        medidas[i].H[j] = __real__ grafo[k].V[fase] / cabs(grafo[k].V[fase]);
                }
                else
                {
                    k = -k;
                    if (medidas[i].fases - 1 == fase)
                        medidas[i].H[j] = __imag__ grafo[k].V[fase] / cabs(grafo[k].V[fase]);
                }
            }
            break;
        case 5: //5: Magnitude de tensão
            for (j = 0; j < medidas[i].nvar; j++)
            {
                k = (int)medidas[i].reguaH[j];
                fase = (int)cabs((medidas[i].reguaH[j] - k) * 10);

                if (medidas[i].reguaH[j] >= 0)
                {
                    if (medidas[i].fases - 1 == fase)
                        medidas[i].H[j] = -1 * __imag__ grafo[k].V[fase] / pow(__real__ grafo[k].V[fase], 2) * 1 / (1 + pow(__imag__ grafo[k].V[fase] / __real__ grafo[k].V[fase], 2));
                }
                else
                {
                    k = -k;
                    if (medidas[i].fases - 1 == fase)
                        medidas[i].H[j] = 1 / __real__ grafo[k].V[fase] * 1 / (1 + pow(__imag__ grafo[k].V[fase] / __real__ grafo[k].V[fase], 2));
                }
            }
            break;
        }
    }
    free(dSaux);
    free(dIaux);
}

//Função que atualiza as variaveis de estado em função do vetor x
void atualiza_estado(GRAFO *grafo, double *x, double *regua, long int nVariaveis)
{
    int i, k, fase;

    //nVariaveis é o tamanho do vetor x
    //regua salva a respectiva barra e fase de cada elemento do vetor x (+ é magnitude de tensão e - é ângulo)
    // Ex: "+1.2" = tensão barra 1 fase 2
    // Ex: "-2.0" = ângulo barra 2 fase 0

    for (i = 0; i < nVariaveis; i++)
    {
        k = (int)regua[i];
        fase = (int)cabs((regua[i] - k) * 10);
        if (regua[i] > 0)
        {                                                                        //magnitude de tensão
            grafo[k].V[fase] = x[i] * grafo[k].V[fase] / cabs(grafo[k].V[fase]); //mantém o ângulo anterior e altera a magnitude
        }
        else
        { //ângulo de tensão
            fase = fase;
            k = -k;
            grafo[k].V[fase] = cabs(grafo[k].V[fase]) * (cos(x[i]) + I * sin(x[i])); //mantém a magnitude e altera o ângulo
        }
    }
}

void add_pseudo(ALIMENTADOR *alimentador,GRAFO *grafo,DRAM *ramos, DMED *medidas, long int nmed)
{
        printf("pseudo medida foi adicionada\n\n");
        int medB=nmed-2,medC=nmed-1;
        int ind;
        int k;
        int j;
        double regua;

        medidas[medB].ligado =1; //se a medida esta ativa
        medidas[medB].tipo = 5;// tipo de medida 
        medidas[medB].DE = grafo[alimentador->noRaiz].barra->ID;// barra DE conforme nomeclatura do DBAR
        medidas[medB].PARA = -1;// barra PARA conforme a nomeclatura do DBAR
        medidas[medB].fases = 2;// De qual fase é a medida
        medidas[medB].id = medB;// id da medida
        medidas[medB].par = -1;// par ???
        medidas[medB].h = 0; //valor do h(x)
        medidas[medB].zmed = -120 * PI / 180;// valor medida
        medidas[medB].sigma = 1; //atof(getfield(dados,7));// sigma, desvio padrão calculado internamente
        medidas[medB].prec = 0.042;// precisão do medidor
        medidas[medB].k=alimentador->noRaiz;
        medidas[medB].m=-1;
        medidas[medB].ramo=-1;

        medidas[medC].ligado =1; //se a medida esta ativa
        medidas[medC].tipo = 5;// tipo de medida 
        medidas[medC].DE = grafo[alimentador->noRaiz].barra->ID;// barra DE conforme nomeclatura do DBAR
        medidas[medC].PARA = -1;// barra PARA conforme a nomeclatura do DBAR
        medidas[medC].fases = 3;// De qual fase é a medida
        medidas[medC].id = medC;// id da medida
        medidas[medC].par = -1;// par ???
        medidas[medC].h = 0; //valor do h(x)
        medidas[medC].zmed = 120 * PI / 180;// valor medida
        medidas[medC].sigma = 1; //atof(getfield(dados,7));// sigma, desvio padrão calculado internamente
        medidas[medC].prec = 0.042;// precisão do medidor
        medidas[medC].k=alimentador->noRaiz;
        medidas[medC].m=-1;
        medidas[medC].ramo=-1;

        k=alimentador->noRaiz;



        // coloca a medida no grafo
        ind = grafo[k].nmed;
        grafo[k].medidores[ind] = &medidas[medB];
        grafo[k].nmed++;

        //aloca o espaço na estrutura medidas
        medidas[medB].nvar = 3;
        medidas[medB].reguaH = (double *)malloc(medidas[medB].nvar * sizeof(double));
        medidas[medB].H = (double *)malloc(medidas[medB].nvar * sizeof(double));
        medidas[medC].nvar = 3;
        medidas[medC].reguaH = (double *)malloc(medidas[medC].nvar * sizeof(double));
        medidas[medC].H = (double *)malloc(medidas[medC].nvar * sizeof(double));
            //inicia os malores com zero
        for (j = 0; j < medidas[medB].nvar; j++)
        {
            medidas[medB].H[j] = 0;
        }

        //preenche a regua

        regua = (double)k;
        regua += 0.01;
        for (j = 0; j < 3; j++)
        {
            medidas[medB].reguaH[j] = -regua;
            regua += 0.1;
        }
        ind = grafo[k].nmed;
        grafo[k].medidores[ind] = &medidas[medC];
        grafo[k].nmed++;

        //aloca o espaço na estrutura medidas
        medidas[medC].nvar = 3;
        medidas[medC].reguaH = (double *)malloc(medidas[medC].nvar * sizeof(double));
        medidas[medC].H = (double *)malloc(medidas[medC].nvar * sizeof(double));
            //inicia os malores com zero
        for (j = 0; j < medidas[medC].nvar; j++)
        {
            medidas[medC].H[j] = 0;
        }

            //preenche a regua

        regua = (double)k;
        regua += 0.01;
        for (j = 0; j < 3; j++)
        {
            medidas[medC].reguaH[j] = -regua;
            regua += 0.1;
        }

    
    
}

void tratamento_referencia_improv(long int *ref_1, long int *ref_2, ALIMENTADOR *alimentador, double *regua, long int nvar,double ***H,GRAFO *grafo, long int numeroBarras,DRAM *ramos, DMED *medidas, long int *nmed,int* pseu)
{
    double **Href;
    int nref;
    int i,j,k;
    *pseu=0;
    
    atualiza_H(grafo, numeroBarras, ramos, medidas, *nmed);
    matrizDinamica(&Href,(int)*nmed,(int)nvar);

    mat_ig(H,*nmed,nvar,Href);
    nref=numref(Href,(int)*nmed,(int)nvar,1,1);
   
    if(nref>1)
    {
        (*nmed)=(*nmed)+2;
        add_pseudo(alimentador,grafo,ramos,medidas,*nmed);
        *pseu=1;
        nref=1;   
    }

    

    for (j = 0; j < nvar; j++)
    {
        if (regua[j] < 0)
        {
            k = (int)regua[j];
            k = -k;
            if (k == alimentador->noRaiz)
            {
                ref_1[0] = j;
                ref_2[0] = j + nref - 1;
                break;
            }
        }
    }
   

}

//Tratamento da referência (por alimentador)
void tratamento_referencia(long int *ref_1, long int *ref_2, ALIMENTADOR *alimentador, double *regua, long int nVariaveis)
{
    int i, j, k, n_ref = 3;

    //Inserir aqui fatoração da matriz H
    //atualizaH
    //fatoração e obtenção da Hdelta
    //número de referências

    //Inicialização com base nas tensões lidas no arquivo Vinicial.csv
    for (j = 0; j < nVariaveis; j++)
    {
        if (regua[j] < 0)
        {
            k = (int)regua[j];
            k = -k;
            if (k == alimentador->noRaiz)
            {
                ref_1[0] = j;
                ref_2[0] = j + n_ref - 1;
                break;
            }
        }
    }
}

//------------------------------------------------------------------------------
//
// ESTIMADOR WLS CONVENCIONAL VIA WLS
//
//------------------------------------------------------------------------------
void estimadorWLS(GRAFO *grafo, long int numeroBarras, DMED *medidas, long int **numeroMedidas, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos, double Sbase,char **argv)
{
    long int nmed, nvar,tratref;
    int i, j, k, r,nref,pseu;
    long int ref_1, ref_2;
    double *z = NULL, **h = NULL, ***H = NULL, **W = NULL, *x = NULL, *regua = NULL,*regua_o=NULL, aux = 0;
    FILE *reguaout;
    reguaout=fopen("regua.csv","w");


    if(argv[1]==NULL || atoi(argv[1])==1)
    {
        tratref=1;
    }
    else if(atoi(argv[1])==2)
    {
        tratref=2;
    }
    else if(atoi(argv[1])==3)
    {
        tratref=3;
    }
   
    
    
    printf("Estimador de Estado WLS Trifásico...\n");
    //--------------------------------------------------------------------------
    //Alocação de memória das variáveis do estimador de estado
    nmed = 0;
    for (i = 0; i < 9; i++)
    {
        for (j = 0; j < 8; j++)
        {
            nmed = nmed + numeroMedidas[i][j];
        }
    }
    nvar = 0;
    //printf("numero barras: %d\n", numeroBarras);
    for (i = 0; i < numeroBarras; i++)
    {
        switch (grafo[i].fases)
        {
        case 1:
            nvar += 2;
            break;
        case 2:
            nvar += 2;
            break;
        case 3:
            nvar += 2;
            break;
        case 4:
            nvar += 4;
            break;
        case 5:
            nvar += 4;
            break;
        case 6:
            nvar += 4;
            break;
        case 7:
            nvar += 6;
            break;
        }
    }

    if ((z = (double *)malloc((nmed+2) * sizeof(double))) == NULL)
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor z!!!!");
        exit(1);
    }
    if ((h = malloc((nmed+2) * sizeof(double *))) == NULL)
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor h!!!!");
        exit(1);
    }
    if ((x = (double *)malloc((nvar) * sizeof(double))) == NULL)
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor x!!!!");
        exit(1);
    }
    if ((regua = (double *)malloc((nvar) * sizeof(double))) == NULL)
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor regua!!!!");
        exit(1);
    }
    if ((regua_o = (double *)malloc((nvar) * sizeof(double))) == NULL)
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor regua!!!!");
        exit(1);
    }

    H = (double ***)malloc((nmed+2) * sizeof(double **));
    for (i = 0; i < nmed+2; i++)
    {
        H[i] = (double **)malloc(nvar * sizeof(double *));
        for (j = 0; j < nvar; j++)
        {
            H[i][j] = &aux;
        }
    }
    //--------------------------------------------------------------------------
    // Direcionamento dos ponteiros que compõem a régua, vetor h e matriz H
    //monta a regua (ordenaçao do vetor x) - V e teta conforme o grafo
    j = 0;
    for (i = 0; i < numeroBarras; i++)
    {
        aux = (double)grafo[i].idNo;
        aux += 0.01;
        switch (grafo[i].fases)
        {
        case 1:
            regua[j] = aux;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            break;
        case 2:
            regua[j] = aux + 0.1;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            break;
        case 3:
            regua[j] = aux + 0.2;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            break;
        case 4:
            regua[j] = aux;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            regua[j] = aux + 0.1;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            break;
        case 5:
            regua[j] = aux;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            regua[j] = aux + 0.2;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            break;
        case 6:
            regua[j] = aux + 0.1;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            regua[j] = aux + 0.2;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            break;
        case 7:
            regua[j] = aux;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            regua[j] = aux + 0.1;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            regua[j] = aux + 0.2;
            regua[j + (int)nvar / 2] = -regua[j];
            j++;
            break;
        }
    }
    for(i=0;i<(int)nvar / 2;i++)
    {
        regua_o[i] = regua[i];
        regua_o[i + (int)nvar / 2] = -regua_o[i];
    }
    aux = 0;

    for (i=0;i<nvar;i++)
    {
        fprintf(reguaout,"%f\n",regua_o[i]);
    }

    


    for (i = 0; i < nmed; i++)
    {
        h[i] = &medidas[i].h;
    }
    //Matriz H aponta para a estrutura de dados das medidas
    for (i = 0; i < nmed; i++)
    {
        for (j = 0; j < medidas[i].nvar; j++)
        {
            for (r = 0; r < nvar; r++)
            {
                if (cabs(medidas[i].reguaH[j] - regua[r]) < EPS)
                {
                    H[i][r] = &medidas[i].H[j];
                    break;
                }
            }
        }
    }

   
    incializa_vetor_x(grafo, numeroBarras, alimentadores, numeroAlimentadores, x, regua, nvar);

    if (tratref==3)
    {
        tratamento_referencia_improv(&ref_1, &ref_2, &alimentadores[0], regua, nvar,H,grafo,numeroBarras,ramos,medidas,&nmed,&pseu);
        if (pseu ==1)
        {
            for (i = 0; i < nmed; i++)
            {
                h[i] = &medidas[i].h;
            }
            //Matriz H aponta para a estrutura de dados das medidas
            for (i = 0; i < nmed; i++)
            {
                for (j = 0; j < medidas[i].nvar; j++)
                {
                    for (r = 0; r < nvar; r++)
                    {
                        if (cabs(medidas[i].reguaH[j] - regua[r]) < EPS)
                        {
                            H[i][r] = &medidas[i].H[j];
                            break;
                        }
                    }
                }
            }

        }
        nref=1;
    }
    else if (tratref==1||tratref==2)
    {
        tratamento_referencia(&ref_1,&ref_2,&alimentadores[0],regua,nvar);
        nref=3;
        if (tratref==2)
        {
            ref_2=ref_1;
            nref=1;    
        }
    }
    
    



    double **Href;
    FILE* JACOBPSEU;
    atualiza_H(grafo, numeroBarras, ramos, medidas, nmed);

    
    matrizDinamica(&Href,nmed,nvar);

    mat_ig(H,nmed,nvar,Href);

    JACOBPSEU=fopen("Hpseudoan.csv","w");
    fimprimirmat(JACOBPSEU,Href,nmed,nvar);

    fclose(JACOBPSEU);
    
    
    //tratamento_referencia(&ref_1, &ref_2, &alimentadores[0], regua, nvar);
    //--------------------------------------------------------------------------
    //Estimação de Estado
    monta_z(z, nmed, medidas);
    monta_W(NULL, nmed, medidas);
    //monta_W_cte(W,nmed,medidas);
    //monta_W_Ident(NULL,nmed,medidas);

    
    //incializa_vetor_x_perturbado(grafo, numeroBarras, alimentadores, numeroAlimentadores, x, regua, nvar);
    FILE *TWLS;
    TWLS=fopen("twls.csv","a");
    double tol = 1e-6;
    clock_t tic  = clock();
    
    (void) otimiza_Gauss_NewtonQR(z, h, H, grafo, numeroBarras, ramos, medidas, nvar, nmed, regua, x, tol, ref_1, ref_2,nref);

    clock_t toc = clock();

    double tempoWLS = (double)(toc - tic) / CLOCKS_PER_SEC;
    printf("\nEstimação WLS: %lf", tempoWLS);
    
    ///fprintf(TWLS,"%.2e\n",tempoWLS);

    exportaCasoReferencia(grafo, numeroBarras, Sbase);
    exportaRamoscorrentes(grafo, numeroBarras, Sbase);
    exportaEstado(grafo, regua, nvar);
    imprimeEstado(grafo, numeroBarras);
    //imprimeInjecoes(grafo, numeroBarras);

    FILE *residuo;
    residuo = fopen("residuo.txt", "wt");
    for (i = 0; i < nmed; i++)
    {
        fprintf(residuo, "%.10lf\n", z[i] - medidas[i].h);
    }
    fclose(residuo);
    fclose(TWLS);
    free(z);
    free(h);
    free(x);
    free(regua);
    for (i = 0; i < nmed; i++)
        free(H[i]);
    free(H);
}

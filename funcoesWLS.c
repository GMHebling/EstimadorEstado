#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <complex.h>

#include "data_structures.h"
#include "funcoesWLS.h"
#include "funcoesTopologia.h"
#include "funcoesCalculoEletrico.h"
#include "funcoesMatematicas.h"

//Imprime na tela o estado atual da rede (tensões complexas nodais em cada fase)
void imprimeEstado(GRAFO *grafo,  long int numeroBarras){
    int i;
    
    printf("\nTensoes Nodais (p.u.):\n");
    for(i=0; i<numeroBarras; i++){ 
        //Retangulares
        //printf("%d\tVa: %.5lf + j%.5lf\tVb: %.5lf + j%.5lf\tVc: %.5lf + j%.5lf\n",grafo[i].barra->ID,__real__ grafo[i].V[0],__imag__ grafo[i].V[0],__real__ grafo[i].V[1],__imag__ grafo[i].V[1],__real__ grafo[i].V[2],__imag__ grafo[i].V[2]);
        //Polares
        switch (grafo[i].fases){
            case 1:
                printf("%d\tVa: %.5lf | %.3lf \tVb:    -    |    -   \tVc:    -    |    -   \n",grafo[i].barra->ID,cabs(grafo[i].V[0]),carg(grafo[i].V[0])*180/PI);
                break;
            case 2:
                printf("%d\tVa:    -    |    -    \tVb: %.5lf | %.3lf\tVc:    -    |    -   \n",grafo[i].barra->ID,cabs(grafo[i].V[1]),carg(grafo[i].V[1])*180/PI);
                break;
            case 3:
                printf("%d\tVa:    -    |    -    \tVb:    -    |    -   \tVc: %.5lf | %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[0]),cabs(grafo[i].V[2]),carg(grafo[i].V[2])*180/PI);
                break;
            case 4:
                printf("%d\tVa: %.5lf | %.3lf \tVb: %.5lf | %.3lf\tVc:    -    |    -   \n",grafo[i].barra->ID,cabs(grafo[i].V[0]),carg(grafo[i].V[0])*180/PI,cabs(grafo[i].V[1]),carg(grafo[i].V[1])*180/PI);
                break;
            case 5:
                printf("%d\tVa: %.5lf | %.3lf \tVb:    -    |    -   \tVc: %.5lf | %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[0]),carg(grafo[i].V[0])*180/PI,cabs(grafo[i].V[2]),carg(grafo[i].V[2])*180/PI);
                break;
            case 6:
                printf("%d\tVa:    -    |    -    \tVb: %.5lf | %.3lf\tVc: %.5lf | %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[1]),carg(grafo[i].V[1])*180/PI,cabs(grafo[i].V[2]),carg(grafo[i].V[2])*180/PI);
                break;
            case 7:
                printf("%d\tVa: %.5lf | %.3lf \tVb: %.5lf | %.3lf\tVc: %.5lf | %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[0]),carg(grafo[i].V[0])*180/PI,cabs(grafo[i].V[1]),carg(grafo[i].V[1])*180/PI,cabs(grafo[i].V[2]),carg(grafo[i].V[2])*180/PI);
                break;    
        }
    }
}

//Salva arquivo com o estado atual da rede (vetor x) conforme régua
void exportaEstado(GRAFO *grafo,  double *regua, long int nvar){
    int i, k, fase;
    FILE *arqout;
    
    arqout = fopen("state.txt","w+");
    for(i=0;i<nvar;i++){
        if (regua[i]>0){
            k = (int) regua[i];
            fase = (int) cabs((regua[i] - k)*10);                    
            fprintf(arqout,"%.2f\t%.15f\n",regua[i],cabs(grafo[k].V[fase]));
        }
        else{
            k = (int) regua[i];
            fase = (int)cabs((regua[i] - k)*10);                    
            k = -k;
            fprintf(arqout,"%.2f\t%.15f\n",regua[i],carg(grafo[k].V[fase]));
        }
    }
    fclose(arqout);
}

void saidaEstado(GRAFO *grafo, long int numeroBarras, int it, double tempoIt, double nFx, double nGx){
    int i,j,k;
    
    FILE *arquivo;
    arquivo = fopen("resultadoEstado.txt","wt");
    
    fprintf(arquivo, "|Dx|_inf =  %.7lf \t |Grad|_inf =  %.7lf \n",nFx,nGx);
    fprintf(arquivo,"\n\n Convergência em %d iteracoes e tempo: %.4lf",it,tempoIt);
    
    fprintf(arquivo,"\nTensoes Nodais: Fase-Terra\n");
    for(i=0; i<numeroBarras; i++){ 
        fprintf(arquivo,"Vbase:%.3lf\t",grafo[i].Vbase);
        switch (grafo[i].fases){
            case 1:
                fprintf(arquivo,"A\t%d\tVa: %.5lf | %.3lf \tVb:    -    |    -   \tVc:    -    |    -   \n",grafo[i].barra->ID,cabs(grafo[i].V[0]),carg(grafo[i].V[0])*180/PI);
                break;
            case 2:
                fprintf(arquivo,"B\t%d\tVa:    -    |    -    \tVb: %.5lf | %.3lf\tVc:    -    |    -   \n",grafo[i].barra->ID,cabs(grafo[i].V[1]),carg(grafo[i].V[1])*180/PI);
                break;
            case 3:
                fprintf(arquivo,"C\t%d\tVa:    -    |    -    \tVb:    -    |    -   \tVc: %.5lf | %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[2]),carg(grafo[i].V[2])*180/PI);
                break;
            case 4:
                fprintf(arquivo,"AB\t%d\tVa: %.5lf | %.3lf \tVb: %.5lf | %.3lf\tVc:    -    |    -   \n",grafo[i].barra->ID,cabs(grafo[i].V[0]),carg(grafo[i].V[0])*180/PI,cabs(grafo[i].V[1]),carg(grafo[i].V[1])*180/PI);
                break;
            case 5:
                fprintf(arquivo,"CA\t%d\tVa: %.5lf | %.3lf \tVb:    -    |    -   \tVc: %.5lf | %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[0]),carg(grafo[i].V[0])*180/PI,cabs(grafo[i].V[2]),carg(grafo[i].V[2])*180/PI);
                break;
            case 6:
                fprintf(arquivo,"BC\t%d\tVa:    -    |    -    \tVb: %.5lf | %.3lf\tVc: %.5lf | %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[1]),carg(grafo[i].V[1])*180/PI,cabs(grafo[i].V[2]),carg(grafo[i].V[2])*180/PI);
                break;
            case 7:
                fprintf(arquivo,"ABC\t%d\tVa: %.5lf | %.3lf \tVb: %.5lf | %.3lf\tVc: %.5lf | %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[0]),carg(grafo[i].V[0])*180/PI,cabs(grafo[i].V[1]),carg(grafo[i].V[1])*180/PI,cabs(grafo[i].V[2]),carg(grafo[i].V[2])*180/PI);
                break;    
        }
    }
    fprintf(arquivo,"\nTensoes Nodais: Fase-Fase\n");
    for(i=0; i<numeroBarras; i++){
        fprintf(arquivo,"Vbase:%.3lf\t",grafo[i].Vbase);
        switch (grafo[i].fases){
            case 1:
                fprintf(arquivo,"A\t%d\tVan: %.5lf | %.3lf \tVbn:    -    |    -   \tVcn:    -    |    -   \n",grafo[i].barra->ID,cabs(grafo[i].V[0]),carg(grafo[i].V[0])*180/PI);
                break;
            case 2:
                fprintf(arquivo,"B\t%d\tVan:    -    |    -    \tVbn: %.5lf | %.3lf\tVcn:    -    |    -   \n",grafo[i].barra->ID,cabs(grafo[i].V[1]),carg(grafo[i].V[1])*180/PI);
                break;
            case 3:
                fprintf(arquivo,"C\t%d\tVan:    -    |    -    \tVbn:    -    |    -   \tVcn: %.5lf | %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[2]),carg(grafo[i].V[2])*180/PI);
                break;
            case 4:
                fprintf(arquivo,"AB\t%d\tVab: %.5lf | %.3lf \tVbc:   -    |   -  \tVca:    -    |    -   \n",grafo[i].barra->ID,cabs(grafo[i].V[0] - grafo[i].V[1]),carg(grafo[i].V[0] - grafo[i].V[1])*180/PI);
                break;
            case 5:
                fprintf(arquivo,"CA\t%d\tVab:   -   |   -  \tVbc:    -    |    -   \tVca: %.5lf | %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[2] - grafo[i].V[0]),carg(grafo[i].V[2] - grafo[i].V[0])*180/PI);
                break;
            case 6:
                fprintf(arquivo,"BC\t%d\tVab:    -    |    -    \tVbc: %.5lf | %.3lf\tVca:   -    |   -  \n",grafo[i].barra->ID,cabs(grafo[i].V[1] - grafo[i].V[2]),carg(grafo[i].V[1]  -grafo[i].V[2])*180/PI);
                break;
            case 7:
                fprintf(arquivo,"ABC\t%d\tVab: %.5lf | %.3lf \tVbc: %.5lf | %.3lf\tVca: %.5lf | %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[0]-grafo[i].V[1]),carg(grafo[i].V[0]- grafo[i].V[1])*180/PI,cabs(grafo[i].V[1] - grafo[i].V[2]),carg(grafo[i].V[1] - grafo[i].V[2])*180/PI,cabs(grafo[i].V[2] - grafo[i].V[0]),carg(grafo[i].V[2] - grafo[i].V[0])*180/PI);
                break;    
        }
    }
    fclose(arquivo);    
}

//Exporta arquivo de texto com as informações da distribuição a Priori SCADA
void exportaCasoReferencia(GRAFO *grafo, long int numeroBarras, double Sbase){
    int i,j,k;
    
    FILE *arquivo;
    arquivo = fopen("referencia.txt","wt");
    
    //----------------------------------------------------------------------
    //Fluxo de potência ativa em kW
    for(i=0;i<numeroBarras;i++){    
        //Percorre os ramos adjacentes
        for(k=0;k<grafo[i].numeroAdjacentes;k++){
            
            switch (grafo[i].adjacentes[k].ramo->fases){
                case 1:
                    fprintf(arquivo,"1,0,%d,%d,1,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__real__ grafo[i].adjacentes[k].S[0]) * Sbase);
                    break;
                case 2:
                    fprintf(arquivo,"1,0,%d,%d,2,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__real__ grafo[i].adjacentes[k].S[1]) * Sbase);
                    break;
                case 3:
                    fprintf(arquivo,"1,0,%d,%d,3,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__real__ grafo[i].adjacentes[k].S[2]) * Sbase);
                    break;
                case 4:
                    fprintf(arquivo,"1,0,%d,%d,1,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__real__ grafo[i].adjacentes[k].S[0]) * Sbase);
                    fprintf(arquivo,"1,0,%d,%d,2,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__real__ grafo[i].adjacentes[k].S[1]) * Sbase);
                    break;
                case 5:
                    fprintf(arquivo,"1,0,%d,%d,1,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__real__ grafo[i].adjacentes[k].S[0]) * Sbase);
                    fprintf(arquivo,"1,0,%d,%d,3,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__real__ grafo[i].adjacentes[k].S[2]) * Sbase);
                    break;
                case 6:
                    fprintf(arquivo,"1,0,%d,%d,2,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__real__ grafo[i].adjacentes[k].S[1]) * Sbase);
                    fprintf(arquivo,"1,0,%d,%d,3,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__real__ grafo[i].adjacentes[k].S[2]) * Sbase);
                    break;
                case 7:
                    fprintf(arquivo,"1,0,%d,%d,1,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__real__ grafo[i].adjacentes[k].S[0]) * Sbase);
                    fprintf(arquivo,"1,0,%d,%d,2,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__real__ grafo[i].adjacentes[k].S[1]) * Sbase);
                    fprintf(arquivo,"1,0,%d,%d,3,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__real__ grafo[i].adjacentes[k].S[2]) * Sbase);
                    break;
            }
        }        
    }
    //----------------------------------------------------------------------
    //Fluxo de potência reativa em kVAr
    for(i=0;i<numeroBarras;i++){    
        //Percorre os ramos adjacentes
        for(k=0;k<grafo[i].numeroAdjacentes;k++){
            
            switch (grafo[i].adjacentes[k].ramo->fases){
                case 1:
                    fprintf(arquivo,"1,1,%d,%d,1,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__imag__ grafo[i].adjacentes[k].S[0]) * Sbase);
                    break;
                case 2:
                    fprintf(arquivo,"1,1,%d,%d,2,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__imag__ grafo[i].adjacentes[k].S[1]) * Sbase);
                    break;
                case 3:
                    fprintf(arquivo,"1,1,%d,%d,3,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__imag__ grafo[i].adjacentes[k].S[2]) * Sbase);
                    break;
                case 4:
                    fprintf(arquivo,"1,1,%d,%d,1,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__imag__ grafo[i].adjacentes[k].S[0]) * Sbase);
                    fprintf(arquivo,"1,1,%d,%d,2,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__imag__ grafo[i].adjacentes[k].S[1]) * Sbase);
                    break;
                case 5:
                    fprintf(arquivo,"1,1,%d,%d,1,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__imag__ grafo[i].adjacentes[k].S[0]) * Sbase);
                    fprintf(arquivo,"1,1,%d,%d,3,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__imag__ grafo[i].adjacentes[k].S[2]) * Sbase);
                    break;
                case 6:
                    fprintf(arquivo,"1,1,%d,%d,2,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__imag__ grafo[i].adjacentes[k].S[1]) * Sbase);
                    fprintf(arquivo,"1,1,%d,%d,3,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__imag__ grafo[i].adjacentes[k].S[2]) * Sbase);
                    break;
                case 7:
                    fprintf(arquivo,"1,1,%d,%d,1,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__imag__ grafo[i].adjacentes[k].S[0]) * Sbase);
                    fprintf(arquivo,"1,1,%d,%d,2,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__imag__ grafo[i].adjacentes[k].S[1]) * Sbase);
                    fprintf(arquivo,"1,1,%d,%d,3,%.8lf,0.020000\n",grafo[i].barra->ID,grafo[grafo[i].adjacentes[k].idNo].barra->ID,(__imag__ grafo[i].adjacentes[k].S[2]) * Sbase);
                    break;
            }
        }        
    }
    
    
    
    //----------------------------------------------------------------------
    //Injeção de potência ativa em kW
    for(i=0;i<numeroBarras;i++){
        switch (grafo[i].fases){
            case 1:
                fprintf(arquivo,"1,2,%d,-1,1,%.8lf,0.020000\n",grafo[i].barra->ID,(__real__ grafo[i].S[0]) * Sbase);
                break;
            case 2:
                fprintf(arquivo,"1,2,%d,-1,2,%.8lf,0.020000\n",grafo[i].barra->ID,(__real__ grafo[i].S[1]) * Sbase);
                break;
            case 3:
                fprintf(arquivo,"1,2,%d,-1,3,%.8lf,0.020000\n",grafo[i].barra->ID,(__real__ grafo[i].S[2]) * Sbase);
                break;
            case 4:
                fprintf(arquivo,"1,2,%d,-1,1,%.8lf,0.020000\n",grafo[i].barra->ID,(__real__ grafo[i].S[0]) * Sbase);
                fprintf(arquivo,"1,2,%d,-1,2,%.8lf,0.020000\n",grafo[i].barra->ID,(__real__ grafo[i].S[1]) * Sbase);
                break;
            case 5:
                fprintf(arquivo,"1,2,%d,-1,1,%.8lf,0.020000\n",grafo[i].barra->ID,(__real__ grafo[i].S[0]) * Sbase);
                fprintf(arquivo,"1,2,%d,-1,3,%.8lf,0.020000\n",grafo[i].barra->ID,(__real__ grafo[i].S[2]) * Sbase);
                break;
            case 6:
                fprintf(arquivo,"1,2,%d,-1,2,%.8lf,0.020000\n",grafo[i].barra->ID,(__real__ grafo[i].S[1]) * Sbase);
                fprintf(arquivo,"1,2,%d,-1,3,%.8lf,0.020000\n",grafo[i].barra->ID,(__real__ grafo[i].S[2]) * Sbase);
                break;
            case 7:
                fprintf(arquivo,"1,2,%d,-1,1,%.8lf,0.020000\n",grafo[i].barra->ID,(__real__ grafo[i].S[0]) * Sbase);
                fprintf(arquivo,"1,2,%d,-1,2,%.8lf,0.020000\n",grafo[i].barra->ID,(__real__ grafo[i].S[1]) * Sbase);
                fprintf(arquivo,"1,2,%d,-1,3,%.8lf,0.020000\n",grafo[i].barra->ID,(__real__ grafo[i].S[2]) * Sbase);
                break;
        }
    }
    //----------------------------------------------------------------------
    //Injeção de potência ativa em kVAr
    for(i=0;i<numeroBarras;i++){
        switch (grafo[i].fases){
            case 1:
                fprintf(arquivo,"1,3,%d,-1,1,%.8lf,0.020000\n",grafo[i].barra->ID,(__imag__ grafo[i].S[0]) * Sbase);
                break;
            case 2:
                fprintf(arquivo,"1,3,%d,-1,2,%.8lf,0.020000\n",grafo[i].barra->ID,(__imag__ grafo[i].S[1]) * Sbase);
                break;
            case 3:
                fprintf(arquivo,"1,3,%d,-1,3,%.8lf,0.020000\n",grafo[i].barra->ID,(__imag__ grafo[i].S[2]) * Sbase);
                break;
            case 4:
                fprintf(arquivo,"1,3,%d,-1,1,%.8lf,0.020000\n",grafo[i].barra->ID,(__imag__ grafo[i].S[0]) * Sbase);
                fprintf(arquivo,"1,3,%d,-1,2,%.8lf,0.020000\n",grafo[i].barra->ID,(__imag__ grafo[i].S[1]) * Sbase);
                break;
            case 5:
                fprintf(arquivo,"1,3,%d,-1,1,%.8lf,0.020000\n",grafo[i].barra->ID,(__imag__ grafo[i].S[0]) * Sbase);
                fprintf(arquivo,"1,3,%d,-1,3,%.8lf,0.020000\n",grafo[i].barra->ID,(__imag__ grafo[i].S[2]) * Sbase);
                break;
            case 6:
                fprintf(arquivo,"1,3,%d,-1,2,%.8lf,0.020000\n",grafo[i].barra->ID,(__imag__ grafo[i].S[1]) * Sbase);
                fprintf(arquivo,"1,3,%d,-1,3,%.8lf,0.020000\n",grafo[i].barra->ID,(__imag__ grafo[i].S[2]) * Sbase);
                break;
            case 7:
                fprintf(arquivo,"1,3,%d,-1,1,%.8lf,0.020000\n",grafo[i].barra->ID,(__imag__ grafo[i].S[0]) * Sbase);
                fprintf(arquivo,"1,3,%d,-1,2,%.8lf,0.020000\n",grafo[i].barra->ID,(__imag__ grafo[i].S[1]) * Sbase);
                fprintf(arquivo,"1,3,%d,-1,3,%.8lf,0.020000\n",grafo[i].barra->ID,(__imag__ grafo[i].S[2]) * Sbase);
                break;
        }
    }
    //----------------------------------------------------------------------
    //Magnitudes de tensão em kV
    for(i=0;i<numeroBarras;i++){
        switch (grafo[i].fases){
            case 1:
                fprintf(arquivo,"1,4,%d,-1,1,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*cabs(grafo[i].V[0]));
                break;
            case 2:
                fprintf(arquivo,"1,4,%d,-1,2,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*cabs(grafo[i].V[1]));
                break;
            case 3:
                fprintf(arquivo,"1,4,%d,-1,3,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*cabs(grafo[i].V[2]));
                break;
            case 4:
                fprintf(arquivo,"1,4,%d,-1,1,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*cabs(grafo[i].V[0]));
                fprintf(arquivo,"1,4,%d,-1,2,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*cabs(grafo[i].V[1]));
                break;
            case 5:
                fprintf(arquivo,"1,4,%d,-1,1,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*cabs(grafo[i].V[0]));
                fprintf(arquivo,"1,4,%d,-1,3,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*cabs(grafo[i].V[2]));
                break;
            case 6:
                fprintf(arquivo,"1,4,%d,-1,2,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*cabs(grafo[i].V[1]));
                fprintf(arquivo,"1,4,%d,-1,3,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*cabs(grafo[i].V[2]));
                break;
            case 7:
                fprintf(arquivo,"1,4,%d,-1,1,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*cabs(grafo[i].V[0]));
                fprintf(arquivo,"1,4,%d,-1,2,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*cabs(grafo[i].V[1]));
                fprintf(arquivo,"1,4,%d,-1,3,%.8lf,0.010000\n",grafo[i].barra->ID,grafo[i].Vbase/1000*cabs(grafo[i].V[2]));
                break;
        }
    }
//    //----------------------------------------------------------------------
//    //Ângulo de tensão em graus
    for(i=0;i<numeroBarras;i++){
        switch (grafo[i].fases){
            case 1:
                fprintf(arquivo,"1,5,%d,-1,1,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*carg(grafo[i].V[0]));
                break;
            case 2:
                fprintf(arquivo,"1,5,%d,-1,2,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*carg(grafo[i].V[1]));
                break;
            case 3:
                fprintf(arquivo,"1,5,%d,-1,3,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*carg(grafo[i].V[2]));
                break;
            case 4:
                fprintf(arquivo,"1,5,%d,-1,1,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*carg(grafo[i].V[0]));
                fprintf(arquivo,"1,5,%d,-1,2,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*carg(grafo[i].V[1]));
                break;
            case 5:
                fprintf(arquivo,"1,5,%d,-1,1,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*carg(grafo[i].V[0]));
                fprintf(arquivo,"1,5,%d,-1,3,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*carg(grafo[i].V[2]));
                break;
            case 6:
                fprintf(arquivo,"1,5,%d,-1,2,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*carg(grafo[i].V[1]));
                fprintf(arquivo,"1,5,%d,-1,3,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*carg(grafo[i].V[2]));
                break;
            case 7:
                fprintf(arquivo,"1,5,%d,-1,1,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*carg(grafo[i].V[0]));
                fprintf(arquivo,"1,5,%d,-1,2,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*carg(grafo[i].V[1]));
                fprintf(arquivo,"1,5,%d,-1,3,%.8lf,0.001000\n",grafo[i].barra->ID,180/PI*carg(grafo[i].V[2]));
                break;
        }
    }
    fclose(arquivo);
}

//Exporta arquivo de texto com as informações da distribuição a Priori SCADA
void exportaPrioriQR(GRAFO *grafo,double *regua, DMED *medidas, long int nvar, long int nmed, long int polar){
    int i,j,k,r, fase;
    double **R = NULL, **H = NULL;
    
    R = aloca_matriz(nmed,nvar);
    H = aloca_matriz(nmed,nvar);
    
    //Monta Matriz W^1/2 H
    for(i=0;i<nmed;i++){
        for(j=0;j<medidas[i].nvar;j++){
            for(r = 0;r<nvar;r++){
                if (cabs(medidas[i].reguaH[j]-regua[r]) < EPS){
                    H[i][r] = medidas[i].H[j]/medidas[i].sigma;
                    break;
                }
            }
        }
    }
    QRfactorization(H, nmed, nvar, R);
    
    //Exporta informações para estimação a posteriori
    FILE *priori;
    priori = fopen("prioriR.txt","w");
    
    fprintf(priori,"%d\t%d",nvar,nmed);
    if (polar == 1){
        for(i=0; i<nvar; i++){ 
            if (regua[i]>0){
                k = (int) regua[i];
                fase = (int) cabs((regua[i] - k)*10);                    
                fprintf(priori,"\n%.14lf",cabs(grafo[k].V[fase]));
            }
            else{
                k = (int) regua[i];
                fase = (int)cabs((regua[i] - k)*10);                    
                k = -k;
                fprintf(priori,"\n%.14lf",carg(grafo[k].V[fase]));
            }
        }
//        for(i=0;i<nvar;i++){
//            fprintf(priori,"\n");
//            for(j=i;j<nvar;j++){
//                fprintf(priori,"%.14f\t",R[i][j]);
//            }
//        }
        fprintf(priori,"\n");
        for(i=0;i<nvar;i++){
            for(j=i;j<nvar;j++){
                if (fabs(R[i][j]) >= 0.0000001) fprintf(priori,"%d\t%d\t%.14f\n",i,j,R[i][j]);
            }
        }
    }
    else{
        for(i=0; i<nvar; i++){ 
            if (regua[i]>0){
                k = (int) regua[i];
                fase = (int) cabs((regua[i] - k)*10);                    
                fprintf(priori,"\n%.14lf",__real__ grafo[k].V[fase]);
            }
            else{
                k = (int) regua[i];
                fase = (int) cabs((regua[i] - k)*10);                    
                k = -k;
                fprintf(priori,"\n%.14lf",__imag__ grafo[k].V[fase]);
            }
        }
        //Montar matriz G em coordenadas retangulares ou fazer rotação da matriz - Motada a matriz H em retangulares
//        for(i=0;i<nvar;i++){
//            fprintf(priori,"\n");
//            for(j=i;j<nvar;j++){
//                fprintf(priori,"%.14f\t",R[i][j]);
//            }
//        }
        fprintf(priori,"\n");
        for(i=0;i<nvar;i++){
            for(j=i;j<nvar;j++){
                if (fabs(R[i][j]) >= 0.0000001) fprintf(priori,"%d\t%d\t%.14f\n",i,j,R[i][j]);
            }
        }
        
    }
    
    
    free(H);
    free(R);
    fclose(priori);
}

//------------------------------------------------------------------------------
//
// FUNÇÕES DO ESTIMADOR WLS 
//
//------------------------------------------------------------------------------
//Função monta vetor z
void monta_z(double *z, long int nmed, DMED *medidas){
    int i;
    
    for(i=0;i<nmed;i++){
        z[i] = medidas[i].zmed;
    }    
}

void monta_z_comVirtuais(double *z, long int nmed, long int nvir, DMED *medidas, DMED *virtuais){
    int i;
    for (i=0; i<nmed; i++){
        z[i] = medidas[i].zmed;
    }
    for (i=0; i<nvir; i++){
        z[i+nmed] = virtuais[i].zmed;
    }
}

//Função monta matriz W
void monta_W(double **W, long int nmed, DMED *medidas){
    int i,j;
    double prec, menorSigma = 1000000;
    double fundoEscala, auxFE;
    
    //Matriz W diagonal - inverso da variância
    for(i=0;i<nmed;i++){
        auxFE = round(medidas[i].zmed * 1.25 * 1000);
        fundoEscala = auxFE/1000;
        switch (medidas[i].tipo){
            case 0:
            case 1:
//                medidas[i].sigma = 0.002;
//                break;
            case 2:    
            case 3:
                //W[i][i] = 40000;
//                medidas[i].sigma = 0.010;
                prec = medidas[i].prec; //0.02; //5% para SCADA de potência
                medidas[i].sigma = 0.33333*prec*cabs(medidas[i].zmed);
//                medidas[i].sigma = 0.33333*prec*cabs(fundoEscala);
                break;
            case 4:
            case 6:
                //W[i][i] = 40000;
//                medidas[i].sigma = 0.001;
                prec =  medidas[i].prec; //0.01; //1% para magnitude de tensão ou corrente
                medidas[i].sigma = 0.33333*prec*cabs(medidas[i].zmed);
//                medidas[i].sigma = 0.33333*prec*cabs(fundoEscala);
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
//                medidas[i].sigma = 0.0001;
                prec = medidas[i].prec; //0.001; //0.1% para PMUs retangulares
                medidas[i].sigma = 0.33333*prec*cabs(medidas[i].zmed);
                break;
        }
        //medidas[i].sigma = 0.33333*prec*cabs(medidas[i].zmed); //Ponderação de acordo com o valor medido (Fórmula B. Pal)
        if (cabs(medidas[i].zmed) > 0.00001){
            if (medidas[i].sigma < menorSigma) menorSigma = medidas[i].sigma;
        }
        medidas[i].sigma = medidas[i].sigma; 
//        W[i][i] = 1/(pow(medidas[i].sigma,2)); 
    }
    for(i=0;i<nmed;i++){ //Tratamento da medida virtual e medidas proximo de zero
        if (cabs(medidas[i].zmed) < 0.00001){
            medidas[i].sigma = menorSigma; 
            //medidas[i].sigma = 0.01*menorSigma; 
//            W[i][i] = 1/(pow(medidas[i].sigma,2));
        }
    }
    
}

//Função monta matriz W
void monta_W_cte(double **W, long int nmed, DMED *medidas){
    int i,j;
    double prec, menorSigma = 1000000;
    double fundoEscala, auxFE;
    
    //Matriz W diagonal - inverso da variância
    for(i=0;i<nmed;i++){
        switch (medidas[i].tipo){
            case 0:
            case 1:
            case 2:    
            case 3:
                if(medidas[i].prec == 0.3) medidas[i].sigma = 0.02; //pseudo medida (20 kVA)
                if(medidas[i].prec == 0.02) medidas[i].sigma = 0.002; //SCADA (2 kVA)
                if(medidas[i].prec == 0.001) medidas[i].sigma = 0.0001; //virtual (0.1 kVA)
                    
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
        if (cabs(medidas[i].zmed) > 0.00001){
            if (medidas[i].sigma < menorSigma) menorSigma = medidas[i].sigma;
        }
        //W[i][i] = 1/(pow(medidas[i].sigma,2));
    }
    for(i=0;i<nmed;i++){ //Tratamento da medida virtual e medidas proximo de zero
        if (cabs(medidas[i].zmed) < 0.00001){
            medidas[i].sigma = menorSigma; 
            //W[i][i] = 1/(pow(medidas[i].sigma,2));
        }
    }
    
}

void monta_W_Ident(double **W, long int nmed, DMED *medidas){
    int i,j;
    
    //Matriz W diagonal - inverso da variância
    for(i=0;i<nmed;i++){
        medidas[i].sigma = 1;
//        W[i][i] = 1;
    }    
}

//Função inicialização do vetor x
void incializa_vetor_x(GRAFO *grafo, long int numeroBarras, ALIMENTADOR *alimentadores, long int numeroAlimentadores, double *x, double *regua, long int nVariaveis){ 
    int i, k, fase;
    BOOL visitado[numeroBarras];
    __complex__ double V0[3], **Yaux;
    
    Yaux = c_matAloca(3);
    
    //Flat start trifásico (Va = Vb = Vc = 1p.u.  Ta = 0  Tb = -120  Tc = 120) - com busca em profundidade para atualizar taps iniciais
    for(i=0; i<numeroBarras; i++){ 
        visitado[i] = false;
    }
    for(i=0; i<numeroAlimentadores; i++)
    {
        //Tensão Inicial da subestação
        V0[0] = grafo[alimentadores[i].noRaiz].barra->Vinicial[0];//1.0*(cos(0) + I*sin(0));
        V0[1] = grafo[alimentadores[i].noRaiz].barra->Vinicial[1];//1.0*(cos(-120*PI/180) + I*sin(-120*PI/180));
        V0[2] = grafo[alimentadores[i].noRaiz].barra->Vinicial[2];//1.0*(cos(120*PI/180) + I*sin(120*PI/180));
        
        FILABARRAS *barraAtual = &alimentadores[i].rnp[0];
        
        int de = barraAtual->idNo;
        grafo[de].V[0] = V0[0];
        grafo[de].V[1] = V0[1];
        grafo[de].V[2] = V0[2];
//    }
//    for(i=0; i<numeroBarras; i++){
//        long int idAlim = grafo[i].idAlim;
//        grafo[i].V[0] = grafo[alimentadores[idAlim].noRaiz].V[0];
//        grafo[i].V[1] = grafo[alimentadores[idAlim].noRaiz].V[1];
//        grafo[i].V[2] = grafo[alimentadores[idAlim].noRaiz].V[2];
//        
//        if (grafo[i].Vbase <= 480){
//            grafo[i].V[0] = cabs(grafo[alimentadores[idAlim].noRaiz].V[0])*(cos(-30*PI/180) + I*sin(-30*PI/180));
//            grafo[i].V[1] = cabs(grafo[alimentadores[idAlim].noRaiz].V[1])*(cos(-150*PI/180) + I*sin(-150*PI/180));
//            grafo[i].V[2] = cabs(grafo[alimentadores[idAlim].noRaiz].V[2])*(cos(90*PI/180) + I*sin(90*PI/180));
//            
////            grafo[para].V[0] = cabs(grafo[de].V[0])*(cos(-30*PI/180) + I*sin(-30*PI/180));
////            grafo[para].V[1] = cabs(grafo[de].V[1])*(cos(-150*PI/180) + I*sin(-150*PI/180));
////            grafo[para].V[2g] = cabs(grafo[de].V[2])*(cos(90*PI/180) + I*sin(90*PI/180)); 
//        }
//            
//    }
        
        while(barraAtual != NULL)
        {
            de = barraAtual->idNo;
            int n_adj = grafo[de].numeroAdjacentes;
            for(k=0;k< n_adj;k++){
                int para = grafo[de].adjacentes[k].idNo;
                if (visitado[para] == false){ 
                    if ((grafo[de].adjacentes[k].tipo == 1)){ //Atualiza o V0 para trafo visto a ligação e tap
                        grafo[para].V[0] = grafo[de].V[0];
                        grafo[para].V[1] = grafo[de].V[1];
                        grafo[para].V[2] = grafo[de].V[2];
                        
                        if ((grafo[de].adjacentes[k].ramo->trafo.lig_pri == 1) && (grafo[de].adjacentes[k].ramo->trafo.lig_sec == 2)){                                
                            grafo[para].V[0] = cabs(grafo[de].V[0])*(cos(-30*PI/180) + I*sin(-30*PI/180));
                            grafo[para].V[1] = cabs(grafo[de].V[1])*(cos(-150*PI/180) + I*sin(-150*PI/180));
                            grafo[para].V[2] = cabs(grafo[de].V[2])*(cos(90*PI/180) + I*sin(90*PI/180));                                    
                        }
                        else if ((grafo[de].adjacentes[k].ramo->trafo.lig_pri == 3) && (grafo[de].adjacentes[k].ramo->trafo.lig_sec == 2)){                                
                            grafo[para].V[0] = cabs(grafo[de].V[0])*(cos(-30*PI/180) + I*sin(-30*PI/180));
                            grafo[para].V[1] = cabs(grafo[de].V[1])*(cos(-150*PI/180) + I*sin(-150*PI/180));
                            grafo[para].V[2] = cabs(grafo[de].V[2])*(cos(90*PI/180) + I*sin(90*PI/180));                                    
                        }
                        else if ((grafo[de].adjacentes[k].ramo->trafo.lig_pri == 2) && (grafo[de].adjacentes[k].ramo->trafo.lig_sec == 1)){
                            if (grafo[de].adjacentes[k].ramo->k == de){
                                grafo[para].V[0] = cabs(grafo[de].V[0])*(cos(-30*PI/180) + I*sin(-30*PI/180));
                                grafo[para].V[1] = cabs(grafo[de].V[1])*(cos(-150*PI/180) + I*sin(-150*PI/180));
                                grafo[para].V[2] = cabs(grafo[de].V[2])*(cos(90*PI/180) + I*sin(90*PI/180));
                            }
                            else{
                                grafo[para].V[0] = cabs(grafo[de].V[0])*(cos(0*PI/180) + I*sin(0*PI/180));
                                grafo[para].V[1] = cabs(grafo[de].V[1])*(cos(-120*PI/180) + I*sin(-120*PI/180));
                                grafo[para].V[2] = cabs(grafo[de].V[2])*(cos(120*PI/180) + I*sin(120*PI/180));
                            }
                        }            

                    }
                    else if (grafo[de].adjacentes[k].tipo == 2){ //Para o caso de regulador de tensão
                        grafo[para].V[0] = grafo[de].V[0]*grafo[de].adjacentes[k].ramo->tap_pri[0]*grafo[de].adjacentes[k].ramo->tap_sec[0];
                        grafo[para].V[1] = grafo[de].V[1]*grafo[de].adjacentes[k].ramo->tap_pri[1]*grafo[de].adjacentes[k].ramo->tap_sec[1];
                        grafo[para].V[2] = grafo[de].V[2]*grafo[de].adjacentes[k].ramo->tap_pri[2]*grafo[de].adjacentes[k].ramo->tap_sec[2];
                    }
                    else{
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
    for(i=0; i<nVariaveis; i++)
    {
        k = (int)regua[i];
        fase = (int) cabs((regua[i] - k)*10);
        if (regua[i] > 0){ //magnitude de tensão
            x[i] = cabs(grafo[k].V[fase]); //mantém o ângulo anterior e altera a magnitude            
        }
        else{ //ângulo de tensão
            fase = fase;
            k = -k;
            x[i] = carg(grafo[k].V[fase]);
        }
    }
}

//Função inicialização do vetor x a partir de arquivo de entrada
void incializa_vetor_x_leitura(GRAFO *grafo, long int numeroBarras, ALIMENTADOR *alimentadores, long int numeroAlimentadores, double *x, double *regua, long int nVariaveis){
    int i, k, fase;
    BOOL visitado[numeroBarras];
    
    
    //Inicialização com base nas tensões lidas no arquivo Vinicial.csv
    for(i=0; i<numeroBarras; i++){ 
        visitado[i] = false;
    }
    for(i=0; i<numeroAlimentadores; i++)
    {
        //Tensão Inicial da subestação
        FILABARRAS *barraAtual = &alimentadores[i].rnp[0];
        while(barraAtual != NULL)
        {
            int de = barraAtual->idNo;
            if (visitado[de] == false){
                grafo[de].V[0] = grafo[de].barra->Vinicial[0];
                grafo[de].V[1] = grafo[de].barra->Vinicial[1];
                grafo[de].V[2] = grafo[de].barra->Vinicial[2];

                visitado[de] = true;
                barraAtual = barraAtual->prox;
            }
        }        
    }
    
    //Montagem do vetor x e da régua
    for(i=0; i<nVariaveis; i++)
    {
        k = (int)regua[i];
        fase = (int) cabs((regua[i] - k)*10);
        if (regua[i] > 0){ //magnitude de tensão
            x[i] = cabs(grafo[k].V[fase]); //mantém o ângulo anterior e altera a magnitude            
        }
        else{ //ângulo de tensão
            fase = fase;
            k = -k;
            x[i] = carg(grafo[k].V[fase]);
        }
    }
    
}

//Função atualiza as grandezas elétricas da rede (fluxos, injeções, correntes)
void atualiza_Rede(GRAFO *grafo, long int numeroBarras){
    int i,j,k,idMed, de, para,ramo,fase;
    __complex__ double *Saux, *Iaux;
    BOOL visitado[numeroBarras];
    
    Saux = (__complex__ double*)malloc(3 * sizeof(__complex__ double));
    Saux[0] = 0;
    Saux[1] = 0;
    Saux[2] = 0;
    Iaux = (__complex__ double*)malloc(3 * sizeof(__complex__ double));
    Iaux[0] = 0;
    Iaux[1] = 0;
    Iaux[2] = 0;
    
   
    //Percorre o grafo atualizando o cálculo de h(x))
    for(i=0;i<numeroBarras;i++){        
        Sk(grafo, i, Saux);
        grafo[i].S[0] = Saux[0];
        grafo[i].S[1] = Saux[1];
        grafo[i].S[2] = Saux[2];

        grafo[i].Cur[0] = conj(Saux[0]/grafo[i].V[0]);
        grafo[i].Cur[1] = conj(Saux[1]/grafo[i].V[1]);
        grafo[i].Cur[2] = conj(Saux[2]/grafo[i].V[2]);
        
        //Percorre os ramos adjacentes
        for(k=0;k<grafo[i].numeroAdjacentes;k++){
            if (i == grafo[i].adjacentes[k].ramo->k){
                Skm(&grafo[i], &grafo[grafo[i].adjacentes[k].idNo], grafo[i].adjacentes[k].ramo, Saux);
                Ikm(&grafo[i], &grafo[grafo[i].adjacentes[k].idNo], grafo[i].adjacentes[k].ramo, Iaux);
            }
            else{
                Smk(&grafo[grafo[i].adjacentes[k].idNo],&grafo[i], grafo[i].adjacentes[k].ramo, Saux);
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
    free(Saux);free(Iaux);
}

//Função atualiza modelo de medição conforme as grandezas elétricas calculadas
void atualiza_Modelo(GRAFO *grafo, long int numeroBarras, long int nmed, DMED *medidas){
    int i,j,k,idMed, de, para,ramo,fase;
    __complex__ double *Saux, *Iaux;
    BOOL visitado[numeroBarras];
    
    for(idMed=0;idMed<nmed;idMed++){
        if (medidas[idMed].PARA == -1){ //Medidor instalado em uma barra
            i = medidas[idMed].k;
            fase = medidas[idMed].fases - 1;//Revisar para trifásico genérico - atual somente A ou B ou C - medida no Delta por exemplo

            switch (medidas[idMed].tipo){
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
        else{ //Medidor instalado em um ramo
            i = medidas[idMed].k;
            fase = medidas[idMed].fases - 1;//Revisar para trifásico genérico - atual somente A ou B ou C - medida no Delta por exemplo
            for(k=0;k<grafo[i].numeroAdjacentes;k++){
                if (grafo[i].adjacentes[k].idNo == medidas[idMed].m){
                    switch (medidas[idMed].tipo){
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
void atualiza_h(GRAFO *grafo, long int numeroBarras, long int nmed, DMED *medidas){
    int i,j,k,idMed, de, para,ramo,fase;
    __complex__ double *Saux, *Iaux;
    BOOL visitado[numeroBarras];
    
    Saux = (__complex__ double*)malloc(3 * sizeof(__complex__ double));
    Saux[0] = 0;
    Saux[1] = 0;
    Saux[2] = 0;
    Iaux = (__complex__ double*)malloc(3 * sizeof(__complex__ double));
    Iaux[0] = 0;
    Iaux[1] = 0;
    Iaux[2] = 0;
    
   
    //Percorre o grafo atualizando o cálculo de h(x))
    for(i=0;i<numeroBarras;i++){
        //Percorre os medidores da barra
        if(grafo[i].nmed>0){
            Sk(grafo, i, Saux);
            grafo[i].S[0] = Saux[0];
            grafo[i].S[1] = Saux[1];
            grafo[i].S[2] = Saux[2];
            
            grafo[i].Cur[0] = conj(Saux[0]/grafo[i].V[0]);
            grafo[i].Cur[1] = conj(Saux[1]/grafo[i].V[1]);
            grafo[i].Cur[2] = conj(Saux[2]/grafo[i].V[2]);
        }
        
        for(j=0;j<grafo[i].nmed;j++){
            idMed = grafo[i].medidores[j]->id;
            fase = grafo[i].medidores[j]->fases - 1; //Revisar para trifásico genérico - atual somente A ou B ou C
            
            switch (medidas[idMed].tipo){
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
        for(k=0;k<grafo[i].numeroAdjacentes;k++){
            if(grafo[i].adjacentes[k].nmed > 0){
                if (i == grafo[i].adjacentes[k].ramo->k){
                    Skm(&grafo[i], &grafo[grafo[i].adjacentes[k].idNo], grafo[i].adjacentes[k].ramo, Saux);
                    Ikm(&grafo[i], &grafo[grafo[i].adjacentes[k].idNo], grafo[i].adjacentes[k].ramo, Iaux);
                }
                else{
                    Smk(&grafo[grafo[i].adjacentes[k].idNo],&grafo[i], grafo[i].adjacentes[k].ramo, Saux);
                    Imk(&grafo[grafo[i].adjacentes[k].idNo], &grafo[i], grafo[i].adjacentes[k].ramo, Iaux);
                }
                grafo[i].adjacentes[k].S[0] = Saux[0];
                grafo[i].adjacentes[k].S[1] = Saux[1];
                grafo[i].adjacentes[k].S[2] = Saux[2];

                grafo[i].adjacentes[k].Cur[0] = Iaux[0];
                grafo[i].adjacentes[k].Cur[1] = Iaux[1];
                grafo[i].adjacentes[k].Cur[2] = Iaux[2];
            }
                        
            for(j=0;j<grafo[i].adjacentes[k].nmed;j++){
                idMed = grafo[i].adjacentes[k].medidores[j]->id;
                fase = grafo[i].adjacentes[k].medidores[j]->fases - 1; //Revisar para trifásico genérico - atual somente A ou B ou C

                switch (medidas[idMed].tipo){
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
    free(Saux);free(Iaux);
}

//Função atualiza matriz H
void atualiza_H(GRAFO *grafo, long int numeroBarras, DRAM *ramos, DMED *medidas, long int numeroMedidas){
    int i,j,k,idMed, de, para,adj,ramo,opt,fase;
    __complex__ double *dSaux, *dIaux;
    int it = 0;
    
    dSaux = (__complex__ double*)malloc(3 * sizeof(__complex__ double));
    dSaux[0] = 0;
    dSaux[1] = 0;
    dSaux[2] = 0;
    dIaux = (__complex__ double*)malloc(3 * sizeof(__complex__ double));
    dIaux[0] = 0;
    dIaux[1] = 0;
    dIaux[2] = 0;
    
    //Atualiza a matriz H de acordo com a régua
    for(i=0;i<numeroMedidas;i++){
        switch (medidas[i].tipo){
            case 0: //0: Fluxo de Potência Ativa - kW
                for(j=0;j<medidas[i].nvar;j++){
                    k = (int) medidas[i].reguaH[j];
                    fase = (int) cabs((medidas[i].reguaH[j] - k)*10);
                    
                    de = medidas[i].k;
                    para = medidas[i].m;
                    ramo = medidas[i].ramo;
                    
                    if (ramos[ramo].k == de){
                        if(medidas[i].reguaH[j]>=0){
                            if (de == k) opt = 0;
                            else opt = 2;
                        }
                        else{
                            if (de == -k) opt = 1;
                            else opt = 3;
                        }
                        dSkm(&grafo[de], &grafo[para], &ramos[ramo], dSaux, opt, fase);
                        medidas[i].H[j] = __real__ dSaux[medidas[i].fases - 1];
                    }
                    else{
                        if(medidas[i].reguaH[j]>=0){
                            if (de == k) opt = 2;
                            else opt = 0;
                        }
                        else{
                            if (de == -k) opt = 3;
                            else opt = 1;
                        }
                        dSmk(&grafo[para], &grafo[de], &ramos[ramo], dSaux, opt, fase);
                        medidas[i].H[j] = __real__ dSaux[medidas[i].fases - 1];
                    }
                }
                break;
            case 1: //1: Fluxo de Potência Reativa - kVAr
                for(j=0;j<medidas[i].nvar;j++){
                    k = (int) medidas[i].reguaH[j];
                    fase = (int) cabs((medidas[i].reguaH[j] - k)*10);
                    
                    de = medidas[i].k;
                    para = medidas[i].m;
                    ramo = medidas[i].ramo;
                    
                    if (ramos[ramo].k == de){
                        if(medidas[i].reguaH[j]>=0){
                            if (de == k) opt = 0;
                            else opt = 2;
                        }
                        else{
                            if (de == -k) opt = 1;
                            else opt = 3;
                        }
                        dSkm(&grafo[de], &grafo[para], &ramos[ramo], dSaux, opt, fase);
                        medidas[i].H[j] = __imag__ dSaux[medidas[i].fases - 1];
                    }
                    else{
                        if(medidas[i].reguaH[j]>=0){
                            if (de == k) opt = 2;
                            else opt = 0;
                        }
                        else{
                            if (de == -k) opt = 3;
                            else opt = 1;
                        }
                        dSmk(&grafo[para], &grafo[de], &ramos[ramo], dSaux, opt, fase);
                        medidas[i].H[j] = __imag__ dSaux[medidas[i].fases - 1];
                    }
                }
                break;    
            case 2: //2: Injeção de Potência Ativa - kW
                for(j=0;j<medidas[i].nvar;j++){
                    k = (int) medidas[i].reguaH[j];
                    fase = (int) cabs((medidas[i].reguaH[j] - k)*10);
                    
                    de = medidas[i].k;
                    
                    if(medidas[i].reguaH[j]>=0){
                        opt = 0;
                    }
                    else{
                        opt = 1;
                    }
                    dSk(grafo, de, dSaux, opt, cabs(k), fase);
                   
                    medidas[i].H[j] = __real__ dSaux[medidas[i].fases - 1];
                }
                break;     
            case 3: //3: Injeção de Potência Reativa - kVAr
                for(j=0;j<medidas[i].nvar;j++){
                    k = (int) medidas[i].reguaH[j];
                    fase = (int) cabs((medidas[i].reguaH[j] - k)*10);
                    
                    de = medidas[i].k;
                    
                    if(medidas[i].reguaH[j]>=0){
                        opt = 0;
                    }
                    else{
                        opt = 1;
                    }
                    dSk(grafo, de, dSaux, opt, cabs(k), fase);
                   
                    medidas[i].H[j] = __imag__ dSaux[medidas[i].fases - 1];
                }break;
            case 4: //4: Magnitude de tensão
                for(j=0;j<medidas[i].nvar;j++){
                    k = (int) medidas[i].reguaH[j];
                    fase = (int) cabs((medidas[i].reguaH[j] - k)*10);
                    
                    if ((medidas[i].fases -1 == fase)&&(medidas[i].reguaH[j]>0))
                        medidas[i].H[j] = 1;
                }
                break;
            case 5: //5: Ângulo de tensão
                for(j=0;j<medidas[i].nvar;j++){
                    k = (int) medidas[i].reguaH[j];
                    fase = (int) cabs((medidas[i].reguaH[j] - k)*10);
                    
                    if ((medidas[i].fases -1 == fase)&&(medidas[i].reguaH[j]<0))
                        medidas[i].H[j] = 1;
                }
                break;
            case 7: //7: Magnitude de Corrente
                for(j=0;j<medidas[i].nvar;j++){
                    k = (int) medidas[i].reguaH[j];
                    fase = (int) cabs((medidas[i].reguaH[j] - k)*10);
                    
                    de = medidas[i].k;
                    para = medidas[i].m;
                    ramo = medidas[i].ramo;
                    
                    if (ramos[ramo].k == de){
                        if(medidas[i].reguaH[j]>=0){
                            if (de == k) opt = 0;
                            else opt = 2;
                        }
                        else{
                            if (de == -k) opt = 1;
                            else opt = 3;
                        }
                        if (it == 0) dSkm(&grafo[de], &grafo[para], &ramos[ramo], dSaux, opt, fase);
                        else dIkm(&grafo[de], &grafo[para], &ramos[ramo], dSaux, opt, fase);
                        medidas[i].H[j] = dSaux[medidas[i].fases - 1];
                    }
                    else{
                        if(medidas[i].reguaH[j]>=0){
                            if (de == k) opt = 2;
                            else opt = 0;
                        }
                        else{
                            if (de == -k) opt = 3;
                            else opt = 1;
                        }
                        dSmk(&grafo[para], &grafo[de], &ramos[ramo], dSaux, opt, fase);
                        medidas[i].H[j] = __real__ dSaux[medidas[i].fases - 1];
                    }
                }
                break;    
                
        }        
    }
    free(dSaux);free(dIaux);
}

//Função atualiza matriz H
void atualiza_H_ret(GRAFO *grafo, long int numeroBarras, DRAM *ramos, DMED *medidas, long int numeroMedidas){
    int i,j,k,idMed, de, para,adj,ramo,opt,fase;
    __complex__ double *dSaux, *dIaux;
    BOOL visitado[numeroBarras];
    
    dSaux = (__complex__ double*)malloc(3 * sizeof(__complex__ double));
    dSaux[0] = 0;
    dSaux[1] = 0;
    dSaux[2] = 0;
    dIaux = (__complex__ double*)malloc(3 * sizeof(__complex__ double));
    dIaux[0] = 0;
    dIaux[1] = 0;
    dIaux[2] = 0;
    
    //Atualiza a matriz H de acordo com a régua
    for(i=0;i<numeroMedidas;i++){
        switch (medidas[i].tipo){
            case 0: //0: Fluxo de Potência Ativa - kW
                for(j=0;j<medidas[i].nvar;j++){
                    k = (int) medidas[i].reguaH[j];
                    fase = (int) cabs((medidas[i].reguaH[j] - k)*10);
                    
                    de = medidas[i].k;
                    para = medidas[i].m;
                    ramo = medidas[i].ramo;
                    
                    if (ramos[ramo].k == de){
                        if(medidas[i].reguaH[j]>=0){
                            if (de == k) opt = 0;
                            else opt = 2;
                        }
                        else{
                            if (de == -k) opt = 1;
                            else opt = 3;
                        }
                        dSkm_ret(&grafo[de], &grafo[para], &ramos[ramo], dSaux, opt, fase);
                        medidas[i].H[j] = __real__ dSaux[medidas[i].fases - 1];
                    }
                    else{
                        if(medidas[i].reguaH[j]>=0){
                            if (de == k) opt = 2;
                            else opt = 0;
                        }
                        else{
                            if (de == -k) opt = 3;
                            else opt = 1;
                        }
                        dSmk_ret(&grafo[para], &grafo[de], &ramos[ramo], dSaux, opt, fase);
                        medidas[i].H[j] = __real__ dSaux[medidas[i].fases - 1];
                    }
                }
                break;
            case 1: //1: Fluxo de Potência Reativa - kVAr
                for(j=0;j<medidas[i].nvar;j++){
                    k = (int) medidas[i].reguaH[j];
                    fase = (int) cabs((medidas[i].reguaH[j] - k)*10);
                    
                    de = medidas[i].k;
                    para = medidas[i].m;
                    ramo = medidas[i].ramo;
                    
                    if (ramos[ramo].k == de){
                        if(medidas[i].reguaH[j]>=0){
                            if (de == k) opt = 0;
                            else opt = 2;
                        }
                        else{
                            if (de == -k) opt = 1;
                            else opt = 3;
                        }
                        dSkm_ret(&grafo[de], &grafo[para], &ramos[ramo], dSaux, opt, fase);
                        medidas[i].H[j] = __imag__ dSaux[medidas[i].fases - 1];
                    }
                    else{
                        if(medidas[i].reguaH[j]>=0){
                            if (de == k) opt = 2;
                            else opt = 0;
                        }
                        else{
                            if (de == -k) opt = 3;
                            else opt = 1;
                        }
                        dSmk_ret(&grafo[para], &grafo[de], &ramos[ramo], dSaux, opt, fase);
                        medidas[i].H[j] = __imag__ dSaux[medidas[i].fases - 1];
                    }
                }
                break;    
            case 2: //2: Injeção de Potência Ativa - kW
                for(j=0;j<medidas[i].nvar;j++){
                    k = (int) medidas[i].reguaH[j];
                    fase = (int) cabs((medidas[i].reguaH[j] - k)*10);
                    
                    de = medidas[i].k;
                    
                    if(medidas[i].reguaH[j]>=0){
                        opt = 0;
                    }
                    else{
                        opt = 1;
                    }
                    dSk_ret(grafo, de, dSaux, opt, cabs(k), fase);
                   
                    medidas[i].H[j] = __real__ dSaux[medidas[i].fases - 1];
                }
                break;     
            case 3: //3: Injeção de Potência Reativa - kVAr
                for(j=0;j<medidas[i].nvar;j++){
                    k = (int) medidas[i].reguaH[j];
                    fase = (int) cabs((medidas[i].reguaH[j] - k)*10);
                    
                    de = medidas[i].k;
                    
                    if(medidas[i].reguaH[j]>=0){
                        opt = 0;
                    }
                    else{
                        opt = 1;
                    }
                    dSk_ret(grafo, de, dSaux, opt, cabs(k), fase);
                   
                    medidas[i].H[j] = __imag__ dSaux[medidas[i].fases - 1];
                }break;
            case 4: //4: Magnitude de tensão 
                for(j=0;j<medidas[i].nvar;j++){
                    k = (int) medidas[i].reguaH[j];
                    fase = (int) cabs((medidas[i].reguaH[j] - k)*10);
                    
                    if(medidas[i].reguaH[j]>=0){
                        if (medidas[i].fases -1 == fase)
                            medidas[i].H[j] = __real__ grafo[k].V[fase]/cabs(grafo[k].V[fase]);
                    }
                    else{
                        k = -k;
                        if (medidas[i].fases -1 == fase)
                            medidas[i].H[j] = __imag__ grafo[k].V[fase]/cabs(grafo[k].V[fase]);
                    }
                }
                break;
            case 5: //5: Magnitude de tensão 
                for(j=0;j<medidas[i].nvar;j++){
                    k = (int) medidas[i].reguaH[j];
                    fase = (int) cabs((medidas[i].reguaH[j] - k)*10);
                    
                    if(medidas[i].reguaH[j]>=0){
                        if (medidas[i].fases -1 == fase)
                            medidas[i].H[j] = -1 * __imag__ grafo[k].V[fase] / pow(__real__ grafo[k].V[fase],2) * 1 / (1 + pow(__imag__ grafo[k].V[fase]/__real__ grafo[k].V[fase],2));
                    }
                    else{
                        k = -k;
                        if (medidas[i].fases -1 == fase)
                            medidas[i].H[j] = 1 / __real__ grafo[k].V[fase] * 1 / (1 + pow(__imag__ grafo[k].V[fase]/__real__ grafo[k].V[fase],2));
                    }
                }
                break;
        }        
    }
    free(dSaux);free(dIaux);
}


//Função que atualiza as variaveis de estado em função do vetor x
void atualiza_estado(GRAFO *grafo, double *x, double *regua, long int nVariaveis){
    int i,k,fase;
    
    //nVariaveis é o tamanho do vetor x
    //regua salva a respectiva barra e fase de cada elemento do vetor x (+ é magnitude de tensão e - é ângulo)
    // Ex: "+1.2" = tensão barra 1 fase 2 
    // Ex: "-2.0" = ângulo barra 2 fase 0 
    
    for(i=0; i<nVariaveis; i++)
    {
        k = (int)regua[i];
        fase = (int) cabs((regua[i] - k)*10);
        if (regua[i] > 0){ //magnitude de tensão
            grafo[k].V[fase] = x[i]*grafo[k].V[fase]/cabs(grafo[k].V[fase]); //mantém o ângulo anterior e altera a magnitude            
        }
        else{ //ângulo de tensão
            fase = fase;
            k = -k;
            grafo[k].V[fase] = cabs(grafo[k].V[fase])*(cos(x[i]) + I*sin(x[i])); //mantém a magnitude e altera o ângulo
        }     
    }
}

//Tratamento da referência (por alimentador)
void tratamento_referencia(long int *ref_1, long int *ref_2, ALIMENTADOR *alimentador, double *regua, long int nVariaveis){ 
    int i, j, k, n_ref = 3;
    
    //Inserir aqui fatoração da matriz H
    //atualizaH
    //fatoração e obtenção da Hdelta
    //número de referências
    
    //Inicialização com base nas tensões lidas no arquivo Vinicial.csv
    for(j=0; j<nVariaveis; j++){
        if (regua[j] < 0){
            k = (int)regua[j];
            k = -k;
            if (k == alimentador->noRaiz){
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
void estimadorWLS(GRAFO *grafo, long int numeroBarras, DMED *medidas, long int **numeroMedidas, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos,double Sbase){
    long int nmed,nvar;
    int i,j,k, r;
    double *z = NULL,**h = NULL,***H = NULL,**W = NULL, *x = NULL, *regua = NULL, aux = 0;
    
    printf("Estimador de Estado WLS Trifásico...\n");
    //--------------------------------------------------------------------------
    //Alocação de memória das variáveis do estimador de estado
    nmed = 0;
    for (i = 0; i < 9; i++){ 
        for (j = 0; j < 8; j++){
            nmed = nmed + numeroMedidas[i][j];
        }
    }
    nvar = 0;
    //printf("numero barras: %d\n", numeroBarras);
    for (i = 0; i < numeroBarras; i++){
        switch (grafo[i].fases){
            case 1:
                nvar +=2;
                break;
            case 2:
                nvar +=2;
                break;
            case 3:
                nvar +=2;
                break;
            case 4:
                nvar +=4;
                break;    
            case 5:
                nvar +=4;
                break;    
            case 6:
                nvar +=4;
                break;    
            case 7:
                nvar +=6;
                break;    
        }
    }    
    //printf("nmed: %d\n", nmed);
    //printf("nvar: %d\n", nvar);
    if ((z = (double *)malloc( (nmed) * sizeof(double)))==NULL){
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor z!!!!");
        exit(1); 
    }
    if ((h = malloc( (nmed) * sizeof(double*)))==NULL){
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor h!!!!");
        exit(1); 
    }
    if ((x = (double *)malloc( (nvar) * sizeof(double)))==NULL){
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor x!!!!");
        exit(1); 
    }
    if ((regua = (double *)malloc( (nvar) * sizeof(double)))==NULL){
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor regua!!!!");
        exit(1); 
    }
    
    H = (double***)malloc(nmed * sizeof(double**)); 
    for (i = 0; i < nmed; i++){ 
         H[i] = (double**) malloc(nvar * sizeof(double*));
         for (j = 0; j < nvar; j++){
              H[i][j] = &aux;
         }
    }
    //--------------------------------------------------------------------------
    // Direcionamento dos ponteiros que compõem a régua, vetor h e matriz H
    //monta a regua (ordenaçao do vetor x) - V e teta conforme o grafo
    j=0;
    for(i=0;i<numeroBarras;i++){
        aux = (double) grafo[i].idNo;
        aux += 0.01;
        switch (grafo[i].fases){
            case 1:
                regua[j] = aux;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 2:
                regua[j] = aux + 0.1;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 3:
                regua[j] = aux + 0.2;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 4:
                regua[j] = aux;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.1;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 5:
                regua[j] = aux;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.2;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 6:
                regua[j] = aux+0.1;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.2;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 7:
                regua[j] = aux;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.1;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.2;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;    
        }
    }
    aux = 0;
    //printf("nmed: %d\n", nmed);
    //printf("nvar: %d\n", nvar);
    //Tratamento da referência
    long int ref_1, ref_2;
    tratamento_referencia(&ref_1, &ref_2, &alimentadores[0], regua, nvar);
    
    tira_refs_regua(nvar, ref_1, ref_2, regua); 
    nvar = nvar - (ref_2 - ref_1 +1);  
    //printf("tira refs\n");
    //vetor h aponta para a estrutura de dados das medidas
    for(i=0;i<nmed;i++){
        h[i] = &medidas[i].h;
    }
    //Matriz H aponta para a estrutura de dados das medidas
    for(i=0;i<nmed;i++){
        for(j=0;j<medidas[i].nvar;j++){
            for(r = 0;r<nvar;r++){
                if (cabs(medidas[i].reguaH[j]-regua[r]) < EPS){
                    H[i][r] = &medidas[i].H[j];
                    break;
                }
            }
        }
    }  
    //--------------------------------------------------------------------------
    //Estimação de Estado    
    monta_z(z,nmed,medidas);
    //monta_W(NULL,nmed,medidas);
    monta_W_cte(W,nmed,medidas);
    //monta_W_Ident(NULL,nmed,medidas);
    
    incializa_vetor_x(grafo, numeroBarras, alimentadores, numeroAlimentadores,x,regua,nvar);
//    incializa_vetor_x_leitura(grafo, numeroBarras, alimentadores, numeroAlimentadores,x,regua,nvar); //Função que le de arquivo externo a condição inicial
    
    double tol = 0.000001;
    clock_t tIni = clock();
    int conv = otimiza_Gauss_NewtonQR(z, h, H, grafo, numeroBarras, ramos, medidas, nvar, nmed, regua, x, tol,ref_1,ref_2);
//    int conv = otimiza_Gauss_Newton_sparsePCGLS(z, h, H, W, grafo, numeroBarras, ramos, medidas, nvar - 3, nmed, regua, x, tol, ref_1, ref_2);
    
    clock_t t1 = clock();
    double tempoWLS = (double)(t1-tIni)/CLOCKS_PER_SEC;
    printf("\nEstimação WLS: %lf",tempoWLS);
    
    exportaCasoReferencia(grafo, numeroBarras, Sbase);
    exportaEstado(grafo,regua,nvar);
    imprimeEstado(grafo,numeroBarras);
    
//    atualiza_H_ret(grafo, numeroBarras, ramos, medidas, nmed);
//    exportaPrioriSCADA(grafo,regua, H, W, nvar, nmed, 0);
//    exportaPrioriQR(grafo,regua,medidas, nvar, nmed, 0);
    
    FILE *residuo;
    residuo = fopen("residuo.txt","wt");
    for(i=0;i<nmed;i++){
        fprintf(residuo,"%.10lf\n",z[i] - medidas[i].h);
    }
    fclose(residuo);
    
    free(z);free(h);free(x);free(regua);
    for (i=0;i<nmed;i++ )free(H[i]);
    free(H);
}


//------------------------------------------------------------------------------
//
// ESTIMADOR WLS CONVENCIONAL VIA WLS com restrições de Igualdade para virtuais e formulação esparsa do resíduo
//
//------------------------------------------------------------------------------
void estimadorWLShachtel(GRAFO *grafo, long int numeroBarras, DMED *medidas, long int **numeroMedidas,DMED *virtuais, long int **numeroVirtuais, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos){
    long int nmed,nvar,nvir;
    int i,j,k, r;
    double *z = NULL,**h = NULL,***H = NULL,**W = NULL, *x = NULL, *regua = NULL, aux = 0;
    double **c = NULL, ***C = NULL;
    
    printf("Estimador de Estado WLS Trifásico Hachtel...\n");
    //--------------------------------------------------------------------------
    //Alocação de memória das variáveis do estimador de estado
    nmed = 0;
    for (i = 0; i < 9; i++){ 
        for (j = 0; j < 8; j++){
            nmed = nmed + numeroMedidas[i][j];
        }
    }
    nvir = 0;
    for (i = 0; i < 9; i++){ 
        for (j = 0; j < 8; j++){
            nvir = nvir + numeroVirtuais[i][j];
        }
    }
    nvar = 0;
    for (i = 0; i < numeroBarras; i++){
        switch (grafo[i].fases){
            case 1:
                nvar +=2;
                break;
            case 2:
                nvar +=2;
                break;
            case 3:
                nvar +=2;
                break;
            case 4:
                nvar +=4;
                break;    
            case 5:
                nvar +=4;
                break;    
            case 6:
                nvar +=4;
                break;    
            case 7:
                nvar +=6;
                break;    
        }
    }
    
    if ((z = (double *)malloc( (nmed) * sizeof(double)))==NULL){
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor z!!!!");
        exit(1); 
    }
    if ((h = malloc( (nmed) * sizeof(double*)))==NULL){
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor h!!!!");
        exit(1); 
    }
    if ((x = (double *)malloc( (nvar) * sizeof(double)))==NULL){
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor x!!!!");
        exit(1); 
    }
    if ((regua = (double *)malloc( (nvar) * sizeof(double)))==NULL){
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor regua!!!!");
        exit(1); 
    }
    if ((c = malloc( (nvir) * sizeof(double*)))==NULL){
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor c!!!!");
        exit(1); 
    }
    
    W = (double**)malloc(nmed * sizeof(double*)); 
    for (i = 0; i < nmed; i++){ 
         W[i] = (double*) malloc(nmed * sizeof(double));
         for (j = 0; j < nmed; j++){
              W[i][j] = 0;
         }
    }
    H = (double***)malloc(nmed * sizeof(double**)); 
    for (i = 0; i < nmed; i++){ 
         H[i] = (double**) malloc(nvar * sizeof(double*));
         for (j = 0; j < nvar; j++){
              H[i][j] = &aux;
         }
    }
    C = (double***)malloc(nvir * sizeof(double**)); 
    for (i = 0; i < nvir; i++){ 
         C[i] = (double**) malloc(nvar * sizeof(double*));
         for (j = 0; j < nvar; j++){
              C[i][j] = &aux;
         }
    }
    
    //--------------------------------------------------------------------------
    // Direcionamento dos ponteiros que compõem a régua, vetor h e matriz H
    //monta a regua (ordenaçao do vetor x) - V e teta conforme o grafo
    j=0;
    for(i=0;i<numeroBarras;i++){
        aux = (double) grafo[i].idNo;
        aux += 0.01;
        switch (grafo[i].fases){
            case 1:
                regua[j] = aux;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 2:
                regua[j] = aux + 0.1;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 3:
                regua[j] = aux + 0.2;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 4:
                regua[j] = aux;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.1;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 5:
                regua[j] = aux;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.2;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 6:
                regua[j] = aux+0.1;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.2;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 7:
                regua[j] = aux;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.1;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.2;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;    
        }
    }
    aux = 0;
    
    
    //vetor h aponta para a estrutura de dados das medidas
    for(i=0;i<nmed;i++){
        h[i] = &medidas[i].h;
    }
    //Matriz H aponta para a estrutura de dados das medidas
    for(i=0;i<nmed;i++){
        for(j=0;j<medidas[i].nvar;j++){
            for(r = 0;r<nvar;r++){
                if (cabs(medidas[i].reguaH[j]-regua[r]) < EPS){
                    //Atualiza a Matriz H
                    H[i][r] = &medidas[i].H[j];
                    break;
                }
            }
        }
    }
    //vetor c aponta para a estrutura de dados das medidas virtuais
    for(i=0;i<nvir;i++){
        c[i] = &virtuais[i].h;
    }
    //Matriz C aponta para a estrutura de dados das medidas virtuais
    for(i=0;i<nvir;i++){
        for(j=0;j<virtuais[i].nvar;j++){
            for(r = 0;r<nvar;r++){
                if (cabs(virtuais[i].reguaH[j]-regua[r]) < EPS){
                    //Atualiza a Matriz H
                    C[i][r] = &virtuais[i].H[j];
                    break;
                }
            }
        }
    }
    
    //--------------------------------------------------------------------------
    //Estimação de Estado    
    monta_z(z,nmed,medidas);
    //monta_W(W,nmed,medidas);
    monta_W_Ident(W,nmed,medidas);
    
    //incializa_vetor_x(grafo, numeroBarras, alimentadores, numeroAlimentadores,x,regua,nvar);
    incializa_vetor_x_leitura(grafo, numeroBarras, alimentadores, numeroAlimentadores,x,regua,nvar); //Função que le de arquivo externo a condição inicial
       
    //Tratamento da referência
    long int ref_1, ref_2;
    tratamento_referencia(&ref_1, &ref_2, &alimentadores[0], regua, nvar);    
    
    
    double tol = 0.00001;
    clock_t tIni = clock();
    int conv = otimiza_Gauss_Newton_sparseHachtel_Virtuais(z, h, c, W, grafo, numeroBarras, ramos, medidas, virtuais, nvar-3, nmed, nvir, regua, x, tol,ref_1,ref_2);
    clock_t t1 = clock();
    double tempoWLS = (double)(t1-tIni)/CLOCKS_PER_SEC;
    printf("\nEstimação WLS: %lf",tempoWLS);
    printf("\nTensoes Nodais (p.u.):\n");
    for(i=0; i<numeroBarras; i++){ 
        //Retangulares
        //printf("%d\tVa: %.5lf + j%.5lf\tVb: %.5lf + j%.5lf\tVc: %.5lf + j%.5lf\n",grafo[i].barra->ID,__real__ grafo[i].V[0],__imag__ grafo[i].V[0],__real__ grafo[i].V[1],__imag__ grafo[i].V[1],__real__ grafo[i].V[2],__imag__ grafo[i].V[2]);
        //Polares
        switch (grafo[i].fases){
            case 1:
                printf("%d\tVa: %.5lf | %.3lf \tVb:    -    |    -   \tVc:    -    |    -   \n",grafo[i].barra->ID,cabs(grafo[i].V[0]),carg(grafo[i].V[0])*180/PI);
                break;
            case 2:
                printf("%d\tVa:    -    |    -    \tVb: %.5lf | %.3lf\tVc:    -    |    -   \n",grafo[i].barra->ID,cabs(grafo[i].V[1]),carg(grafo[i].V[1])*180/PI);
                break;
            case 3:
                printf("%d\tVa:    -    |    -    \tVb:    -    |    -   \tVc: %.5lf | %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[2]),carg(grafo[i].V[2])*180/PI);
                break;
            case 4:
                printf("%d\tVa: %.5lf | %.3lf \tVb: %.5lf | %.3lf\tVc:    -    |    -   \n",grafo[i].barra->ID,cabs(grafo[i].V[0]),carg(grafo[i].V[0])*180/PI,cabs(grafo[i].V[1]),carg(grafo[i].V[1])*180/PI);
                break;
            case 5:
                printf("%d\tVa: %.5lf | %.3lf \tVb:    -    |    -   \tVc: %.5lf | %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[0]),carg(grafo[i].V[0])*180/PI,cabs(grafo[i].V[2]),carg(grafo[i].V[2])*180/PI);
                break;
            case 6:
                printf("%d\tVa:    -    |    -    \tVb: %.5lf | %.3lf\tVc: %.5lf | %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[1]),carg(grafo[i].V[1])*180/PI,cabs(grafo[i].V[2]),carg(grafo[i].V[2])*180/PI);
                break;
            case 7:
                printf("%d\tVa: %.5lf | %.3lf \tVb: %.5lf | %.3lf\tVc: %.5lf | %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[0]),carg(grafo[i].V[0])*180/PI,cabs(grafo[i].V[1]),carg(grafo[i].V[1])*180/PI,cabs(grafo[i].V[2]),carg(grafo[i].V[2])*180/PI);
                break;    
        }
    }
    atualiza_H_ret(grafo, numeroBarras, ramos, medidas, nmed);
    
    free(z);free(h);free(H);free(W);free(x);
}


//------------------------------------------------------------------------------
//
// FLUXO DE POTÊNCIA VIA MINIMIZAÇÃO QR
//
//------------------------------------------------------------------------------
void fluxoPotencia_NRQR(GRAFO *grafo, long int numeroBarras, DMED *medidas, long int **numeroMedidas, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos,double Sbase){
    long int nmed,nvar;
    int i,j,k, r;
    double *z = NULL,**h = NULL,***H = NULL, *x = NULL, *regua = NULL, aux = 0;
    
    printf("Calculo de Fluxo de Potencia via Minimizacao de Residuo via NR-QR...\n");
    //--------------------------------------------------------------------------
    //Alocação de memória das variáveis do estimador de estado
    nmed = 0;
    for (i = 0; i < 9; i++){ 
        for (j = 0; j < 8; j++){
            nmed = nmed + numeroMedidas[i][j];
        }
    }
    nvar = 0;
    for (i = 0; i < numeroBarras; i++){
        switch (grafo[i].fases){
            case 1:
                nvar +=2;
                break;
            case 2:
                nvar +=2;
                break;
            case 3:
                nvar +=2;
                break;
            case 4:
                nvar +=4;
                break;    
            case 5:
                nvar +=4;
                break;    
            case 6:
                nvar +=4;
                break;    
            case 7:
                nvar +=6;
                break;    
        }
    }    
    
    if ((z = (double *)malloc( (nmed) * sizeof(double)))==NULL){
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor z!!!!");
        exit(1); 
    }
    if ((h = malloc( (nmed) * sizeof(double*)))==NULL){
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor h!!!!");
        exit(1); 
    }
    if ((x = (double *)malloc( (nvar) * sizeof(double)))==NULL){
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor x!!!!");
        exit(1); 
    }
    if ((regua = (double *)malloc( (nvar) * sizeof(double)))==NULL){
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor regua!!!!");
        exit(1); 
    }
    
    H = (double***)malloc(nmed * sizeof(double**)); 
    for (i = 0; i < nmed; i++){ 
         H[i] = (double**) malloc(nvar * sizeof(double*));
         for (j = 0; j < nvar; j++){
              H[i][j] = &aux;
         }
    }
    //--------------------------------------------------------------------------
    // Direcionamento dos ponteiros que compõem a régua, vetor h e matriz H
    //monta a regua (ordenaçao do vetor x) - V e teta conforme o grafo
    j=0;
    for(i=0;i<numeroBarras;i++){
        aux = (double) grafo[i].idNo;
        aux += 0.01;
        switch (grafo[i].fases){
            case 1:
                regua[j] = aux;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 2:
                regua[j] = aux + 0.1;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 3:
                regua[j] = aux + 0.2;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 4:
                regua[j] = aux;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.1;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 5:
                regua[j] = aux;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.2;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 6:
                regua[j] = aux+0.1;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.2;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 7:
                regua[j] = aux;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.1;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.2;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;    
        }
    }
    aux = 0;
    
    //Tratamento da referência
    long int ref_1, ref_2;
    tratamento_referencia(&ref_1, &ref_2, &alimentadores[0], regua, nvar); 
      
    
    //vetor h aponta para a estrutura de dados das medidas
    for(i=0;i<nmed;i++){
        h[i] = &medidas[i].h;
    }
    //Matriz H aponta para a estrutura de dados das medidas
    for(i=0;i<nmed;i++){
        for(j=0;j<medidas[i].nvar;j++){
            for(r = 0;r<nvar;r++){
                if (cabs(medidas[i].reguaH[j]-regua[r]) < EPS){
                    H[i][r] = &medidas[i].H[j];
                    break;
                }
            }
        }
    }  
    //--------------------------------------------------------------------------
    //Estimação de Estado - Método de Newton Raphson  
    monta_z(z,nmed,medidas); // para o fluxo somente valores das barras PQ, PV e VTeta devem compor o vetor de medidas
    monta_W_Ident(NULL,nmed,medidas);
    
    incializa_vetor_x(grafo, numeroBarras, alimentadores, numeroAlimentadores,x,regua,nvar);
    //incializa_vetor_x_leitura(grafo, numeroBarras, alimentadores, numeroAlimentadores,x,regua,nvar); //Função que le de arquivo externo a condição inicial
    
    double tol = 0.000001;
    clock_t tIni = clock();
    int conv = otimiza_Gauss_NewtonQR(z, h, H, grafo, numeroBarras, ramos, medidas, nvar-3, nmed, regua, x, tol,ref_1,ref_2);
   
    clock_t t1 = clock();
    double tempoWLS = (double)(t1-tIni)/CLOCKS_PER_SEC;
    printf("\nFluxo de Potencia NR-QR: %lf",tempoWLS);
    
    exportaCasoReferencia(grafo, numeroBarras, Sbase);
    exportaEstado(grafo,regua,nvar);
    imprimeEstado(grafo,numeroBarras);
        
    FILE *residuo;
    residuo = fopen("residuo.txt","wt");
    for(i=0;i<nmed;i++){
        fprintf(residuo,"%.10lf\n",z[i] - medidas[i].h);
    }
    fclose(residuo);
    
    free(z);free(h);free(x);free(regua);
    for (i=0;i<nmed;i++ )free(H[i]);
    free(H);
}
void estimadorNEC(GRAFO *grafo, long int numeroBarras, DMED *medidas, DMED *virtuais, long int **numeroMedidas, long int **numeroVirtuais, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos,double Sbase){
   long int nmed,nvar, nvir;
    int i,j,k, r;
    double *z = NULL,**h = NULL,***H = NULL, ***C = NULL, **W = NULL, *x = NULL, *regua = NULL, aux = 0;
    
    printf("Estimador de Estado NEC Trifásico...\n");
    //--------------------------------------------------------------------------
    //Alocação de memória das variáveis do estimador de estado
    nmed = 0;
    for (i = 0; i < 9; i++){ 
        for (j = 0; j < 8; j++){
            nmed = nmed + numeroMedidas[i][j];
        }
    }
    
    nvir = 0;
    for (i = 0; i < 9; i++){
        for (j = 0; j < 8; j++){
            nvir = nvir + numeroVirtuais[i][j];
        }
    }
    
    nvar = 0;
    //printf("numero barras: %d\n", numeroBarras);
    for (i = 0; i < numeroBarras; i++){
        switch (grafo[i].fases){
            case 1:
                nvar +=2;
                break;
            case 2:
                nvar +=2;
                break;
            case 3:
                nvar +=2;
                break;
            case 4:
                nvar +=4;
                break;    
            case 5:
                nvar +=4;
                break;    
            case 6:
                nvar +=4;
                break;    
            case 7:
                nvar +=6;
                break;    
        }
    } 
     
    //printf("nmed: %d\n", nmed);
    //printf("nvar: %d\n", nvar);
    if ((z = (double *)malloc( (nmed) * sizeof(double)))==NULL){
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor z!!!!");
        exit(1); 
    }
    if ((h = malloc( (nmed) * sizeof(double*)))==NULL){
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor h!!!!");
        exit(1); 
    }
    if ((x = (double *)malloc( (nvar+nvir) * sizeof(double)))==NULL){
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor x!!!!");
        exit(1); 
    }
    if ((regua = (double *)malloc( (nvar) * sizeof(double)))==NULL){
        printf("Erro -- Nao foi possivel alocar espaco de memoria para o vetor regua!!!!");
        exit(1); 
    }
    
    H = (double***)malloc(nmed * sizeof(double**)); 
    for (i = 0; i < nmed; i++){ 
         H[i] = (double**) malloc(nvar * sizeof(double*));
         for (j = 0; j < nvar; j++){
              H[i][j] = &aux;
         }
    }
    
    C = (double***)malloc(nvir * sizeof(double**));
    for (i = 0; i < nvir; i ++){
        C[i] = (double**)malloc(nvar * sizeof(double*));
        for (j = 0; j < nvar; j++){
            C[i][j] = &aux;
        }
    }
    
    j=0;
    for(i=0;i<numeroBarras;i++){
        aux = (double) grafo[i].idNo;
        aux += 0.01;
        switch (grafo[i].fases){
            case 1:
                regua[j] = aux;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 2:
                regua[j] = aux + 0.1;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 3:
                regua[j] = aux + 0.2;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 4:
                regua[j] = aux;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.1;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 5:
                regua[j] = aux;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.2;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 6:
                regua[j] = aux+0.1;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.2;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;
            case 7:
                regua[j] = aux;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.1;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                regua[j] = aux + 0.2;
                regua[j + (int) nvar/2] = -regua[j];
                j++;
                break;    
        }
    }
    aux = 0;
    
    //printf("nmed: %d\n", nmed);
    //printf("nvar: %d\n", nvar);
    //Tratamento da referência
    long int ref_1, ref_2;
    tratamento_referencia(&ref_1, &ref_2, &alimentadores[0], regua, nvar);
   
    
    tira_refs_regua(nvar, ref_1, ref_2, regua); 
   
    nvar = nvar - (ref_2 - ref_1 +1);  
    //printf("tira refs\n");
    //vetor h aponta para a estrutura de dados das medidas
    for(i=0;i<nmed;i++){
        h[i] = &medidas[i].h;
    }
    
    //Matriz H aponta para a estrutura de dados das medidas
    for(i=0;i<nmed;i++){
        for(j=0;j<medidas[i].nvar;j++){
            for(r = 0;r<nvar;r++){
                if (cabs(medidas[i].reguaH[j]-regua[r]) < EPS){
                    H[i][r] = &medidas[i].H[j];
                    break;
                }
            }
        }
    }
    

    for(i=0;i<nvir;i++){
        for(j=0;j<virtuais[i].nvar;j++){
            for(r = 0;r<nvar;r++){
                if (cabs(virtuais[i].reguaH[j]-regua[r]) < EPS){
                    C[i][r] = &virtuais[i].H[j];
                    break;
                }
            }
        }
    }
    
    monta_z_comVirtuais(z, nmed, nvir, medidas, virtuais);
    //monta W -> H'WH
    //monta_W_Ident(NULL, nmed, medidas);
    monta_W(NULL, nmed, medidas);
    //monta_W_cte(NULL, nmed, medidas);
    //inicializa primeiros nvar valores do vetor x

    
//    int conv = otimiza_Gauss_Newton_sparsePCGLS(z, h, H, W, grafo, numeroBarras, ramos, medidas, nvar - 3, nmed, regua, x, tol, ref_1, ref_2)
    

    incializa_vetor_x(grafo, numeroBarras, alimentadores, numeroAlimentadores,x,regua,nvar);
    double tol = 0.000001;
    clock_t tIni = clock();
    int conv = otimizaNEC(z, h, H, C, grafo, numeroBarras, ramos, medidas, virtuais, nvir, nvar, nmed, regua, x, tol, ref_1, ref_2);
    clock_t t1 = clock();
    double tempoWLS = (double)(t1-tIni)/CLOCKS_PER_SEC;
    printf("\nEstimação NEC: %lf",tempoWLS);
    
    BOOL printV = true;
    if (printV){
        printf("\nTensoes Nodais (p.u.):\n");
    for(i=0; i<numeroBarras; i++){ 
        //Retangulares
        //printf("%d\tVa: %.5lf + j%.5lf\tVb: %.5lf + j%.5lf\tVc: %.5lf + j%.5lf\n",grafo[i].barra->ID,__real__ grafo[i].V[0],__imag__ grafo[i].V[0],__real__ grafo[i].V[1],__imag__ grafo[i].V[1],__real__ grafo[i].V[2],__imag__ grafo[i].V[2]);
        //Polares
        switch (grafo[i].fases){
            case 1:
                printf("%d\tVa: %.5lf | %.3lf \tVb:    -    |    -   \tVc:    -    |    -   \n",grafo[i].barra->ID,cabs(grafo[i].V[0]),carg(grafo[i].V[0])*180/PI);
                break;
            case 2:
                printf("%d\tVa:    -    |    -    \tVb: %.5lf | %.3lf\tVc:    -    |    -   \n",grafo[i].barra->ID,cabs(grafo[i].V[1]),carg(grafo[i].V[1])*180/PI);
                break;
            case 3:
                printf("%d\tVa:    -    |    -    \tVb:    -    |    -   \tVc: %.5lf | %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[2]),carg(grafo[i].V[2])*180/PI);
                break;
            case 4:
                printf("%d\tVa: %.5lf | %.3lf \tVb: %.5lf | %.3lf\tVc:    -    |    -   \n",grafo[i].barra->ID,cabs(grafo[i].V[0]),carg(grafo[i].V[0])*180/PI,cabs(grafo[i].V[1]),carg(grafo[i].V[1])*180/PI);
                break;
            case 5:
                printf("%d\tVa: %.5lf | %.3lf \tVb:    -    |    -   \tVc: %.5lf | %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[0]),carg(grafo[i].V[0])*180/PI,cabs(grafo[i].V[2]),carg(grafo[i].V[2])*180/PI);
                break;
            case 6:
                printf("%d\tVa:    -    |    -    \tVb: %.5lf | %.3lf\tVc: %.5lf | %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[1]),carg(grafo[i].V[1])*180/PI,cabs(grafo[i].V[2]),carg(grafo[i].V[2])*180/PI);
                break;
            case 7:
                printf("%d\tVa: %.5lf | %.3lf \tVb: %.5lf | %.3lf\tVc: %.5lf | %.3lf\n",grafo[i].barra->ID,cabs(grafo[i].V[0]),carg(grafo[i].V[0])*180/PI,cabs(grafo[i].V[1]),carg(grafo[i].V[1])*180/PI,cabs(grafo[i].V[2]),carg(grafo[i].V[2])*180/PI);
                break;    
        }
    }
    }
    
    atualiza_H_ret(grafo, numeroBarras, ramos, medidas, nmed);
    atualiza_H_ret(grafo, numeroBarras, ramos, virtuais, nvir);
    
}
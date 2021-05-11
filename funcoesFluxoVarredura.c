#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <complex.h>

#include "data_structures.h"
#include "funcoesTopologia.h"
#include "funcoesFluxoVarredura.h"
#include "funcoesCalculoEletrico.h"
#include "funcoesMatematicas.h"

//Imprime na tela as tensões nodais trifásicas
void imprimeTensoesNodais(GRAFO *grafo){
    int i;
    
//    printf("Vbase:%.3lf\t",grafo->Vbase);
    switch (grafo->fases){
        case 1:
            printf("%d\tVa: %.5lf | %.3lf \tVb:    -    |    -   \tVc:    -    |    -   \n",grafo->barra->ID,cabs(grafo->V[0]),carg(grafo->V[0])*180/PI);
            break;
        case 2:
            printf("%d\tVa:    -    |    -    \tVb: %.5lf | %.3lf\tVc:    -    |    -   \n",grafo->barra->ID,cabs(grafo->V[1]),carg(grafo->V[1])*180/PI);
            break;
        case 3:
            printf("%d\tVa:    -    |    -    \tVb:    -    |    -   \tVc: %.5lf | %.3lf\n",grafo->barra->ID,cabs(grafo->V[2]),carg(grafo->V[2])*180/PI);
            break;
        case 4:
            printf("%d\tVa: %.5lf | %.3lf \tVb: %.5lf | %.3lf\tVc:    -    |    -   \n",grafo->barra->ID,cabs(grafo->V[0]),carg(grafo->V[0])*180/PI,cabs(grafo->V[1]),carg(grafo->V[1])*180/PI);
            break;
        case 5:
            printf("%d\tVa: %.5lf | %.3lf \tVb:    -    |    -   \tVc: %.5lf | %.3lf\n",grafo->barra->ID,cabs(grafo->V[0]),carg(grafo->V[0])*180/PI,cabs(grafo->V[2]),carg(grafo->V[2])*180/PI);
            break;
        case 6:
            printf("%d\tVa:    -    |    -    \tVb: %.5lf | %.3lf\tVc: %.5lf | %.3lf\n",grafo->barra->ID,cabs(grafo->V[1]),carg(grafo->V[1])*180/PI,cabs(grafo->V[2]),carg(grafo->V[2])*180/PI);
            break;
        case 7:
            printf("%d\tVa: %.5lf | %.3lf \tVb: %.5lf | %.3lf\tVc: %.5lf | %.3lf\n",grafo->barra->ID,cabs(grafo->V[0]),carg(grafo->V[0])*180/PI,cabs(grafo->V[1]),carg(grafo->V[1])*180/PI,cabs(grafo->V[2]),carg(grafo->V[2])*180/PI);
            break;    
    }
    
}
//Imprime na tela os taps dos reguladores
void imprimeTaps(GRAFO *grafo){
    int i;
    for (i=0;i < grafo->numeroAdjacentes;i++){
        if((grafo->adjacentes[i].tipo  == regulador) && (grafo->adjacentes[i].ramo->k == grafo->idNo)){
            int idram = grafo->adjacentes[i].idram;
            printf("Regulador Barra %d:    ",grafo->adjacentes[i].ramo->DE);
            printf("%.3lf kVA\t%.2lf kV\tControl {%d}\t%.0lf\t%.0lf\t%.0lf\n",grafo->adjacentes[i].ramo->regulador.Snominal,grafo->adjacentes[i].ramo->regulador.Vnom/1000,grafo->adjacentes[i].ramo->regulador.controle,grafo->adjacentes[i].ramo->regulador.tap[0],grafo->adjacentes[i].ramo->regulador.tap[1],grafo->adjacentes[i].ramo->regulador.tap[2]);
        }
    }
}
//Imprime na tela as tensões nodais trifásicas
void imprimeiInjecoesCorrentes(GRAFO *grafo){
    int i;
    
//    printf("Vbase:%.3lf\t",grafo->Vbase);
    switch (grafo->fases){
        case 1:
            printf("%d\tIa: %.5lf | %.3lf \tIb:    -    |    -   \tIc:    -    |    -   \n",grafo->barra->ID,cabs(grafo->Cur[0]),carg(grafo->Cur[0])*180/PI);
            break;
        case 2:
            printf("%d\tIa:    -    |    -    \tIb: %.5lf | %.3lf\tIc:    -    |    -   \n",grafo->barra->ID,cabs(grafo->Cur[1]),carg(grafo->Cur[1])*180/PI);
            break;
        case 3:
            printf("%d\tIa:    -    |    -    \tIb:    -    |    -   \tIc: %.5lf | %.3lf\n",grafo->barra->ID,cabs(grafo->Cur[2]),carg(grafo->Cur[2])*180/PI);
            break;
        case 4:
            printf("%d\tIa: %.5lf | %.3lf \tIb: %.5lf | %.3lf\tIc:    -    |    -   \n",grafo->barra->ID,cabs(grafo->Cur[0]),carg(grafo->Cur[0])*180/PI,cabs(grafo->Cur[1]),carg(grafo->Cur[1])*180/PI);
            break;
        case 5:
            printf("%d\tIa: %.5lf | %.3lf \tIb:    -    |    -   \tIc: %.5lf | %.3lf\n",grafo->barra->ID,cabs(grafo->Cur[0]),carg(grafo->Cur[0])*180/PI,cabs(grafo->Cur[2]),carg(grafo->Cur[2])*180/PI);
            break;
        case 6:
            printf("%d\tIa:    -    |    -    \tIb: %.5lf | %.3lf\tIc: %.5lf | %.3lf\n",grafo->barra->ID,cabs(grafo->Cur[1]),carg(grafo->Cur[1])*180/PI,cabs(grafo->Cur[2]),carg(grafo->Cur[2])*180/PI);
            break;
        case 7:
            printf("%d\tIa: %.5lf | %.3lf \tIb: %.5lf | %.3lf\tIc: %.5lf | %.3lf\n",grafo->barra->ID,cabs(grafo->Cur[0]),carg(grafo->Cur[0])*180/PI,cabs(grafo->Cur[1]),carg(grafo->Cur[1])*180/PI,cabs(grafo->Cur[2]),carg(grafo->Cur[2])*180/PI);
            break;    
    }
    
}

//Imprime na tela as tensões nodais trifásicas
void imprimeCorrentes(NOADJACENTE *noAdj){
   int i;
   
//    printf("Vbase:%.3lf\t",grafo->Vbase);
   switch (noAdj->ramo->fases){
       case 1:
           printf("%d - %d\tIa: %.5lf | %.3lf \tIb:    -    |    -   \tIc:    -    |    -   \n",noAdj->ramo->DE,noAdj->ramo->PARA,cabs(noAdj->Cur[0]),carg(noAdj->Cur[0])*180/PI);
           break;
       case 2:
           printf("%d - %d\tIa:    -    |    -    \tIb: %.5lf | %.3lf\tIc:    -    |    -   \n",noAdj->ramo->DE,noAdj->ramo->PARA,cabs(noAdj->Cur[1]),carg(noAdj->Cur[1])*180/PI);
           break;
       case 3:
           printf("%d - %d\tIa:    -    |    -    \tIb:    -    |    -   \tIc: %.5lf | %.3lf\n",noAdj->ramo->DE,noAdj->ramo->PARA,cabs(noAdj->Cur[2]),carg(noAdj->Cur[2])*180/PI);
           break;
       case 4:
           printf("%d - %d\tIa: %.5lf | %.3lf \tIb: %.5lf | %.3lf\tIc:    -    |    -   \n",noAdj->ramo->DE,noAdj->ramo->PARA,cabs(noAdj->Cur[0]),carg(noAdj->Cur[0])*180/PI,cabs(noAdj->Cur[1]),carg(noAdj->Cur[1])*180/PI);
           break;
       case 5:
           printf("%d - %d\tIa: %.5lf | %.3lf \tIb:    -    |    -   \tIc: %.5lf | %.3lf\n",noAdj->ramo->DE,noAdj->ramo->PARA,cabs(noAdj->Cur[0]),carg(noAdj->Cur[0])*180/PI,cabs(noAdj->Cur[2]),carg(noAdj->Cur[2])*180/PI);
           break;
       case 6:
           printf("%d - %d\tIa:    -    |    -    \tIb: %.5lf | %.3lf\tIc: %.5lf | %.3lf\n",noAdj->ramo->DE,noAdj->ramo->PARA,cabs(noAdj->Cur[1]),carg(noAdj->Cur[1])*180/PI,cabs(noAdj->Cur[2]),carg(noAdj->Cur[2])*180/PI);
           break;
       case 7:
           printf("%d - %d\tIa: %.5lf | %.3lf \tIb: %.5lf | %.3lf\tIc: %.5lf | %.3lf\n",noAdj->ramo->DE,noAdj->ramo->PARA,cabs(noAdj->Cur[0]),carg(noAdj->Cur[0])*180/PI,cabs(noAdj->Cur[1]),carg(noAdj->Cur[1])*180/PI,cabs(noAdj->Cur[2]),carg(noAdj->Cur[2])*180/PI);
           break;    
   }
   
}


//Função inicialização do vetor x para sistemas radiais (sem o arquivo Vinicial)
void incializaTensoesRaiz(GRAFO *grafo, long int numeroBarras, ALIMENTADOR *alimentadores, long int numeroAlimentadores){ 
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
        while(barraAtual != NULL)
        {
            grafo[barraAtual->idNo].V[0] = V0[0];
            grafo[barraAtual->idNo].V[1] = V0[1];
            grafo[barraAtual->idNo].V[2] = V0[2];
            
            
            visitado[barraAtual->idNo] = true;
            barraAtual = barraAtual->prox;
        }        
    }
}

//Função inicialização do vetor x para sistemas radiais (sem o arquivo Vinicial)
void incializaTensoesVarredura(GRAFO *grafo, long int numeroBarras, ALIMENTADOR *alimentadores, long int numeroAlimentadores){ 
    int i, k, fase;
    BOOL visitado[numeroBarras];
    __complex__ double V0[3], **Yaux;
    BOOL control_REG_OPT;

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
        visitado[de] = true;
        
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
                            grafo[para].V[0] = cabs(grafo[de].V[0])*(cos(-30*PI/180) + I*sin(-30*PI/180));
                            grafo[para].V[1] = cabs(grafo[de].V[1])*(cos(-150*PI/180) + I*sin(-150*PI/180));
                            grafo[para].V[2] = cabs(grafo[de].V[2])*(cos(90*PI/180) + I*sin(90*PI/180));                                    
                        }            
                                    
                
                    }
                    else if (grafo[de].adjacentes[k].tipo == 2){ //Para o caso de regulador de tensão
                        if (control_REG_OPT == 1){
                            grafo[de].adjacentes[k].ramo->regulador.tap[0] = 0;
                            grafo[de].adjacentes[k].ramo->regulador.tap[1] = 0;
                            grafo[de].adjacentes[k].ramo->regulador.tap[2] = 0;

                            atualizaTapRegulador(grafo[de].adjacentes[k].ramo);
                        }

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
    for (i=0;i<3;i++) free(Yaux[i]);
    free(Yaux);
    
}

void atualizaTapRegulador(DRAM *ramo){
    if (ramo->tipo == 2){
        // Tira efeito dos antigos valores do tap nas matrizes de quadripolo
            ramo->Ypp[0][0] = 1/pow(ramo->tap_pri[0],2)*ramo->Ypp[0][0];
            ramo->Ypp[1][1] = 1/pow(ramo->tap_pri[1],2)*ramo->Ypp[1][1];
            ramo->Ypp[2][2] = 1/pow(ramo->tap_pri[2],2)*ramo->Ypp[2][2];
            
            ramo->Yps[0][0] = 1/ramo->tap_sec[0]*1/ramo->tap_pri[0]*ramo->Yps[0][0];
            ramo->Yps[1][1] = 1/ramo->tap_sec[1]*1/ramo->tap_pri[1]*ramo->Yps[1][1];
            ramo->Yps[2][2] = 1/ramo->tap_sec[2]*1/ramo->tap_pri[2]*ramo->Yps[2][2];
            
            ramo->Ysp[0][0] = 1/ramo->tap_sec[0]*1/ramo->tap_pri[0]*ramo->Ysp[0][0];
            ramo->Ysp[1][1] = 1/ramo->tap_sec[1]*1/ramo->tap_pri[1]*ramo->Ysp[1][1];
            ramo->Ysp[2][2] = 1/ramo->tap_sec[2]*1/ramo->tap_pri[2]*ramo->Ysp[2][2];
            
            ramo->Yss[0][0] = 1/pow(ramo->tap_sec[0],2)*ramo->Yss[0][0];
            ramo->Yss[1][1] = 1/pow(ramo->tap_sec[1],2)*ramo->Yss[1][1];
            ramo->Yss[2][2] = 1/pow(ramo->tap_sec[2],2)*ramo->Yss[2][2];
    
    // Atualiaza os taps
            ramo->tap_pri[0] = (1 + ramo->regulador.tap[0]*ramo->regulador.regulacao/ramo->regulador.ntaps);
            ramo->tap_pri[1] = (1 + ramo->regulador.tap[1]*ramo->regulador.regulacao/ramo->regulador.ntaps);
            ramo->tap_pri[2] = (1 + ramo->regulador.tap[2]*ramo->regulador.regulacao/ramo->regulador.ntaps);
            
            ramo->tap_sec[0] = 1;
            ramo->tap_sec[1] = 1;
            ramo->tap_sec[2] = 1;
            
            ramo->Ypp[0][0] = pow(ramo->tap_pri[0],2)*ramo->Ypp[0][0];
            ramo->Ypp[1][1] = pow(ramo->tap_pri[1],2)*ramo->Ypp[1][1];
            ramo->Ypp[2][2] = pow(ramo->tap_pri[2],2)*ramo->Ypp[2][2];
            
            ramo->Yps[0][0] = ramo->tap_sec[0]*ramo->tap_pri[0]*ramo->Yps[0][0];
            ramo->Yps[1][1] = ramo->tap_sec[1]*ramo->tap_pri[1]*ramo->Yps[1][1];
            ramo->Yps[2][2] = ramo->tap_sec[2]*ramo->tap_pri[2]*ramo->Yps[2][2];
            
            ramo->Ysp[0][0] = ramo->tap_sec[0]*ramo->tap_pri[0]*ramo->Ysp[0][0];
            ramo->Ysp[1][1] = ramo->tap_sec[1]*ramo->tap_pri[1]*ramo->Ysp[1][1];
            ramo->Ysp[2][2] = ramo->tap_sec[2]*ramo->tap_pri[2]*ramo->Ysp[2][2];
            
            ramo->Yss[0][0] = pow(ramo->tap_sec[0],2)*ramo->Yss[0][0];
            ramo->Yss[1][1] = pow(ramo->tap_sec[1],2)*ramo->Yss[1][1];
            ramo->Yss[2][2] = pow(ramo->tap_sec[2],2)*ramo->Yss[2][2];
        }
}

//Função que atualiza as injeções de corrente nas barras da rede conforme as cargas e e respectivos modelos de carga, shunts e geradores distribuidos
void atualizaInjecoes(GRAFO *no){
    int i, j;
    __complex__ double Saux[3], Iaux[3], IauxDelta[3], Vl[3], V0;

    no->Cur[0] = 0;
    no->Cur[1] = 0;
    no->Cur[2] = 0;
    no->S[0] = 0;
    no->S[1] = 0;
    no->S[2] = 0;
    
    // Atualiza Cargas
    for (i = 0; i<no->barra->nloads;i++){
        switch (no->barra->loads[i].lig){
            case YN:
                Saux[0] = (no->barra->loads[i].Pnom[0] + I * no->barra->loads[i].Qnom[0]) * pow(cabs(no->V[0]) , no->barra->loads[i].ZIP);
                Saux[1] = (no->barra->loads[i].Pnom[1] + I * no->barra->loads[i].Qnom[1]) * pow(cabs(no->V[1]) , no->barra->loads[i].ZIP);
                Saux[2] = (no->barra->loads[i].Pnom[2] + I * no->barra->loads[i].Qnom[2]) * pow(cabs(no->V[2]) , no->barra->loads[i].ZIP);
                
                Iaux[0] = conj(Saux[0]/no->V[0]);
                Iaux[1] = conj(Saux[1]/no->V[1]);
                Iaux[2] = conj(Saux[2]/no->V[2]);
                
                break;
            case Y:
                V0 = no->V[0] + no->V[1] + no->V[2]; //tensão de neutro
                Saux[0] = (no->barra->loads[i].Pnom[0] + I * no->barra->loads[i].Qnom[0]) * pow(cabs(no->V[0] - V0) , no->barra->loads[i].ZIP);
                Saux[1] = (no->barra->loads[i].Pnom[1] + I * no->barra->loads[i].Qnom[1]) * pow(cabs(no->V[1] - V0) , no->barra->loads[i].ZIP);
                Saux[2] = (no->barra->loads[i].Pnom[2] + I * no->barra->loads[i].Qnom[2]) * pow(cabs(no->V[2] - V0) , no->barra->loads[i].ZIP);
                
                Iaux[0] = conj(Saux[0]/(no->V[0] - V0));
                Iaux[1] = conj(Saux[1]/(no->V[1] - V0));
                Iaux[2] = conj(Saux[2]/(no->V[2] - V0));
                
                break;
            case D:
                tensaoDelta(no->V, Vl);
                Saux[0] = (no->barra->loads[i].Pnom[0] + I * no->barra->loads[i].Qnom[0]) * pow(cabs(no->V[0]) , no->barra->loads[i].ZIP);
                Saux[1] = (no->barra->loads[i].Pnom[1] + I * no->barra->loads[i].Qnom[1]) * pow(cabs(no->V[1]) , no->barra->loads[i].ZIP);
                Saux[2] = (no->barra->loads[i].Pnom[2] + I * no->barra->loads[i].Qnom[2]) * pow(cabs(no->V[2]) , no->barra->loads[i].ZIP);
                
                IauxDelta[0] = conj(Saux[0]/Vl[0]);
                IauxDelta[1] = conj(Saux[1]/Vl[1]);
                IauxDelta[2] = conj(Saux[2]/Vl[2]);
                
                Iaux[0] = IauxDelta[0] - IauxDelta[2];
                Iaux[1] = IauxDelta[1] - IauxDelta[0];
                Iaux[2] = IauxDelta[2] - IauxDelta[1];
                
                break;
        }
        no->Cur[0] += Iaux[0];
        no->Cur[1] += Iaux[1];
        no->Cur[2] += Iaux[2];
        
        no->S[0] += no->V[0] * conj(Iaux[0]);
        no->S[1] += no->V[1] * conj(Iaux[1]);
        no->S[2] += no->V[2] * conj(Iaux[2]);
    } 
    
    
    // Atualiza Shunts
    for (i = 0; i<no->barra->nshunts;i++){
        switch (no->barra->shunts[i].lig){
            case YN:
                Saux[0] = ( I * no->barra->shunts[i].Qnom[0]) * pow(cabs(no->V[0]) , 2);
                Saux[1] = ( I * no->barra->shunts[i].Qnom[1]) * pow(cabs(no->V[1]) , 2);
                Saux[2] = ( I * no->barra->shunts[i].Qnom[2]) * pow(cabs(no->V[2]) , 2);
                
                Iaux[0] = conj(Saux[0]/no->V[0]);
                Iaux[1] = conj(Saux[1]/no->V[1]);
                Iaux[2] = conj(Saux[2]/no->V[2]);
                
                break;
            case Y:
                V0 = no->V[0] + no->V[1] + no->V[2]; //tensão de neutro
                Saux[0] = ( I * no->barra->shunts[i].Qnom[0]) * pow(cabs(no->V[0] - V0) , 2);
                Saux[1] = ( I * no->barra->shunts[i].Qnom[1]) * pow(cabs(no->V[1] - V0) , 2);
                Saux[2] = ( I * no->barra->shunts[i].Qnom[2]) * pow(cabs(no->V[2] - V0) , 2);
                
                Iaux[0] = conj(Saux[0]/(no->V[0] - V0));
                Iaux[1] = conj(Saux[1]/(no->V[1] - V0));
                Iaux[2] = conj(Saux[2]/(no->V[2] - V0));
                
                break;
            case D:
                tensaoDelta(no->V, Vl);
                Saux[0] = ( I * no->barra->shunts[i].Qnom[0]) * pow(cabs(no->V[0]) , 2);
                Saux[1] = ( I * no->barra->shunts[i].Qnom[1]) * pow(cabs(no->V[1]) , 2);
                Saux[2] = ( I * no->barra->shunts[i].Qnom[2]) * pow(cabs(no->V[2]) , 2);
                
                IauxDelta[0] = conj(Saux[0]/Vl[0]);
                IauxDelta[1] = conj(Saux[1]/Vl[1]);
                IauxDelta[2] = conj(Saux[2]/Vl[2]);
                
                Iaux[0] = IauxDelta[0] - IauxDelta[2];
                Iaux[1] = IauxDelta[1] - IauxDelta[0];
                Iaux[2] = IauxDelta[2] - IauxDelta[1];
                
                break;
        }
        no->Cur[0] += Iaux[0];
        no->Cur[1] += Iaux[1];
        no->Cur[2] += Iaux[2];
        
        no->S[0] += no->V[0] * conj(Iaux[0]);
        no->S[1] += no->V[1] * conj(Iaux[1]);
        no->S[2] += no->V[2] * conj(Iaux[2]);
    } 
    
    // Atualiza Geradores Distribuídos
    for (i = 0; i<no->barra->ngds;i++){
        //Futuro para geração distribuída------------------------------------------****************************
        Iaux[0] = 0;
        Iaux[1] = 0;
        Iaux[2] = 0;
        
        no->Cur[0] += Iaux[0];
        no->Cur[1] += Iaux[1];
        no->Cur[2] += Iaux[2];
        
        no->S[0] += no->V[0] * conj(Iaux[0]);
        no->S[1] += no->V[0] * conj(Iaux[1]);
        no->S[2] += no->V[0] * conj(Iaux[2]);
    }
}



//cálculo de corrente a montante no quadripolo
void calculaCorrenteMontante(){
    
}

//cálculo de tensão a jusante no quadripolo
void calculaTensaoJusante(complex double *Vp, complex double *Vs, complex double *Ips, DRAM *ramo){
    
    int i,j;
    __complex__ double Vaux[3], Aux[3], **Y, **Yp, **Ys, Iaux[3], V0;
    BOOL singular1 = false, singular2 = false;
    
    Y = c_matAloca(3);
    Yp = c_matAloca(3);
    Ys = c_matAloca(3);
    
    Vaux[0]=0;
    Vaux[1]=0;
    Vaux[2]=0;
    
    //Verifica posto das submatrizes do quadripolo de acordo com a ligacao
    if(ramo->tipo == 1){
        switch (ramo->trafo.lig_pri){
            case 3:
            case 2:
                singular1 = true;
                break;
        }
        switch (ramo->trafo.lig_sec){
            case 2:
                singular2 = true;
                break;
        }

    }
    

    switch (ramo->tipo){
        case 0: //Linha
            //Calcula queda de tensão na linha através da impedância
            for(i=0;i<3;i++){
                for(j=0;j<3;j++){
                    Vaux[i] += ramo->B[i][j] * Vp[j];
                }
                Aux[i] = Ips[i] - Vaux[i];
            }
            for(i=0;i<3;i++){
                Vaux[i] = 0;
                for(j=0;j<3;j++){
                    Vaux[i] += ramo->Z[i][j] * Aux[j];
                }
                Vs[i] = Vp[i] - Vaux[i];
            }


            break;
        //----------------------------------------------------------------------
        case 1: //Trafo - está para o YN-YN
            V0 = 0;
            Iaux[0] = Ips[0];
            Iaux[1] = Ips[1];
            Iaux[2] = Ips[2];

            c_matIgual(Y, ramo->Yps, 3);
            c_matIgual(Yp, ramo->Ypp, 3);
            
            //tratamento de excessões das ligações de trafos com matriz de admitância singulares
            //--------------------------------------
            //Caso a matriz de quadripolo é singular 
            if (singular1) {
                //Calcula a tensão de neutro do secundário
                c_matIgual(Ys, ramo->Yss, 3);
                c_matInversaZ(Ys, 3);
                for(i=0;i<3;i++){
                    Vaux[i] = 0;
                    for(j=0;j<3;j++){
                        Vaux[i] += ramo->Ysp[i][j] * Vp[j];
                    }
                    Aux[i] = - Vaux[i];
                    //  Aux[i] = Isp[i] - Vaux[i]; // Falta pensar na corrente do secundário!
                }
                for(i=0;i<3;i++){
                    Vs[i] = 0;
                    for(j=0;j<3;j++){
                        Vs[i] += Ys[i][j] * Aux[j];
                    }
                    V0 += Vs[i]/3;
                }
                                
                for (i=0;i<3;i++){
                    Y[2][i] = 1;
                    Yp[2][i] = 0;
                }
                Iaux[2] = 0;
            }
            else if (singular2){
                V0 = 0; //Aproximação para caso de ligação em Delta no secundário - ver referência (Tratativa específica para calcular V0 no secundario - sistemas multiaterrados)
                for (i=0;i<3;i++){
                    Y[2][i] = 1;
                    Yp[2][i] = 0;
                }
                Iaux[2] = 0;
            }
                
            //---------------------------------------
            
            //Calcula Vs = Yps ^-1 * (Ips - Ypp * Vp)
            //Calcula a tensão nas três fases
            for(i=0;i<3;i++){
                for(j=0;j<3;j++){
                    Vaux[i] += Yp[i][j]*Vp[j];
                }
                Aux[i] = Iaux[i] - Vaux[i];
            }
            
            c_matInversaZ(Y, 3);
            for(i=0;i<3;i++){
                Vs[i] = 0;
                for(j=0;j<3;j++){
                    Vs[i] += Y[i][j] * Aux[j];
                }
                Vs[i] += V0;
                if (cabs(Vs[i]) == 0) Vs[i] = 1; //Caso em que não tem a fase resulta em V = 0 - coloca 1 para evitar NaN na hora de atualizar a corrente
            }

            break;
        case 2: //Regulador de Tensão
            //Calcula Vs = Yps ^-1 * (Ips - Ypp * Vp)
            //Calcula a tensão nas três fases
            for(i=0;i<3;i++){
                for(j=0;j<3;j++){
                    Vaux[i] += ramo->Ypp[i][j]*Vp[j];
                }
                Aux[i] = Ips[i] - Vaux[i];
            }
            //Inverte matriz Yps - dá para melhorar este ponto computacionalmente
            c_matIgual(Y, ramo->Yps, 3);

            //tratamento de excessões das ligações de trafos com matriz de admitância singulares

            c_matInversaZ(Y, 3);
            for(i=0;i<3;i++){
                Vs[i] = 0;
                for(j=0;j<3;j++){
                    Vs[i] += Y[i][j] * Aux[j];
                }
                if (cabs(Vs[i]) == 0) Vs[i] = 1; //Caso em que não tem a fase resulta em V = 0 - coloca 1 para evitar NaN na hora de atualizar a corrente
            }
            
            
            break;
        case 3: //Chave
            Vs[0] = Vp[0];
            Vs[1] = Vp[1];
            Vs[2] = Vp[2];
            break; 
    }
    for(i = 0;i<3;i++) free(Y[i]);
    free(Y);

    for(i = 0;i<3;i++) free(Yp[i]);
    free(Yp);

    for(i = 0;i<3;i++) free(Ys[i]);
    free(Ys);
    
}


//Função que efetua a etapa backward
void backward(GRAFO *noP, GRAFO *grafo){
    long int i,j,noAdj,noMont = -1, auxNoMont = -1;
    complex double Iacc[3], Iaux[3], Vaux[3], **Yaux;
    BOOL singular = false;
    complex double V0;
    
    Yaux = c_matAloca(3);
    
    //Soma correntes à jusante e identifica o nó a montante
    Iacc[0] = 0;
    Iacc[1] = 0;
    Iacc[2] = 0;
    
    for (i = 0; i < noP->numeroAdjacentes;i ++){
        noAdj = noP->adjacentes[i].idNo;
        
        if (noP->profundidade < grafo[noAdj].profundidade){
            Iacc[0] += noP->adjacentes[i].Cur[0];
            Iacc[1] += noP->adjacentes[i].Cur[1];
            Iacc[2] += noP->adjacentes[i].Cur[2];
        }
        else{
            noMont = noAdj;
            auxNoMont = i;
        }
    }
    if (noMont == -1){ //Nó raiz - corrente é a soma dos jusantes somente
        noP->Cur[0] = -Iacc[0];
        noP->Cur[1] = -Iacc[1];
        noP->Cur[2] = -Iacc[2];
    }
    //Soma injeção de corrente na barra
    Iacc[0] += noP->Cur[0];
    Iacc[1] += noP->Cur[1];
    Iacc[2] += noP->Cur[2];
    
    //--------------------------------------------------------------------------
    //Atualiza a corrente no nó a montante
    if (noMont != -1){
        noP->adjacentes[auxNoMont].Cur[0] = -Iacc[0];
        noP->adjacentes[auxNoMont].Cur[1] = -Iacc[1];
        noP->adjacentes[auxNoMont].Cur[2] = -Iacc[2]; 
        
        //Verifica posto das submatrizes do quadripolo de acordo com a ligacao
        if(noP->adjacentes[auxNoMont].tipo == 1){
            switch (noP->adjacentes[auxNoMont].ramo->trafo.lig_sec){
                case 2:
                case 3:
                    singular = true;
                    break;
            }

        }
        //Calcula a corrente montante de acordo com o tipo de ramo
        switch (noP->adjacentes[auxNoMont].tipo){
            //------------------------------------------------------------------
            case 0: //Linha

                for (i=0;i<3;i++){
                    Iaux[i] = 0;
                    for (j=0;j<3;j++){
                        Iaux[i] += noP->adjacentes[auxNoMont].ramo->B[i][j] * noP->V[j] + noP->adjacentes[auxNoMont].ramo->B[i][j] * grafo[noMont].V[j];
                    }
                }
                
                Iacc[0] = Iacc[0] + Iaux[0];
                Iacc[1] = Iacc[1] + Iaux[1];
                Iacc[2] = Iacc[2] + Iaux[2];
                break;
            //------------------------------------------------------------------
            case 1: //Trafo 
                V0 = 0;
                Iacc[0] = Iacc[0];
                Iacc[1] = Iacc[1];
                Iacc[2] = Iacc[2];
                
                //Injeção de corrente no secundário
                //I = (-Iacc - traf.Yss*Vsec); 
                for (i=0;i<3;i++){
                    Iaux[i] = 0;
                    for (j=0;j<3;j++){
                        Iaux[i] += noP->adjacentes[auxNoMont].ramo->Yss[i][j] * noP->V[j];
                    }
                    Iaux[i] = -Iacc[i] - Iaux[i];
                }
                c_matIgual(Yaux, noP->adjacentes[auxNoMont].ramo->Ysp, 3);
                
                //--------------------------------------
                //Caso a matriz de quadripolo é singular 
                if (singular) {
                    for (i=0;i<3;i++){
                        V0 += grafo[noMont].V[i]/3;
                        Yaux[2][i] = 1;
                    }
                    Iaux[2] = 0;
                }
                //---------------------------------------
                
                //Vpm = traf.Ysp\(-Iacc - traf.Yss*Vsec);
                c_matInversaZ(Yaux, 3);
                for (i=0;i<3;i++){
                    Vaux[i] = 0;
                    for (j=0;j<3;j++){
                        Vaux[i] += Yaux[i][j] * Iaux[j];
                    }
                    Vaux[i] += V0; //quando é Ysp singular soma a tensão de sequência zero - se não for singular V0 = 0
                }
                
                // traf.Ikm = (traf.Ypp*Vpm + traf.Yps*Vsec);
                for (i=0;i<3;i++){
                    Iacc[i] = 0;
                    for (j=0;j<3;j++){
                        Iacc[i] += noP->adjacentes[auxNoMont].ramo->Ypp[i][j] * Vaux[j] + noP->adjacentes[auxNoMont].ramo->Yps[i][j] * noP->V[j];
                    }
                }
                break;
            //------------------------------------------------------------------
            case 2: //Regulador de Tensão
                Iacc[0] = noP->adjacentes[auxNoMont].ramo->tap_pri[0]*Iacc[0];
                Iacc[1] = noP->adjacentes[auxNoMont].ramo->tap_pri[1]*Iacc[1];
                Iacc[2] = noP->adjacentes[auxNoMont].ramo->tap_pri[2]*Iacc[2]; 
                
                break;
            //------------------------------------------------------------------
            case 3: //Chave
                Iacc[0] = Iacc[0];
                Iacc[1] = Iacc[1];
                Iacc[2] = Iacc[2];                
                break; 
        } 
        //Atualiza a corrente a montante no grafo
        for (i = 0; i < grafo[noMont].numeroAdjacentes;i ++){
            if (grafo[noMont].adjacentes[i].idNo == noP->idNo){
                grafo[noMont].adjacentes[i].Cur[0] = Iacc[0];
                grafo[noMont].adjacentes[i].Cur[1] = Iacc[1];
                grafo[noMont].adjacentes[i].Cur[2] = Iacc[2];
            }
            // imprimeCorrentes(&grafo[noMont].adjacentes[i]);
        }
    }

    for(i = 0;i<3;i++) free(Yaux[i]);
    free(Yaux);
}

//Função que atualiza controle de reguladores de tensão de acordo com parametrização do LDC
BOOL controleReguladorTensao_LDC(double Vbase, double Ibase, complex double *Vp, complex double *Vs, complex double *Ips, complex double *Isp, DRAM *ramo){
    int i;
    complex double   Zcont[3], Zcont_r[3];
    BOOL ctr_act;
    double Vdiff[3], taps_0[3], Vset[3], Vset_r[3];
    int Direction[3];

    // printf("\n Regulador %d - %d \n",ramo->DE, ramo->PARA);

    //Parâmetros do controlador LDC
    double reg = ramo->regulador.regulacao / ramo->regulador.ntaps;

    Zcont[0] = ramo->regulador.R1 + I* ramo->regulador.X1;
    Zcont[1] = ramo->regulador.R2 + I* ramo->regulador.X2;
    Zcont[2] = ramo->regulador.R3 + I* ramo->regulador.X3;

    Zcont_r[0] = ramo->regulador.R1r + I* ramo->regulador.X1r;
    Zcont_r[1] = ramo->regulador.R2r + I* ramo->regulador.X2r;
    Zcont_r[2] = ramo->regulador.R3r + I* ramo->regulador.X3r;

    Vset[0] = ramo->regulador.V1;
    Vset[1] = ramo->regulador.V2;
    Vset[2] = ramo->regulador.V3;
    
    Vset_r[0] = ramo->regulador.V1r;
    Vset_r[1] = ramo->regulador.V2r;
    Vset_r[2] = ramo->regulador.V3r;

    for (i=0;i<3;i++) taps_0[i] = ramo->regulador.tap[i];

    //Calcula direção do fluxo
    for (i=0;i<3;i++) {
        if ((creal(Vs[i] * conj(Isp[i])) && (cimag(Vs[i] * conj(Isp[i])) >= 0)) >= 0) Direction[i] = 1;
        else Direction[i] = 0;
    }
    // printf("\n\n");
    //Malha de Controle
    switch (ramo->regulador.controle){
        //------------------------------------------------------
        // CONTROLE SEM RESTRIÇÃO
        case 0 :
            for (i=0;i<3;i++) Vdiff[i] = Vset[i] - cabs(Vbase*Vs[i]/ramo->regulador.TP - Zcont[i] * Ibase / sqrt(3)* ( - Isp[i])/ramo->regulador.TC);

            for (i=0;i<3;i++){
                // printf("\t Vdiff[%d]= %lf",i, Vdiff[i]);
                // printf("\t tap[%d]= %lf",i, ramo->regulador.tap[i]);
                if (Vdiff[i] >= ramo->regulador.deltaV) 
                    ramo->regulador.tap[i] = ramo->regulador.tap[i] + ceil((Vdiff[i] - ramo->regulador.deltaV)/(reg*ramo->regulador.Vnom/pow(3,0.5)/ramo->regulador.TP));
                else if (Vdiff[i] <= -ramo->regulador.deltaV) 
                    ramo->regulador.tap[i] = ramo->regulador.tap[i] + floor((Vdiff[i] + ramo->regulador.deltaV)/(reg*ramo->regulador.Vnom/pow(3,0.5)/ramo->regulador.TP));
                // printf("\t tap[%d]= %lf",i, ramo->regulador.tap[i]);
            }

            break;
        //------------------------------------------------------
        // CONTROLE LOCKED forward_sweep
        case 1 :
            for (i=0;i<3;i++) Vdiff[i] = Vset[i] - cabs(Vbase*Vs[i]/ramo->regulador.TP - Zcont[i] * Ibase * ( - Isp[i])/ramo->regulador.TC);

            for (i=0;i<3;i++){
                if (Direction[i] ==1){
                    if (Vdiff[i] >= ramo->regulador.deltaV) 
                        ramo->regulador.tap[i] = ramo->regulador.tap[i] + ceil((Vdiff[i] - ramo->regulador.deltaV)*(reg*ramo->regulador.Vnom/pow(3,0.5)/ramo->regulador.TP));
                    else if (Vdiff[i] <= -ramo->regulador.deltaV) 
                        ramo->regulador.tap[i] = ramo->regulador.tap[i] + floor((Vdiff[i] + ramo->regulador.deltaV)*(reg*ramo->regulador.Vnom/pow(3,0.5)/ramo->regulador.TP));
                }
                else{
                    ramo->regulador.tap[i] = ramo->regulador.tap_ini[i];
                }
            }
            break;
        //------------------------------------------------------
        // CONTROLE LOCKED REVERSE
        case 2 :
            for (i=0;i<3;i++) Vdiff[i] = Vset_r[i] - cabs(Vbase*Vp[i]/ramo->regulador.TP - Zcont_r[i] * Ibase * ( -Ips[i])/ramo->regulador.TC);

            for (i=0;i<3;i++){
                if (Direction[i] ==1){
                    ramo->regulador.tap[i] = ramo->regulador.tap_ini[i];
                }
                else{
                    if (Vdiff[i] >= ramo->regulador.deltaVr) 
                        ramo->regulador.tap[i] = ramo->regulador.tap[i] - ceil((Vdiff[i] - ramo->regulador.deltaVr)*(reg*ramo->regulador.Vnom/pow(3,0.5)/ramo->regulador.TP));
                    else if (Vdiff[i] <= -ramo->regulador.deltaVr) 
                        ramo->regulador.tap[i] = ramo->regulador.tap[i] - floor((Vdiff[i] + ramo->regulador.deltaVr)*(reg*ramo->regulador.Vnom/pow(3,0.5)/ramo->regulador.TP));
                }
            }
            break;
        //------------------------------------------------------
        // CONTROLE BIDIRECTIONAL 
        case 3:
            for (i=0;i<3;i++){
                if (Direction[i] == 1)
                    Vdiff[i] = Vset[i] - cabs(Vs[i]/ramo->regulador.TP - Zcont[i] * Isp[i]/ramo->regulador.TC);
                else
                    Vdiff[i] = Vset_r[i] - cabs(Vp[i]/ramo->regulador.TP - Zcont_r[i] * Ips[i]/ramo->regulador.TC);                    
            }

            for (i=0;i<3;i++){
                if (Direction[i] ==1){
                    if (Vdiff[i] >= ramo->regulador.deltaV) 
                        ramo->regulador.tap[i] = ramo->regulador.tap[i] + ceil((Vdiff[i] - ramo->regulador.deltaV)*(reg*ramo->regulador.Vnom/pow(3,0.5)/ramo->regulador.TP));
                    else if (Vdiff[i] <= -ramo->regulador.deltaV) 
                        ramo->regulador.tap[i] = ramo->regulador.tap[i] + floor((Vdiff[i] + ramo->regulador.deltaV)*(reg*ramo->regulador.Vnom/pow(3,0.5)/ramo->regulador.TP));
                }
                else{
                    if (Vdiff[i] >= ramo->regulador.deltaVr) 
                        ramo->regulador.tap[i] = ramo->regulador.tap[i] - ceil((Vdiff[i] - ramo->regulador.deltaVr)*(reg*ramo->regulador.Vnom/pow(3,0.5)/ramo->regulador.TP));
                    else if (Vdiff[i] <= -ramo->regulador.deltaVr) 
                        ramo->regulador.tap[i] = ramo->regulador.tap[i] - floor((Vdiff[i] + ramo->regulador.deltaVr)*(reg*ramo->regulador.Vnom/pow(3,0.5)/ramo->regulador.TP));
                
                }
            }
            break;
        //------------------------------------------------------
        // CONTROLE IDLE
        case 4 :
            for (i=0;i<3;i++) Vdiff[i] = Vset[i] - cabs(Vs[i]/ramo->regulador.TP - Zcont[i] * Isp[i]/ramo->regulador.TC);

            for (i=0;i<3;i++){
                if (Direction[i] ==1){
                    if (Vdiff[i] >= ramo->regulador.deltaV) 
                        ramo->regulador.tap[i] = ramo->regulador.tap[i] + ceil((Vdiff[i] - ramo->regulador.deltaV)*(reg*ramo->regulador.Vnom/pow(3,0.5)/ramo->regulador.TP));
                    else if (Vdiff[i] <= -ramo->regulador.deltaV) 
                        ramo->regulador.tap[i] = ramo->regulador.tap[i] + floor((Vdiff[i] + ramo->regulador.deltaV)*(reg*ramo->regulador.Vnom/pow(3,0.5)/ramo->regulador.TP));
                }
                else{
                    ramo->regulador.tap[i] = ramo->regulador.tap_ini[i];
                }
            }
            break;
        //------------------------------------------------------
        // CONTROLE NEUTRAL REVERSE  
        case 5:
            for (i=0;i<3;i++) Vdiff[i] = Vset[i] - cabs(Vs[i]/ramo->regulador.TP - Zcont[i] * Isp[i]/ramo->regulador.TC);

            for (i=0;i<3;i++){
                if (Direction[i] ==1){
                    if (Vdiff[i] >= ramo->regulador.deltaV) 
                        ramo->regulador.tap[i] = ramo->regulador.tap[i] + ceil((Vdiff[i] - ramo->regulador.deltaV)*(reg*ramo->regulador.Vnom/pow(3,0.5)/ramo->regulador.TP));
                    else if (Vdiff[i] <= -ramo->regulador.deltaV) 
                        ramo->regulador.tap[i] = ramo->regulador.tap[i] + floor((Vdiff[i] + ramo->regulador.deltaV)*(reg*ramo->regulador.Vnom/pow(3,0.5)/ramo->regulador.TP));
                }
                else{
                    ramo->regulador.tap[i] = 0;
                }
            }
            break;
        //------------------------------------------------------
        // CONTROLE COGENERATION 
        case 6:
            for (i=0;i<3;i++){
                if (Direction[i] == 1)
                    Vdiff[i] = Vset[i] - cabs(Vs[i]/ramo->regulador.TP - Zcont[i] * Isp[i]/ramo->regulador.TC);
                else
                    Vdiff[i] = Vset[i] - cabs(Vs[i]/ramo->regulador.TP - Zcont_r[i] * Isp[i]/ramo->regulador.TC);                    
            }

            for (i=0;i<3;i++){                
                if (Vdiff[i] >= ramo->regulador.deltaV) 
                    ramo->regulador.tap[i] = ramo->regulador.tap[i] + ceil((Vdiff[i] - ramo->regulador.deltaV)*(reg*ramo->regulador.Vnom/pow(3,0.5)/ramo->regulador.TP));
                else if (Vdiff[i] <= -ramo->regulador.deltaV) 
                    ramo->regulador.tap[i] = ramo->regulador.tap[i] + floor((Vdiff[i] + ramo->regulador.deltaV)*(reg*ramo->regulador.Vnom/pow(3,0.5)/ramo->regulador.TP));
            }
            break;
    }
    // Verifica limites de corrente
    
    
    // Verifica limites de tap
    for (i=0;i<3;i++){
        if ((ramo->regulador.tap[i] > ramo->regulador.ntaps) && (Direction[i] == 1))
            ramo->regulador.tap[i] = ramo->regulador.ntaps;
        else if ((ramo->regulador.tap[i] < -ramo->regulador.ntaps) && (Direction[i] == 1))
            ramo->regulador.tap[i] = -ramo->regulador.ntaps;
        else if ((ramo->regulador.tap[i] > ramo->regulador.ntaps) && (Direction[i] == -1))
            ramo->regulador.tap[i] = -ramo->regulador.ntaps;
        else if ((ramo->regulador.tap[i] < ramo->regulador.ntaps) && (Direction[i] == -1))
            ramo->regulador.tap[i] = ramo->regulador.ntaps;
    }

    // Indica mudança de tap
    ctr_act = 0;
    for (i=0;i<3;i++){
        if (ramo->regulador.tap[i] != taps_0[i]) ctr_act = 1;
    }
    
    // Atualiza akm
    atualizaTapRegulador(ramo);

    return(ctr_act);
}

//Função que efetua a etapa forward_sweep
BOOL forward_sweep(GRAFO *noP, GRAFO *grafo){
    int i, noAdj, idx;
    complex double Vaux[3];
    BOOL control_action = 0; //variável para indicar se houve transição por parte dos controladores
    double Sbase;
    double Ibase;

    BOOL control_REG_OPT;
    BOOL control_CAP_OPT;
    
    for (i = 0; i < noP->numeroAdjacentes;i ++){
        noAdj = noP->adjacentes[i].idNo;
        Vaux[0] = grafo[noAdj].V[0];
        Vaux[1] = grafo[noAdj].V[1];
        Vaux[2] = grafo[noAdj].V[2];
        
        if (noP->profundidade < grafo[noAdj].profundidade){
            
            calculaTensaoJusante(noP->V, Vaux, noP->adjacentes[i].Cur, noP->adjacentes[i].ramo);
            
            //printf("Vaux: %f + i*%f\n", creal(Vaux[0]), cimag(Vaux[0]));

            grafo[noAdj].V[0] = Vaux[0];
            grafo[noAdj].V[1] = Vaux[1];
            grafo[noAdj].V[2] = Vaux[2]; 

            //------------------------------------------------------------------
            //Atualiza TAPs de acordo com o controlador LDC do regulador
            if (control_REG_OPT == 1){
                if(noP->adjacentes[i].tipo == 2){
                    //Atualiza controle de taps do regulador de tensão
                    for (int j = 0; j < grafo[noAdj].numeroAdjacentes;j ++){
                        if (grafo[noAdj].adjacentes[j].idNo == noP->idNo){
                            idx = j;
                        }
                    }
                    control_action = controleReguladorTensao_LDC(noP->Vbase, Sbase/noP->Vbase, noP->V, Vaux, noP->adjacentes[i].Cur, grafo[noAdj].adjacentes[idx].Cur, noP->adjacentes[i].ramo);
                    // printf("\nLDC %d\n",control_action);
                }
            }
            //------------------------------------------------------------------
            //Atualiza Bancos de Capacitor Chaveados com os controladores de Bancos de Capacitor
            if (control_CAP_OPT == 1){

            }
            // //------------------------------------------------------------------
            // //Atualiza Geradores Distribuídos com os controladores de Tensão (Futuro)
            // if (control_GD_OPT == 1){

            // }

        }
    }
    return (control_action);
}

int fluxoPotencia_BFS_Alimentador(GRAFO *grafo, long int numeroBarras, ALIMENTADOR alimentador, DRAM *ramos, double Sbase){
    int it, nvar, k = 0,i, conv = 0;
    int MAXIT = 30;
    double tolerance = 0.00001;
    double *DV, nDV;
    long int *RNP;
    BOOL control_action = 0;

    FILABARRAS *barraAtual, *barraProx;
    NOADJACENTE *ramoAdj;
    
    //----------------------------------------------------------------------
    //RNP no formato de vetor
    //
    RNP = aloca_vetor_int(alimentador.numeroNos+1);
    barraAtual = &alimentador.rnp[0];
    while(barraAtual != NULL){
        RNP[k] = barraAtual->idNo;
        k++;
        barraAtual = barraAtual->prox;
    }
    nvar = alimentador.numeroNos*6;
    DV = aloca_vetor(nvar);
        
    ////----------------------------------------------------------------------
    // Fluxo de  Potência por Varredura Direta/Inversa via Soma de Correntes
    for (it = 0;it < MAXIT; it ++){
        //----------------------------------------------------------------------
        //Atualiza Correntes Nodais
        //
        barraAtual = &alimentador.rnp[0];
        k=0;
        while(barraAtual != NULL){
            atualizaInjecoes(&grafo[barraAtual->idNo]);
            for(i=0;i<3;i++){
                DV[k] = cabs(grafo[barraAtual->idNo].V[i]);
                k++;
                DV[k] = carg(grafo[barraAtual->idNo].V[i]);
                k++;                      
            }
            barraAtual = barraAtual->prox;
        }

        //---------------------------------------------------------------------- 
        //Backward Sweep
        for(k = alimentador.numeroNos-1; k >= 0; k--){
            backward(&grafo[RNP[k]], grafo);
        }

        //---------------------------------------------------------------------- 
        //forward_sweep Sweep
        for(k = 0; k < alimentador.numeroNos; k++){
            control_action = forward_sweep(&grafo[RNP[k]], grafo);
        }

        //Critério de Convergência
        barraAtual = &alimentador.rnp[0];
        k=0;
        while(barraAtual != NULL){
            for(i=0;i<3;i++){
                DV[k] = DV[k] - cabs(grafo[barraAtual->idNo].V[i]);
                k++;
                DV[k] = DV[k] - carg(grafo[barraAtual->idNo].V[i]);
                k++;                        
            }
            barraAtual = barraAtual->prox;
        }
        
        nDV = norma_inf(DV,nvar);
        printf("\nDv = %lf",nDV);
        if ((fabs(nDV) < tolerance) && !control_action){
            conv = 1;
            break;
        }
    }
    free(DV); free(RNP);
    return(it);
}

//Método Baseado na Varredura Direta/Inversa Trifasica
void fluxoPotencia_BFS_Multiplos(GRAFO *grafo, long int numeroBarras, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos,double Sbase){
    long int nmed,nvar,nmedTotal;
    int MAXIT = 30;
    int i,j, idAlim, it;
    FILABARRAS *barraAtual;
    GRAFO *no;
    
    // printf("Calculo de Fluxo de Potencia via Varredura Direta/Inversa...\n");
    incializaTensoesVarredura(grafo, numeroBarras, alimentadores, numeroAlimentadores);
//    incializaTensoesRaiz(grafo, numeroBarras, alimentadores, numeroAlimentadores);
    
    //Calculo de fluxo de potência para todos os alimentadores individualmente
    double tol = 0.000001;
    // clock_t tIni = clock();
    for (idAlim = 0; idAlim < numeroAlimentadores; idAlim++){
        it = fluxoPotencia_BFS_Alimentador(grafo, numeroBarras, alimentadores[idAlim], ramos, Sbase);
        
    }//Fim do loop dos alimentadores  
    // clock_t t1 = clock();
    // double tempoBFS = (double)(t1-tIni)/CLOCKS_PER_SEC;
    // printf("\nTempo computacional BFS: %lf\n\n",tempoBFS);   

    
    
    //Impressão de resultados
    for (idAlim = 0; idAlim < numeroAlimentadores; idAlim++){
        printf("\n\n");
        printf("Alimentador: %d\n",idAlim +1);
        printf("...\n");
        
        if (it < MAXIT) printf("\nConvergencia em %d iteracoes \n",it);
        else printf("\nNumero maximo de iteracoes atingido %d\n", MAXIT);
        
        //Imprime as tensões nodais
       printf("\nTensoes Nodais:\n");
       barraAtual = &alimentadores[idAlim].rnp[0];
       while(barraAtual != NULL){
           imprimeTensoesNodais(&grafo[barraAtual->idNo]);
           barraAtual = barraAtual->prox;
       }
        
    //    printf("\nInjecoes de Corrente\n");
    //    barraAtual = &alimentadores[idAlim].rnp[0];
    //    while(barraAtual != NULL){
    //        imprimeiInjecoesCorrentes(&grafo[barraAtual->idNo]);
    //        barraAtual = barraAtual->prox;
    //    }
        printf("\nTaps de Reguladores:\n");
        barraAtual = &alimentadores[idAlim].rnp[0];
        while(barraAtual != NULL){
            imprimeTaps(&grafo[barraAtual->idNo]);
            barraAtual = barraAtual->prox;
        }
    }

    
        
        
    free(barraAtual);
   
}







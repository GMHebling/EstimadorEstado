#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <complex.h>

#include "data_structures.h"
#include "funcoesTopologia.h"
#include "funcoesMatematicas.h"

void adicionaNo(FILABARRAS **setor, long int idNo)
{
    FILABARRAS *novoSetor;
    FILABARRAS *aux = NULL;
    
    novoSetor = (FILABARRAS *) malloc(sizeof(FILABARRAS));
    
    if(novoSetor == NULL) exit(EXIT_FAILURE);
    
    novoSetor->idNo = idNo;
    //novoSetor->prox = *setor;
    novoSetor->prox = NULL;
    aux = *setor;
    aux->prox = novoSetor;
    
    *setor = novoSetor;
}

void adicionaNoNaFila(FILABARRAS ** fila, long int idNo) {
    FILABARRAS *novoVertice = NULL;
    FILABARRAS *aux = NULL;
 
    novoVertice = (FILABARRAS *)malloc(sizeof(FILABARRAS));
 
    if(novoVertice == NULL) {
        printf("erro insere_fila\n");
        exit(EXIT_FAILURE);
    }
    
    novoVertice->idNo = idNo;
    novoVertice->prox = NULL;
    
    if(*fila == NULL)
        *fila = novoVertice;
    else {
        aux = *fila;
        while(aux->prox !=NULL) aux = aux->prox;
        aux->prox = novoVertice;
    }
}
void apontaProxNoNaFila(FILABARRAS ** fila) {
    FILABARRAS *novoVertice = NULL;
    FILABARRAS *aux = NULL;
    
    aux = *fila;
    while(aux->prox !=NULL) aux = aux->prox;
    *fila = aux;
    
}

int retiraNoDaFila(FILABARRAS ** fila) {
    FILABARRAS *aux = NULL;
    
    int idNo = -1;
    
    aux = *fila;
    
    if(aux != NULL) *fila = aux->prox;
    
    idNo = aux->idNo;
    
    free(aux);
 
    return idNo;
}

BOOL filaNaoVazia(FILABARRAS * fila) {
     if(fila == NULL)
         return false;
     return true;
}

BOOL estaLista(FILABARRAS *setor, int idNo) {
    FILABARRAS *p;
    p = setor;
    while(p != NULL && p->idNo != idNo)
        p = p->prox;
    if(p != NULL && p->idNo == idNo)
        return true;
    else
        return false;

}

BOOL ramoLigado(NOADJACENTE adjacente) {
    if(adjacente.estado == fechado) return true;
    return false;
}

BOOL estaListaAdjacencias(GRAFO *grafo, long int idNoRaiz, long int idNoAdj)
{
    int contador;
    for(contador =0; contador < grafo[idNoRaiz].numeroAdjacentes; contador++)
    {
        if(grafo[idNoRaiz].adjacentes[contador].idNo == idNoAdj)
            return true;
    }
    return false;
}

void buscaProfundidade(FILABARRAS *barraAtual, long int idNo, int profundidade,  BOOL *visitado, GRAFO * grafo, long int idAlim)
{
    //Depth-Search Algorithm - busca no e a sua profundidade (gera RNP))
    long int barraAdj,i = 0;
    
    visitado[idNo] = true;
    barraAtual->profundidade = profundidade;
    GRAFO * no = &grafo[idNo];
    grafo[idNo].idAlim = idAlim;
    grafo[idNo].profundidade = profundidade;
    //printf("\nidNo: %d  -  %d", grafo[idNo].barra->ID, profundidade);
    profundidade++;
    for(i = 0; i < no->numeroAdjacentes; i++)
    {   
        barraAdj = no->adjacentes[i].idNo;
        if ((visitado[barraAdj]== false))
            {
                idNo= barraAdj;
                //adicionaNo(&barraAtual, idNo);
                adicionaNoNaFila(&barraAtual, idNo);
                apontaProxNoNaFila(&barraAtual);
                buscaProfundidade(barraAtual, idNo, profundidade, visitado, grafo, idAlim);
            } 
    }
}

void buscaLargura(GRAFO * grafo, ALIMENTADOR *alimentador, long int idAlim, long int idNoRaiz, BOOL * visitado) {
    long int i, barraAdj, idNo,profundidade = 0;    
    
    FILABARRAS *barraAtual = NULL;
    FILABARRAS *filaProf = NULL;
    FILABARRAS *filaBarras = NULL;
    
    //adicionaNo(&alimentador[idAlim].rnp, idNoRaiz);
    barraAtual = &alimentador[idAlim].rnp[0];
    adicionaNoNaFila(&filaBarras, idNoRaiz);
    adicionaNoNaFila(&filaProf, profundidade);
    
    while(filaNaoVazia(filaBarras))
    {
        idNo = retiraNoDaFila(&filaBarras);
        if (idNo != idNoRaiz){
            adicionaNo(&barraAtual, idNo);
        }
        profundidade = retiraNoDaFila(&filaProf);
        barraAtual->profundidade = profundidade;
        GRAFO * no = &grafo[idNo];
        grafo[idNo].idAlim = idAlim;
        for(i = 0; i < no->numeroAdjacentes; i++)
        {
            barraAdj = no->adjacentes[i].idNo;
            if(ramoLigado(no->adjacentes[i]))
            {
                if(!visitado[barraAdj])
                {
                    /*if(estaLista(barraAtual, barraAdj))
                        printf("CICLO barra %ld idAlim %ld\n", barraAdj, idAlim);
                    else
                    {*/


                    adicionaNoNaFila(&filaBarras, barraAdj);
                    adicionaNoNaFila(&filaProf, profundidade + 1);
                    alimentador[idAlim].numeroNos++;

                    //}
                }
            }
        }
        visitado[idNo] = true;        
    }   
}


void buscaProfundidadeAlimentadores(GRAFO *grafo, long int numeroBarras, ALIMENTADOR **alimentadores, long int numeroAlimentadores) {
    int i, idAlim = 0;
    FILABARRAS * lista_barras = NULL;
    FILABARRAS *barraAtual = NULL;
    
    BOOL visitado[numeroBarras];
    
    if (((*alimentadores)= (ALIMENTADOR *)malloc( numeroAlimentadores * sizeof(ALIMENTADOR)))==NULL)
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para alimentadores !!!!");
        exit(1); 
    }
    
    for(i=0; i<numeroBarras; i++){ 
        visitado[i] = false;
        if(grafo[i].tipo == 2){
            (*alimentadores)[idAlim].idAlim = idAlim;
            (*alimentadores)[idAlim].noRaiz = i;
            (*alimentadores)[idAlim].numeroNos = 1;
            (*alimentadores)[idAlim].idRaiz = grafo[i].barra->ID;
            (*alimentadores)[idAlim].rnp[0].idNo = i;
            (*alimentadores)[idAlim].rnp[0].profundidade = 0;
            (*alimentadores)[idAlim].rnp[0].prox = NULL;
            idAlim++;
        }
    }
    
    for(i=0; i<numeroAlimentadores; i++)
    {
        //buscaLargura(grafo, (*alimentadores), i, (*alimentadores)[i].noRaiz, visitado);
        barraAtual = &(*alimentadores)[i].rnp[0];
        buscaProfundidade(barraAtual,(*alimentadores)[i].noRaiz,0,visitado,grafo,i);
        //printf("\n alimentador %d Raiz: %d   Nos: %d",i,(*alimentadores)[i].idRaiz,(*alimentadores)[i].numeroNos);
    }
    
}

//------------------------------------------------------------------------------
//
// ROTINAS DE TRATAMENTO DOS DADOS ELÉTRICOS
//
//------------------------------------------------------------------------------


//Converte a matriz de impedância em admitância e monta quadripolo dos ramos do grafo
void montaQuadripoloLinha(DRAM *ramo, DLIN *linha){
    int aux = 1;
    __complex__ double y, **Zl,**B;
    
    //Aloca Matrizes de Quadripolos
    ramo->Ypp = c_matAloca(3);
    ramo->Yps = c_matAloca(3);
    ramo->Ysp = c_matAloca(3);
    ramo->Yss = c_matAloca(3);
    
    //Aloca Matrizes de Impedância e Admitância
    Zl = c_matAloca(3);
    B = c_matAloca(3);
    
    //Matriz Impedância da linha
    Zl[0][0] = linha->Zaa;
    Zl[0][1] = linha->Zab;
    Zl[0][2] = linha->Zac;
    Zl[1][0] = linha->Zab;
    Zl[1][1] = linha->Zbb;
    Zl[1][2] = linha->Zbc;
    Zl[2][0] = linha->Zac;
    Zl[2][1] = linha->Zbc;
    Zl[2][2] = linha->Zcc;
    
    //Matriz Susceptãncia Shunt da linha
    B[0][0] = I * linha->Baa/2;
    B[0][1] = I * linha->Bab/2;
    B[0][2] = I * linha->Bac/2;
    B[1][0] = I * linha->Bab/2;
    B[1][1] = I * linha->Bbb/2;
    B[1][2] = I * linha->Bbc/2;
    B[2][0] = I * linha->Bac/2;
    B[2][1] = I * linha->Bbc/2;
    B[2][2] = I * linha->Bcc/2;
    
    //Inversa de Z - salva na variável Zl
    c_matInversaZ(Zl, 3);
    
    c_matIgual(ramo->Ypp, Zl, 3);
    c_matIgual(ramo->Yss, Zl, 3);
    c_matMultEsc(Zl, -1, 3);
    c_matIgual(ramo->Yps, Zl, 3);
    c_matIgual(ramo->Ysp, Zl, 3);
    
    
    ramo->Ypp[0][0] = ramo->Ypp[0][0] + B[0][0];
    ramo->Ypp[0][1] = ramo->Ypp[0][1] + B[0][1];
    ramo->Ypp[0][2] = ramo->Ypp[0][2] + B[0][2];
    ramo->Ypp[1][0] = ramo->Ypp[1][0] + B[1][0];
    ramo->Ypp[1][1] = ramo->Ypp[1][1] + B[1][1];
    ramo->Ypp[1][2] = ramo->Ypp[1][2] + B[1][2];
    ramo->Ypp[2][0] = ramo->Ypp[2][0] + B[2][0];
    ramo->Ypp[2][1] = ramo->Ypp[2][1] + B[2][1];
    ramo->Ypp[2][2] = ramo->Ypp[2][2] + B[2][2];
    
    ramo->Yss[0][0] = ramo->Yss[0][0] + B[0][0];
    ramo->Yss[0][1] = ramo->Yss[0][1] + B[0][1];
    ramo->Yss[0][2] = ramo->Yss[0][2] + B[0][2];
    ramo->Yss[1][0] = ramo->Yss[1][0] + B[1][0];
    ramo->Yss[1][1] = ramo->Yss[1][1] + B[1][1];
    ramo->Yss[1][2] = ramo->Yss[1][2] + B[1][2];
    ramo->Yss[2][0] = ramo->Yss[2][0] + B[2][0];
    ramo->Yss[2][1] = ramo->Yss[2][1] + B[2][1];
    ramo->Yss[2][2] = ramo->Yss[2][2] + B[2][2];
    
//    printf("\n\Quadripolo barra para %d\n\n", noadj->idNo);
//    c_matImprime(ramo->Ypp,3);
//    printf("\n");
//    c_matImprime(ramo->Yps,3);
//    printf("\n");
//    c_matImprime(ramo->Ysp,3);
//    printf("\n");
//    c_matImprime(ramo->Yss,3);
//    printf("\n");
    
}

//Converte a matriz de impedância em admitância e monta quadripolo dos ramos do grafo
void montaQuadripoloTrafo(DRAM *ramo, DTRF *trafo){
    int aux = 1;
    __complex__ double y, **Yi, **Yii,**Yiii;
    
    //Aloca Matrizes de Quadripolos
    ramo->Ypp = c_matAloca(3);
    ramo->Yps = c_matAloca(3);
    ramo->Ysp = c_matAloca(3);
    ramo->Yss = c_matAloca(3);
    
    //Aloca Matrizes de Ligação D-Y-YN
    Yi = c_matAloca(3);
    Yii = c_matAloca(3);
    Yiii = c_matAloca(3);
    
    //Admitância do trafo
    y = 1/(trafo->R + I*trafo->X);
    
    //Yi
    Yi[0][0] = 1;
    Yi[0][1] = 0;
    Yi[0][2] = 0;
    Yi[1][0] = 0;
    Yi[1][1] = 1;
    Yi[1][2] = 0;
    Yi[2][0] = 0;
    Yi[2][1] = 0;
    Yi[2][2] = 1;
    //Yii
    Yii[0][0] = 2;
    Yii[0][1] = -1;
    Yii[0][2] = -1;
    Yii[1][0] = -1;
    Yii[1][1] = 2;
    Yii[1][2] = -1;
    Yii[2][0] = -1;
    Yii[2][1] = -1;
    Yii[2][2] = 2;
    //Yiii
    Yiii[0][0] = -1/(pow(3,0.5));
    Yiii[0][1] = 1/(pow(3,0.5));
    Yiii[0][2] = 0;
    Yiii[1][0] = 0;
    Yiii[1][1] = -1/(pow(3,0.5));
    Yiii[1][2] = 1/(pow(3,0.5));
    Yiii[2][0] = 1/(pow(3,0.5));
    Yiii[2][1] = 0;
    Yiii[2][2] = -1/(pow(3,0.5));
    
    c_matMultEsc(Yi, y, 3);
    c_matMultEsc(Yii, y/3, 3);
    c_matMultEsc(Yiii, y, 3);
        
    if((trafo->lig_pri == 1)&& (trafo->lig_sec == 1)){ //Ligação YN-YN
        c_matIgual(ramo->Ypp, Yi, 3);
        c_matIgual(ramo->Yps, Yi, 3);
        c_matMultEsc(ramo->Yps, -1, 3);
        c_matIgual(ramo->Ysp, Yi, 3);
        c_matMultEsc(ramo->Ysp, -1, 3);
        c_matIgual(ramo->Yss, Yi, 3);
    }
    else if((trafo->lig_pri == 1)&& (trafo->lig_sec == 3)){ //Ligação YN-Y
        c_matIgual(ramo->Ypp, Yii, 3);
        c_matIgual(ramo->Yps, Yii, 3);
        c_matMultEsc(ramo->Yps, -1, 3);
        c_matIgual(ramo->Ysp, Yii, 3);
        c_matMultEsc(ramo->Ysp, -1, 3);
        c_matIgual(ramo->Yss, Yii, 3);
        ramo->Ypp[2][0] = ramo->Ypp[2][0] + 100*cabs(y); //Sequência zero para manter o quadripólo com posto completo
        ramo->Ypp[2][1] = ramo->Ypp[2][1] + 100*cabs(y);
        ramo->Ypp[2][2] = ramo->Ypp[2][2] + 100*cabs(y);
    }
    else if((trafo->lig_pri == 1)&& (trafo->lig_sec == 2)){ //Ligação YN-D
        c_matIgual(ramo->Ypp, Yi, 3);
        c_matIgual(ramo->Yps, Yiii, 3);
        c_matIgual(ramo->Ysp, Yiii, 3);
        c_matTransp(ramo->Ysp, 3);
        c_matIgual(ramo->Yss, Yii, 3);
        ramo->Yss[2][0] = ramo->Yss[2][0] + 100*cabs(y); //Sequência zero para manter o quadripólo com posto completo
        ramo->Yss[2][1] = ramo->Yss[2][1] + 100*cabs(y);
        ramo->Yss[2][2] = ramo->Yss[2][2] + 100*cabs(y);
    }
    else if((trafo->lig_pri == 3)&& (trafo->lig_sec == 1)){ //Ligação Y-YN
        c_matIgual(ramo->Ypp, Yii, 3);
        c_matIgual(ramo->Yps, Yii, 3);
        c_matMultEsc(ramo->Yps, -1, 3);
        c_matIgual(ramo->Ysp, Yii, 3);
        c_matMultEsc(ramo->Ysp, -1, 3);
        c_matIgual(ramo->Yss, Yii, 3);
        ramo->Yss[2][0] = ramo->Yss[2][0] + 100*cabs(y); //Sequência zero para manter o quadripólo com posto completo
        ramo->Yss[2][1] = ramo->Yss[2][1] + 100*cabs(y);
        ramo->Yss[2][2] = ramo->Yss[2][2] + 100*cabs(y);
    }
    else if((trafo->lig_pri == 3)&& (trafo->lig_sec == 3)){ //Ligação Y-Y
        c_matIgual(ramo->Ypp, Yii, 3);
        c_matIgual(ramo->Yps, Yii, 3);
        c_matMultEsc(ramo->Yps, -1, 3);
        c_matIgual(ramo->Ysp, Yii, 3);
        c_matMultEsc(ramo->Ysp, -1, 3);
        c_matIgual(ramo->Yss, Yii, 3);
        ramo->Yss[2][0] = ramo->Yss[2][0] + 100*cabs(y); //Sequência zero para manter o quadripólo com posto completo
        ramo->Yss[2][1] = ramo->Yss[2][1] + 100*cabs(y); //Pensar se trata-se de simplificação, pois pode existir tensão de neutro neste caso, mas não irá te corrente de neutro.
        ramo->Yss[2][2] = ramo->Yss[2][2] + 100*cabs(y);
    }
    else if((trafo->lig_pri == 3)&& (trafo->lig_sec == 2)){ //Ligação Y-D
        c_matIgual(ramo->Ypp, Yii, 3);
        c_matIgual(ramo->Yps, Yiii, 3);
        c_matIgual(ramo->Ysp, Yiii, 3);
        c_matTransp(ramo->Ysp, 3);
        c_matIgual(ramo->Yss, Yii, 3);
        ramo->Yss[2][0] = ramo->Yss[2][0] + 100*cabs(y); //Sequência zero para manter o quadripólo com posto completo
        ramo->Yss[2][1] = ramo->Yss[2][1] + 100*cabs(y);
        ramo->Yss[2][2] = ramo->Yss[2][2] + 100*cabs(y);
    }
    else if((trafo->lig_pri == 2)&& (trafo->lig_sec == 1)){ //Ligação D-YN
        c_matIgual(ramo->Ypp, Yii, 3);
        c_matIgual(ramo->Yps, Yiii, 3);
        c_matIgual(ramo->Ysp, Yiii, 3);
        c_matTransp(ramo->Ysp, 3);
        c_matIgual(ramo->Yss, Yi, 3);
        ramo->Ypp[2][0] = ramo->Ypp[2][0] + 100*cabs(y); //Sequência zero para manter o quadripólo com posto completo
        ramo->Ypp[2][1] = ramo->Ypp[2][1] + 100*cabs(y);
        ramo->Ypp[2][2] = ramo->Ypp[2][2] + 100*cabs(y);
        ramo->Yps[2][0] = ramo->Yps[2][0] + 100*cabs(y); //Sequência zero para manter o quadripólo com posto completo
        ramo->Yps[2][1] = ramo->Yps[2][1] + 100*cabs(y);
        ramo->Yps[2][2] = ramo->Yps[2][2] + 100*cabs(y);
    }
    else if((trafo->lig_pri == 2)&& (trafo->lig_sec == 3)){ //Ligação D-Y
        c_matIgual(ramo->Ypp, Yii, 3);
        c_matIgual(ramo->Yps, Yiii, 3);
        c_matTransp(ramo->Yps, 3);
        c_matIgual(ramo->Ysp, Yiii, 3);
        c_matIgual(ramo->Yss, Yii, 3);
        ramo->Ypp[2][0] = ramo->Ypp[2][0] + 100*cabs(y); //Sequência zero para manter o quadripólo com posto completo
        ramo->Ypp[2][1] = ramo->Ypp[2][1] + 100*cabs(y);
        ramo->Ypp[2][2] = ramo->Ypp[2][2] + 100*cabs(y);
        
        // Aproximação para o caso do Y não aterrado - tensão de sequência é igual a zero no secundário - fica igual DYn
        ramo->Yps[2][0] = ramo->Yps[2][0] + 100*cabs(y); //Sequência zero para manter o quadripólo com posto completo
        ramo->Yps[2][1] = ramo->Yps[2][1] + 100*cabs(y);
        ramo->Yps[2][2] = ramo->Yps[2][2] + 100*cabs(y);
    }
    else if((trafo->lig_pri == 2)&& (trafo->lig_sec == 2)){ //Ligação D-D
        c_matIgual(ramo->Ypp, Yii, 3);
        c_matIgual(ramo->Yps, Yii, 3);
        c_matMultEsc(ramo->Yps, -1, 3);
        c_matIgual(ramo->Ysp, Yii, 3);
        c_matMultEsc(ramo->Ysp, -1, 3);
        c_matIgual(ramo->Yss, Yii, 3);
//        ramo->Ypp[2][0] = ramo->Ypp[2][0] + 1; //Sequência zero para manter o quadripólo com posto completo
//        ramo->Ypp[2][1] = ramo->Ypp[2][1] + 1;
//        ramo->Ypp[2][2] = ramo->Ypp[2][2] + 1;
//        ramo->Yps[2][0] = ramo->Yps[2][0] + 1; //Sequência zero para manter o quadripólo com posto completo
//        ramo->Yps[2][1] = ramo->Yps[2][1] + 1;
//        ramo->Yps[2][2] = ramo->Yps[2][2] + 1;
        ramo->Yss[2][0] = ramo->Yss[2][0] + 100*cabs(y); //Sequência zero para manter o quadripólo com posto completo
        ramo->Yss[2][1] = ramo->Yss[2][1] + 100*cabs(y);
        ramo->Yss[2][2] = ramo->Yss[2][2] + 100*cabs(y);
//        ramo->Ysp[2][0] = ramo->Ysp[2][0] + 1; //Sequência zero para manter o quadripólo com posto completo
//        ramo->Ysp[2][1] = ramo->Ysp[2][1] + 1;
//        ramo->Ysp[2][2] = ramo->Ysp[2][2] + 1;
        
    }
//    printf("\n\Quadripolo\n\n");
//    c_matImprime(ramo->Ypp,3);
//    printf("\n");
//    c_matImprime(ramo->Yps,3);
//    printf("\n");
//    c_matImprime(ramo->Ysp,3);
//    printf("\n");
//    c_matImprime(ramo->Yss,3);
//    printf("\n");
    
}

//Converte a matriz de impedância em admitância e monta quadripolo dos ramos do grafo
void montaQuadripoloRegulador(DRAM *ramo, DREG *reg){
    int aux = 1;
    __complex__ double y, **Yi, **Yii,**Yiii;
    
    //Aloca Matrizes de Quadripolos
    ramo->Ypp = c_matAloca(3);
    ramo->Yps = c_matAloca(3);
    ramo->Ysp = c_matAloca(3);
    ramo->Yss = c_matAloca(3);
    
    //Aloca Matrizes de Ligação D-Y-YN
    Yi = c_matAloca(3);
    Yii = c_matAloca(3);
    Yiii = c_matAloca(3);
    
    //Admitância do trafo
    y = 1/(reg->R + I*reg->X);
    
    //Yi
    Yi[0][0] = 1;
    Yi[0][1] = 0;
    Yi[0][2] = 0;
    Yi[1][0] = 0;
    Yi[1][1] = 1;
    Yi[1][2] = 0;
    Yi[2][0] = 0;
    Yi[2][1] = 0;
    Yi[2][2] = 1;
    //Yii
    Yii[0][0] = 2;
    Yii[0][1] = -1;
    Yii[0][2] = -1;
    Yii[1][0] = -1;
    Yii[1][1] = 2;
    Yii[1][2] = -1;
    Yii[2][0] = -1;
    Yii[2][1] = -1;
    Yii[2][2] = 2;
    //Yiii
    Yiii[0][0] = -1/(pow(3,0.5));
    Yiii[0][1] = 1/(pow(3,0.5));
    Yiii[0][2] = 0;
    Yiii[1][0] = 0;
    Yiii[1][1] = -1/(pow(3,0.5));
    Yiii[1][2] = 1/(pow(3,0.5));
    Yiii[2][0] = 1/(pow(3,0.5));
    Yiii[2][1] = 0;
    Yiii[2][2] = -1/(pow(3,0.5));
    
    
    //Feito somente para Yn
    switch (ramo->fases){
        case 1:
            Yi[1][1] = 0;
            Yi[2][2] = 0;
            break;    
        case 2:
            Yi[0][0] = 0;
            Yi[2][2] = 0;
            break;
        case 3:
            Yi[1][1] = 0;
            Yi[0][0] = 0;
            break;
        case 4:
            Yi[2][2] = 0;
            break;    
        case 5:
            Yi[1][1] = 0;
            break;
        case 6:
            Yi[0][0] = 0;
            break; 
    }   
       
    c_matMultEsc(Yi, y, 3);
    c_matMultEsc(Yii, y/3, 3);
    c_matMultEsc(Yiii, y, 3);
    
    //Feito Somente para YN!!!!!!!!!!!!!!!!!!!!!!!!!
    if(reg->lig == 1){ //Ligação YN
        c_matIgual(ramo->Ypp, Yi, 3);
        c_matIgual(ramo->Yps, Yi, 3);
        c_matMultEsc(ramo->Yps, -1, 3);
        c_matIgual(ramo->Ysp, Yi, 3);
        c_matMultEsc(ramo->Ysp, -1, 3);
        c_matIgual(ramo->Yss, Yi, 3);
    }
    else if(reg->lig == 2){ //Ligação D
        c_matIgual(ramo->Ypp, Yii, 3);
        c_matIgual(ramo->Yps, Yii, 3);
        c_matMultEsc(ramo->Yps, -1, 3);
        c_matIgual(ramo->Ysp, Yii, 3);
        c_matMultEsc(ramo->Ysp, -1, 3);
        c_matIgual(ramo->Yss, Yii, 3);
    }
    else if(reg->lig == 3){ //Ligação Y
        c_matIgual(ramo->Ypp, Yi, 3);
        c_matIgual(ramo->Yps, Yiii, 3);
        c_matIgual(ramo->Ysp, Yiii, 3);
        c_matTransp(ramo->Ysp, 3);
        c_matIgual(ramo->Yss, Yii, 3);
    }
    else if(reg->lig == 4){ //Ligação OYn
        c_matIgual(ramo->Ypp, Yii, 3);
        c_matIgual(ramo->Yps, Yii, 3);
        c_matMultEsc(ramo->Yps, -1, 3);
        c_matIgual(ramo->Ysp, Yii, 3);
        c_matMultEsc(ramo->Ysp, -1, 3);
        c_matIgual(ramo->Yss, Yii, 3);
    }
    else if(reg->lig == 5){ //Ligação OD
        c_matIgual(ramo->Ypp, Yii, 3);
        c_matIgual(ramo->Yps, Yii, 3);
        c_matMultEsc(ramo->Yps, -1, 3);
        c_matIgual(ramo->Ysp, Yii, 3);
        c_matMultEsc(ramo->Ysp, -1, 3);
        c_matIgual(ramo->Yss, Yii, 3);
    }
    
    
//    printf("\n\Quadripolo\n\n");
//    c_matImprime(ramo->Ypp,3);
//    printf("\n");
//    c_matImprime(ramo->Yps,3);
//    printf("\n");
//    c_matImprime(ramo->Ysp,3);
//    printf("\n");
//    c_matImprime(ramo->Yss,3);
//    printf("\n");
}

//Monta quadripolo dos shunts do grafo
void montaQuadripoloShunt(GRAFO *no, DSHNT *shunt){
    int aux = 1, i, j;
    __complex__ double ya,yb,yc, **Yi, **Yii,**Yiii,**Ysh;
    
    //Aloca Matrizes de Quadripolos
    Ysh = c_matAloca(3);
    
    
    //Aloca Matrizes de Ligação D-Y-YN
    Yi = c_matAloca(3);
    Yii = c_matAloca(3);
    Yiii = c_matAloca(3);
    
    //Admitância do trafo
    
    ya = I*shunt->Qnom[0]/(pow(shunt->Vbase/sqrt(3),2)/pow(no->Vbase,2));
    yb = I*shunt->Qnom[1]/(pow(shunt->Vbase/sqrt(3),2)/pow(no->Vbase,2));
    yc = I*shunt->Qnom[2]/(pow(shunt->Vbase/sqrt(3),2)/pow(no->Vbase,2));
    
    //Yi
    Yi[0][0] = 1*ya;
    Yi[0][1] = 0;
    Yi[0][2] = 0;
    Yi[1][0] = 0;
    Yi[1][1] = 1*yb;
    Yi[1][2] = 0;
    Yi[2][0] = 0;
    Yi[2][1] = 0;
    Yi[2][2] = 1*yc;
    //Yii
    Yii[0][0] = 2*ya;
    Yii[0][1] = -1*yb;
    Yii[0][2] = -1*yc;
    Yii[1][0] = -1*ya;
    Yii[1][1] = 2*yb;
    Yii[1][2] = -1*yc;
    Yii[2][0] = -1*ya;
    Yii[2][1] = -1*yb;
    Yii[2][2] = 2*yc;
    //Yiii
    Yiii[0][0] = -1/(pow(3,0.5))*ya;
    Yiii[0][1] = 1/(pow(3,0.5))*yb;
    Yiii[0][2] = 0;
    Yiii[1][0] = 0;
    Yiii[1][1] = -1/(pow(3,0.5))*yb;
    Yiii[1][2] = 1/(pow(3,0.5))*yc;
    Yiii[2][0] = 1/(pow(3,0.5))*ya;
    Yiii[2][1] = 0;
    Yiii[2][2] = -1/(pow(3,0.5))*yc;
    
    if(shunt->lig == 1){ //Ligação YN
        c_matIgual(Ysh, Yi, 3);
    }
    else if(shunt->lig == 2){ //Ligação D trifásico
        c_matIgual(Ysh, Yiii, 3);
        Ysh[2][0] = Ysh[2][0] + 1; //Sequência zero para manter o quadripólo com posto completo
        Ysh[2][1] = Ysh[2][1] + 1;
        Ysh[2][2] = Ysh[2][2] + 1;
    }
    else if(shunt->lig == 3){ //Ligação Y
        c_matIgual(Ysh, Yii, 3);
    }
    for(i=0;i<3;i++){
        for(j=0;j<3;j++){
            no->Ysh[i][j] += Ysh[i][j];
        }
    }
    
//    printf("\n\Quadripolo\n\n");
//    c_matImprime(ramo->Ypp,3);
//    printf("\n");
//    c_matImprime(ramo->Yps,3);
//    printf("\n");
//    c_matImprime(ramo->Ysp,3);
//    printf("\n");
//    c_matImprime(ramo->Yss,3);
//    printf("\n");
    
}

//Função que transforma os dados do grafo em pu e monta as matrizes de admitância dos elementos do grafo - 
// sem busca em profundidade no grafo - Monta Elemento a Elemento de acordo com dados de entrada
void calculaPU(GRAFO *grafo, long int numeroBarras, DRAM *ramos, long int numeroRamos, double Sbase) {
    int i, idNo, idRam;
    FILABARRAS * lista_barras = NULL;
    long int barraAdj = 0;
    double Vbase;
    
    
    //Transforma em PU informações dos nós do grafo
    for (idNo=0;idNo<numeroBarras;idNo++){
        //Transforma em pu os dados de carga - Futuro para fluxo de carga
        /*for(i=0;i<grafo[idNo].barra->nloads;i++){

        }*/ 
        
        //Transforma em PU os shunts da rede elétrica
        for(i=0;i<grafo[idNo].barra->nshunts;i++){
            grafo[idNo].barra->shunts[i].Qnom[0] = 1000*grafo[idNo].barra->shunts[i].Qnom[0]/Sbase; 
            grafo[idNo].barra->shunts[i].Qnom[1] = 1000*grafo[idNo].barra->shunts[i].Qnom[1]/Sbase;
            grafo[idNo].barra->shunts[i].Qnom[2] = 1000*grafo[idNo].barra->shunts[i].Qnom[2]/Sbase;
            montaQuadripoloShunt(&grafo[idNo],&grafo[idNo].barra->shunts[i]);
        }
        
        /*//Transforma em pu os dados de GDs - Futuro para fluxo de carga
        for(i=0;i<grafo[idNo].barra->ngds;i++){

        }*/
    }
    
    //Transforma em PU informações dos ramos do grafo
    for (idRam=0;idRam<numeroRamos;idRam++){
        Vbase = grafo[ramos[idRam].m].Vbase;
        //Transforma as impedâncias em pu
        switch(ramos[idRam].tipo){
            case 0:
                ramos[idRam].linha.Zaa = ramos[idRam].linha.Zaa/((pow(Vbase,2))/Sbase);
                ramos[idRam].linha.Zab = ramos[idRam].linha.Zab/((pow(Vbase,2))/Sbase);
                ramos[idRam].linha.Zac = ramos[idRam].linha.Zac/((pow(Vbase,2))/Sbase);
                ramos[idRam].linha.Zbb = ramos[idRam].linha.Zbb/((pow(Vbase,2))/Sbase);
                ramos[idRam].linha.Zbc = ramos[idRam].linha.Zbc/((pow(Vbase,2))/Sbase);
                ramos[idRam].linha.Zcc = ramos[idRam].linha.Zcc/((pow(Vbase,2))/Sbase);
                
                ramos[idRam].linha.Baa = ramos[idRam].linha.Baa/((pow(Vbase,2))/Sbase);
                ramos[idRam].linha.Bab = ramos[idRam].linha.Bab/((pow(Vbase,2))/Sbase);
                ramos[idRam].linha.Bac = ramos[idRam].linha.Bac/((pow(Vbase,2))/Sbase);
                ramos[idRam].linha.Bbb = ramos[idRam].linha.Bbb/((pow(Vbase,2))/Sbase);
                ramos[idRam].linha.Bbc = ramos[idRam].linha.Bbc/((pow(Vbase,2))/Sbase);
                ramos[idRam].linha.Bcc = ramos[idRam].linha.Bcc/((pow(Vbase,2))/Sbase);
                
                
                montaQuadripoloLinha(&ramos[idRam], &ramos[idRam].linha);
                break;
            case 1:
                ramos[idRam].trafo.R = 3*ramos[idRam].trafo.R/((pow(Vbase,2))/Sbase);
                ramos[idRam].trafo.X = 3*ramos[idRam].trafo.X/((pow(Vbase,2))/Sbase);
                
                montaQuadripoloTrafo(&ramos[idRam], &ramos[idRam].trafo);
                break;
            case 2:
                ramos[idRam].regulador.R = 3*ramos[idRam].regulador.R/((pow(Vbase,2))/Sbase);
                ramos[idRam].regulador.X = 3*ramos[idRam].regulador.X/((pow(Vbase,2))/Sbase);
                
                montaQuadripoloRegulador(&ramos[idRam], &ramos[idRam].regulador);
                break;    
        }
    }   
}

//Função que transforma os dados do grafo em pu e monta as matrizes de admitância dos elementos do grafo
void buscaPU(long int idNo, double Vbase,  BOOL *visitado, GRAFO * grafo, double Sbase)
{
    //Depth-Search Algorithm - busca no e a sua profundidade (gera RNP))
    long int barraAdj,i = 0;
    
    visitado[idNo] = true;
    
    GRAFO * no = &grafo[idNo];
    grafo[idNo].Vbase = Vbase;
    
    //Transforma em pu os dados de carga - Futuro para fluxo de carga
    /*for(i=0;i<grafo[idNo].barra->nloads;i++){
                
    }*/    
    //Transforma em pu os dados de shunts
    for(i=0;i<grafo[idNo].barra->nshunts;i++){
        grafo[idNo].barra->shunts[i].Qnom[0] = 1000*grafo[idNo].barra->shunts[i].Qnom[0]/Sbase; 
        grafo[idNo].barra->shunts[i].Qnom[1] = 1000*grafo[idNo].barra->shunts[i].Qnom[1]/Sbase;
        grafo[idNo].barra->shunts[i].Qnom[2] = 1000*grafo[idNo].barra->shunts[i].Qnom[2]/Sbase;
        montaQuadripoloShunt(&grafo[idNo],&grafo[idNo].barra->shunts[i]);
    }    
    /*//Transforma em pu os dados de GDs - Futuro para fluxo de carga
    for(i=0;i<grafo[idNo].barra->ngds;i++){
                
    }*/
    
    
    //Tratamento dos Ramos - Recursivamente
    for(i = 0; i < no->numeroAdjacentes; i++)
    {   
        barraAdj = no->adjacentes[i].idNo;
        
        if ((visitado[barraAdj]== false)){
            idNo= barraAdj;
            //Atualiza a Base para a próxima barra caso seja um trafo               
            if(no->adjacentes[i].tipo == 1){
                Vbase = Vbase*no->adjacentes[i].relacao;
            }
            
            //Transforma as impedâncias em pu
            switch(no->adjacentes[i].tipo){
                case 0:
                    no->adjacentes[i].ramo->linha.Zaa = no->adjacentes[i].ramo->linha.Zaa/((pow(Vbase,2))/Sbase);
                    no->adjacentes[i].ramo->linha.Zab = no->adjacentes[i].ramo->linha.Zab/((pow(Vbase,2))/Sbase);
                    no->adjacentes[i].ramo->linha.Zac = no->adjacentes[i].ramo->linha.Zac/((pow(Vbase,2))/Sbase);
                    no->adjacentes[i].ramo->linha.Zbb = no->adjacentes[i].ramo->linha.Zbb/((pow(Vbase,2))/Sbase);
                    no->adjacentes[i].ramo->linha.Zbc = no->adjacentes[i].ramo->linha.Zbc/((pow(Vbase,2))/Sbase);
                    no->adjacentes[i].ramo->linha.Zcc = no->adjacentes[i].ramo->linha.Zcc/((pow(Vbase,2))/Sbase);
                    
                    no->adjacentes[i].ramo->linha.Baa = no->adjacentes[i].ramo->linha.Baa*((pow(Vbase,2))/Sbase);
                    no->adjacentes[i].ramo->linha.Bab = no->adjacentes[i].ramo->linha.Bab*((pow(Vbase,2))/Sbase);
                    no->adjacentes[i].ramo->linha.Bac = no->adjacentes[i].ramo->linha.Bac*((pow(Vbase,2))/Sbase);
                    no->adjacentes[i].ramo->linha.Bbb = no->adjacentes[i].ramo->linha.Bbb*((pow(Vbase,2))/Sbase);
                    no->adjacentes[i].ramo->linha.Bbc = no->adjacentes[i].ramo->linha.Bbc*((pow(Vbase,2))/Sbase);
                    no->adjacentes[i].ramo->linha.Bcc = no->adjacentes[i].ramo->linha.Bcc*((pow(Vbase,2))/Sbase);
                    
                    montaQuadripoloLinha(no->adjacentes[i].ramo, &no->adjacentes[i].ramo->linha);
                    break;
                case 1:
                    no->adjacentes[i].ramo->trafo.R = 3*no->adjacentes[i].ramo->trafo.R/((pow(Vbase,2))/Sbase);
                    no->adjacentes[i].ramo->trafo.X = 3*no->adjacentes[i].ramo->trafo.X/((pow(Vbase,2))/Sbase);
                    montaQuadripoloTrafo(no->adjacentes[i].ramo, &no->adjacentes[i].ramo->trafo);
                    break;
                case 2:
                    no->adjacentes[i].ramo->regulador.R = 3*no->adjacentes[i].ramo->regulador.R/((pow(Vbase,2))/Sbase);
                    no->adjacentes[i].ramo->regulador.X = 3*no->adjacentes[i].ramo->regulador.X/((pow(Vbase,2))/Sbase);
                    montaQuadripoloRegulador(no->adjacentes[i].ramo, &no->adjacentes[i].ramo->regulador);
                    break;    
            }
            
            buscaPU(idNo, Vbase, visitado, grafo, Sbase); //recursão
        }
    }
}

//Transforma os dados em pu e termina de montar os dados do grado (matrizes Y)
void tranformaPU(GRAFO *grafo, long int numeroBarras, double Sbase, ALIMENTADOR *alimentadores, long int numeroAlimentadores) {
    int i;
    FILABARRAS * lista_barras = NULL;
    
    BOOL visitado[numeroBarras];
    for(i=0; i<numeroBarras; i++){ 
        visitado[i] = false;
    }
    
    //Transforma dados lidos em pu
    for(i=0; i<numeroAlimentadores; i++)
    {
        buscaPU(alimentadores[i].noRaiz, grafo[alimentadores[i].noRaiz].Vbase, visitado, grafo, Sbase);
    }
    
    /*for(i=0; i<numeroAlimentadores; i++)
    {
        FILABARRAS *barraAtual = &alimentadores[i].rnp[0];
        while(barraAtual != NULL)
        {
            printf("\nBarra[%d]: %d  Prof: %d  Vbase = %.5lf",barraAtual->idNo,grafo[barraAtual->idNo].barra->ID,barraAtual->profundidade,grafo[barraAtual->idNo].Vbase);
            barraAtual = barraAtual->prox;
        }
    }*/
}

//Atualiza as matrizes Y de acordo com os taps dos transformadores
void atualizaTaps(DRAM *ramos, long int numeroRamos){
    int i;
    
    for(i=0;i<numeroRamos;i++){
        if (ramos[i].tipo == 2){
            ramos[i].tap_pri[0] = (1 + ramos[i].regulador.tap[0]*ramos[i].regulador.regulacao/ramos[i].regulador.ntaps);
            ramos[i].tap_pri[1] = (1 + ramos[i].regulador.tap[1]*ramos[i].regulador.regulacao/ramos[i].regulador.ntaps);
            ramos[i].tap_pri[2] = (1 + ramos[i].regulador.tap[2]*ramos[i].regulador.regulacao/ramos[i].regulador.ntaps);
            
            ramos[i].tap_sec[0] = 1;
            ramos[i].tap_sec[1] = 1;
            ramos[i].tap_sec[2] = 1;
            
            ramos[i].Ypp[0][0] = pow(ramos[i].tap_pri[0],2)*ramos[i].Ypp[0][0];
            ramos[i].Ypp[1][1] = pow(ramos[i].tap_pri[1],2)*ramos[i].Ypp[1][1];
            ramos[i].Ypp[2][2] = pow(ramos[i].tap_pri[2],2)*ramos[i].Ypp[2][2];
            
            ramos[i].Yps[0][0] = ramos[i].tap_sec[0]*ramos[i].tap_pri[0]*ramos[i].Yps[0][0];
            ramos[i].Yps[1][1] = ramos[i].tap_sec[1]*ramos[i].tap_pri[1]*ramos[i].Yps[1][1];
            ramos[i].Yps[2][2] = ramos[i].tap_sec[2]*ramos[i].tap_pri[2]*ramos[i].Yps[2][2];
            
            ramos[i].Ysp[0][0] = ramos[i].tap_sec[0]*ramos[i].tap_pri[0]*ramos[i].Ysp[0][0];
            ramos[i].Ysp[1][1] = ramos[i].tap_sec[1]*ramos[i].tap_pri[1]*ramos[i].Ysp[1][1];
            ramos[i].Ysp[2][2] = ramos[i].tap_sec[2]*ramos[i].tap_pri[2]*ramos[i].Ysp[2][2];
            
            ramos[i].Yss[0][0] = pow(ramos[i].tap_sec[0],2)*ramos[i].Yss[0][0];
            ramos[i].Yss[1][1] = pow(ramos[i].tap_sec[1],2)*ramos[i].Yss[1][1];
            ramos[i].Yss[2][2] = pow(ramos[i].tap_sec[2],2)*ramos[i].Yss[2][2];
        }
        else if (ramos[i].tipo == 1){ //REVISAR PARA TAP FORA DA NOMINAL EM TRAFOS
            ramos[i].tap_pri[0] = 1;
            ramos[i].tap_pri[1] = 1;
            ramos[i].tap_pri[2] = 1;
            
            ramos[i].tap_sec[0] = 1;
            ramos[i].tap_sec[1] = 1;
            ramos[i].tap_sec[2] = 1;
        }
                
    
    }
    
    
    
}

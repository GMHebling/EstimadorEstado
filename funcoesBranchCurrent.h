/* 
 * File:   funcoesBranchCurrent.h
 * Author: Gustavo Hebling
 *
 * Created on 15 de Abril de 2021, 18:20
 */

#ifndef FUNCOESBRANCHCURRENT_H
#define FUNCOESBRANCHCURRENT_H


void estimadorBC_RECT(GRAFO *grafo, long int numeroBarras, DMED *medidas, long int **numeroMedidas, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos,double Sbase);

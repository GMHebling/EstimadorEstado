/* 
 * File:   funcoesBranchCurrent.h
 * Author: Gustavo Hebling
 *
 * Created on 15 de Abril de 2021, 18:20
 */

#ifndef funcoesBranchCurrent_H
#define funcoesBranchCurrent_H

void estimadorBC_RECT(GRAFO *grafo, long int numeroRamos, long int numeroBarras, DMED *medidas, long int **numeroMedidas, ALIMENTADOR *alimentadores, long int numeroAlimentadores, DRAM *ramos,double Sbase);
int *montaRNP(ALIMENTADOR alimentadores)
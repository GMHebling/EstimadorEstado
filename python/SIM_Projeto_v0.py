#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SIMULADOR MONTE CARLO - ESTIMADOR DE ESTADO TRIFÁSICO

Script para análise estatística do estimador de estado trifásico por simulação de Monte Carlo

1. Monta caso de referência com base em cálculo de fluxo de potência

2. Simulação quasi-estacionária para avaliação temporal

3. Emulação de planos de medição

4. Amostragem de ruído aleatório e Simulação de Monte Carlo 

Created on Sat Mar 21 12:05:29 2020

@author: Julio Massignan
"""
import pandas as pd

def LeituraDados(foldername):
    #----------------------------------------------------------------
    #Função de leitura de dados da rede elétrica e préprocessamento
    # foldername: pasta com arquivos a serem lidos
    
    
    df = pd.DataFrame(columns = ['ID'])
    
    return df

def PowerFlow(network_model, metodo):
    #----------------------------------------------------------------
    #Função que calcula o fluxo de potência para a rede em determinado instnate de tempo
    # Chama rotina externa para o cálculo
    # network_model: modelo da rede elétrica em dataframe
    # metodo: esolha do método 1= Newthon Rapshon e 2= Varredura Direta/Inversa (futuro)
    
           
    #Exporta condicao de carga como plano de medição sem redundância
    
    
    #Calula fluxo de potência pelo método de newton
    
    
    #Leitura do resultado
    
    #Salva resultado
    
    return 0


#-----------------------------------------------------------------------------
# LEITURA DE DADOS
#
#-----------------------------------------------------------------------------
md = '/home/julio/projetos'
sd = '/IEEE342SIM/' 

# Valores de Precisão dos tipos de medidores
prec_PSEUDO = 0.30
prec_SMeter = 0.05
prec_SCADA = 0.02
prec_SCADA_V = 0.01
prec_PMU = 0.001
prec_VIRTUAL = 0.0001

# Parâmetros das simulações
Namostras = 1
Dt = 1
tem_pmu = 0
pmu_polar = 0

network_model = LeituraDados(md + sd)

#-----------------------------------------------------------------------------
# SELEÇÃO DE EVENTOS PARA SIMULAÇÃO QUASI-ESTACIONÁRIA
#
#----------------------------------------------------------------------------- 


#-----------------------------------------------------------------------------
# CASO DE REFERÊNCIA
#
#-----------------------------------------------------------------------------   





#-----------------------------------------------------------------------------
# SIMULAÇÃO DE MONTE CARLO
#
#-----------------------------------------------------------------------------

# Montagem do plano de medição

# Amostragem de ruído aleatório

# Exporta medidas


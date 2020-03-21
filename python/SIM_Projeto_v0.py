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
import subprocess
import numpy as np

def LeituraDados(foldername):
    #----------------------------------------------------------------
    #Função de leitura de dados da rede elétrica e préprocessamento
    # foldername: pasta com arquivos a serem lidos
    
    # Dados de Barras
    filename = '/DBAR.csv'
    try:
        df_DBAR = pd.read_csv(foldername + filename, sep = ',',header=None)
    except:
        df_DBAR = pd.DataFrame(columns = range(1,18))
    df_DBAR.columns = ["ID", "LIGACAO", "FASES", "TENSAO_NOM", "PNOM_A", "PNOM_B", "PNOM_C", "QNOM_A", "QNOM_B", "QNOM_C", \
            "ZIP", "VNOM_A", "VNOM_B", "VNOM_C", "ANG_VNOM_A", "ANG_VNOM_B", "ANG_VNOM_C"]
    
    # Dados de Ramais e Circuitos
    filename = '/DLIN.csv'
    try:
        df_DLIN = pd.read_csv(foldername + filename, sep = ',',header=None)
    except:
        df_DLIN = pd.DataFrame(columns = range(1,23))
    df_DLIN.columns = ["DE", "PARA", "FASES", "COMPRIMENTO", "Raa","Xaa", "Rab","Xab", "Rac","Xac", "Rbb","Xbb",\
            "Rbc","Xbc", "Rcc","Xcc",\
            "Baa", "Bab", "Bac", "Bbb", "Bbc", "Bcc" ]
    
    # Dados de Reguladores de Tensão
    filename = '/DREG.csv'
    try:
        df_DREG = pd.read_csv(foldername + filename, sep = ',',header=None)
    except:
        df_DREG = pd.DataFrame(columns = range(1,27))    
    df_DREG.columns = ["DE", "PARA", "FASES", "VNOM", "REGULACAO", "NTAPS", "S_NOMINAL",
                        "RESISTENCIA","REATANCIA","LIGACAO","RELACAO_TP","RELACAO_TC","DELTA_V","R1","X1","R2","X2","R3",
                        "X3","V1","V2","V3","TIPOCONT","TAP1","TAP2","TAP3"]
    
    # Dados de Bancos de Capacitor
    filename = '/DSHNT.csv'
    try:
        df_DSHNT = pd.read_csv(foldername + filename, sep = ',',header=None)
    except:
        df_DSHNT = pd.DataFrame(columns = range(1,9))
    df_DSHNT.columns = ["ID", "LIGACAO", "FASES","TENSAO_NOM","QNOM_A","QNOM_B","QNOM_C","TipoCont"]
    
    # Dados de Transformadores
    filename = '/DTRF.csv'
    try:
        df_DTRF = pd.read_csv(foldername + filename, sep = ',',header=None)
    except:
        df_DTRF = pd.DataFrame(columns = range(1,14))
    df_DTRF.columns = ["DE", "PARA", "FASES", "VPRI", "VSEC", "S_NOMINAL", "RESISTENCIA", "REATANCIA",\
            "LIGACAO_PRI", "LIGACAO_SEC", "DEFASAMENTO", "TAP_PRI", "TAP_SEC"]
    
        
    return df_DBAR, df_DLIN, df_DREG, df_DSHNT, df_DTRF

def ExportMeasurementSet(md, sd, filename, df_DMED):
    #----------------------------------------------------------------
    #Função que exporta em arquvio DMED.csv o plano de medição e respecitvos valores medidos
    df_DMED.to_csv(md + sd + filename, sep = ',', index=False, header=False, float_format='%.15f')


def PowerFlow(md, sd, network_model, loading, method):
    #----------------------------------------------------------------
    #Função que calcula o fluxo de potência para a rede elétrica
    # Chama rotina externa para o cálculo
    # network_model: modelo da rede elétrica em dataframe
    # loading: cenário de carga em dataframe
    # method: esolha do método 1= Newthon Rapshon e 2= Varredura Direta/Inversa (futuro)
    
    
    
           
    #Exporta condicao de carga como plano de medição sem redundância
    # ExportMeasurementSet(md, sd, "\DMED.csv")
    
    #Calula fluxo de potência pelo método de newton
    subprocess.check_call('./powerflow', cwd = md, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    #Leitura do resultado
    filename = '/referencia.txt'
    df_DREF = pd.read_csv(md + filename, sep = ',',header=None)
    df_DREF.columns = ['Estado','Tipo','De','Para','Fases','Zmed','Sigma']
    
    return df_DREF    


def SampleMeasurementsMC():
    #----------------------------------------------------------------
    #Função que exporta em arquvio DMED.csv o plano de medição e respecitvos valores medidos
    
    
    
    return 0

def StateEstimation(network_model, measurement_set, method):
    #----------------------------------------------------------------
    #Função que roda o estimador de estado para um
    # Chama rotina externa para o cálculo
    # network_model: modelo da rede elétrica em dataframe
    # measurement_set: plano de medição em dataframe
    # method: esolha do método 1= Newthon Rapshon
    
           
    #Exporta medidas
    
    
    #Roda estimador de estado
    
    
    #Leitura do resultado
    
    
    #Salva resultado
    
    
    return 0


#-----------------------------------------------------------------------------
# LEITURA DE DADOS
#
#-----------------------------------------------------------------------------
md = '/home/julio/projetos/EstimadorEstado'
sd = '/IEEE342SIM' 

# Valores de Precisão dos tipos de medidores
precision = {'PSEUDO': 0.30,
             'SMeter': 0.05,
             'SCADA': 0.02,
             'PMU': 0.001,
             'VIRTUAL': 0.0001}

# Parâmetros das simulações
Namostras = 1
Dt = 1
tem_pmu = 0
pmu_polar = 0

df_DBAR, df_DLIN, df_DREG, df_DSHNT, df_DTRF = LeituraDados(md + sd)

# imprime config.txt


#-----------------------------------------------------------------------------
# SELEÇÃO DE EVENTOS PARA SIMULAÇÃO QUASI-ESTACIONÁRIA
#
#----------------------------------------------------------------------------- 

# Define o cenário de carregamento




#-----------------------------------------------------------------------------
# CASO DE REFERÊNCIA
#
#-----------------------------------------------------------------------------   
for t in range(1,Dt+1):
    df_DREF = PowerFlow(md, sd, network_model = 1, loading = 1, method = 1)
    
    print(t)


#-----------------------------------------------------------------------------
# SIMULAÇÃO DE MONTE CARLO
#
#-----------------------------------------------------------------------------
np.random.seed(100)


# Montagem do plano de medição

# Amostragem de ruído aleatório

# Exporta medidas






















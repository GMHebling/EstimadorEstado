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

class sim_data:
    # Classe para armazenar resultados de cálculo do estimador de estado    
    nvar = 0
    nmed = 0
    
    def __init__(self, df_DREF, x, Cov_X, z, residuo):
        self.df_DREF = df_DREF
        self.x = x
        self.Cov_X = Cov_X
        self.z = z
        self.residuo = residuo

class network_data:
    # Classe para salvar dados da rede elétrica no formato de data_frame de leitura
    n_bus = 0
    n_branch = 0
        
    def __init__(self, df_DBAR, df_DLIN, df_DREG, df_DSHNT, df_DTRF):
        self.df_DBAR = df_DBAR
        self.df_DLIN = df_DLIN
        self.df_DREG = df_DREG
        self.df_DSHNT = df_DSHNT
        self.df_DTRF = df_DTRF


def PrintConfig(md,sd):
    # Imprime arquivo de configuração do estimador trifásico
    # md: pasta principal com executável e arquivos de saída
    # sd: pasta do sistema a ser simulado
    file = open(md + '/config.tx','w')
    file.write(md) 
    file.write(sd)
    file.close()
    

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
    
    network_model = network_data(df_DBAR, df_DLIN, df_DREG, df_DSHNT, df_DTRF)
    network_model.n_bus = len(network_model.df_DBAR)
    network_model.n_branch = len(network_model.df_DLIN) + len(network_model.df_DREG) + len(network_model.df_DTRF)   
  
    
    return network_model

    
    return df_DBAR, df_DLIN, df_DREG, df_DSHNT, df_DTRF

def ExportMeasurementSet(md, sd, filename, df_DREF, locMed):
    #----------------------------------------------------------------
    #Função que exporta em arquvio DMED.csv o plano de medição e respecitvos valores medidos
    # 
    # filename: nome do arquivo a ser salvo
    # df_DMED: Valores de referência do plano de medição
    # locMed: plano de medição, local de instalação de cada medidor e respectivo tipo
    
    df_DMED = df_DREF
    
    df_DMED.to_csv(md + sd + filename, sep = ',', index=False, header=False, float_format='%.15f')

def SampleMeasurementsMC():
    #----------------------------------------------------------------
    #Função que exporta em arquvio DMED.csv o plano de medição e respecitvos valores medidos
    
    
    
    return 0

def PowerFlow(md, sd, network_model, loading, method):
    #----------------------------------------------------------------
    #Função que calcula o fluxo de potência para a rede elétrica
    # Chama rotina externa para o cálculo
    # network_model: modelo da rede elétrica em dataframe
    # loading: cenário de carga em dataframe
    # method: esolha do método 1= Newthon Rapshon e 2= Varredura Direta/Inversa (futuro)
    
    # Plano de medição para fluxo de carga - sem redundância seguindo modelo de tipo de barra PQ, PV ou VTeta
    locMed_PF = {'IPQ': network_model.df_DBAR['ID'][1:network_model.n_bus].to_list(),
                'FPQ': [],
                'V': [1]} 
           
    #Exporta condicao de carga como plano de medição sem redundância
    df_DREF = 1
    
    # ExportMeasurementSet(md, sd, "\DMED.csv", df_DREF, locMed_PF)
    
    #Calula fluxo de potência pelo método de newton
    subprocess.check_call('./powerflow', cwd = md, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    #Leitura do resultado
    filename = '/referencia.txt'
    df_DSIM = pd.read_csv(md + filename, sep = ',',header=None)
    df_DSIM.columns = ['Estado','Tipo','De','Para','Fases','Zmed','Sigma']
    
    filename = '/state.txt'
    x = pd.read_csv(md + filename, sep = '\t',header=None)
    x.columns = ['regua','val']
    
    Cov_X = []
    
    residuo = []
    
    z = []
        
    #Salva resultado
    sim_simul = sim_data(df_DREF, x, Cov_X, z, residuo)
    sim_simul.nvar = len(x)
    sim_simul.nmed = len(z)
    
    return sim_simul   


def StateEstimation(md, sd, network_model, measurement_set, method):
    #----------------------------------------------------------------
    #Função que roda o estimador de estado para um
    # Chama rotina externa para o cálculo
    # network_model: modelo da rede elétrica em dataframe
    # measurement_set: plano de medição em dataframe
    # method: esolha do método 1= Newthon Rapshon
    
           
    #Exporta medidas
    
    
    #Roda estimador de estado
    subprocess.check_call('./ss', cwd = md, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    #Leitura do resultado
    filename = '/referencia.txt'
    df_DSIM = pd.read_csv(md + filename, sep = ',',header=None)
    df_DSIM.columns = ['Estado','Tipo','De','Para','Fases','Zmed','Sigma']
    
    filename = '/state.txt'
    x = pd.read_csv(md + filename, sep = ',',header=None)
    x.columns = ['regua','val']
    
    Cov_X = []
    
    filename = '/residuoNormalizado.txt'
    residuo = pd.read_csv(md + filename, sep = ',',header=None)
    residuo.columns = ['id','r','rN','ec','UI']
    
    z = []    
    
    #Salva resultado
    sim_simul = sim_data(df_DSIM, x, Cov_X, z, residuo)
    sim_simul.nvar = len(x)
    sim_simul.nmed = len(z)
    
    return sim_simul


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

# Leitura de Dados
network_model = LeituraDados(md + sd)

# imprime config.txt
PrintConfig(md,sd)

#-----------------------------------------------------------------------------
# SELEÇÃO DE EVENTOS PARA SIMULAÇÃO QUASI-ESTACIONÁRIA
#
#----------------------------------------------------------------------------- 

# Define o cenário de carregamento




#-----------------------------------------------------------------------------
# CASO DE REFERÊNCIA
#
#-----------------------------------------------------------------------------   
sim_ref = []
for t in range(1,Dt+1):
    print(t)
    sim_ref.append(PowerFlow(md, sd, network_model, loading = 1, method = 1))

#-----------------------------------------------------------------------------
# SIMULAÇÃO DE MONTE CARLO
#
#-----------------------------------------------------------------------------
np.random.seed(100)


# Montagem do plano de medição
locMed_SCADA = {'IPQ': [],
                'FPQ': [],
                'V': []}

locMed_PMU = {'ICur': [],
              'FCur': [],
              'V': []}

locMed_Virtual = {'IPQ': [],
                'FPQ': [],
                'V': []}

locMed_SM = {'IPQ': [1001 , 1002 , 1003 , 1005 , 1006 , 1007 , 1008 , 1010 , 1011 , 1012 , 
                     1013 , 1015 , 1016 , 1017 , 1018 , 1019 , 1020 , 1021 , 1022 , 1023 ,
                     1024 , 1025 , 1026 , 1031 , 1032 , 1037 , 1038 , 1039 , 1040 , 1041 ,
                     1042 , 1043 , 1046 , 1047 , 1048 , 1051 , 1052 , 1053 , 1056 , 1057 ,
                     1058 , 1061 , 1062 , 1063 , 1064 , 1065 , 1066 , 1067 , 1072 , 1073 ,
                     1078 , 1079 , 1080 , 1081 , 1082 , 1083 , 1084 , 1087 , 1088 , 1089 ,
                     1092 , 1093 , 1094 , 1097 , 1098 , 1099 , 1102 , 1103 , 1104 , 1105 ,
                     1106 , 1107 , 1108 , 1113 , 1114 , 1120 , 1121 , 1122 , 1123 , 1124 , 
                     1125 , 1127 , 1128 , 1129 , 1130 , 1132 , 1133 , 1134 , 1135 , 1137 , 
                     1138 , 1139 , 1140 , 1142 , 1143 , 1144],
             'FPQ': [],
             'V': []}




# Início da Simulação de Monte Carlo
sim_MC = []
for amostra in range(1,Namostras+1):
    
    # Amostragem de ruído aleatório
    
    
    
    
    
    sim_simul = []
    for t in range(1,Dt+1):
        print(t)
        
        # Exporta medidas
        
        # Executa o Estimador
        # sim_simul.append(StateEstimation(md, sd, network_model, measurement_set, method))
    
    
    sim_MC.append(sim_simul)
    

# Análise do Erro
    
    
# Figuras de Resultados    






















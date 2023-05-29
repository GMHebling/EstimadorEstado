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
import matplotlib.pyplot as plt
from pylab import *

pmu_polar = 0

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
    file = open(md + '/config.txt','w')
    file.write(md + '\n') 
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
        print("DBAR não encontrado!!!")
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


def ExportSystemData(md, sd, network_model):
    #----------------------------------------------------------------
    #Função que exporta arquivos da rede elétrica
    # network_model: classe com dataframes da rede elétrica
    
    network_model.df_DBAR.to_csv(md + sd + "\DBAR.csv", sep = ',', index=False, header=False, float_format='%.15f')
    network_model.df_DSHNT.to_csv(md + sd + "\DSHNT.csv", sep = ',', index=False, header=False, float_format='%.15f')
    network_model.df_DLIN.to_csv(md + sd + "\DLIN.csv", sep = ',', index=False, header=False, float_format='%.15f')
    network_model.df_DREG.to_csv(md + sd + "\DREG.csv", sep = ',', index=False, header=False, float_format='%.15f')
    network_model.df_DSWTC.to_csv(md + sd + "\DSWTC.csv", sep = ',', index=False, header=False, float_format='%.15f')
    network_model.df_DTRF.to_csv(md + sd + "\DTRF.csv", sep = ',', index=False, header=False, float_format='%.15f')
    

def ExportMeasurementSet(md, sd, filename, mode, df_DREF, locMed):
    #----------------------------------------------------------------
    #Função que exporta em arquvio DMED.csv o plano de medição e respecitvos valores medidos
    # 
    # filename: nome do arquivo a ser salvo
    # mode: modo de escrita do arquivo
    # df_DMED: Valores de referência do plano de medição
    # locMed: plano de medição, local de instalação de cada medidor e respectivo tipo
    
    if pmu_polar == 1:
        # Tipos para PMUs em coordenadas polares
        tipos = {'FPQ': [0 , 1],    
                  'IPQ': [2 , 3],
                  'V': [4],
                  'ICur': [],
                  'FCur': [6 , 7],
                  'Vp': [4 , 5]} 
    elif pmu_polar == 0:
        # Tipos para PMUs em coordenadas retangulares
        tipos = {'FPQ': [0 , 1],    
                 'IPQ': [2 , 3],
                 'V': [4],
                 'ICur': [10 , 11],
                 'FCur': [12 , 13],
                 'Vp': [8, 9]} 
    
    
    df_DMED = pd.DataFrame()
    for key in locMed:
        if len(locMed[key]) > 0:
            
            # Testa se valor do localizador é tupla (medidas em ramos) ou int (medidas em barras)
            if type(locMed[key][0]) == tuple:
                df_Aux = df_DREF[df_DREF[['De', 'Para']].apply(tuple, axis=1).isin(locMed[key])]
                df_DMED = df_DMED.append(df_Aux[df_Aux.Tipo.isin(tipos[key])], ignore_index=True)
                
            else: 
                
                df_Aux = df_DREF[df_DREF.De.isin(locMed[key])]
                df_DMED = df_DMED.append(df_Aux[df_Aux.Tipo.isin(tipos[key])], ignore_index=True)
    
    df_DMED.to_csv(md + sd + filename, sep = ',', index=False, header=False, float_format='%.15f', mode=mode)


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
                'V': [network_model.df_DBAR['ID'][0]]} 
         
    #Exporta condicao de carga como plano de medição sem redundância
    ExportMeasurementSet(md, sd, "/DMED.csv", 'w', loading, locMed_PF)
    
    #Calula fluxo de potência pelo método de newton
    #subprocess.check_call('./powerflow', cwd = md, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.check_output('../powerflow')
    
    #Leitura do resultado
    filename = 'referencia.txt'
    #df_DSIM = pd.read_csv(md + filename, sep = ',',header=None)
    df_DSIM = pd.read_csv(filename, sep = ',',header=None)
    df_DSIM.columns = ['Estado','Tipo','De','Para','Fases','Zmed','Sigma']
    
    filename = 'state.txt'
    x = pd.read_csv(filename, sep = '\t',header=None)
    x.columns = ['regua','val']
    
    #x = x[(x['regua'] != -0.01) & (x['regua'] != -0.11) & (x['regua'] != -0.21)]
    
    Cov_X = []
    
    residuo = []
    
    z = []
        
    #Salva resultado
    sim_simul = sim_data(df_DSIM, x, Cov_X, z, residuo)
    sim_simul.nvar = len(x)
    sim_simul.nmed = len(z)
    
    return sim_simul   


def SampleMeasurementsMC(md, sd, filename, mode, df_DREF, locMed, precision):
    #----------------------------------------------------------------
    #Função que exporta em arquvio DMED.csv o plano de medição e respecitvos valores medidos com inserção de ruído
    # 
    # filename: nome do arquivo a ser salvo
    # mode: modo de escrita do arquivo
    # df_DMED: Valores de referência do plano de medição
    # locMed: plano de medição, local de instalação de cada medidor e respectivo tipo
    
    nmed = len(df_DREF)
    
    # Inseri ruído aleatório nos valores de referência
    df_DMED = df_DREF.copy()
    df_DMED['Zmed'] = df_DREF['Zmed'].values + precision * abs(df_DREF['Zmed'].values) / 3 * np.random.randn(nmed)
    
    if precision == 0:
        precision = 0.00001
    df_DMED['Sigma'] = precision
    
    #Exporta condicao de carga como plano de medição sem redundância
    ExportMeasurementSet(md, sd, filename, mode, df_DMED, locMed)
    

def StateEstimation(md, sd, network_model, measurement_set, method):
    #----------------------------------------------------------------
    #Função que roda o estimador de estado para um
    # Chama rotina externa para o cálculo
    # network_model: modelo da rede elétrica em dataframe
    # measurement_set: vetor de medidas em dataframe
    # method: esolha do método 1= Newthon Rapshon
    
           
    
    #Roda estimador de estado
    #subprocess.check_call('./estimator', cwd = md, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.check_output('../estimator')
    #Leitura do resultado
    filename = 'referencia.txt'
    #df_DSIM = pd.read_csv(md + filename, sep = ',',header=None)
    df_DSIM = pd.read_csv(filename, sep = ',',header=None)
    df_DSIM.columns = ['Estado','Tipo','De','Para','Fases','Zmed','Sigma']
    
    filename = 'state.txt'
    x = pd.read_csv(filename, sep = '\t',header=None)
    x.columns = ['regua','val']
    print('state estimation', len(x))
    Cov_X = []
    
    # filename = '/residuoNormalizado.txt'
    # residuo = pd.read_csv(md + filename, sep = ',',header=None)
    # residuo.columns = ['id','r','rN','ec','UI']
    residuo = []
    
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
md = '..'
# sd = '/IEEE342SIM' 
# sd = '/IEEE123' 
# sd = '/IEEE34'
# sd = '/IEEE906'
# sd = '/IEEE906'
sd = '/REAL1058'

# Valores de Precisão dos tipos de medidores
precision = {'PSEUDO': 0.30,
             'SMeter': 0.05,
             'SCADA': 0.02,
             'PMU': 0.001,
             'VIRTUAL': 0.001}

# Parâmetros das simulações
Namostras = 100
Dt = 1
tem_pmu = 0

# Leitura de Dados
network_model = LeituraDados(md + sd)
print("")
# imprime config.txt
#PrintConfig(md,sd)




#-----------------------------------------------------------------------------
# SELEÇÃO DE EVENTOS PARA SIMULAÇÃO QUASI-ESTACIONÁRIA
#
#----------------------------------------------------------------------------- 

# Define o cenário base de carregamento
try:
    filename = '/DMED_fp.csv'
    base_loading =  pd.read_csv(md + sd + filename, sep = ',',header=None)
    base_loading.columns = ['Estado','Tipo','De','Para','Fases','Zmed','Sigma']
    base_loading.Sigma = 1
except: 
    # Não encontrou arquivo com carregamento base para o fluxo de potência
    # Varre o DBAR para preencher o carregamento base - Por enquanto limitado ao modelo PQ constante
    base_loading = pd.DataFrame(columns = ['Estado','Tipo','De','Para','Fases','Zmed','Sigma'])
    
    # Tem que montar um DREF com os valores de P e Q do DBAR
    

# Preenche a estrutrua de carregamento
loading = []
for t in range(1,Dt+1):
    loading.append(base_loading)

# Possibilita inserir contingências, transições, variação de carga e geração, etc.
    
  
    

#-----------------------------------------------------------------------------
# CASO DE REFERÊNCIA
#
#-----------------------------------------------------------------------------   

print('------- Caso de Referência com Fluxo de Potência ----------')
sim_ref = []
for t in range(1,Dt+1):
    print('t: [' + str(t) + '] / ['+ str(Dt) + ']')
    sim_ref.append(PowerFlow(md, sd, network_model, loading[t-1], method = 1))



#-----------------------------------------------------------------------------
# SIMULAÇÃO DE MONTE CARLO
#
#-----------------------------------------------------------------------------
np.random.seed(100)


# Montagem do plano de medição 
# # Plano de Medição para o IEEE34
# locMed_SCADA = {'IPQ': [800,802,806,808,810,812,814,816,818,820,822,824,826,828,830,832,834,836,838,840,842,844,846,848,850,852,854,856,858,860,862,864,888,890,814,852],
#                 'FPQ': [(800, 802), (802, 800), (802, 806), (806, 802), (806, 808), (808, 806),
#                         (808, 810), (808, 812), (810, 808), (812, 808), (812, 814), (814, 812),
#                         (814, 8140), (816, 818), (816, 824), (816, 850), (818, 816), (818, 820),
#                         (820, 818), (820, 822), (822, 820), (824, 816), (824, 826), (824, 828),
#                         (828, 824), (828, 830), (830, 828), (830, 854), (832, 858), (832, 8520),
#                         (832, 888), (834, 860), (834, 842), (834, 858), (835, 840), (836, 862),
#                         (836, 860), (838, 862), (840, 836), (842, 834), (842, 844), (844, 842),
#                         (844, 846), (846, 844), (846, 848), (848, 846), (850, 8140), (850, 816), 
#                         (852, 854), (852,8520), (854, 830), (854,856), (854,852), (856,854), 
#                         (858,832), (858,864), (858,834), (860,834), (860,836), (862,836), 
#                         (864,858), (888,890), (888,832)],
#                 'V': [800,802,806,808,810,812,814,816,818,820,822,824,826,828,830,832,834,836,838,840,842,844,846,848,850,852,854,856,858,860,862,864,888,890,8140,8520]}
# #
# locMed_PMU = {'ICur': [],
#               'FCur': [],
#               'Vp': []}
# locMed_Pseudo = {'IPQ': [],
#                  'FPQ': [],
#                  'V': []}
# locMed_SM = {'IPQ': [],
#              'FPQ': [],
#              'V': []}



#Plano de Medição para o IEEE123
# locMed_SCADA = {'IPQ': [150, 83, 88, 90, 92],
#                'FPQ': [(150, 149), (149, 150), (9, 9000), (25, 2500), (60, 6000), (61, 610),
#                        (18, 35), (13, 52), (97, 101)],
#                'V': [150, 149, 6000]}

# locMed_PMU = {'ICur': [],
#              'FCur': [(150, 149), (60, 6000), (18, 35)],
#              'Vp': [150, 149, 60, 18, 83]}
# locMed_Pseudo = {'IPQ': [1, 2, 4, 5, 6, 7, 9, 12, 10, 11, 16, 17,
#                         19, 20, 22, 24, 28, 29, 30, 31, 32, 33, 34,
#                         35, 37, 38, 39, 41, 42, 43, 45, 46, 47, 48,
#                         49, 50, 51, 52, 53, 55, 56, 58, 59, 60, 62,
#                         63, 64, 65, 66, 68, 69, 70, 71, 73, 74, 75,
#                         76, 77, 79, 80, 82, 84, 85, 86, 87, 94, 95,
#                         96, 98, 99, 100, 102, 103, 104, 106, 107, 109,
#                         111, 112, 113, 114],
#                 'FPQ': [],
#                 'V': []}
# locMed_SM = {'IPQ': [],
#             'FPQ': [],
#             'V': []}

# Plano de Medição para o IEEE342
# locMed_SCADA = {'IPQ': [1],
#                'FPQ': [(1, 2), (2, 1), (2, 3), (2, 7), (3, 5), (5, 3), (7, 9), (9, 7),
#                        (5, 44), (9, 116), (5, 11), (5, 27), (44,
#                                                              45), (44, 62), (9, 81), (9, 97),
#                        (116, 117), (116, 134), (11, 13), (27,
#                                                           29), (45, 47), (62, 64), (81, 83),
#                        (97, 99), (117, 119), (134, 136)],
#                'V': [1, 3, 5, 7, 9,
#                      11, 27, 45, 62, 81, 97, 117, 134,
#                      1001, 1196, 1201, 1207, 1214, 1221, 1228, 1234, 1238]}

# locMed_PMU = {'ICur': [],
#              'FCur': [(1, 2), (3, 5), (5, 3), (7, 9), (9, 7)],
#              'Vp': [1, 3, 5, 7, 9]}

# locMed_Pseudo = {'IPQ': [1193, 1198, 1203, 1210, 1217, 1224, 1231, 1236,
#                         1001, 1002, 1003, 1005, 1006, 1007, 1008, 1010, 1011, 1012,
#                         1013, 1015, 1016, 1017, 1018, 1019, 1020, 1021, 1022, 1023,
#                         1024, 1025, 1026, 1031, 1032, 1037, 1038, 1039, 1040, 1041,
#                         1042, 1043, 1046, 1047, 1048, 1051, 1052, 1053, 1056, 1057,
#                         1058, 1061, 1062, 1063, 1064, 1065, 1066, 1067, 1072, 1073,
#                         1078, 1079, 1080, 1081, 1082, 1083, 1084, 1087, 1088, 1089,
#                         1092, 1093, 1094, 1097, 1098, 1099, 1102, 1103, 1104, 1105,
#                         1106, 1107, 1108, 1113, 1114, 1120, 1121, 1122, 1123, 1124,
#                         1125, 1127, 1128, 1129, 1130, 1132, 1133, 1134, 1135, 1137,
#                         1138, 1139, 1140, 1142, 1143, 1144],
#                 'FPQ': [],
#                 'V': []}

# locMed_SM = {'IPQ': [],
#             'FPQ': [],
#             'V': []}

# Plano de Medição para o IEEE906
# locMed_SCADA = {'IPQ': [1,2,3,5,7,9,11,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150],
#                'FPQ': [(1, 2), (2, 3), (15, 16)],
#                'V': [1, 3, 5, 7, 9, 11, 34, 66, 123]}

#ajuste para o hatchel
# locMed_SCADA = {'IPQ': [34,70,225,289,349,387,388,502,562,563,611,629,817,860,861,896,898,900,906,47,83,178,248,249,276,314,406,522,639,676,682,688,702,755,785,813,886,899,208,264,320,327,337,342,458,539,556,614,619,701,778,780,835,73,74],
#                'FPQ': [],
#                'V': [0]}

# locMed_PMU = {'ICur': [],
#              'FCur': [],
#              'Vp': []}

# locMed_Pseudo = {'IPQ': [],
#                 'FPQ': [],
#                 'V': []}

# locMed_SM = {'IPQ': [],
#             'FPQ': [],
#             'V': []}

### 
# Plano de Medição para o IEEE906
# locMed_SCADA = {'IPQ': [1,2,3,5,7,9,11,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150],
#                'FPQ': [(1, 2), (2, 3), (15, 16)],
#                'V': [1, 3, 5, 7, 9, 11, 34, 66, 123]}

#ajuste para o hatchel
locMed_SCADA = {'IPQ': [0],
               'FPQ': [(1,0), (17,18), (20,21), (41,42), (0,1)],
               'V': [1, 0]}

locMed_PMU = {'ICur': [],
             'FCur': [],
             'Vp': []}

locMed_Pseudo = {'IPQ': [803, 799, 798, 797, 784, 779, 774, 766, 765, 809, 925, 908, 898, 884, 880, 865, 856, 838, 824, 761, 705, 697, 695, 689, 681, 675, 671, 664, 661, 714, 718, 758, 748, 745, 731, 729, 728, 725, 722, 973, 1248, 1239, 1232, 1224, 1221, 1202, 1194, 1189, 1180, 1258, 1265, 1315, 1311, 1306, 1303, 1297, 1294, 1287, 1281, 657, 1174, 1167, 1050, 1030, 1025, 1022, 1018, 1005, 998, 989, 984, 1055, 1158, 1154, 1097, 1092, 1082, 1079, 1069, 1066, 27, 656, 310, 304, 297, 289, 281, 279, 255, 231, 222, 217, 319, 381, 375, 366, 361, 354, 348, 343, 330, 188, 85, 79, 56, 52, 48, 42, 38, 34, 91, 393, 100, 179, 176, 174, 162, 159, 137, 123, 117, 112, 405, 595, 586, 576, 562, 559, 555, 547, 545, 536, 599, 609, 646, 644, 640, 634, 631, 620, 618, 527, 479, 476, 472, 464, 445, 436, 435, 430, 408, 483, 522, 518, 515, 510, 507, 502, 498, 488],
                'FPQ': [],
                'V': []}

locMed_SM = {'IPQ': [],
            'FPQ': [],
            'V': []}

# Injecoes virtuais - barras sem carga que não estao no vetor de medidas SCADA, Pseudo e SMeter
injecoes = locMed_Pseudo['IPQ'] + locMed_SCADA['IPQ'] + locMed_SM['IPQ']
aux_inj = network_model.df_DBAR['ID'].values
# locMed_Virtual = {'IPQ': list(set(aux_inj)-set(injecoes)),
#                   'FPQ': [],
#                   'V': []}
# {'IPQ': [1090, 1313, 1310, 1309, 1305, 1304, 1302, 1301, 1314, 5, 4, 3, 2, 1, 1316, 1300, 1289, 1286, 1285, 1284, 1283, 1280, 1279, 1290, 1299, 1296, 1293, 1292, 1291, 6, 28, 26, 25, 24, 23, 22, 21, 29, 39, 37, 36, 35, 33, 31, 20, 13, 12, 11, 10, 9, 7, 14, 19, 18, 17, 16, 15, 1278, 1215, 1214, 1213, 1212, 1211, 1210, 1216, 1217, 1223, 1222, 1220, 1219, 1218, 1209, 1200, 1199, 1198, 1197, 1196, 1195, 1201, 1208, 1207, 1206, 1225, 1205, 1204, 1203, 1228, 1259, 1257, 1256, 1249, 1247, 1246, 1245, 1271, 1277, 1276, 1275, 1274, 1272, 1244, 1235, 1233, 1231, 1230, 1229, 1236, 1237, 1243, 1242, 1241, 1240, 1238, 1192, 40, 145, 144, 143, 142, 140, 139, 136, 146, 153, 152, 151, 149, 148, 147, 133, 122, 121, 120, 119, 116, 115, 114, 124, 130, 129, 128, 127, 125, 154, 190, 187, 186, 185, 184, 183, 191, 192, 199, 198, 196, 194, 193, 182, 168, 167, 165, 164, 156, 155, 169, 181, 178, 173, 172, 171, 170, 113, 64, 63, 62, 61, 60, 59, 66, 68, 73, 72, 71, 70, 69, 58, 47, 46, 45, 44, 43, 41, 49, 57, 55, 54, 74, 53, 51, 50, 75, 104, 103, 101, 99, 98, 96, 95, 106, 111, 110, 109, 108, 107, 93, 82, 80, 78, 77, 76, 83, 84, 92, 90, 89, 88, 86, 1191, 975, 974, 970, 969, 968, 967, 978, 980, 986, 985, 983, 982, 981, 987, 966, 957, 956, 955, 954, 952, 951, 959, 960, 965, 964, 963, 962, 961, 988, 1019, 1017, 1016, 1015, 1014, 1012, 1011, 1020, 1034, 1033, 1032, 1031, 1029, 1021, 1010, 996, 995, 994, 993, 992, 1000, 1001, 1009, 1008, 1007, 1003, 1002, 950, 911, 910, 909, 907, 906, 905, 904, 912, 918, 917, 916, 915, 914, 913, 903, 895, 894, 893, 891, 890, 889, 888, 896, 902, 901, 900, 899, 897, 919, 941, 940, 939, 938, 937, 936, 942, 944, 949, 948, 947, 946, 945, 935, 926, 924, 923, 922, 921, 920, 927, 933, 932, 931, 930, 929, 928, 1035, 1138, 1137, 1136, 1135, 1132, 1131, 1130, 1139, 1146, 1145, 1144, 1143, 1141, 1140, 1129, 1121, 1120, 1119, 1118, 1117, 1116, 1122, 1123, 1128, 1127, 1126, 1125, 1124, 1147, 1183, 1178, 1176, 1175, 1173, 1172, 1170, 1184, 1190, 1188, 1187, 1186, 1185, 1169, 1159, 1157, 1155, 1153, 1152, 1151, 1161, 1166, 1165, 1164, 1163, 1162, 200, 1115, 1064, 1062, 1061, 1060, 1059, 1058, 1052, 1065, 1073, 1072, 1071, 1070, 1068, 1067, 1049, 1042, 1041, 1039, 1038, 1037, 1036, 1043, 1048, 1047, 1046, 1045, 1044, 1074, 1106, 1105, 1098, 1096, 1095, 1094, 1107, 1108, 1114, 1113, 1112, 1111, 1110, 1093, 1083, 1080, 1078, 1077, 1076, 1075, 1085, 1091, 1089, 1088, 1087, 1086, 201, 643, 641, 639, 638, 637, 635, 645, 648, 655, 654, 653, 652, 649, 660, 633, 622, 621, 619, 617, 616, 615, 623, 624, 630, 629, 628, 627, 626, 663, 694, 693, 692, 691, 690, 688, 687, 699, 708, 707, 706, 704, 703, 702, 686, 674, 673, 670, 669, 666, 677, 678, 685, 684, 683, 682, 680, 614, 569, 568, 567, 566, 565, 564, 563, 570, 577, 575, 574, 573, 572, 571, 561, 549, 548, 546, 544, 543, 542, 541, 550, 558, 557, 554, 553, 551, 578, 604, 602, 601, 600, 598, 596, 605, 606, 613, 612, 611, 610, 608, 594, 585, 584, 583, 581, 580, 579, 587, 593, 592, 591, 590, 589, 588, 709, 855, 853, 852, 851, 850, 849, 857, 859, 864, 863, 862, 861, 860, 866, 848, 836, 835, 833, 832, 831, 830, 837, 839, 847, 846, 844, 842, 841, 867, 887, 885, 883, 872, 871, 870, 869, 868, 873, 874, 881, 879, 878, 876, 875, 829, 755, 754, 753, 752, 751, 742, 741, 756, 775, 773, 772, 771, 770, 764, 739, 720, 715, 713, 712, 711, 710, 721, 737, 730, 727, 724, 723, 777, 820, 819, 818, 817, 816, 815, 821, 822, 828, 827, 826, 825, 823, 813, 801, 800, 796, 794, 783, 782, 802, 812, 810, 808, 807, 806, 805, 540, 308, 307, 303, 300, 299, 296, 311, 295, 312, 317, 316, 315, 314, 313, 318, 294, 285, 284, 283, 278, 277, 276, 286, 293, 292, 291, 290, 288, 320, 347, 346, 345, 342, 341, 340, 339, 350, 359, 358, 357, 355, 353, 352, 338, 328, 327, 326, 323, 321, 329, 332, 337, 336, 335, 275, 334, 333, 274, 227, 226, 225, 224, 221, 220, 219, 228, 236, 235, 234, 233, 230, 229, 215, 206, 205, 204, 203, 202, 207, 208, 213, 212, 211, 210, 209, 237, 262, 261, 260, 259, 258, 254, 253, 263, 273, 272, 269, 267, 265, 252, 243, 242, 241, 240, 239, 238, 244, 251, 249, 247, 246, 245, 360, 480, 478, 477, 475, 474, 473, 471, 481, 494, 487, 486, 485, 484, 482, 469, 460, 459, 458, 457, 455, 454, 461, 462, 468, 467, 466, 465, 463, 495, 525, 524, 521, 520, 519, 517, 516, 526, 539, 538, 537, 535, 534, 514, 504, 501, 500, 499, 497, 496, 505, 513, 512, 511, 509, 508, 506, 453, 390, 389, 388, 387, 386, 385, 382, 391, 399, 398, 397, 396, 394, 392, 380, 370, 367, 365, 364, 363, 362, 371, 379, 378, 376, 374, 373, 400, 433, 432, 431, 429, 428, 426, 434, 437, 452, 451, 450, 449, 448, 425, 407, 406, 404, 403, 402, 401, 409, 417, 415, 414, 413, 412, 411, 1315, 671, 714, 695, 884, 765, 304, 289, 179, 162, 137, 27, 393, 436, 408, 405, 366, 595, 1174, 1221, 1154, 52],
locMed_Virtual = {'IPQ': [1090, 1313, 1310, 1309, 1305, 1304, 1302, 1301, 1314, 5, 4, 3, 2, 1, 1316, 1300, 1289, 1286, 1285, 1284, 1283, 1280, 1279, 1290, 1299, 1296, 1293, 1292, 1291, 6, 28, 26, 25, 24, 23, 22, 21, 29, 39, 37, 36, 35, 33, 31, 20, 13, 12, 11, 10, 9, 7, 14, 19, 18, 17, 16, 15, 1278, 1215, 1214, 1213, 1212, 1211, 1210, 1216, 1217, 1223, 1222, 1220, 1219, 1218, 1209, 1200, 1199, 1198, 1197, 1196, 1195, 1201, 1208, 1207, 1206, 1225, 1205, 1204, 1203, 1228, 1259, 1257, 1256, 1249, 1247, 1246, 1245, 1271, 1277, 1276, 1275, 1274, 1272, 1244, 1235, 1233, 1231, 1230, 1229, 1236, 1237, 1243, 1242, 1241, 1240, 1238, 1192, 40, 145, 144, 143, 142, 140, 139, 136, 146, 153, 152, 151, 149, 148, 147, 133, 122, 121, 120, 119, 116, 115, 114, 124, 130, 129, 128, 127, 125, 154, 190, 187, 186, 185, 184, 183, 191, 192, 199, 198, 196, 194, 193, 182, 168, 167, 165, 164, 156, 155, 169, 181, 178, 173, 172, 171, 170, 113, 64, 63, 62, 61, 60, 59, 66, 68, 73, 72, 71, 70, 69, 58, 47, 46, 45, 44, 43, 41, 49, 57, 55, 54, 74, 53, 51, 50, 75, 104, 103, 101, 99, 98, 96, 95, 106, 111, 110, 109, 108, 107, 93, 82, 80, 78, 77, 76, 83, 84, 92, 90, 89, 88, 86, 1191, 975, 974, 970, 969, 968, 967, 978, 980, 986, 985, 983, 982, 981, 987, 966, 957, 956, 955, 954, 952, 951, 959, 960, 965, 964, 963, 962, 961, 988, 1019, 1017, 1016, 1015, 1014, 1012, 1011, 1020, 1034, 1033, 1032, 1031, 1029, 1021, 1010, 996, 995, 994, 993, 992, 1000, 1001, 1009, 1008, 1007, 1003, 1002, 950, 911, 910, 909, 907, 906, 905, 904, 912, 918, 917, 916, 915, 914, 913, 903, 895, 894, 893, 891, 890, 889, 888, 896, 902, 901, 900, 899, 897, 919, 941, 940, 939, 938, 937, 936, 942, 944, 949, 948, 947, 946, 945, 935, 926, 924, 923, 922, 921, 920, 927, 933, 932, 931, 930, 929, 928, 1035, 1138, 1137, 1136, 1135, 1132, 1131, 1130, 1139, 1146, 1145, 1144, 1143, 1141, 1140, 1129, 1121, 1120, 1119, 1118, 1117, 1116, 1122, 1123, 1128, 1127, 1126, 1125, 1124, 1147, 1183, 1178, 1176, 1175, 1173, 1172, 1170, 1184, 1190, 1188, 1187, 1186, 1185, 1169, 1159, 1157, 1155, 1153, 1152, 1151, 1161, 1166, 1165, 1164, 1163, 1162, 200, 1115, 1064, 1062, 1061, 1060, 1059, 1058, 1052, 1065, 1073, 1072, 1071, 1070, 1068, 1067, 1049, 1042, 1041, 1039, 1038, 1037, 1036, 1043, 1048, 1047, 1046, 1045, 1044, 1074, 1106, 1105, 1098, 1096, 1095, 1094, 1107, 1108, 1114, 1113, 1112, 1111, 1110, 1093, 1083, 1080, 1078, 1077, 1076, 1075, 1085, 1091, 1089, 1088, 1087, 1086, 201, 643, 641, 639, 638, 637, 635, 645, 648, 655, 654, 653, 652, 649, 660, 633, 622, 621, 619, 617, 616, 615, 623, 624, 630, 629, 628, 627, 626, 663, 694, 693, 692, 691, 690, 688, 687, 699, 708, 707, 706, 704, 703, 702, 686, 674, 673, 670, 669, 666, 677, 678, 685, 684, 683, 682, 680, 614, 569, 568, 567, 566, 565, 564, 563, 570, 577, 575, 574, 573, 572, 571, 561, 549, 548, 546, 544, 543, 542, 541, 550, 558, 557, 554, 553, 551, 578, 604, 602, 601, 600, 598, 596, 605, 606, 613, 612, 611, 610, 608, 594, 585, 584, 583, 581, 580, 579, 587, 593, 592, 591, 590, 589, 588, 709, 855, 853, 852, 851, 850, 849, 857, 859, 864, 863, 862, 861, 860, 866, 848, 836, 835, 833, 832, 831, 830, 837, 839, 847, 846, 844, 842, 841, 867, 887, 885, 883, 872, 871, 870, 869, 868, 873, 874, 881, 879, 878, 876, 875, 829, 755, 754, 753, 752, 751, 742, 741, 756, 775, 773, 772, 771, 770, 764, 739, 720, 715, 713, 712, 711, 710, 721, 737, 730, 727, 724, 723, 777, 820, 819, 818, 817, 816, 815, 821, 822, 828, 827, 826, 825, 823, 813, 801, 800, 796, 794, 783, 782, 802, 812, 810, 808, 807, 806, 805, 540, 308, 307, 303, 300, 299, 296, 311, 295, 312, 317, 316, 315, 314, 313, 318, 294, 285, 284, 283, 278, 277, 276, 286, 293, 292, 291, 290, 288, 320, 347, 346, 345, 342, 341, 340, 339, 350, 359, 358, 357, 355, 353, 352, 338, 328, 327, 326, 323, 321, 329, 332, 337, 336, 335, 275, 334, 333, 274, 227, 226, 225, 224, 221, 220, 219, 228, 236, 235, 234, 233, 230, 229, 215, 206, 205, 204, 203, 202, 207, 208, 213, 212, 211, 210, 209, 237, 262, 261, 260, 259, 258, 254, 253, 263, 273, 272, 269, 267, 265, 252, 243, 242, 241, 240, 239, 238, 244, 251, 249, 247, 246, 245, 360, 480, 478, 477, 475, 474, 473, 471, 481, 494, 487, 486, 485, 484, 482, 469, 460, 459, 458, 457, 455, 454, 461, 462, 468, 467, 466, 465, 463, 495, 525, 524, 521, 520, 519, 517, 516, 526, 539, 538, 537, 535, 534, 514, 504, 501, 500, 499, 497, 496, 505, 513, 512, 511, 509, 508, 506, 453, 390, 389, 388, 387, 386, 385, 382, 391, 399, 398, 397, 396, 394, 392, 380, 370, 367, 365, 364, 363, 362, 371, 379, 378, 376, 374, 373, 400, 433, 432, 431, 429, 428, 426, 434, 437, 452, 451, 450, 449, 448, 425, 407, 406, 404, 403, 402, 401, 409, 417, 415, 414, 413, 412, 411, 1315, 671, 714, 695, 884, 765, 304, 289, 179, 162, 137, 27, 393, 436, 408, 405, 366, 595, 1174, 1221, 1154, 52],
                  'FPQ': [],
                  'V': []}

# --------------------------------------
# Início da Simulação de Monte Carlo

print(sim_ref[0].df_DREF)

print('------- Simulação de Monte Carlo ----------')
sim_MC = []
for amostra in range(1,Namostras+1):
    print('Amostra: [' + str(amostra) + '] / ['+ str(Namostras) + ']')
    sim_simul = []
    for t in range(1,Dt+1):
        print('t: [' + str(t) + '] / ['+ str(Dt) + ']')
                
        # Exporta medidas
        SampleMeasurementsMC(md, sd, "/DMED.csv", 'w', sim_ref[t-1].df_DREF, locMed_Pseudo, precision['PSEUDO'])
        SampleMeasurementsMC(md, sd, "/DMED.csv", 'a', sim_ref[t-1].df_DREF, locMed_SM, precision['SMeter'])
        SampleMeasurementsMC(md, sd, "/DMED.csv", 'a', sim_ref[t-1].df_DREF, locMed_Virtual, precision['VIRTUAL'])
        SampleMeasurementsMC(md, sd, "/DMED.csv", 'a', sim_ref[t-1].df_DREF, locMed_SCADA, precision['SCADA'])        
        # SampleMeasurementsMC(md, sd, "/DMED.csv", 'a', sim_ref[t-1].df_DREF, locMed_PMU, precision['PMU'])
        
        # Executa o Estimador
        sim_simul.append(StateEstimation(md, sd, network_model, measurement_set = 1, method = 1))
    
    
    sim_MC.append(sim_simul)
    


# -----------------------------------------------------------------------------
# ANÁLISES DOS RESULTADOS
# 
# 
# -----------------------------------------------------------------------------
# Análise do Erro
t = 0    
erro_x = []
for amostra in range(0,Namostras):
    erro_x.append(sim_MC[amostra][t].x.val.values - sim_ref[t].x.val.values[:len(sim_ref[t].x.val.values)])

    
# MAE
MAE_t = abs(np.nanmean(erro_x,0))
MAE_x = abs(np.nanmean(erro_x,1))

print('MAE: ', max(MAE_x))

# RMSE
aux = [i ** 2 for i in erro_x]
RMSE_t = np.sqrt(np.nanmean(aux,0))
print('MAE: ', max(RMSE_t))

# KEMA







# -----------------------------------------------------------------------------
# Figuras de Resultados    
params = {
   'axes.labelsize': 8,
   'font.size': 8,
   'legend.fontsize': 10,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'text.usetex': False,
   'figure.figsize': [8.5, 4.5]
   }
rcParams.update(params)


fig = plt.figure()
axes(frameon=0)
grid()

plot(MAE_x, linewidth=2, linestyle='--', color='#B22400', label = 'MAE performance index')
#plot(MAE_x, linewidth=2, linestyle='--', color='#B22400', label = 'MAE performance index')
#plot(RMSE_t, linewidth=2, linestyle='--', color='#006BB2', label = 'RMSE performance index')
# xlim(-5, 400)
# ylim(-5000, 300)
plt.legend(frameon=False)
leg = plt.legend()
frame = leg.get_frame()
frame.set_facecolor('0.9')
frame.set_edgecolor('0.9')

plt.show()






















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
    filename = '/referencia.txt'
    # df_DSIM = pd.read_csv(md + filename, sep = ',',header=None)
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
    filename = '/referencia.txt'
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
sd = '/IEEE123' 
#sd = '/IEEE34'

# Valores de Precisão dos tipos de medidores
precision = {'PSEUDO': 0.30,
             'SMeter': 0.05,
             'SCADA': 0.02,
             'PMU': 0.001,
             'VIRTUAL': 0.000}

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
locMed_SCADA = {'IPQ': [150, 83, 88, 90, 92],
               'FPQ': [(150, 149), (149, 150), (9, 9000), (25, 2500), (60, 6000), (61, 610),
                       (18, 35), (13, 52), (97, 101)],
               'V': [150, 149, 6000]}

locMed_PMU = {'ICur': [],
             'FCur': [(150, 149), (60, 6000), (18, 35)],
             'Vp': [150, 149, 60, 18, 83]}
locMed_Pseudo = {'IPQ': [1, 2, 4, 5, 6, 7, 9, 12, 10, 11, 16, 17,
                        19, 20, 22, 24, 28, 29, 30, 31, 32, 33, 34,
                        35, 37, 38, 39, 41, 42, 43, 45, 46, 47, 48,
                        49, 50, 51, 52, 53, 55, 56, 58, 59, 60, 62,
                        63, 64, 65, 66, 68, 69, 70, 71, 73, 74, 75,
                        76, 77, 79, 80, 82, 84, 85, 86, 87, 94, 95,
                        96, 98, 99, 100, 102, 103, 104, 106, 107, 109,
                        111, 112, 113, 114],
                'FPQ': [],
                'V': []}
locMed_SM = {'IPQ': [],
            'FPQ': [],
            'V': []}

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

#
# Injecoes virtuais - barras sem carga que não estao no vetor de medidas SCADA, Pseudo e SMeter
injecoes = locMed_Pseudo['IPQ'] + locMed_SCADA['IPQ'] + locMed_SM['IPQ']
aux_inj = network_model.df_DBAR['ID'].values
locMed_Virtual = {'IPQ': list(set(aux_inj)-set(injecoes)),
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
    erro_x.append(sim_MC[amostra][t].x.val.values - sim_ref[t].x.val.values)

    
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






















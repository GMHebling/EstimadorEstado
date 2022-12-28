CC=gcc
CFLAGS= -g -o -W
LDFLAGS= -lcholmod -lspqr -lsuitesparseconfig -lm -lstdc++
OBJFILES = main.o funcoesAMB.o funcoesBranchCurrent.o funcoesBadData.o funcoesCalculoEletrico.o funcoesLeitura.o funcoesMatematicas.o funcoesOtimizacao.o funcoesTopologia.o funcoesWLS.o funcoesFluxoVarredura.o 
OBJFILESFP = main_fp.o funcoesBadData.o funcoesCalculoEletrico.o funcoesLeitura.o funcoesMatematicas.o funcoesOtimizacao.o funcoesTopologia.o funcoesWLS.o
TARGET = estimator powerflow

all: $(TARGET)

estimator: $(OBJFILES)
		$(CC) $(CFLAGS) -o estimatorFIXED $(OBJFILES) $(LDFLAGS)

powerflow: $(OBJFILESFP)
		$(CC) $(CFLAGS) -o powerflow $(OBJFILESFP) $(LDFLAGS)
clean:
	rm -f $(OBJFILES) $(OBJFILESFP) *~ 

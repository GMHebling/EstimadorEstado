CC=gcc
CFLAGS= -o -W 
LDFLAGS= -lcholmod -lspqr -lsuitesparseconfig -lm -lstdc++
OBJFILES = main.o funcoesBranchCurrent.o funcoesBadData.o funcoesCalculoEletrico.o funcoesLeitura.o funcoesMatematicas.o funcoesOtimizacao.o funcoesTopologia.o funcoesWLS.o
OBJFILESFP = main_fp.o funcoesBadData.o funcoesCalculoEletrico.o funcoesLeitura.o funcoesMatematicas.o funcoesOtimizacao.o funcoesTopologia.o funcoesWLS.o
TARGET = estimator powerflow clean

all: $(TARGET)

estimator: $(OBJFILES)
		$(CC) $(CFLAGS) -o estimator $(OBJFILES) $(LDFLAGS)

powerflow: $(OBJFILESFP)
		$(CC) $(CFLAGS) -o powerflow $(OBJFILESFP) $(LDFLAGS)
clean:
	rm -f $(OBJFILES) $(OBJFILESFP) *~ 
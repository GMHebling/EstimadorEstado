CC=gcc
CFLAGS= -o -W 
LDFLAGS= -lcholmod -lspqr -lsuitesparseconfig -lm -lstdc++
OBJFILES = main.o funcoesBadData.o funcoesCalculoEletrico.o funcoesLeitura.o funcoesMatematicas.o funcoesOtimizacao.o funcoesTopologia.o funcoesWLS.o
OBJFILESFP = main_fp.o funcoesBadData.o funcoesCalculoEletrico.o funcoesLeitura.o funcoesMatematicas.o funcoesOtimizacao.o funcoesTopologia.o funcoesWLS.o
TARGET = necfixed powerflow 

all: $(TARGET)

necfixed: $(OBJFILES)
		$(CC) $(CFLAGS) -o necfixed $(OBJFILES) $(LDFLAGS)

powerflow: $(OBJFILESFP)
		$(CC) $(CFLAGS) -o powerflow $(OBJFILESFP) $(LDFLAGS)
clean:
	rm -f $(OBJFILES) $(OBJFILESFP) *~ 
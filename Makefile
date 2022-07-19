CC=gcc
CFLAGS= -o -W -g -Wall
LDFLAGS= -lcholmod -lspqr -lsuitesparseconfig -lm -lstdc++ -lumfpack
OBJFILES = main.o funcoesCalculoEletrico.o funcoesLeitura.o funcoesMatematicas.o funcoesOtimizacao.o funcoesTopologia.o funcoesWLS.o Observabilidade.o numref.o matriz.o leitura.o fluxoNRQR.o 
TARGET = ss

all: $(TARGET)

$(TARGET): $(OBJFILES)
		$(CC) $(CFLAGS) -o $(TARGET) $(OBJFILES) $(LDFLAGS)

clean:
	rm -f $(OBJFILES) $(TARGET) *~ 
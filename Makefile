CC=gcc
CFLAGS= -o -W 
LDFLAGS= -lcholmod -lspqr -lsuitesparseconfig -lm -lstdc++
OBJFILES = main.o funcoesBadData.o funcoesCalculoEletrico.o funcoesLeitura.o funcoesMatematicas.o funcoesOtimizacao.o funcoesTopologia.o funcoesWLS.o
TARGET = ss

all: build clean 
	

build: $(TARGET)

$(TARGET): $(OBJFILES)
		$(CC) $(CFLAGS) -o $(TARGET) $(OBJFILES) $(LDFLAGS)

clean:
	clear && rm -f $(OBJFILES) *~ 
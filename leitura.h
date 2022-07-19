#ifndef leitura_H
#define leitura_H

char *config();
int r_dmatfcsv(double ***mat,char path[100],char sep,int *nlin,int *ncol);
void imprimirmat(double** mat, int linha, int col);
void fimprimirmat(FILE* arquivo, double** mat, int linha, int col);
void imprimirvet_double(double* vet, int linha);
void fimprimirvet_double(FILE* arquivo, double* vet, int linha);
void imprimirvet_int(int* vet, int linha);
void fimprimirvet_int(FILE* arquivo, int* vet, int linha);


#endif
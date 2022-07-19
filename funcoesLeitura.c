#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <ctype.h>

#include "data_structures.h"
#include "funcoesLeitura.h"

//------------------------------------------------------------------------------
//
// FUNÇÕES DE LEITURA DE DADOS
//
//------------------------------------------------------------------------------
// Troca vírgula por ponto em string
char *replace(char *st)
{
    int i = 0;
    while (st[i] != '\0')
    {
        if (st[i] == ',')
            st[i] = '.';
        i++;
    }
    return st;
}

// Lê cada campo de uma string separada por vírgula
char *getfield(char *lin, int num)
{
    char *tok;
    char *line;
    line = strdup(lin);
    //printf("\nteste -  %s",line);
    tok = strtok(line, ",\t\n\r");
    while (tok != NULL)
    {
        if (!--num)
        {
            return tok;
        }
        tok = strtok(NULL, ",\t\n\r");
    }
    return NULL;
}

char *charLigacao(LIGACAO num)
{
    char *tok = (char *)malloc(5 * sizeof(char));
    switch (num)
    {
    case O:
        strcpy(tok, "O");
        return (tok);
        break;
    case YN:
        strcpy(tok, "Yn");
        return (tok);
        break;
    case Y:
        strcpy(tok, "Y");
        return (tok);
        break;
    case D:
        strcpy(tok, "D");
        return (tok);
        break;
    case OY:
        strcpy(tok, "OY");
        return (tok);
        break;
    case OD:
        strcpy(tok, "OD");
        return (tok);
        break;
    }
}
char *charFases(FASES num)
{
    char *tok = (char *)malloc(5 * sizeof(char));
    switch (num)
    {
    case N:
        strcpy(tok, "N");
        return (tok);
        break;
    case A:
        strcpy(tok, "A");
        return (tok);
        break;
    case B:
        strcpy(tok, "B");
        return (tok);
        break;
    case C:
        strcpy(tok, "C");
        return (tok);
        break;
    case AB:
        strcpy(tok, "AB");
        return (tok);
        break;
    case BC:
        strcpy(tok, "BC");
        return (tok);
        break;
    case CA:
        strcpy(tok, "CA");
        return (tok);
        break;
    case ABC:
        strcpy(tok, "ABC");
        return (tok);
        break;
    }
}

char *charMedidor(long int num)
{
    char *tok = (char *)malloc(10 * sizeof(char));
    switch (num)
    {
    case 0:
        strcpy(tok, "Pkm");
        return (tok);
        break;
    case 1:
        strcpy(tok, "Qkm");
        return (tok);
        break;
    case 2:
        strcpy(tok, "Pk");
        return (tok);
        break;
    case 3:
        strcpy(tok, "Qk");
        return (tok);
        break;
    case 4:
        strcpy(tok, "Vk");
        return (tok);
        break;
    case 5:
        strcpy(tok, "Tk");
        return (tok);
        break;
    case 6:
        strcpy(tok, "Ikm");
        return (tok);
        break;
    case 7:
        strcpy(tok, "Akm");
        return (tok);
        break;
    case 8:
        strcpy(tok, "Vk_real");
        return (tok);
        break;
    case 9:
        strcpy(tok, "Vk_imag");
        return (tok);
        break;
    case 10:
        strcpy(tok, "Ik_real");
        return (tok);
        break;
    case 11:
        strcpy(tok, "Ik_imag");
        return (tok);
        break;
    case 12:
        strcpy(tok, "Ikm_real");
        return (tok);
        break;
    case 13:
        strcpy(tok, "Ikm_imag");
        return (tok);
        break;
    }
}

//------------------------------------------------------------------------------
//
// FUNÇÕES DE LEITURA DE DADOS
//
//------------------------------------------------------------------------------
// Leitura de dados da rede elétrica
char *leituraDados(DBAR **barra, DRAM **ramo, long int *numeroBarras, long int *numeroRamos, long int *numeroAlimentadores)
{

    // Funcao que lê os dados da rede elétrica, recebe as estruturas DBAR e DRAM (dados barras e dados dos ramos)
    // Recebe também  o número de barras, o número de ramos e o número de alimentadores
    // Todas as variáveis são recebidas como ponteiros

    FILE *arquivo = NULL;// Declara arquivo que vai receber os dados de leitura 
    char linha[1000], *pasta, *folder,aux[200];//declara strings para receber nome do diretório de arquivos

    pasta = (char *)malloc(200);//aloca memória dinâmica para a str pasta
    folder = (char *)malloc(200); //aloca memória dinamica para o nome do arquivo de entrada 
    printf("Leitura de dados da rede elétrica...\n"); // printa no console que a leitura de dados da rede começou 

    //Recebe o nome da pasta com os dados a serem lidos - arquivo config.txt
    FILE *config = NULL; // arquivo config.txt contendo o nome dos diretórios 
    config = fopen("config.txt", "r");
    if (config == NULL) // se não houver arquivo de dados
    {
        printf("Erro ao abrir arquivo config.txt !!!\n");
        exit(1);
    }
    fgets(linha, 1000, config);//lê a primeira linha contendo o endereço do arquivo principal
    folder = getfield(linha, 1);//Salva na string folder
    fgets(linha, 1000, config);//lê a segunda linha do arquivo com a sub pasta dos dados da rede
    pasta = getfield(linha, 1);//Salva na string pasta
    printf("Main directory: \n %s \n", folder);// printa o nome dos arquivos no console
    printf("Data sub-folder: \n %s \n", pasta);// printa o nome dos arquivos no console
    fclose(config);//fecha o arquivo


     strcpy(aux,pasta);// string aux recebe o string em pasta
    //strcpy(aux2, aux);// utiliza o aux 2 pra concatenar nomes convenientemente no aux 
    // Leitura dos dados de barras

    arquivo = fopen(strcat(aux, "DBAR.csv"), "r"); //abre o arquivo /pasta/DBAR.csv

    if (arquivo != NULL)
    {
        leituraDBAR(arquivo, barra, numeroBarras, numeroAlimentadores);//lê o DBAR
        fclose(arquivo);//
    }
    else
    {
        printf("Erro ao abrir arquivo DBAR.csv !!!\n");
        exit(1);
    }

    printf("DBAR ok\n");
    strcpy(aux, pasta);
    arquivo = fopen(strcat(aux, "DSHNT.csv"), "r"); //Le somente se existir o arquivo
    if (arquivo != NULL)
    {
        leituraDSHNT(arquivo, barra, numeroBarras);
        //printf("DSHNT ok\n");
        fclose(arquivo);
    }
    
    strcpy(aux, pasta);
    arquivo = fopen(strcat(aux, "DGD.csv"), "r"); //Le somente se existir o arquivo
    if (arquivo != NULL)
    {
        leituraDGD(arquivo, barra, numeroBarras);
        //printf("DGD ok\n");
        fclose(arquivo);
    }

    // Leitura dos dados de ramos

    strcpy(aux, pasta);
    arquivo = fopen(strcat(aux, "DLIN.csv"), "r"); //Le somente se existir o arquivo
    if (arquivo != NULL)
    {
        leituraDLIN(arquivo, ramo, numeroRamos, barra, numeroBarras);
        //printf("DLIN ok\n");
        fclose(arquivo);
    }
    //Leitura dos dados dos trafos, que entram nos ramos
    strcpy(aux, pasta);
    arquivo = fopen(strcat(aux, "DTRF.csv"), "r"); //Le somente se existir o arquivo
    if (arquivo != NULL)
    {
        leituraDTRF(arquivo, ramo, numeroRamos, barra, numeroBarras);
        //printf("DTRF ok\n");
        fclose(arquivo);
    }
    //Leitura dos dados dos reguladores, que entram nos ramos
    strcpy(aux, pasta);
    arquivo = fopen(strcat(aux, "DREG.csv"), "r"); //Le somente se existir o arquivo
    if (arquivo != NULL)
    {
        leituraDREG(arquivo, ramo, numeroRamos, barra, numeroBarras);
        //printf("DREG ok\n");
        fclose(arquivo);
    }
    // Dados das chaves 
    strcpy(aux, pasta);
    arquivo = fopen(strcat(aux, "DSWTC.csv"), "r"); //Le somente se existir o arquivo
    if (arquivo != NULL)
    {
        leituraDSWTC(arquivo, ramo, numeroRamos, barra, numeroBarras);
        //printf("DSWTC ok\n");
        fclose(arquivo);
    }

    //Leitura do arquivo de tensões iniciais
    strcpy(aux, pasta);
    arquivo = fopen(strcat(aux, "Vinicial.csv"), "r"); //Le somente se existir o arquivo
    if (arquivo != NULL)
    {
        leituraVinicial(arquivo, barra, numeroBarras);
        //printf("Vinicial ok\n");
        fclose(arquivo);
    }

    folder = getfield(pasta, 1);
    return (folder);
}

//------------------------------------------------------------------------------
// Leitura de dados de BARRA
void leituraDBAR(FILE *arquivo, DBAR **barras, long int *numeroBarras, long int *numeroAlimentadores)
{
    char blocoLeitura[2000];     /* Variável para realizar a leitura do bloco de caracteres do arquivo. */
    char *dados;                 /* Variável do tipo ponteiro para char, utilizada para alterar o ponteiro da string lida do arquivo de forma a realizar o loop no sscanf. */
    int contador = 0, i, aux, k; /* Variáveis contadores para percorrer o arquivo e a string de leitura. */
    int carac, numLinhas = 0;    /* Variável com o número de linhas do arquivo a serem lidas. */
    double PA, PB, PC, QA, QB, QC;
    dados = (char *)malloc(100);
    //Aloca na memória espaço para as barras
    while ((carac = fgetc(arquivo)) != EOF)
    {
        if (carac == '\n')
            numLinhas++;
    }//conta o número de linhas no arquivo
    rewind(arquivo);// retorna ao topo do arquivo
    if (((*barras) = (DBAR *)malloc((numLinhas + 1) * sizeof(DBAR))) == NULL)//aloca a memória para a estrutura de dados do DBAR
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para as barras !!!!");
        exit(1);
    }
    //printf("numero linhas: %d\n", numLinhas);
    // Le o arquivo de curva de cargas até o fim
    while ((fgets(blocoLeitura, 2000, arquivo)) != NULL)//lê linha por linha do arquivo DBAR
    {

        dados = blocoLeitura;
        //printf("%s\n", blocoLeitura);

        //Verifica se a barra já foi criada
        aux = -1;
        for (i = 0; i < contador; i++)
        {
            if ((*barras)[i].ID == atoi(getfield(dados, 1)))//lê campo 1 da string e vê se ela já exite na esturtura de o DBAR
            {
                aux = i;
            }
        }

        if (aux == -1)
        { //Criando novo DBAR
            (*barras)[contador].ID = atoi(getfield(dados, 1));//salva ID da barra, numero conforme o dbar
            (*barras)[contador].i = contador;//salva o número da barra de acordo com a ordem que ela aparece no dbar
            (*barras)[contador].ligacao = atoi(getfield(dados, 2)); //Tipo da ligacao carga (1 = YN / 2 = Delta / 3 = Y) 
            (*barras)[contador].fases = atoi(getfield(dados, 3)); //Número de fases (1=A; 2=B; 3=C; 4=AB; 5=CA; 6=BC; 7=ABC);
            (*barras)[contador].Vbase = atof(getfield(dados, 4)) / pow(3, 0.5);//Tensão nominal da barra em V dividido por raiz de 3 para encontrar de fase
            (*barras)[contador].tipo = 0;//tipo recebe 0

            //(*barras)[contador].loads = (DLOAD *)malloc( 1 * sizeof(DLOAD));
            (*barras)[contador].nloads = 0;//numero de cargas na barra
            //(*barras)[contador].shunts = (DSHNT *)malloc( 1 * sizeof(DSHNT));
            (*barras)[contador].nshunts = 0;//numero de shunts
            //(*barras)[contador].gds = (DGD *)malloc( 1 * sizeof(DGD));
            (*barras)[contador].ngds = 0;//numero de geração distribuida

            (*barras)[contador].Vinicial[0] = 1 * (cos(0) + I * sin(0));//flat start sem compensar a defazagem do trafo
            (*barras)[contador].Vinicial[1] = 1 * (cos(-120 * PI / 180) + I * sin(-120 * PI / 180));//flat start sem compensar a defazagem do trafo
            (*barras)[contador].Vinicial[2] = 1 * (cos(120 * PI / 180) + I * sin(120 * PI / 180));//flat start sem compensar a defazagem do trafo

            //Leitura das Cargas
            PA = atof(getfield(dados, 5));// potencia ativa da carga trifásica em kW
            PB = atof(getfield(dados, 6));// potencia ativa da carga trifásica em kW
            PC = atof(getfield(dados, 7));// potencia ativa da carga trifásica em kW
            QA = atof(getfield(dados, 8));// potência reativa da carga trifásica em kVAr
            QB = atof(getfield(dados, 9));// potência reativa da carga trifásica em kVAr
            QC = atof(getfield(dados, 10));// potência reativa da carga trifásica em kVAr

            if ((PA != 0) || (PB != 0) || (PC != 0) || (QA != 0) || (QB != 0) || (QC != 0))
            {
                // caso algum dos campos de carga do dbar seja diferente de 0 ele insere isto DBAR
                //Inseri um valor de load
                (*barras)[contador].nloads++;
                /*if (((*barras)[contador].loads = (DLOAD *)realloc((*barras)[contador].loads, ((*barras)[contador].nloads +1) * sizeof(DLOAD)))==NULL)
                {
                    printf("Erro -- Nao foi possivel alocar espaco de memoria para cargas !!!!");
                    exit(1); 
                }*/
                k = (*barras)[contador].nloads - 1;

                (*barras)[contador].loads[k].ID = (*barras)[contador].ID;
                (*barras)[contador].loads[k].Vbase = (*barras)[contador].Vbase;
                (*barras)[contador].loads[k].fases = (*barras)[contador].fases;
                (*barras)[contador].loads[k].lig = (*barras)[contador].ligacao;
                (*barras)[contador].loads[k].ZIP = atof(getfield(dados, 11));

                (*barras)[contador].loads[k].Pnom[0] = PA;
                (*barras)[contador].loads[k].Pnom[1] = PB;
                (*barras)[contador].loads[k].Pnom[2] = PC;
                (*barras)[contador].loads[k].Qnom[0] = QA;
                (*barras)[contador].loads[k].Qnom[1] = QB;
                (*barras)[contador].loads[k].Qnom[2] = QC;
            }

            //Leitura da Barra de Referência definida por ser aquela que tem valors de tensão e angulo 
            if (getfield(dados, 13) != NULL)
            {
                (*barras)[contador].tipo = 2;//tipo da barra de referencia
                numeroAlimentadores[0]++;//incrementa o número de alimentadores

                double VA = atof(getfield(dados, 12));
                double VB = atof(getfield(dados, 13));
                double VC = atof(getfield(dados, 14));
                double TA = atof(getfield(dados, 15));
                double TB = atof(getfield(dados, 16));
                double TC = atof(getfield(dados, 17));

                __real__(*barras)[contador].Vref[0] = VA * cos(TA * M_PI / 180);//atribui os vaores reais e imaginarios das tensões inciaias obs o operador __real__/__imag__ permite voce acessar as partes reais e imaginárias de um numero complexo, para escrever ou ler, mas só funcina em compiladores GCC
                __imag__(*barras)[contador].Vref[0] = VA * sin(TA * M_PI / 180);//atribui os vaores reais e imaginarios das tensões inciaias obs o operador __real__/__imag__ permite voce acessar as partes reais e imaginárias de um numero complexo, para escrever ou ler, mas só funcina em compiladores GCC
                __real__(*barras)[contador].Vref[1] = VB * cos(TB * M_PI / 180);//atribui os vaores reais e imaginarios das tensões inciaias obs o operador __real__/__imag__ permite voce acessar as partes reais e imaginárias de um numero complexo, para escrever ou ler, mas só funcina em compiladores GCC
                __imag__(*barras)[contador].Vref[1] = VB * sin(TB * M_PI / 180);//atribui os vaores reais e imaginarios das tensões inciaias obs o operador __real__/__imag__ permite voce acessar as partes reais e imaginárias de um numero complexo, para escrever ou ler, mas só funcina em compiladores GCC
                __real__(*barras)[contador].Vref[2] = VC * cos(TC * M_PI / 180);//atribui os vaores reais e imaginarios das tensões inciaias obs o operador __real__/__imag__ permite voce acessar as partes reais e imaginárias de um numero complexo, para escrever ou ler, mas só funcina em compiladores GCC
                __imag__(*barras)[contador].Vref[2] = VC * sin(TC * M_PI / 180);//atribui os vaores reais e imaginarios das tensões inciaias obs o operador __real__/__imag__ permite voce acessar as partes reais e imaginárias de um numero complexo, para escrever ou ler, mas só funcina em compiladores GCC

                (*barras)[contador].Vinicial[0] = VA * (cos(TA * M_PI / 180) + I * sin(TA * M_PI / 180));// atribui o numero complexo a tensão inicial da barra de referencia
                (*barras)[contador].Vinicial[1] = VB * (cos(TB * M_PI / 180) + I * sin(TB * M_PI / 180));// atribui o numero complexo a tensão inicial da barra de referencia
                (*barras)[contador].Vinicial[2] = VC * (cos(TC * M_PI / 180) + I * sin(TC * M_PI / 180));// atribui o numero complexo a tensão inicial da barra de referencia
            }
            contador++;
        }
        else
        { // Inserindo nova carga em DBAR já existente

            // caso a barra esteja repitida ele insere ela como uma carga nova
            //Leitura das Cargas
            PA = atof(getfield(dados, 5));
            PB = atof(getfield(dados, 6));
            PC = atof(getfield(dados, 7));
            QA = atof(getfield(dados, 8));
            QB = atof(getfield(dados, 9));
            QC = atof(getfield(dados, 10));

            if ((PA != 0) || (PB != 0) || (PC != 0) || (QA != 0) || (QB != 0) || (QC != 0))
            {
                //Inseri um valor de load
                (*barras)[aux].nloads++;
                /*if (((*barras)[aux].loads = (DLOAD *)realloc((*barras)[aux].loads, (*barras)[aux].nloads * sizeof(DLOAD)))==NULL)
                {
                    printf("Erro -- Nao foi possivel alocar espaco de memoria para cargas !!!!");
                    exit(1); 
                }*/
                k = (*barras)[aux].nloads - 1;
                (*barras)[aux].loads[k].ID = (*barras)[aux].ID;
                (*barras)[aux].loads[k].lig = atoi(getfield(dados, 2));
                (*barras)[aux].loads[k].fases = atoi(getfield(dados, 3));
                (*barras)[aux].loads[k].Vbase = atof(getfield(dados, 4));
                (*barras)[aux].loads[k].ZIP = atof(getfield(dados, 11));

                (*barras)[aux].loads[k].Pnom[0] = PA;
                (*barras)[aux].loads[k].Pnom[1] = PB;
                (*barras)[aux].loads[k].Pnom[2] = PC;
                (*barras)[aux].loads[k].Qnom[0] = QA;
                (*barras)[aux].loads[k].Qnom[1] = QB;
                (*barras)[aux].loads[k].Qnom[2] = QC;
            }
        }
        //printf("contador: %d\n", contador);
    }
    //printf("contador: %d, ", contador);
    numeroBarras[0] = contador;
}

//------------------------------------------------------------------------------
// Leitura de dados de SHUNTS
void leituraDSHNT(FILE *arquivo, DBAR **barras, long int *numeroBarras)
{
    char blocoLeitura[2000]; /* Variável para realizar a leitura do bloco de caracteres do arquivo. */
    char *dados;             /* Variável do tipo ponteiro para char, utilizada para alterar o ponteiro da string lida do arquivo de forma a realizar o loop no sscanf. */
    int i, aux, k;           /* Variáveis contadores para percorrer o arquivo e a string de leitura. */
    dados = (char *)malloc(100);
    // Le o arquivo de curva de cargas até o fim
    while ((fgets(blocoLeitura, 2000, arquivo)) != NULL)
    {
        dados = blocoLeitura;

        //Verifica se a barra existe
        aux = -1;
        for (i = 0; i < numeroBarras[0]; i++)
        {
            if ((*barras)[i].ID == atoi(getfield(dados, 1)))
            {
                aux = i;
            }
        }
        if (aux != -1)//se existe ele insere o shunt na barra
        {
            (*barras)[aux].nshunts++;//complementa o numero de shunts naquela barra
            if (((*barras)[aux].shunts = (DSHNT *)realloc((*barras)[aux].shunts, (*barras)[aux].nshunts * sizeof(DSHNT))) == NULL)//aloca ou realoca memória
            {
                printf("Erro -- Nao foi possivel alocar espaco de memoria para shunts !!!!");
                exit(1);
            }
            k = (*barras)[aux].nshunts - 1;//indica correntamente o contador do shunt

            (*barras)[aux].shunts[k].ID = (*barras)[aux].ID;//ID conforme o nome no DBAR
            (*barras)[aux].shunts[k].lig = atoi(getfield(dados, 2));//Tipo de ligação do banco shunt (1 = YN / 2 = Delta / 3 = Y);
            (*barras)[aux].shunts[k].fases = atoi(getfield(dados, 3));//Número de fases (1=A; 2=B; 3=C; 4=AB; 5=CA; 6=BC; 7=ABC);
            (*barras)[aux].shunts[k].Vbase = atof(getfield(dados, 4));//Tensão base de linha

            (*barras)[aux].shunts[k].Qnom[0] = atof(getfield(dados, 5));//Potencia nominal por fase
            (*barras)[aux].shunts[k].Qnom[1] = atof(getfield(dados, 6));//Potencia nominal por fase
            (*barras)[aux].shunts[k].Qnom[2] = atof(getfield(dados, 7));//Potencia nominal por fase
            (*barras)[aux].shunts[k].controle = atoi(getfield(dados, 8));//Tipo de controle (0= sem controle / 1=Controle de tensão da barra)

            if (getfield(dados, 10) != NULL)
            {
                (*barras)[aux].shunts[k].num = atoi(getfield(dados, 9));// Numero de bancos Shunt
                (*barras)[aux].shunts[k].DV = atof(getfield(dados, 10));//Intervalo de tensão de controle em PU
                (*barras)[aux].shunts[k].Vset[0] = atof(getfield(dados, 11));//Set point da tensão (em pu)
                (*barras)[aux].shunts[k].Vset[1] = atof(getfield(dados, 12));//Set point da tensão (em pu)
                (*barras)[aux].shunts[k].Vset[2] = atof(getfield(dados, 13));//Set point da tensão (em pu)
            }
        }
    }
}

//------------------------------------------------------------------------------
// Leitura de dados de GDs
void leituraDGD(FILE *arquivo, DBAR **barras, long int *numeroBarras)
{
    //as GDS ebtram na estrutura DBAR
    char blocoLeitura[2000]; /* Variável para realizar a leitura do bloco de caracteres do arquivo. */
    char *dados;             /* Variável do tipo ponteiro para char, utilizada para alterar o ponteiro da string lida do arquivo de forma a realizar o loop no sscanf. */
    int i, aux, k;           /* Variáveis contadores para percorrer o arquivo e a string de leitura. */
    dados = (char *)malloc(100);// aloca memória para a string q vai receber cada linha do arquivo
    // Le o arquivo de curva de cargas até o fim
    while ((fgets(blocoLeitura, 2000, arquivo)) != NULL)//lê linha por linha do arquivo DGD.csv
    {
        dados = blocoLeitura;//pega linha lida

        //Verifica se a barra já foi existe
        aux = -1;
        for (i = 0; i < numeroBarras[0]; i++)
        {
            if ((*barras)[i].ID == atoi(getfield(dados, 1)))
            {
                aux = i;
            }
        }
        if (aux != -1)
        {
            // se a barra existe ele insere o GD
            (*barras)[aux].ngds++;//incrementa o contador de GDS
            if (((*barras)[aux].gds = (DGD *)realloc((*barras)[aux].gds, (*barras)[aux].ngds * sizeof(DGD))) == NULL)//Aloca memória para a estrutura GD na barra
            {
                printf("Erro -- Nao foi possivel alocar espaco de memoria para gds !!!!");
                exit(1);
            }
            k = (*barras)[aux].ngds - 1;

            (*barras)[aux].gds[k].ID = (*barras)[aux].ID;//salva o id da barra do GD igual no DBAR
            (*barras)[aux].gds[k].lig = atoi(getfield(dados, 2));//salva o tipo de ligação
            (*barras)[aux].gds[k].fases = atoi(getfield(dados, 3));// salva as fases conectadas
            (*barras)[aux].gds[k].Vbase = atof(getfield(dados, 4));//tensão base, diferente do DBAR n é dividido por raiz de 3
            (*barras)[aux].gds[k].Snominal = atof(getfield(dados, 5));//potência nominal

            (*barras)[aux].gds[k].Pnom[0] = atof(getfield(dados, 6));//Potência ativa nominal por fase 0 == a,1==b ,2 == c
            (*barras)[aux].gds[k].Pnom[1] = atof(getfield(dados, 7));//Potência ativa nominal por fase 0 == a,1==b ,2 == c
            (*barras)[aux].gds[k].Pnom[2] = atof(getfield(dados, 8));//Potência ativa nominal por fase 0 == a,1==b ,2 == c
            (*barras)[aux].gds[k].Qnom[0] = atof(getfield(dados, 9));//Potência reativa nominal por fase  0 == a,1==b ,2 == c
            (*barras)[aux].gds[k].Qnom[1] = atof(getfield(dados, 10));//Potência reativa nominal por fase 0 == a,1==b ,2 == c
            (*barras)[aux].gds[k].Qnom[2] = atof(getfield(dados, 11));//Potência reativa nominal por fase 0 == a,1==b ,2 == c
            (*barras)[aux].gds[k].controle = atoi(getfield(dados, 12));// Controle da excitatriz do gerador (0 - Máquina como PQ constante / 1 = Máquina como PV constante / 2 = Corrente constante);

            if (getfield(dados, 14) != NULL)//se existe dados após o campo 14
            {
                (*barras)[aux].gds[k].Qmin = atof(getfield(dados, 13));//Limites mínimo e máximo de geração de potência reativa;
                (*barras)[aux].gds[k].Qmax = atof(getfield(dados, 14));//Limites mínimo e máximo de geração de potência reativa;
                (*barras)[aux].gds[k].Vset[0] = atof(getfield(dados, 15));//Setpoint de tensão para máquina operando como PV constante
                (*barras)[aux].gds[k].Vset[1] = atof(getfield(dados, 16));//Setpoint de tensão para máquina operando como PV constante
                (*barras)[aux].gds[k].Vset[2] = atof(getfield(dados, 17));//Setpoint de tensão para máquina operando como PV constante
                (*barras)[aux].gds[k].controlePV = atoi(getfield(dados, 18));//0 Sequência positiva / 2 por fase
            }
        }
    }
}

//------------------------------------------------------------------------------
// Leitura de dados de LINHAS
void leituraDLIN(FILE *arquivo, DRAM **ramos, long int *numeroRamos, DBAR **barras, long int *numeroBarras)
{
    char blocoLeitura[2000];       /* Variável para realizar a leitura do bloco de caracteres do arquivo. */
    char *dados;                   /* Variável do tipo ponteiro para char, utilizada para alterar o ponteiro da string lida do arquivo de forma a realizar o loop no sscanf. */
    long int contador = 0, i, aux; /* Variáveis contadores para percorrer o arquivo e a string de leitura. */
    int carac, numLinhas = 0;      /* Variável com o número de linhas do arquivo a serem lidas. */
    dados = (char *)malloc(100); //string q vai receber linha por linha do arquivo de dados

    //Aloca na memória espaço para as linhas
    while ((carac = fgetc(arquivo)) != EOF)
    {
        if (carac == '\n')
            numLinhas++;
    }
    rewind(arquivo);
    if (((*ramos) = (DRAM *)malloc((numLinhas) * sizeof(DRAM))) == NULL)
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para as linhas !!!!");
        exit(1);
    }

    // Le o arquivo de curva de cargas até o fim
    while ((fgets(blocoLeitura, 2000, arquivo)) != NULL)
    {
        dados = blocoLeitura;

        (*ramos)[contador].DE = atoi(getfield(dados, 1));//nós terminais de cada ramo   
        (*ramos)[contador].PARA = atoi(getfield(dados, 2));//nós terminais de cada ramo
        (*ramos)[contador].tipo = 0;//Tipo fixo em 0
        (*ramos)[contador].estado = 1;//estado fixo em 1
        (*ramos)[contador].fases = atoi(getfield(dados, 3));//Número de fases (1=A; 2=B; 3=C; 4=AB; 5=CA; 6=BC; 7=ABC); 

        for (i = 0; i < numeroBarras[0]; i++)//atribui o nome de DE para conforme a ordem do DBAR
        {
            if ((*ramos)[contador].DE == (*barras)[i].ID)
            {
                (*ramos)[contador].k = (*barras)[i].i;
            }
            if ((*ramos)[contador].PARA == (*barras)[i].ID)
            {
                (*ramos)[contador].m = (*barras)[i].i;
            }
        }

        //Preenche o (*ramos)[contador].linha
        (*ramos)[contador].linha.fases = atoi(getfield(dados, 3));
        (*ramos)[contador].linha.comprimento = atof(getfield(dados, 4));
        double comp = atof(getfield(dados, 4));
        __real__(*ramos)[contador].linha.Zaa = comp * atof(getfield(dados, 5));//lê parte real e imaginária do CSV 
        __imag__(*ramos)[contador].linha.Zaa = comp * atof(getfield(dados, 6));//lê parte real e imaginária do CSV 
        __real__(*ramos)[contador].linha.Zab = comp * atof(getfield(dados, 7));//lê parte real e imaginária do CSV 
        __imag__(*ramos)[contador].linha.Zab = comp * atof(getfield(dados, 8));//lê parte real e imaginária do CSV 
        __real__(*ramos)[contador].linha.Zac = comp * atof(getfield(dados, 9));//lê parte real e imaginária do CSV 
        __imag__(*ramos)[contador].linha.Zac = comp * atof(getfield(dados, 10));//lê parte real e imaginária do csv 
        __real__(*ramos)[contador].linha.Zbb = comp * atof(getfield(dados, 11));//lê parte real e imaginária do csv 
        __imag__(*ramos)[contador].linha.Zbb = comp * atof(getfield(dados, 12));//lê parte real e imaginária do csv 
        __real__(*ramos)[contador].linha.Zbc = comp * atof(getfield(dados, 13));//lê parte real e imaginária do csv 
        __imag__(*ramos)[contador].linha.Zbc = comp * atof(getfield(dados, 14));//lê parte real e imaginária do csv 
        __real__(*ramos)[contador].linha.Zcc = comp * atof(getfield(dados, 15));//lê parte real e imaginária do csv 
        __imag__(*ramos)[contador].linha.Zcc = comp * atof(getfield(dados, 16));//lê parte real e imaginária do csv 
        (*ramos)[contador].linha.Baa = comp * atof(getfield(dados, 17)) /1000000;//lê os dados do shunt e divide por 1000000????
        (*ramos)[contador].linha.Bab = comp * atof(getfield(dados, 18)) /1000000;//lê os dados do shunt e divide por 1000000????
        (*ramos)[contador].linha.Bac = comp * atof(getfield(dados, 19)) /1000000;//lê os dados do shunt e divide por 1000000????
        (*ramos)[contador].linha.Bbb = comp * atof(getfield(dados, 20)) /1000000;//lê os dados do shunt e divide por 1000000????
        (*ramos)[contador].linha.Bbc = comp * atof(getfield(dados, 21)) /1000000;//lê os dados do shunt e divide por 1000000????
        (*ramos)[contador].linha.Bcc = comp * atof(getfield(dados, 22)) /1000000;//lê os dados do shunt e divide por 1000000????

        contador++;//inclementa contador de linhas
    }
    numeroRamos[0] = contador;
}

//------------------------------------------------------------------------------
// Leitura de dados de TRAFOS
void leituraDTRF(FILE *arquivo, DRAM **ramos, long int *numeroRamos, DBAR **barras, long int *numeroBarras)
{
    char blocoLeitura[2000];  /* Variável para realizar a leitura do bloco de caracteres do arquivo. */
    char *dados;              /* Variável do tipo ponteiro para char, utilizada para alterar o ponteiro da string lida do arquivo de forma a realizar o loop no sscanf. */
    int contador = 0, i, aux; /* Variáveis contadores para percorrer o arquivo e a string de leitura. */
    int carac, numLinhas = 0; /* Variável com o número de linhas do arquivo a serem lidas. */
    dados = (char *)malloc(100); //variavel com cada linha

    //Aloca na memória espaço para os trafos
    while ((carac = fgetc(arquivo)) != EOF)
    {
        if (carac == '\n')
            numLinhas++;// conta numero de linhas do arquivo
    }
    rewind(arquivo);//reinicia a leitura do arquivo

    if (((*ramos) = (DRAM *)realloc((*ramos), (numeroRamos[0] + numLinhas) * sizeof(DRAM))) == NULL)//aloca a memória para os ramos
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para os trafos !!!!");
        exit(1);
    }

    contador = numeroRamos[0];//recebe o número atual de ramos, contando os arquivos anteriores lidos
    // Le o arquivo até o fim
    while ((fgets(blocoLeitura, 2000, arquivo)) != NULL)
    {
        dados = blocoLeitura;//lê cada linha

        (*ramos)[contador].DE = atoi(getfield(dados, 1));//lê o de PARA conforme o nome da barra no DBAR
        (*ramos)[contador].PARA = atoi(getfield(dados, 2));//lê o de PARA conforme o nome da barra no DBAR
        (*ramos)[contador].tipo = 1;//Tipo 1 == trafo
        (*ramos)[contador].estado = 1;// estado ==1
        (*ramos)[contador].fases = atoi(getfield(dados, 3));//Número de fases (1=A; 2=B; 3=C; 4=AB; 5=CA; 6=BC; 7=ABC); 

        for (i = 0; i < numeroBarras[0]; i++)//atribui o k m (numeração do DBAR de acordo com a ordem no contador)
        {
            if ((*ramos)[contador].DE == (*barras)[i].ID)
            {
                (*ramos)[contador].k = (*barras)[i].i;
            }
            if ((*ramos)[contador].PARA == (*barras)[i].ID)
            {
                (*ramos)[contador].m = (*barras)[i].i;
            }
        }

        //Preencher o (*ramos)[contador].trafo
        (*ramos)[contador].trafo.fases = atoi(getfield(dados, 3));//Número de fases (1=A; 2=B; 3=C; 4=AB; 5=CA; 6=BC; 7=ABC); 
        (*ramos)[contador].trafo.Vpri = atof(getfield(dados, 4));//Tensão do primario em V
        (*ramos)[contador].trafo.Vsec = atof(getfield(dados, 5));//Tensão do sencundario em V
        (*ramos)[contador].trafo.Snominal = atof(getfield(dados, 6));// Potencia aparente nominak
        (*ramos)[contador].trafo.R = atof(getfield(dados, 7)) * pow((*ramos)[contador].trafo.Vsec, 2) / ((*ramos)[contador].trafo.Snominal * 1000) / 3; //calcula o R e o x já transformando pra pu R*(Vpri^2/(S*1000))
        (*ramos)[contador].trafo.X = atof(getfield(dados, 8)) * pow((*ramos)[contador].trafo.Vsec, 2) / ((*ramos)[contador].trafo.Snominal * 1000) / 3;//calcula o R e o x já transformando pra pu x*(Vpri^2/(S*1000))
        (*ramos)[contador].trafo.lig_pri = atoi(getfield(dados, 9));//Tipo de ligação do primário (1 = YN / 2 = Delta / 3 = Y / 4 - Estrela aberto / 5 - Delta aberto);
        (*ramos)[contador].trafo.lig_sec = atoi(getfield(dados, 10));//Tipo de ligação do Secundário (1 = YN / 2 = Delta / 3 = Y / 4 - Estrela aberto / 5 - Delta aberto);
        (*ramos)[contador].trafo.defasamento = atoi(getfield(dados, 11));//Defasamento angular (Ex: DYn0 ou DYn1);
        (*ramos)[contador].trafo.tap_pri = atof(getfield(dados, 12));//Tap do primário fora do nominal (Ex: 1.025);
        (*ramos)[contador].trafo.tap_sec = atof(getfield(dados, 13));//Tap do Secundário fora do nominal (Ex: 1.025);

        contador++;//complementa o contador de ramos
    }
    numeroRamos[0] = numeroRamos[0] + numLinhas;//complementa o número de ramos
}

//------------------------------------------------------------------------------
// Leitura de dados de REGULADORES
void leituraDREG(FILE *arquivo, DRAM **ramos, long int *numeroRamos, DBAR **barras, long int *numeroBarras)
{
    char blocoLeitura[2000];  /* Variável para realizar a leitura do bloco de caracteres do arquivo. */
    char *dados;              /* Variável do tipo ponteiro para char, utilizada para alterar o ponteiro da string lida do arquivo de forma a realizar o loop no sscanf. */
    int contador = 0, i, aux; /* Variáveis contadores para percorrer o arquivo e a string de leitura. */
    int carac, numLinhas = 0; /* Variável com o número de linhas do arquivo a serem lidas. */
    dados = (char *)malloc(100);

    //Aloca na memória espaço para os reguladores
    while ((carac = fgetc(arquivo)) != EOF)//lê o arquivo char a char e conta quantas linhas tem
    {
        if (carac == '\n')
            numLinhas++;
    }
    rewind(arquivo);//retorna ao topo do arquivo
    if (((*ramos) = (DRAM *)realloc((*ramos), (numeroRamos[0] + numLinhas) * sizeof(DRAM))) == NULL)//aumenta o espaço alocado na memória
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para os reguladores de tensão !!!!");
        exit(1);
    }
    contador = numeroRamos[0];
    // Le o arquivo de curva de cargas até o fim
    while ((fgets(blocoLeitura, 2000, arquivo)) != NULL)
    {
        dados = blocoLeitura;

        (*ramos)[contador].DE = atoi(getfield(dados, 1));//Preenche o DE conforme o nome dado no DREG
        (*ramos)[contador].PARA = atoi(getfield(dados, 2));//Preenche o PARA conforme o nome dado no DREG
        (*ramos)[contador].tipo = 2;//Tipo 2 == Regulador
        (*ramos)[contador].estado = 1;//Estado == 1
        (*ramos)[contador].fases = atoi(getfield(dados, 3));//Número de fases (1=A; 2=B; 3=C; 4=AB; 5=CA; 6=BC; 7=ABC); 

        for (i = 0; i < numeroBarras[0]; i++)//Atribui o km conforme a numeração da barra do DBAR, segundo a ordem
        {
            if ((*ramos)[contador].DE == (*barras)[i].ID)
            {
                (*ramos)[contador].k = (*barras)[i].i;
            }
            if ((*ramos)[contador].PARA == (*barras)[i].ID)
            {
                (*ramos)[contador].m = (*barras)[i].i;
            }
        }

        //Preencher o (*ramos)[contador].regulador
        (*ramos)[contador].regulador.fases = atoi(getfield(dados, 3));//Número de fases (1=A; 2=B; 3=C; 4=AB; 5=CA; 6=BC; 7=ABC); 
        (*ramos)[contador].regulador.Vnom = atof(getfield(dados, 4));//Tensão nominal
        (*ramos)[contador].regulador.regulacao = atof(getfield(dados, 5));// Regulação de tensão Ex: 10% entrar como 0.1
        (*ramos)[contador].regulador.ntaps = atoi(getfield(dados, 6));//Número de taps (Ex: 16 para cima e 16 para baixo);
        (*ramos)[contador].regulador.Snominal = atof(getfield(dados, 7));//Potência nominal do regulador de tensão;
        (*ramos)[contador].regulador.R = atof(getfield(dados, 8)) * pow((*ramos)[contador].regulador.Vnom, 2) / ((*ramos)[contador].regulador.Snominal * 1000) / 3;//Resistência série do regulador o programa já transorma em pu
        (*ramos)[contador].regulador.X = atof(getfield(dados, 9)) * pow((*ramos)[contador].regulador.Vnom, 2) / ((*ramos)[contador].regulador.Snominal * 1000) / 3;//Reatância série do regulador o programa já transorma em pu
        (*ramos)[contador].regulador.lig = atoi(getfield(dados, 10));// Tipo de ligação do regulador (1 = YN / 2 = Delta / 3 = Y / 4 - Estrela aberto / 5 - Delta aberto);
        (*ramos)[contador].regulador.TP = atof(getfield(dados, 11));// Relação do TP do controlador do regulador;
        (*ramos)[contador].regulador.TC = atof(getfield(dados, 12));//Relação do TC do controlador do regulador;
        (*ramos)[contador].regulador.deltaV = atof(getfield(dados, 13));//Intervalo de tensão para atuação em V secundários (Ex: 2 V significa +1 ou -1V em relação ao ajuste V1, V2 e V3 para atuação);
        (*ramos)[contador].regulador.R1 = atof(getfield(dados, 14));//Parâmetros do controlador LDC por fase em ohms;
        (*ramos)[contador].regulador.X1 = atof(getfield(dados, 15));//Parâmetros do controlador LDC por fase em ohms;
        (*ramos)[contador].regulador.R2 = atof(getfield(dados, 16));//Parâmetros do controlador LDC por fase em ohms;
        (*ramos)[contador].regulador.X2 = atof(getfield(dados, 17));//Parâmetros do controlador LDC por fase em ohms;
        (*ramos)[contador].regulador.R3 = atof(getfield(dados, 18));//Parâmetros do controlador LDC por fase em ohms;
        (*ramos)[contador].regulador.X3 = atof(getfield(dados, 19));//Parâmetros do controlador LDC por fase em ohms;
        (*ramos)[contador].regulador.V1 = atof(getfield(dados, 20));//Ajuste de tensão do controlador LDC em V secundários;
        (*ramos)[contador].regulador.V2 = atof(getfield(dados, 21));//Ajuste de tensão do controlador LDC em V secundários;
        (*ramos)[contador].regulador.V3 = atof(getfield(dados, 22));//Ajuste de tensão do controlador LDC em V secundários;
        (*ramos)[contador].regulador.controle = atoi(getfield(dados, 23));//Tipo de Controle (0 = Somente Forward sem Restrição / 1 = Locked Forward / 2 = Locked Reverse / 3 = Bidirectional / 4 = Idle / 5 = Neutral Reverse / 6 = Cogenartion / Ver manual Siemens MJ4A)
        (*ramos)[contador].regulador.tap[0] = atof(getfield(dados, 24));//Parâmetro de tap inicial (opcional - serve como tap do controle reverse locked)
        (*ramos)[contador].regulador.tap[1] = atof(getfield(dados, 25));//Parâmetro de tap inicial (opcional - serve como tap do controle reverse locked)
        (*ramos)[contador].regulador.tap[2] = atof(getfield(dados, 26));//Parâmetro de tap inicial (opcional - serve como tap do controle reverse locked)

        //Se tiver parametros de controle reverso
        /*
        (*ramos)[contador].regulador.deltaVr = atof(getfield(dados,27));
        (*ramos)[contador].regulador.R1r = atof(getfield(dados,28));
        (*ramos)[contador].regulador.X1r = atof(getfield(dados,29));
        (*ramos)[contador].regulador.R2r = atof(getfield(dados,30));
        (*ramos)[contador].regulador.X2r = atof(getfield(dados,31));
        (*ramos)[contador].regulador.R3r = atof(getfield(dados,32));
        (*ramos)[contador].regulador.X3r = atof(getfield(dados,33));
        (*ramos)[contador].regulador.V1r = atof(getfield(dados,34));
        (*ramos)[contador].regulador.V2r = atof(getfield(dados,35));
        (*ramos)[contador].regulador.V3r = atof(getfield(dados,36));
        */
        contador++;//incrementa o contador
    }
    numeroRamos[0] = numeroRamos[0] + numLinhas;//incrementa o número de ramos com os taps
}

//------------------------------------------------------------------------------
// Leitura de dados de CHAVES
void leituraDSWTC(FILE *arquivo, DRAM **ramos, long int *numeroRamos, DBAR **barras, long int *numeroBarras)
{
    char blocoLeitura[2000];  /* Variável para realizar a leitura do bloco de caracteres do arquivo. */
    char *dados;              /* Variável do tipo ponteiro para char, utilizada para alterar o ponteiro da string lida do arquivo de forma a realizar o loop no sscanf. */
    int contador = 0, i, aux; /* Variáveis contadores para percorrer o arquivo e a string de leitura. */
    int carac, numLinhas = 0; /* Variável com o número de linhas do arquivo a serem lidas. */
    dados = (char *)malloc(100);

    //Aloca na memória espaço para os trafos
    while ((carac = fgetc(arquivo)) != EOF)//Percorre caractere a caractere 
    {
        if (carac == '\n')
            numLinhas++;//conta o número de linhas
    }
    rewind(arquivo);//retorna ao início do arquivo

    if (((*ramos) = (DRAM *)realloc((*ramos), (numeroRamos[0] + numLinhas) * sizeof(DRAM))) == NULL) // aumenta o espaço alocado na estrutura dos ramos
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para as chaves !!!!");
        exit(1);
    }

    contador = numeroRamos[0];//contador aponta para a primeira posição vazia da estrutura dos ramos
    // Le o arquivo até o fim
    while ((fgets(blocoLeitura, 2000, arquivo)) != NULL)
    {
        dados = blocoLeitura;

        (*ramos)[contador].DE = atoi(getfield(dados, 1));//Recebe nome barra DE 
        (*ramos)[contador].PARA = atoi(getfield(dados, 2));//Recebe nome barra PARA
        (*ramos)[contador].tipo = 3;//TIPO 3 == switch 

        for (i = 0; i < numeroBarras[0]; i++)
        {
            if ((*ramos)[contador].DE == (*barras)[i].ID)//atribui o k e o m de acordo com a ordem do Dbar e a numeração de indices das barras
            {
                (*ramos)[contador].k = (*barras)[i].i;
            }
            if ((*ramos)[contador].PARA == (*barras)[i].ID)
            {
                (*ramos)[contador].m = (*barras)[i].i;
            }
        }
        (*ramos)[contador].fases = atoi(getfield(dados, 3));//Numero de fases
        (*ramos)[contador].estado = atoi(getfield(dados, 4));//Estado 0 aberto 1 fechado

        contador++;//incrementa o contador 
    }
    numeroRamos[0] = numeroRamos[0] + numLinhas;//atualiza o número de ramos
}

//------------------------------------------------------------------------------
// Leitura de valores de tensão de inicialização do estimador
void leituraVinicial(FILE *arquivo, DBAR **barras, long int *numeroBarras)
{
    char blocoLeitura[2000]; /* Variável para realizar a leitura do bloco de caracteres do arquivo. */
    char *dados;             /* Variável do tipo ponteiro para char, utilizada para alterar o ponteiro da string lida do arquivo de forma a realizar o loop no sscanf. */
    int i, aux, k;           /* Variáveis contadores para percorrer o arquivo e a string de leitura. */
    double Va, Vb, Vc, Ta, Tb, Tc;
    dados = (char *)malloc(100);
    // Le o arquivo de curva de cargas até o fim
    while ((fgets(blocoLeitura, 2000, arquivo)) != NULL)//lê arquivo linha por linha
    {
        dados = blocoLeitura;

        //Verifica se a barra existe
        aux = -1;
        for (i = 0; i < numeroBarras[0]; i++)
        {
            if ((*barras)[i].ID == atoi(getfield(dados, 1)))
            {
                aux = i;
            }
        }
        if (aux != -1)
        {
            Va = atof(getfield(dados, 2));//Recebe os Os moduos de tensão por fases
            Vb = atof(getfield(dados, 3));//Recebe os Os moduos de tensão por fases
            Vc = atof(getfield(dados, 4));//Recebe os Os moduos de tensão por fases
            Ta = atof(getfield(dados, 5)) * PI / 180;//recebe os angulos por fase em graus e transforma pra radianos
            Tb = atof(getfield(dados, 6)) * PI / 180;//recebe os angulos por fase em graus e transforma pra radianos
            Tc = atof(getfield(dados, 7)) * PI / 180;//recebe os angulos por fase em graus e transforma pra radianos

            if (Va != Va)/////????????????????????
                Va = 1;/////????????????????????
            if (Vb != Vb)//??????????????????????
                Vb = 1;//??????????????????????
            if (Vc != Vc)//????????????????????
                Vc = 1;//????????????????????
            if (Ta != Ta)/////????????????????
                Ta = 0 * PI / 180;/////????????????????
            if (Tb != Tb)//////////??????????????????
                Tb = -120 * PI / 180;//////////??????/////???????????????
            if (Tc != Tc)//////////??????????????????/////???
                Tc = 120 * PI / 180;/////???

            (*barras)[aux].Vinicial[0] = Va * (cos(Ta) + I * sin(Ta));//Transforma para complexo
            (*barras)[aux].Vinicial[1] = Vb * (cos(Tb) + I * sin(Tb));//Transforma para complexo
            (*barras)[aux].Vinicial[2] = Vc * (cos(Tc) + I * sin(Tc));//Transforma para complexo
        }
    }
}

//------------------------------------------------------------------------------
// Leitura de dados de MEDIDORES
long int **leituraMedidas(char *folder, char *file, DMED **medidas, DRAM *ramos, long int numeroRamos, DBAR *barras, long int numeroBarras, GRAFO *grafo, double Sbase)
{
    char blocoLeitura[2000];                             /* Variável para realizar a leitura do bloco de caracteres do arquivo. */
    char *dados;                                         /* Variável do tipo ponteiro para char, utilizada para alterar o ponteiro da string lida do arquivo de forma a realizar o loop no sscanf. */
    long int contador = 0, i, j, k, l, m, ind, aux, adj; /* Variáveis contadores para percorrer o arquivo e a string de leitura. */
    int carac, numLinhas = 0;                            /* Variável com o número de linhas do arquivo a serem lidas. */
    FILE *arquivo;
    long int **numeroMedidas, nmed;
    double Vbase = 1;
    double regua;
    char text_aux[500];
    dados = (char *)malloc(100);
    // Leitura dos dados de medidores
    strcpy(text_aux, folder);
    arquivo = fopen(strcat(text_aux, file), "r");
    

    //arquivo = fopen(folder,"r");

    // Verifica se O DMED foi lido
    if (arquivo == NULL)
    {
        printf("Erro ao abrir arquivo %s !!!\n", strcat(text_aux, file));
        exit(1);
    }

    //????
    numeroMedidas = (long int **)malloc(14 * sizeof(long int *));//matriz com o tipo de medidor em cada fase
    for (i = 0; i < 14; i++)
    {
        numeroMedidas[i] = (long int *)malloc(8 * sizeof(long int));
        for (j = 0; j < 8; j++)
        {
            numeroMedidas[i][j] = 0;
        }
    }

    //Aloca na memória espaço para a estrutura do DMED conforme o numero de linhas no arquivo
    while ((carac = fgetc(arquivo)) != EOF)
    {
        if (carac == '\n')
            numLinhas++;
    }
    rewind(arquivo);
    if (((*medidas) = (DMED *)malloc((numLinhas+2) * sizeof(DMED))) == NULL)
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para as medidas !!!!");
        exit(1);
    }

    // Le o arquivo de curva de cargas até o fim
    // Preenche a estrutura DMED
    while ((fgets(blocoLeitura, 2000, arquivo)) != NULL)
    {
        dados = blocoLeitura;

        (*medidas)[contador].ligado = atoi(getfield(dados, 1)); //se a medida esta ativa
        (*medidas)[contador].tipo = atoi(getfield(dados, 2));// tipo de medida 
        (*medidas)[contador].DE = atoi(getfield(dados, 3));// barra DE conforme nomeclatura do DBAR
        (*medidas)[contador].PARA = atoi(getfield(dados, 4));// barra PARA conforme a nomeclatura do DBAR
        (*medidas)[contador].fases = atoi(getfield(dados, 5));// De qual fase é a medida
        (*medidas)[contador].id = contador;// id da medida
        (*medidas)[contador].par = -1;// par ???

        (*medidas)[contador].h = 0; //valor do h(x)
        (*medidas)[contador].zmed = atof(getfield(dados, 6));// valor medida
        (*medidas)[contador].sigma = 1; //atof(getfield(dados,7));// sigma, desvio padrão calculado internamente
        (*medidas)[contador].prec = atof(getfield(dados, 7));// precisão do medidor

        numeroMedidas[(*medidas)[contador].tipo][(*medidas)[contador].fases - 1]++;//tipo de medidor em cada fase
        switch ((*medidas)[contador].tipo)// medidas que são em ramos ele atribui a numeração que serve de indice tanto no grafo quanto no dbar e no dram
        {
        case 0: // medidas que são em ramos ele atribui a numeração que serve de indice tanto no grafo quanto no dbar e no dram
        case 1:
        case 6:
        case 7:
        case 12: //PMU de corrente retangular real
        case 13: //PMU de corrente retangular imaginário
            for (i = 0; i < numeroRamos; i++)
            {
                if (((*medidas)[contador].DE == ramos[i].DE) && ((*medidas)[contador].PARA == ramos[i].PARA))
                {
                    (*medidas)[contador].k = ramos[i].k;
                    (*medidas)[contador].m = ramos[i].m;
                    (*medidas)[contador].ramo = i;
                }
                if (((*medidas)[contador].DE == ramos[i].PARA) && ((*medidas)[contador].PARA == ramos[i].DE))
                {
                    (*medidas)[contador].k = ramos[i].m;
                    (*medidas)[contador].m = ramos[i].k;
                    (*medidas)[contador].ramo = i;
                }
            }
            break;
        case 2: //Medidas nas barras ele atribui a numeração de DE para semelhante ao grafo  
        case 3:
        case 4:
        case 5:
        case 8:  //PMU de tensão retangular real
        case 9:  //PMU de tensão retangular imaginário
        case 10: //PMU de injeção de corrente retangular real
        case 11: //PMU de injeção de corrente retangular imaginário
            for (i = 0; i < numeroBarras; i++)
            {
                if ((*medidas)[contador].DE == barras[i].ID)
                {
                    (*medidas)[contador].k = barras[i].i;
                }
                (*medidas)[contador].m = -1;
                (*medidas)[contador].ramo = -1;
            }
            break;
        }
        contador++;
    }
    fclose(arquivo);

    //Associa as medidas ao grafo e transforma em pu os dados medidos

    //conta o numero de medidas
    nmed = 0;
    for (i = 0; i < 14; i++)
    {
        for (j = 0; j < 8; j++)
        {
            nmed = nmed + numeroMedidas[i][j];
        }
    }

    //--------------------------------------------------------------------------
    //Associa os medidores ao grafo e transforma medidas em pu
    for (i = 0; i < nmed; i++)
    {
        k = (*medidas)[i].k;
        m = (*medidas)[i].m;

        (*medidas)[i].idAlim = grafo[k].idAlim;

        //Associa a medida ao grafo e transforma em pu o valor medido e sigma
        switch ((*medidas)[i].tipo)
        {
        case 0: //Medida de Fluxo de Potência Ativa em kW
            (*medidas)[i].zmed = (*medidas)[i].zmed / (Sbase / 1000);//Transforma em PU
            (*medidas)[i].sigma = (*medidas)[i].sigma / (Sbase / 1000);

            for (j = 0; j < grafo[k].numeroAdjacentes; j++)// encontra o no adjacente ao qual aquela medida de fluxo está 
            {
                if (grafo[k].adjacentes[j].idNo == m)
                {
                    adj = j;
                }
            }
            ind = grafo[k].adjacentes[adj].nmed;//adiciona a medida de fluxo naquele nó para aquele adjacente
            grafo[k].adjacentes[adj].medidores[ind] = &(*medidas)[i];
            grafo[k].adjacentes[adj].nmed++;

            (*medidas)[i].nvar = 12;//considera que se tem 12 variáveis que aquela medida relaciona
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));//regua da Jacobiana
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));//aloca o espaço na jacobiana para todas as derivadas daquela medida
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;//inicia com valor zero
            }


            //coloca na régua as variáveis tanto do De qaunto do para
            regua = (double)k;//barra da medida 
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                (*medidas)[i].reguaH[j + 3] = -regua;
                regua += 0.1;
            }
            regua = (double)m;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j + 6] = regua;
                (*medidas)[i].reguaH[j + 9] = -regua;
                regua += 0.1;
            }
            break;
        case 1: //Medida de Fluxo de Potência Reativa em kVAr semelhante ao ativo
            (*medidas)[i].zmed = (*medidas)[i].zmed / (Sbase / 1000);
            (*medidas)[i].sigma = (*medidas)[i].sigma / (Sbase / 1000);

            for (j = 0; j < grafo[k].numeroAdjacentes; j++)
            {
                if (grafo[k].adjacentes[j].idNo == m)
                {
                    adj = j;
                }
            }
            ind = grafo[k].adjacentes[adj].nmed;
            grafo[k].adjacentes[adj].medidores[ind] = &(*medidas)[i];
            grafo[k].adjacentes[adj].nmed++;

            (*medidas)[i].nvar = 12;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                (*medidas)[i].reguaH[j + 3] = -regua;
                regua += 0.1;
            }
            regua = (double)m;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j + 6] = regua;
                (*medidas)[i].reguaH[j + 9] = -regua;
                regua += 0.1;
            }
            break;
        case 2: //Medida de Injeção de Potência Ativa em kW
            (*medidas)[i].zmed = (*medidas)[i].zmed / (Sbase / 1000);//passa medida pra p.u
            (*medidas)[i].sigma = (*medidas)[i].sigma / (Sbase / 1000);

            ind = grafo[k].nmed;//coloca a medida na barra do grafo
            grafo[k].medidores[ind] = &(*medidas)[i];
            grafo[k].nmed++;

            (*medidas)[i].nvar = 6 + 6 * grafo[k].numeroAdjacentes;//na estrutura medidas adiciona as variáveis das barras adjacentes
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));//aloca a régua para aquela medida
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));// aloca o espaço na jacobiana para aquela medida
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            //preenche a regua para aquela barra em questao
            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                (*medidas)[i].reguaH[j + 3] = -regua;
                regua += 0.1;
            }
            for (l = 0; l < grafo[k].numeroAdjacentes; l++)
            {
                regua = grafo[k].adjacentes[l].idNo;
                regua += 0.01;
                for (j = 0; j < 3; j++)
                {
                    (*medidas)[i].reguaH[6 + 6 * l + j] = regua;
                    (*medidas)[i].reguaH[6 + 6 * l + j + 3] = -regua;
                    regua += 0.1;
                }
            }
            break;
        case 3: //Medida de Injeção de Potência Reativa em kVAr | mesma coisa que a injeção ativa
            (*medidas)[i].zmed = (*medidas)[i].zmed / (Sbase / 1000);
            (*medidas)[i].sigma = (*medidas)[i].sigma / (Sbase / 1000);

            ind = grafo[k].nmed;
            grafo[k].medidores[ind] = &(*medidas)[i];
            grafo[k].nmed++;

            (*medidas)[i].nvar = 6 + 6 * grafo[k].numeroAdjacentes;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                (*medidas)[i].reguaH[j + 3] = -regua;
                regua += 0.1;
            }
            for (l = 0; l < grafo[k].numeroAdjacentes; l++)
            {
                regua = grafo[k].adjacentes[l].idNo;
                regua += 0.01;
                for (j = 0; j < 3; j++)
                {
                    (*medidas)[i].reguaH[6 + 6 * l + j] = regua;
                    (*medidas)[i].reguaH[6 + 6 * l + j + 3] = -regua;
                    regua += 0.1;
                }
            }
            break;
        case 4: //Medida de Magnitude de Tensão - kV
            switch ((*medidas)[i].fases)
            {
            case 1:
            case 2:
            case 3:
                Vbase = grafo[k].Vbase;//tensao base monofasica
                break;
            case 4:
            case 5:
            case 6:
                Vbase = grafo[k].Vbase * (pow(3, 0.5)); //tensao de base de linha
                break;
            }
            (*medidas)[i].zmed = (*medidas)[i].zmed / (Vbase / 1000);//passa pra p.u
            (*medidas)[i].sigma = (*medidas)[i].sigma / (Vbase / 1000);//passa pra p.u

            //add o medidor no grafo e na barra
            ind = grafo[k].nmed;
            grafo[k].medidores[ind] = &(*medidas)[i];
            grafo[k].nmed++;

            //seta o numero de variáveis 
            (*medidas)[i].nvar = 6;
            //aloca o espaco na regua
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            // aloca o espaço na linha da jacobiana
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            //preenche a regua

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                regua += 0.1;
            }
            regua = (double)k;
            regua += 0.01;
            for (j = 3; j < 6; j++)
            {
                (*medidas)[i].reguaH[j] = -regua;
                regua += 0.1;
            }
            break;
        case 5: //Medida de Ângulo de tensão - graus
            (*medidas)[i].zmed = (*medidas)[i].zmed * PI / 180;//passa o valor pra radianos
            (*medidas)[i].sigma = (*medidas)[i].sigma * PI / 180;
            ;

            // coloca a medida no grafo
            ind = grafo[k].nmed;
            grafo[k].medidores[ind] = &(*medidas)[i];
            grafo[k].nmed++;

            //aloca o espaço na estrutura medidas
            (*medidas)[i].nvar = 3;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            //inicia os malores com zero
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            //preenche a regua

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = -regua;
                regua += 0.1;
            }
            break;
        case 6: //Medida de Magnitude de Corrente em A - Fluxo
            Vbase = grafo[k].Vbase;
            (*medidas)[i].zmed = (*medidas)[i].zmed / ((Sbase / 1000) / (pow(3, 0.5) * (Vbase / 1000)));
            (*medidas)[i].sigma = (*medidas)[i].sigma / ((Sbase / 1000) / (pow(3, 0.5) * (Vbase / 1000)));

            for (j = 0; j < grafo[k].numeroAdjacentes; j++)
            {
                if (grafo[k].adjacentes[j].idNo == m)
                {
                    adj = j;
                }
            }
            ind = grafo[k].adjacentes[adj].nmed;
            grafo[k].adjacentes[adj].medidores[ind] = &(*medidas)[i];
            grafo[k].adjacentes[adj].nmed++;

            (*medidas)[i].nvar = 12;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                (*medidas)[i].reguaH[j + 3] = -regua;
                regua += 0.1;
            }
            regua = (double)m;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j + 6] = regua;
                (*medidas)[i].reguaH[j + 9] = -regua;
                regua += 0.1;
            }
            break;
        case 7: //Medida de Ângulo de Corrente em graus
            (*medidas)[i].zmed = (*medidas)[i].zmed * PI / 180;
            (*medidas)[i].sigma = (*medidas)[i].sigma * PI / 180;

            for (j = 0; j < grafo[k].numeroAdjacentes; j++)
            {
                if (grafo[k].adjacentes[j].idNo == m)
                {
                    adj = j;
                }
            }
            ind = grafo[k].adjacentes[adj].nmed;
            grafo[k].adjacentes[adj].medidores[ind] = &(*medidas)[i];
            grafo[k].adjacentes[adj].nmed++;

            (*medidas)[i].nvar = 12;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                (*medidas)[i].reguaH[j + 3] = -regua;
                regua += 0.1;
            }
            regua = (double)m;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j + 6] = regua;
                (*medidas)[i].reguaH[j + 9] = -regua;
                regua += 0.1;
            }
            break;
        case 8: //Medida PMU de tensão retangular Real
            switch ((*medidas)[i].fases)
            {
            case 1:
            case 2:
            case 3:
                Vbase = grafo[k].Vbase;
                break;
            case 4:
            case 5:
            case 6:
                Vbase = grafo[k].Vbase * (pow(3, 0.5));
                break;
            }
            (*medidas)[i].zmed = (*medidas)[i].zmed / (Vbase / 1000);
            (*medidas)[i].sigma = (*medidas)[i].sigma / (Vbase / 1000);

            ind = grafo[k].nmed;
            grafo[k].medidores[ind] = &(*medidas)[i];
            grafo[k].nmed++;

            (*medidas)[i].nvar = 3;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                regua += 0.1;
            }
            break;
        case 9: //Medida PMU de tensão retangular Imaginário
            switch ((*medidas)[i].fases)
            {
            case 1:
            case 2:
            case 3:
                Vbase = grafo[k].Vbase;
                break;
            case 4:
            case 5:
            case 6:
                Vbase = grafo[k].Vbase * (pow(3, 0.5));
                break;
            }
            (*medidas)[i].zmed = (*medidas)[i].zmed / (Vbase / 1000);
            (*medidas)[i].sigma = (*medidas)[i].sigma / (Vbase / 1000);

            ind = grafo[k].nmed;
            grafo[k].medidores[ind] = &(*medidas)[i];
            grafo[k].nmed++;

            (*medidas)[i].nvar = 3;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = -regua;
                regua += 0.1;
            }
            break;
        case 10: //Medida PMU de injeção de corrente retangular real
            Vbase = grafo[k].Vbase;
            (*medidas)[i].zmed = (*medidas)[i].zmed / ((Sbase / 1000) / ((Vbase / 1000)));
            (*medidas)[i].sigma = (*medidas)[i].sigma / ((Sbase / 1000) / ((Vbase / 1000)));

            ind = grafo[k].nmed;
            grafo[k].medidores[ind] = &(*medidas)[i];
            grafo[k].nmed++;

            (*medidas)[i].nvar = 6 + 6 * grafo[k].numeroAdjacentes;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                (*medidas)[i].reguaH[j + 3] = -regua;
                regua += 0.1;
            }
            for (l = 0; l < grafo[k].numeroAdjacentes; l++)
            {
                regua = grafo[k].adjacentes[l].idNo;
                regua += 0.01;
                for (j = 0; j < 3; j++)
                {
                    (*medidas)[i].reguaH[6 + 6 * l + j] = regua;
                    (*medidas)[i].reguaH[6 + 6 * l + j + 3] = -regua;
                    regua += 0.1;
                }
            }
            break;
        case 11: //Medida PMU de injeção de corrente retangular imaginário
            Vbase = grafo[k].Vbase;
            (*medidas)[i].zmed = (*medidas)[i].zmed / ((Sbase / 1000) / ((Vbase / 1000)));
            (*medidas)[i].sigma = (*medidas)[i].sigma / ((Sbase / 1000) / ((Vbase / 1000)));

            ind = grafo[k].nmed;
            grafo[k].medidores[ind] = &(*medidas)[i];
            grafo[k].nmed++;

            (*medidas)[i].nvar = 6 + 6 * grafo[k].numeroAdjacentes;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                (*medidas)[i].reguaH[j + 3] = -regua;
                regua += 0.1;
            }
            for (l = 0; l < grafo[k].numeroAdjacentes; l++)
            {
                regua = grafo[k].adjacentes[l].idNo;
                regua += 0.01;
                for (j = 0; j < 3; j++)
                {
                    (*medidas)[i].reguaH[6 + 6 * l + j] = regua;
                    (*medidas)[i].reguaH[6 + 6 * l + j + 3] = -regua;
                    regua += 0.1;
                }
            }

            break;
        case 12: //Medida PMU de corrente retangular real
            Vbase = grafo[k].Vbase;
            (*medidas)[i].zmed = (*medidas)[i].zmed / ((Sbase / 1000) / ((Vbase / 1000)));
            (*medidas)[i].sigma = (*medidas)[i].sigma / ((Sbase / 1000) / ((Vbase / 1000)));

            for (j = 0; j < grafo[k].numeroAdjacentes; j++)
            {
                if (grafo[k].adjacentes[j].idNo == m)
                {
                    adj = j;
                }
            }
            ind = grafo[k].adjacentes[adj].nmed;
            grafo[k].adjacentes[adj].medidores[ind] = &(*medidas)[i];
            grafo[k].adjacentes[adj].nmed++;

            (*medidas)[i].nvar = 12;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                (*medidas)[i].reguaH[j + 3] = -regua;
                regua += 0.1;
            }
            regua = (double)m;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j + 6] = regua;
                (*medidas)[i].reguaH[j + 9] = -regua;
                regua += 0.1;
            }
            break;
        case 13: //Medida PMU de corrente retangular Imaginário
            Vbase = grafo[k].Vbase;
            (*medidas)[i].zmed = (*medidas)[i].zmed / ((Sbase / 1000) / ((Vbase / 1000)));
            (*medidas)[i].sigma = (*medidas)[i].sigma / ((Sbase / 1000) / ((Vbase / 1000)));

            for (j = 0; j < grafo[k].numeroAdjacentes; j++)
            {
                if (grafo[k].adjacentes[j].idNo == m)
                {
                    adj = j;
                }
            }
            ind = grafo[k].adjacentes[adj].nmed;
            grafo[k].adjacentes[adj].medidores[ind] = &(*medidas)[i];
            grafo[k].adjacentes[adj].nmed++;

            (*medidas)[i].nvar = 12;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                (*medidas)[i].reguaH[j + 3] = -regua;
                regua += 0.1;
            }
            regua = (double)m;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j + 6] = regua;
                (*medidas)[i].reguaH[j + 9] = -regua;
                regua += 0.1;
            }
            break;
        }
    }   

    //--------------------------------------------------------------------------
    //Associa os medidores ao grafo e transforma medidas em pu (equivalente de corrente para o AMB e Baran precisa do par de medida)
    for (i = 0; i < nmed; i++)
    {
        k = (*medidas)[i].k;
        m = (*medidas)[i].m;

        //Associa a medida ao grafo e transforma em pu o valor medido e sigma
        switch ((*medidas)[i].tipo)
        {
        case 0: //Medida de Fluxo de Potência Ativa em kW
            for (j = 0; j < nmed; j++)
            {
                if (((*medidas)[j].tipo == 1) && ((*medidas)[j].k == k) && ((*medidas)[j].m == m))
                {
                    (*medidas)[i].par = j;
                    (*medidas)[j].par = i;
                    j = nmed;
                }
            }
            break;
        case 2: //Medida de Injeção de Potência Ativa em kW
            for (j = 0; j < nmed; j++)
            {
                if (((*medidas)[j].tipo == 3) && ((*medidas)[j].k == k) && ((*medidas)[j].m == m))
                {
                    (*medidas)[i].par = j;
                    (*medidas)[j].par = i;
                    j = nmed;
                }
            }
            break;
        case 4: //Medida de Magnitude de Tensão - kV
            for (j = 0; j < nmed; j++)
            {
                if (((*medidas)[j].tipo == 5) && ((*medidas)[j].k == k) && ((*medidas)[j].m == m))
                {
                    (*medidas)[i].par = j;
                    (*medidas)[j].par = i;
                    j = nmed;
                }
            }
            break;
        case 6: //Medida de Magnitude de Corrente em A
            for (j = 0; j < nmed; i++)
            {
                if (((*medidas)[j].tipo == 7) && ((*medidas)[j].k == k) && ((*medidas)[j].m == m))
                {
                    (*medidas)[i].par = j;
                    (*medidas)[j].par = i;
                    j = nmed;
                }
            }
            break;
        }
    }

    return (numeroMedidas);
}

long int **LeituraFP(char *folder, DMED **medidas, DRAM *ramos, long int numeroRamos, DBAR *barras, long int numeroBarras, GRAFO *grafo, double Sbase)
{
    char blocoLeitura[2000];                             /* Variável para realizar a leitura do bloco de caracteres do arquivo. */
    char *dados;
    int tipoBFP;                                         /* Variável do tipo ponteiro para char, utilizada para alterar o ponteiro da string lida do arquivo de forma a realizar o loop no sscanf. */
    long int contador = 0, i, j, k, l, m, ind, aux, adj;; /* Variáveis contadores para percorrer o arquivo e a string de leitura. */
    int carac, numLinhas = 0;                            /* Variável com o número de linhas do arquivo a serem lidas. */
    FILE *arquivo;
    long int **numeroMedidas, nmed,ID_FP;
    double Vbase = 1;
    double regua;
    double Va,Vb,Vc,Ta,Tb,Tc;
    char text_aux[500];
    char file[100]="loads.csv";
    char file2[100]="DBTYPE.csv";
    dados = (char *)malloc(100);
    // Leitura dos dados de medidores
    strcpy(text_aux, folder);
    arquivo = fopen(strcat(text_aux, file), "r");
    

    //arquivo = fopen(folder,"r");

    // Verifica se O DMED foi lido
    if (arquivo == NULL)
    {
        printf("Erro ao abrir arquivo %s !!!\n", strcat(text_aux, file));
        exit(1);
    }

    //????
    numeroMedidas = (long int **)malloc(14 * sizeof(long int *));//matriz com o tipo de medidor em cada fase
    for (i = 0; i < 14; i++)
    {
        numeroMedidas[i] = (long int *)malloc(8 * sizeof(long int));
        for (j = 0; j < 8; j++)
        {
            numeroMedidas[i][j] = 0;
        }
    }

    //Aloca na memória espaço para a estrutura do DMED conforme o numero de linhas no arquivo
    while ((carac = fgetc(arquivo)) != EOF)
    {
        if (carac == '\n')
            numLinhas++;
    }
    rewind(arquivo);
    if (((*medidas) = (DMED *)malloc((numLinhas+2) * sizeof(DMED))) == NULL)
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para as medidas !!!!");
        exit(1);
    }

    // Le o arquivo de curva de cargas até o fim
    // Preenche a estrutura DMED
    while ((fgets(blocoLeitura, 2000, arquivo)) != NULL)
    {
        dados = blocoLeitura;

        (*medidas)[contador].ligado = atoi(getfield(dados, 1)); //se a medida esta ativa
        (*medidas)[contador].tipo = atoi(getfield(dados, 2));// tipo de medida 
        (*medidas)[contador].DE = atoi(getfield(dados, 3));// barra DE conforme nomeclatura do DBAR
        (*medidas)[contador].PARA = atoi(getfield(dados, 4));// barra PARA conforme a nomeclatura do DBAR
        (*medidas)[contador].fases = atoi(getfield(dados, 5));// De qual fase é a medida
        (*medidas)[contador].id = contador;// id da medida
        (*medidas)[contador].par = -1;// par ???

        (*medidas)[contador].h = 0; //valor do h(x)
        (*medidas)[contador].zmed = atof(getfield(dados, 6));// valor medida
        (*medidas)[contador].sigma = 1; //atof(getfield(dados,7));// sigma, desvio padrão calculado internamente
        (*medidas)[contador].prec = atof(getfield(dados, 7));// precisão do medidor

        numeroMedidas[(*medidas)[contador].tipo][(*medidas)[contador].fases - 1]++;//tipo de medidor em cada fase
        switch ((*medidas)[contador].tipo)// medidas que são em ramos ele atribui a numeração que serve de indice tanto no grafo quanto no dbar e no dram
        {
        case 0: // medidas que são em ramos ele atribui a numeração que serve de indice tanto no grafo quanto no dbar e no dram
        case 1:
        case 6:
        case 7:
        case 12: //PMU de corrente retangular real
        case 13: //PMU de corrente retangular imaginário
            for (i = 0; i < numeroRamos; i++)
            {
                if (((*medidas)[contador].DE == ramos[i].DE) && ((*medidas)[contador].PARA == ramos[i].PARA))
                {
                    (*medidas)[contador].k = ramos[i].k;
                    (*medidas)[contador].m = ramos[i].m;
                    (*medidas)[contador].ramo = i;
                }
                if (((*medidas)[contador].DE == ramos[i].PARA) && ((*medidas)[contador].PARA == ramos[i].DE))
                {
                    (*medidas)[contador].k = ramos[i].m;
                    (*medidas)[contador].m = ramos[i].k;
                    (*medidas)[contador].ramo = i;
                }
            }
            break;
        case 2: //Medidas nas barras ele atribui a numeração de DE para semelhante ao grafo  
        case 3:
        case 4:
        case 5:
        case 8:  //PMU de tensão retangular real
        case 9:  //PMU de tensão retangular imaginário
        case 10: //PMU de injeção de corrente retangular real
        case 11: //PMU de injeção de corrente retangular imaginário
            for (i = 0; i < numeroBarras; i++)
            {
                if ((*medidas)[contador].DE == barras[i].ID)
                {
                    (*medidas)[contador].k = barras[i].i;
                }
                (*medidas)[contador].m = -1;
                (*medidas)[contador].ramo = -1;
            }
            break;
        }
        contador++;
    }
    fclose(arquivo);

    //Associa as medidas ao grafo e transforma em pu os dados medidos

    //conta o numero de medidas
    nmed = 0;
    for (i = 0; i < 14; i++)
    {
        for (j = 0; j < 8; j++)
        {
            nmed = nmed + numeroMedidas[i][j];
        }
    }

    //--------------------------------------------------------------------------
    //Associa os medidores ao grafo e transforma medidas em pu
    for (i = 0; i < nmed; i++)
    {
        k = (*medidas)[i].k;
        m = (*medidas)[i].m;

        (*medidas)[i].idAlim = grafo[k].idAlim;

        //Associa a medida ao grafo e transforma em pu o valor medido e sigma
        switch ((*medidas)[i].tipo)
        {
        case 0: //Medida de Fluxo de Potência Ativa em kW
            (*medidas)[i].zmed = (*medidas)[i].zmed / (Sbase / 1000);//Transforma em PU
            (*medidas)[i].sigma = (*medidas)[i].sigma / (Sbase / 1000);

            for (j = 0; j < grafo[k].numeroAdjacentes; j++)// encontra o no adjacente ao qual aquela medida de fluxo está 
            {
                if (grafo[k].adjacentes[j].idNo == m)
                {
                    adj = j;
                }
            }
            ind = grafo[k].adjacentes[adj].nmed;//adiciona a medida de fluxo naquele nó para aquele adjacente
            grafo[k].adjacentes[adj].medidores[ind] = &(*medidas)[i];
            grafo[k].adjacentes[adj].nmed++;

            (*medidas)[i].nvar = 12;//considera que se tem 12 variáveis que aquela medida relaciona
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));//regua da Jacobiana
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));//aloca o espaço na jacobiana para todas as derivadas daquela medida
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;//inicia com valor zero
            }


            //coloca na régua as variáveis tanto do De qaunto do para
            regua = (double)k;//barra da medida 
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                (*medidas)[i].reguaH[j + 3] = -regua;
                regua += 0.1;
            }
            regua = (double)m;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j + 6] = regua;
                (*medidas)[i].reguaH[j + 9] = -regua;
                regua += 0.1;
            }
            break;
        case 1: //Medida de Fluxo de Potência Reativa em kVAr semelhante ao ativo
            (*medidas)[i].zmed = (*medidas)[i].zmed / (Sbase / 1000);
            (*medidas)[i].sigma = (*medidas)[i].sigma / (Sbase / 1000);

            for (j = 0; j < grafo[k].numeroAdjacentes; j++)
            {
                if (grafo[k].adjacentes[j].idNo == m)
                {
                    adj = j;
                }
            }
            ind = grafo[k].adjacentes[adj].nmed;
            grafo[k].adjacentes[adj].medidores[ind] = &(*medidas)[i];
            grafo[k].adjacentes[adj].nmed++;

            (*medidas)[i].nvar = 12;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                (*medidas)[i].reguaH[j + 3] = -regua;
                regua += 0.1;
            }
            regua = (double)m;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j + 6] = regua;
                (*medidas)[i].reguaH[j + 9] = -regua;
                regua += 0.1;
            }
            break;
        case 2: //Medida de Injeção de Potência Ativa em kW
            (*medidas)[i].zmed = (*medidas)[i].zmed / (Sbase / 1000);//passa medida pra p.u
            (*medidas)[i].sigma = (*medidas)[i].sigma / (Sbase / 1000);

            ind = grafo[k].nmed;//coloca a medida na barra do grafo
            grafo[k].medidores[ind] = &(*medidas)[i];
            grafo[k].nmed++;

            (*medidas)[i].nvar = 6 + 6 * grafo[k].numeroAdjacentes;//na estrutura medidas adiciona as variáveis das barras adjacentes
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));//aloca a régua para aquela medida
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));// aloca o espaço na jacobiana para aquela medida
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            //preenche a regua para aquela barra em questao
            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                (*medidas)[i].reguaH[j + 3] = -regua;
                regua += 0.1;
            }
            for (l = 0; l < grafo[k].numeroAdjacentes; l++)
            {
                regua = grafo[k].adjacentes[l].idNo;
                regua += 0.01;
                for (j = 0; j < 3; j++)
                {
                    (*medidas)[i].reguaH[6 + 6 * l + j] = regua;
                    (*medidas)[i].reguaH[6 + 6 * l + j + 3] = -regua;
                    regua += 0.1;
                }
            }
            break;
        case 3: //Medida de Injeção de Potência Reativa em kVAr | mesma coisa que a injeção ativa
            (*medidas)[i].zmed = (*medidas)[i].zmed / (Sbase / 1000);
            (*medidas)[i].sigma = (*medidas)[i].sigma / (Sbase / 1000);

            ind = grafo[k].nmed;
            grafo[k].medidores[ind] = &(*medidas)[i];
            grafo[k].nmed++;

            (*medidas)[i].nvar = 6 + 6 * grafo[k].numeroAdjacentes;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                (*medidas)[i].reguaH[j + 3] = -regua;
                regua += 0.1;
            }
            for (l = 0; l < grafo[k].numeroAdjacentes; l++)
            {
                regua = grafo[k].adjacentes[l].idNo;
                regua += 0.01;
                for (j = 0; j < 3; j++)
                {
                    (*medidas)[i].reguaH[6 + 6 * l + j] = regua;
                    (*medidas)[i].reguaH[6 + 6 * l + j + 3] = -regua;
                    regua += 0.1;
                }
            }
            break;
        case 4: //Medida de Magnitude de Tensão - kV
            switch ((*medidas)[i].fases)
            {
            case 1:
            case 2:
            case 3:
                Vbase = grafo[k].Vbase;//tensao base monofasica
                break;
            case 4:
            case 5:
            case 6:
                Vbase = grafo[k].Vbase * (pow(3, 0.5)); //tensao de base de linha
                break;
            }
            (*medidas)[i].zmed = (*medidas)[i].zmed / (Vbase / 1000);//passa pra p.u
            (*medidas)[i].sigma = (*medidas)[i].sigma / (Vbase / 1000);//passa pra p.u

            //add o medidor no grafo e na barra
            ind = grafo[k].nmed;
            grafo[k].medidores[ind] = &(*medidas)[i];
            grafo[k].nmed++;

            //seta o numero de variáveis 
            (*medidas)[i].nvar = 6;
            //aloca o espaco na regua
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            // aloca o espaço na linha da jacobiana
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            //preenche a regua

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                regua += 0.1;
            }
            regua = (double)k;
            regua += 0.01;
            for (j = 3; j < 6; j++)
            {
                (*medidas)[i].reguaH[j] = -regua;
                regua += 0.1;
            }
            break;
        case 5: //Medida de Ângulo de tensão - graus
            (*medidas)[i].zmed = (*medidas)[i].zmed * PI / 180;//passa o valor pra radianos
            (*medidas)[i].sigma = (*medidas)[i].sigma * PI / 180;
            ;

            // coloca a medida no grafo
            ind = grafo[k].nmed;
            grafo[k].medidores[ind] = &(*medidas)[i];
            grafo[k].nmed++;

            //aloca o espaço na estrutura medidas
            (*medidas)[i].nvar = 3;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            //inicia os malores com zero
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            //preenche a regua

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = -regua;
                regua += 0.1;
            }
            break;
        case 6: //Medida de Magnitude de Corrente em A - Fluxo
            Vbase = grafo[k].Vbase;
            (*medidas)[i].zmed = (*medidas)[i].zmed / ((Sbase / 1000) / (pow(3, 0.5) * (Vbase / 1000)));
            (*medidas)[i].sigma = (*medidas)[i].sigma / ((Sbase / 1000) / (pow(3, 0.5) * (Vbase / 1000)));

            for (j = 0; j < grafo[k].numeroAdjacentes; j++)
            {
                if (grafo[k].adjacentes[j].idNo == m)
                {
                    adj = j;
                }
            }
            ind = grafo[k].adjacentes[adj].nmed;
            grafo[k].adjacentes[adj].medidores[ind] = &(*medidas)[i];
            grafo[k].adjacentes[adj].nmed++;

            (*medidas)[i].nvar = 12;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                (*medidas)[i].reguaH[j + 3] = -regua;
                regua += 0.1;
            }
            regua = (double)m;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j + 6] = regua;
                (*medidas)[i].reguaH[j + 9] = -regua;
                regua += 0.1;
            }
            break;
        case 7: //Medida de Ângulo de Corrente em graus
            (*medidas)[i].zmed = (*medidas)[i].zmed * PI / 180;
            (*medidas)[i].sigma = (*medidas)[i].sigma * PI / 180;

            for (j = 0; j < grafo[k].numeroAdjacentes; j++)
            {
                if (grafo[k].adjacentes[j].idNo == m)
                {
                    adj = j;
                }
            }
            ind = grafo[k].adjacentes[adj].nmed;
            grafo[k].adjacentes[adj].medidores[ind] = &(*medidas)[i];
            grafo[k].adjacentes[adj].nmed++;

            (*medidas)[i].nvar = 12;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                (*medidas)[i].reguaH[j + 3] = -regua;
                regua += 0.1;
            }
            regua = (double)m;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j + 6] = regua;
                (*medidas)[i].reguaH[j + 9] = -regua;
                regua += 0.1;
            }
            break;
        case 8: //Medida PMU de tensão retangular Real
            switch ((*medidas)[i].fases)
            {
            case 1:
            case 2:
            case 3:
                Vbase = grafo[k].Vbase;
                break;
            case 4:
            case 5:
            case 6:
                Vbase = grafo[k].Vbase * (pow(3, 0.5));
                break;
            }
            (*medidas)[i].zmed = (*medidas)[i].zmed / (Vbase / 1000);
            (*medidas)[i].sigma = (*medidas)[i].sigma / (Vbase / 1000);

            ind = grafo[k].nmed;
            grafo[k].medidores[ind] = &(*medidas)[i];
            grafo[k].nmed++;

            (*medidas)[i].nvar = 3;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                regua += 0.1;
            }
            break;
        case 9: //Medida PMU de tensão retangular Imaginário
            switch ((*medidas)[i].fases)
            {
            case 1:
            case 2:
            case 3:
                Vbase = grafo[k].Vbase;
                break;
            case 4:
            case 5:
            case 6:
                Vbase = grafo[k].Vbase * (pow(3, 0.5));
                break;
            }
            (*medidas)[i].zmed = (*medidas)[i].zmed / (Vbase / 1000);
            (*medidas)[i].sigma = (*medidas)[i].sigma / (Vbase / 1000);

            ind = grafo[k].nmed;
            grafo[k].medidores[ind] = &(*medidas)[i];
            grafo[k].nmed++;

            (*medidas)[i].nvar = 3;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = -regua;
                regua += 0.1;
            }
            break;
        case 10: //Medida PMU de injeção de corrente retangular real
            Vbase = grafo[k].Vbase;
            (*medidas)[i].zmed = (*medidas)[i].zmed / ((Sbase / 1000) / ((Vbase / 1000)));
            (*medidas)[i].sigma = (*medidas)[i].sigma / ((Sbase / 1000) / ((Vbase / 1000)));

            ind = grafo[k].nmed;
            grafo[k].medidores[ind] = &(*medidas)[i];
            grafo[k].nmed++;

            (*medidas)[i].nvar = 6 + 6 * grafo[k].numeroAdjacentes;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                (*medidas)[i].reguaH[j + 3] = -regua;
                regua += 0.1;
            }
            for (l = 0; l < grafo[k].numeroAdjacentes; l++)
            {
                regua = grafo[k].adjacentes[l].idNo;
                regua += 0.01;
                for (j = 0; j < 3; j++)
                {
                    (*medidas)[i].reguaH[6 + 6 * l + j] = regua;
                    (*medidas)[i].reguaH[6 + 6 * l + j + 3] = -regua;
                    regua += 0.1;
                }
            }
            break;
        case 11: //Medida PMU de injeção de corrente retangular imaginário
            Vbase = grafo[k].Vbase;
            (*medidas)[i].zmed = (*medidas)[i].zmed / ((Sbase / 1000) / ((Vbase / 1000)));
            (*medidas)[i].sigma = (*medidas)[i].sigma / ((Sbase / 1000) / ((Vbase / 1000)));

            ind = grafo[k].nmed;
            grafo[k].medidores[ind] = &(*medidas)[i];
            grafo[k].nmed++;

            (*medidas)[i].nvar = 6 + 6 * grafo[k].numeroAdjacentes;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                (*medidas)[i].reguaH[j + 3] = -regua;
                regua += 0.1;
            }
            for (l = 0; l < grafo[k].numeroAdjacentes; l++)
            {
                regua = grafo[k].adjacentes[l].idNo;
                regua += 0.01;
                for (j = 0; j < 3; j++)
                {
                    (*medidas)[i].reguaH[6 + 6 * l + j] = regua;
                    (*medidas)[i].reguaH[6 + 6 * l + j + 3] = -regua;
                    regua += 0.1;
                }
            }

            break;
        case 12: //Medida PMU de corrente retangular real
            Vbase = grafo[k].Vbase;
            (*medidas)[i].zmed = (*medidas)[i].zmed / ((Sbase / 1000) / ((Vbase / 1000)));
            (*medidas)[i].sigma = (*medidas)[i].sigma / ((Sbase / 1000) / ((Vbase / 1000)));

            for (j = 0; j < grafo[k].numeroAdjacentes; j++)
            {
                if (grafo[k].adjacentes[j].idNo == m)
                {
                    adj = j;
                }
            }
            ind = grafo[k].adjacentes[adj].nmed;
            grafo[k].adjacentes[adj].medidores[ind] = &(*medidas)[i];
            grafo[k].adjacentes[adj].nmed++;

            (*medidas)[i].nvar = 12;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                (*medidas)[i].reguaH[j + 3] = -regua;
                regua += 0.1;
            }
            regua = (double)m;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j + 6] = regua;
                (*medidas)[i].reguaH[j + 9] = -regua;
                regua += 0.1;
            }
            break;
        case 13: //Medida PMU de corrente retangular Imaginário
            Vbase = grafo[k].Vbase;
            (*medidas)[i].zmed = (*medidas)[i].zmed / ((Sbase / 1000) / ((Vbase / 1000)));
            (*medidas)[i].sigma = (*medidas)[i].sigma / ((Sbase / 1000) / ((Vbase / 1000)));

            for (j = 0; j < grafo[k].numeroAdjacentes; j++)
            {
                if (grafo[k].adjacentes[j].idNo == m)
                {
                    adj = j;
                }
            }
            ind = grafo[k].adjacentes[adj].nmed;
            grafo[k].adjacentes[adj].medidores[ind] = &(*medidas)[i];
            grafo[k].adjacentes[adj].nmed++;

            (*medidas)[i].nvar = 12;
            (*medidas)[i].reguaH = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            (*medidas)[i].H = (double *)malloc((*medidas)[i].nvar * sizeof(double));
            for (j = 0; j < (*medidas)[i].nvar; j++)
            {
                (*medidas)[i].H[j] = 0;
            }

            regua = (double)k;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j] = regua;
                (*medidas)[i].reguaH[j + 3] = -regua;
                regua += 0.1;
            }
            regua = (double)m;
            regua += 0.01;
            for (j = 0; j < 3; j++)
            {
                (*medidas)[i].reguaH[j + 6] = regua;
                (*medidas)[i].reguaH[j + 9] = -regua;
                regua += 0.1;
            }
            break;
        }
    }   

    //--------------------------------------------------------------------------
    //Associa os medidores ao grafo e transforma medidas em pu (equivalente de corrente para o AMB e Baran precisa do par de medida)
    for (i = 0; i < nmed; i++)
    {
        k = (*medidas)[i].k;
        m = (*medidas)[i].m;

        //Associa a medida ao grafo e transforma em pu o valor medido e sigma
        switch ((*medidas)[i].tipo)
        {
        case 0: //Medida de Fluxo de Potência Ativa em kW
            for (j = 0; j < nmed; j++)
            {
                if (((*medidas)[j].tipo == 1) && ((*medidas)[j].k == k) && ((*medidas)[j].m == m))
                {
                    (*medidas)[i].par = j;
                    (*medidas)[j].par = i;
                    j = nmed;
                }
            }
            break;
        case 2: //Medida de Injeção de Potência Ativa em kW
            for (j = 0; j < nmed; j++)
            {
                if (((*medidas)[j].tipo == 3) && ((*medidas)[j].k == k) && ((*medidas)[j].m == m))
                {
                    (*medidas)[i].par = j;
                    (*medidas)[j].par = i;
                    j = nmed;
                }
            }
            break;
        case 4: //Medida de Magnitude de Tensão - kV
            for (j = 0; j < nmed; j++)
            {
                if (((*medidas)[j].tipo == 5) && ((*medidas)[j].k == k) && ((*medidas)[j].m == m))
                {
                    (*medidas)[i].par = j;
                    (*medidas)[j].par = i;
                    j = nmed;
                }
            }
            break;
        case 6: //Medida de Magnitude de Corrente em A
            for (j = 0; j < nmed; i++)
            {
                if (((*medidas)[j].tipo == 7) && ((*medidas)[j].k == k) && ((*medidas)[j].m == m))
                {
                    (*medidas)[i].par = j;
                    (*medidas)[j].par = i;
                    j = nmed;
                }
            }
            break;
        }
    }

    // Atritbui o tipo das barras
    strcpy(text_aux, folder);
    arquivo = fopen(strcat(text_aux, file2), "r");

    if (arquivo == NULL)
    {
        printf("Erro ao abrir arquivo de Tipos de Barras !!!\n");
        exit(1);
    }

     while ((fgets(blocoLeitura, 2000, arquivo)) != NULL)
    {
        dados=blocoLeitura;
        ID_FP =atoi(getfield(dados,1));
        tipoBFP = atoi(getfield(dados,2));

        Va=atof(getfield(dados,3));
        Vb=atof(getfield(dados,4));
        Vc=atof(getfield(dados,5));
        Ta=atof(getfield(dados,6));
        Tb=atof(getfield(dados,7));
        Tc=atof(getfield(dados,8));

        for (k=0;k<numeroBarras;k++)
        {
           if(grafo[k].barra->ID==ID_FP)
           {
               grafo[k].tipo=tipoBFP;
               if (tipoBFP==1)
               {
                   switch ( grafo[k].fases)
                   {
                    case 1:
                        grafo[k].barra->Vinicial[0]=Va;
                        break;
                    case 2:
                        grafo[k].barra->Vinicial[1]=Vb*(cos(-120*PI/180) + I*sin(-120*PI/180));
                        break;
                    case 3:
                        grafo[k].barra->Vinicial[2]=Vc*(cos(120*PI/180) + I*sin(120*PI/180));
                        j++;
                        break;
                    case 4:
                        grafo[k].barra->Vinicial[0]=Va;
                        grafo[k].barra->Vinicial[1]=Vb*(cos(-120*PI/180) + I*sin(120*PI/180));
                        break;
                    case 5:
                        grafo[k].barra->Vinicial[0]=Va;
                        grafo[k].barra->Vinicial[2]=Vc*(cos(120*PI/180) + I*sin(120*PI/180));
                        break;
                    case 6:
                        grafo[k].barra->Vinicial[1]=Vb*(cos(-120*PI/180) + I*sin(-120*PI/180));
                        grafo[k].barra->Vinicial[2]=Vc*(cos(120*PI/180) + I*sin(120*PI/180));
                        break;
                    case 7:
                        grafo[k].barra->Vinicial[0]=Va;
                        grafo[k].barra->Vinicial[1]=Vb*(cos(-120*PI/180) + I*sin(-120*PI/180));
                        grafo[k].barra->Vinicial[2]=Vc*(cos(120*PI/180) + I*sin(120*PI/180));
                        break;
                   }
               }
                if (tipoBFP==2)
               {
                   switch ( grafo[k].fases)
                   {
                    case 1:
                        grafo[k].barra->Vinicial[0]=Va;
                        break;
                    case 2:
                        grafo[k].barra->Vinicial[1]=Vb*(cos(Tb*PI/180) + I*sin(Tb*PI/180));
                        break;
                    case 3:
                        grafo[k].barra->Vinicial[2]=Vc*(cos(Tc*PI/180) + I*sin(Tc*PI/180));
                        j++;
                        break;
                    case 4:
                        grafo[k].barra->Vinicial[0]=Va;
                        grafo[k].barra->Vinicial[1]=Vb*(cos(Tb*PI/180) + I*sin(Tc*PI/180));
                        break;
                    case 5:
                        grafo[k].barra->Vinicial[0]=Va;
                        grafo[k].barra->Vinicial[2]=Vc*(cos(Tc*PI/180) + I*sin(Tc*PI/180));
                        break;
                    case 6:
                        grafo[k].barra->Vinicial[1]=Vb*(cos(Tb*PI/180) + I*sin(Tb*PI/180));
                        grafo[k].barra->Vinicial[2]=Vc*(cos(Tc*PI/180) + I*sin(Tc*PI/180));
                        break;
                    case 7:
                        grafo[k].barra->Vinicial[0]=Va;
                        grafo[k].barra->Vinicial[1]=Vb*(cos(Tb*PI/180) + I*sin(Tb*PI/180));
                        grafo[k].barra->Vinicial[2]=Vc*(cos(Tc*PI/180) + I*sin(Tc*PI/180));
                        break;
                   }
               }
               break;
           }
        }

    }


    return (numeroMedidas);
}

// Gera estrutura de dados do grafo que representa a rede elétricas
void geraGrafo(GRAFO **grafo, DBAR *barras, long int numeroBarras, DRAM *ramos, long int numeroRamos)
{
    long int i, j, k;
    long int barraDe, barraPara, nadj;

    if (((*grafo) = (GRAFO *)malloc((numeroBarras) * sizeof(GRAFO))) == NULL)//aloca a memória para o grafo da rede
    {
        printf("Erro -- Nao foi possivel alocar espaco de memoria para as barras !!!!");
        exit(1);
    }

    for (i = 0; i < numeroBarras; i++)
    {
        (*grafo)[i].idNo = i; //Cada barra do grafo é um nó e recebe a identificação como a ordem do DBAR
        (*grafo)[i].tipo = barras[i].tipo;//tipo recebe 0 hard coded
        (*grafo)[i].fases = barras[i].fases;//Número de fases (1=A; 2=B; 3=C; 4=AB; 5=CA; 6=BC; 7=ABC);
        (*grafo)[i].Vbase = barras[i].Vbase;//Tensão nominal da barra em V dividido por raiz de 3 para encontrar de fase
        (*grafo)[i].barra = &barras[i];//Passa a estrutura do DBAR inteira do DBAR para aquela barra do grafo
        (*grafo)[i].nmed = 0;//numero de medidas na barra
        (*grafo)[i].numeroAdjacentes = 0;//numero de barras adjacentes

        (*grafo)[i].medidores = (DMED **)malloc(30 * sizeof(DMED *));//Aloca o espaço para 30 medidores "hard coded"
        for (j = 0; j < 3; j++)//Atribui valor nulo para os shunts da barra
        {
            for (k = 0; k < 3; k++)
            {
                (*grafo)[i].Ysh[j][k] = 0;
            }
        }
    }

    for (i = 0; i < numeroRamos; i++)//percorre todos os ramos
    {
        barraDe = ramos[i].k;//Salva o ID da barra De ou k
        barraPara = ramos[i].m;//Salva o ID da barra Para ou m
        nadj = (*grafo)[barraDe].numeroAdjacentes;//Recebe o número de adjacentes da barra De

        (*grafo)[barraDe].adjacentes[nadj].ramo = &ramos[i];//Salva o ramo como um adjacente daquela barra do Grafo, salva a estrutura do DRAM
        (*grafo)[barraDe].adjacentes[nadj].idNo = barraPara;//Salva o ID da barra para
        (*grafo)[barraDe].adjacentes[nadj].tipo = ramos[i].tipo;//Salva o tipo
        (*grafo)[barraDe].adjacentes[nadj].estado = ramos[i].estado;//Salva o estado
        (*grafo)[barraDe].adjacentes[nadj].nmed = 0;//Salva o número de medidas
        (*grafo)[barraDe].adjacentes[nadj].idram = i;//atribui o Id do ramo de acordo com a estrutura do DRAM
        (*grafo)[barraDe].numeroAdjacentes++;//incrementa o número de adjacentes
        if (ramos[i].tipo == 1)//Se for transaformador, guarda a elação de transformação
        {
            (*grafo)[barraDe].adjacentes[nadj].relacao = ramos[i].trafo.Vsec / ramos[i].trafo.Vpri;
        }
        (*grafo)[barraDe].adjacentes[nadj].medidores = (DMED **)malloc(30 * sizeof(DMED *));//aloca uma estrutura para 30 medidores dentro 
        //faz a mesma coisa que fez para a barra DE, mas para a barra para
        nadj = (*grafo)[barraPara].numeroAdjacentes;//Recebe o número de adjacentes da barra Para
        (*grafo)[barraPara].adjacentes[nadj].ramo = &ramos[i];//Salva o ramo como um adjacente daquela barra do Grafo, salva a estrutura do DRAM
        (*grafo)[barraPara].adjacentes[nadj].idNo = barraDe;//Salva o ID da narra DE
        (*grafo)[barraPara].adjacentes[nadj].tipo = ramos[i].tipo;//Salva o tipo 0 para linha 1 trafo...
        (*grafo)[barraPara].adjacentes[nadj].estado = ramos[i].estado;//Salva o estado
        (*grafo)[barraPara].adjacentes[nadj].nmed = 0;//Inicia o numero de medidas
        (*grafo)[barraPara].adjacentes[nadj].idram = i;//atribui o Id do ramo de acordo com a estrutura do DRAM
        (*grafo)[barraPara].numeroAdjacentes++;//incrementa o número de adjacentes
        if (ramos[i].tipo == 1)//Se for transaformador, guarda a elação de transformação
        {
            (*grafo)[barraPara].adjacentes[nadj].relacao = ramos[i].trafo.Vpri / ramos[i].trafo.Vsec;
        }
        (*grafo)[barraPara].adjacentes[nadj].medidores = (DMED **)malloc(30 * sizeof(DMED *));
    }
}
//------------------------------------------------------------------------------
//
// FUNÇÕES DE IMPRESSÃO DE DADOS
//
//------------------------------------------------------------------------------
// Imprime arquivo com elementos lidos
void salvaDadosRedeEletrica(DBAR *barras, long int numeroBarras, DRAM *ramos, long int numeroRamos, DMED *medidas, long int **numeroMedidas)
{
    long int i, j, nmed;
    FILE *arquivo;

    nmed = 0;
    for (i = 0; i < 14; i++)
    {
        for (j = 0; j < 8; j++)
        {
            nmed = nmed + numeroMedidas[i][j];
        }
    }

    // Leitura dos dados de barras
    arquivo = fopen("dadosRedeEletrica.dad", "w");
    if (arquivo == NULL)
    {
        printf("Erro ao abrir arquivo dadosRedeEletrica.dad !!!\n");
        exit(1);
    }

    fprintf(arquivo, "LISTA DE BARRAS:\n");
    for (i = 0; i < numeroBarras; i++)
    {
        fprintf(arquivo, "Barra: %d\t%d\t%s\t%.2lf", barras[i].i, barras[i].ID, charFases(barras[i].fases), barras[i].Vbase);
        if (barras[i].nloads > 0)
        {
            for (j = 0; j < barras[i].nloads; j++)
                fprintf(arquivo, "\n\tCarga %d:\t%s\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.1lf", j, charLigacao(barras[i].loads[j].lig), barras[i].loads[j].Pnom[0], barras[i].loads[j].Pnom[1], barras[i].loads[j].Pnom[2], barras[i].loads[j].Qnom[0], barras[i].loads[j].Qnom[1], barras[i].loads[j].Qnom[2], barras[i].loads[j].ZIP);
        }
        if (barras[i].nshunts > 0)
        {
            for (j = 0; j < barras[i].nshunts; j++)
            {
                fprintf(arquivo, "\n\tShunt %d:\t%s\t%.5lf\t%.5lf\t%.5lf\t%d", j, charLigacao(barras[i].shunts[j].lig), barras[i].shunts[j].Qnom[0], barras[i].shunts[j].Qnom[1], barras[i].shunts[j].Qnom[2], barras[i].shunts[j].controle);
                if (barras[i].shunts[j].controle != 0)
                    fprintf(arquivo, "\t%.5lf\t%.5lf\t%.5lf\t%.5lf", barras[i].shunts[j].DV, barras[i].shunts[j].Vset[0], barras[i].shunts[j].Vset[1], barras[i].shunts[j].Vset[2]);
            }
        }
        if (barras[i].ngds > 0)
        {
            for (j = 0; j < barras[i].ngds; j++)
            {
                fprintf(arquivo, "\n\tGD %d:\t%s\t%.5lf kVA\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.1lf", j, charLigacao(barras[i].gds[j].lig), barras[i].gds[j].Snominal, barras[i].gds[j].Pnom[0], barras[i].gds[j].Pnom[1], barras[i].gds[j].Pnom[2], barras[i].gds[j].Qnom[0], barras[i].gds[j].Qnom[1], barras[i].gds[j].Qnom[2], barras[i].gds[j].controle);
                if (barras[i].gds[j].controle != 0)
                    fprintf(arquivo, "\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf", barras[i].gds[j].Qmin, barras[i].gds[j].Qmax, barras[i].gds[j].Vset[0], barras[i].gds[j].Vset[1], barras[i].gds[j].Vset[2], barras[i].gds[j].controlePV);
            }
        }
        fprintf(arquivo, "\n");
    }
    fprintf(arquivo, "\nLISTA DE RAMOS:\n");
    for (i = 0; i < numeroRamos; i++)
    {
        switch (ramos[i].tipo)
        {
        case ramal:
            fprintf(arquivo, "Linha %d (%d)\t%d (%d):\t%s\t%d", ramos[i].DE, ramos[i].k, ramos[i].PARA, ramos[i].m, charFases(ramos[i].fases), ramos[i].estado);
            fprintf(arquivo, "\n\t%.4lf  %.4lf\t\t%.4lf  %.4lf\t\t%.4lf  %.4lf", __real__ ramos[i].linha.Zaa, __imag__ ramos[i].linha.Zaa, __real__ ramos[i].linha.Zab, __imag__ ramos[i].linha.Zab, __real__ ramos[i].linha.Zac, __imag__ ramos[i].linha.Zac);
            fprintf(arquivo, "\n\t\t\t\t\t%.4lf  %.4lf\t\t%.4lf  %.4lf", __real__ ramos[i].linha.Zbb, __imag__ ramos[i].linha.Zbb, __real__ ramos[i].linha.Zbc, __imag__ ramos[i].linha.Zbc);
            fprintf(arquivo, "\n\t\t\t\t\t\t\t\t\t%.4lf  %.4lf", __real__ ramos[i].linha.Zcc, __imag__ ramos[i].linha.Zcc);
            break;
        case trafo:
            fprintf(arquivo, "Trafo %d (%d)\t%d (%d):\t%s\t%d", ramos[i].DE, ramos[i].k, ramos[i].PARA, ramos[i].m, charFases(ramos[i].fases), ramos[i].estado);
            fprintf(arquivo, "\n\t%.4lf  %.4lf\t%.3lf\t%.2lf / %.2lf\t%s%s%d\t%.4lf\t%.4lf", ramos[i].trafo.R, ramos[i].trafo.X, ramos[i].trafo.Snominal, ramos[i].trafo.Vpri, ramos[i].trafo.Vsec, charLigacao(ramos[i].trafo.lig_pri), charLigacao(ramos[i].trafo.lig_sec), ramos[i].trafo.defasamento, ramos[i].trafo.tap_pri, ramos[i].trafo.tap_sec);
            break;
        case regulador:
            fprintf(arquivo, "Regulador %d (%d)\t%d (%d):\t%s\t%d", ramos[i].DE, ramos[i].k, ramos[i].PARA, ramos[i].m, charFases(ramos[i].fases), ramos[i].estado);
            fprintf(arquivo, "\n\t%.4lf  %.4lf\t%.3lf\t%.2lf\t%s\t%d\t%.4lf\t%.4lf\t%.4lf", ramos[i].regulador.R, ramos[i].regulador.X, ramos[i].regulador.Snominal, ramos[i].regulador.Vnom, charLigacao(ramos[i].regulador.lig), ramos[i].regulador.controle, ramos[i].regulador.tap[0], ramos[i].regulador.tap[1], ramos[i].regulador.tap[2]);
            break;
        case chave:
            fprintf(arquivo, "Chave %d (%d)\t%d (%d):\t%s\t%d", ramos[i].DE, ramos[i].k, ramos[i].PARA, ramos[i].m, charFases(ramos[i].fases), ramos[i].estado);
            break;
        }
        fprintf(arquivo, "\n");
    }

    fprintf(arquivo, "\nLISTA DE MEDIDORES:\n");
    for (i = 0; i < nmed; i++)
    {
        fprintf(arquivo, "Medidor %d %s %d (%d)\t%d (%d):\t%s\t%d", i, charMedidor(medidas[i].tipo), medidas[i].DE, medidas[i].k, medidas[i].PARA, medidas[i].m, charFases(medidas[i].fases), medidas[i].ligado);
        fprintf(arquivo, "\n\tNumero de variaveis: %d\n", medidas[i].nvar);
        for (j = 0; j < medidas[i].nvar; j++)
        {
            fprintf(arquivo, "\n\t\t%.1lf", medidas[i].reguaH[j]);
        }
        fprintf(arquivo, "\n");
    }

    fclose(arquivo);
}

// Imprime arquivo com elementos lidos
void salvaMedidasRedeEletrica(DMED *medidas, long int **numeroMedidas)
{
    long int i, j;
    FILE *arquivo;

    // Leitura dos dados de barras
    arquivo = fopen("medidasRedeEletrica.dad", "w");
    if (arquivo == NULL)
    {
        printf("Erro ao abrir arquivo dadosRedeEletrica.dad !!!\n");
        exit(1);
    }

    fclose(arquivo);
}

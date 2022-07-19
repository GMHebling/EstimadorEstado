#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <complex.h>

#include "data_structures.h"
#include "funcoesMatematicas.h"

long int maxIT = 10000;

//------------------------------------------------------------------------------
//
// Funções com Matrizes esparas reais
//
//------------------------------------------------------------------------------

//Insere mais um elemento na lista encadeada esparsa
void sparseAdd(SPARSE **A, long int i, long int j, double *val)
{
    SPARSE *novoVertice = NULL;
    SPARSE *aux = NULL;

    novoVertice = (SPARSE *)malloc(sizeof(SPARSE));

    if (novoVertice == NULL)
    {
        printf("erro insere_SPARSE\n");
        exit(EXIT_FAILURE);
    }

    novoVertice->i = i;
    novoVertice->j = j;
    novoVertice->ij = val;
    novoVertice->x = val[0];
    novoVertice->prox = NULL;

    if (*A == NULL)
        *A = novoVertice;
    else
    {
        aux = *A;
        while (aux->prox != NULL)
            aux = aux->prox;
        aux->prox = novoVertice;
    }
}
void sparseIncCholesky(SPARSE *A, SPARSE *L)
{
    SPARSE *tmp = A;

    while (tmp != NULL)
    {
        tmp = tmp->prox;
    }
}
void sparseJacobiPrecon(SPARSE *A, SPARSE **M)
{
    SPARSE *tmp = A;

    while (tmp != NULL)
    {
        if (tmp->i == tmp->j)
            sparseAdd(M, tmp->i, tmp->i, tmp->ij);
        tmp = tmp->prox;
    }
}

void sparseFree(SPARSE *A)
{
    SPARSE *tmp = A;

    while (tmp != NULL)
    {
        tmp = A->prox;
        free(A);
        A = tmp;
    }
    A = NULL;
}
void sparsePrint(SPARSE *A)
{
    SPARSE *atual = A;
    while (atual != NULL)
    {
        printf("\n[%d , %d]: %.5lf", atual->i, atual->j, atual->x);
        atual = atual->prox;
    }
}
void sparseTiraCol(SPARSE *A, long int n)
{
    SPARSE *atual = A, *anterior = NULL, *next = NULL;
    while (atual != NULL)
    {
        if (atual->j > n)
        {
            atual->j = atual->j - 1;
            anterior = atual;
            atual = atual->prox;
        }
        else if (atual->j == n)
        {
            next = atual->prox;
            free(atual);
            if (anterior != NULL)
                anterior->prox = next;
            atual = next;
        }
        else
        {
            anterior = atual;
            atual = atual->prox;
        }
    }
}
void sparseTiraLin(SPARSE *A, long int n)
{
    SPARSE *atual = A, *anterior = NULL, *next = NULL;
    while (atual != NULL)
    {
        if (atual->i > n)
        {
            atual->j = atual->j - 1;
            anterior = atual;
            atual = atual->prox;
        }
        else if (atual->i == n)
        {
            next = atual->prox;
            free(atual);
            if (anterior != NULL)
                anterior->prox = next;
            atual = next;
        }
        else
        {
            anterior = atual;
            atual = atual->prox;
        }
    }
}
//Multiplicação da transposta de matriz esparsa por vetor
void sparseMultMatTransVet(SPARSE *A, double *b, double *c)
{
    SPARSE *next;
    int i, j;

    next = A;
    if (next != NULL)
    {
        j = next->i;
        i = next->j;
        c[i] += *next->ij * b[j];
        sparseMultMatTransVet(next->prox, b, c);
    }
}
//Multiplicação matriz esparsa por vetor
void sparseMultMatVet(SPARSE *A, double *b, double *c)
{
    SPARSE *next;
    int i, j;

    next = A;
    if (next != NULL)
    {
        i = next->i;
        j = next->j;
        c[i] += *next->ij * b[j];
        sparseMultMatVet(next->prox, b, c);
    }
}
//Matriz esparsa transposta
void sparseTrans(SPARSE *A)
{
    SPARSE *next;
    int i, j;

    next = A;
    if (next != NULL)
    {
        i = next->i;
        j = next->j;
        next->i = j;
        next->j = i;
        sparseTrans(next->prox);
    }
}
//Soluciona Lx = b e salva a resposta no vetor b (L no formato compressed column)
void sparseForward(CCSPARSE *L, double *b, int n)
{
    int i, j;
    SPARSE *tmp;

    for (i = 0; i < n; i++)
    {
        b[i] = b[i] / L[i].val->x;

        tmp = L[i].val;
        tmp = tmp->prox;
        while (tmp != NULL)
        {
            b[tmp->i] = b[tmp->i] - tmp->x * b[i];
            tmp = tmp->prox;
        }
    }
}
//Soluciona Lx = b e salva a resposta no vetor b (L no formato compressed column)
void sparseBackward(CCSPARSE *L, double *b, int n)
{
    int i, j;
    double Lii;
    SPARSE *tmp;

    for (i = (n - 1); i >= 0; i--)
    {
        tmp = L[i].val;
        Lii = L[i].val->x;
        tmp = tmp->prox;
        while (tmp != NULL)
        {
            b[i] = b[i] - tmp->x * b[tmp->i];
            tmp = tmp->prox;
        }
        b[i] = b[i] / Lii;
    }
}

//Solução de sistema linear esparso via gradientes conjugados
//recebe ponto inicial x0
//http://www2.math.ethz.ch/education/bachelor/seminars/hs2010/wave/mueller
//https://www-users.cs.umn.edu/~saad/IterMethBook_2ndEd.pdf
void sparseCG(SPARSE *A, double *b, double *x, long int n)
{
    int i, j;
    double r[n], aux[n], p[n], Dx[n], aux2[n];
    double alpha, beta;
    double direc;
    double tol = 0.0000001;
    //Computa r0 = b - A*x0 e p0 = r0
    for (i = 0; i < n; i++)
        aux[i] = 0;
    sparseMultMatVet(A, x, aux);
    for (i = 0; i < n; i++)
    {
        r[i] = b[i] - aux[i];
        p[i] = r[i];
    }

    //Loop do método
    for (j = 0; j < maxIT; j++)
    {
        for (i = 0; i < n; i++)
            aux[i] = 0;
        sparseMultMatVet(A, p, aux);
        double alpha1 = prod_escalar(p, r, n);
        alpha = alpha1 / prod_escalar(p, aux, n);

        for (i = 0; i < n; i++)
        {
            Dx[i] = alpha * p[i];
            x[i] += Dx[i];
            r[i] = r[i] - alpha * aux[i];
        }
        beta = prod_escalar(aux, r, n) / prod_escalar(aux, p, n);
        for (i = 0; i < n; i++)
        {
            p[i] = r[i] - beta * p[i];
        }

        //Convergência com tolerância numérica para o Newton CG
        //if (norma_inf(Dx,n) < tol) break;
        if (norma_inf(Dx, n) < tol)
        {
            break;
        }
        direc = -prod_escalar(b, x, n);

        if ((direc <= -0.0001) && (j > 10000))
            break;
        /*else{
            
            
        }*/
    }
    printf("\nIteracoees CG: %d \t%.4lf ", j, direc);
}

//Gradiente Conjugado Precondicionado Cholesky
void sparseCG_preCholesky(CCSPARSE *L, SPARSE *A, double *b, double *x, long int n)
{
    int i, j;
    double Mr[n], r[n], aux[n], p[n], Dx[n], aux2[n];
    double alpha, beta;
    double direc;
    double tol = 0.0000001;
    //Computa r0 = b - A*x0 e p0 = r0
    for (i = 0; i < n; i++)
        aux[i] = 0;
    sparseMultMatVet(A, x, aux);
    for (i = 0; i < n; i++)
    {
        r[i] = b[i] - aux[i];
        Mr[i] = r[i];
    }
    sparseForward(L, Mr, n);
    sparseBackward(L, Mr, n);
    for (i = 0; i < n; i++)
        p[i] = Mr[i];

    //Loop do método
    for (j = 0; j < maxIT; j++)
    {
        for (i = 0; i < n; i++)
            aux[i] = 0;
        sparseMultMatVet(A, p, aux);
        double alpha1 = prod_escalar(r, Mr, n);
        alpha = alpha1 / prod_escalar(p, aux, n);

        for (i = 0; i < n; i++)
        {
            Dx[i] = alpha * p[i];
            x[i] += Dx[i];
            r[i] = r[i] - alpha * aux[i];
        }
        for (i = 0; i < n; i++)
            Mr[i] = r[i];
        sparseForward(L, Mr, n);
        sparseBackward(L, Mr, n);

        beta = prod_escalar(r, Mr, n) / alpha1;
        for (i = 0; i < n; i++)
        {
            p[i] = Mr[i] + beta * p[i];
        }

        //Convergência com tolerância numérica para o Newton CG
        //if (norma_inf(Dx,n) < tol) break;
        if (norma_inf(Dx, n) < tol)
        {
            break;
        }
        direc = -prod_escalar(b, x, n);

        if ((direc <= -0.0001) && (j > 10000))
            break;
        /*else{
            
            
        }*/
    }
    printf("\nIteracoees CG: %d \t%.4lf ", j, direc);
}

//Gradiente conjugado para equaçao normal com precondicionador a esquerda (trinagular CCSPARSE)
void sparsePCGLS(CCSPARSE *L, SPARSE *A, double *b, double *x, long int m, long int n)
{
    int i, j;
    double z[n], rt[n], p[n], Dx[n]; //,grad[n];
    double r[m], w[m], aux[m];
    double alpha, beta;
    double tol = 1E-14; //0.00000000000000001;

    //Computa r0 = b - A*x0, rt=Atr0 e p0 = M-1*r0
    for (i = 0; i < m; i++)
        w[i] = 0;
    sparseMultMatVet(A, x, w);
    for (i = 0; i < m; i++)
        r[i] = b[i] - w[i];

    for (i = 0; i < n; i++)
        rt[i] = 0;
    sparseMultMatTransVet(A, b, rt); //rt = At*r0

    for (i = 0; i < n; i++)
        z[i] = rt[i];
    sparseForward(L, z, n);
    sparseBackward(L, z, n); // z = M^-1*rt
    for (i = 0; i < n; i++)
        p[i] = z[i];

    //    for (i=0;i<n;i++) grad[i] = 0;
    //    sparseMultMatTransVet(A, b,grad);
    //
    //Loop do método
    for (j = 0; j < maxIT; j++)
    {
        for (i = 0; i < m; i++)
            w[i] = 0;
        sparseMultMatVet(A, p, w);              //w = A*p
        double alpha1 = prod_escalar(z, rt, n); //(z . rt)
        alpha = alpha1 / prod_escalar(w, w, m); //apha = (z.rt)/||w||^2

        for (i = 0; i < n; i++)
        {
            Dx[i] = alpha * p[i];
            x[i] += Dx[i]; //x = x + alpha*p
            rt[i] = 0;
        }
        for (i = 0; i < m; i++)
            r[i] = r[i] - alpha * w[i];  //r = r - alpha*w
        sparseMultMatTransVet(A, r, rt); //rt = At*r0

        for (i = 0; i < n; i++)
            z[i] = rt[i];
        sparseForward(L, z, n);
        sparseBackward(L, z, n); // z = M-1*rt

        beta = prod_escalar(z, rt, n) / alpha1;
        for (i = 0; i < n; i++)
            p[i] = z[i] + beta * p[i]; //p = z + beta*p

        //Convergência com tolerância numérica para o Newton CG
        //if (norma_inf(Dx,n) < tol) break;
        if (norma_inf(Dx, n) < tol)
        {
            break;
        }
        //        double direc = -prod_escalar(x,grad,n);
        //        printf("%lf\n",direc);

        //        if ((direc <= -0.0001)||(j>10000)) break;
    }
    printf("\nIteracoees PCGLS: %d ", j);
}

//Gradiente conjugado para equaçao normal
void sparseCGLS(SPARSE *A, SPARSE *At, double *b, double *x, long int m, long int n)
{
    int i, j;
    double r[n], p[n], Dx[n];
    double d[m], t[m];
    double alpha, beta;
    double tol = 0.00000000000000001;

    //Computa r0 = b - A*x0 e p0 = r0
    for (i = 0; i < n; i++)
        r[i] = 0;
    sparseMultMatVet(At, b, r);

    for (i = 0; i < n; i++)
        p[i] = r[i];

    for (i = 0; i < m; i++)
    {
        d[i] = b[i];
        t[i] = 0;
    }
    sparseMultMatVet(A, p, t);

    //Loop do método
    for (j = 1; j < maxIT; j++)
    {
        double alpha1 = prod_escalar(r, r, n);
        alpha = alpha1 / prod_escalar(t, t, m);

        for (i = 0; i < n; i++)
        {
            Dx[i] = alpha * p[i];
            x[i] += Dx[i];
        }
        for (i = 0; i < m; i++)
        {
            d[i] = d[i] - alpha * t[i];
        }
        for (i = 0; i < n; i++)
            r[i] = 0;
        sparseMultMatVet(At, d, r);
        beta = prod_escalar(r, r, n) / alpha1;
        for (i = 0; i < n; i++)
        {
            p[i] = r[i] + beta * p[i];
        }
        for (i = 0; i < m; i++)
            t[i] = 0;
        sparseMultMatVet(A, p, t);
        //Convergência com tolerância numérica para o Newton CG
        //if (norma_inf(Dx,n) < tol) break;
        if (norma_inf(Dx, n) < tol)
        {
            break;
        }
    }

    printf("\nIteracoees CGLS: %d ", j);
}

//Gradiente conjugado para equaçao normal - précondicionado - Jacobi preconditioner
void sparseCGLS_preJacobi(double *M, SPARSE *A, SPARSE *At, double *b, double *x, long int m, long int n)
{
    int i, j;
    double Mr[n], r[n], p[n], Dx[n];
    double d[m], t[m];
    double alpha, beta;
    double tol = 0.00000000000001;

    //Computa r0 = b - A*x0 e p0 = r0
    for (i = 0; i < n; i++)
        r[i] = 0;
    sparseMultMatVet(At, b, r);
    for (i = 0; i < n; i++)
        p[i] = M[i] * r[i];

    //    for (i=0;i<n;i++) p[i] = r[i];

    for (i = 0; i < m; i++)
    {
        d[i] = b[i];
        t[i] = 0;
    }
    sparseMultMatVet(A, p, t);

    //Loop do método
    for (j = 1; j < maxIT; j++)
    {
        for (i = 0; i < n; i++)
            Mr[i] = M[i] * r[i];
        double alpha1 = prod_escalar(r, Mr, n);
        alpha = alpha1 / prod_escalar(t, t, m);

        for (i = 0; i < n; i++)
        {
            Dx[i] = alpha * p[i];
            x[i] += Dx[i];
        }
        for (i = 0; i < m; i++)
        {
            d[i] = d[i] - alpha * t[i];
        }
        for (i = 0; i < n; i++)
            r[i] = 0;
        sparseMultMatVet(At, d, r);
        for (i = 0; i < n; i++)
            Mr[i] = M[i] * r[i];

        beta = prod_escalar(r, Mr, n) / alpha1;
        for (i = 0; i < n; i++)
        {
            p[i] = Mr[i] + beta * p[i];
        }
        for (i = 0; i < m; i++)
            t[i] = 0;
        sparseMultMatVet(A, p, t);
        //Convergência com tolerância numérica para o Newton CG
        //if (norma_inf(Dx,n) < tol) break;
        if (norma_inf(Dx, n) < tol)
        {
            break;
        }
    }

    printf("\nIteracoees CGLS preJacobi: %d ", j);
}

//Gradiente conjugado para equaçao normal - précondicionado - Cholesky preconditioner
void sparseCGLS_preCholesky(CCSPARSE *L, SPARSE *A, SPARSE *At, double *b, double *x, long int m, long int n)
{
    int i, j;
    double Mr[n], r[n], p[n], Dx[n];
    double d[m], t[m];
    double alpha, beta;
    double tol = 0.00000000000001;

    //Computa r0 = b - A*x0 e p0 = r0
    for (i = 0; i < n; i++)
        r[i] = 0;
    sparseMultMatVet(At, b, r);
    for (i = 0; i < n; i++)
        p[i] = r[i];
    sparseForward(L, p, n);
    sparseBackward(L, p, n);
    for (i = 0; i < n; i++)
        Mr[i] = p[i];

    for (i = 0; i < m; i++)
    {
        d[i] = b[i];
        t[i] = 0;
    }
    sparseMultMatVet(A, p, t);

    //Loop do método
    for (j = 1; j < maxIT; j++)
    {
        double alpha1 = prod_escalar(r, Mr, n);
        alpha = alpha1 / prod_escalar(t, t, m);

        for (i = 0; i < n; i++)
        {
            Dx[i] = alpha * p[i];
            x[i] += Dx[i];
        }
        for (i = 0; i < m; i++)
        {
            d[i] = d[i] - alpha * t[i];
        }
        for (i = 0; i < n; i++)
            r[i] = 0;
        sparseMultMatVet(At, d, r);
        for (i = 0; i < n; i++)
            Mr[i] = r[i];
        sparseForward(L, Mr, n);
        sparseBackward(L, Mr, n);

        beta = prod_escalar(r, Mr, n) / alpha1;
        for (i = 0; i < n; i++)
        {
            p[i] = Mr[i] + beta * p[i];
        }
        for (i = 0; i < m; i++)
            t[i] = 0;
        sparseMultMatVet(A, p, t);
        //Convergência com tolerância numérica para o Newton CG
        //if (norma_inf(Dx,n) < tol) break;
        if (norma_inf(Dx, n) < tol)
        {
            break;
        }
    }

    printf("\nIteracoees CGLS preCholesky: %d ", j);
}

//Solução de sistema linear esparso via gradientes conjugados quadrático
//recebe ponto inicial x0
//http://utminers.utep.edu/xzeng/2017spring_math5330/MATH_5330_Computational_Methods_of_Linear_Algebra_files/ln07.pdf
void sparseCGS(SPARSE *A, double *b, double *x, long int n)
{
    int i, j;
    double r[n], rs[n], aux[n], aux2[n], p[n], q[n], u[n], Dx[n], tmp[n];
    double alpha, beta, tol = 0.00001;

    //Computa r0 = b - A*x0 e p0 = r0
    for (i = 0; i < n; i++)
        aux[i] = 0;
    sparseMultMatVet(A, x, aux);
    for (i = 0; i < n; i++)
    {
        r[i] = b[i] - aux[i];
        p[i] = r[i];
        u[i] = r[i];
        rs[i] = r[i] + 1;
    }

    //Loop do método
    for (j = 0; j < maxIT; j++)
    {
        for (i = 0; i < n; i++)
            aux[i] = 0;
        sparseMultMatVet(A, p, aux);
        alpha = prod_escalar(r, rs, n) / prod_escalar(aux, rs, n);
        double alpha1 = prod_escalar(r, rs, n);

        for (i = 0; i < n; i++)
        {
            q[i] = u[i] - alpha * aux[i];
            Dx[i] = alpha * (u[i] + q[i]);
            x[i] += Dx[i];
            tmp[i] = u[i] + q[i];
        }
        for (i = 0; i < n; i++)
            aux2[i] = 0;
        sparseMultMatVet(A, tmp, aux2);
        for (i = 0; i < n; i++)
        {
            r[i] = r[i] - alpha * aux2[i];
        }
        beta = prod_escalar(r, rs, n) / alpha1;
        for (i = 0; i < n; i++)
        {
            u[i] = r[i] + beta * q[i];
            p[i] = u[i] + beta * (q[i] + beta * p[i]);
        }
        //Convergência com tolerância numérica
        if (norma_inf(Dx, n) < tol)
            break;
    }

    printf("\nIteracoees CGS: %d", j);
}

//Le uma matriz no formato compressed column de arquivo
void sparseCCread(FILE *arquivo, CCSPARSE *L)
{
    int i, j;
    double val, *aux = NULL;

    while ((fscanf(arquivo, "%d %d %lf\n", &i, &j, &val)) != EOF)
    {
        L[j].n++;
        aux = &val;
        sparseAdd(&L[j].val, i, j, aux);
    }
    aux = NULL;
}
//Fatoração Cholesky da matriz A (SPD) salva a Lower como compressed column
void sparseCCcholesky(CCSPARSE *A, CCSPARSE *L, int n)
{
}

//------------------------------------------------------------------------------
//
// Funções com Matrizes complexas
//
//------------------------------------------------------------------------------
//Aloca Vetor
__complex__ double *c_vetAloca(int n)
{
    int j;
    __complex__ double *m;

    m = (__complex__ double *)malloc(n * sizeof(__complex__ double));

    for (j = 0; j < n; j++)
    {
        m[j] = 0; //Inicializa com 0.
    }

    return m;
}
//Aloca Matriz
__complex__ double **c_matAloca(int n)
{
    int i, j;
    __complex__ double **m = (__complex__ double **)malloc(n * sizeof(__complex__ double *));

    for (i = 0; i < n; i++)
    {
        m[i] = (__complex__ double *)malloc(n * sizeof(__complex__ double));
        for (j = 0; j < n; j++)
        {
            m[i][j] = 0; //Inicializa com 0.
        }
    }
    return m;
}
//Conjugado de matrizes
void c_matConj(__complex__ double **A, int n)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            A[i][j] = conj(A[i][j]);
        }
    }
}

//Igualdade de matrizes
void c_matIgual(__complex__ double **A, __complex__ double **B, int n)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            A[i][j] = B[i][j];
        }
    }
}
//Inversa da matriz método Decomposição LU (Para calcular a matriz admitância de linhas)
void c_matInversaZ(__complex__ double **A, int n)
{
    int i, j, k;
    __complex__ double **Zaux, piv;

    //Aloca Matriz
    Zaux = c_matAloca(3);
    Zaux[0][0] = 1;
    Zaux[1][1] = 1;
    Zaux[2][2] = 1;

    //Etapa Forward
    for (k = 0; k < 3 - 1; k++)
    {
        for (i = (k + 1); i < 3; i++)
        {
            if (A[k][k] != (0 + I * 0) && A[i][k] != (0 + I * 0))
            {
                piv = -((A[i][k]) / (A[k][k]));
                //printf("%lf %lf\n", real(piv), imag(piv));
                //getchar();
                for (j = 0; j < 3; j++)
                {
                    A[i][j] = (A[i][j] + (piv * A[k][j]));
                    Zaux[i][j] = (Zaux[i][j] + (piv * Zaux[k][j]));
                }
            }
        }
    }
    //Etapa Diagonalização
    for (i = 0; i < 3; i++)
    {
        if (A[i][i] != (0 + I * 0))
        {
            piv = (1 / (A[i][i]));
            // printf("%lf %lf\n", real(piv), imag(piv));
            for (j = 0; j < 3; j++)
            {
                A[i][j] = (piv * A[i][j]);
                Zaux[i][j] = (piv * Zaux[i][j]);
            }
        }
        else
        {
            Zaux[i][i] = (0 + I * 0);
        }
    }
    //Etapa Backward
    for (k = 0; k < 3 - 1; k++)
    {
        for (i = (k + 1); i < 3; i++)
        {
            piv = -(A[k][i]);
            //printf("%lf %lf\n", real(piv), imag(piv));
            for (j = 0; j < 3; j++)
            {
                A[k][j] = (A[k][j] + (piv * A[i][j]));
                Zaux[k][j] = (Zaux[k][j] + (piv * Zaux[i][j]));
            }
        }
    }
    c_matIgual(A, Zaux, 3);
}

//Transposta
void c_matTransp(__complex__ double **A, int n)
{
    int i, j;
    __complex__ double Aux[n][n];

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            Aux[i][j] = A[i][j];
        }
    }
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            A[i][j] = Aux[j][i];
        }
    }
}
//multiplica por escalar
void c_matMultEsc(__complex__ double **A, __complex__ double b, int n)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            A[i][j] = b * A[i][j];
        }
    }
}

////multiplica matrizes quadradas e retorna em A
void c_multMatMat(__complex__ double **A, __complex__ double **B, int n)
{
    int i, j, k;
    __complex__ double aux, laux[n];

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            aux = 0;
            for (k = 0; k < n; k++)
            {
                aux += A[i][k] * B[k][j];
            }
            laux[j] = aux;
        }
        for (j = 0; j < n; j++)
        {
            A[i][j] = laux[j];
        }
    }
}
////multiplica matrizes quadradas e retorna em A
void c_multMatVet(__complex__ double **A, __complex__ double *B, int n)
{
    int i, j, k;
    __complex__ double aux, laux[n];

    for (i = 0; i < n; i++)
    {
        aux = 0;
        for (j = 0; j < n; j++)
        {
            aux += A[i][j] * B[j];
        }
        A[i][0] = aux;
    }
}

//Imprime Matriz na Tela
void c_matImprime(__complex__ double **A, int n)
{
    int i, j;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("%.3lf +i %.3lf \t\t", __real__ A[i][j], __imag__ A[i][j]);
        }
        printf("\n");
    }
}

//------------------------------------------------------------------------------
//
// Funções com Matrizes Reais
//
//------------------------------------------------------------------------------
//ALoca matriz de double
double **aloca_matriz(int m, int n)
{
    int i, j;
    double **A;

    A = (double **)malloc(m * sizeof(double *));
    for (i = 0; i < m; i++)
        A[i] = (double *)malloc(n * sizeof(double));

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            A[i][j] = 0;
        }
    }
    return (A);
}
//Aloca vetor de double
double *aloca_vetor(int m)
{
    int i;
    double *A;

    A = (double *)malloc(m * sizeof(double));

    for (i = 0; i < m; i++)
    {
        A[i] = 0;
    }
    return (A);
}
//Função que retorna o produto escalar (interno) de dois vetores c = <a,b>
double prod_escalar(double *a, double *b, int n)
{
    int i;
    double c;

    c = 0;
    for (i = 0; i < n; i++)
    {
        c += a[i] * b[i];
    }
    return (c);
}
//Iguala Matrizes B = A
long int mat_ig_sparse(double ***A, int m, int n, long int *i_sparse, long int *j_sparse, double *x_sparse)
{
    int i, j;
    int index = 0;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (*A[i][j] != 0)
            {
                i_sparse[index] = i;
                j_sparse[index] = j;
                x_sparse[index] = *A[i][j];
                index += 1;
            }
        }
    }
    return index;
}
void mat_ig(double ***A, int m, int n, double **B)
{
    int i, j;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            B[i][j] = *A[i][j];
        }
    }
}
//Função que retorna a norma infinita de um vetor
//c = max( abs(a))
double norma_inf(double *a, int n)
{
    int i;
    double c;

    c = cabs(a[0]);
    for (i = 0; i < n; i++)
    {
        if (cabs(a[i]) > c)
            c = cabs(a[i]);
    }

    return (c);
}
//Função que retorna a norma euclidiana de um vetor
//c = sum( abs(a))
double norma_euc(double *a, int n)
{
    int i;
    double c = 0;

    for (i = 0; i < n; i++)
    {
        c = c + pow(a[i], 2);
    }
    c = sqrt(c);

    return (c);
}
void matTransp(double **A, int m, int n, double **At)
{
    int i, j;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            At[j][i] = A[i][j];
        }
    }
}

long int tira_refs_sparse(DMED *medidas, double ***A, int m, int n, int col1, int col2, long int *sparse_i, long int *sparse_j, double *sparse_x, double *regua, double *x, long int it)
{
    int i, j, n_cols;
    int index_i = 0;
    int index_j = 0;
    long int nzeros = 0;

    n_cols = col2 - col1 + 1;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            if (j < col1)
            {
                if (it == 0)
                {
                    regua[j] = regua[j];
                    x[j] = x[j];
                }
                if (*A[i][j] != 0)
                {
                    sparse_i[index_i] = i;
                    sparse_j[index_i] = j;
                    sparse_x[index_i] = *A[i][j];
                    index_i += 1;
                    nzeros += 1;
                }
            }
            else if (j > col2)
            {
                if (it == 0)
                {
                    regua[(j - n_cols)] = regua[j];
                    x[(j - n_cols)] = x[j];
                }
                if (*A[i][j] != 0)
                {
                    sparse_i[index_i] = i;
                    sparse_x[index_i] = *A[i][j];
                    sparse_j[index_i] = j;
                    index_i += 1;
                    nzeros += 1;
                }
            }
        }
    }

    return nzeros;
}
void tira_refs_regua(int m, int col1, int col2, double *regua)
{
    int i, j, n_cols;
    n_cols = col2 - col1 + 1;
    for (j = 0; j < m; j++)
    {
        if (j < col1)
        {
            regua[j] = regua[j];
        }
        else if (j > col2)
        {
            regua[(j - n_cols)] = regua[j];
        }
    }
}

//Tira coluna de matriz da regua e do vetor x
void tira_refs(double ***A, int m, int n, int col1, int col2, double ***temp, double *regua, double *x, long int it)
{
    int i, j, n_cols;
    int index = 0;

    n_cols = col2 - col1 + 1;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            if (j < col1)
            {
                temp[i][j] = A[i][j];
            }

            else if (j > col2)
                temp[i][(j - n_cols)] = A[i][j];
        }
    }
    if (it == 0)
    {
        for (j = 0; j < m; j++)
        {
            if (j < col1)
            {
                regua[j] = regua[j];
                x[j] = x[j];
            }
            else if (j > col2)
            {
                regua[(j - n_cols)] = regua[j];
                x[(j - n_cols)] = x[j];
            }
        }
    }
}


void tira_varsFP(double ***A, int m, int n,double *regua_rem,int nrem ,double ***temp, double *regua, double *x)
{
    int i, j, n_cols;
    int index = 0;
    int k,l;

    n_cols=m;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {

            temp[i][j] = A[i][j];
        }
    }

   

    for (i=0;i<nrem;i++)
    {
        for(j=0;j<n_cols;j++)
        {
            if(regua_rem[i]==regua[j])
            {
                for(k=j;k<n_cols-1;k++)
                {
                    regua[k]=regua[k+1];
                    x[k] = x[k+1];
                    for(l=0;l<n;l++)
                    {
                        temp[l][k] = temp[l][k+1];
                    }
                }
                n_cols--;
                break;
            }

        }
    }


}

void tira_refs2(double ***A, int m, int n, int col1, int col2, double ***temp, double *regua, double *x, long int it)
{
    int i, j, n_cols;
    int index = 0;

    n_cols = col2 - col1 + 1;
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            if (j < col1)
            {
                temp[i][j] = A[i][j];
            }

            else if (j > col2)
                temp[i][(j - n_cols)] = A[i][j];
        }
    }
    if (it == 1)
    {
        for (j = 0; j < m; j++)
        {
            if (j < col1)
            {
                regua[j] = regua[j];
                x[j] = x[j];
            }
            else if (j > col2)
            {
                regua[(j - n_cols)] = regua[j];
                x[(j - n_cols)] = x[j];
            }
        }
    }
}
//Concatena horizontalmente duas matrizes
void cat_hor(double **A, int m1, int n1, double **B, int m2, int n2, double **temp)
{
    int i, j;

    for (j = 0; j < m1; j++)
    {
        for (i = 0; i < n1; i++)
        {
            temp[j][i] = A[j][i];
        }
        for (i = 0; i < n2; i++)
        {
            temp[j][i + n1] = B[j][i];
        }
    }
}
//Fun��o que multiplica Matrizes
//C = A*B  sizeA[m1 x n1]   sizeB[m2 x n2]
void mult_mat(double **A, int m1, int n1, double **B, int m2, int n2, double **C)
{
    int i, j, k;
    double aux_sum;

    for (i = 0; i < m1; i++)
    {
        for (j = 0; j < n2; j++)
        {
            aux_sum = 0;
            for (k = 0; k < m2; k++)
                aux_sum = aux_sum + A[i][k] * B[k][j];

            //if (fabs(aux_sum) <precs) C[i][j] = 0;
            //            else
            C[i][j] = aux_sum;
        }
    }
}
void multMatVet(double **A, double *b, long int m, long int n, double *c)
{
    int i, j;
    for (i = 0; i < m; i++)
    {
        c[i] = 0;
        for (j = 0; j < n; j++)
        {
            if (A[i][j] != 0)
                c[i] += A[i][j] * b[j];
        }
    }
}
void somaMatMat(double **A, double **B, long int m, long int n)
{
    int i, j;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            A[i][j] = A[i][j] + B[i][j];
        }
    }
}
//Concatena horizontalmente duas matrizes
void cat_vert(double **A, int m1, int n1, double **B, int m2, int n2, double **temp)
{
    int i, j;

    for (j = 0; j < m1; j++)
    {
        for (i = 0; i < n1; i++)
        {
            temp[j][i] = A[j][i];
        }
    }
    for (j = 0; j < m2; j++)
    {
        for (i = 0; i < n2; i++)
        {
            temp[j + m1][i] = B[j][i];
        }
    }
}
//Concatena horizontalmente dois vetores
void cat_vert_vet(double *A, int m1, double *B, int m2, double *temp)
{
    int j;

    for (j = 0; j < m1; j++)
    {
        temp[j] = A[j];
    }
    for (j = 0; j < m2; j++)
    {
        temp[j + m1] = B[j];
    }
}

//------------------------------------------------------------------------------
//
// Funções de Algebra Linear
//
//------------------------------------------------------------------------------
//Retorna o maior autovalor da matriz A caculado pelo método das potências
double eigenvalue_largest(double **A, long int m, long int n)
{
    int i, j;
    double eigen = 0, nDx = 1, tol = 0.00001;
    double *x, *x0, *Dx;

    x = aloca_vetor(m);
    x0 = aloca_vetor(m);
    Dx = aloca_vetor(m);
    for (i = 0; i < m; i++)
        x0[i] = 1;
    while (nDx > tol)
    {
        double nx0 = norma_euc(x0, m);

        for (i = 0; i < m; i++)
        {
            x[i] = 0;
            for (j = 0; j < n; j++)
            {
                if (cabs(A[i][j]) >= 0.00000001)
                    x[i] += A[i][j] * x0[i] / nx0;
            }
        }
        for (i = 0; i < m; i++)
        {
            Dx[i] = x[i] - x0[i];
            x0[i] = x[i];
        }
        nDx = norma_inf(Dx, m);
    }

    for (i = 0; i < m; i++)
    {
        x0[i] = 0;
        for (j = 0; j < n; j++)
        {
            if (cabs(A[i][j]) >= 0.00000001)
                x0[i] += A[i][j] * x[i];
        }
    }
    eigen = prod_escalar(x0, x, m) / prod_escalar(x, x, m);

    return (eigen);
}
//Realiza a fatoração QR e retorna a matriz R
void QRfactorization(double **A, int m, int n, double **R)
{
    int i, j, k;
    double *x, *v, *col, **HA;
    double d, w, f, sum;

    x = (double *)calloc(n, sizeof(double));
    v = (double *)calloc(m, sizeof(double));
    col = (double *)calloc(m, sizeof(double));
    HA = (double **)calloc(m, sizeof(double *));
    for (i = 0; i < m; i++)
    {
        HA[i] = (double *)calloc(n, sizeof(double));
    }

    //Inicialização de variáveis
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            HA[i][j] = A[i][j];
        }
    }
    //Aplicação das transformações de Householder
    for (j = 0; j < n; j++)
    {
        for (i = j; i < m; i++)
        {
            col[i - j] = HA[i][j];
        }
        if (HA[j][j] > 0)
            d = -sqrt(prod_escalar(col, col, m - j));
        else
            d = sqrt(prod_escalar(col, col, m - j));

        w = HA[j][j] - d;
        f = sqrt(-2 * w * d);
        for (i = j; i < m; i++)
        {
            v[i - j] = HA[i][j] / f;
        }
        v[0] = w / f;
        for (i = j + 1; i < m; i++)
        {
            HA[i][j] = 0;
        }
        HA[j][j] = d;

        for (k = j + 1; k < n; k++)
        {
            for (i = j; i < m; i++)
            {
                col[i - j] = HA[i][k];
            }
            f = 2 * prod_escalar(v, col, m - j);
            for (i = j; i < m; i++)
            {
                HA[i][k] = HA[i][k] - f * v[i - j];
            }
        }
    }
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            R[i][j] = HA[i][j];
        }
    }
}

//Efetua a backward substituion em matriz triangular superior e retorna o vetor x no lugar de b
void backwardSubs(double **A, long int m, long int n, double *b)
{
    int i, j;
    double sum;

    //Substituição Backward
    for (i = (n - 1); i >= 0; i--)
    {
        sum = 0;
        for (j = (n - 1); j > i; j--)
        {
            sum = sum + A[i][j] * b[j];
        }
        b[i] = (b[i] - sum) / A[i][i];
    }
}
//Efetua a forward substituion em matriz triangular inferior e retorna o vetor x no lugar de b
void forwardSubs(double **A, long int m, long int n, double *b)
{
    int i, j;
    double sum;

    //Substituição Forward
    for (i = 0; i < m; i++)
    {
        sum = 0;
        for (j = 0; j < i; j++)
        {
            sum = sum + A[i][j] * b[j];
        }
        b[i] = (b[i] - sum) / A[i][i];
    }
}

//Householder refelctions com precondicionador Jacobi
double *solve_PreconHouseholder(double **A, int m, int n, double *b)
{
    int i, j, k;
    double *x, *v, *col, **HA, *Hb;
    double d, w, f, sum;

    double **Precon;
    Precon = aloca_matriz(m, n);
    for (i = 0; i < m; i++)
        Precon[i][i] = 1 / sqrt(A[i][i]);

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            A[i][j] = A[i][j] * Precon[i][i];
        }
    }

    x = (double *)calloc(m, sizeof(double));
    v = (double *)calloc(m, sizeof(double));
    col = (double *)calloc(m, sizeof(double));
    Hb = (double *)calloc(m, sizeof(double));
    HA = (double **)calloc(m, sizeof(double *));
    for (i = 0; i < m; i++)
    {
        HA[i] = (double *)calloc(n, sizeof(double));
    }

    //Inicialização de variáveis
    for (i = 0; i < m; i++)
    {
        Hb[i] = b[i];
        for (j = 0; j < n; j++)
        {
            HA[i][j] = A[i][j];
        }
    }
    //Aplicação das transformações de Householder
    for (j = 0; j < n - 1; j++)
    {
        for (i = j; i < m; i++)
        {
            col[i - j] = HA[i][j];
        }
        if (HA[j][j] > 0)
            d = -sqrt(prod_escalar(col, col, m - j));
        else
            d = sqrt(prod_escalar(col, col, m - j));

        w = HA[j][j] - d;
        f = sqrt(-2 * w * d);
        for (i = j; i < m; i++)
        {
            v[i - j] = HA[i][j] / f;
        }
        v[0] = w / f;
        for (i = j + 1; i < m; i++)
        {
            HA[i][j] = 0;
        }
        HA[j][j] = d;

        for (k = j + 1; k < n; k++)
        {
            for (i = j; i < m; i++)
            {
                col[i - j] = HA[i][k];
            }
            f = 2 * prod_escalar(v, col, m - j);
            for (i = j; i < m; i++)
            {
                HA[i][k] = HA[i][k] - f * v[i - j];
            }
        }
        //Aplica as tranformações de Householder no vetor b
        for (i = j; i < m; i++)
        {
            col[i - j] = Hb[i];
        }
        f = 2 * prod_escalar(v, col, m - j);
        for (i = j; i < m; i++)
        {
            Hb[i] = Hb[i] - f * v[i - j];
        }
    }
    //Substituição Backward
    for (i = (n - 1); i >= 0; i--)
    {
        sum = 0;
        for (j = (n - 1); j > i; j--)
        {
            sum = sum + HA[i][j] * x[j];
        }
        x[i] = (Hb[i] - sum) / HA[i][i];
    }
    /*for (i=0;i<m;i++)
    {
        for (j=0;j<n;j++)
        {
            printf("%lf\t", HA[i][j]);
        }
        printf("\n");
    }*/

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            A[i][j] = A[i][j] / Precon[i][i];
        }
    }
    for (j = 0; j < n; j++)
    {
        x[j] = x[j] / Precon[j][j];
    }

    return (x);
}

double *solve_Householder_LS(double **A, int m, int n, double *b)
{
    int i, j, k;
    double *x, *v, *col, **HA, *Hb;
    double d, w, f, sum;

    x = (double *)calloc(n, sizeof(double));
    v = (double *)calloc(m, sizeof(double));
    col = (double *)calloc(m, sizeof(double));
    Hb = (double *)calloc(m, sizeof(double));
    HA = (double **)calloc(m, sizeof(double *));
    for (i = 0; i < m; i++)
    {
        HA[i] = (double *)calloc(n, sizeof(double));
    }

    //Inicialização de variáveis
    for (i = 0; i < m; i++)
    {
        Hb[i] = b[i];
        for (j = 0; j < n; j++)
        {
            HA[i][j] = A[i][j];
        }
    }
    //Aplicação das transformações de Householder
    for (j = 0; j < n; j++)
    {
        for (i = j; i < m; i++)
        {
            col[i - j] = HA[i][j];
        }
        if (HA[j][j] > 0)
            d = -sqrt(prod_escalar(col, col, m - j));
        else
            d = sqrt(prod_escalar(col, col, m - j));

        w = HA[j][j] - d;
        f = sqrt(-2 * w * d);
        for (i = j; i < m; i++)
        {
            v[i - j] = HA[i][j] / f;
        }
        v[0] = w / f;
        for (i = j + 1; i < m; i++)
        {
            HA[i][j] = 0;
        }
        HA[j][j] = d;

        for (k = j + 1; k < n; k++)
        {
            for (i = j; i < m; i++)
            {
                col[i - j] = HA[i][k];
            }
            f = 2 * prod_escalar(v, col, m - j);
            for (i = j; i < m; i++)
            {
                HA[i][k] = HA[i][k] - f * v[i - j];
            }
        }
        //Aplica as tranformações de Householder no vetor b
        for (i = j; i < m; i++)
        {
            col[i - j] = Hb[i];
        }
        f = 2 * prod_escalar(v, col, m - j);
        for (i = j; i < m; i++)
        {
            Hb[i] = Hb[i] - f * v[i - j];
        }
    }
    //Substituição Backward
    for (i = (n - 1); i >= 0; i--)
    {
        sum = 0;
        for (j = (n - 1); j > i; j--)
        {
            sum = sum + HA[i][j] * x[j];
        }
        x[i] = (Hb[i] - sum) / HA[i][i];
    }
    /*for (i=0;i<m;i++){
        for (j=0;j<n;j++){
            printf("%lf\t", HA[i][j]);
        }
        printf("\n");
    }*/
    free(v);
    free(col);
    free(Hb);
    for (i = 0; i < m; i++)
        free(HA[i]);
    free(HA);
    return (x);
}

double *solve_Householder(double **A, int m, int n, double *b)
{
    int i, j, k;
    double *x, *v, *col, **HA, *Hb;
    double d, w, f, sum;

    x = (double *)calloc(m, sizeof(double));
    v = (double *)calloc(m, sizeof(double));
    col = (double *)calloc(m, sizeof(double));
    Hb = (double *)calloc(m, sizeof(double));
    HA = (double **)calloc(m, sizeof(double *));
    for (i = 0; i < m; i++)
    {
        HA[i] = (double *)calloc(n, sizeof(double));
    }

    //Inicialização de variáveis
    for (i = 0; i < m; i++)
    {
        Hb[i] = b[i];
        for (j = 0; j < n; j++)
        {
            HA[i][j] = A[i][j];
        }
    }
    //Aplicação das transformações de Householder
    for (j = 0; j < n - 1; j++)
    {
        for (i = j; i < m; i++)
        {
            col[i - j] = HA[i][j];
        }
        if (HA[j][j] > 0)
            d = -sqrt(prod_escalar(col, col, m - j));
        else
            d = sqrt(prod_escalar(col, col, m - j));

        w = HA[j][j] - d;
        f = sqrt(-2 * w * d);
        for (i = j; i < m; i++)
        {
            v[i - j] = HA[i][j] / f;
        }
        v[0] = w / f;
        for (i = j + 1; i < m; i++)
        {
            HA[i][j] = 0;
        }
        HA[j][j] = d;

        for (k = j + 1; k < n; k++)
        {
            for (i = j; i < m; i++)
            {
                col[i - j] = HA[i][k];
            }
            f = 2 * prod_escalar(v, col, m - j);
            for (i = j; i < m; i++)
            {
                HA[i][k] = HA[i][k] - f * v[i - j];
            }
        }
        //Aplica as tranformações de Householder no vetor b
        for (i = j; i < m; i++)
        {
            col[i - j] = Hb[i];
        }
        f = 2 * prod_escalar(v, col, m - j);
        for (i = j; i < m; i++)
        {
            Hb[i] = Hb[i] - f * v[i - j];
        }
    }
    //Substituição Backward
    for (i = (n - 1); i >= 0; i--)
    {
        sum = 0;
        for (j = (n - 1); j > i; j--)
        {
            sum = sum + HA[i][j] * x[j];
        }
        x[i] = (Hb[i] - sum) / HA[i][i];
    }
    return (x);
}

//Algoritmo de Crout
double *solve_Crout(double **A, int m, int n, double *b)
{
    int i, j, k;
    double **F, **L, **U, *x, *y;
    double sum;

    x = (double *)calloc(m, sizeof(double));
    y = (double *)calloc(m, sizeof(double));
    F = (double **)calloc(m, sizeof(double *));
    L = (double **)calloc(m, sizeof(double *));
    U = (double **)calloc(m, sizeof(double *));
    for (i = 0; i < m; i++)
    {
        F[i] = (double *)calloc(n, sizeof(double));
        L[i] = (double *)calloc(n, sizeof(double));
        U[i] = (double *)calloc(n, sizeof(double));
    }

    //----------------------------------------------------------
    //Fatoração pelo algoritmo de Crout
    for (j = 0; j < n; j++)
    {
        F[j][0] = A[j][0];
    }
    for (j = 0; j < n; j++)
    {
        F[0][j] = A[0][j] / A[0][0];
    }

    for (j = 0; j < n; j++)
    {
        for (k = j; k < n; k++)
        {
            sum = 0;
            for (i = 0; i < j; i++)
            {
                sum = sum + F[k][i] * F[i][j];
            }
            F[k][j] = A[k][j] - sum;
            if (fabs(F[k][j]) < EPS)
                F[k][j] = 0;
        }
        if (j != n)
        {
            for (k = (j + 1); k < n; k++)
            {
                sum = 0;
                for (i = 0; i < j; i++)
                {
                    sum = sum + F[j][i] * F[i][k];
                }
                F[j][k] = (A[j][k] - sum) / F[j][j];
                if (fabs(F[j][k]) < EPS)
                    F[j][k] = 0;
            }
        }
    }
    //Fim da Fatoração
    //----------------------------------------------------------

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (j <= i)
                L[i][j] = F[i][j];
            else
                U[i][j] = F[i][j];
        }
        U[i][i] = 1;
    }

    //---------------------------------------------------
    //Resolve Ly = b
    for (i = 0; i < n; i++)
    {
        sum = 0;
        for (j = 0; j < i; j++)
        {
            sum = sum + L[i][j] * y[j];
        }
        if (L[i][i] != 0)
            y[i] = (b[i] - sum) / L[i][i];
    }
    //---------------------------------------------------
    //Resolve Ux = y
    for (i = (n - 1); i >= 0; i--)
    {
        sum = 0;
        for (j = (n - 1); j > i; j--)
        {
            sum = sum + U[i][j] * x[j];
        }
        x[i] = (y[i] - sum);
    }

    return (x);
}


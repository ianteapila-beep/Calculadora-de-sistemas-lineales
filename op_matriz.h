#ifndef OP_MATRIZ_H
#define OP_MATRIZ_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <gmp.h>

#define MAX 1024

typedef struct
{
    int r;
    int c;
    mpq_t *V;
} Matriz;

typedef struct
{
    Matriz *A;
    Matriz *V;
} Sistema;

// Entrada segura y utilidades
void limpiabuffer(void);
void pausa(void);
int vint(void);
double vdouble(void);
char vchar(void);

// Gestión de memoria, indexación y submatrices
int I(int r, int c, Matriz X);
Matriz *crea_MD(int r, int c);
void libera_matriz(Matriz *V);
Matriz *copia_matriz(Matriz B);
Matriz *sis_matriz(Matriz A, Matriz V);
Matriz *extraer_submatriz(Matriz A, int r0, int c0, int filas, int cols);
void insertar_submatriz(Matriz *dest, Matriz src, int r0, int c0);

// Operaciones elementales y pivoteo
void inter_ren(int des, int ori, Matriz *A);
void inter_col(int des, int ori, Matriz *A);
void div_ren(int r, int col_inicio, mpq_t divisor, Matriz *A);
void comb_ren(int r_dest, int r_piv, int col_inicio, mpq_t factor, Matriz *A);
int pivoteo_parcial_idx(Matriz *A, int k, int *cam, int *fila_piv);
int pivoteo_parcial(Matriz *A, int k, int *cam);
int pivoteo_total(Matriz *A, int k, int lim_col, int *fila_piv, int *col_piv);

// Propiedades de las matrices
int cuadrada(Matriz A);
int simetrica(Matriz A);
int vect_null(Matriz A);
mpq_t *determinante(Matriz B);
mpq_t *traza(Matriz A);
int rango_matriz(Matriz A);
int es_diagonalmente_dominante(Matriz A);
int hacer_diagonalmente_dominante(Matriz *M);

// Operaciones con matrices y sistemas triangulares
Matriz *suma_matriz(Matriz A, Matriz B);
Matriz *resta_matriz(Matriz A, Matriz B);
Matriz *producto_matriz(Matriz A, Matriz B);
Matriz *producto_e_matriz(mpq_t k, Matriz A);
Matriz *transpuesta(Matriz A);
Matriz *potencia_matriz(Matriz A, int exp);
Matriz *inversa(Matriz A);
Matriz diag_sup(Matriz A, int *cam);
Matriz diag_inf(Matriz A, int *cam);
Matriz *permutar_vector(Matriz b, const int *p);
Matriz *sustitucion_adelante(Matriz L, Matriz b, int diag_unitaria);
Matriz *sustitucion_atras(Matriz U, Matriz y, int diag_unitaria);

// Lectura e impresión alineada
void imprime_matriz(Matriz *V);
void imprime_matriz_dec(Matriz *V, int dec);
void imprime_sistema(Sistema S);
void imprime_sistema_dec(Sistema S, int dec);
void imprime_aumentada(Matriz *M, int col_sep, int dec, char modo_vista);
void leer_decimal_o_fraccion(char *tok, mpq_t V);
void ingresa_fraccion(char *f, Matriz *V, int r);
void crea_matriz(Matriz *V);
void crea_matriz_rand(Matriz *V);

// Normas de matrices y vectores
mpq_t *norma_1(Matriz A);
mpq_t *norma_inf(Matriz A);
double norma_euclidiana(Matriz A);
double norma_p(Matriz A, double p);
double norma_espectral(Matriz A);

#endif
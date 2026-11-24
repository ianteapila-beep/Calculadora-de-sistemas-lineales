#ifndef METODOS_H
#define METODOS_H

#include "op_matriz.h"

typedef struct
{
    Matriz *L;
    Matriz *U;
    int *p;
} LU;

typedef struct
{
    Matriz *sol;
    int iteraciones;
    int convergio;
    double error_final;
} ResultadoIterativo;

// 1. Factorizaciones LU (Doolittle, Doolittle con Pivoteo, Crout) y Cholesky
void libera_lu(LU *lu);
LU *doolittle(Matriz A);
Matriz *resuelve_lu(Sistema S);

LU *doolittle_pivoteo(Matriz A);
Matriz *resuelve_lu_pivoteo(Sistema S);

LU *crout(Matriz A);
Matriz *resuelve_crout(Sistema S);

Matriz *ejecutar_cholesky_res(Sistema S, int dec);

// 2. Eliminación de Gauss, Gauss-Jordan y Particionado
Matriz *ejecutar_gauss_res(Sistema S, int tipo_pivoteo, int dec, char pasos, char modo_vista);
Matriz *ejecutar_gauss_jordan_res(Sistema S, int tipo_pivoteo, int dec, char pasos, char modo_vista);
Matriz *ejecutar_gauss_particionado_res(Sistema S, int tam_bloque, int dec, char pasos);

// 3. Inversa por Gauss-Jordan con pasos
Matriz *ejecutar_inversa_pasos_res(Matriz A, int dec, char pasos);

// 4. Métodos Iterativos (Jacobi, Gauss-Seidel, SOR)
ResultadoIterativo ejecutar_jacobi_res(Sistema S, int dec, double tol, int max_iter);
ResultadoIterativo ejecutar_gauss_seidel_res(Sistema S, int dec, double tol, int max_iter);
ResultadoIterativo ejecutar_sor_res(Sistema S, double omega, int dec, double tol, int max_iter);

// 5. Análisis Matricial (Número de Condición y Eigenvalores)
void ejecutar_numero_condicion(Matriz A, int dec);
void ejecutar_metodo_potencias(Matriz A, int inverso, int dec, double tol, int max_iter);

#endif
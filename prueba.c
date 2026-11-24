#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include "op_matriz.h"
#include "metodos.h"

#define LIMPIAR_PANTALLA() do { printf("\033[H\033[J"); fflush(stdout); } while (0)
#define Sal putchar('\n')

#define RES_SISTEMA 1
#define RES_MATRIZ  2

// ============================================================================
// ESTRUCTURAS PARA BANCO DE MATRICES E HISTORIAL DE RESULTADOS
// ============================================================================

typedef struct
{
    int id;
    char nombre[64];
    char origen[64];
    Sistema sis;
    int tiene_vector;
} EntradaMatriz;

typedef struct
{
    EntradaMatriz *items;
    int count;
    int capacity;
    int next_id;
    int activo_idx;
} BancoMatrices;

typedef struct
{
    int id;
    char metodo[96];
    char sistema_nombre[64];
    int tipo;               // RES_SISTEMA o RES_MATRIZ
    Matriz *A_orig;         // Copia de A usada en el calculo
    Matriz *b_orig;         // Copia de b usada (si tipo == RES_SISTEMA)
    Matriz *resultado;      // Vector solucion x (n x 1) o matriz resultado (r x c)
    int es_exacto;          // 1 = exacto (GMP), 0 = aproximado/iterativo
    int iteraciones;        // Iteraciones realizadas
    int convergio;          // 1 = convergio o exacto, 0 = no convergio
    double error_iter;      // Error final de parada
    double tiempo_ms;       // Tiempo de ejecucion en ms
} ResultadoGuardado;

typedef struct
{
    ResultadoGuardado *items;
    int count;
    int capacity;
    int next_id;
} HistorialResultados;

// ============================================================================
// UTILIDADES DE ENTRADA Y GESTIÓN DEL BANCO DE MATRICES
// ============================================================================

static void leer_linea_texto(char *dest, size_t max_len)
{
    if (fgets(dest, (int)max_len, stdin) != NULL)
    {
        dest[strcspn(dest, "\r\n")] = '\0';
    }
    else
    {
        dest[0] = '\0';
    }
}

static void iniciar_banco(BancoMatrices *banco)
{
    banco->capacity = 8;
    banco->count = 0;
    banco->next_id = 1;
    banco->activo_idx = 0;
    banco->items = (EntradaMatriz *)malloc(sizeof(EntradaMatriz) * banco->capacity);
    if (banco->items == NULL)
    {
        printf("NO HAY MEMORIA DISPONIBLE\n");
        exit(1);
    }
}

static int banco_agregar(BancoMatrices *banco, const char *nombre, const char *origen,
                         Matriz *A, Matriz *b, int tiene_vector, int activar)
{
    if (banco->count >= banco->capacity)
    {
        banco->capacity *= 2;
        EntradaMatriz *nuevo = (EntradaMatriz *)realloc(banco->items, sizeof(EntradaMatriz) * banco->capacity);
        if (nuevo == NULL)
        {
            printf("NO HAY MEMORIA DISPONIBLE\n");
            exit(1);
        }
        banco->items = nuevo;
    }

    int idx = banco->count++;
    EntradaMatriz *e = &banco->items[idx];
    e->id = banco->next_id++;

    if (nombre != NULL && nombre[0] != '\0')
        snprintf(e->nombre, sizeof(e->nombre), "%s", nombre);
    else
        snprintf(e->nombre, sizeof(e->nombre), "%s_%d", tiene_vector ? "Sistema" : "Matriz", e->id);

    snprintf(e->origen, sizeof(e->origen), "%s", (origen != NULL) ? origen : "Manual");
    e->sis.A = A;
    e->sis.V = (b != NULL) ? b : crea_MD(A->r, 1);
    e->tiene_vector = tiene_vector;

    if (activar)
        banco->activo_idx = idx;

    return e->id;
}

static void banco_eliminar(BancoMatrices *banco, int idx)
{
    if (idx < 0 || idx >= banco->count) return;

    libera_matriz(banco->items[idx].sis.A);
    libera_matriz(banco->items[idx].sis.V);

    for (int i = idx; i < banco->count - 1; ++i)
        banco->items[i] = banco->items[i + 1];

    banco->count--;
    if (banco->activo_idx >= banco->count)
        banco->activo_idx = (banco->count > 0) ? banco->count - 1 : 0;
}

static void liberar_banco(BancoMatrices *banco)
{
    if (banco == NULL || banco->items == NULL) return;
    for (int i = 0; i < banco->count; ++i)
    {
        libera_matriz(banco->items[i].sis.A);
        libera_matriz(banco->items[i].sis.V);
    }
    free(banco->items);
    banco->items = NULL;
    banco->count = 0;
}

static EntradaMatriz *banco_actual(BancoMatrices *banco)
{
    return &banco->items[banco->activo_idx];
}

static int banco_buscar_por_id(const BancoMatrices *banco, int id)
{
    for (int i = 0; i < banco->count; ++i)
    {
        if (banco->items[i].id == id)
            return i;
    }
    return -1;
}

static void mostrar_tabla_banco(const BancoMatrices *banco)
{
    printf("%-8s %-5s %-22s %-11s %-12s %-20s\n",
           "Activa", "ID", "Nombre", "Dimension", "Tipo", "Origen");
    printf("--------------------------------------------------------------------------------\n");
    for (int i = 0; i < banco->count; ++i)
    {
        const EntradaMatriz *e = &banco->items[i];
        char dim[32];
        snprintf(dim, sizeof(dim), "%d x %d", e->sis.A->r, e->sis.A->c);
        printf("%-8s %-5d %-22s %-11s %-12s %-20s\n",
               (i == banco->activo_idx) ? "  >>>" : "",
               e->id,
               e->nombre,
               dim,
               e->tiene_vector ? "Sistema [A|b]" : "Matriz A",
               e->origen);
    }
    Sal;
}

static void ofrecer_guardar_en_banco(BancoMatrices *banco, Matriz *M_res, const char *prefijo_nombre, const char *origen)
{
    if (M_res == NULL) return;

    printf("\nDeseas guardar este resultado en el Banco de Matrices? (s/n): ");
    char g = vchar();
    if (g == 's' || g == 'S')
    {
        char nombre[64];
        printf("Nombre para guardar (ENTER para usar '%s_%d'): ", prefijo_nombre, banco->next_id);
        leer_linea_texto(nombre, sizeof(nombre));
        if (nombre[0] == '\0')
            snprintf(nombre, sizeof(nombre), "%s_%d", prefijo_nombre, banco->next_id);

        printf("Establecer tambien como matriz activa? (s/n): ");
        char act = vchar();
        int activar = (act == 's' || act == 'S') ? 1 : 0;

        int id = banco_agregar(banco, nombre, origen, copia_matriz(*M_res), NULL, 0, activar);
        printf("Guardado en el Banco de Matrices con ID #%d ('%s').\n", id, nombre);
    }
}

// ============================================================================
// GESTIÓN DEL HISTORIAL DE RESULTADOS
// ============================================================================

static void iniciar_historial(HistorialResultados *hist)
{
    hist->capacity = 16;
    hist->count = 0;
    hist->next_id = 1;
    hist->items = (ResultadoGuardado *)malloc(sizeof(ResultadoGuardado) * hist->capacity);
    if (hist->items == NULL)
    {
        printf("NO HAY MEMORIA DISPONIBLE\n");
        exit(1);
    }
}

static int historial_agregar_sistema(HistorialResultados *hist, const char *metodo, const char *sis_nombre,
                                     Sistema S_orig, Matriz *sol, int es_exacto,
                                     int iteraciones, int convergio, double error_iter, double tiempo_ms)
{
    if (sol == NULL) return -1;

    if (hist->count >= hist->capacity)
    {
        hist->capacity *= 2;
        ResultadoGuardado *nuevo = (ResultadoGuardado *)realloc(hist->items, sizeof(ResultadoGuardado) * hist->capacity);
        if (nuevo == NULL)
        {
            printf("NO HAY MEMORIA DISPONIBLE\n");
            exit(1);
        }
        hist->items = nuevo;
    }

    ResultadoGuardado *r = &hist->items[hist->count++];
    r->id = hist->next_id++;
    snprintf(r->metodo, sizeof(r->metodo), "%s", metodo);
    snprintf(r->sistema_nombre, sizeof(r->sistema_nombre), "%s", sis_nombre);
    r->tipo = RES_SISTEMA;
    r->A_orig = copia_matriz(*(S_orig.A));
    r->b_orig = copia_matriz(*(S_orig.V));
    r->resultado = copia_matriz(*sol);
    r->es_exacto = es_exacto;
    r->iteraciones = iteraciones;
    r->convergio = convergio;
    r->error_iter = error_iter;
    r->tiempo_ms = tiempo_ms;

    return r->id;
}

static int historial_agregar_matriz(HistorialResultados *hist, const char *metodo, const char *sis_nombre,
                                    Matriz *A_orig, Matriz *M_res, double tiempo_ms)
{
    if (M_res == NULL) return -1;

    if (hist->count >= hist->capacity)
    {
        hist->capacity *= 2;
        ResultadoGuardado *nuevo = (ResultadoGuardado *)realloc(hist->items, sizeof(ResultadoGuardado) * hist->capacity);
        if (nuevo == NULL)
        {
            printf("NO HAY MEMORIA DISPONIBLE\n");
            exit(1);
        }
        hist->items = nuevo;
    }

    ResultadoGuardado *r = &hist->items[hist->count++];
    r->id = hist->next_id++;
    snprintf(r->metodo, sizeof(r->metodo), "%s", metodo);
    snprintf(r->sistema_nombre, sizeof(r->sistema_nombre), "%s", sis_nombre);
    r->tipo = RES_MATRIZ;
    r->A_orig = (A_orig != NULL) ? copia_matriz(*A_orig) : NULL;
    r->b_orig = NULL;
    r->resultado = copia_matriz(*M_res);
    r->es_exacto = 1;
    r->iteraciones = 0;
    r->convergio = 1;
    r->error_iter = 0.0;
    r->tiempo_ms = tiempo_ms;

    return r->id;
}

static void liberar_historial(HistorialResultados *hist)
{
    if (hist == NULL || hist->items == NULL) return;
    for (int i = 0; i < hist->count; ++i)
    {
        libera_matriz(hist->items[i].A_orig);
        libera_matriz(hist->items[i].b_orig);
        libera_matriz(hist->items[i].resultado);
    }
    free(hist->items);
    hist->items = NULL;
    hist->count = 0;
}

static int historial_buscar_por_id(const HistorialResultados *hist, int id)
{
    for (int i = 0; i < hist->count; ++i)
    {
        if (hist->items[i].id == id)
            return i;
    }
    return -1;
}

static void mostrar_tabla_historial(const HistorialResultados *hist)
{
    printf("%-5s %-11s %-18s %-30s %-18s %-10s\n",
           "ID", "Tipo", "Origen", "Metodo / Operacion", "Precision/Iter", "Tiempo(ms)");
    printf("----------------------------------------------------------------------------------------------------\n");
    for (int i = 0; i < hist->count; ++i)
    {
        const ResultadoGuardado *r = &hist->items[i];
        char dim_tipo[24];
        char info_iter[32];

        if (r->tipo == RES_SISTEMA)
            snprintf(dim_tipo, sizeof(dim_tipo), "Sol (%dx1)", r->resultado->r);
        else
            snprintf(dim_tipo, sizeof(dim_tipo), "Mat (%dx%d)", r->resultado->r, r->resultado->c);

        if (r->es_exacto)
            snprintf(info_iter, sizeof(info_iter), "Exacto (GMP)");
        else
            snprintf(info_iter, sizeof(info_iter), "%s (%d it)", r->convergio ? "Conv." : "NoConv", r->iteraciones);

        printf("%-5d %-11s %-18s %-30s %-18s %-10.3f\n",
               r->id, dim_tipo, r->sistema_nombre, r->metodo, info_iter, r->tiempo_ms);
    }
    Sal;
}

// ============================================================================
// MENU 1: GESTIONAR MÚLTIPLES MATRICES Y SISTEMAS
// ============================================================================

static void pedir_dimensiones(int *r, int *c, int solo_cuadrada)
{
    do
    {
        printf("Ingresa el numero de renglones: ");
        *r = vint();
        if (*r <= 1)
            printf("Dato invalido (debe ser mayor a 1).\n");
    } while (*r <= 1);

    if (solo_cuadrada)
    {
        *c = *r;
    }
    else
    {
        do
        {
            printf("Ingresa el numero de columnas: ");
            *c = vint();
            if (*c <= 0)
                printf("Dato invalido (debe ser mayor a 0).\n");
        } while (*c <= 0);
    }
}

static void agregar_matriz_a_banco(BancoMatrices *banco)
{
    LIMPIAR_PANTALLA();
    printf("AGREGAR NUEVA MATRIZ / SISTEMA AL BANCO\n\n");
    printf("1.- Matriz cuadrada (n x n) + vector b (para sistemas, det, inversa, LU)\n");
    printf("2.- Sistema rectangular (r x c) + vector b (para minimos cuadrados, rango, normas)\n");
    printf("3.- Matriz sola (r x c) (para operaciones basicas, transpuesta y normas)\n");

    int tipo = 0;
    do
    {
        printf("Selecciona tipo: ");
        tipo = vint();
        if (tipo < 1 || tipo > 3)
            printf("Opcion invalida (elige 1, 2 o 3).\n");
    } while (tipo < 1 || tipo > 3);

    char nombre[64];
    printf("\nNombre identificador (ENTER para automatico): ");
    leer_linea_texto(nombre, sizeof(nombre));

    int r = 3, c = 3;
    pedir_dimensiones(&r, &c, (tipo == 1));

    int modo = 0;
    do
    {
        printf("\nModo de llenado (0 = Aleatoria, 1 = Manual): ");
        modo = vint();
        if (modo != 0 && modo != 1)
            printf("Opcion invalida (elige 0 o 1).\n");
    } while (modo != 0 && modo != 1);

    Matriz *A = crea_MD(r, c);
    Matriz *b = crea_MD(r, 1);
    int con_vector = (tipo == 1 || tipo == 2) ? 1 : 0;

    if (modo == 0)
    {
        crea_matriz_rand(A);
        crea_matriz_rand(b);
    }
    else
    {
        printf("\nIngresa la matriz A (%d x %d) por renglon:\n", r, c);
        crea_matriz(A);
        if (con_vector)
        {
            printf("\nIngresa el vector b (%d x 1):\n", r);
            crea_matriz(b);
        }
    }

    int id = banco_agregar(banco, nombre, (modo == 0) ? "Aleatoria" : "Manual",
                           A, b, con_vector, 1);
    printf("\nMatriz/Sistema agregado con ID #%d y establecido como activo.\n", id);
    pausa();
}

static void menu_gestionar_matrices(BancoMatrices *banco)
{
    int op = 0;
    do
    {
        LIMPIAR_PANTALLA();
        printf("BANCO DE MATRICES Y SISTEMAS (%d guardadas)\n\n", banco->count);
        mostrar_tabla_banco(banco);

        printf("1.- Agregar nueva matriz o sistema\n");
        printf("2.- Seleccionar matriz/sistema activo\n");
        printf("3.- Ver contenido de una matriz del banco\n");
        printf("4.- Duplicar matriz/sistema activo\n");
        printf("5.- Eliminar una matriz del banco\n");
        printf("6.- Regresar al menu principal\n");
        printf("Selecciona una opcion: ");
        op = vint();

        switch (op)
        {
            case 1:
                agregar_matriz_a_banco(banco);
                break;

            case 2:
            {
                printf("Ingresa el ID de la matriz/sistema a activar: ");
                int id = vint();
                int idx = banco_buscar_por_id(banco, id);
                if (idx == -1)
                {
                    printf("No se encontro ninguna matriz con ID #%d.\n", id);
                    pausa();
                }
                else
                {
                    banco->activo_idx = idx;
                }
                break;
            }

            case 3:
            {
                printf("Ingresa el ID de la matriz/sistema a inspeccionar: ");
                int id = vint();
                int idx = banco_buscar_por_id(banco, id);
                if (idx == -1)
                {
                    printf("No se encontro ninguna matriz con ID #%d.\n", id);
                }
                else
                {
                    LIMPIAR_PANTALLA();
                    EntradaMatriz *e = &banco->items[idx];
                    printf("Detalle de #%d '%s' (%d x %d) [Origen: %s]:\n\n",
                           e->id, e->nombre, e->sis.A->r, e->sis.A->c, e->origen);
                    if (e->tiene_vector)
                        imprime_sistema(e->sis);
                    else
                        imprime_matriz(e->sis.A);
                }
                pausa();
                break;
            }

            case 4:
            {
                EntradaMatriz *act = banco_actual(banco);
                char nombre_copia[64];
                snprintf(nombre_copia, sizeof(nombre_copia), "%.57s_copia", act->nombre);
                int id = banco_agregar(banco, nombre_copia, "Duplicada",
                                       copia_matriz(*(act->sis.A)),
                                       copia_matriz(*(act->sis.V)),
                                       act->tiene_vector, 1);
                printf("Duplicada como #%d '%s' y marcada como activa.\n", id, nombre_copia);
                pausa();
                break;
            }

            case 5:
            {
                if (banco->count <= 1)
                {
                    printf("Debe permanecer al menos una matriz en el banco.\n");
                    pausa();
                    break;
                }
                printf("Ingresa el ID de la matriz a eliminar: ");
                int id = vint();
                int idx = banco_buscar_por_id(banco, id);
                if (idx == -1)
                {
                    printf("No se encontro ninguna matriz con ID #%d.\n", id);
                }
                else
                {
                    banco_eliminar(banco, idx);
                    printf("Matriz #%d eliminada.\n", id);
                }
                pausa();
                break;
            }

            case 6:
                break;

            default:
                printf("Opcion no valida\n");
                pausa();
                break;
        }
    } while (op != 6);
}

// ============================================================================
// MENÚ 2: MÉTODOS DIRECTOS Y FACTORIZACIONES PARA SISTEMAS LINEALES
// ============================================================================

static int seleccionar_pivoteo(void)
{
    int piv = 0;
    do
    {
        printf("\nTipo de pivoteo:\n");
        printf("1.- Parcial\n");
        printf("2.- Total\n");
        printf("3.- Ninguno\n");
        printf("Selecciona una opcion: ");
        piv = vint();
    } while (piv < 1 || piv > 3);
    return piv;
}

static const char *nombre_pivoteo(int piv)
{
    if (piv == 1) return "Piv. Parcial";
    if (piv == 2) return "Piv. Total";
    return "Sin Pivoteo";
}

static void menu_sistemas_directos(BancoMatrices *banco, HistorialResultados *hist, int dec, char pasos, char modo_vista)
{
    EntradaMatriz *act = banco_actual(banco);
    Sistema *M = &act->sis;

    int op = 0;
    do
    {
        LIMPIAR_PANTALLA();
        printf("METODOS DIRECTOS Y FACTORIZACIONES PARA SISTEMAS LINEALES\n\n");
        printf("Sistema activo: #%d '%s' (%d x %d):\n", act->id, act->nombre, M->A->r, M->A->c);
        imprime_sistema(*M);
        Sal;
        printf("1.- Eliminacion de Gauss\n");
        printf("2.- Metodo de Gauss-Jordan\n");
        printf("3.- Gauss-Jordan Particionado (por bloques)\n");
        printf("4.- Factorizacion LU (Doolittle sin pivoteo)\n");
        printf("5.- Factorizacion LU con Pivoteo Parcial (PA = LU)\n");
        printf("6.- Factorizacion de Crout\n");
        printf("7.- Factorizacion de Cholesky\n");
        printf("8.- Volver al menu principal\n");
        printf("Selecciona una opcion: ");
        op = vint();

        if (op >= 1 && op <= 7 && !cuadrada(*(M->A)))
        {
            printf("\nERROR: La opcion %d requiere una matriz cuadrada (actual: %d x %d).\n",
                   op, M->A->r, M->A->c);
            pausa();
            continue;
        }

        switch (op)
        {
            case 1:
            {
                int piv = seleccionar_pivoteo();
                LIMPIAR_PANTALLA();
                printf("METODO DE GAUSS (%s)\n", nombre_pivoteo(piv));

                clock_t t0 = clock();
                Matriz *sol = ejecutar_gauss_res(*M, piv, dec, pasos, modo_vista);
                double ms = (double)(clock() - t0) * 1000.0 / CLOCKS_PER_SEC;

                if (sol != NULL)
                {
                    char etiqueta[96];
                    snprintf(etiqueta, sizeof(etiqueta), "Gauss (%s)", nombre_pivoteo(piv));
                    int id_res = historial_agregar_sistema(hist, etiqueta, act->nombre, *M,
                                                           sol, 1, 0, 1, 0.0, ms);
                    printf("\n[Resultado #%d guardado en historial (%.3f ms)]\n", id_res, ms);
                    ofrecer_guardar_en_banco(banco, sol, "Sol_Gauss", etiqueta);
                    libera_matriz(sol);
                }
                pausa();
                break;
            }
            case 2:
            {
                int piv = seleccionar_pivoteo();
                LIMPIAR_PANTALLA();
                printf("METODO DE GAUSS-JORDAN (%s)\n", nombre_pivoteo(piv));

                clock_t t0 = clock();
                Matriz *sol = ejecutar_gauss_jordan_res(*M, piv, dec, pasos, modo_vista);
                double ms = (double)(clock() - t0) * 1000.0 / CLOCKS_PER_SEC;

                if (sol != NULL)
                {
                    char etiqueta[96];
                    snprintf(etiqueta, sizeof(etiqueta), "Gauss-Jordan (%s)", nombre_pivoteo(piv));
                    int id_res = historial_agregar_sistema(hist, etiqueta, act->nombre, *M,
                                                           sol, 1, 0, 1, 0.0, ms);
                    printf("\n[Resultado #%d guardado en historial (%.3f ms)]\n", id_res, ms);
                    ofrecer_guardar_en_banco(banco, sol, "Sol_GJ", etiqueta);
                    libera_matriz(sol);
                }
                pausa();
                break;
            }
            case 3:
            {
                int n = M->A->r;
                int tam_bloque = 0;
                do
                {
                    printf("\nIngresa el tamano de la particion (divisor exacto de %d, menor que %d): ", n, n);
                    tam_bloque = vint();
                    if (tam_bloque < 1 || tam_bloque >= n || n % tam_bloque != 0)
                    {
                        printf("Tamano no valido. Debe dividir exactamente a %d.\n", n);
                        tam_bloque = -1;
                    }
                } while (tam_bloque < 1);

                LIMPIAR_PANTALLA();
                printf("METODO DE GAUSS-JORDAN PARTICIONADO (bloque = %d)\n", tam_bloque);

                clock_t t0 = clock();
                Matriz *sol = ejecutar_gauss_particionado_res(*M, tam_bloque, dec, pasos);
                double ms = (double)(clock() - t0) * 1000.0 / CLOCKS_PER_SEC;

                if (sol != NULL)
                {
                    char etiqueta[96];
                    snprintf(etiqueta, sizeof(etiqueta), "GJ Particionado (b=%d)", tam_bloque);
                    int id_res = historial_agregar_sistema(hist, etiqueta, act->nombre, *M,
                                                           sol, 1, 0, 1, 0.0, ms);
                    printf("\n[Resultado #%d guardado en historial (%.3f ms)]\n", id_res, ms);
                    ofrecer_guardar_en_banco(banco, sol, "Sol_Part", etiqueta);
                    libera_matriz(sol);
                }
                pausa();
                break;
            }
            case 4:
            {
                LIMPIAR_PANTALLA();
                printf("METODO DE DOOLITTLE (LU SIN PIVOTEO)\n\n");

                clock_t t0 = clock();
                Matriz *sol = resuelve_lu(*M);
                double ms = (double)(clock() - t0) * 1000.0 / CLOCKS_PER_SEC;

                if (sol != NULL)
                {
                    printf("Solucion del sistema:\n");
                    imprime_matriz(sol);
                    int id_res = historial_agregar_sistema(hist, "Doolittle LU", act->nombre, *M,
                                                           sol, 1, 0, 1, 0.0, ms);
                    printf("\n[Resultado #%d guardado en historial (%.3f ms)]\n", id_res, ms);
                    ofrecer_guardar_en_banco(banco, sol, "Sol_LU", "Doolittle LU");
                    libera_matriz(sol);
                }
                pausa();
                break;
            }
            case 5:
            {
                LIMPIAR_PANTALLA();
                printf("FACTORIZACION LU CON PIVOTEO PARCIAL (PA = LU)\n\n");

                clock_t t0 = clock();
                Matriz *sol = resuelve_lu_pivoteo(*M);
                double ms = (double)(clock() - t0) * 1000.0 / CLOCKS_PER_SEC;

                if (sol != NULL)
                {
                    printf("Solucion del sistema:\n");
                    imprime_matriz(sol);
                    int id_res = historial_agregar_sistema(hist, "LU Pivoteo Parcial", act->nombre, *M,
                                                           sol, 1, 0, 1, 0.0, ms);
                    printf("\n[Resultado #%d guardado en historial (%.3f ms)]\n", id_res, ms);
                    ofrecer_guardar_en_banco(banco, sol, "Sol_LUPiv", "LU Pivoteo Parcial");
                    libera_matriz(sol);
                }
                pausa();
                break;
            }
            case 6:
            {
                LIMPIAR_PANTALLA();
                printf("FACTORIZACION DE CROUT (LU)\n\n");

                clock_t t0 = clock();
                Matriz *sol = resuelve_crout(*M);
                double ms = (double)(clock() - t0) * 1000.0 / CLOCKS_PER_SEC;

                if (sol != NULL)
                {
                    printf("Solucion del sistema:\n");
                    imprime_matriz(sol);
                    int id_res = historial_agregar_sistema(hist, "Crout LU", act->nombre, *M,
                                                           sol, 1, 0, 1, 0.0, ms);
                    printf("\n[Resultado #%d guardado en historial (%.3f ms)]\n", id_res, ms);
                    ofrecer_guardar_en_banco(banco, sol, "Sol_Crout", "Crout LU");
                    libera_matriz(sol);
                }
                pausa();
                break;
            }
            case 7:
            {
                LIMPIAR_PANTALLA();
                printf("FACTORIZACION DE CHOLESKY (A = L * L^T)\n\n");

                clock_t t0 = clock();
                Matriz *sol = ejecutar_cholesky_res(*M, dec);
                double ms = (double)(clock() - t0) * 1000.0 / CLOCKS_PER_SEC;

                if (sol != NULL)
                {
                    int id_res = historial_agregar_sistema(hist, "Cholesky (LL^T)", act->nombre, *M,
                                                           sol, 0, 0, 1, 0.0, ms);
                    printf("\n[Resultado #%d guardado en historial (%.3f ms)]\n", id_res, ms);
                    ofrecer_guardar_en_banco(banco, sol, "Sol_Chol", "Cholesky");
                    libera_matriz(sol);
                }
                pausa();
                break;
            }
            case 8:
                break;
            default:
                printf("Opcion no valida\n");
                pausa();
                break;
        }
    } while (op != 8);
}

// ============================================================================
// MENÚ 3: MÉTODOS ITERATIVOS (JACOBI, GAUSS-SEIDEL, SOR)
// ============================================================================

static void pedir_criterio_paro_iterativo(double *tol, int *max_iter)
{
    *tol = 1e-10;
    *max_iter = 10000;

    printf("Selecciona criterio de paro:\n");
    printf("1.- Ingresar tolerancia (Max iteraciones por defecto: 10000)\n");
    printf("2.- Ingresar maximo de iteraciones (Tolerancia por defecto: 1e-10)\n");

    int op = 0;
    do
    {
        printf("Selecciona una opcion: ");
        op = vint();
        if (op != 1 && op != 2)
            printf("Opcion invalida (elige 1 o 2).\n");
    } while (op != 1 && op != 2);

    if (op == 1)
    {
        do
        {
            printf("Tolerancia del error (ej. 1e-8): ");
            *tol = vdouble();
        } while (*tol <= 0.0);
    }
    else
    {
        do
        {
            printf("Maximo de iteraciones permitidas: ");
            *max_iter = vint();
        } while (*max_iter <= 0);
    }
}

static void menu_metodos_iterativos(BancoMatrices *banco, HistorialResultados *hist, int dec)
{
    EntradaMatriz *act = banco_actual(banco);
    Sistema *M = &act->sis;

    if (!cuadrada(*(M->A)))
    {
        printf("\nERROR: La matriz actual '%s' (%d x %d) no es cuadrada.\n",
               act->nombre, M->A->r, M->A->c);
        pausa();
        return;
    }

    int op_m = 0;
    do
    {
        LIMPIAR_PANTALLA();
        printf("METODOS ITERATIVOS PARA SISTEMAS LINEALES\n\n");
        printf("Sistema activo: #%d '%s' (%d x %d):\n", act->id, act->nombre, M->A->r, M->A->c);
        imprime_sistema(*M);
        Sal;
        printf("1.- Metodo de Jacobi\n");
        printf("2.- Metodo de Gauss-Seidel\n");
        printf("3.- Metodo SOR (Sobre-relajacion Sucesiva)\n");
        printf("4.- Volver al menu principal\n");
        printf("Selecciona una opcion: ");
        op_m = vint();

        if (op_m >= 1 && op_m <= 3)
        {
            double omega = 1.0;
            if (op_m == 3)
            {
                do
                {
                    printf("\nIngresa el factor de relajacion omega en (0, 2) [ej. 1.25]: ");
                    omega = vdouble();
                    if (omega <= 0.0 || omega >= 2.0)
                        printf("Valor invalido. Debe cumplir 0 < omega < 2.\n");
                } while (omega <= 0.0 || omega >= 2.0);
            }

            double tol = 1e-10;
            int max_iter = 10000;
            Sal;
            pedir_criterio_paro_iterativo(&tol, &max_iter);

            LIMPIAR_PANTALLA();
            clock_t t0 = clock();
            ResultadoIterativo res = {NULL, 0, 0, 0.0};
            char etiqueta[96];

            if (op_m == 1)
            {
                printf("METODO ITERATIVO DE JACOBI\n");
                res = ejecutar_jacobi_res(*M, dec, tol, max_iter);
                snprintf(etiqueta, sizeof(etiqueta), "Jacobi (tol=%.1e, %d it)", tol, res.iteraciones);
            }
            else if (op_m == 2)
            {
                printf("METODO ITERATIVO DE GAUSS-SEIDEL\n");
                res = ejecutar_gauss_seidel_res(*M, dec, tol, max_iter);
                snprintf(etiqueta, sizeof(etiqueta), "Gauss-Seidel (tol=%.1e, %d it)", tol, res.iteraciones);
            }
            else
            {
                printf("METODO ITERATIVO SOR (omega = %.4f)\n", omega);
                res = ejecutar_sor_res(*M, omega, dec, tol, max_iter);
                snprintf(etiqueta, sizeof(etiqueta), "SOR w=%.2f (%d it)", omega, res.iteraciones);
            }
            double ms = (double)(clock() - t0) * 1000.0 / CLOCKS_PER_SEC;

            if (res.sol != NULL)
            {
                int id_res = historial_agregar_sistema(hist, etiqueta, act->nombre, *M,
                                                       res.sol, 0, res.iteraciones, res.convergio,
                                                       res.error_final, ms);
                printf("\n[Resultado #%d guardado en historial (%.3f ms)]\n", id_res, ms);
                ofrecer_guardar_en_banco(banco, res.sol, "Sol_Iter", etiqueta);
                libera_matriz(res.sol);
            }
            pausa();
        }
        else if (op_m != 4)
        {
            printf("Opcion no valida\n");
            pausa();
        }
    } while (op_m != 4);
}

// ============================================================================
// MENÚ 5: COMPARACIÓN DE RESULTADOS DE MÉTODOS Y OPERACIONES
// ============================================================================

static void calcular_residuo_sistema(const Matriz *A, const Matriz *b, const Matriz *x,
                                     double *res_inf, double *res_eucl)
{
    *res_inf = -1.0;
    *res_eucl = -1.0;
    if (A == NULL || b == NULL || x == NULL) return;

    Matriz *Ax = producto_matriz(*A, *x);
    if (Ax == NULL) return;

    Matriz *diff = resta_matriz(*b, *Ax);
    libera_matriz(Ax);
    if (diff == NULL) return;

    mpq_t *n_inf = norma_inf(*diff);
    if (n_inf != NULL)
    {
        *res_inf = mpq_get_d(*n_inf);
        mpq_clear(*n_inf);
        free(n_inf);
    }
    *res_eucl = norma_euclidiana(*diff);
    libera_matriz(diff);
}

static void comparar_soluciones_sistemas(const HistorialResultados *hist, int dec)
{
    int disponibles = 0;
    for (int i = 0; i < hist->count; ++i)
    {
        if (hist->items[i].tipo == RES_SISTEMA)
            disponibles++;
    }

    if (disponibles < 2)
    {
        printf("\nSe necesitan al menos 2 soluciones de sistemas guardadas para comparar (actualmente hay %d).\n", disponibles);
        pausa();
        return;
    }

    int k = 0;
    do
    {
        printf("Cuantas soluciones deseas comparar? (2 - %d): ", disponibles);
        k = vint();
        if (k < 2 || k > disponibles)
            printf("Cantidad invalida.\n");
    } while (k < 2 || k > disponibles);

    int *sel_idx = (int *)malloc(sizeof(int) * k);
    int n_dim = -1;

    for (int s = 0; s < k; ++s)
    {
        while (1)
        {
            printf("Ingresa el ID de la solucion #%d: ", s + 1);
            int id = vint();
            int idx = historial_buscar_por_id(hist, id);

            if (idx == -1 || hist->items[idx].tipo != RES_SISTEMA)
            {
                printf("ID invalido o no corresponde a una solucion de sistema.\n");
                continue;
            }

            int repetido = 0;
            for (int prev = 0; prev < s; ++prev)
            {
                if (sel_idx[prev] == idx)
                {
                    repetido = 1;
                    break;
                }
            }
            if (repetido)
            {
                printf("Ese ID ya fue seleccionado.\n");
                continue;
            }

            if (n_dim == -1)
            {
                n_dim = hist->items[idx].resultado->r;
            }
            else if (hist->items[idx].resultado->r != n_dim)
            {
                printf("Dimension incompatible: se esperaba vector de %d incognitas y #%d tiene %d.\n",
                       n_dim, id, hist->items[idx].resultado->r);
                continue;
            }

            sel_idx[s] = idx;
            break;
        }
    }

    LIMPIAR_PANTALLA();
    printf("====================================================================================================\n");
    printf("                    COMPARACION DE RESULTADOS ENTRE METODOS (%d incognitas)\n", n_dim);
    printf("====================================================================================================\n\n");

    int ref_exacta = -1;
    for (int s = 0; s < k; ++s)
    {
        if (hist->items[sel_idx[s]].es_exacto)
        {
            ref_exacta = s;
            break;
        }
    }

    int col_w = (dec + 10 < 18) ? 18 : dec + 10;

    printf("1. VALORES OBTENIDOS POR INCOGNITA:\n\n");
    printf("%-10s", "Variable");
    for (int s = 0; s < k; ++s)
    {
        char hdr[48];
        snprintf(hdr, sizeof(hdr), "[#%d] %.16s", hist->items[sel_idx[s]].id, hist->items[sel_idx[s]].metodo);
        printf("%-*s", col_w, hdr);
    }
    if (ref_exacta != -1)
        printf("%-18s", "Fraccion Exacta");
    Sal;

    for (int s = 0; s < k + 1; ++s)
        printf("------------------");
    Sal;

    for (int r = 0; r < n_dim; ++r)
    {
        char var_lbl[24];
        snprintf(var_lbl, sizeof(var_lbl), "x_%d", r + 1);
        printf("%-10s", var_lbl);

        for (int s = 0; s < k; ++s)
        {
            Matriz *sol = hist->items[sel_idx[s]].resultado;
            double val = mpq_get_d(sol->V[I(r, 0, *sol)]);
            printf("%-*.*f", col_w, dec, val);
        }

        if (ref_exacta != -1)
        {
            Matriz *sol_ref = hist->items[sel_idx[ref_exacta]].resultado;
            gmp_printf("%Qd", sol_ref->V[I(r, 0, *sol_ref)]);
        }
        Sal;
    }
    Sal;

    printf("2. METRICAS DE EXACTITUD (RESIDUO r = b - Ax) Y RENDIMIENTO:\n\n");
    printf("%-6s %-26s %-14s %-10s %-14s %-14s %-10s\n",
           "ID", "Metodo", "Tipo", "Iter.", "||b-Ax||_inf", "||b-Ax||_2", "Tiempo(ms)");
    printf("----------------------------------------------------------------------------------------------------\n");

    for (int s = 0; s < k; ++s)
    {
        const ResultadoGuardado *rg = &hist->items[sel_idx[s]];
        double res_inf = 0.0, res_eucl = 0.0;
        calcular_residuo_sistema(rg->A_orig, rg->b_orig, rg->resultado, &res_inf, &res_eucl);

        char iter_buf[16];
        if (rg->es_exacto)
            snprintf(iter_buf, sizeof(iter_buf), "N/A");
        else
            snprintf(iter_buf, sizeof(iter_buf), "%d", rg->iteraciones);

        printf("%-6d %-26.26s %-14s %-10s %-14.3e %-14.3e %-10.3f\n",
               rg->id,
               rg->metodo,
               rg->es_exacto ? "Exacto (GMP)" : "Iterativo",
               iter_buf,
               res_inf,
               res_eucl,
               rg->tiempo_ms);
    }
    Sal;

    printf("3. DIFERENCIAS DIRECTAS ENTRE METODOS (||x_A - x_B||):\n\n");
    printf("%-28s %-18s %-14s %-14s %-14s\n",
           "Par comparado", "Igualdad Exacta", "Dif Max (inf)", "Dif Eucl (L2)", "Err Rel Max");
    printf("----------------------------------------------------------------------------------------------------\n");

    for (int a = 0; a < k; ++a)
    {
        for (int b = a + 1; b < k; ++b)
        {
            const ResultadoGuardado *ra = &hist->items[sel_idx[a]];
            const ResultadoGuardado *rb = &hist->items[sel_idx[b]];

            int identicos = 1;
            for (int r = 0; r < n_dim; ++r)
            {
                if (!mpq_equal(ra->resultado->V[I(r, 0, *(ra->resultado))],
                               rb->resultado->V[I(r, 0, *(rb->resultado))]))
                {
                    identicos = 0;
                    break;
                }
            }

            Matriz *diff = resta_matriz(*(ra->resultado), *(rb->resultado));
            double d_inf = 0.0, d_eucl = 0.0, err_rel_max = 0.0;

            if (diff != NULL)
            {
                mpq_t *ni = norma_inf(*diff);
                if (ni != NULL)
                {
                    d_inf = mpq_get_d(*ni);
                    mpq_clear(*ni);
                    free(ni);
                }
                d_eucl = norma_euclidiana(*diff);

                const Matriz *ref_vec = rb->es_exacto ? rb->resultado : ra->resultado;
                for (int r = 0; r < n_dim; ++r)
                {
                    double d_val = fabs(mpq_get_d(diff->V[I(r, 0, *diff)]));
                    double ref_val = fabs(mpq_get_d(ref_vec->V[I(r, 0, *ref_vec)]));
                    double rel = (ref_val > 1e-15) ? (d_val / ref_val) : d_val;
                    if (rel > err_rel_max)
                        err_rel_max = rel;
                }
                libera_matriz(diff);
            }

            char par_lbl[64];
            snprintf(par_lbl, sizeof(par_lbl), "#%d vs #%d", ra->id, rb->id);
            printf("%-28s %-18s %-14.3e %-14.3e %-14.3e\n",
                   par_lbl,
                   identicos ? "IDENTICOS (100%)" : "Diferentes",
                   d_inf,
                   d_eucl,
                   err_rel_max);
        }
    }

    free(sel_idx);
    pausa();
}

static void comparar_dos_matrices(const HistorialResultados *hist, int dec)
{
    if (hist->count < 2)
    {
        printf("\nSe necesitan al menos 2 resultados guardados en el historial.\n");
        pausa();
        return;
    }

    printf("Ingresa el ID del primer resultado: ");
    int id1 = vint();
    printf("Ingresa el ID del segundo resultado: ");
    int id2 = vint();

    int idx1 = historial_buscar_por_id(hist, id1);
    int idx2 = historial_buscar_por_id(hist, id2);

    if (idx1 == -1 || idx2 == -1)
    {
        printf("Uno o ambos IDs no existen en el historial.\n");
        pausa();
        return;
    }

    const ResultadoGuardado *r1 = &hist->items[idx1];
    const ResultadoGuardado *r2 = &hist->items[idx2];

    if (r1->resultado->r != r2->resultado->r || r1->resultado->c != r2->resultado->c)
    {
        printf("Dimensiones incompatibles: #%d es %dx%d y #%d es %dx%d.\n",
               id1, r1->resultado->r, r1->resultado->c,
               id2, r2->resultado->r, r2->resultado->c);
        pausa();
        return;
    }

    LIMPIAR_PANTALLA();
    printf("COMPARACION DIRECTA: #%d (%s) vs #%d (%s)\n\n",
           r1->id, r1->metodo, r2->id, r2->metodo);

    int identicas = 1;
    int total = r1->resultado->r * r1->resultado->c;
    for (int i = 0; i < total; ++i)
    {
        if (!mpq_equal(r1->resultado->V[i], r2->resultado->V[i]))
        {
            identicas = 0;
            break;
        }
    }

    printf("Resultado #%d (%s):\n", r1->id, r1->metodo);
    imprime_matriz(r1->resultado);
    Sal;

    printf("Resultado #%d (%s):\n", r2->id, r2->metodo);
    imprime_matriz(r2->resultado);
    Sal;

    if (identicas)
    {
        printf(">>> AMBOS RESULTADOS SON ALGEBRAICAMENTE IDENTICOS (Exactitud GMP) <<<\n");
    }
    else
    {
        Matriz *diff = resta_matriz(*(r1->resultado), *(r2->resultado));
        printf("Matriz Diferencia (Resultado #%d - Resultado #%d):\n", r1->id, r2->id);
        imprime_matriz(diff);
        Sal;

        mpq_t *n1 = norma_1(*diff);
        mpq_t *ninf = norma_inf(*diff);
        double neucl = norma_euclidiana(*diff);

        if (n1 != NULL)
        {
            gmp_printf("- Norma L1 de la diferencia:          %Qd (%.*f)\n", *n1, dec, mpq_get_d(*n1));
            mpq_clear(*n1);
            free(n1);
        }
        if (ninf != NULL)
        {
            gmp_printf("- Norma L_inf de la diferencia:       %Qd (%.*f)\n", *ninf, dec, mpq_get_d(*ninf));
            mpq_clear(*ninf);
            free(ninf);
        }
        printf("- Norma Euclidiana de la diferencia:  %.*f\n", dec, neucl);
        libera_matriz(diff);
    }

    pausa();
}

static void menu_comparar_resultados(HistorialResultados *hist, int dec)
{
    int op = 0;
    do
    {
        LIMPIAR_PANTALLA();
        printf("HISTORIAL Y COMPARACION DE RESULTADOS (%d guardados)\n\n", hist->count);

        if (hist->count == 0)
        {
            printf("Aun no hay resultados guardados. Ejecuta algun metodo o operacion primero.\n\n");
            pausa();
            return;
        }

        mostrar_tabla_historial(hist);

        printf("1.- Comparar soluciones de sistemas (x) entre diferentes metodos\n");
        printf("2.- Comparar dos resultados matriciales cualesquiera (M1 vs M2)\n");
        printf("3.- Ver detalle completo de un resultado guardado\n");
        printf("4.- Limpiar historial de resultados\n");
        printf("5.- Regresar al menu principal\n");
        printf("Selecciona una opcion: ");
        op = vint();

        switch (op)
        {
            case 1:
                comparar_soluciones_sistemas(hist, dec);
                break;

            case 2:
                comparar_dos_matrices(hist, dec);
                break;

            case 3:
            {
                printf("Ingresa el ID del resultado a inspeccionar: ");
                int id = vint();
                int idx = historial_buscar_por_id(hist, id);
                if (idx == -1)
                {
                    printf("No se encontro el resultado #%d.\n", id);
                }
                else
                {
                    LIMPIAR_PANTALLA();
                    const ResultadoGuardado *r = &hist->items[idx];
                    printf("Detalle Resultado #%d\n", r->id);
                    printf("- Metodo / Operacion: %s\n", r->metodo);
                    printf("- Matriz/Sistema:     %s\n", r->sistema_nombre);
                    printf("- Tiempo ejecucion:   %.3f ms\n\n", r->tiempo_ms);
                    imprime_matriz(r->resultado);
                }
                pausa();
                break;
            }

            case 4:
            {
                liberar_historial(hist);
                iniciar_historial(hist);
                printf("Historial limpiado.\n");
                pausa();
                return;
            }

            case 5:
                break;

            default:
                printf("Opcion no valida\n");
                pausa();
                break;
        }
    } while (op != 5);
}

// ============================================================================
// MENÚ 6: OPCIONES EXTRA (OPERACIONES BÁSICAS, PROPIEDADES, INVERSA, ANÁLISIS)
// ============================================================================

static void menu_operaciones_basicas(BancoMatrices *banco, HistorialResultados *hist)
{
    int op = 0;
    do
    {
        EntradaMatriz *act = banco_actual(banco);
        Sistema *M = &act->sis;

        LIMPIAR_PANTALLA();
        printf("OPERACIONES BASICAS CON MATRICES\n\n");
        printf("Matriz A activa: #%d '%s' (%d x %d):\n", act->id, act->nombre, M->A->r, M->A->c);
        imprime_matriz(M->A);
        Sal;
        printf("1.- Suma (A + B)\n");
        printf("2.- Resta (A - B)\n");
        printf("3.- Producto matricial (A * B)\n");
        printf("4.- Producto por escalar (k * A)\n");
        printf("5.- Matriz Transpuesta (A^T)\n");
        printf("6.- Potencia de Matriz (A^k)\n");
        printf("7.- Regresar\n");
        printf("Selecciona una opcion: ");
        op = vint();

        if (op >= 1 && op <= 3)
        {
            int rb = M->A->r, cb = M->A->c;
            if (op == 3)
            {
                printf("\nPara A * B, B debe tener %d renglones.\n", M->A->c);
                rb = M->A->c;
                do
                {
                    printf("Ingresa el numero de columnas de B: ");
                    cb = vint();
                } while (cb <= 0);
            }
            else
            {
                printf("\nLa matriz B tendra dimension %d x %d.\n", rb, cb);
            }

            int modo = 0;
            do
            {
                printf("Origen de B (0 = Aleatoria, 1 = Manual, 2 = Elegir del Banco de Matrices): ");
                modo = vint();
            } while (modo < 0 || modo > 2);

            Matriz *B = NULL;
            if (modo == 2)
            {
                Sal;
                mostrar_tabla_banco(banco);
                printf("Ingresa el ID de la matriz B (%d x %d): ", rb, cb);
                int id_b = vint();
                int idx_b = banco_buscar_por_id(banco, id_b);
                if (idx_b == -1)
                {
                    printf("ID no encontrado.\n");
                    pausa();
                    continue;
                }
                Matriz *cand = banco->items[idx_b].sis.A;
                if (cand->r != rb || cand->c != cb)
                {
                    printf("Dimensiones incompatibles: #%d es %d x %d (se requiere %d x %d).\n",
                           id_b, cand->r, cand->c, rb, cb);
                    pausa();
                    continue;
                }
                B = copia_matriz(*cand);
            }
            else
            {
                B = crea_MD(rb, cb);
                if (modo == 0) crea_matriz_rand(B);
                else crea_matriz(B);
            }

            LIMPIAR_PANTALLA();
            printf("Matriz A:\n");
            imprime_matriz(M->A);
            printf("\nMatriz B:\n");
            imprime_matriz(B);
            Sal;

            clock_t t0 = clock();
            Matriz *C = NULL;
            const char *op_nombre = "";
            if (op == 1) { C = suma_matriz(*(M->A), *B); op_nombre = "Suma (A + B)"; }
            else if (op == 2) { C = resta_matriz(*(M->A), *B); op_nombre = "Resta (A - B)"; }
            else if (op == 3) { C = producto_matriz(*(M->A), *B); op_nombre = "Producto (A * B)"; }
            double ms = (double)(clock() - t0) * 1000.0 / CLOCKS_PER_SEC;

            if (C != NULL)
            {
                printf("Resultado:\n");
                imprime_matriz(C);

                int id_res = historial_agregar_matriz(hist, op_nombre, act->nombre, M->A, C, ms);
                printf("\n[Resultado #%d guardado en historial]\n", id_res);
                ofrecer_guardar_en_banco(banco, C, "Res_Op", op_nombre);
                libera_matriz(C);
            }
            libera_matriz(B);
            pausa();
        }
        else if (op == 4)
        {
            char buf[128], clean_buf[128];
            mpq_t k;
            mpq_init(k);

            printf("\nIngresa el escalar k (entero, fraccion o decimal): ");
            if (fgets(buf, sizeof(buf), stdin) != NULL)
            {
                int j = 0;
                for (int i = 0; buf[i] != '\0'; ++i)
                {
                    if (!isspace((unsigned char)buf[i]))
                        clean_buf[j++] = buf[i];
                }
                clean_buf[j] = '\0';
                leer_decimal_o_fraccion(clean_buf, k);
            }

            clock_t t0 = clock();
            Matriz *C = producto_e_matriz(k, *(M->A));
            double ms = (double)(clock() - t0) * 1000.0 / CLOCKS_PER_SEC;

            LIMPIAR_PANTALLA();
            gmp_printf("Resultado de (%Qd) * A:\n", k);
            imprime_matriz(C);

            if (C != NULL)
            {
                int id_res = historial_agregar_matriz(hist, "Producto Escalar (k*A)", act->nombre, M->A, C, ms);
                printf("\n[Resultado #%d guardado en historial]\n", id_res);
                ofrecer_guardar_en_banco(banco, C, "Res_Escalar", "Escalar k*A");
                libera_matriz(C);
            }

            mpq_clear(k);
            pausa();
        }
        else if (op == 5)
        {
            clock_t t0 = clock();
            Matriz *T = transpuesta(*(M->A));
            double ms = (double)(clock() - t0) * 1000.0 / CLOCKS_PER_SEC;

            LIMPIAR_PANTALLA();
            printf("Matriz Transpuesta A^T (%d x %d):\n\n", T->r, T->c);
            imprime_matriz(T);

            int id_res = historial_agregar_matriz(hist, "Transpuesta (A^T)", act->nombre, M->A, T, ms);
            printf("\n[Resultado #%d guardado en historial (%.3f ms)]\n", id_res, ms);
            ofrecer_guardar_en_banco(banco, T, "Transpuesta", "Transpuesta A^T");
            libera_matriz(T);
            pausa();
        }
        else if (op == 6)
        {
            printf("\nIngresa el exponente entero k (puede ser negativo si A es invertible): ");
            int exp = vint();

            LIMPIAR_PANTALLA();
            clock_t t0 = clock();
            Matriz *P = potencia_matriz(*(M->A), exp);
            double ms = (double)(clock() - t0) * 1000.0 / CLOCKS_PER_SEC;

            if (P != NULL)
            {
                printf("Resultado de A^%d:\n\n", exp);
                imprime_matriz(P);

                char etiqueta[64];
                snprintf(etiqueta, sizeof(etiqueta), "Potencia (A^%d)", exp);
                int id_res = historial_agregar_matriz(hist, etiqueta, act->nombre, M->A, P, ms);
                printf("\n[Resultado #%d guardado en historial (%.3f ms)]\n", id_res, ms);
                ofrecer_guardar_en_banco(banco, P, "Potencia", etiqueta);
                libera_matriz(P);
            }
            pausa();
        }
        else if (op != 7)
        {
            printf("Opcion no valida\n");
            pausa();
        }
    } while (op != 7);
}

static void menu_propiedades_y_normas(BancoMatrices *banco, HistorialResultados *hist, int dec, char pasos)
{
    int op = 0;
    do
    {
        EntradaMatriz *act = banco_actual(banco);
        Sistema *M = &act->sis;

        LIMPIAR_PANTALLA();
        printf("PROPIEDADES, INVERSA, NORMAS Y ANALISIS ESPECTRAL\n\n");
        printf("Matriz A activa: #%d '%s' (%d x %d):\n", act->id, act->nombre, M->A->r, M->A->c);
        imprime_matriz(M->A);
        Sal;
        printf("1.- Determinante\n");
        printf("2.- Matriz Inversa (Gauss-Jordan)\n");
        printf("3.- Propiedades, Rango, Nulidad y Traza\n");
        printf("4.- Calcular Normas (L1, L_inf, Euclidiana, Espectral y Norma-p)\n");
        printf("5.- Numero de Condicion kappa(A)\n");
        printf("6.- Eigenvalor y Eigenvector (Metodo de Potencias)\n");
        printf("7.- Regresar\n");
        printf("Selecciona una opcion: ");
        op = vint();

        switch (op)
        {
            case 1:
            {
                LIMPIAR_PANTALLA();
                mpq_t *det = determinante(*(M->A));
                if (det != NULL)
                {
                    gmp_printf("Determinante exacto: %Qd\n", *det);
                    printf("Determinante decimal: %.*f\n", dec, mpq_get_d(*det));
                    mpq_clear(*det);
                    free(det);
                }
                pausa();
                break;
            }
            case 2:
            {
                LIMPIAR_PANTALLA();
                clock_t t0 = clock();
                Matriz *inv = ejecutar_inversa_pasos_res(*(M->A), dec, pasos);
                double ms = (double)(clock() - t0) * 1000.0 / CLOCKS_PER_SEC;

                if (inv != NULL)
                {
                    int id_res = historial_agregar_matriz(hist, "Inversa (Gauss-Jordan)", act->nombre, M->A, inv, ms);
                    printf("\n[Inversa guardada en historial con ID #%d (%.3f ms)]\n", id_res, ms);
                    ofrecer_guardar_en_banco(banco, inv, "Inversa", "Inversa Gauss-Jordan");
                    libera_matriz(inv);
                }
                pausa();
                break;
            }
            case 3:
            {
                LIMPIAR_PANTALLA();
                int rk = rango_matriz(*(M->A));
                int nul = M->A->c - rk;

                printf("Propiedades de la matriz A (%d x %d):\n\n", M->A->r, M->A->c);
                printf("- Es cuadrada:                 %s\n", cuadrada(*(M->A)) ? "SI" : "NO");
                printf("- Es simetrica:                %s\n", simetrica(*(M->A)) ? "SI" : "NO");
                printf("- Sin renglones/cols nulos:    %s\n", vect_null(*(M->A)) ? "SI" : "NO");
                printf("- Rango rank(A):               %d\n", rk);
                printf("- Nulidad nullity(A):          %d\n", nul);

                if (cuadrada(*(M->A)))
                {
                    mpq_t *tr = traza(*(M->A));
                    if (tr != NULL)
                    {
                        gmp_printf("- Traza tr(A):                 %Qd  (%.*f)\n", *tr, dec, mpq_get_d(*tr));
                        mpq_clear(*tr);
                        free(tr);
                    }
                }
                pausa();
                break;
            }
            case 4:
            {
                printf("\nIngresa el valor de p >= 1 para la Norma-p (ej. 3): ");
                double p = vdouble();

                LIMPIAR_PANTALLA();
                printf("Calculando normas de A (%d x %d)...\n\n", M->A->r, M->A->c);
                mpq_t *n1 = norma_1(*(M->A));
                mpq_t *ninf = norma_inf(*(M->A));
                double neucl = norma_euclidiana(*(M->A));
                double nesp = norma_espectral(*(M->A));
                double np = norma_p(*(M->A), p);

                if (n1 != NULL)
                {
                    gmp_printf("- Norma L1 (max columna):   %Qd  (%.*f)\n", *n1, dec, mpq_get_d(*n1));
                    mpq_clear(*n1);
                    free(n1);
                }
                if (ninf != NULL)
                {
                    gmp_printf("- Norma L_inf (max fila):   %Qd  (%.*f)\n", *ninf, dec, mpq_get_d(*ninf));
                    mpq_clear(*ninf);
                    free(ninf);
                }
                printf("- Norma Euclidiana (Frob.): %.*f\n", dec, neucl);
                printf("- Norma Espectral (sigma):  %.*f\n", dec, nesp);
                if (np >= 0.0)
                    printf("- Norma-p (p = %.2f):       %.*f\n", p, dec, np);

                pausa();
                break;
            }
            case 5:
            {
                LIMPIAR_PANTALLA();
                ejecutar_numero_condicion(*(M->A), dec);
                pausa();
                break;
            }
            case 6:
            {
                printf("\n1.- Eigenvalor dominante (Metodo de Potencias Directo)\n");
                printf("2.- Eigenvalor de menor magnitud (Metodo de Potencias Inverso)\n");
                printf("Selecciona variante: ");
                int var = vint();

                LIMPIAR_PANTALLA();
                ejecutar_metodo_potencias(*(M->A), (var == 2) ? 1 : 0, dec, 1e-12, 5000);
                pausa();
                break;
            }
            case 7:
                break;
            default:
                printf("Opcion no valida\n");
                pausa();
                break;
        }
    } while (op != 7);
}

static void menu_extras(BancoMatrices *banco, HistorialResultados *hist, int dec, char pasos)
{
    int op = 0;
    do
    {
        LIMPIAR_PANTALLA();
        printf("OPCIONES EXTRA:\n\n");
        printf("1.- Operaciones Basicas (Suma, Resta, Producto, Escalar, Transpuesta, Potencia)\n");
        printf("2.- Propiedades, Determinante, Inversa, Normas, Condicion y Eigenvalores\n");
        printf("3.- Regresar al menu principal\n");
        printf("\nSelecciona una opcion: ");
        op = vint();

        switch (op)
        {
            case 1:
                menu_operaciones_basicas(banco, hist);
                break;
            case 2:
                menu_propiedades_y_normas(banco, hist, dec, pasos);
                break;
            case 3:
                break;
            default:
                printf("Opcion no valida\n");
                pausa();
                break;
        }
    } while (op != 3);
}

// ============================================================================
// MENÚ 4: CONFIGURACIÓN
// ============================================================================

static void menu_configuracion(int *dec, char *pasos, char *modo_vista)
{
    LIMPIAR_PANTALLA();
    printf("CONFIGURACION DE VISUALIZACION\n\n");

    do
    {
        printf("Formato de matriz en Gauss y Gauss-Jordan (f = Fracciones, d = Punto flotante) [actual: %c]: ", *modo_vista);
        *modo_vista = (char)tolower((unsigned char)vchar());
    } while (*modo_vista != 'f' && *modo_vista != 'd');

    do
    {
        printf("Cantidad de decimales a mostrar (1 - 15) [actual: %d]: ", *dec);
        *dec = vint();
    } while (*dec < 1 || *dec > 15);

    do
    {
        printf("Mostrar pasos intermedios en Gauss/Gauss-Jordan/Inversa? (s/n) [actual: %c]: ", *pasos);
        *pasos = (char)tolower((unsigned char)vchar());
    } while (*pasos != 's' && *pasos != 'n');
}

// ============================================================================
// MAIN
// ============================================================================

int main(void)
{
    int op = 0;
    int dec = 5;
    char pasos = 's';
    char modo_vista = 'f'; // 'f' = Fracciones por defecto, 'd' = Punto flotante

    srand((unsigned int)time(NULL));

    BancoMatrices banco;
    HistorialResultados hist;
    iniciar_banco(&banco);
    iniciar_historial(&hist);

    Matriz *A_ini = crea_MD(3, 3);
    Matriz *b_ini = crea_MD(3, 1);
    crea_matriz_rand(A_ini);
    crea_matriz_rand(b_ini);
    banco_agregar(&banco, "Sistema_1", "Aleatoria", A_ini, b_ini, 1, 1);

    do
    {
        EntradaMatriz *act = banco_actual(&banco);

        LIMPIAR_PANTALLA();
        printf("CALCULADORA DE MATRICES Y SISTEMAS LINEALES\n\n");
        printf("Activa: #%d '%s' | Matrices en banco: %d | Resultados guardados: %d\n",
               act->id, act->nombre, banco.count, hist.count);

        if (cuadrada(*(act->sis.A)) || act->tiene_vector)
        {
            printf("Sistema actual (%d x %d):\n", act->sis.A->r, act->sis.A->c);
            imprime_sistema(act->sis);
        }
        else
        {
            printf("Matriz actual (%d x %d):\n", act->sis.A->r, act->sis.A->c);
            imprime_matriz(act->sis.A);
        }
        Sal;
        printf("1.- Gestionar Matrices (Agregar / Cambiar / Ver / Eliminar)\n");
        printf("2.- Resolver Sistema: Metodos Directos y Factorizaciones\n");
        printf("3.- Resolver Sistema: Metodos Iterativos (Jacobi, Gauss-Seidel, SOR)\n");
        printf("4.- Configuracion (Decimales: %d | Pasos: %c | Gauss/GJ: %s)\n",
               dec, pasos, (modo_vista == 'd') ? "Punto flotante" : "Fracciones");
        printf("5.- Comparar Resultados de Metodos (%d guardados)\n", hist.count);
        printf("6.- Extras (Operaciones Basicas, Propiedades, Inversa, Normas, Eigenvalores)\n");
        printf("7.- Salir\n");

        printf("Selecciona una opcion: ");
        op = vint();

        switch (op)
        {
            case 1:
                menu_gestionar_matrices(&banco);
                break;
            case 2:
                menu_sistemas_directos(&banco, &hist, dec, pasos, modo_vista);
                break;
            case 3:
                menu_metodos_iterativos(&banco, &hist, dec);
                break;
            case 4:
                menu_configuracion(&dec, &pasos, &modo_vista);
                break;
            case 5:
                menu_comparar_resultados(&hist, dec);
                break;
            case 6:
                menu_extras(&banco, &hist, dec, pasos);
                break;
            case 7:
                break;
            default:
                printf("Opcion no valida\n");
                pausa();
                break;
        }
    } while (op != 7);

    liberar_historial(&hist);
    liberar_banco(&banco);
    return 0;
}
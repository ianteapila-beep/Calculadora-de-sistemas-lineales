#include "metodos.h"

// ============================================================================
// 1. FACTORIZACIONES LU (DOOLITTLE, DOOLITTLE CON PIVOTEO, CROUT) Y CHOLESKY
// ============================================================================

void libera_lu(LU *lu)
{
    if (lu == NULL) return;
    if (lu->L != NULL) libera_matriz(lu->L);
    if (lu->U != NULL) libera_matriz(lu->U);
    if (lu->p != NULL) free(lu->p);
    free(lu);
}

LU *doolittle(Matriz A)
{
    if (!cuadrada(A))
    {
        printf("LA MATRIZ DEBE SER CUADRADA\n");
        return NULL;
    }

    int n = A.r;
    LU *res = (LU *)malloc(sizeof(LU));
    if (res == NULL)
    {
        printf("NO HAY MEMORIA DISPONIBLE\n");
        exit(1);
    }

    res->L = crea_MD(n, n);
    res->U = crea_MD(n, n);
    res->p = (int *)malloc(sizeof(int) * n);
    for (int k = 0; k < n; ++k)
    {
        res->p[k] = k;
        mpq_set_si(res->L->V[I(k, k, *(res->L))], 1, 1);
    }

    mpq_t suma, prod;
    mpq_inits(suma, prod, NULL);

    for (int i = 0; i < n; ++i)
    {
        for (int j = i; j < n; ++j)
        {
            mpq_set_si(suma, 0, 1);
            for (int k = 0; k < i; ++k)
            {
                mpq_mul(prod, res->L->V[I(i, k, *(res->L))], res->U->V[I(k, j, *(res->U))]);
                mpq_add(suma, suma, prod);
            }
            mpq_sub(res->U->V[I(i, j, *(res->U))], A.V[I(i, j, A)], suma);
        }

        if (mpq_sgn(res->U->V[I(i, i, *(res->U))]) == 0)
        {
            printf("ERROR: Pivote cero en U[%d][%d]. Doolittle sin pivoteo no puede continuar.\n\n", i + 1, i + 1);
            printf("Matriz L obtenida hasta el momento:\n");
            imprime_matriz(res->L);
            putchar('\n');
            printf("Matriz U obtenida hasta el momento:\n");
            imprime_matriz(res->U);
            putchar('\n');

            mpq_clears(suma, prod, NULL);
            libera_lu(res);
            return NULL;
        }

        for (int j = i + 1; j < n; ++j)
        {
            mpq_set_si(suma, 0, 1);
            for (int k = 0; k < i; ++k)
            {
                mpq_mul(prod, res->L->V[I(j, k, *(res->L))], res->U->V[I(k, i, *(res->U))]);
                mpq_add(suma, suma, prod);
            }
            mpq_sub(prod, A.V[I(j, i, A)], suma);
            mpq_div(res->L->V[I(j, i, *(res->L))], prod, res->U->V[I(i, i, *(res->U))]);
        }
    }

    mpq_clears(suma, prod, NULL);
    return res;
}

Matriz *resuelve_lu(Sistema S)
{
    if (S.A == NULL || S.V == NULL || !cuadrada(*(S.A)) || S.V->r != S.A->r || S.V->c != 1)
    {
        printf("ERROR: Dimensiones incompatibles. A debe ser n x n y b debe ser n x 1.\n");
        return NULL;
    }

    printf("Sistema original [A | b]:\n");
    imprime_sistema(S);
    putchar('\n');

    LU *lu = doolittle(*(S.A));
    if (lu == NULL)
        return NULL;

    printf("Matriz L:\n");
    imprime_matriz(lu->L);
    putchar('\n');

    printf("Matriz U:\n");
    imprime_matriz(lu->U);
    putchar('\n');

    Matriz *y = sustitucion_adelante(*(lu->L), *(S.V), 1);
    Matriz *x = (y != NULL) ? sustitucion_atras(*(lu->U), *y, 0) : NULL;

    libera_matriz(y);
    libera_lu(lu);
    return x;
}

LU *doolittle_pivoteo(Matriz A)
{
    if (!cuadrada(A))
    {
        printf("LA MATRIZ DEBE SER CUADRADA\n");
        return NULL;
    }

    int n = A.r;
    LU *res = (LU *)malloc(sizeof(LU));
    if (res == NULL)
    {
        printf("NO HAY MEMORIA DISPONIBLE\n");
        exit(1);
    }

    res->L = crea_MD(n, n);
    res->U = crea_MD(n, n);
    res->p = (int *)malloc(sizeof(int) * n);
    for (int k = 0; k < n; ++k) res->p[k] = k;

    Matriz *W = copia_matriz(A);
    mpq_t mult;
    mpq_init(mult);

    for (int k = 0; k < n; ++k)
    {
        int piv = k;
        if (!pivoteo_parcial_idx(W, k, NULL, &piv))
        {
            printf("ERROR: La matriz es singular (pivote nulo en columna %d). No tiene solucion unica.\n", k + 1);
            mpq_clear(mult);
            libera_matriz(W);
            libera_lu(res);
            return NULL;
        }

        if (piv != k)
        {
            int tmp_p = res->p[k];
            res->p[k] = res->p[piv];
            res->p[piv] = tmp_p;

            for (int c = 0; c < k; ++c)
                mpq_swap(res->L->V[I(k, c, *(res->L))], res->L->V[I(piv, c, *(res->L))]);
        }

        mpq_set_si(res->L->V[I(k, k, *(res->L))], 1, 1);

        for (int i = k + 1; i < n; ++i)
        {
            mpq_div(mult, W->V[I(i, k, *W)], W->V[I(k, k, *W)]);
            mpq_set(res->L->V[I(i, k, *(res->L))], mult);
            comb_ren(i, k, k, mult, W);
        }
    }

    for (int i = 0; i < n; ++i)
        for (int j = i; j < n; ++j)
            mpq_set(res->U->V[I(i, j, *(res->U))], W->V[I(i, j, *W)]);

    mpq_clear(mult);
    libera_matriz(W);
    return res;
}

Matriz *resuelve_lu_pivoteo(Sistema S)
{
    if (S.A == NULL || S.V == NULL || !cuadrada(*(S.A)) || S.V->r != S.A->r || S.V->c != 1)
    {
        printf("ERROR: Dimensiones incompatibles.\n");
        return NULL;
    }

    printf("Sistema original [A | b]:\n");
    imprime_sistema(S);
    putchar('\n');

    LU *lu = doolittle_pivoteo(*(S.A));
    if (lu == NULL)
        return NULL;

    int n = S.A->r;
    printf("Vector de permutacion P (orden de filas): [ ");
    for (int i = 0; i < n; ++i)
        printf("%d ", lu->p[i] + 1);
    printf("]\n\n");

    printf("Matriz L:\n");
    imprime_matriz(lu->L);
    putchar('\n');

    printf("Matriz U:\n");
    imprime_matriz(lu->U);
    putchar('\n');

    Matriz *Pb = permutar_vector(*(S.V), lu->p);
    Matriz *y = (Pb != NULL) ? sustitucion_adelante(*(lu->L), *Pb, 1) : NULL;
    Matriz *x = (y != NULL) ? sustitucion_atras(*(lu->U), *y, 0) : NULL;

    libera_matriz(Pb);
    libera_matriz(y);
    libera_lu(lu);
    return x;
}

LU *crout(Matriz A)
{
    if (!cuadrada(A))
    {
        printf("LA MATRIZ DEBE SER CUADRADA\n");
        return NULL;
    }

    int n = A.r;
    LU *res = (LU *)malloc(sizeof(LU));
    if (res == NULL)
    {
        printf("NO HAY MEMORIA DISPONIBLE\n");
        exit(1);
    }

    res->L = crea_MD(n, n);
    res->U = crea_MD(n, n);
    res->p = (int *)malloc(sizeof(int) * n);
    for (int k = 0; k < n; ++k)
    {
        res->p[k] = k;
        mpq_set_si(res->U->V[I(k, k, *(res->U))], 1, 1);
    }

    mpq_t suma, prod;
    mpq_inits(suma, prod, NULL);

    for (int j = 0; j < n; ++j)
    {
        for (int i = j; i < n; ++i)
        {
            mpq_set_si(suma, 0, 1);
            for (int k = 0; k < j; ++k)
            {
                mpq_mul(prod, res->L->V[I(i, k, *(res->L))], res->U->V[I(k, j, *(res->U))]);
                mpq_add(suma, suma, prod);
            }
            mpq_sub(res->L->V[I(i, j, *(res->L))], A.V[I(i, j, A)], suma);
        }

        if (mpq_sgn(res->L->V[I(j, j, *(res->L))]) == 0)
        {
            printf("ERROR: Pivote cero en L[%d][%d]. Crout sin pivoteo no puede continuar.\n", j + 1, j + 1);
            mpq_clears(suma, prod, NULL);
            libera_lu(res);
            return NULL;
        }

        for (int i = j + 1; i < n; ++i)
        {
            mpq_set_si(suma, 0, 1);
            for (int k = 0; k < j; ++k)
            {
                mpq_mul(prod, res->L->V[I(j, k, *(res->L))], res->U->V[I(k, i, *(res->U))]);
                mpq_add(suma, suma, prod);
            }
            mpq_sub(prod, A.V[I(j, i, A)], suma);
            mpq_div(res->U->V[I(j, i, *(res->U))], prod, res->L->V[I(j, j, *(res->L))]);
        }
    }

    mpq_clears(suma, prod, NULL);
    return res;
}

Matriz *resuelve_crout(Sistema S)
{
    if (S.A == NULL || S.V == NULL || !cuadrada(*(S.A)) || S.V->r != S.A->r || S.V->c != 1)
    {
        printf("ERROR: Dimensiones incompatibles.\n");
        return NULL;
    }

    printf("Sistema original [A | b]:\n");
    imprime_sistema(S);
    putchar('\n');

    LU *lu = crout(*(S.A));
    if (lu == NULL)
        return NULL;

    printf("Matriz L (Crout):\n");
    imprime_matriz(lu->L);
    putchar('\n');

    printf("Matriz U (Crout):\n");
    imprime_matriz(lu->U);
    putchar('\n');

    Matriz *y = sustitucion_adelante(*(lu->L), *(S.V), 0);
    Matriz *x = (y != NULL) ? sustitucion_atras(*(lu->U), *y, 1) : NULL;

    libera_matriz(y);
    libera_lu(lu);
    return x;
}

Matriz *ejecutar_cholesky_res(Sistema S, int dec)
{
    if (S.A == NULL || S.V == NULL || !cuadrada(*(S.A)) || S.V->r != S.A->r || S.V->c != 1)
    {
        printf("ERROR: El sistema debe ser cuadrado (A de n x n y b de n x 1).\n");
        return NULL;
    }

    if (!simetrica(*(S.A)))
    {
        printf("ERROR: La matriz A no es simetrica (A != A^T).\n");
        printf("La factorizacion de Cholesky (A = L * L^T) requiere una matriz simetrica definida positiva.\n");
        return NULL;
    }

    int n = S.A->r;
    double *L = (double *)calloc(n * n, sizeof(double));
    if (!L)
    {
        printf("NO HAY MEMORIA DISPONIBLE\n");
        return NULL;
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            double suma = 0.0;
            for (int k = 0; k < j; ++k)
                suma += L[i * n + k] * L[j * n + k];

            double a_ij = mpq_get_d(S.A->V[I(i, j, *(S.A))]);
            if (i == j)
            {
                double diag = a_ij - suma;
                if (diag <= 1e-14)
                {
                    printf("ERROR: La matriz no es definida positiva (pivote <= 0 en renglon %d).\n", i + 1);
                    free(L);
                    return NULL;
                }
                L[i * n + j] = sqrt(diag);
            }
            else
            {
                L[i * n + j] = (a_ij - suma) / L[j * n + j];
            }
        }
    }

    Matriz *L_mat = crea_MD(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j)
            mpq_set_d(L_mat->V[I(i, j, *L_mat)], L[i * n + j]);
    free(L);

    Matriz *LT_mat = transpuesta(*L_mat);

    printf("Matriz Triangular Inferior L:\n");
    imprime_matriz_dec(L_mat, dec);
    printf("\nMatriz Transpuesta L^T:\n");
    imprime_matriz_dec(LT_mat, dec);

    Matriz *y = sustitucion_adelante(*L_mat, *(S.V), 0);
    Matriz *sol = (y != NULL) ? sustitucion_atras(*LT_mat, *y, 0) : NULL;

    if (sol != NULL)
    {
        printf("\n--- SOLUCION POR CHOLESKY ---\n");
        for (int i = 0; i < n; ++i)
            printf("x_%d = %.*f\n", i + 1, dec, mpq_get_d(sol->V[I(i, 0, *sol)]));
    }

    libera_matriz(y);
    libera_matriz(L_mat);
    libera_matriz(LT_mat);
    return sol;
}

// ============================================================================
// 2. ELIMINACION DE GAUSS Y GAUSS-JORDAN
// ============================================================================

static int *iniciar_orden_variables(int t)
{
    int *x = (int *)malloc(sizeof(int) * t);
    if (x == NULL)
    {
        printf("NO HAY MEMORIA DISPONIBLE\n");
        exit(1);
    }
    for (int i = 0; i < t; ++i)
        x[i] = i + 1;
    return x;
}

static void imprimir_titulo_vista(const char *titulo_frac, const char *titulo_dec, char modo_vista)
{
    if (modo_vista == 'd')
        printf("\n **** %s ****\n", titulo_dec);
    else
        printf("\n **** %s ****\n", titulo_frac);
}

static void imprimir_raices_sistema(Matriz *sol_ordenada, const int *x, int dec, char modo_vista)
{
    printf("\n--- RAICES DEL SISTEMA ---\n");
    for (int i = 0; i < sol_ordenada->r; ++i)
    {
        if (modo_vista == 'd')
            printf("x_%d = %.*f\n", x[i], dec, mpq_get_d(sol_ordenada->V[I(i, 0, *sol_ordenada)]));
        else
            gmp_printf("x_%d = %Qd\n", x[i], sol_ordenada->V[I(i, 0, *sol_ordenada)]);
    }
}

static Matriz *reordenar_solucion(Matriz *sol_perm, const int *x)
{
    int t = sol_perm->r;
    Matriz *sol = crea_MD(t, 1);
    for (int i = 0; i < t; ++i)
        mpq_set(sol->V[I(x[i] - 1, 0, *sol)], sol_perm->V[I(i, 0, *sol_perm)]);
    return sol;
}

Matriz *ejecutar_gauss_res(Sistema S, int tipo_pivoteo, int dec, char pasos, char modo_vista)
{
    if (S.A == NULL || S.V == NULL || !cuadrada(*(S.A)) || S.V->r != S.A->r || S.V->c != 1)
    {
        printf("ERROR: El sistema debe ser cuadrado (A de n x n y b de n x 1).\n");
        return NULL;
    }

    int t = S.A->r;
    int *x = iniciar_orden_variables(t);
    Matriz *M = sis_matriz(*(S.A), *(S.V));
    Matriz *sol = NULL;

    imprimir_titulo_vista("MATRIZ DE FRACCIONES", "MATRIZ EN PUNTO FLOTANTE", modo_vista);
    imprime_aumentada(M, t, dec, modo_vista);
    putchar('\n');

    mpq_t pivote, factor;
    mpq_inits(pivote, factor, NULL);
    int exito = 1;

    for (int r = 0; r < t; ++r)
    {
        if (tipo_pivoteo == 1)
        {
            pivoteo_parcial(M, r, NULL);
        }
        else if (tipo_pivoteo == 2)
        {
            int fr = r, fc = r;
            pivoteo_total(M, r, t, &fr, &fc);
            if (fc != r)
            {
                int tmp = x[r];
                x[r] = x[fc];
                x[fc] = tmp;
            }
        }

        if (mpq_sgn(M->V[I(r, r, *M)]) == 0)
        {
            printf("\n!ERROR: EL SISTEMA NO TIENE SOLUCION UNICA!\n");
            printf("PIVOTE NULO EN LA FILA: %d\n", r + 1);
            exito = 0;
            break;
        }

        mpq_set(pivote, M->V[I(r, r, *M)]);
        div_ren(r, r, pivote, M);

        for (int r_2 = r + 1; r_2 < t; ++r_2)
        {
            mpq_set(factor, M->V[I(r_2, r, *M)]);
            comb_ren(r_2, r, r, factor, M);
        }

        if (pasos == 's')
        {
            printf("---> PASO %d <---\n", r + 1);
            imprime_aumentada(M, t, dec, modo_vista);
            putchar('\n');
        }
    }

    if (exito)
    {
        imprimir_titulo_vista("SISTEMA TRIANGULADO (FRACCIONES)", "SISTEMA TRIANGULADO (PUNTO FLOTANTE)", modo_vista);
        imprime_aumentada(M, t, dec, modo_vista);

        Matriz *U = extraer_submatriz(*M, 0, 0, t, t);
        Matriz *y = extraer_submatriz(*M, 0, t, t, 1);
        Matriz *raices = sustitucion_atras(*U, *y, 1);

        if (raices != NULL)
        {
            imprimir_raices_sistema(raices, x, dec, modo_vista);
            sol = reordenar_solucion(raices, x);
            libera_matriz(raices);
        }
        libera_matriz(U);
        libera_matriz(y);
    }

    mpq_clears(pivote, factor, NULL);
    libera_matriz(M);
    free(x);
    return sol;
}

Matriz *ejecutar_gauss_jordan_res(Sistema S, int tipo_pivoteo, int dec, char pasos, char modo_vista)
{
    if (S.A == NULL || S.V == NULL || !cuadrada(*(S.A)) || S.V->r != S.A->r || S.V->c != 1)
    {
        printf("ERROR: El sistema debe ser cuadrado (A de n x n y b de n x 1).\n");
        return NULL;
    }

    int t = S.A->r;
    int *x = iniciar_orden_variables(t);
    Matriz *M = sis_matriz(*(S.A), *(S.V));
    Matriz *sol = NULL;

    imprimir_titulo_vista("MATRIZ DE FRACCIONES", "MATRIZ EN PUNTO FLOTANTE", modo_vista);
    imprime_aumentada(M, t, dec, modo_vista);
    putchar('\n');

    mpq_t pivote, factor;
    mpq_inits(pivote, factor, NULL);
    int exito = 1;

    for (int k = 0; k < t; ++k)
    {
        if (tipo_pivoteo == 1)
        {
            pivoteo_parcial(M, k, NULL);
        }
        else if (tipo_pivoteo == 2)
        {
            int fr = k, fc = k;
            pivoteo_total(M, k, t, &fr, &fc);
            if (fc != k)
            {
                int tmp = x[k];
                x[k] = x[fc];
                x[fc] = tmp;
            }
        }

        if (mpq_sgn(M->V[I(k, k, *M)]) == 0)
        {
            printf("\n!ERROR: EL SISTEMA NO TIENE SOLUCION UNICA!\n");
            printf("PIVOTE NULO EN LA FILA/COLUMNA: %d\n", k + 1);
            exito = 0;
            break;
        }

        mpq_set(pivote, M->V[I(k, k, *M)]);
        div_ren(k, k, pivote, M);

        for (int r = 0; r < t; ++r)
        {
            if (r != k)
            {
                mpq_set(factor, M->V[I(r, k, *M)]);
                comb_ren(r, k, k, factor, M);
            }
        }

        if (pasos == 's')
        {
            printf("---> PASO %d <---\n", k + 1);
            imprime_aumentada(M, t, dec, modo_vista);
            putchar('\n');
        }
    }

    if (exito)
    {
        imprimir_titulo_vista("SISTEMA RESUELTO (FRACCIONES)", "SISTEMA RESUELTO (PUNTO FLOTANTE)", modo_vista);
        imprime_aumentada(M, t, dec, modo_vista);

        Matriz *raices = extraer_submatriz(*M, 0, t, t, 1);
        if (raices != NULL)
        {
            imprimir_raices_sistema(raices, x, dec, modo_vista);
            sol = reordenar_solucion(raices, x);
            libera_matriz(raices);
        }
    }

    mpq_clears(pivote, factor, NULL);
    libera_matriz(M);
    free(x);
    return sol;
}

// ============================================================================
// 3. GAUSS-JORDAN PARTICIONADO (POR BLOQUES)
// ============================================================================

typedef struct
{
    int tam_bloque;
    int num_bloques;
    Matriz ***A;
    Matriz **B;
} SistemaBloques;

static void imprimir_sistema_bloques(const SistemaBloques *sb)
{
    int N = sb->num_bloques;
    int b = sb->tam_bloque;
    char buffer[256];

    printf("\n=======================================================================\n");
    for (int i = 0; i < N; ++i)
    {
        for (int r = 0; r < b; ++r)
        {
            for (int j = 0; j < N; ++j)
            {
                printf("[ ");
                for (int c = 0; c < b; ++c)
                {
                    Matriz *blk = sb->A[i][j];
                    gmp_snprintf(buffer, sizeof(buffer), "%Qd", blk->V[I(r, c, *blk)]);
                    printf("%-8s ", buffer);
                }
                printf("]  ");
            }

            Matriz *vec = sb->B[i];
            printf("|  [ ");
            gmp_snprintf(buffer, sizeof(buffer), "%Qd", vec->V[I(r, 0, *vec)]);
            printf("%-8s ]\n", buffer);
        }
        putchar('\n');
    }
}

static SistemaBloques particionar_sistema_general(Matriz *M, int tam_bloque)
{
    SistemaBloques sb;
    sb.tam_bloque = tam_bloque;
    sb.num_bloques = M->r / tam_bloque;

    int N = sb.num_bloques;
    int b = sb.tam_bloque;

    sb.A = (Matriz ***)malloc(sizeof(Matriz **) * N);
    sb.B = (Matriz **)malloc(sizeof(Matriz *) * N);

    for (int i = 0; i < N; ++i)
    {
        sb.A[i] = (Matriz **)malloc(sizeof(Matriz *) * N);
        for (int j = 0; j < N; ++j)
            sb.A[i][j] = extraer_submatriz(*M, i * b, j * b, b, b);

        sb.B[i] = extraer_submatriz(*M, i * b, M->r, b, 1);
    }

    return sb;
}

static int resolver_gauss_jordan_bloques(SistemaBloques *sb, char pasos)
{
    int N = sb->num_bloques;

    if (pasos == 's')
    {
        printf("\n---> ESTADO INICIAL PARTICIONADO <---");
        imprimir_sistema_bloques(sb);
    }

    for (int k = 0; k < N; ++k)
    {
        Matriz *A_kk_inv = inversa(*(sb->A[k][k]));
        if (A_kk_inv == NULL)
        {
            printf("\n!ERROR: El bloque diagonal A_%d,%d no es invertible!\n", k + 1, k + 1);
            return 0;
        }

        for (int j = 0; j < N; ++j)
        {
            Matriz *temp_A = producto_matriz(*A_kk_inv, *(sb->A[k][j]));
            libera_matriz(sb->A[k][j]);
            sb->A[k][j] = temp_A;
        }

        Matriz *temp_B = producto_matriz(*A_kk_inv, *(sb->B[k]));
        libera_matriz(sb->B[k]);
        sb->B[k] = temp_B;
        libera_matriz(A_kk_inv);

        for (int i = 0; i < N; ++i)
        {
            if (i != k)
            {
                Matriz *Factor = copia_matriz(*(sb->A[i][k]));

                for (int j = 0; j < N; ++j)
                {
                    Matriz *prod = producto_matriz(*Factor, *(sb->A[k][j]));
                    Matriz *nueva_A = resta_matriz(*(sb->A[i][j]), *prod);
                    libera_matriz(sb->A[i][j]);
                    libera_matriz(prod);
                    sb->A[i][j] = nueva_A;
                }

                Matriz *prod_B = producto_matriz(*Factor, *(sb->B[k]));
                Matriz *nueva_B = resta_matriz(*(sb->B[i]), *prod_B);
                libera_matriz(sb->B[i]);
                libera_matriz(prod_B);
                sb->B[i] = nueva_B;

                libera_matriz(Factor);
            }
        }

        if (pasos == 's')
        {
            printf("\n---> PASO %d <---", k + 1);
            imprimir_sistema_bloques(sb);
        }
    }
    return 1;
}

static void reconstruir_matriz_desde_bloques(Matriz *M, const SistemaBloques *sb)
{
    int N = sb->num_bloques;
    int b = sb->tam_bloque;
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
            insertar_submatriz(M, *(sb->A[i][j]), i * b, j * b);
        insertar_submatriz(M, *(sb->B[i]), i * b, M->r);
    }
}

static void liberar_sistema_bloques(SistemaBloques *sb)
{
    int N = sb->num_bloques;
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
            libera_matriz(sb->A[i][j]);
        free(sb->A[i]);
        libera_matriz(sb->B[i]);
    }
    free(sb->A);
    free(sb->B);
}

Matriz *ejecutar_gauss_particionado_res(Sistema S, int tam_bloque, int dec, char pasos)
{
    if (S.A == NULL || S.V == NULL || !cuadrada(*(S.A)) || S.V->r != S.A->r || S.V->c != 1)
    {
        printf("ERROR: El sistema debe ser cuadrado (A de n x n y b de n x 1).\n");
        return NULL;
    }

    int t = S.A->r;
    if (tam_bloque < 1 || tam_bloque >= t || t % tam_bloque != 0)
    {
        printf("ERROR: El tamano de bloque (%d) debe ser un divisor exacto de %d (menor que %d).\n", tam_bloque, t, t);
        return NULL;
    }

    int *x = iniciar_orden_variables(t);
    Matriz *M = sis_matriz(*(S.A), *(S.V));
    Matriz *sol = NULL;

    printf("\n **** MATRIZ DECIMAL vs MATRIZ DE FRACCIONES ****\n");
    imprime_aumentada(M, t, dec, 'a');

    SistemaBloques sb = particionar_sistema_general(M, tam_bloque);
    if (resolver_gauss_jordan_bloques(&sb, pasos))
    {
        reconstruir_matriz_desde_bloques(M, &sb);
        printf("\n **** SISTEMA RESUELTO ****\n");
        imprime_aumentada(M, t, dec, 'a');

        Matriz *raices = extraer_submatriz(*M, 0, t, t, 1);
        if (raices != NULL)
        {
            printf("\n--- RAICES DEL SISTEMA ---\n");
            for (int i = 0; i < t; ++i)
            {
                printf("x_%d = %-.*f = ", x[i], dec, mpq_get_d(raices->V[I(i, 0, *raices)]));
                gmp_printf("%Qd\n", raices->V[I(i, 0, *raices)]);
            }
            sol = reordenar_solucion(raices, x);
            libera_matriz(raices);
        }
    }

    liberar_sistema_bloques(&sb);
    libera_matriz(M);
    free(x);
    return sol;
}

// ============================================================================
// 4. INVERSA POR GAUSS-JORDAN CON PASOS
// ============================================================================

Matriz *ejecutar_inversa_pasos_res(Matriz A, int dec, char pasos)
{
    if (!cuadrada(A))
    {
        printf("ERROR: La matriz debe ser cuadrada.\n");
        return NULL;
    }

    int t = A.r;
    Matriz *M_inv = crea_MD(t, 2 * t);
    insertar_submatriz(M_inv, A, 0, 0);
    for (int r = 0; r < t; ++r)
        mpq_set_si(M_inv->V[I(r, t + r, *M_inv)], 1, 1);

    mpq_t pivote, factor;
    mpq_inits(pivote, factor, NULL);
    int exito = 1;

    for (int k = 0; k < t; ++k)
    {
        if (!pivoteo_parcial(M_inv, k, NULL))
        {
            printf("\n!LA MATRIZ NO ES INVERTIBLE (Determinante = 0)!\n");
            exito = 0;
            break;
        }

        mpq_set(pivote, M_inv->V[I(k, k, *M_inv)]);
        div_ren(k, k, pivote, M_inv);

        for (int r = 0; r < t; ++r)
        {
            if (r != k)
            {
                mpq_set(factor, M_inv->V[I(r, k, *M_inv)]);
                comb_ren(r, k, k, factor, M_inv);
            }
        }

        if (pasos == 's')
        {
            printf("\n---> PASO %d <---\n", k + 1);
            imprime_aumentada(M_inv, t, dec, 'a');
        }
    }

    Matriz *inv_res = NULL;
    if (exito)
    {
        inv_res = extraer_submatriz(*M_inv, 0, t, t, t);
        printf("\n **** MATRIZ INVERSA (FRACCIONES) ****\n");
        imprime_matriz(inv_res);
        printf("\n **** MATRIZ INVERSA (DECIMAL) ****\n");
        imprime_matriz_dec(inv_res, dec);
    }

    mpq_clears(pivote, factor, NULL);
    libera_matriz(M_inv);
    return inv_res;
}

// ============================================================================
// 5. METODOS ITERATIVOS (JACOBI, GAUSS-SEIDEL, SOR)
// ============================================================================

static void preparar_sistema_iterativo(Matriz *M, int dec, const char *nombre_metodo)
{
    int t = M->r;
    printf("\n **** MATRIZ DECIMAL vs MATRIZ DE FRACCIONES ****\n");
    imprime_aumentada(M, t, dec, 'a');

    if (!es_diagonalmente_dominante(*M))
    {
        printf("\nADVERTENCIA: La matriz no es estrictamente diagonal dominante.\n");
        printf("El metodo de %s podria divergir o requerir mas iteraciones.\n", nombre_metodo);

        char opcion;
        do
        {
            printf("\nDeseas intentar reorganizar las filas para hacerlo diagonal dominante? (s/n): ");
            opcion = vchar();
        } while (opcion != 's' && opcion != 'S' && opcion != 'n' && opcion != 'N');

        if (opcion == 's' || opcion == 'S')
        {
            if (hacer_diagonalmente_dominante(M))
            {
                printf("\nExito! Se reorganizaron las filas para lograr dominancia diagonal.\n");
                printf("\n **** NUEVA MATRIZ REORGANIZADA ****\n");
                imprime_aumentada(M, t, dec, 'a');
            }
            else
            {
                printf("\nNo fue posible lograr dominancia diagonal intercambiando filas.\n");
            }
        }
    }
    else
    {
        printf("\nLa matriz es estrictamente diagonal dominante. La convergencia esta garantizada.\n");
    }
    pausa();
}

static ResultadoIterativo jacobi_iterar(Matriz *M, int dec, double tol, int max_iter)
{
    ResultadoIterativo res = {NULL, 0, 0, 0.0};
    int t = M->r;
    double *x_sol = (double *)calloc(t, sizeof(double));
    double *x_nuevo = (double *)malloc(sizeof(double) * t);
    double *rel_err = (double *)malloc(sizeof(double) * t);

    if (!x_sol || !x_nuevo || !rel_err)
    {
        printf("Error: Memoria insuficiente.\n");
        free(x_sol); free(x_nuevo); free(rel_err);
        return res;
    }

    int iter = 0;
    int convergio = 0;
    double error_max = 0.0;
    int col_w = (dec + 14 < 16) ? 16 : dec + 14;

    printf("\n=======================================================\n");
    printf("              ITERACIONES DEL METODO DE JACOBI          \n");
    printf("=======================================================\n\n");

    printf("%-12s", "Iteracion");
    for (int i = 0; i < t; i++)
    {
        char var_name[32];
        snprintf(var_name, sizeof(var_name), "X_%d", i + 1);
        printf("%-*s", col_w, var_name);
    }
    printf("\n");

    printf("%-12s", "");
    for (int i = 0; i < t; i++)
        printf("%-*.*f", col_w, dec, x_sol[i]);
    printf("\n");

    printf("%-12d", 0);
    for (int i = 0; i < t; i++)
    {
        char err_str[64];
        snprintf(err_str, sizeof(err_str), "e(x_%d) = N/A", i + 1);
        printf("%-*s", col_w, err_str);
    }
    printf("\n");

    for (iter = 1; iter <= max_iter; iter++)
    {
        error_max = 0.0;

        for (int i = 0; i < t; i++)
        {
            double suma = 0.0;
            for (int j = 0; j < t; j++)
            {
                if (i != j)
                    suma += mpq_get_d(M->V[I(i, j, *M)]) * x_sol[j];
            }

            double a_ii = mpq_get_d(M->V[I(i, i, *M)]);
            if (fabs(a_ii) < 1e-12)
            {
                printf("\n!Error: Elemento diagonal A[%d][%d] es cero o muy cercano a cero!\n", i + 1, i + 1);
                free(x_sol); free(x_nuevo); free(rel_err);
                return res;
            }

            double b_i = mpq_get_d(M->V[I(i, t, *M)]);
            x_nuevo[i] = (b_i - suma) / a_ii;

            if (!isfinite(x_nuevo[i]))
            {
                printf("\n!ERROR DE DIVERGENCIA DETECTADO EN LA ITERACION %d!\n", iter);
                printf("El metodo de Jacobi esta divergiendo (Infinito o NaN).\n");
                free(x_sol); free(x_nuevo); free(rel_err);
                return res;
            }

            double diff = fabs(x_nuevo[i] - x_sol[i]);
            rel_err[i] = (fabs(x_nuevo[i]) > 1e-12) ? (diff / fabs(x_nuevo[i])) : diff;

            if (rel_err[i] > error_max)
                error_max = rel_err[i];
        }

        if (!isfinite(error_max))
        {
            printf("\n!ERROR DE DIVERGENCIA DETECTADO EN LA ITERACION %d!\n", iter);
            free(x_sol); free(x_nuevo); free(rel_err);
            return res;
        }

        for (int i = 0; i < t; i++)
            x_sol[i] = x_nuevo[i];

        printf("%-12s", "");
        for (int i = 0; i < t; i++)
            printf("%-*.*f", col_w, dec, x_sol[i]);
        printf("\n");

        printf("%-12d", iter);
        for (int i = 0; i < t; i++)
        {
            char err_str[64];
            snprintf(err_str, sizeof(err_str), "e(x_%d) = %.*f", i + 1, dec, rel_err[i]);
            printf("%-*s", col_w, err_str);
        }
        printf("\n");

        if (error_max <= tol)
        {
            convergio = 1;
            break;
        }
    }

    int iter_realizadas = convergio ? iter : max_iter;
    if (!convergio)
        printf("\nAdvertencia: Alcanzado el limite maximo de iteraciones (%d) sin lograr convergencia.\n", max_iter);
    else
        printf("\nConvergencia lograda exitosamente en %d iteraciones.\n", iter_realizadas);

    printf("\n--- SOLUCION FINAL ---\n");
    res.sol = crea_MD(t, 1);
    res.iteraciones = iter_realizadas;
    res.convergio = convergio;
    res.error_final = error_max;

    for (int i = 0; i < t; i++)
    {
        printf("x_%d = %.*f\n", i + 1, dec, x_sol[i]);
        mpq_set_d(res.sol->V[I(i, 0, *(res.sol))], x_sol[i]);
    }

    free(x_sol);
    free(x_nuevo);
    free(rel_err);
    return res;
}

ResultadoIterativo ejecutar_jacobi_res(Sistema S, int dec, double tol, int max_iter)
{
    ResultadoIterativo res = {NULL, 0, 0, 0.0};
    if (S.A == NULL || S.V == NULL || !cuadrada(*(S.A)) || S.V->r != S.A->r || S.V->c != 1)
    {
        printf("ERROR: El sistema debe ser cuadrado (A de n x n y b de n x 1).\n");
        return res;
    }

    Matriz *M = sis_matriz(*(S.A), *(S.V));
    preparar_sistema_iterativo(M, dec, "Jacobi");
    res = jacobi_iterar(M, dec, tol, max_iter);
    libera_matriz(M);
    return res;
}

static ResultadoIterativo gauss_seidel_sor_iterar(Matriz *M, double omega, int dec, double tol, int max_iter)
{
    ResultadoIterativo res = {NULL, 0, 0, 0.0};
    int t = M->r;
    double *x_sol = (double *)calloc(t, sizeof(double));
    double *rel_err = (double *)malloc(sizeof(double) * t);

    if (!x_sol || !rel_err)
    {
        printf("Error: Memoria insuficiente.\n");
        free(x_sol); free(rel_err);
        return res;
    }

    int iter = 0;
    int convergio = 0;
    double error_max = 0.0;
    int col_w = (dec + 14 < 16) ? 16 : dec + 14;

    printf("\n=======================================================\n");
    if (fabs(omega - 1.0) < 1e-12)
        printf("           ITERACIONES DEL METODO DE GAUSS-SEIDEL       \n");
    else
        printf("      ITERACIONES DEL METODO SOR (omega = %.4f)       \n", omega);
    printf("=======================================================\n\n");

    printf("%-12s", "Iteracion");
    for (int i = 0; i < t; i++)
    {
        char var_name[32];
        snprintf(var_name, sizeof(var_name), "X_%d", i + 1);
        printf("%-*s", col_w, var_name);
    }
    printf("\n");

    printf("%-12s", "");
    for (int i = 0; i < t; i++)
        printf("%-*.*f", col_w, dec, x_sol[i]);
    printf("\n");

    printf("%-12d", 0);
    for (int i = 0; i < t; i++)
    {
        char err_str[64];
        snprintf(err_str, sizeof(err_str), "e(x_%d) = N/A", i + 1);
        printf("%-*s", col_w, err_str);
    }
    printf("\n");

    for (iter = 1; iter <= max_iter; iter++)
    {
        error_max = 0.0;

        for (int i = 0; i < t; i++)
        {
            double suma = 0.0;
            for (int j = 0; j < t; j++)
            {
                if (i != j)
                    suma += mpq_get_d(M->V[I(i, j, *M)]) * x_sol[j];
            }

            double a_ii = mpq_get_d(M->V[I(i, i, *M)]);
            if (fabs(a_ii) < 1e-12)
            {
                printf("\n!Error: Elemento diagonal A[%d][%d] es cero o muy cercano a cero!\n", i + 1, i + 1);
                free(x_sol); free(rel_err);
                return res;
            }

            double b_i = mpq_get_d(M->V[I(i, t, *M)]);
            double x_prev = x_sol[i];
            double x_gs = (b_i - suma) / a_ii;
            x_sol[i] = (1.0 - omega) * x_prev + omega * x_gs;

            if (!isfinite(x_sol[i]))
            {
                printf("\n!ERROR DE DIVERGENCIA DETECTADO EN LA ITERACION %d!\n", iter);
                free(x_sol); free(rel_err);
                return res;
            }

            double diff = fabs(x_sol[i] - x_prev);
            rel_err[i] = (fabs(x_sol[i]) > 1e-12) ? (diff / fabs(x_sol[i])) : diff;

            if (rel_err[i] > error_max)
                error_max = rel_err[i];
        }

        if (!isfinite(error_max))
        {
            printf("\n!ERROR DE DIVERGENCIA DETECTADO EN LA ITERACION %d!\n", iter);
            free(x_sol); free(rel_err);
            return res;
        }

        printf("%-12s", "");
        for (int i = 0; i < t; i++)
            printf("%-*.*f", col_w, dec, x_sol[i]);
        printf("\n");

        printf("%-12d", iter);
        for (int i = 0; i < t; i++)
        {
            char err_str[64];
            snprintf(err_str, sizeof(err_str), "e(x_%d) = %.*f", i + 1, dec, rel_err[i]);
            printf("%-*s", col_w, err_str);
        }
        printf("\n");

        if (error_max <= tol)
        {
            convergio = 1;
            break;
        }
    }

    int iter_realizadas = convergio ? iter : max_iter;
    if (!convergio)
        printf("\nAdvertencia: Alcanzado el limite maximo de iteraciones (%d) sin lograr convergencia.\n", max_iter);
    else
        printf("\nConvergencia lograda exitosamente en %d iteraciones.\n", iter_realizadas);

    printf("\n--- SOLUCION FINAL ---\n");
    res.sol = crea_MD(t, 1);
    res.iteraciones = iter_realizadas;
    res.convergio = convergio;
    res.error_final = error_max;

    for (int i = 0; i < t; i++)
    {
        printf("x_%d = %.*f\n", i + 1, dec, x_sol[i]);
        mpq_set_d(res.sol->V[I(i, 0, *(res.sol))], x_sol[i]);
    }

    free(x_sol);
    free(rel_err);
    return res;
}

ResultadoIterativo ejecutar_gauss_seidel_res(Sistema S, int dec, double tol, int max_iter)
{
    ResultadoIterativo res = {NULL, 0, 0, 0.0};
    if (S.A == NULL || S.V == NULL || !cuadrada(*(S.A)) || S.V->r != S.A->r || S.V->c != 1)
    {
        printf("ERROR: El sistema debe ser cuadrado (A de n x n y b de n x 1).\n");
        return res;
    }

    Matriz *M = sis_matriz(*(S.A), *(S.V));
    preparar_sistema_iterativo(M, dec, "Gauss-Seidel");
    res = gauss_seidel_sor_iterar(M, 1.0, dec, tol, max_iter);
    libera_matriz(M);
    return res;
}

ResultadoIterativo ejecutar_sor_res(Sistema S, double omega, int dec, double tol, int max_iter)
{
    ResultadoIterativo res = {NULL, 0, 0, 0.0};
    if (S.A == NULL || S.V == NULL || !cuadrada(*(S.A)) || S.V->r != S.A->r || S.V->c != 1)
    {
        printf("ERROR: El sistema debe ser cuadrado (A de n x n y b de n x 1).\n");
        return res;
    }

    if (omega <= 0.0 || omega >= 2.0)
    {
        printf("ERROR: El factor de relajacion omega debe estar en (0, 2).\n");
        return res;
    }

    Matriz *M = sis_matriz(*(S.A), *(S.V));
    preparar_sistema_iterativo(M, dec, "SOR");
    res = gauss_seidel_sor_iterar(M, omega, dec, tol, max_iter);
    libera_matriz(M);
    return res;
}

// ============================================================================
// 6. ANALISIS MATRICIAL: NUMERO DE CONDICION Y METODO DE LAS POTENCIAS
// ============================================================================

void ejecutar_numero_condicion(Matriz A, int dec)
{
    if (!cuadrada(A))
    {
        printf("ERROR: El numero de condicion requiere una matriz cuadrada.\n");
        return;
    }

    Matriz *A_inv = inversa(A);
    if (A_inv == NULL)
    {
        printf("Numero de condicion kappa(A) = INFINITO (matriz singular).\n");
        return;
    }

    mpq_t *n1_A = norma_1(A);
    mpq_t *n1_inv = norma_1(*A_inv);
    mpq_t *ninf_A = norma_inf(A);
    mpq_t *ninf_inv = norma_inf(*A_inv);
    double n2_A = norma_espectral(A);
    double n2_inv = norma_espectral(*A_inv);

    mpq_t k1, kinf;
    mpq_inits(k1, kinf, NULL);
    mpq_mul(k1, *n1_A, *n1_inv);
    mpq_mul(kinf, *ninf_A, *ninf_inv);
    double k2 = n2_A * n2_inv;

    printf("NUMERO DE CONDICION DE LA MATRIZ A (%d x %d)\n", A.r, A.c);
    printf("kappa(A) = ||A|| * ||A^-1||\n\n");

    gmp_printf("- Norma L1:        ||A||_1 = %Qd,  ||A^-1||_1 = %Qd\n", *n1_A, *n1_inv);
    gmp_printf("  -> kappa_1(A)   = %Qd  (%.*f)\n\n", k1, dec, mpq_get_d(k1));

    gmp_printf("- Norma L_inf:     ||A||_inf = %Qd,  ||A^-1||_inf = %Qd\n", *ninf_A, *ninf_inv);
    gmp_printf("  -> kappa_inf(A) = %Qd  (%.*f)\n\n", kinf, dec, mpq_get_d(kinf));

    printf("- Norma Espectral: ||A||_2 = %.*f,  ||A^-1||_2 = %.*f\n", dec, n2_A, dec, n2_inv);
    printf("  -> kappa_2(A)   = %.*f\n\n", dec, k2);

    if (k2 < 10.0)
        printf("Diagnostico: Sistema muy bien condicionado (kappa cercano a 1).\n");
    else if (k2 < 1000.0)
        printf("Diagnostico: Condicionamiento moderado.\n");
    else
        printf("Diagnostico: Sistema mal condicionado (alta sensibilidad a perturbaciones).\n");

    mpq_clears(k1, kinf, NULL);
    mpq_clear(*n1_A); free(n1_A);
    mpq_clear(*n1_inv); free(n1_inv);
    mpq_clear(*ninf_A); free(ninf_A);
    mpq_clear(*ninf_inv); free(ninf_inv);
    libera_matriz(A_inv);
}

void ejecutar_metodo_potencias(Matriz A, int inverso, int dec, double tol, int max_iter)
{
    if (!cuadrada(A))
    {
        printf("ERROR: El calculo de eigenvalores requiere una matriz cuadrada.\n");
        return;
    }

    int n = A.r;
    double *A_d = NULL;
    double *L_d = NULL;
    double *U_d = NULL;
    int *perm = NULL;
    double *y_tmp = NULL;

    if (inverso)
    {
        LU *lu = doolittle_pivoteo(A);
        if (lu == NULL)
        {
            printf("ERROR: La matriz es singular; lambda = 0 es eigenvalor y no se puede aplicar potencia inversa.\n");
            return;
        }

        L_d = (double *)malloc(sizeof(double) * n * n);
        U_d = (double *)malloc(sizeof(double) * n * n);
        perm = (int *)malloc(sizeof(int) * n);
        y_tmp = (double *)malloc(sizeof(double) * n);

        if (!L_d || !U_d || !perm || !y_tmp)
        {
            printf("NO HAY MEMORIA DISPONIBLE\n");
            free(L_d); free(U_d); free(perm); free(y_tmp);
            libera_lu(lu);
            return;
        }

        for (int i = 0; i < n; ++i)
        {
            perm[i] = lu->p[i];
            for (int j = 0; j < n; ++j)
            {
                L_d[i * n + j] = mpq_get_d(lu->L->V[I(i, j, *(lu->L))]);
                U_d[i * n + j] = mpq_get_d(lu->U->V[I(i, j, *(lu->U))]);
            }
        }
        libera_lu(lu);
    }
    else
    {
        A_d = (double *)malloc(sizeof(double) * n * n);
        if (!A_d)
        {
            printf("NO HAY MEMORIA DISPONIBLE\n");
            return;
        }

        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                A_d[i * n + j] = mpq_get_d(A.V[I(i, j, A)]);
    }

    double *v = (double *)malloc(sizeof(double) * n);
    double *Mv = (double *)malloc(sizeof(double) * n);
    if (!v || !Mv)
    {
        printf("NO HAY MEMORIA DISPONIBLE\n");
        free(A_d); free(L_d); free(U_d); free(perm); free(y_tmp);
        free(v); free(Mv);
        return;
    }

    double norm0 = 0.0;
    for (int i = 0; i < n; ++i)
    {
        v[i] = 1.0 + 0.1 * (double)i;
        norm0 += v[i] * v[i];
    }
    norm0 = sqrt(norm0);
    for (int i = 0; i < n; ++i)
        v[i] /= norm0;

    double mu = 0.0, mu_prev = 0.0;
    int iter = 0;
    int iter_realizadas = 0;
    int convergio = 0;

    for (iter = 1; iter <= max_iter; ++iter)
    {
        iter_realizadas = iter;

        if (inverso)
        {
            for (int i = 0; i < n; ++i)
            {
                double s = 0.0;
                for (int j = 0; j < i; ++j)
                    s += L_d[i * n + j] * y_tmp[j];
                y_tmp[i] = v[perm[i]] - s;
            }

            for (int i = n - 1; i >= 0; --i)
            {
                double s = 0.0;
                for (int j = i + 1; j < n; ++j)
                    s += U_d[i * n + j] * Mv[j];
                Mv[i] = (y_tmp[i] - s) / U_d[i * n + i];
            }
        }
        else
        {
            for (int i = 0; i < n; ++i)
            {
                double s = 0.0;
                for (int j = 0; j < n; ++j)
                    s += A_d[i * n + j] * v[j];
                Mv[i] = s;
            }
        }

        mu = 0.0;
        double norm_Mv = 0.0;
        for (int i = 0; i < n; ++i)
        {
            mu += v[i] * Mv[i];
            norm_Mv += Mv[i] * Mv[i];
        }
        norm_Mv = sqrt(norm_Mv);

        if (norm_Mv < 1e-15)
        {
            mu = 0.0;
            convergio = 1;
            break;
        }

        int max_idx = 0;
        for (int i = 1; i < n; ++i)
            if (fabs(Mv[i]) > fabs(Mv[max_idx]))
                max_idx = i;
        double signo = (Mv[max_idx] < 0.0) ? -1.0 : 1.0;

        for (int i = 0; i < n; ++i)
            v[i] = signo * Mv[i] / norm_Mv;

        if (iter > 1 && fabs(mu - mu_prev) <= tol * (1.0 + fabs(mu)))
        {
            convergio = 1;
            break;
        }

        mu_prev = mu;
    }

    double lambda = inverso ? ((fabs(mu) > 1e-15) ? (1.0 / mu) : 0.0) : mu;

    printf("RESULTADO DEL METODO DE LAS POTENCIAS %s\n\n",
           inverso ? "(INVERSO - EIGENVALOR MENOR)" : "(DIRECTO - EIGENVALOR DOMINANTE)");

    if (!convergio)
    {
        printf("ADVERTENCIA: Se alcanzo el limite maximo de iteraciones (%d) sin lograr convergencia dentro de la tolerancia.\n\n",
               max_iter);
    }

    printf("- Iteraciones realizadas: %d\n", iter_realizadas);
    printf("- Eigenvalor aproximado:  lambda = %.*f\n\n", dec, lambda);
    printf("- Eigenvector unitario asociado v:\n");
    for (int i = 0; i < n; ++i)
        printf("  v_%d = %.*f\n", i + 1, dec, v[i]);

    // Verificar residuo ||A*v - lambda*v||_2 usando op_matriz
    Matriz *v_mat = crea_MD(n, 1);
    for (int i = 0; i < n; ++i)
        mpq_set_d(v_mat->V[I(i, 0, *v_mat)], v[i]);

    mpq_t lam_q;
    mpq_init(lam_q);
    mpq_set_d(lam_q, lambda);

    Matriz *Av = producto_matriz(A, *v_mat);
    Matriz *lv = producto_e_matriz(lam_q, *v_mat);
    Matriz *res_v = resta_matriz(*Av, *lv);
    printf("\n- Error de eigenpar ||A*v - lambda*v||_2 = %.*e\n", dec, norma_euclidiana(*res_v));

    mpq_clear(lam_q);
    libera_matriz(v_mat);
    libera_matriz(Av);
    libera_matriz(lv);
    libera_matriz(res_v);

    free(A_d);
    free(L_d);
    free(U_d);
    free(perm);
    free(y_tmp);
    free(v);
    free(Mv);
}
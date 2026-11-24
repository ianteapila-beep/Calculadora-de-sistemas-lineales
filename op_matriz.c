#include "op_matriz.h"

// ============================================================================
// ENTRADA SEGURA Y UTILIDADES
// ============================================================================

void limpiabuffer(void)
{
    int a;
    while ((a = getchar()) != '\n' && a != EOF);
}

void pausa(void)
{
    printf("\nPresiona ENTER para continuar...");
    limpiabuffer();
}

int vint(void)
{
    int num = 0, res;
    while ((res = scanf("%d", &num)) != 1)
    {
        if (res == EOF) exit(0);
        printf("Dato invalido, intenta de nuevo: ");
        limpiabuffer();
    }
    limpiabuffer();
    return num;
}

double vdouble(void)
{
    double num = 0.0;
    int res;
    while ((res = scanf("%lf", &num)) != 1)
    {
        if (res == EOF) exit(0);
        printf("Dato invalido, intenta de nuevo: ");
        limpiabuffer();
    }
    limpiabuffer();
    return num;
}

char vchar(void)
{
    char r;
    int res;
    while ((res = scanf(" %c", &r)) != 1)
    {
        if (res == EOF) exit(0);
        printf("Dato invalido, intenta de nuevo: ");
        limpiabuffer();
    }
    limpiabuffer();
    return r;
}

// ============================================================================
// GESTION DE MEMORIA, INDEXACION Y SUBMATRICES
// ============================================================================

int I(int r, int c, Matriz X)
{
    if (r < 0 || c < 0 || r >= X.r || c >= X.c)
    {
        printf("INDICE FUERA DEL RANGO DE LA MATRIZ\n");
        return -1;
    }
    return r * X.c + c;
}

Matriz *crea_MD(int r, int c)
{
    Matriz *V = (Matriz *)malloc(sizeof(Matriz));
    if (V == NULL)
    {
        printf("NO HAY MEMORIA DISPONIBLE\n");
        exit(1);
    }
    V->r = r;
    V->c = c;
    V->V = (mpq_t *)malloc(sizeof(mpq_t) * c * r);
    if (V->V == NULL)
    {
        printf("NO HAY MEMORIA DISPONIBLE\n");
        exit(1);
    }
    for (int i = 0; i < r * c; ++i)
        mpq_init(V->V[i]);
    return V;
}

void libera_matriz(Matriz *V)
{
    if (V == NULL)
        return;

    if (V->V != NULL)
    {
        for (int i = 0; i < V->r * V->c; ++i)
            mpq_clear(V->V[i]);
        free(V->V);
        V->V = NULL;
    }
    free(V);
}

Matriz *copia_matriz(Matriz B)
{
    Matriz *A = crea_MD(B.r, B.c);
    if (A == NULL) return NULL;

    for (int i = 0; i < A->r * A->c; ++i)
        mpq_set(A->V[i], B.V[i]);

    return A;
}

Matriz *sis_matriz(Matriz A, Matriz V)
{
    if (A.r != V.r || V.c != 1)
    {
        printf("Matriz y vector no compatibles\n");
        return NULL;
    }
    Matriz *S = crea_MD(A.r, A.c + V.c);
    for (int i = 0; i < S->r; ++i)
    {
        for (int j = 0; j < S->c; ++j)
        {
            if (j < A.c)
                mpq_set(S->V[I(i, j, *S)], A.V[I(i, j, A)]);
            else
                mpq_set(S->V[I(i, j, *S)], V.V[i]);
        }
    }
    return S;
}

Matriz *extraer_submatriz(Matriz A, int r0, int c0, int filas, int cols)
{
    if (r0 < 0 || c0 < 0 || filas <= 0 || cols <= 0 ||
        r0 + filas > A.r || c0 + cols > A.c)
    {
        printf("Submatriz fuera de rango.\n");
        return NULL;
    }

    Matriz *sub = crea_MD(filas, cols);
    for (int i = 0; i < filas; ++i)
        for (int j = 0; j < cols; ++j)
            mpq_set(sub->V[I(i, j, *sub)], A.V[I(r0 + i, c0 + j, A)]);

    return sub;
}

void insertar_submatriz(Matriz *dest, Matriz src, int r0, int c0)
{
    if (dest == NULL || r0 < 0 || c0 < 0 ||
        r0 + src.r > dest->r || c0 + src.c > dest->c)
        return;

    for (int i = 0; i < src.r; ++i)
        for (int j = 0; j < src.c; ++j)
            mpq_set(dest->V[I(r0 + i, c0 + j, *dest)], src.V[I(i, j, src)]);
}

// ============================================================================
// OPERACIONES ELEMENTALES Y PIVOTEO
// ============================================================================

void inter_ren(int des, int ori, Matriz *A)
{
    if (des == ori) return;
    for (int i = 0; i < A->c; ++i)
        mpq_swap(A->V[I(des, i, *A)], A->V[I(ori, i, *A)]);
}

void inter_col(int des, int ori, Matriz *A)
{
    if (des == ori) return;
    for (int i = 0; i < A->r; ++i)
        mpq_swap(A->V[I(i, des, *A)], A->V[I(i, ori, *A)]);
}

void div_ren(int r, int col_inicio, mpq_t divisor, Matriz *A)
{
    if (A == NULL || r < 0 || r >= A->r || mpq_sgn(divisor) == 0)
        return;

    for (int c = col_inicio; c < A->c; ++c)
        mpq_div(A->V[I(r, c, *A)], A->V[I(r, c, *A)], divisor);
}

void comb_ren(int r_dest, int r_piv, int col_inicio, mpq_t factor, Matriz *A)
{
    if (A == NULL || r_dest < 0 || r_dest >= A->r || r_piv < 0 || r_piv >= A->r)
        return;
    if (mpq_sgn(factor) == 0)
        return;

    mpq_t aux;
    mpq_init(aux);
    for (int c = col_inicio; c < A->c; ++c)
    {
        mpq_mul(aux, factor, A->V[I(r_piv, c, *A)]);
        mpq_sub(A->V[I(r_dest, c, *A)], A->V[I(r_dest, c, *A)], aux);
    }
    mpq_clear(aux);
}

int pivoteo_parcial_idx(Matriz *A, int k, int *cam, int *fila_piv)
{
    if (A == NULL || k < 0 || k >= A->r || k >= A->c)
        return 0;

    int fila_max = k;
    mpq_t max_val, abs_val;
    mpq_init(max_val);
    mpq_init(abs_val);

    mpq_abs(max_val, A->V[I(k, k, *A)]);

    for (int r = k + 1; r < A->r; ++r)
    {
        mpq_abs(abs_val, A->V[I(r, k, *A)]);
        if (mpq_cmp(abs_val, max_val) > 0)
        {
            mpq_set(max_val, abs_val);
            fila_max = r;
        }
    }

    if (mpq_sgn(max_val) == 0)
    {
        mpq_clear(max_val);
        mpq_clear(abs_val);
        return 0;
    }

    if (fila_max != k)
    {
        inter_ren(k, fila_max, A);
        if (cam != NULL)
            ++(*cam);
    }

    if (fila_piv != NULL)
        *fila_piv = fila_max;

    mpq_clear(max_val);
    mpq_clear(abs_val);
    return 1;
}

int pivoteo_parcial(Matriz *A, int k, int *cam)
{
    return pivoteo_parcial_idx(A, k, cam, NULL);
}

int pivoteo_total(Matriz *A, int k, int lim_col, int *fila_piv, int *col_piv)
{
    if (A == NULL || k < 0 || k >= A->r || k >= lim_col || lim_col > A->c)
        return 0;

    int max_r = k, max_c = k;
    mpq_t max_val, abs_val;
    mpq_init(max_val);
    mpq_init(abs_val);

    mpq_abs(max_val, A->V[I(k, k, *A)]);

    for (int r = k; r < A->r; ++r)
    {
        for (int c = k; c < lim_col; ++c)
        {
            mpq_abs(abs_val, A->V[I(r, c, *A)]);
            if (mpq_cmp(abs_val, max_val) > 0)
            {
                mpq_set(max_val, abs_val);
                max_r = r;
                max_c = c;
            }
        }
    }

    if (mpq_sgn(max_val) == 0)
    {
        mpq_clear(max_val);
        mpq_clear(abs_val);
        return 0;
    }

    if (max_r != k)
        inter_ren(k, max_r, A);
    if (max_c != k)
        inter_col(k, max_c, A);

    if (fila_piv != NULL) *fila_piv = max_r;
    if (col_piv != NULL)  *col_piv = max_c;

    mpq_clear(max_val);
    mpq_clear(abs_val);
    return 1;
}

// ============================================================================
// PROPIEDADES DE LAS MATRICES Y DOMINANCIA DIAGONAL
// ============================================================================

int cuadrada(Matriz A)
{
    return (A.c == A.r) ? 1 : 0;
}

int simetrica(Matriz A)
{
    if (!cuadrada(A))
        return 0;

    for (int i = 0; i < A.r; ++i)
    {
        for (int j = i + 1; j < A.c; ++j)
        {
            if (mpq_equal(A.V[I(i, j, A)], A.V[I(j, i, A)]) == 0)
                return 0;
        }
    }
    return 1;
}

int vect_null(Matriz A)
{
    for (int i = 0; i < A.r; ++i)
    {
        int nz = 0;
        for (int j = 0; j < A.c; ++j)
        {
            if (mpq_sgn(A.V[I(i, j, A)]) != 0) { nz = 1; break; }
        }
        if (!nz) return 0;
    }

    for (int j = 0; j < A.c; ++j)
    {
        int nz = 0;
        for (int i = 0; i < A.r; ++i)
        {
            if (mpq_sgn(A.V[I(i, j, A)]) != 0) { nz = 1; break; }
        }
        if (!nz) return 0;
    }
    return 1;
}

mpq_t *determinante(Matriz B)
{
    if (!cuadrada(B))
    {
        printf("La matriz no es cuadrada.\n");
        return NULL;
    }
    Matriz *A = copia_matriz(B);
    mpq_t *det = (mpq_t *)malloc(sizeof(mpq_t));
    if (det == NULL)
    {
        printf("No hay memoria disponible\n");
        exit(1);
    }
    int cambios_signo = 0;
    mpq_init(*det);
    *A = diag_sup(*A, &cambios_signo);
    mpq_set(*det, A->V[0]);

    for (int i = 1; i < A->c; ++i)
    {
        if (mpq_sgn(A->V[I(i, i, *A)]) == 0)
        {
            mpq_set_si(*det, 0, 1);
            libera_matriz(A);
            return det;
        }
        mpq_mul(*det, *det, A->V[I(i, i, *A)]);
    }

    if ((cambios_signo % 2) == 1)
        mpq_neg(*det, *det);

    libera_matriz(A);
    return det;
}

mpq_t *traza(Matriz A)
{
    if (!cuadrada(A))
    {
        printf("La traza solo esta definida para matrices cuadradas.\n");
        return NULL;
    }

    mpq_t *tr = (mpq_t *)malloc(sizeof(mpq_t));
    if (tr == NULL)
    {
        printf("No hay memoria disponible\n");
        exit(1);
    }
    mpq_init(*tr);
    mpq_set_si(*tr, 0, 1);

    for (int i = 0; i < A.r; ++i)
        mpq_add(*tr, *tr, A.V[I(i, i, A)]);

    return tr;
}

int rango_matriz(Matriz A)
{
    if (A.r <= 0 || A.c <= 0 || A.V == NULL)
        return 0;

    Matriz *M = copia_matriz(A);
    mpq_t factor;
    mpq_init(factor);

    int fila = 0;
    for (int col = 0; col < M->c && fila < M->r; ++col)
    {
        int piv = -1;
        for (int r = fila; r < M->r; ++r)
        {
            if (mpq_sgn(M->V[I(r, col, *M)]) != 0)
            {
                piv = r;
                break;
            }
        }

        if (piv == -1)
            continue;

        if (piv != fila)
            inter_ren(fila, piv, M);

        for (int r = fila + 1; r < M->r; ++r)
        {
            if (mpq_sgn(M->V[I(r, col, *M)]) != 0)
            {
                mpq_div(factor, M->V[I(r, col, *M)], M->V[I(fila, col, *M)]);
                comb_ren(r, fila, col, factor, M);
            }
        }
        fila++;
    }

    mpq_clear(factor);
    libera_matriz(M);
    return fila;
}

static int es_fila_dominante_mat(Matriz M, int r, int col_diag, int n)
{
    mpq_t suma, abs_val, diag_abs;
    mpq_inits(suma, abs_val, diag_abs, NULL);

    for (int c = 0; c < n; ++c)
    {
        if (c != col_diag)
        {
            mpq_abs(abs_val, M.V[I(r, c, M)]);
            mpq_add(suma, suma, abs_val);
        }
    }
    mpq_abs(diag_abs, M.V[I(r, col_diag, M)]);
    int dom = (mpq_cmp(diag_abs, suma) > 0) ? 1 : 0;

    mpq_clears(suma, abs_val, diag_abs, NULL);
    return dom;
}

int es_diagonalmente_dominante(Matriz A)
{
    int n = A.r;
    if (A.c < n) return 0;

    for (int i = 0; i < n; ++i)
    {
        if (!es_fila_dominante_mat(A, i, i, n))
            return 0;
    }
    return 1;
}

static int buscar_perm_dom_mat(Matriz M, int n, int pos, int *usados, int *orden)
{
    if (pos == n) return 1;
    for (int r = 0; r < n; ++r)
    {
        if (!usados[r] && es_fila_dominante_mat(M, r, pos, n))
        {
            usados[r] = 1;
            orden[pos] = r;
            if (buscar_perm_dom_mat(M, n, pos + 1, usados, orden))
                return 1;
            usados[r] = 0;
        }
    }
    return 0;
}

int hacer_diagonalmente_dominante(Matriz *M)
{
    if (M == NULL || M->c < M->r) return 0;
    int n = M->r;

    int *usados = (int *)calloc(n, sizeof(int));
    int *orden = (int *)malloc(sizeof(int) * n);
    if (!usados || !orden)
    {
        free(usados);
        free(orden);
        return 0;
    }

    int exito = 0;
    if (buscar_perm_dom_mat(*M, n, 0, usados, orden))
    {
        Matriz *copia = copia_matriz(*M);
        for (int r = 0; r < n; ++r)
        {
            int r_orig = orden[r];
            for (int c = 0; c < M->c; ++c)
                mpq_set(M->V[I(r, c, *M)], copia->V[I(r_orig, c, *copia)]);
        }
        libera_matriz(copia);
        exito = 1;
    }

    free(usados);
    free(orden);
    return exito;
}

// ============================================================================
// OPERACIONES CON MATRICES Y SISTEMAS TRIANGULARES
// ============================================================================

Matriz *suma_matriz(Matriz A, Matriz B)
{
    if (A.r != B.r || A.c != B.c)
    {
        printf("NO ES POSIBLE REALIZAR LA SUMA DE LAS MATRICES (%dx%d vs %dx%d)\n", A.r, A.c, B.r, B.c);
        return NULL;
    }
    Matriz *C = crea_MD(A.r, A.c);
    for (int i = 0; i < A.r * A.c; ++i)
        mpq_add(C->V[i], A.V[i], B.V[i]);
    return C;
}

Matriz *resta_matriz(Matriz A, Matriz B)
{
    if (A.r != B.r || A.c != B.c)
    {
        printf("NO ES POSIBLE REALIZAR LA RESTA DE LAS MATRICES (%dx%d vs %dx%d)\n", A.r, A.c, B.r, B.c);
        return NULL;
    }
    Matriz *C = crea_MD(A.r, A.c);
    for (int i = 0; i < A.r; ++i)
        for (int j = 0; j < A.c; ++j)
            mpq_sub(C->V[I(i, j, *C)], A.V[I(i, j, A)], B.V[I(i, j, B)]);
    return C;
}

Matriz *producto_matriz(Matriz A, Matriz B)
{
    if (A.c != B.r)
    {
        printf("NO ES POSIBLE EL PRODUCTO DE UNA MATRIZ DE %d x %d CON UNA DE %d x %d\n", A.r, A.c, B.r, B.c);
        return NULL;
    }
    Matriz *C = crea_MD(A.r, B.c);
    mpq_t aux;
    mpq_init(aux);
    for (int i = 0; i < A.r; ++i)
    {
        for (int j = 0; j < B.c; ++j)
        {
            for (int k = 0; k < A.c; ++k)
            {
                mpq_mul(aux, A.V[I(i, k, A)], B.V[I(k, j, B)]);
                mpq_add(C->V[I(i, j, *C)], C->V[I(i, j, *C)], aux);
            }
        }
    }
    mpq_clear(aux);
    return C;
}

Matriz *producto_e_matriz(mpq_t k, Matriz A)
{
    Matriz *C = crea_MD(A.r, A.c);
    for (int i = 0; i < A.r * A.c; ++i)
        mpq_mul(C->V[i], k, A.V[i]);
    return C;
}

Matriz *transpuesta(Matriz A)
{
    Matriz *T = crea_MD(A.c, A.r);
    for (int i = 0; i < A.r; ++i)
    {
        for (int j = 0; j < A.c; ++j)
            mpq_set(T->V[I(j, i, *T)], A.V[I(i, j, A)]);
    }
    return T;
}

Matriz *potencia_matriz(Matriz A, int exp)
{
    if (!cuadrada(A))
    {
        printf("La matriz debe ser cuadrada para calcular A^k.\n");
        return NULL;
    }

    int n = A.r;
    Matriz *base = NULL;
    int e = exp;

    if (exp < 0)
    {
        base = inversa(A);
        if (base == NULL)
            return NULL;
        e = -exp;
    }
    else
    {
        base = copia_matriz(A);
    }

    Matriz *res = crea_MD(n, n);
    for (int i = 0; i < n; ++i)
        mpq_set_si(res->V[I(i, i, *res)], 1, 1);

    while (e > 0)
    {
        if (e % 2 == 1)
        {
            Matriz *tmp = producto_matriz(*res, *base);
            libera_matriz(res);
            res = tmp;
        }
        e /= 2;
        if (e > 0)
        {
            Matriz *tmp_b = producto_matriz(*base, *base);
            libera_matriz(base);
            base = tmp_b;
        }
    }

    libera_matriz(base);
    return res;
}

Matriz *inversa(Matriz A)
{
    if (!cuadrada(A))
    {
        printf("La matriz debe ser cuadrada para calcular su inversa.\n");
        return NULL;
    }

    int n = A.r;
    Matriz *aug = crea_MD(n, 2 * n);
    insertar_submatriz(aug, A, 0, 0);
    for (int r = 0; r < n; ++r)
        mpq_set_si(aug->V[I(r, n + r, *aug)], 1, 1);

    mpq_t pivote, factor;
    mpq_inits(pivote, factor, NULL);

    for (int k = 0; k < n; ++k)
    {
        if (!pivoteo_parcial(aug, k, NULL))
        {
            printf("La matriz no es invertible (Determinante = 0).\n");
            mpq_clears(pivote, factor, NULL);
            libera_matriz(aug);
            return NULL;
        }

        mpq_set(pivote, aug->V[I(k, k, *aug)]);
        div_ren(k, k, pivote, aug);

        for (int r = 0; r < n; ++r)
        {
            if (r != k)
            {
                mpq_set(factor, aug->V[I(r, k, *aug)]);
                comb_ren(r, k, k, factor, aug);
            }
        }
    }

    Matriz *inv = extraer_submatriz(*aug, 0, n, n, n);
    mpq_clears(pivote, factor, NULL);
    libera_matriz(aug);
    return inv;
}

Matriz diag_sup(Matriz A, int *cam)
{
    mpq_t piv, sig, op;
    mpq_inits(piv, sig, op, NULL);

    for (int i = 0; i < A.r - 1; ++i)
    {
        mpq_set(piv, A.V[I(i, i, A)]);

        if (mpq_sgn(piv) == 0)
        {
            int aux = A.r;
            while ((--aux > i) && mpq_sgn(A.V[I(aux, i, A)]) == 0);

            if (aux == i)
                continue;

            if (cam != NULL) ++(*cam);
            inter_ren(i, aux, &A);
            mpq_set(piv, A.V[I(i, i, A)]);
        }

        for (int r = i + 1; r < A.r; ++r)
        {
            mpq_set(sig, A.V[I(r, i, A)]);
            if (mpq_sgn(sig) == 0) continue;
            mpq_div(op, sig, piv);
            comb_ren(r, i, i, op, &A);
        }
    }

    mpq_clears(piv, sig, op, NULL);
    return A;
}

Matriz diag_inf(Matriz A, int *cam)
{
    mpq_t piv, sig, op;
    mpq_inits(piv, sig, op, NULL);

    for (int i = A.r - 1; i > 0; --i)
    {
        mpq_set(piv, A.V[I(i, i, A)]);

        if (mpq_sgn(piv) == 0)
        {
            int aux = i;
            while ((--aux >= 0) && mpq_sgn(A.V[I(aux, i, A)]) == 0);

            if (aux < 0)
                continue;

            if (cam != NULL) ++(*cam);
            inter_ren(i, aux, &A);
            mpq_set(piv, A.V[I(i, i, A)]);
        }

        for (int r = i - 1; r >= 0; --r)
        {
            mpq_set(sig, A.V[I(r, i, A)]);
            if (mpq_sgn(sig) == 0) continue;
            mpq_div(op, sig, piv);
            comb_ren(r, i, 0, op, &A);
        }
    }

    mpq_clears(piv, sig, op, NULL);
    return A;
}

Matriz *permutar_vector(Matriz b, const int *p)
{
    if (b.c != 1 || p == NULL) return NULL;
    Matriz *pb = crea_MD(b.r, 1);
    for (int i = 0; i < b.r; ++i)
        mpq_set(pb->V[I(i, 0, *pb)], b.V[I(p[i], 0, b)]);
    return pb;
}

Matriz *sustitucion_adelante(Matriz L, Matriz b, int diag_unitaria)
{
    if (!cuadrada(L) || b.r != L.r || b.c != 1)
        return NULL;

    int n = L.r;
    Matriz *y = crea_MD(n, 1);
    mpq_t suma, prod, num;
    mpq_inits(suma, prod, num, NULL);

    for (int i = 0; i < n; ++i)
    {
        mpq_set_si(suma, 0, 1);
        for (int j = 0; j < i; ++j)
        {
            mpq_mul(prod, L.V[I(i, j, L)], y->V[I(j, 0, *y)]);
            mpq_add(suma, suma, prod);
        }
        mpq_sub(num, b.V[I(i, 0, b)], suma);

        if (diag_unitaria)
        {
            mpq_set(y->V[I(i, 0, *y)], num);
        }
        else
        {
            if (mpq_sgn(L.V[I(i, i, L)]) == 0)
            {
                mpq_clears(suma, prod, num, NULL);
                libera_matriz(y);
                return NULL;
            }
            mpq_div(y->V[I(i, 0, *y)], num, L.V[I(i, i, L)]);
        }
    }

    mpq_clears(suma, prod, num, NULL);
    return y;
}

Matriz *sustitucion_atras(Matriz U, Matriz y, int diag_unitaria)
{
    if (!cuadrada(U) || y.r != U.r || y.c != 1)
        return NULL;

    int n = U.r;
    Matriz *x = crea_MD(n, 1);
    mpq_t suma, prod, num;
    mpq_inits(suma, prod, num, NULL);

    for (int i = n - 1; i >= 0; --i)
    {
        mpq_set_si(suma, 0, 1);
        for (int j = i + 1; j < n; ++j)
        {
            mpq_mul(prod, U.V[I(i, j, U)], x->V[I(j, 0, *x)]);
            mpq_add(suma, suma, prod);
        }
        mpq_sub(num, y.V[I(i, 0, y)], suma);

        if (diag_unitaria)
        {
            mpq_set(x->V[I(i, 0, *x)], num);
        }
        else
        {
            if (mpq_sgn(U.V[I(i, i, U)]) == 0)
            {
                mpq_clears(suma, prod, num, NULL);
                libera_matriz(x);
                return NULL;
            }
            mpq_div(x->V[I(i, 0, *x)], num, U.V[I(i, i, U)]);
        }
    }

    mpq_clears(suma, prod, num, NULL);
    return x;
}

// ============================================================================
// LECTURA E IMPRESION ALINEADA
// ============================================================================

static int ancho_mpq(mpq_t val)
{
    return gmp_snprintf(NULL, 0, "% Qd", val);
}

static int ancho_double(double val, int dec)
{
    return snprintf(NULL, 0, "%.*f", dec, val);
}

void imprime_matriz(Matriz *V)
{
    if (V == NULL || V->V == NULL)
    {
        printf("Matriz vacia o no inicializada.\n");
        return;
    }

    int *anchos = (int *)calloc(V->c, sizeof(int));
    if (anchos == NULL)
    {
        printf("NO HAY MEMORIA DISPONIBLE\n");
        exit(1);
    }

    for (int i = 0; i < V->r; ++i)
    {
        for (int j = 0; j < V->c; ++j)
        {
            int len = ancho_mpq(V->V[I(i, j, *V)]);
            if (len > anchos[j])
                anchos[j] = len;
        }
    }

    for (int i = 0; i < V->r; ++i)
    {
        printf("[ ");
        for (int j = 0; j < V->c; ++j)
        {
            int idx = I(i, j, *V);
            int len = ancho_mpq(V->V[idx]);

            for (int s = 0; s < anchos[j] - len; ++s)
                putchar(' ');

            gmp_printf("% Qd", V->V[idx]);

            if (j < V->c - 1)
                printf("  ");
        }
        printf(" ]\n");
    }

    free(anchos);
}

void imprime_matriz_dec(Matriz *V, int dec)
{
    if (V == NULL || V->V == NULL)
    {
        printf("Matriz vacia o no inicializada.\n");
        return;
    }

    int *anchos = (int *)calloc(V->c, sizeof(int));
    if (anchos == NULL)
    {
        printf("NO HAY MEMORIA DISPONIBLE\n");
        exit(1);
    }

    for (int i = 0; i < V->r; ++i)
    {
        for (int j = 0; j < V->c; ++j)
        {
            double val = mpq_get_d(V->V[I(i, j, *V)]);
            int len = ancho_double(val, dec);
            if (len > anchos[j])
                anchos[j] = len;
        }
    }

    for (int i = 0; i < V->r; ++i)
    {
        printf("[ ");
        for (int j = 0; j < V->c; ++j)
        {
            double val = mpq_get_d(V->V[I(i, j, *V)]);
            printf("%*.*f", anchos[j], dec, val);
            if (j < V->c - 1)
                printf("  ");
        }
        printf(" ]\n");
    }

    free(anchos);
}

void imprime_sistema(Sistema S)
{
    if (S.A == NULL || S.V == NULL || S.A->V == NULL || S.V->V == NULL)
    {
        printf("Sistema vacio o no inicializado.\n");
        return;
    }

    if (S.A->r != S.V->r)
    {
        printf("ERROR: La matriz A (%d renglones) y el vector b (%d renglones) no coinciden.\n", S.A->r, S.V->r);
        return;
    }

    Matriz *aug = sis_matriz(*(S.A), *(S.V));
    imprime_aumentada(aug, S.A->c, 5, 'f');
    libera_matriz(aug);
}

void imprime_sistema_dec(Sistema S, int dec)
{
    if (S.A == NULL || S.V == NULL || S.A->V == NULL || S.V->V == NULL)
    {
        printf("Sistema vacio o no inicializado.\n");
        return;
    }

    if (S.A->r != S.V->r)
    {
        printf("ERROR: La matriz A (%d renglones) y el vector b (%d renglones) no coinciden.\n", S.A->r, S.V->r);
        return;
    }

    Matriz *aug = sis_matriz(*(S.A), *(S.V));
    imprime_aumentada(aug, S.A->c, dec, 'd');
    libera_matriz(aug);
}

void imprime_aumentada(Matriz *M, int col_sep, int dec, char modo_vista)
{
    if (M == NULL || M->V == NULL)
    {
        printf("Matriz vacia o no inicializada.\n");
        return;
    }

    int *ancho_dec = (int *)calloc(M->c, sizeof(int));
    int *ancho_frac = (int *)calloc(M->c, sizeof(int));
    if (ancho_dec == NULL || ancho_frac == NULL)
    {
        printf("NO HAY MEMORIA DISPONIBLE\n");
        exit(1);
    }

    for (int r = 0; r < M->r; ++r)
    {
        for (int c = 0; c < M->c; ++c)
        {
            int ld = ancho_double(mpq_get_d(M->V[I(r, c, *M)]), dec);
            if (ld > ancho_dec[c]) ancho_dec[c] = ld;

            int lf = ancho_mpq(M->V[I(r, c, *M)]);
            if (lf > ancho_frac[c]) ancho_frac[c] = lf;
        }
    }

    int sep_entre_vistas = (M->r <= 3) ? 16 : 8;

    for (int r = 0; r < M->r; ++r)
    {
        if (modo_vista == 'd' || modo_vista == 'a')
        {
            printf("[ ");
            for (int c = 0; c < M->c; ++c)
            {
                if (c == col_sep) printf("| ");
                double val = mpq_get_d(M->V[I(r, c, *M)]);
                printf("%*.*f", ancho_dec[c], dec, val);
                if (c < M->c - 1) printf("  ");
            }
            printf(" ]");
        }

        if (modo_vista == 'a')
            printf("%*c", sep_entre_vistas, ' ');

        if (modo_vista == 'f' || modo_vista == 'a')
        {
            printf("[ ");
            for (int c = 0; c < M->c; ++c)
            {
                if (c == col_sep) printf("| ");
                int idx = I(r, c, *M);
                int len = ancho_mpq(M->V[idx]);
                for (int s = 0; s < ancho_frac[c] - len; ++s)
                    putchar(' ');
                gmp_printf("% Qd", M->V[idx]);
                if (c < M->c - 1) printf("  ");
            }
            printf(" ]");
        }
        putchar('\n');
    }

    free(ancho_dec);
    free(ancho_frac);
}

void leer_decimal_o_fraccion(char *tok, mpq_t V)
{
    char *exp_ptr = strpbrk(tok, "eE");
    char *punto = strchr(tok, '.');

    if (punto != NULL || exp_ptr != NULL)
    {
        char mantisa[MAX];
        long exp10 = 0;

        if (exp_ptr != NULL)
        {
            size_t len_m = (size_t)(exp_ptr - tok);
            if (len_m == 0 || len_m >= MAX)
            {
                printf("Error leyendo elemento '%s'. Asignando 0.\n", tok);
                mpq_set_si(V, 0, 1);
                return;
            }
            memcpy(mantisa, tok, len_m);
            mantisa[len_m] = '\0';

            char *endptr = NULL;
            exp10 = strtol(exp_ptr + 1, &endptr, 10);
            if (endptr == exp_ptr + 1 || *endptr != '\0')
            {
                printf("Error leyendo elemento '%s'. Asignando 0.\n", tok);
                mpq_set_si(V, 0, 1);
                return;
            }
        }
        else
        {
            strncpy(mantisa, tok, MAX - 1);
            mantisa[MAX - 1] = '\0';
        }

        int decimales = 0;
        char temp_str[MAX];
        int i = 0, j = 0, puntos = 0;
        while (mantisa[i] != '\0')
        {
            if (mantisa[i] == '.')
            {
                puntos++;
                decimales = (int)strlen(mantisa) - i - 1;
            }
            else
            {
                temp_str[j++] = mantisa[i];
            }
            i++;
        }
        temp_str[j] = '\0';

        if (puntos > 1)
        {
            printf("Error leyendo elemento '%s'. Asignando 0.\n", tok);
            mpq_set_si(V, 0, 1);
            return;
        }

        mpz_t num, den, pow10;
        mpz_init(num);
        mpz_init(den);
        mpz_init(pow10);

        if (mpz_set_str(num, temp_str, 10) == 0)
        {
            long exp_neto = exp10 - (long)decimales;
            if (exp_neto >= 0)
            {
                mpz_ui_pow_ui(pow10, 10, (unsigned long)exp_neto);
                mpz_mul(num, num, pow10);
                mpz_set_ui(den, 1);
            }
            else
            {
                mpz_ui_pow_ui(den, 10, (unsigned long)(-exp_neto));
            }
            mpq_set_num(V, num);
            mpq_set_den(V, den);
        }
        else
        {
            printf("Error leyendo elemento '%s'. Asignando 0.\n", tok);
            mpq_set_si(V, 0, 1);
        }

        mpz_clear(num);
        mpz_clear(den);
        mpz_clear(pow10);
    }
    else
    {
        if (mpq_set_str(V, tok, 10) != 0)
        {
            printf("Error leyendo elemento '%s'. Asignando 0.\n", tok);
            mpq_set_si(V, 0, 1);
        }
    }
    mpq_canonicalize(V);
}

void ingresa_fraccion(char *f, Matriz *V, int r)
{
    char clean_f[MAX];
    int j = 0;

    for (int i = 0; f[i] != '\0'; i++)
        if (!isspace((unsigned char)f[i]))
            clean_f[j++] = f[i];

    clean_f[j] = '\0';

    const char *del = ",";
    char *tok = strtok(clean_f, del);

    int c = 0;
    while (c < V->c && tok != NULL)
    {
        leer_decimal_o_fraccion(tok, V->V[I(r, c, *V)]);
        ++c;
        tok = strtok(NULL, del);
    }
    while (c < V->c)
    {
        mpq_set_si(V->V[I(r, c, *V)], 0, 1);
        ++c;
    }
}

void crea_matriz(Matriz *V)
{
    char fila[MAX];
    printf("(Separa los %d elementos de cada renglon con comas, ej: 1/2, -3, 0.25)\n", V->c);
    for (int i = 0; i < V->r; ++i)
    {
        printf("Renglon %d: ", i + 1);
        if (fgets(fila, MAX, stdin) != NULL)
            ingresa_fraccion(fila, V, i);
    }
}

void crea_matriz_rand(Matriz *V)
{
    for (int i = 0; i < V->r * V->c; ++i)
    {
        long num = (rand() % 19) - 9;
        if (num == 0) num = 1;
        mpq_set_si(V->V[i], num, 1);
    }
}

// ============================================================================
// NORMAS DE MATRICES Y VECTORES
// ============================================================================

mpq_t *norma_1(Matriz A)
{
    if (A.r <= 0 || A.c <= 0 || A.V == NULL)
    {
        printf("Matriz vacia o no inicializada.\n");
        return NULL;
    }

    mpq_t *max_norma = (mpq_t *)malloc(sizeof(mpq_t));
    if (max_norma == NULL)
    {
        printf("No hay memoria disponible\n");
        exit(1);
    }
    mpq_init(*max_norma);

    mpq_t suma_col, abs_val;
    mpq_init(suma_col);
    mpq_init(abs_val);

    for (int j = 0; j < A.c; ++j)
    {
        mpq_set_ui(suma_col, 0, 1);
        for (int i = 0; i < A.r; ++i)
        {
            mpq_abs(abs_val, A.V[I(i, j, A)]);
            mpq_add(suma_col, suma_col, abs_val);
        }

        if (j == 0 || mpq_cmp(suma_col, *max_norma) > 0)
            mpq_set(*max_norma, suma_col);
    }

    mpq_clear(suma_col);
    mpq_clear(abs_val);
    return max_norma;
}

mpq_t *norma_inf(Matriz A)
{
    if (A.r <= 0 || A.c <= 0 || A.V == NULL)
    {
        printf("Matriz vacia o no inicializada.\n");
        return NULL;
    }

    mpq_t *max_norma = (mpq_t *)malloc(sizeof(mpq_t));
    if (max_norma == NULL)
    {
        printf("No hay memoria disponible\n");
        exit(1);
    }
    mpq_init(*max_norma);

    mpq_t suma_ren, abs_val;
    mpq_init(suma_ren);
    mpq_init(abs_val);

    for (int i = 0; i < A.r; ++i)
    {
        mpq_set_ui(suma_ren, 0, 1);
        for (int j = 0; j < A.c; ++j)
        {
            mpq_abs(abs_val, A.V[I(i, j, A)]);
            mpq_add(suma_ren, suma_ren, abs_val);
        }

        if (i == 0 || mpq_cmp(suma_ren, *max_norma) > 0)
            mpq_set(*max_norma, suma_ren);
    }

    mpq_clear(suma_ren);
    mpq_clear(abs_val);
    return max_norma;
}

double norma_euclidiana(Matriz A)
{
    if (A.r <= 0 || A.c <= 0 || A.V == NULL)
    {
        printf("Matriz vacia o no inicializada.\n");
        return -1.0;
    }

    mpq_t suma_cuad, cuad;
    mpq_init(suma_cuad);
    mpq_init(cuad);

    for (int i = 0; i < A.r * A.c; ++i)
    {
        mpq_mul(cuad, A.V[i], A.V[i]);
        mpq_add(suma_cuad, suma_cuad, cuad);
    }

    double res = sqrt(mpq_get_d(suma_cuad));

    mpq_clear(suma_cuad);
    mpq_clear(cuad);
    return res;
}

double norma_p(Matriz A, double p)
{
    if (A.r <= 0 || A.c <= 0 || A.V == NULL)
    {
        printf("Matriz vacia o no inicializada.\n");
        return -1.0;
    }
    if (p < 1.0)
    {
        printf("El parametro p debe ser mayor o igual a 1.\n");
        return -1.0;
    }

    double suma = 0.0;
    for (int i = 0; i < A.r * A.c; ++i)
    {
        double val = fabs(mpq_get_d(A.V[i]));
        suma += pow(val, p);
    }

    return pow(suma, 1.0 / p);
}

double norma_espectral(Matriz A)
{
    if (A.r <= 0 || A.c <= 0 || A.V == NULL)
    {
        printf("Matriz vacia o no inicializada.\n");
        return -1.0;
    }

    int n = A.c;
    double *B = (double *)calloc(n * n, sizeof(double));
    double *v = (double *)malloc(sizeof(double) * n);
    double *Bv = (double *)malloc(sizeof(double) * n);

    if (B == NULL || v == NULL || Bv == NULL)
    {
        printf("No hay memoria disponible\n");
        exit(1);
    }

    mpq_t acc, prod;
    mpq_init(acc);
    mpq_init(prod);

    int es_nula = 1;
    for (int i = 0; i < n; ++i)
    {
        for (int j = i; j < n; ++j)
        {
            mpq_set_ui(acc, 0, 1);
            for (int k = 0; k < A.r; ++k)
            {
                mpq_mul(prod, A.V[I(k, i, A)], A.V[I(k, j, A)]);
                mpq_add(acc, acc, prod);
            }
            double val = mpq_get_d(acc);
            B[i * n + j] = val;
            B[j * n + i] = val;
            if (val != 0.0)
                es_nula = 0;
        }
    }

    mpq_clear(acc);
    mpq_clear(prod);

    if (es_nula)
    {
        free(B);
        free(v);
        free(Bv);
        return 0.0;
    }

    double norm_v = 0.0;
    for (int i = 0; i < n; ++i)
    {
        v[i] = 1.0 + (double)i / (double)n;
        norm_v += v[i] * v[i];
    }
    norm_v = sqrt(norm_v);
    for (int i = 0; i < n; ++i)
        v[i] /= norm_v;

    double lambda = 0.0, lambda_prev = 0.0;
    const int max_iter = 1000;
    const double tol = 1e-12;

    for (int iter = 0; iter < max_iter; ++iter)
    {
        for (int i = 0; i < n; ++i)
        {
            double suma = 0.0;
            for (int j = 0; j < n; ++j)
                suma += B[i * n + j] * v[j];
            Bv[i] = suma;
        }

        lambda = 0.0;
        double norm_Bv = 0.0;
        for (int i = 0; i < n; ++i)
        {
            lambda += v[i] * Bv[i];
            norm_Bv += Bv[i] * Bv[i];
        }
        norm_Bv = sqrt(norm_Bv);

        if (norm_Bv == 0.0)
        {
            lambda = 0.0;
            break;
        }

        for (int i = 0; i < n; ++i)
            v[i] = Bv[i] / norm_Bv;

        if (fabs(lambda - lambda_prev) <= tol * (1.0 + fabs(lambda)))
            break;

        lambda_prev = lambda;
    }

    free(B);
    free(v);
    free(Bv);

    if (lambda < 0.0)
        lambda = 0.0;

    return sqrt(lambda);
}
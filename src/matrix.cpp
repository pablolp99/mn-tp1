#include <math.h>
#include <cstring>
#include <cstdio>
#include "matrix.hpp"
#include "furnace.hpp"

void print_matrix(double** A, int n, int m){
    for (int i = 0; i < n; ++i){
        for (int j = 0; j < m; ++j){
            printf("%f ", A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_vector(double* v, int n){
    for (int i = 0; i < n; ++i){
        printf("%f ", v[i]);
        printf("\n");
    }
    printf("\n");
}

double** create_2d_array(int rows, int cols){
    double** M = 0;
    M = new double*[rows];

    for (int i = 0; i < rows; ++i){
        M[i] = new double[cols];
        for (int j = 0; j < cols; ++j){
            M[i][j] = 0;
        }
    }

    return M;
}

double** calculate_coefficient_matrix(double r_i, double r_e, int m_1, int n){
    int size = m_1 * n;
    double delta_r = (r_e - r_i) / m_1;
    double delta_theta = 2 * M_PI / n;

    double** M = create_2d_array(size, size);

    // El for que esta en python
    for (int i = 1; i < m_1-1; ++i){
        for (int j = 0; j < n; ++j){
            double** v = create_2d_array(m_1, n);
            double r = r_i + delta_r * i;
            double tmp = calculate_t_j_kpn(r, delta_theta);

            v[i-1][j] = calculate_t_jp_k(delta_r, r);
            v[i][j] = calculate_t_ja_k(delta_r, r, delta_theta);
            v[i+1][j] = calculate_t_jn_k(delta_r);
            v[i][(j-1) % n] = tmp;
            v[i][(j+1) % n] = tmp;
            // Reshape v to be a vector
            double* reshaped_v = reshape_1d_array(v, size, m_1, n);
            // Make M[i*n + j] = reshaped_v
            M[i*n + j] = reshaped_v;
        }
    }

    double** M_r = create_2d_array(size - 2 * n, size - 2 * n);
    for (int i = n; i < (m_1 - 1) * n; ++i){
        for (int j = n; j < (m_1 - 1) * n; ++j){
            M_r[i - n][j - n] = M[i][j];
        }
    }

    return M_r;
}

double* reshape_1d_array(double** v, int size, int m, int n){
    double* reshaped;
    reshaped = new double[size];
    for (int i = 0; i < m; ++i){
        for (int j = 0; j < n; ++j){
            reshaped[i*n + j] = v[i][j];
        }
    }
    return reshaped;
}

double* calculate_b_vector(double* t_i, double* t_e, double** A, double r_i, double r_e, int n, int m_1){
    int size = n * (m_1 - 2);
    double* b = new double[size];
    double delta_r = (r_e - r_i) / m_1;

    for (int i = 0; i < n; ++i){
        b[i] = -A[n][0] * t_i[i];
        b[(size-1) - i] = -calculate_t_jn_k(delta_r) * t_e[(n-1) - i];
    }

    return b;
}

// Este metodo realiza una iteracion de la triangulacion (La iteracion i que pone en cero los elementos de la columna
// i desde la fila i+1 y actualiza el resto)
double* gaussian_elimination_step(double** A, double* b, int i, int n){
    // Vector de coeficientes (para posterior uso con la L)
    auto* coefs = new double[n];
    // Seteamos todos los valores iniciales en 0
    memset(coefs, 0, sizeof(coefs));
    for (int j = i+1; j < n; ++j) {
        // De coefs[0] a coefs[j-1] quedan con 0's
        coefs[j] = A[j][i] / A[i][i];
        A[j][i] = 0;
        for (int k = i+1; k < n; ++k) {
            A[j][k] = A[j][k] - coefs[j] * A[i][k];
        }
        b[j] = b[j] - coefs[j] * b[i];
    }
    return coefs;
}

void gaussian_elimination(double** A, double* b, int n){
    for (int i = 0; i < n; ++i){
        gaussian_elimination_step(A, b, i, n);
    }
}

double** LU_factorization(double** A, double* b, int n) {
    double** L = create_2d_array(n, n);
    double* coef;

    for (int i = 0; i < n; ++i){
        // Saves coeficients for L
        coef = gaussian_elimination_step(A, b, i, n);
        L[i][i] = 1;
        for (int j = i+1; j < n; ++j) {
            L[j][i] = coef[j];
        }
    }

    return L;
}

double* upper_triangular_system_solver(double** A, const double* b, int n) {
    auto* x = new double[n];

    memset(x, 0, n);
    for (int i = n-1; i >= 0; i--) {
        // x_i = b_i
        x[i] = b[i];
        for (int j = n-1; j > i; j--) { // x_i = b_i - sum_{j=i+1}^{n} a_{i,j} x_j
            x[i] = x[i] - A[i][j]*x[j];
        }
        // x_i = (b_i - sum_{j=i+1}^{n} a_{i,j} x_j) / a_{i,i}
        x[i] = x[i] / A[i][i];
    }
    return x;
}

double* lower_triangular_system_solver(double** A, const double* b, int n) {
    auto* x = new double[n];

    memset(x, 0, n);
    for (int i = 0; i < n; i++) {
        // x_i = b_i
        x[i] = b[i];
        for (int j = 0; j < i; j++) {
            x[i] = x[i] - A[i][j]*x[j];
        }
        // x_i = b_i - sum_{j=0}^{i} a_{i,j} x_j
        x[i] = x[i] / A[i][i];
        // x_i = (b_i - sum_{j=i+1}^{n} a_{i,j} x_j) / a_{i,i}
    }

    return x;
}

double* LU_resolver(double** L, double** U, double* b, int n) {
    double* y = lower_triangular_system_solver(L, b, n);

    return upper_triangular_system_solver(U, y, n);
}
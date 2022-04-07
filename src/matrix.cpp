#include <bits/stdc++.h>
#include <math.h>
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
            int index = j-1;
            if (j-1 < 0){
                index = n-1;
            }
            double r = r_i + (r_e - r_i) / (m_1 - 1);
            double tmp = calculate_t_j_kpn(r, delta_theta);

            v[i-1][j] = calculate_t_jp_k(delta_r, r);
            v[i][j] = calculate_t_ja_k(delta_r, r, delta_theta);
            v[i+1][j] = calculate_t_jn_k(delta_r);
            v[i][index] = tmp;
            v[i][(j+1) % n] = tmp;
            // Reshape v to be a vector
            double* reshaped_v = reshape_1d_array(v, size, m_1, n);
            // Make M[i*n + j] = reshaped_v
            M[i*n + j] = reshaped_v;
            // Habra que deletear los arrays?
            // delete reshaped_v;
            // delete v;
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

double* calculate_b_vector(double* t_i, double* t_e, double** A, double r_i, double r_e, int n, int m){
    int size = (n - 1) * (m - 1);
    double* b = new double[size];
    double delta_r = (r_e - r_i) / m;

    for (int i = 0; i < n; ++i){
        b[i] = -A[n][0] * t_i[i];
        b[(size-1) - i] = -calculate_t_jn_k(delta_r) * t_e[(n-1) - i];
    }

    return b;
}

double* gaussian_elimination_step(double** A, int i, int n){
    // Vector de coeficientes (para posterior uso con la L)
    double* coefs = new double[n];
    // Seteamos todos los valores iniciales en 0
    memset(coefs, 0, sizeof(coefs));
    for (int j = i+1; j < n; ++j) {
        // De coefs[0] a coefs[j-1] quedan con 0's
        coefs[j] = A[j][i] / A[i][i];
        for (int k = i; k < n; ++k) {
            A[j][k] = A[j][k] - coefs[j] * A[i][k];
        }
    }
    return coefs;
}

void gaussian_elimination(double** A, int n){
    for (int i = 0; i < n; ++i){
        gaussian_elimination_step(A, i, n);
    }
}
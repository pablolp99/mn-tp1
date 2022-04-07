#include <math.h>
#include "matrix.hpp"
#include "furnace.hpp"

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
    // 5 6 5 4 500 1
    int total = 0;
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
            total += 1;
        }
    }

    return M;
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
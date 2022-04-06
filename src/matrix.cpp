#include "matrix.hpp"
#include "furnace.hpp"

double** create_empty_matrix(int rows, int cols){
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

double** calculate_coefficient_matrix(double r, double delta_r, double delta_theta, int m_1, int n){
    int size = m_1 * n;
    double** M = create_empty_matrix(size, size);

    // El for que esta en python
    for (int i = 1; i < m_1-1; ++i){
        for (int j = 0; j < n; ++j){
            double** v = create_empty_matrix(m_1, n);

            v[i-1][j] = calculate_t_jp_k(delta_r, r);
            v[i][j] = calculate_t_ja_k(delta_r, r, delta_theta);
            v[i+1][j] = calculate_t_jn_k(delta_r);
            double tmp = calculate_t_j_kpn(r, delta_theta);
            v[i][j-1] = tmp;
            v[i][(j+1) % n] = tmp;
            // Reshape v to be a vector

            // Make M[i*n + j] = reshaped(v)

            delete v;
        }
    }

    return M;
}
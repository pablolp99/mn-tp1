#include <iostream>
#include "matrix.hpp"

using namespace std;

int main(){

    int m_1, n, n_inst;
    double r_i, r_e, iso;


    cin >> r_i >> r_e >> m_1 >> n >> iso >> n_inst;

    double t_i[n];
    double t_e[n];

    for (int i = 0; i < n_inst; ++i){
        double buff;
        for (int j = 0; j < n; ++j){
            cin >> buff;
            t_i[j] = buff;
        }
        for (int j = 0; j < n; ++j){
            cin >> buff;
            t_e[j] = buff;
        }
    }

    // The first approach is to use just one instance
    double** A = calculate_coefficient_matrix(r_i, r_e, m_1, n);
    // Calculate b vector
    double* b = calculate_b_vector(t_i, t_e, A, r_i, r_e, n, m_1);
    // Triangulate
    print_matrix(A, 12, 12);
    gaussian_elimination(A, 12);
    printf("\n");
    print_matrix(A, 12, 12);
    // Return x
    printf("Here\n");

    return 0;
}
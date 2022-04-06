#include <iostream>
#include "matrix.hpp"

using namespace std;

int main(){

    int m_1, n, n_inst;
    double r_i, r_e, iso;

    cin >> r_i >> r_e >> m_1 >> n >> iso >> n_inst;
    // The first approach is to use just one instance


    double** A = calculate_coefficient_matrix(r_i, r_e, m_1, n);
    // Calculate b vector
    // Triangulate
    // Return x
    printf("Here\n");

    return 0;
}
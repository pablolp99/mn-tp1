#include <iostream>
#include "matrix.hpp"
#include <math.h>

using namespace std;

int main(){

    int m_1, n, n_inst;
    double r_i, r_e, iso;

    cin >> r_i >> r_e >> m_1 >> n >> iso >> n_inst;
    // The first approach is to use just one instance
    double delta_r = (r_e - r_i) / m_1;
    double delta_theta = 2 * M_PI / n;

//    double A = calculate_coefficient_matrix();

    return 0;
}
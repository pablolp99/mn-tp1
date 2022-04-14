#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "matrix.hpp"

using namespace std;

#define GAUSSIAN 0
#define LU 1

int main(int argc, char *argv[]){
    ifstream infile("/home/pablo/UBA/comp2022/MN/mn-tp1/inputs/example_input.txt");
    string line;

    int m_1, n, n_inst;
    double r_i, r_e, iso;

    getline(infile, line);
    istringstream iss(line);
    iss >> r_i >> r_e >> m_1 >> n >> iso >> n_inst;

//    Arreglar esto para que soporte multiples instancias
    double t_i[n];
    double t_e[n];

    for (int i = 0; i < n_inst; ++i){
        getline(infile, line);
        double buff;
        istringstream iss(line);
        for (int j = 0; j < n; ++j){
            iss >> buff;
            t_i[j] = buff;
        }
        for (int j = 0; j < n; ++j){
            iss >> buff;
            t_e[j] = buff;
        }
    }

    int temperatures_amount = n * (m_1 - 2);

    // The first approach is to use just one instance
    double** U = calculate_coefficient_matrix(r_i, r_e, m_1, n);
    // Calculate b vector
    double* b = calculate_b_vector(t_i, t_e, U, r_i, r_e, n, m_1);

    // Read argv config
    // Triangulate
    double* x;
    int elimination = stoi(argv[1]);
    if (elimination == GAUSSIAN){
        gaussian_elimination(U, b, temperatures_amount);
        // Calculate vector X
        x = upper_triangular_system_solver(U, b, temperatures_amount);
    } else if (elimination == LU){
        double** L = LU_factorization(U, b, 12);
        // Calculate vector X
        x = LU_resolver(L, U, b, n);
    }

    // Interpolate results
    double delta_r = (r_e - r_i) / m_1;
    double** result = interpolate_results(x, n, n, m_1, r_i, delta_r, iso);

    print_matrix(result, n, 2);

    return 0;
}
#include <fstream>
#include <iostream>
#include <cmath>
#include <sstream>
#include <cstdio>
#include <string>
#include <chrono>

#include "matrix.hpp"

using namespace std;

#define GAUSSIAN 0
#define LU 1

int main(int argc, char *argv[]){
    // Read input file
    FILE * pFile;
    pFile = fopen(argv[2],"w");

    ifstream infile(argv[1]);
    string line;

    int m_1, n, n_inst;
    double r_i, r_e, iso;

    getline(infile, line);
    istringstream iss(line);
    iss >> r_i >> r_e >> m_1 >> n >> iso >> n_inst;

//    Arreglar esto para que soporte multiples instancias
    double** t_i;
    double** t_e;
    t_i = new double*[n_inst];
    t_e = new double*[n_inst];

    for (int i = 0; i < n_inst; ++i){
        getline(infile, line);
        istringstream iss(line);

        double buff;
        t_i[i] = new double[n];
        t_e[i] = new double[n];

        for (int j = 0; j < n; ++j){
            iss >> buff;
            t_i[i][j] = buff;
        }
        for (int j = 0; j < n; ++j){
            iss >> buff;
            t_e[i][j] = buff;
        }
    }

    int temperatures_amount = n * (m_1 - 2);

    // The first approach is to use just one instance
    double** U = calculate_coefficient_matrix(r_i, r_e, m_1, n);
    double** A = calculate_coefficient_matrix(r_i, r_e, m_1, n);
    for (int instance = 0; instance < n_inst; ++instance){
//        printf("Instance %d\n", instance);
        auto start = std::chrono::high_resolution_clock::now();

        // Calculate b vector
        double* b = calculate_b_vector(t_i[instance], t_e[instance], A, r_i, r_e, n, m_1);
        // Read argv config
        // Triangulate
        double* x;
        int elimination = stoi(argv[3]);

        if (elimination == GAUSSIAN){
            gaussian_elimination(U, b, temperatures_amount);
            // Calculate vector X
            x = upper_triangular_system_solver(U, b, temperatures_amount);
        } else if (elimination == LU){
            double **L;
            if(instance == 0) { //LU factorization has not been calculated
                L = LU_factorization(U, temperatures_amount);
            }
            // Calculate vector X
            x = LU_resolver(L, U, b, temperatures_amount);
        } else {
            printf("Elimination method not supported. Please choose 0 (Gaussian Elmination) or 1 (LU Factorization)\n");
            return 1;
        }
        auto stop = std::chrono::high_resolution_clock::now();

        // Interpolate results
        double delta_theta = 2 * M_PI / n;
        double delta_r = (r_e - r_i) / m_1;
        double** result = interpolate_results(x, n, n, m_1, r_i, delta_r, delta_theta, iso);
        print_matrix(result, n, 3, pFile);

        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);

//        cout <<  << endl;
        printf("%lu\n", duration.count());
    }

    return 0;
}
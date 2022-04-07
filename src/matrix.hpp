double** create_2d_array(int rows, int cols);
double** calculate_coefficient_matrix(double r_i, double r_e, int m_1, int n);
double* reshape_1d_array(double** v, int size, int m, int n);
double* calculate_b_vector(double* t_i, double* t_e, double** A, double r_i, double r_e, int n, int m);
void lu_factorization(double** a, double** l, double** u, int n);
double* gaussian_elimination_step(double** A, int i, int n);
void gaussian_elimination(double** A, int n);
#include <math.h>
#include "furnace.hpp"

double calculate_t_jp_k(double delta_r, double r){
    // Calculates t_{j-1, k} coeficients
    return pow(delta_r, -2) - pow((r * delta_r), -1);
}

double calculate_t_ja_k(double delta_r, double r, double delta_theta){
    // Calculates t_{j, k} coeficients
    return pow((r * delta_r), -1) - 2 * pow(delta_r, -2) - 2 * pow(pow( (r * delta_theta), 2 ), -1);
}

double calculate_t_jn_k(double delta_r){
    // Calculates t_{j+1, k} coeficients
    return pow(delta_r, -2);
}

double calculate_t_j_kpn(double r, double delta_theta){
    // Calculates t_{j, k-1} and t_{j, k+1} coeficients
    return pow(pow(r * delta_theta, 2), -1);
}

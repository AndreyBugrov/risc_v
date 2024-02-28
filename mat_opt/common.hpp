#pragma once
#include <random>  // random matrix generation
#include <cstring> // memset in zero matrix generating

typedef void (*mult_func)(double*, double*, double*, int);

void generate_rand_matrix(double* matr, int n, double min, double max);
void generate_zero_matrix(double* matr, int n);
void generate_identity_matrix(double* matr, int n);

#pragma once
#include "../open_blas/include/cblas.h"
void simple_matrix_mult(double* a, double* b, double* c, int n_a, int m_b, int elements_in_vector);
void matrix_mult_second_transposed(double* a, double* b, double* c, int n_a, int m_b, int elements_in_vector);
void block_matrix_mult(double* a, double* b, double* c, int n_a, int m_b, int elements_in_vector);
void block_matrix_mult_second_transposed(double* a, double* b, double* c, int n_a, int m_b, int elements_in_vector);
void transpose_matrix(double* matr, int n, int m);
void block_transpose_matrix(double* matr, int n, int m);
#pragma once
#include <omp.h>    // parallel versions 
#include <iostream> // std::swap for

void base_matrix_mult(double* a, double* b, double* c, int n_a, int m_b, int elements_in_vector);
void base_matrix_mult_omp(double* a, double* b, double* c, int n_a, int m_b, int elements_in_vector);

void b_transposed_matrix_mult(double* a, double* b, double* c, int n_a, int m_b, int elements_in_vector);
void b_transposed_matrix_mult_omp(double* a, double* b, double* c, int n_a, int m_b, int elements_in_vector);

void base_block_matrix_mult(double* a, double* b, double* c, int n_a, int m_b, int elements_in_vector, int block_size);
void optimal_block_matrix_mult(double* a, double* b, double* c, int n_a, int m_b, int elements_in_vector);
void b_transposed_block_matrix_mult(double* a, double* b, double* c, int n_a, int m_b, int elements_in_vector);

void strassen_matrix_mult(double* a, double* b, double* c, int n_a, int m_b, int elements_in_vector);

void transpose_square_matrix(double* matr, int n, int m);
void transpose_common_matrix(double* matr, int n, int m);
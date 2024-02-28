#pragma once

#include <omp.h>    // parallel versions

#include <utility> // std::swap in transposision
#include <cstring>  // std::memcpy in transposition

void base_matrix_mult(double* a, double* b, double* c, int n);
void base_matrix_mult_omp(double* a, double* b, double* c, int n);

void b_transposed_matrix_mult(double* a, double* b, double* c, int n);
void b_transposed_matrix_mult_omp(double* a, double* b, double* c, int n);

void transposed_matrix_mult(double* a, double* b, double* c, int n);
void transposed_matrix_mult_omp(double* a, double* b, double* c, int n);

void base_block_matrix_mult(double* a, double* b, double* c, int n, int block_size);
void optimal_block_matrix_mult(double* a, double* b, double* c, int n);
void b_transposed_block_matrix_mult(double* a, double* b, double* c, int n);

void strassen_matrix_mult(double* a, double* b, double* c, int n);

void transpose_matrix_in_place(double* matr, int n);
double* get_transposed_matrix(double* matr, int n);
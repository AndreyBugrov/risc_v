#pragma once

#include <omp.h>    // parallel versions

#include <utility>  // std::swap in transposision
#include <cstring>  // std::memcpy in transposition & increase_matrixes
#include <cmath>    // log2 in base_block_matrix_mult
#include <vector>   // std:vector in increase_matrixes
#include <numeric>  // std::inner_product in row_matrix_mult

#include "common.hpp" // generate_zero_matrix in increase_matrixes

void base_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);
void base_matrix_mult_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);

void row_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);
void row_matrix_mult_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);

void b_transposed_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);
void b_transposed_matrix_mult_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);
void b_transposed_matrix_mult_omp_simd(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);

void transposed_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);
void transposed_matrix_mult_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);
void transposed_matrix_mult_omp_simd(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);

// unsafe for matrices with odd sides 
void base_block_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);
void optimal_block_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);
void b_transposed_block_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);

void strassen_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);
void recursive_strassen_part(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);

void transpose_matrix_in_place(double* matr, int n);
double* get_transposed_matrix(double* matr, int n);

// n - length of one matrix dimension, block_n - length of one block dimension 
void blockcpy(double* __restrict__ src, double* __restrict__ dest, int n, int block_n);

//make matrix dimension lengths power of 2, return new matrix dimension length
int increase_matrices(double* __restrict__ a, double* __restrict__ inc_a, double* __restrict__ b, double* __restrict__ inc_b, int n);
void split_matrices(double* __restrict__ a, double* __restrict__ a11, double* __restrict__ a12,double* __restrict__ a21, double* __restrict__ a22, int n);
void collect_matrices(double* __restrict__ a, double* __restrict__ a11, double* __restrict__ a12,double* __restrict__ a21, double* __restrict__ a22, int n);
void matrix_add(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);
void matrix_sub(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);
#pragma once

#include <omp.h>    // parallel versions

#include <utility>  // std::swap in transposision
#include <cstring>  // std::memcpy in transposition
#include <cmath>    // log2 in base_block_matrix_mult

void base_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);
void base_matrix_mult_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);
void base_matrix_mult_omp_simd(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);

void b_transposed_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);
void b_transposed_matrix_mult_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);

void transposed_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);
void transposed_matrix_mult_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);

// unsafe for matrices with odd sides 
void base_block_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);
void optimal_block_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);
void b_transposed_block_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);

void strassen_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n);

void transpose_matrix_in_place(double* matr, int n);
double* get_transposed_matrix(double* matr, int n);

// n - length of one matrix dimension, block_n - length of one block dimension 
void blockcpy(double* __restrict__ src, double* __restrict__ dest, int n, int block_n);

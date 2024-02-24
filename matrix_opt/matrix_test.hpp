#pragma once

#include <float.h> // DBL_MAX

#include <iostream> // test result printing
#include <chrono> // time measurment
#include <vector> // get_unequal_elements
#include <random> // random tests

#include "../open_blas/include/cblas.h" // reference version

#include "all.hpp" // multiplication and transpose from "mult_types.hpp", matrix generation from "common.hpp"


// void cblas_domatcopy(OPENBLAS_CONST enum CBLAS_ORDER CORDER, OPENBLAS_CONST enum CBLAS_TRANSPOSE CTRANS, OPENBLAS_CONST blasint crows, OPENBLAS_CONST blasint ccols, OPENBLAS_CONST double calpha, OPENBLAS_CONST double *a,
// 		     OPENBLAS_CONST blasint clda, double *b, OPENBLAS_CONST blasint cldb);

void fixed_multiplication_test(void (*matrix_mult_function)(double*, double*, double*, int, int, int)); // fixed data test for matrix multiplication
void open_blas_multiplication_test(void (*matrix_mult_function)(double*, double*, double*, int, int, int)); // random big data test for matrix multiplication

void fixed_transpose_matrix_test(void (*transpose_matrix_function)(double*, int, int)); // fixed data test for matrix transpose
void open_blas_transpose_matrix_test(void (*transpose_matrix_function)(double*, int, int)); // random big data test for matrix transpose

void fill_in_matrixes_to_multiply(double* a, double* b, double* c);
void fill_in_matrix_to_transpose(double* a);
// return two pairs: 1 - pair of unequal element number and first met different unequal element index 2 - pair of minimum and maximum difference
// return percentage of unequal elements, minimum and maximum difference
std::vector<double> get_unequal_elements(double* base, double* current, int n);
void print_test_result(std::vector<double> comparison_output, double seconds);

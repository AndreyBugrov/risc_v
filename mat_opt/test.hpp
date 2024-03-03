#pragma once

#include <float.h>                      // DBL_MAX in get_unequal_elements

#include <iostream>                     // print_test_result
#include <chrono>                       // time measurment in tests
#include <vector>                       // get_unequal_elements
#include <cmath>                        // abs in print_test_result
#include <string>                       // argument parsing in main
#include <map>                          // argument parsing in main

#include "../open_blas/include/cblas.h" // reference version

#include "multiplication.hpp"           // multiplication functions
#include "common.hpp"                   // matrix generation in tests


// void cblas_domatcopy(OPENBLAS_CONST enum CBLAS_ORDER CORDER, OPENBLAS_CONST enum CBLAS_TRANSPOSE CTRANS, OPENBLAS_CONST blasint crows, OPENBLAS_CONST blasint ccols, OPENBLAS_CONST double calpha, OPENBLAS_CONST double *a,
// 		     OPENBLAS_CONST blasint clda, double *b, OPENBLAS_CONST blasint cldb);

enum class test_type{
    zero,
    identity,
    equal,
    random,
    big
};

bool multiplication_test(std::string test_name, mult_func matrix_mult_function, test_type type);

// return two pairs: 1 - pair of unequal element number and first met different unequal element index 2 - pair of minimum and maximum difference
// return percentage of unequal elements, minimum and maximum difference
std::vector<double> get_unequal_elements(double* base, double* current, int size);
// return is test passed or not
bool print_test_result(std::vector<double> comparison_output, double seconds);
void print_test_statistics(int passed_num, int failed_num, double seconds);

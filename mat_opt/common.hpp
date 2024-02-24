#pragma once
#include <iostream>  // save_result
#include <random>    // random matrix generation
#include <cmath>     // abs
#include <string>    // function sets
#include <map>       // function name matching
#include <stdexcept> // exceptions in map.at()

#include "all.hpp"

typedef void (*mult_func)(double*, double*, double*, int, int, int);
typedef void (*tr_func)(double*, int, int);

void generate_rand_matrix(double* matr, int n, int m, double min, double max);
void generate_zero_matrix(double* matr, int n, int m);
int get_max_value_index(double* vec, int n);
void save_result(double* total_seconds, int exp_num);

mult_func set_multiplication_function(std::string function_name);
tr_func set_transpose_function(std::string function_name);

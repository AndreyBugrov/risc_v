#include <iostream>                     // std::cout in print_result
#include <map>                          // function name matching in get_multiplication_function
#include <stdexcept>                    // exceptions in map.at() in get_multiplication_function
#include <chrono>                       // time measurment in main

#include "../open_blas/include/cblas.h" // reference version for main

#include "common.hpp"                   // multiplication types in get_multiplication_function
#include "multiplication.hpp"           // multiplication functions in get_multiplication_function

int get_max_value_index(double* vec, int n);
void print_result(double* total_seconds, int exp_num);

mult_func get_multiplication_function(std::string function_name);
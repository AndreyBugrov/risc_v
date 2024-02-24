#pragma once
#include <iostream>
#include <random> // random matrix generation
#include <cmath> // abs

void generate_rand_matrix(double* matr, int n, int m, double min, double max);
void generate_zero_matrix(double* matr, int n, int m);
void print_result(double seconds, bool is_automatic);
int get_max_value_index(double* vec, int n);
void save_result(double* total_seconds, int exp_num, bool is_automatic);

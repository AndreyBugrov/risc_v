#pragma once
#include <string>
#include <iostream>
#include <random>
#include <cmath> // abs
using std::string;
void transpose_matrix(double* matr, int n, int m);
void generate_rand_matrix(double* matr, int n, int m, double min, double max);
void generate_zero_matrix(double* matr, int n, int m);
void generate_test_matrixes(double* a, double* b, double* c);
bool check_test_result(double* result_matrix);
std::string print_test_result(double* result_matrix);
void print_result(double seconds, bool is_automatic);
int get_max_value_index(double* vec, int n);
void save_result(double* total_seconds, int exp_num, bool is_automatic);
bool are_vectors_equal(double* a, double* b, int n);
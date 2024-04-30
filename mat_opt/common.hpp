#pragma once
#include <random>  // random matrix generation
#include <cstring> // memset in zero matrix generating
#include <fstream> // file reading in get_cache_sizes
#include <string>  // file reading in get_cache_sizes

const int kCacheLevelNumber = 15;

typedef void (*mult_func)(double*, double*, double*, int);

void generate_rand_matrix(double* matr, int n, double min, double max);
void generate_zero_matrix(double* matr, int n);
void generate_identity_matrix(double* matr, int n);

void get_cache_sizes(std::string filename, int* cache_sizes);

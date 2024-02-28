#include "multiplication.hpp"

void base_block_matrix_mult(double* a, double* b, double* c, int n){}
void optimal_block_matrix_mult(double* a, double* b, double* c, int n){}
void b_transposed_block_matrix_mult(double* a, double* b, double* c, int n){}

void strassen_matrix_mult(double* a, double* b, double* c, int n){}

void base_matrix_mult(double* a, double* b, double* c, int n){
    for(int i=0;i<n;i++){ // i-th row in a
        for(int j=0;j<n;j++){ // j-th column in b
            for(int k=0;k<n; k++){ // k-th element in vector
                c[i*n+j]+=a[i*n+k]*b[k*n+j];
            }
        }
    }
}
void base_matrix_mult_omp(double* a, double* b, double* c, int n){
    #pragma omp parallel for shared(a, b, c, n)
        for(int i=0;i<n;i++){ // i-th row in a
            for(int j=0;j<n;j++){ // j-th column in b
                for(int k=0;k<n; k++){ // k-th element in vector
                    c[i*n+j]+=a[i*n+k]*b[k*n+j];
                }
            }
        }
}

void b_transposed_matrix_mult(double* a, double* b, double* c, int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            for(int k=0;k<n;k++){
                c[i*n+j] += a[i*n+k]*b[j*n+k];
            }
        }
    }
}
void b_transposed_matrix_mult_omp(double* a, double* b, double* c, int n){
    #pragma omp parallel for shared(a, b, c, n)
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                for(int k=0;k<n;k++){
                    c[i*n+j] += a[i*n+k]*b[j*n+k];
                }
            }
        }
}
void transposed_matrix_mult(double* a, double* b, double* c, int n){
    double* bT = get_transposed_matrix(b, n);
    b_transposed_matrix_mult(a, bT, c, n);
}
void transposed_matrix_mult_omp(double* a, double* b, double* c, int n){
    double* bT = get_transposed_matrix(b, n);
    b_transposed_matrix_mult_omp(a, bT, c, n);
}

void transpose_matrix_in_place(double* matr, int n){
    for(int i=0;i<n;i++){
        for(int j=i;j<n;j++){
            std::swap(matr[i * n + j], matr[j * n + i]);
        }
    }
}
double* get_transposed_matrix(double* matrix, int n){
    double* Tmatrix = new double[n*n];
    std::memcpy(Tmatrix, matrix, n*n*sizeof(double));
    transpose_matrix_in_place(Tmatrix, n);
    return Tmatrix;
}
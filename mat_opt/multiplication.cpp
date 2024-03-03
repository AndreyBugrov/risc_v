#include "multiplication.hpp"

void base_block_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    int block_n;
    int tmp_block_n = n;
    while(tmp_block_n%2==0){
        tmp_block_n /= 2;
    }
    if(tmp_block_n>1){
        block_n = tmp_block_n;
    }else{
        block_n = 2; // unsafe part
    }

    double* __restrict__ a_block = new double[n*n];
    double* __restrict__ b_block = new double[n*n];
    blockcpy(a, a_block, n, block_n);
    blockcpy(b, b_block, n, block_n);
    //////////////// multiplication itself
}
void optimal_block_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){}
void b_transposed_block_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){}

void strassen_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){}

void base_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    for(int i=0;i<n;i++){ // i-th row in a
        for(int j=0;j<n;j++){ // j-th column in b
            for(int k=0;k<n; k++){ // k-th element in vector
                c[i*n+j]+=a[i*n+k]*b[k*n+j];
            }
        }
    }
}
void base_matrix_mult_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    #pragma omp parallel for shared(a, b, c, n)
        for(int i=0;i<n;i++){ // i-th row in a
            for(int j=0;j<n;j++){ // j-th column in b
                for(int k=0;k<n; k++){ // k-th element in vector
                    c[i*n+j]+=a[i*n+k]*b[k*n+j];
                }
            }
        }
}
void base_matrix_mult_omp_simd(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    #pragma omp parallel for shared(a, b, c, n)
    for(int i=0;i<n;i++){ // i-th row in a
        for(int j=0;j<n;j++){ // j-th column in b
            #pragma omp simd
                for(int k=0;k<n; k++){ // k-th element in vector
                    c[i*n+j]+=a[i*n+k]*b[k*n+j];
                }
        }
    }
}

void b_transposed_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            for(int k=0;k<n;k++){
                c[i*n+j] += a[i*n+k]*b[j*n+k];
            }
        }
    }
}
void b_transposed_matrix_mult_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    #pragma omp parallel for shared(a, b, c, n)
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                for(int k=0;k<n;k++){
                    c[i*n+j] += a[i*n+k]*b[j*n+k];
                }
            }
        }
}
void transposed_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    double* bT = get_transposed_matrix(b, n);
    b_transposed_matrix_mult(a, bT, c, n);
}
void transposed_matrix_mult_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
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
void blockcpy(double* __restrict__ src, double* __restrict__ dest, int n, int block_n){

    const int block_line = n / block_n; // block number in src matrix line
    const int block_size = block_n * block_n; // element number in block
    const int block_line_size = block_n * block_n * block_line; // element number in block line

    for(int i=0;i<block_line;i++){ // i-th block row
        for(int j=0;j<block_line;j++){ // j-th block column
            for(int k=0;k<block_n;k++){ // k-th row in block
                std::memcpy(dest+i*block_line_size+j*block_size+k*block_n, src+i*n*block_n+k*n+j*block_n, block_n*sizeof(double));
            }
        }
    }
}
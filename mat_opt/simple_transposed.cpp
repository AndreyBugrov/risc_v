#include "mult_types.hpp"

void transpose_matrix(double* matr, int n, int m){
    double* tmp_matr = new double[n*m];
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            tmp_matr[i*m+j]=matr[j*m+i];
        }
    }
    for(int i=0;i<n*m;i++){
        matr[i]=tmp_matr[i];
    }
    delete[] tmp_matr;
}

void matrix_mult_second_transposed(double* a, double* b, double* c, int n_a, int n_b, int elements_in_vector){
    for(int i=0;i<n_a;i++){
        for(int j=0;j<n_b;j++){
            for(int k=0;k<elements_in_vector;k++){
                c[i*n_b+j]=a[i*elements_in_vector+k]*b[j*elements_in_vector+k];
            }
        }
    }
}
#include "mult_types.hpp"

void transpose_square_matrix(double* matr, int n, int m){
    for(int i=0;i<n;i++){
        for(int j=i;j<m;j++){
            std::swap(matr[i * m + j], matr[j * n + i]);
        }
    }
}

void transpose_common_matrix(double* matr, int n, int m){
    double* tmp_matr = new double[m*n];
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            tmp_matr[j*n+i]=matr[i*m+j];
        }
    }
    for(int i=0;i<n*m;i++){
        matr[i]=tmp_matr[i];
    }
    delete[] tmp_matr;
}
#include "mult_types.hpp"

void base_matrix_mult(double* a, double* b, double* c, int n_a, int m_b, int elements_in_vector){
    for(int i=0;i<n_a;i++){ // i-th row in a
        for(int j=0;j<m_b;j++){ // j-th column in b
            for(int k=0;k<elements_in_vector; k++){ // k-th element in vector
                c[i*m_b+j]+=a[i*elements_in_vector+k]*b[k*m_b+j];
            }
        }
    }
}

void base_matrix_mult_omp(double* a, double* b, double* c, int n_a, int m_b, int elements_in_vector){
    #pragma omp parallel for shared(a, b, c, n_a, m_b, elements_in_vector)
        for(int i=0;i<n_a;i++){ // i-th row in a
            for(int j=0;j<m_b;j++){ // j-th column in b
                for(int k=0;k<elements_in_vector; k++){ // k-th element in vector
                    c[i*m_b+j]+=a[i*elements_in_vector+k]*b[k*m_b+j];
                }
            }
        }
}
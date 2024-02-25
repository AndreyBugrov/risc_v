#include "mult_types.hpp"

void b_transposed_matrix_mult(double* a, double* b, double* c, int n_a, int m_b, int elements_in_vector){
    for(int i=0;i<n_a;i++){
        for(int j=0;j<m_b;j++){
            for(int k=0;k<elements_in_vector;k++){
                c[i*m_b+j] += a[i*elements_in_vector+k]*b[j*elements_in_vector+k];
            }
        }
    }
}

void b_transposed_matrix_mult_omp(double* a, double* b, double* c, int n_a, int m_b, int elements_in_vector){
    #pragma omp parallel for shared(a, b, c, n_a, m_b, elements_in_vector)
        for(int i=0;i<n_a;i++){
            for(int j=0;j<m_b;j++){
                for(int k=0;k<elements_in_vector;k++){
                    c[i*m_b+j] += a[i*elements_in_vector+k]*b[j*elements_in_vector+k];
                }
            }
        }
}

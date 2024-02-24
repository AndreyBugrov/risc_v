#include <chrono> // time measurment

#include "all.hpp" // matrix multiplication and generation
#include "../open_blas/include/cblas.h"

int main(int argc, char* argv[]){
    std::string function_name = std::string(argv[1]);

    // n x m matrixes:
    // a: n_a x elements_in_vector
    // b: elements_in_vector x m_b
    // c: n_a * m_b
    const int matrix_size = argc>2? std::stoi(argv[2]) : 1000;
    const int n_a = matrix_size;
    const int elements_in_vector = matrix_size;
    const int m_b = matrix_size;
    const int exp_num = argc>3? std::stoi(argv[3]) : 1;
    double* a, *b, *c;
    double* total_seconds = new double[exp_num];
    double* dgemm_seconds = new double[exp_num];

    a = new double[n_a*elements_in_vector];
    b = new double[elements_in_vector*m_b];
    c = new double[n_a*m_b];
    double* c_blas = new double[n_a*m_b];

    generate_rand_matrix(a, n_a, elements_in_vector, 0.0, 100.0);
    generate_rand_matrix(b, elements_in_vector, m_b, -50.0, 50.0);

    mult_func matrix_mult_function = set_multiplication_function(function_name);

    for(int i=0;i<exp_num;i++){
        generate_zero_matrix(c, n_a, m_b);

        const auto start_my_mult{std::chrono::steady_clock::now()};
        matrix_mult_function(a, b, c, n_a, m_b, elements_in_vector);
        const auto end_my_mult{std::chrono::steady_clock::now()};

        const auto start_dgemm{std::chrono::steady_clock::now()};
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, n_a, m_b, elements_in_vector, 1.0, a, elements_in_vector, b, m_b, 0.0, c_blas, m_b);
        const auto end_dgemm{std::chrono::steady_clock::now()};
        
        std::chrono::duration<double> elapsed_seconds = end_my_mult - start_my_mult;
        total_seconds[i] = elapsed_seconds.count();
        std::chrono::duration<double> elapsed_dgemm = end_dgemm - start_dgemm;
        dgemm_seconds[i] = elapsed_dgemm.count();
    }
    save_result(total_seconds, exp_num);
    save_result(dgemm_seconds, exp_num);
    
    delete[] a; 
    delete[] b;
    delete[] c;
    return 0;
}

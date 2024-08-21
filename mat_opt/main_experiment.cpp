#include "experiment.hpp"


int main(int argc, char* argv[]){
    const bool is_automatic = std::string(argv[1]) == std::string("a") ? true : false; // it is strange but argv[1] != "a" && argv[1] != "a\0"
    std::string function_name = std::string(argv[2]);

    // n x m matrixes:
    // a: n_a x elements_in_vector
    // b: elements_in_vector x m_b
    // c: n_a * m_b
    const int n = argc>3 ? std::stoi(argv[3]) : 1000;
    const int exp_num = argc>4 ? std::stoi(argv[4]) : 1;
    double* a, *b, *c;
    double* total_seconds = new double[exp_num];
    double* dgemm_seconds = new double[exp_num];

    a = new double[n*n];
    b = new double[n*n];
    c = new double[n*n];
    double* c_blas = new double[n*n];

    generate_rand_matrix(a, n, 0.0, 100.0);
    generate_rand_matrix(b, n, -50.0, 50.0);

    mult_func matrix_mult_function = get_multiplication_function(function_name);

    for(int i=0;i<exp_num;i++){
        generate_zero_matrix(c, n);
        generate_zero_matrix(c_blas, n);

        const auto start_my_mult{std::chrono::steady_clock::now()};
        matrix_mult_function(a, b, c, n);
        const auto end_my_mult{std::chrono::steady_clock::now()};

        const auto start_dgemm{std::chrono::steady_clock::now()};
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, n, n, n, 1.0, a, n, b, n, 0.0, c_blas, n);
        const auto end_dgemm{std::chrono::steady_clock::now()};
        
        std::chrono::duration<double> elapsed_seconds = end_my_mult - start_my_mult;
        total_seconds[i] = elapsed_seconds.count();

        std::chrono::duration<double> elapsed_dgemm = end_dgemm - start_dgemm;
        dgemm_seconds[i] = elapsed_dgemm.count();
    }

    print_experiment_result(dgemm_seconds, total_seconds, exp_num, c_blas, c, n*n, is_automatic);

    delete[] a;
    delete[] b;
    delete[] c;

    return 0;
}

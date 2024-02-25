#include <iostream>
#include <random>
#include <string> // stoi
#include <chrono> // time
#include <cmath> // abs

#include "open_blas/include/cblas.h"

const int test_element_num = 12;
const double etalon[test_element_num]={    
    5, -1, 4, -1, 
    -1, 2, 1, 2, 
    3, 0, 3, 0};


void simple_matrix_mult(double* a, double* b, double* c, int n_a, int m_b, int elements_in_vector){
    for(int i=0;i<n_a;i++){ // i-th row in a
        for(int j=0;j<m_b;j++){ // j-th column in b
            for(int k=0;k<elements_in_vector; k++){ // k-th element in vector
                c[i*m_b+j]+=a[i*elements_in_vector+k]*b[k*m_b+j];
            }
        }
    }
}

void matrix_mult_second_transposed(double* a, double* b, double* c, int n_a, int m_b, int elements_in_vector){
        #pragma omp parallel for shared(a, b, c, n_a, m_b, elements_in_vector)
    for(int i=0;i<n_a;i++){
        for(int j=0;j<m_b;j++){
            for(int k=0;k<elements_in_vector;k++){
                c[i*m_b+j] += a[i*elements_in_vector+k]*b[j*elements_in_vector+k];
                // std::cout<<"c="<<c[i*n_b+j]<<"\n";
            }
        }
        // std::cout<<"i="<<i<<"\n";
    }
}

void transpose_matrix(double* matr, int n, int m){
    // double* tmp_matr = new double[m*n];
    // for(int i=0;i<n;i++){
    //     for(int j=0;j<m;j++){
    //         tmp_matr[j*n+i]=matr[i*m+j];
    //     }
    // }
    // for(int i=0;i<n*m;i++){
    //     matr[i]=tmp_matr[i];
    // }
    // delete[] tmp_matr;

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

    // double i_j, j_i;
    // for(int i=0;i<n;i++){
    //     for(int j=i+1;j<m;j++){
    //         i_j = matr[i*m+j];
    //         j_i = matr[j*n+i];
    //         matr[i*m+j]=j_i;
    //         matr[j*n+i]=i_j;
    //     }
    // }

    // for(int i=0;i<n;i++){
    //     for(int j=i;j<m;j++){
    //         std::swap(matr[i * m + j], matr[j * n + i]);
    //     }
    // }
}

void memory_transpose_square_matrix(double* matr, int n, int m){
    for(int i=0;i<n;i++){
        for(int j=i;j<m;j++){
            std::swap(matr[i * m + j], matr[j * n + i]);
        }
    }
}

void generate_rand_matrix(double* matr, int n, int m, double min, double max){
    std::random_device rd;
    std::mt19937 engine(rd());
    std::uniform_real_distribution<double> gen(min, max);
    for(int i=0;i<n*m;i++){
        matr[i] = gen(engine);
    }
}
void generate_zero_matrix(double* matr, int n, int m){
    for(int i=0;i<n*m;i++){
        matr[i]=0.0;
    }
}
void generate_test_matrixes(double* a, double* b, double* c){
    int n_a = 3;
    int elements_in_vector = 2;
    int m_b = 4;
    a[0]=1;
    a[1] = 2;
    a[2] = 1;
    a[3] = -1;
    a[4] = a[5]= 1;

    b[0] = 1;
    b[1] = 1;
    b[2] = 2;
    b[3] = 1;
    b[4] = 2;
    b[5] = -1;
    b[6] = 1;
    b[7] = -1;
    for(int i=0;i<n_a*m_b;i++){
        c[i]=0.0;
    }
    /*
    c= 
    5 -1 4 -1 
    -1 2 1 2 
    3 0 3 0
    */
}
bool check_test_result(double* result_matrix){
    for(int i=0;i<test_element_num;i++){
        if(result_matrix[i]!=etalon[i]){
            return false;
        }
    }
    return true;
}
std::string print_test_result(double* result_matrix){
    return check_test_result(result_matrix)? "Test passed": "Test failed";
}
void print_result(double seconds, bool is_automatic){
    // if(is_automatic){
    //     std::cout<<seconds;
    // }
    // else{
    //     std::cout<<seconds<<"\n";
    // }
    std::cout<<seconds<<"\n";
}

int get_max_value_index(double* vec, int n){
    double max_value = vec[0];
    int max_index = 0;
    for(int i=1;i<n;i++){
        if(max_value<vec[i]){
            max_value = vec[i];
            max_index = i;
        }
    }
    return max_index;
}

void save_result(double* total_seconds, int exp_num, bool is_automatic){
  if(exp_num==1){
    print_result(total_seconds[0], is_automatic);
  }
  else{
    int max_index = get_max_value_index(total_seconds, exp_num);
    double seconds_without_outliers = 0;
    for(int i=0;i<exp_num;i++){
        if(i!=max_index){
            seconds_without_outliers+=total_seconds[i];
        }
    }
    print_result(seconds_without_outliers / (exp_num-1), is_automatic);
  }
}

bool are_vectors_equal(double* a, double* b, int n){
    for(int i=0; i<n; i++){
        if(abs(a[i]-b[i])>0.0001)
            {
                std::cout<<"abs = "<<abs(a[i]-b[i])<<"\n";
            return false;
            }
    }
    return true;
}


int main(int argc, char* argv[]){
    const bool is_automatic = std::string(argv[1]) == std::string("a") ? true : false; // it is strange but argv[1] != "a" && argv[1] != "a\0"
    // n x m matrixes:
    //a: n_a x elements_in_vector
    //b: elements_in_vector x m_b
    const int matrix_size = argc>2? std::stoi(argv[2]) : 1000;
    const int n_a = matrix_size;
    const int elements_in_vector = matrix_size;
    const int m_b = matrix_size;
    const int exp_num = argc>3? std::stoi(argv[3]) : 1;
    double* a, *b, *c;
    double* total_seconds = new double[exp_num];
    double* dgemm_seconds = new double[exp_num];

    double* base_b = new double[elements_in_vector*m_b];

    a = new double[n_a*elements_in_vector];
    b = new double[elements_in_vector*m_b];
    c = new double[n_a*m_b];
    double* c_blas = new double[n_a*m_b];

    generate_rand_matrix(a,n_a,elements_in_vector,0.0,10);
    generate_rand_matrix(b,elements_in_vector,m_b,-5.0,5.0);

    for(int i=0;i<elements_in_vector*m_b;i++){
        base_b[i]=b[i];
    }

    for(int i=0;i<exp_num;i++){
        generate_zero_matrix(c,n_a,m_b);
        memory_transpose_square_matrix(b, elements_in_vector, m_b);
        const auto start_simple{std::chrono::steady_clock::now()};
        // simple_matrix_mult(a,b,c,n_a,m_b,elements_in_vector);
        matrix_mult_second_transposed(a,b,c,n_a,m_b,elements_in_vector);
        const auto end_simple{std::chrono::steady_clock::now()};
        const auto start_dgemm{std::chrono::steady_clock::now()};

        memory_transpose_square_matrix(b, elements_in_vector, m_b);
        cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, n_a, m_b, elements_in_vector, 1.0, a, elements_in_vector, b, m_b,0.0,c_blas,m_b);
        const auto end_dgemm{std::chrono::steady_clock::now()};

        std::cout<<are_vectors_equal(c, c_blas, n_a*m_b)<<"\n";

        std::cout<<"b and base_b: "<<are_vectors_equal(base_b, b, elements_in_vector*m_b)<<"\n";


        std::chrono::duration<double> elapsed_seconds = end_simple - start_simple;
        total_seconds[i] = elapsed_seconds.count();
        std::chrono::duration<double> elapsed_dgemm = end_dgemm-start_dgemm;
        dgemm_seconds[i] = elapsed_dgemm.count();
    }
    save_result(total_seconds, exp_num, is_automatic);
    save_result(dgemm_seconds, exp_num, is_automatic);
    
    delete[] a; 
    delete[] b;
    delete[] c;
    // return 0;

    a = new double[3*2];
    b = new double [2*4];
    c = new double [3*4];
    generate_test_matrixes(a,b,c);
    simple_matrix_mult(a,b,c,3,4,2);
    for(int i=0;i<3;i++){
        for(int j=0;j<4;j++){
            std::cout<<c[i*4+j]<<" ";
        }
        std::cout<<"\n";
    }

    std::cout<<print_test_result(c)<<"\n";

    for(int i=0;i<3*4;i++){
        c[i]=0.0;
    }

    transpose_matrix(b,2,4);
    matrix_mult_second_transposed(a,b,c,3,4,2);
    std::cout<<"\n";
    for(int i=0;i<3;i++){
        for(int j=0;j<4;j++){
            std::cout<<c[i*4+j]<<" ";
        }
        std::cout<<"\n";
    }
    std::cout<<print_test_result(c)<<"\n";
    return 0;
}

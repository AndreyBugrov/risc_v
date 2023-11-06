#include <iostream>
#include <random>
#include <string> // stoi
#include <chrono> // time

// functor matrix generating
// class matrix
// functor matrix multiplication

const int test_element_num = 12;
const double etalon[test_element_num]={    
    5, -1, 4, -1, 
    -1, 2, 1, 2, 
    3, 0, 3, 0};


void simple_matrix_mult(double* a, double* b, double* c, int n_a, int m_b, int elements_in_vector){
    
    /*
        The i-th row is multiplied by the j-th column of the matrix:
            1. Number of elements (k) in a.row and in b.column are equal: m(a) = n(b)
            2. j = [0, m(b)-1] - column number of b
            3. i = [0, n(a)-1] - row number of a
            4. n(c) = n(a), m(c) = m(b) 
    */

    // for(int i=0;i<n_a;i++){
    //     for(int j=0;j<n_b;j++){
    //         for(int k=0;k<n_b; k++){
    //             std::cout<<"i: "<<i<<", " << "j: "<< j <<", k:"<< k <<"\n";
    //             c[i*n_a+j]+=a[i*n_a+k]*b[j*n_b+k];
    //         }          
    //     }    
    // }

    //    std::cout<<"------------------\n";
    for(int i=0;i<n_a;i++){ // i-th row in a
        for(int j=0;j<m_b;j++){ // j-th column in b
            // std::cout<<"c="<<c[i*m_b+j]<<"\n";
            for(int k=0;k<elements_in_vector; k++){ // k-th element in vector
                // std::cout<<"\n";
                // std::cout<<"a["<<i<<"]["<<k<<"] * b["<<k<<"]["<<j<<"]\n";
                // std::cout<<a[i*elements_in_vector+k]<<" * "<< b[k*m_b+j]<<"\n";
                c[i*m_b+j]+=a[i*elements_in_vector+k]*b[k*m_b+j];
                // c[i*n_a+j]+=a[i*n_a+k]*b[j*m_b+k]; // for b^T or for storing b by columns
            }
            // std::cout<<"------------------";
            // std::cout<<"c["<<i<<"]["<<j<<"]="<<c[i*m_b+j]<<"\n";
            // std::cout<<"------------------"<<"etalon="<<etalon[i*m_b+j]<<"\n";
        }
    }
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
void print_result(const std::chrono::duration<double>& seconds){
    std::cout<<seconds.count()<<"\n";
}

int main(int argc, char* argv[]){
    // n x m matrixes:
    //a: n_a x elements_in_vector
    //b: elements_in_vector x m_b
    const int n_a = argc>1? std::stoi(argv[1]) : 333;
    const int elements_in_vector = argc > 2? std::stoi(argv[2]): 997;
    const int m_b = argc>3? std::stoi(argv[3]): 100;
    double* a, *b, *c;

    a = new double[n_a*elements_in_vector];
    b = new double[elements_in_vector*m_b];
    c = new double[n_a*m_b];

    generate_rand_matrix(a,n_a,elements_in_vector,0.0,10);
    generate_rand_matrix(b,elements_in_vector,m_b,-5.0,5.0);
    generate_zero_matrix(c,n_a,m_b);

    const auto start_simple{std::chrono::steady_clock::now()};
    simple_matrix_mult(a,b,c,n_a,m_b,elements_in_vector);
    const auto end_simple{std::chrono::steady_clock::now()};
    std::chrono::duration<double> elapsed_seconds = end_simple - start_simple;
    print_result(elapsed_seconds);
    
    delete[] a,b,c;

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

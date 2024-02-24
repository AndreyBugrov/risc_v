#include "common.hpp"
const int test_element_num = 12;
const double etalon[test_element_num]={    
    5, -1, 4, -1, 
    -1, 2, 1, 2, 
    3, 0, 3, 0};

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
    if(is_automatic){
        std::cout<<seconds;
    }
    else{
        std::cout<<seconds<<"\n";
    }
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
            return false;
    }
    return true;
}

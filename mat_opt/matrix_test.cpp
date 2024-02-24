// #include "all.hpp"
#include "common.hpp"
#include "matrix_test.hpp"

void print_test_result(std::vector<double> comparing_output, double seconds){
    if(comparing_output[0]>0){
        std::cout<<"FAILED!\n";
        std::cout<<"\tNumber of unequal elements: "<<comparing_output[0]<<"\n";
        std::cout<<"\tMinimum difference: "<<comparing_output[1]<<"\n";
        std::cout<<"\tMaximum difference: "<<comparing_output[2]<<"\n";
    }else{
        std::cout<<"PASSED!\n";
    }
    std::cout<<"\tTime: "<<seconds<<" seconds\n";
}

void fixed_multiplication_test(mult_func matrix_mult_function){
    const auto start_my_mult{std::chrono::steady_clock::now()};

    int n_a = 3;
    int elements_in_vector = 2;
    int m_b = 4;
    double* a, *b, *c;

    a = new double[n_a*elements_in_vector];
    b = new double[elements_in_vector*m_b];
    c = new double[n_a*m_b];

    double multiply_etalon[n_a*m_b]={    
    5, -1, 4, -1, 
    -1, 2, 1, 2, 
    3, 0, 3, 0};

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

    matrix_mult_function(a,b,c,n_a,m_b,elements_in_vector);

    delete[] a;
    delete[] b;
    delete[] c;

    std::vector<double> comparison_output = get_unequal_elements(multiply_etalon, c, n_a*m_b);
    const auto end_my_mult{std::chrono::steady_clock::now()};
    std::chrono::duration<double> elapsed_seconds = end_my_mult - start_my_mult;
    print_test_result(comparison_output, elapsed_seconds.count());
}
void fixed_transpose_matrix_test(tr_func transpose_matrix_function){
    const auto start_my_transpose{std::chrono::steady_clock::now()};

    int n = 103;
    int m = 29;
    double* matrix, *transpose_etalon;
    matrix = new double[n*m];
    transpose_etalon = new double[n*m];

    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            matrix[i*m+j]=i*m+j;
        }
    }
    for(int j=0;j<m;j++){
        for(int i=0;i<n;i++){
            transpose_etalon[i*m+j]=i*m+j;
        }
    }

    transpose_matrix_function(matrix, n, m);

    std::vector<double> comparison_output = get_unequal_elements(transpose_etalon, matrix, n*m);
    const auto end_my_transpose{std::chrono::steady_clock::now()};
    std::chrono::duration<double> elapsed_seconds = end_my_transpose - start_my_transpose;
    print_test_result(comparison_output, elapsed_seconds.count());
}

void open_blas_multiplication_test(mult_func matrix_mult_function){
    const auto start_my_mult{std::chrono::steady_clock::now()};

    int min_random_length = 17;
    int max_random_length = 180;

    std::random_device rd;
    std::mt19937 engine(rd());
    std::uniform_int_distribution<std::mt19937::result_type> gen(min_random_length, max_random_length);

    int n_a = gen(engine);
    int elements_in_vector = gen(engine);
    int m_b = gen(engine);
    double* a, *b, *c, *base_c;

    a = new double[n_a*elements_in_vector];
    b = new double[elements_in_vector*m_b];
    c = new double[n_a*m_b];
    base_c = new double[n_a*m_b];

    generate_rand_matrix(a, n_a, elements_in_vector, -100.0, 100.0);
    generate_rand_matrix(b, elements_in_vector, m_b, -100.0, 100.0);
    generate_zero_matrix(c, n_a, m_b);
    generate_zero_matrix(base_c, n_a, m_b);

    matrix_mult_function(a,b,c,n_a,m_b,elements_in_vector);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, n_a, m_b, elements_in_vector, 1.0, a, elements_in_vector, b, m_b, 0.0, base_c, m_b);

    delete[] a;
    delete[] b;
    delete[] c;

    std::vector<double> comparison_output = get_unequal_elements(base_c, c, n_a*m_b);
    const auto end_my_mult{std::chrono::steady_clock::now()};
    std::chrono::duration<double> elapsed_seconds = end_my_mult - start_my_mult;
    print_test_result(comparison_output, elapsed_seconds.count());
}

void open_blas_transpose_matrix_test(tr_func transpose_matrix_function){
    const auto start_my_transpose{std::chrono::steady_clock::now()};

    int min_random_length = 17;
    int max_random_length = 180;

    std::random_device rd;
    std::mt19937 engine(rd());
    std::uniform_int_distribution<std::mt19937::result_type> gen(min_random_length, max_random_length);

    int n = gen(engine);
    int m = gen(engine);
    double* matrix, *base_matrix;
    matrix = new double[n*m];
    base_matrix = new double[n*m];

    generate_rand_matrix(matrix, n, m, -100.0, 100.0);

    for(int i=0;i<n*m;i++){
        base_matrix[i] = matrix[i];
    }

    transpose_matrix_function(matrix, n, m);

     ///////////////////////////////////////////
    // TO DO: use cblas function and compare //
   ///////////////////////////////////////////

    std::vector<double> comparison_output = get_unequal_elements(base_matrix, matrix, n*m);
    const auto end_my_transpose{std::chrono::steady_clock::now()};
    std::chrono::duration<double> elapsed_seconds = end_my_transpose - start_my_transpose;
    print_test_result(comparison_output, elapsed_seconds.count());
}

std::vector<double> get_unequal_elements(double* base, double* current, int n){
    const double eps = 0.0001;
    double fault_counter = 0;
    double min_difference = DBL_MAX;
    double max_difference = 0.0;
    double difference;
    for(int i=0; i<n; i++){
        difference = abs(base[i]-current[i]);
        if(difference>eps){
            fault_counter++;
            if(min_difference>difference){
                min_difference = difference;
            }
            if(max_difference<difference){
                max_difference = difference;
            }
        }
    }
    double fault_percentage = fault_counter / n * 100.0;
    std::vector<double> out;
    out.push_back(fault_percentage); 
    out.push_back(min_difference);
    out.push_back(max_difference);
    return out;
}

int main(int argc, char* argv[]){
    std::map<std::string, int> all_test_types={{"all", 0}, {"all_mult", 1}, {"all_transpose", 2}, {"special_mult", 3}};
    std::string test_type = std::string(argv[1]);
    int test_type_num = all_test_types[test_type];
    switch (test_type_num)
    {
    case 0:
        break;
    case 1:
        break;
    case 2:
        break;
    case 3:
        break;
    
    default:
        return 0;
    }
}

#include "test.hpp"

bool print_test_result(std::vector<double> comparison_output, double seconds){
    bool passed;
    if(comparison_output[0]>0){
        std::cout<<"FAILED!\n";
        std::cout<<"\tNumber of unequal elements:     "<<comparison_output[0]<<"\n";
        std::cout<<"\tPercentage of unequal elements: "<<comparison_output[1]<<"\n";
        std::cout<<"\tMinimum difference:             "<<comparison_output[2]<<"\n";
        std::cout<<"\tMaximum difference:             "<<comparison_output[3]<<"\n";
        passed = false;
    }else{
        std::cout<<"PASSED!\n";
        passed = true;
    }
    std::cout<<"\tTime: "<<seconds<<" seconds\n";
    std::cout<<"---------------------------------\n";
    return passed;
}
void print_test_statistics(int passed_num, int failed_num, double seconds, std::vector<std::string> failed_test_names){
    int all = passed_num + failed_num;
    double passed_percent = passed_num / double(all) * 100.0;
    double failed_percent = failed_num / double(all) * 100.0;
    std::cout<<"SUMMARY\n";
    std::cout<<"\tTest number: "<<all<<"\n";
    std::cout<<"\tPassed:      "<<passed_num<<" ("<<passed_percent<<"%)\n";
    std::cout<<"\tFailed:      "<<failed_num<<"  ("<<failed_percent<<"%)\n";
    std::cout<<"\tFull time:   "<<seconds<<" seconds\n";
    if(!failed_test_names.empty()){
        std::cout<<"FAILED TESTS:"<<"\n";
        for(auto failed_test_name : failed_test_names){
            std::cout<<"\t"<<failed_test_name<<"\n";
        }
    }
}

bool multiplication_test(std::string test_name, mult_func matrix_mult_function, test_type type){
    std::cout<<"TEST:\t"<<test_name;
    const auto start_my_mult{std::chrono::steady_clock::now()};

    int min_random_length;
    int max_random_length;
    if(type==test_type::big){
        // min_random_length = 128;
        // max_random_length = 256; // right end is included
        min_random_length = kRecursiveStrassenMultLimit<<2;
        max_random_length = kRecursiveStrassenMultLimit<<3; // right end is included
    }else{
        // min_random_length = 16;
        // max_random_length = 64; // right end is included
        min_random_length = kRecursiveStrassenMultLimit;
        max_random_length = kRecursiveStrassenMultLimit<<1; // right end is included
    }

    std::random_device rd;
    std::mt19937 engine(rd());
    std::uniform_int_distribution<std::mt19937::result_type> gen(min_random_length, max_random_length);

    int n = gen(engine);

    std::cout<<" (n = "<<n<<")\n";

    double* a, *b, *c, *base_c;

    a = new double[n*n];
    b = new double[n*n];
    c = new double[n*n];
    base_c = new double[n*n];

    if(type==test_type::identity || type==test_type::zero){
        if(gen(engine)%2){
            if(type==test_type::identity){
                generate_identity_matrix(b, n);
            }else{
                generate_zero_matrix(b, n);
            }
            generate_rand_matrix(a, n, -100.0, 100.0);
        }else{
            if(type==test_type::identity){
                generate_identity_matrix(a, n);
            }else{
                generate_zero_matrix(a, n);
            }
            generate_rand_matrix(b, n, -100.0, 100.0);
        }
    }else{
        generate_rand_matrix(a, n, -100.0, 100.0);
        generate_rand_matrix(b, n, -100.0, 100.0);
    }
    generate_zero_matrix(c, n);
    generate_zero_matrix(base_c, n);

    matrix_mult_function(a, b, c, n);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans, n, n, n, 1.0, a, n, b, n, 0.0, base_c, n);
    std::vector<double> comparison_output = get_unequal_elements(base_c, c, n*n);

    delete[] a;
    delete[] b;
    delete[] c;

    const auto end_my_mult{std::chrono::steady_clock::now()};
    std::chrono::duration<double> elapsed_seconds = end_my_mult - start_my_mult;
    return print_test_result(comparison_output, elapsed_seconds.count());
}
bool split_matrices_test(test_type type){
    std::string test_name = "split matrices";
    std::cout<<"TEST:\t"<<test_name;
    const auto start_my_mult{std::chrono::steady_clock::now()};

    int min_random_deg;
    int max_random_deg;
    if(type==test_type::big){
        min_random_deg = 10;
        max_random_deg = 12; // right end is included
    }else{
        min_random_deg = 6;
        max_random_deg = 9; // right end is included
    }

    std::random_device rd;
    std::mt19937 engine(rd());
    std::uniform_int_distribution<std::mt19937::result_type> gen(min_random_deg, max_random_deg);

    double* a, *a11, *a12, *a21, *a22, *zero_a_ij, *identity_a_ij;
    int n = 1 << gen(engine);
    a = new double[n*n];
    generate_identity_matrix(a, n);
    std::cout<<" (n = "<<n<<")\n";
    n >>= 1;

    a11 = new double[n*n];
    a12 = new double[n*n];
    a21 = new double[n*n];
    a22 = new double[n*n];
    zero_a_ij = new double[n*n];
    identity_a_ij = new double[n*n];

    generate_identity_matrix(identity_a_ij, n);
    generate_zero_matrix(zero_a_ij, n);


    split_matrices(a, a11, a12, a21, a22, 2*n);
    std::vector<double> comparison_output = get_unequal_elements(identity_a_ij, a11, n*n);

    const auto end_my_mult{std::chrono::steady_clock::now()};
    std::chrono::duration<double> elapsed_seconds = end_my_mult - start_my_mult;
    return print_test_result(comparison_output, elapsed_seconds.count());
}
bool collect_matrices_test(test_type type){
    std::string test_name = "collect matrices";
    std::cout<<"TEST:\t"<<test_name;
    const auto start_my_mult{std::chrono::steady_clock::now()};

    int min_random_deg;
    int max_random_deg;
    if(type==test_type::big){
        min_random_deg = 10;
        max_random_deg = 12; // right end is included
    }else{
        min_random_deg = 6;
        max_random_deg = 9; // right end is included
    }

    std::random_device rd;
    std::mt19937 engine(rd());
    std::uniform_int_distribution<std::mt19937::result_type> gen(min_random_deg, max_random_deg);

    double* a, *a11, *a12, *a21, *a22, *zero_a_ij, *identity_a_ij;
    int n = 1 << gen(engine);
    a = new double[n*n];
    std::cout<<" (n = "<<n<<")\n";
    n >>= 1;

    a11 = new double[n*n];
    a12 = new double[n*n];
    a21 = new double[n*n];
    a22 = new double[n*n];
    zero_a_ij = new double[n*n];
    identity_a_ij = new double[n*n];

    generate_identity_matrix(identity_a_ij, n);
    generate_zero_matrix(zero_a_ij, n);
    generate_identity_matrix(a11, n);
    generate_identity_matrix(a22, n);
    generate_zero_matrix(a12, n);
    generate_zero_matrix(a21, n);


    collect_matrices(a, a11, a12, a21, a22, 2*n);
    std::vector<double> comparison_output = get_unequal_elements(identity_a_ij, a11, n*n);

    const auto end_my_mult{std::chrono::steady_clock::now()};
    std::chrono::duration<double> elapsed_seconds = end_my_mult - start_my_mult;
    return print_test_result(comparison_output, elapsed_seconds.count());
}
bool split_and_collect_matrices_test(std::string test_name, test_type type){
    std::cout<<"TEST:\t"<<test_name;
    const auto start_my_mult{std::chrono::steady_clock::now()};

    int min_random_deg;
    int max_random_deg;
    if(type==test_type::big){
        min_random_deg = 10;
        max_random_deg = 12; // right end is included
    }else{
        min_random_deg = 6;
        max_random_deg = 9; // right end is included
    }

    std::random_device rd;
    std::mt19937 engine(rd());
    std::uniform_int_distribution<std::mt19937::result_type> gen(min_random_deg, max_random_deg);

    double* base_a, *a11, *a12, *a21, *a22, *a;
    int n = 1 << gen(engine);
    base_a = new double[n*n];
    a = new double[n*n];
    std::cout<<" (n = "<<n<<")\n";
    
    switch (type)
    {
        case test_type::zero:
            generate_zero_matrix(base_a, n);
            break;
        case test_type::identity:
            generate_identity_matrix(base_a, n);
            break;
        default:
            generate_rand_matrix(base_a, n, -100, 100);
            break;
    }
    
    n >>= 1;
    a11 = new double[n*n];
    a12 = new double[n*n];
    a21 = new double[n*n];
    a22 = new double[n*n];
    generate_zero_matrix(a11, n);
    generate_zero_matrix(a12, n);
    generate_zero_matrix(a21, n);
    generate_zero_matrix(a22, n);
    n <<=1;

    split_matrices(base_a, a11, a12, a21, a22, n);
    collect_matrices(a, a11, a12, a21, a22, n);
    std::vector<double> comparison_output = get_unequal_elements(base_a, a, n*n);

    const auto end_my_mult{std::chrono::steady_clock::now()};
    std::chrono::duration<double> elapsed_seconds = end_my_mult - start_my_mult;
    return print_test_result(comparison_output, elapsed_seconds.count());
}
bool increase_and_decrease_matrices(std::string test_name){
    std::cout<<"TEST:\t"<<test_name;
    const auto start_my_mult{std::chrono::steady_clock::now()};

    const int min_random_length = 512;
    const int max_random_length = 2048; // right end is included

    std::random_device rd;
    std::mt19937 engine(rd());
    std::uniform_int_distribution<std::mt19937::result_type> gen(min_random_length, max_random_length);

    double* base_a, *a;
    int n = gen(engine);
    base_a = new double[n*n];
    a = new double[n*n];
    std::cout<<" (n = "<<n<<")\n";

    int inc_dim = 1;
    int offset = log2(n);
    if(offset<log2(n)){
        offset++;
    }
    inc_dim <<= offset;

    generate_rand_matrix(base_a, n, -100, 100);

    double* inc_a = new double[inc_dim*inc_dim];

    for(int i=0;i<n;i++){
        std::memcpy(inc_a+i*inc_dim, base_a+i*n, n*sizeof(double));
    }
    decrease_matrix(a, inc_a, n, inc_dim);
    std::vector<double> comparison_output = get_unequal_elements(base_a, a, n*n);

    const auto end_my_mult{std::chrono::steady_clock::now()};
    std::chrono::duration<double> elapsed_seconds = end_my_mult - start_my_mult;
    return print_test_result(comparison_output, elapsed_seconds.count());
}
bool matrix_alg_sum_test(std::string test_name, test_type type, bool is_add){
    std::cout<<"TEST:\t"<<test_name;
    const auto start_my_mult{std::chrono::steady_clock::now()};

    int min_random_length=512;
    int max_random_length=2048;

    std::random_device rd;
    std::mt19937 engine(rd());
    std::uniform_int_distribution<std::mt19937::result_type> gen(min_random_length, max_random_length);

    int n = gen(engine);

    std::cout<<" (n = "<<n<<")\n";

    double* a, *b, *c, *base_c;

    a = new double[n*n];
    b = new double[n*n];
    c = new double[n*n];
    base_c = new double[n*n];

    if(type==test_type::identity || type==test_type::zero){
        if(type==test_type::identity){
            generate_identity_matrix(b, n);
            generate_identity_matrix(a, n);
            generate_zero_matrix(base_c, n);
            if(is_add){
                for(int i=0;i<n;i++){
                    base_c[i*n+i]=2.0;
                }
            }else{
                generate_zero_matrix(base_c, n);
            }

        }
        else{
            generate_zero_matrix(a, n);
            generate_zero_matrix(b, n);
            generate_zero_matrix(base_c, n);
        }
    }else{
        return true;
    }
    if(is_add){
        matrix_add(a, b, c, n);
    }else{
        matrix_sub(a, b, c, n);
    }

    std::vector<double> comparison_output = get_unequal_elements(base_c, c, n*n);

    delete[] a;
    delete[] b;
    delete[] c;

    const auto end_my_mult{std::chrono::steady_clock::now()};
    std::chrono::duration<double> elapsed_seconds = end_my_mult - start_my_mult;
    return print_test_result(comparison_output, elapsed_seconds.count());
}

std::vector<double> get_unequal_elements(double* base, double* current, int size){
    const double eps = 0.0001;
    double fault_counter = 0;
    double min_difference = DBL_MAX;
    double max_difference = 0.0;
    double difference;
    for(int i=0; i<size; i++){
        difference = abs(base[i]-current[i]);
        if(difference>eps){
            // std::cout<<"size="<<size<<"\n";
            // std::cout<<i<<"\n";
            fault_counter++;
            if(min_difference>difference){
                min_difference = difference;
            }
            if(max_difference<difference){
                max_difference = difference;
            }
        }
    }
    double fault_percentage = fault_counter / size * 100.0;
    std::vector<double> out;
    out.push_back(fault_counter);
    out.push_back(fault_percentage); 
    out.push_back(min_difference);
    out.push_back(max_difference);
    return out;
}

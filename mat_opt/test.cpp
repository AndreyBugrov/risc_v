#include "test.hpp"

bool print_test_result(std::vector<double> comparing_output, double seconds){
    bool passed;
    if(comparing_output[0]>0){
        std::cout<<"FAILED!\n";
        std::cout<<"\tPercentage of unequal elements: "<<comparing_output[0]<<"\n";
        std::cout<<"\tMinimum difference: "<<comparing_output[1]<<"\n";
        std::cout<<"\tMaximum difference: "<<comparing_output[2]<<"\n";
        passed = false;
    }else{
        std::cout<<"PASSED!\n";
        passed = true;
    }
    std::cout<<"\tTime: "<<seconds<<" seconds\n";
    std::cout<<"---------------------------------\n";
    return passed;
}
void print_test_statistics(int passed_num, int failed_num, double seconds){
    int all = passed_num + failed_num;
    double pased_percent = passed_num / double(all) * 100.0;
    double failed_percent = failed_num / double(all) * 100.0;
    std::cout<<"SUMMARY\n";
    std::cout<<"\tTest number: "<<all<<"\n";
    std::cout<<"\tPassed:      "<<passed_num<<" ("<<pased_percent<<"%)\n";
    std::cout<<"\tFailed:      "<<failed_num<<"  ("<<failed_percent<<"%)\n";
    std::cout<<"\tFull time:   "<<seconds<<"\n";
}

bool multiplication_test(std::string test_name, mult_func matrix_mult_function, test_type type){
    std::cout<<"TEST:\t"<<test_name<<"\n";
    const auto start_my_mult{std::chrono::steady_clock::now()};

    int min_random_length;
    int max_random_length;
    if(type==test_type::big){
        min_random_length = 128;
        max_random_length = 256+1;
    }else{
        min_random_length = 16;
        max_random_length = 64+1;
    }

    std::random_device rd;
    std::mt19937 engine(rd());
    std::uniform_int_distribution<std::mt19937::result_type> gen(min_random_length, max_random_length);

    int n = gen(engine);

    double* a, *b, *c, *base_c;

    a = new double[n*n];
    b = new double[n*n];
    c = new double[n*n];
    base_c = new double[n*n];

    if(type==test_type::identity || type==test_type::zero){
        if(gen(engine)%2){
            if(type==test_type::identity){
                generate_zero_matrix(b, n);
            }else{
                generate_identity_matrix(b, n);
            }
            generate_rand_matrix(a, n, -100.0, 100.0);
        }else{
            if(type==test_type::identity){
                generate_zero_matrix(a, n);
            }else{
                generate_identity_matrix(a, n);
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
    std::vector<double> comparison_output = get_unequal_elements(base_c, c, n);

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
    // std::cout<<fault_counter<<" "<<size<<"\n";
    std::vector<double> out;
    out.push_back(fault_percentage); 
    out.push_back(min_difference);
    out.push_back(max_difference);
    return out;
}

#include "test.hpp"

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

void multiplication_test(std::string test_name, mult_func matrix_mult_function, test_type type){
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

    delete[] a;
    delete[] b;
    delete[] c;

    std::vector<double> comparison_output = get_unequal_elements(base_c, c, n);
    const auto end_my_mult{std::chrono::steady_clock::now()};
    std::chrono::duration<double> elapsed_seconds = end_my_mult - start_my_mult;
    print_test_result(comparison_output, elapsed_seconds.count());
}

void transpose_test(std::string test_name, test_type type){
    std::cout<<"TEST:\t"<<test_name<<"\n";
    const auto start_my_transpose{std::chrono::steady_clock::now()};

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
    double* matrix, *base_matrix;
    matrix = new double[n*n];
    base_matrix = new double[n*n];

    switch (type)
    {
    case test_type::zero:
        generate_zero_matrix(matrix, n);
        break;
    case test_type::identity:
        generate_identity_matrix(matrix, n);
    default:
        generate_rand_matrix(matrix, n, -100.0, 100.0);
        break;
    }

    std::memcpy(base_matrix, matrix, n*n*sizeof(double));

    matrix = get_transposed_matrix(matrix, n);

     ///////////////////////////////////////////
    // TO DO: use cblas function and compare //
   ///////////////////////////////////////////

    std::vector<double> comparison_output = get_unequal_elements(base_matrix, matrix, n*n);
    const auto end_my_transpose{std::chrono::steady_clock::now()};
    std::chrono::duration<double> elapsed_seconds = end_my_transpose - start_my_transpose;
    print_test_result(comparison_output, elapsed_seconds.count());
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
    std::vector<double> out;
    out.push_back(fault_percentage); 
    out.push_back(min_difference);
    out.push_back(max_difference);
    return out;
}

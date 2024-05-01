#include "test.hpp"

int main(int argc, char* argv[]){
    const auto start_tests{std::chrono::steady_clock::now()};
    std::vector<std::pair<std::string, mult_func>> test_mult_funcs={
        {"base multiplication", base_matrix_mult}, 
        {"base multiplication (omp)", base_matrix_mult_omp},
        {"row multiplication", row_matrix_mult},
        {"row multiplication (simd)", row_matrix_mult_simd},
        {"row multiplication (optimized)", row_matrix_mult_opt},
        {"row muliplication (optimized+simd)", row_matrix_mult_opt_simd},
        {"row multiplication (omp)", row_matrix_mult_omp},
        {"row multiplication (optimized+omp)", row_matrix_mult_opt_omp},
        {"row multiplication (omp+simd)", row_matrix_mult_omp_simd},
        {"row multiplication (optimized+omp+simd)", row_matrix_mult_opt_omp_simd},
        {"multiplication with transposed matrix", transposed_matrix_mult},
        {"multiplication with transposed matrix (omp)", transposed_matrix_mult_omp},
        {"multiplication with transposed matrix (omp+simd)", transposed_matrix_mult_omp_simd},
        {"strassen multiplication", strassen_matrix_mult},
        {"strassen multiplication (naive omp)", strassen_matrix_mult_omp},
        {"strassen multiplication (recurcive omp)", strassen_matrix_mult_rec_omp}};
    std::map<int, std::string> test_type_names={{0, "zero"}, {1, "identity"}, {2, "equal"}, {3, "random"}, {4, "big"}};
    int passed_counter = 0;
    int failed_counter = 0;
    std::vector<std::string> failed_test_names;
    std::string full_test_name;

    std::string filename = "cache.txt";
    int* cache_sizes = new int[kCacheLevelNumber];
    get_cache_sizes(filename, cache_sizes);
    // for(int i=0;i<kCacheLevelNumber;i++){
    //     std::cout<<i<<": "<<cache_sizes[i]<<"\n";
    // }

    if (split_matrices_test(test_type::big)){
        passed_counter++;
    }else{
        failed_counter++;
        full_test_name = "collect matrices";
        failed_test_names.push_back(full_test_name);
    }
    if (collect_matrices_test(test_type::big)){
        passed_counter++;
    }else{
        failed_counter++;
        full_test_name = "collect matrices";
        failed_test_names.push_back(full_test_name);
    }
    std::string short_test_name = "split + collect matrices";

    for(int i=0; i<=static_cast<int>(test_type::big); i++){
        if(i==static_cast<int>(test_type::equal)){
            continue;
        }
        full_test_name = short_test_name + std::string(" : "+test_type_names[i]);
        if(split_and_collect_matrices_test(full_test_name, static_cast<test_type>(i))){
            passed_counter++;
        }else{
            failed_counter++;
            failed_test_names.push_back(full_test_name);
        }
    }
    short_test_name = "increase and decrease matrices";
    if(increase_and_decrease_matrices(short_test_name)){
        passed_counter++;
    }else{
        failed_counter++;
        failed_test_names.push_back(full_test_name);
    }

    short_test_name = "matrix addition test";

    for(int i=0; i<=static_cast<int>(test_type::identity); i++){
        full_test_name = short_test_name + std::string(" : "+test_type_names[i]);
        if(matrix_alg_sum_test(full_test_name, static_cast<test_type>(i), true)){
            passed_counter++;
        }else{
            failed_counter++;
            failed_test_names.push_back(full_test_name);
        }
    }
    short_test_name = "matrix subtraction test";
    for(int i=0; i<=static_cast<int>(test_type::identity); i++){
        full_test_name = short_test_name + std::string(" : "+test_type_names[i]);
        if(matrix_alg_sum_test(full_test_name, static_cast<test_type>(i), false)){
            passed_counter++;
        }else{
            failed_counter++;
            failed_test_names.push_back(full_test_name);
        }
        
    }
    for(auto test_arguments : test_mult_funcs){
        for(int i=0; i<=static_cast<int>(test_type::big); i++){
            full_test_name = test_arguments.first + std::string(" : "+test_type_names[i]);
            if(multiplication_test(full_test_name, test_arguments.second, static_cast<test_type>(i))){
                passed_counter++;
            }else{
                failed_counter++;
                failed_test_names.push_back(full_test_name);
            }
        }
    }
    const auto end_tests{std::chrono::steady_clock::now()};
    std::chrono::duration<double> elapsed_seconds = end_tests - start_tests;
    print_test_statistics(passed_counter, failed_counter, elapsed_seconds.count(), failed_test_names);
    std::cout<<sizeof(double)<<"\n";
}

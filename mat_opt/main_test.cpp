#include "test.hpp"

int main(int argc, char* argv[]){
    std::vector<std::pair<std::string, mult_func>> test_mult_funcs={
        {"base multiplication", base_matrix_mult}, 
        {"base multiplication (omp)", base_matrix_mult_omp},
        {"multiplication with transposed matrix", transposed_matrix_mult},
        {"multiplication with transposed matrix (omp)", transposed_matrix_mult_omp}};
    std::map<int, std::string> test_type_names={{0, "zero"}, {1, "identity"}, {2, "equal"}, {3, "random"}, {4, "big"}};
    int passed_counter = 0;
    int failed_counter = 0;
    for(auto test_arguments : test_mult_funcs){
        for(int i=0; i<=static_cast<int>(test_type::big); i++){
            if(multiplication_test(test_arguments.first + std::string(" : "+test_type_names[i]), test_arguments.second, static_cast<test_type>(i))){
                passed_counter++;
            }else{
                failed_counter++;
            }
        }
    }
    print_test_statistics(passed_counter, failed_counter);
}
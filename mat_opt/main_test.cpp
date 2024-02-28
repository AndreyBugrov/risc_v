#include "test.hpp"

int main(int argc, char* argv[]){
    std::map<std::string, int> all_test_types={{"all", 0}, {"all_mult", 1}, {"all_transpose", 2}};
    std::string test_type = argc>1? std::string(argv[1]): "all";
    const auto found = all_test_types.find(test_type);
    int test_type_num = 999;
    if(found == all_test_types.cend()){
        test_type_num = 3;
    }else{
        test_type_num = all_test_types[test_type];
    }
    switch (test_type_num)
    {
    case 0:
    case 1:{
        std::vector<std::pair<std::string, mult_func>> test_mult_funcs={{"base multiplication", base_matrix_mult}, 
        {"base multiplication (omp included)", base_matrix_mult_omp},
        {"multiplication with second transposed matrix", b_transposed_matrix_mult},
        {"multiplication with second transposed matrix (omp included)", b_transposed_matrix_mult_omp}};
    }
    if(!test_type_num){
            break;
    }
    case 2:{

    }
        break;
    case 3:
        // specified test type parcing
        break;
    
    default:
        std::cout<<"Error: incorrect type name\n";
        return 0;
    }
}

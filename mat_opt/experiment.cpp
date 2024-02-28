#include "experiment.hpp"

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

void print_result(double* total_seconds, int exp_num){
  if(exp_num==1){
    std::cout<<total_seconds[0]<<"\n";
  }
  else{
    int max_index = get_max_value_index(total_seconds, exp_num);
    double seconds_without_outliers = 0;
    for(int i=0;i<exp_num;i++){
        if(i!=max_index){
            seconds_without_outliers+=total_seconds[i];
        }
    }
    std::cout<<seconds_without_outliers / (exp_num-1)<<"\n";
  }
}

mult_func get_multiplication_function(std::string function_name){
    const int names_num = 5;
    std::string all_function_names[names_num]={"simple", "transposed", "opt_block", "block_transposed"};
    const std::map<std::string, mult_func> func_map={{all_function_names[0], base_matrix_mult},
    {all_function_names[1], b_transposed_matrix_mult},
    {all_function_names[2], optimal_block_matrix_mult},
    {all_function_names[3], b_transposed_block_matrix_mult}
    };
    try{
        return func_map.at(function_name);
    }
    catch(std::out_of_range ex){
        std::string msg = "Unknown function name: \"" + function_name + "\"\nUse theese: \"" + all_function_names[0]+"\"";
        for(int i=1;i<names_num;i++){
            msg+=", \""+all_function_names[i]+"\"";
        }
        msg += "\n";
        throw msg;
    }
}

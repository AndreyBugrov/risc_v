#include "experiment.hpp"

void print_experiment_result(double* cblas_seconds, double* current_seconds, int experiment_num, double* base, double* current, int size, bool is_automatic) {
  double cblas_seconds_mean = std::accumulate(cblas_seconds, cblas_seconds+experiment_num, 0.0, std::plus<double>()) / experiment_num;
  double current_seconds_mean = std::accumulate(current_seconds, current_seconds+experiment_num, 0.0, std::plus<double>()) / experiment_num;
  double max_difference = 0.0;
  double difference;
  for(int i=0;i<size;i++){
    difference = abs(base[i]-current[i]);
    if(difference>max_difference){
        max_difference = difference;
    }
  }
  if(is_automatic){
      std::cout << cblas_seconds_mean << "\n" << current_seconds_mean << "\n" << current_seconds_mean / cblas_seconds_mean << "\n" << max_difference << "\n";
  }else{
    std::cout<<"OpenBLAS time:  " << cblas_seconds_mean<<"\n";
    std::cout<<"Current time:   " << current_seconds_mean<<"\n";
    std::cout<<"Ratio:          " << current_seconds_mean / cblas_seconds_mean<<"\n";
    std::cout<<"Max inaccuracy: " << max_difference <<"\n";
  }
}

mult_func get_multiplication_function(std::string function_name){
    const int names_num = 16;
    std::string all_function_names[names_num]={"base", "base_omp", "row", "row_omp", "row_omp_simd", "tr", "tr_omp", "tr_omp_simd", "strassen", 
    "strassen_omp", "strassen_rec_omp", "row_opt", "row_opt_simd","row_opt_omp", "row_opt_omp_simd", "row_simd"};
    const std::map<std::string, mult_func> func_map={
    {all_function_names[0], base_matrix_mult},
    {all_function_names[1], base_matrix_mult_omp},
    {all_function_names[2], row_matrix_mult},
    {all_function_names[3], row_matrix_mult_omp},
    {all_function_names[4], row_matrix_mult_omp_simd},
    {all_function_names[5], transposed_matrix_mult},
    {all_function_names[6], transposed_matrix_mult_omp},
    {all_function_names[7], transposed_matrix_mult_omp_simd}, // do not use it! it is slow function
    {all_function_names[8], strassen_matrix_mult},
    {all_function_names[9], strassen_matrix_mult_omp},
    {all_function_names[10], strassen_matrix_mult_rec_omp},
    {all_function_names[11], row_matrix_mult_opt},
    {all_function_names[12], row_matrix_mult_opt_simd},
    {all_function_names[13], row_matrix_mult_opt_omp},
    {all_function_names[14], row_matrix_mult_opt_omp_simd},
    {all_function_names[15], row_matrix_mult_simd}
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
        std::cout<<msg;
        throw msg;
    }
}

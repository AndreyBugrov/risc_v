#define _USE_MATH_DEFINES  // for C++
#include <iostream> // std::cin and std::cout
#include <iomanip> // std::fixed and std::setprecision
#include <cmath> // M_PI
#include <string> // std::stoi
#include <chrono> // time
#include <omp.h> //Open MP
double pi_rectangle(double x) {
  return (1.0 / (1.0 + x * x)); // only height
}

void print_sum_and_pi(double sum, double seconds, bool is_automatic) {
  if (is_automatic){
    std::cout<<fabs(M_PI - sum)<<" "<<seconds;
  }else{
    std::cout <<"Result:     "<< sum << "\n";
    std::cout <<"Inaccuracy: "<< fabs(M_PI - sum) << "\n";
    std::cout <<"Seconds:    " << seconds<<"\n";
    std::cout<<"\n";
  }
}

void zero_vector(double* vec, int n){
  for(int i=0;i<n;i++){
    vec[i]=0.0;
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

void save_result(double* total_seconds, int exp_num, double& sum, bool is_automatic){
  if(exp_num==1){
    print_sum_and_pi(sum, total_seconds[0], is_automatic);
  }
  else{
    int max_index = get_max_value_index(total_seconds, exp_num);
    double seconds_without_outliers = 0;
    for(int i=0;i<exp_num;i++){
      if(i!=max_index){
        seconds_without_outliers+=total_seconds[i];
      }
    }
    print_sum_and_pi(sum, seconds_without_outliers / (exp_num-1), is_automatic);
  }
}

void merge_sum(double* thread_sum, double& sum, double step, int proc){
  for(int i=0;i<proc;i++){
    sum+=thread_sum[i];
  }            
  sum *= 4.0 * step;
}

enum class counting_type{
  left,
  middle,
  right,
  all
};

int main(int argc, char* argv[]) {
  const bool is_automatic = std::string(argv[1]) == std::string("a") ? true : false; // it is strange but argv[1] != "a" && argv[1] != "a\0"
  const int N = argc > 1 ? std::stoi(argv[2]) : 1000;
  counting_type type;
  if(argc>3){
    int mode = std::stoi(std::string(argv[3]));
    switch (mode)
    {
    case 1:
      type = counting_type::left;
      break;
    case 2:
      type = counting_type::middle;
      break;
    case 3:
      type = counting_type::right;
      break;
    case 0:
      type = counting_type::all;
      break;
    default:
      break;
    }
  }
  const int exp_num = (argc>4)? std::stoi(argv[4]) : 1;

  double step = 1.0 / N;
  double sum = 0.0;
  const int precision = 13;
  std::cout << std::fixed; // to set precision to every value
  std::cout << std::setprecision(precision);
  double* total_seconds = new double[exp_num];
  double average_seconds;

  int proc = omp_get_num_procs();
  double* thread_sum = new double[proc];

  zero_vector(total_seconds, exp_num);
  if(!is_automatic||is_automatic&&(type==counting_type::all||type==counting_type::left)){
    for(int n=0; n<exp_num; n++){
      sum = 0.0; /// zeroing every exp_num not to bring garbage to the next experiment.
                 /// do not do it at the end of iteration to prevent saving null sum
      const auto start_left{std::chrono::steady_clock::now()};
      zero_vector(thread_sum, proc);
      #pragma omp parallel num_threads(proc) default(none) shared(step, N, thread_sum)
      {
        int thread_num = omp_get_thread_num();
        #pragma omp for 
        for (int i = 0; i < N; i++) {
          thread_sum[thread_num] += pi_rectangle(i * step);
        }
      }
      merge_sum(thread_sum, sum, step, proc);
      const auto end_left{std::chrono::steady_clock::now()};
      std::chrono::duration<double> elapsed_seconds = end_left - start_left;
      total_seconds[n] = elapsed_seconds.count();
    }
    save_result(total_seconds, exp_num, sum, is_automatic);
  }

  if(!is_automatic||is_automatic&&(type==counting_type::all||type==counting_type::middle)){
  for(int n=0; n<exp_num; n++){
    sum = 0.0; /// zeroing every exp_num not to bring garbage to the next experiment.
                /// do not do it at the end of iteration to prevent saving null sum
    zero_vector(thread_sum, proc);
    const auto start_middle{std::chrono::steady_clock::now()};
    #pragma omp parallel num_threads(proc) default(none) shared(step, N, thread_sum)
    {
      int thread_num = omp_get_thread_num();
      #pragma omp for
      for (int i=0; i < N; i++) {
        thread_sum[thread_num] += pi_rectangle((i * step + (i + 1) * step) * 0.5);
      }
    }
    merge_sum(thread_sum, sum, step, proc);
    const auto end_middle{std::chrono::steady_clock::now()};
    std::chrono::duration<double> elapsed_seconds = end_middle - start_middle;
    total_seconds[n] = elapsed_seconds.count();
  }
  save_result(total_seconds, exp_num, sum, is_automatic);
  }

  if(!is_automatic||is_automatic&&(type==counting_type::all||type==counting_type::right)){
    for(int n=0; n<exp_num; n++){
      sum = 0.0; /// zeroing every exp_num not to bring garbage to the next experiment.
                 /// do not do it at the end of iteration to prevent saving null sum
      const auto start_right{std::chrono::steady_clock::now()};
      zero_vector(thread_sum, proc);
      #pragma omp parallel num_threads(proc) default(none) shared(step, N, thread_sum)
      {
        int thread_num = omp_get_thread_num();
        #pragma omp for
        for (int i = 0; i < N; i++) {
          thread_sum[thread_num] += pi_rectangle((i + 1) * step);
        }
      }
      merge_sum(thread_sum, sum, step, proc);
      const auto end_right{std::chrono::steady_clock::now()};
      std::chrono::duration<double> elapsed_seconds = end_right - start_right;
      total_seconds[n] = elapsed_seconds.count();
    }
    save_result(total_seconds, exp_num, sum, is_automatic);
  }

  delete[] total_seconds;
  return 0;
}

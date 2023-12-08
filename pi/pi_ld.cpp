#define _USE_MATH_DEFINES  // for C++
#include <iostream> // std::cin and std::cout
#include <iomanip> // std::fixed and std::setprecision
#include <cmath> // M_PI
#include <string> // std::stoi
#include <chrono> // time
long double pi_rectangle(long double x) {
  return (1.0 / (1.0 + x * x)); // only height
}

void print_sum_and_pi(double sum, double seconds, bool is_automatic) {
  if (is_automatic){
    std::cout<<fabsl(M_PI - sum)<<" "<<seconds;
  }else{
    std::cout <<"Result:     "<< sum << "\n";
    std::cout <<"Inaccuracy: "<< fabs(M_PI - sum) << "\n";
    std::cout <<"Seconds:    " << seconds<<"\n";
    std::cout<<"\n";
  }
}

// zero seconds vector therefore type of vec is double*
void zero_vector_d(double* vec, int n){
  for(int i=0;i<n;i++){
    vec[i]=static_cast<long double>(0.0);
  }
}

// compare seconds therefore type of vec is double*
int get_max_value_index_d(double* vec, int n){
    long double max_value = vec[0];
    int max_index = 0;
    for(int i=1;i<n;i++){
        if(max_value<vec[i]){
        max_value = vec[i];
        max_index = i;
        }
    }
    return max_index;
}

void save_result(double* total_seconds, int exp_num, long double& sum, bool is_automatic){
  if(exp_num==1){
    print_sum_and_pi(sum, total_seconds[0], is_automatic);
  }
  else{
      int max_index = get_max_value_index_d(total_seconds, exp_num);
    double seconds_without_outliers;
    for(int i=0;i<exp_num;i++){
      if(i!=max_index){
        seconds_without_outliers+=total_seconds[i];
      }
    }
    print_sum_and_pi(sum, seconds_without_outliers / (exp_num-1), is_automatic);
  }
  sum = static_cast<long double> (0.0);
  zero_vector_d(total_seconds, exp_num);
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

  long double step = static_cast<long double>(1.0 / N);
  long double sum = static_cast<long double> (0.0);
  const int precision = 13;
  std::cout << std::fixed; // to set precision to every value
  std::cout << std::setprecision(precision);
  double* total_seconds = new double[exp_num];
  double average_seconds;

  zero_vector_d(total_seconds, exp_num);
  if(!is_automatic||is_automatic&&(type==counting_type::all||type==counting_type::left)){
    for(int n=0; n<exp_num; n++){
      const auto start_left{std::chrono::steady_clock::now()};
      for (int i = 0; i < N; i++) {
        sum += pi_rectangle(i * step);
      }
      sum *= 4.0 * step;
      const auto end_left{std::chrono::steady_clock::now()};
      std::chrono::duration<double> elapsed_seconds = end_left - start_left;
      total_seconds[n] = elapsed_seconds.count();
    }
    save_result(total_seconds, exp_num, sum, is_automatic);
  }

  if(!is_automatic||is_automatic&&(type==counting_type::all||type==counting_type::middle)){
    for(int n=0; n<exp_num; n++){
      const auto start_middle{std::chrono::steady_clock::now()};
      for (int i = 0; i < N; i++) {
        sum += pi_rectangle((i * step + (i + 1) * step) * 0.5);
      }
      sum *= 4.0 * step;
      const auto end_middle{std::chrono::steady_clock::now()};
      std::chrono::duration<double> elapsed_seconds = end_middle - start_middle;
      total_seconds[n] = elapsed_seconds.count();
    }
    save_result(total_seconds, exp_num, sum, is_automatic);
  }

  if(!is_automatic||is_automatic&&(type==counting_type::all||type==counting_type::right)){
    for(int n=0; n<exp_num; n++){
      const auto start_right{std::chrono::steady_clock::now()};
      for (int i = 0; i < N; i++) {
        sum += pi_rectangle((i + 1) * step);
      }
      sum *= 4.0 * step;
      const auto end_right{std::chrono::steady_clock::now()};
      std::chrono::duration<double> elapsed_seconds = end_right - start_right;
      total_seconds[n] = elapsed_seconds.count();
    }
    save_result(total_seconds, exp_num, sum, is_automatic);
  }

  delete[] total_seconds;
  return 0;
}

#define _USE_MATH_DEFINES  // for C++
#include <iostream> // std::cin and std::cout
#include <iomanip> // std::fixed and std::setprecision
#include <cmath> // M_PI
#include <string> // std::stoi
#include <chrono> // time
double pi_rectangle(double x) {
  return (1.0 / (1.0 + x * x)); // only height
}

void print_sum_and_pi(double sum, const std::chrono::duration<double>& seconds, bool is_automatic) {
  if (is_automatic){
    std::cout<<seconds.count();
  }else{
    std::cout <<"Result:     "<< sum << "\n";
    std::cout <<"Inaccuracy: "<< fabs(M_PI - sum) << "\n";
    std::cout <<"Seconds:    " << seconds.count()<<"\n";
    std::cout<<"\n";
  }
}

/// 1 000 000 000 operations make inaccuracy 0 in every test
/// 100 000 operations make inaccuracy 0 in middle rectangles

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
      /* code */
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

  double step = static_cast<double>(1.0) / N;
  double sum = 0.0;
  const int precision = 10;
  std::cout << std::fixed; // to set precision to every value
  std::cout << std::setprecision(precision);
  std::chrono::duration<double> elapsed_seconds;

  if(!is_automatic||is_automatic&&(type==counting_type::all||type==counting_type::left)){
    const auto start_left{std::chrono::steady_clock::now()};
    for (int i = 0; i < N; i++) {
      sum += pi_rectangle(static_cast<double>(i) * step);
    }
    sum *= 4.0 * step;
    const auto end_left{std::chrono::steady_clock::now()};
    elapsed_seconds = end_left - start_left;
    print_sum_and_pi(sum, elapsed_seconds, is_automatic);
  }

  sum = 0.0;
    if(!is_automatic||is_automatic&&(type==counting_type::all||type==counting_type::middle)){
    const auto start_middle{std::chrono::steady_clock::now()};
    for (int i = 0; i < N; i++) {
      sum += pi_rectangle((i * step + (i + 1) * step) * 0.5);
    }
    sum *= 4.0 * step;
    const auto end_middle{std::chrono::steady_clock::now()};
    elapsed_seconds = end_middle - start_middle;
    print_sum_and_pi(sum, elapsed_seconds, is_automatic);
  }

  sum = 0.0;
  if(!is_automatic||is_automatic&&(type==counting_type::all||type==counting_type::right)){
    const auto start_right{std::chrono::steady_clock::now()};
    for (int i = 0; i < N; i++) {
      sum += pi_rectangle((static_cast<double>(i) + 1) * step);
    }
    sum *= 4.0 * step;
    const auto end_right{std::chrono::steady_clock::now()};
    elapsed_seconds = end_right - start_right;
    print_sum_and_pi(sum, elapsed_seconds, is_automatic);
  }

  return 0;
}

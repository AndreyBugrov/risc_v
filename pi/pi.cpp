#define _USE_MATH_DEFINES  // for C++
#include <iostream> // std::cin and std::cout
#include <iomanip> // std::fixed and std::setprecision
#include <cmath> // M_PI
#include <string> // std::stoi
#include <chrono> // time
long double pi_rectangle(long double x, long double width) {
  const long double one = 1.0;
  return (one / (one + x * x)) * width;
}

void print_sum_and_pi(long double sum, const std::chrono::duration<double>& seconds, bool is_automatic) {
  if (is_automatic){
    std::cout<<static_cast<long double>(4.0) * sum<<" "<<fabsl(M_PI - 4.0 * sum)<<" "<<seconds.count()<<"\n";
  }else{
    std::cout <<"Result:     "<< static_cast<long double>(4.0) * sum << "\n";
    std::cout <<"Inaccuracy: "<< fabsl(M_PI - 4.0 * sum) << "\n";
    std::cout <<"Seconds:    " << seconds.count()<<"\n";
    std::cout<<"\n";
  }
}

/// 1 000 000 000 operations make inaccuracy 0 in every test
/// 1 000 000 operations make inaccuracy 0 in middle rectangles

enum class counting_type{
  left,
  middle,
  right,
  all
};

int main(int argc, char* argv[]) {
  const int N = argc > 1 ? std::stoi(argv[1]) : 1000;
  bool is_automatic = false;
  counting_type type;
  if(argc>2){
    is_automatic = std::string(argv[2]) == std::string("a") ? true : false; // it is strange but argv[2] != "a" && argv[2] != "a\0"
    int mode = std::stoi(std::string(argv[3]));
    switch (mode)
    {
    case 0:
      type = counting_type::left;
      /* code */
      break;
    case 1:
      type = counting_type::middle;
          break;
    case 2:
      type = counting_type::right;
          break;
    case 3:
    type = counting_type::all;
          break;
    default:
      break;
    }
  }
  long double step = static_cast<long double>(1.0) / N;
  long double sum = 0.0;
  const int precision = 10;
  std::cout << std::fixed; // to set precision to every value
  std::cout << std::setprecision(precision);
  std::chrono::duration<double> elapsed_seconds;

  if(!is_automatic||is_automatic&&(type==counting_type::all||type==counting_type::left)){
    const auto start_left{std::chrono::steady_clock::now()};
    for (int i = 0; i < N; i++) {
      sum += pi_rectangle(static_cast<long double>(i) * step, step);
    }
    const auto end_left{std::chrono::steady_clock::now()};
    elapsed_seconds = end_left - start_left;
    print_sum_and_pi(sum, elapsed_seconds, is_automatic);
  }


  sum = 0.0;
    if(!is_automatic||is_automatic&&(type==counting_type::all||type==counting_type::middle)){
    const auto start_middle{std::chrono::steady_clock::now()};
    for (int i = 0; i < N; i++) {
      sum += pi_rectangle((i * step + (i + 1) * step) * 0.5, step);
    }
      const auto end_middle{std::chrono::steady_clock::now()};
      elapsed_seconds = end_middle - start_middle;
      print_sum_and_pi(sum, elapsed_seconds, is_automatic);
  }


  sum = 0.0;
  if(!is_automatic||is_automatic&&(type==counting_type::all||type==counting_type::right)){
    const auto start_right{std::chrono::steady_clock::now()};
    for (int i = 0; i < N; i++) {
      sum += pi_rectangle((static_cast<long double>(i) + 1) * step, step);
  }
    const auto end_right{std::chrono::steady_clock::now()};
    elapsed_seconds = end_right - start_right;
    print_sum_and_pi(sum, elapsed_seconds, is_automatic);
  }

  return 0;
}

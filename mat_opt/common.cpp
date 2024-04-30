#include "common.hpp"

using std::string;

void generate_rand_matrix(double* matr, int n, double min, double max){
    std::random_device rd;
    std::mt19937 engine(rd());
    std::uniform_real_distribution<double> gen(min, max);
    const int n2 = n*n;
    for(int i=0;i<n2;i++){
        matr[i] = gen(engine);
    }
}
void generate_zero_matrix(double* matr, int n){
    memset(matr, 0, sizeof(double)*n*n);
}
void generate_identity_matrix(double* matr, int n){
    generate_zero_matrix(matr, n);
    for(int i=0;i<n;i++){
        matr[i*n+i]=1.0;
    }
}

void get_cache_sizes(std::string filename, int* cache_sizes){
    std::ifstream input_file(filename);
    if(input_file.is_open()){
        std::string line;
        int i=0;
        while(std::getline(input_file, line)){
            cache_sizes[i] = std::stoi(line);
            i++;
        }
    }else{
        throw(std::string("File \"")+filename+std::string("\" has not been opened\n"));
    }
}

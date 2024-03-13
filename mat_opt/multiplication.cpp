#include "multiplication.hpp"
#include <iostream>

void base_block_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    int block_n;
    int tmp_block_n = n;
    while(tmp_block_n%2==0){
        tmp_block_n /= 2;
    }
    if(tmp_block_n>1){
        block_n = tmp_block_n;
    }else{
        block_n = 2; // unsafe part
    }

    double* __restrict__ a_block = new double[n*n];
    double* __restrict__ b_block = new double[n*n];
    blockcpy(a, a_block, n, block_n);
    blockcpy(b, b_block, n, block_n);
    //////////////// multiplication itself
}
void optimal_block_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){}
void b_transposed_block_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){}

void base_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    for(int i=0;i<n;i++){ // i-th row in a
        for(int j=0;j<n;j++){ // j-th column in b
            for(int k=0;k<n; k++){ // k-th element in vector
                c[i*n+j]+=a[i*n+k]*b[k*n+j];
            }
        }
    }
}
void base_matrix_mult_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    #pragma omp parallel for shared(a, b, c, n)
        for(int i=0;i<n;i++){ // i-th row in a
            for(int j=0;j<n;j++){ // j-th column in b
                for(int k=0;k<n; k++){ // k-th element in vector
                    c[i*n+j]+=a[i*n+k]*b[k*n+j];
                }
            }
        }
}
void row_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            for(int k=0;k<n;k++){
                c[i*n+k]+=a[i*n+j]*b[j*n+k];
            }
        }
    }
}
void row_matrix_mult_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    #pragma omp parallel for shared(a, b, c, n)
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                for(int k=0;k<n;k++){
                    c[i*n+k]+=a[i*n+j]*b[j*n+k];
                }
            }
        }
}
void b_transposed_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            for(int k=0;k<n;k++){
                c[i*n+j] += a[i*n+k]*b[j*n+k];
            }
        }
    }
}
void b_transposed_matrix_mult_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    #pragma omp parallel for shared(a, b, c, n)
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                for(int k=0;k<n;k++){
                    c[i*n+j] += a[i*n+k]*b[j*n+k];
                }
            }
        }
}
// void b_transposed_matrix_mult_omp_simd(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
//     #pragma omp parallel for shared(a, b, c, n)
//         for(int i=0;i<n;i++){
//             for(int j=0;j<n;j++){
//                 c[i*n+j]=std::inner_product(a, a+i*n, b, 0, std::multiplies<double>(), std::plus<double>());
//                 // for(int k=0;k<n;k++){
//                 //     c[i*n+j] += a[i*n+k]*b[j*n+k];
//                 // }
//             }
//         }
// }
void b_transposed_matrix_mult_omp_simd(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    #pragma omp parallel
    {
        #pragma omp for simd
            for(int i=0;i<n;i++){
                for(int k=0;k<n;k++){
                    for(int j=0;j<n;j++){
                        c[i*n+j] += a[i*n+k]*b[j*n+k];
                    }
                }
            }
    }
}
void transposed_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    double* bT = get_transposed_matrix(b, n);
    b_transposed_matrix_mult(a, bT, c, n);
}
void transposed_matrix_mult_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    double* bT = get_transposed_matrix(b, n);
    b_transposed_matrix_mult_omp(a, bT, c, n);
}
void transposed_matrix_mult_omp_simd(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    double* bT = get_transposed_matrix(b, n);
    b_transposed_matrix_mult_omp_simd(a, bT, c, n);
}

void strassen_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){  
    double* inc_a, *inc_b, *inc_c;

    // int inc_n = increase_matrices(a, inc_a, b, inc_b, c, inc_c, n);
        int inc_dim = 1;
    int offset = log2(n);
    if(offset<log2(n)){
        offset++;
    }
    inc_dim <<= offset;

    inc_a = new double[inc_dim*inc_dim];
    inc_b = new double[inc_dim*inc_dim];
    inc_c = new double[inc_dim*inc_dim];
    
    generate_zero_matrix(inc_a, inc_dim);
    generate_zero_matrix(inc_b, inc_dim);
    generate_zero_matrix(inc_c, inc_dim);

    for(int i=0;i<n;i++){
        std::memcpy(inc_a+i*inc_dim, a+i*n, n*sizeof(double));
        std::memcpy(inc_b+i*inc_dim, b+i*n, n*sizeof(double));
        // there is no memcpy for inc_c because c is zero matrix
    }
    recursive_strassen_part(inc_a, inc_b, inc_c, inc_dim);
    decrease_matrix(c, inc_c, n, inc_dim);
}
void recursive_strassen_part(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    // std::cout<<"n_start = "<<n<<"\n";
    if(n<=64){
        row_matrix_mult_omp(a, b, c, n);
        // for(int i=0;i<n;i++){
        //     for(int j=0;j<n;j++){
        //         std::cout<<c[i*n+j]<<" ";
        //     }
        //     std::cout<<"\n";
        // }
        return;
    }
    n >>= 1;
    double *a11 = new double[n*n], *a22 = new double[n*n], *a12 = new double[n*n], *a21 = new double[n*n];
    double *b11 = new double[n*n], *b12 = new double[n*n], *b21 = new double[n*n], *b22 = new double[n*n];
    split_matrices(a, a11, a12, a21, a22, n);
    split_matrices(b, b11, b12, b21, b22, n);
    //double *p1, *p2, *p3, *p4, *p5, *p6, *p7;
    double** p = new double*[7];
    for(int i=0;i<7;i++){
        p[i] = new double[n*n];
        generate_zero_matrix(p[i], n);
    }

    double* a11_add_a22, *b11_add_b22, *a21_add_a22, *b12_sub_b22, *b21_sub_b11, *a11_add_a12, *a21_sub_a11, *b11_add_b12, *a12_sub_a22, *b21_add_b22;
    a11_add_a22 = new double[n*n];
    b11_add_b22 = new double[n*n];
    a21_add_a22 = new double[n*n];
    b12_sub_b22 = new double[n*n];
    b21_sub_b11 = new double[n*n];
    a11_add_a12 = new double[n*n];
    a21_sub_a11 = new double[n*n];
    b11_add_b12 = new double[n*n];
    a12_sub_a22 = new double[n*n];
    b21_add_b22 = new double[n*n];

    // std::cout<<"n = "<<n<<"\n";
    matrix_add(a11, a22, a11_add_a22, n);
    matrix_add(a21, a22, a21_add_a22, n);
    matrix_sub(b12, b22, b12_sub_b22, n);
    matrix_sub(b21, b11, b21_sub_b11, n);
    matrix_add(a11, a12, a11_add_a12, n);
    matrix_sub(a21, a11, a21_sub_a11, n);
    matrix_add(b11, b12, b11_add_b12, n);
    matrix_sub(a12, a22, a12_sub_a22, n);
    matrix_add(b21, b22, b21_add_b22, n);
    // std::cout<<"n_recursive = "<<n<<"\n";

    

    // recursive_strassen_part(a11_add_a22, b11_add_b22, p1, n);
    // recursive_strassen_part(a21_add_a22, b11, p2, n);
    // recursive_strassen_part(a11, b12_sub_b22, p3, n);
    // recursive_strassen_part(a22, b21_sub_b11, p4, n);
    // recursive_strassen_part(a11_add_a12, b22, p5, n);
    // recursive_strassen_part(a21_sub_a11, b11_add_b12, p6, n);
    // recursive_strassen_part(a12_sub_a22, b21_add_b22, p7, n);
    recursive_strassen_part(a11_add_a22, b11_add_b22, p[0], n);
    recursive_strassen_part(a21_add_a22, b11, p[1], n);
    recursive_strassen_part(a11, b12_sub_b22, p[2], n);
    recursive_strassen_part(a22, b21_sub_b11, p[3], n);
    recursive_strassen_part(a11_add_a12, b22, p[4], n);
    recursive_strassen_part(a21_sub_a11, b11_add_b12, p[5], n);
    recursive_strassen_part(a12_sub_a22, b21_add_b22, p[6], n);
    // std::cout<<"n1 = "<<n<<"\n";
    double *c11 = new double[n*n], *c12 = new double[n*n], *c21 = new double[n*n], *c22 = new double[n*n];
    // double** c_ij = new double*[4]; // c11, c12, c21, c22
    // for(int i=0;i<4;i++){
    //     c_ij[i] = new double[n*n];
    //     generate_zero_matrix(c_ij[i],n);
    // }

    double* p1_sum_p4 = new double[n*n], *p7_sub_p5 = new double[n*n], *p1_sub_p2 = new double[n*n], *p3_sum_p6 = new double[n*n]; // can be reduced to reusing c_ij
    // matrix_add(p1_sum_p4, p7_sub_p5, c11, n);
    // matrix_add(p3, p5, c11, n);
    // matrix_add(p2, p4, c21, n);
    // matrix_add(p1_sub_p2, p3_sum_p6, c22, n);
    matrix_add(p[0], p[3], p1_sum_p4, n);
    matrix_sub(p[0], p[1], p1_sub_p2, n);
    matrix_add(p[2], p[5], p3_sum_p6, n);

    matrix_add(p1_sum_p4, p7_sub_p5, c11, n);
    matrix_add(p[2], p[4], c12, n);
    matrix_add(p[1], p[3], c21, n);
    matrix_add(p1_sub_p2, p3_sum_p6, c22, n);
    // matrix_add(p1_sum_p4, p7_sub_p5, c_ij[0], n);
    // matrix_add(p[2], p[4], c_ij[1], n);
    // matrix_add(p[1], p[3], c_ij[2], n);
    // matrix_add(p1_sub_p2, p3_sum_p6, c_ij[3], n);

    // std::cout<<"n2 = "<<n<<"\n";
    collect_matrices(c, c11, c12, c21, c22, n<<1); // may be delete c_ij here
    // collect_matrices(c, c_ij[0], c_ij[1], c_ij[2], c_ij[3], n<<1); // may be delete c_ij here

}

//     private static int[][] multiStrassen(int[][] a, int[][] b, int n) {
//     if (n <= 64) {
//         return multiply(a, b);
//     }

//     n = n >> 1;

//     int[][] a11 = new int[n][n];
//     int[][] a12 = new int[n][n];
//     int[][] a21 = new int[n][n];
//     int[][] a22 = new int[n][n];

//     int[][] b11 = new int[n][n];
//     int[][] b12 = new int[n][n];
//     int[][] b21 = new int[n][n];
//     int[][] b22 = new int[n][n];

//     splitMatrix(a, a11, a12, a21, a22);
//     splitMatrix(b, b11, b12, b21, b22);

//     int[][] p1 = multiStrassen(summation(a11, a22), summation(b11, b22), n);
//     int[][] p2 = multiStrassen(summation(a21, a22), b11, n);
//     int[][] p3 = multiStrassen(a11, subtraction(b12, b22), n);
//     int[][] p4 = multiStrassen(a22, subtraction(b21, b11), n);
//     int[][] p5 = multiStrassen(summation(a11, a12), b22, n);
//     int[][] p6 = multiStrassen(subtraction(a21, a11), summation(b11, b12), n);
//     int[][] p7 = multiStrassen(subtraction(a12, a22), summation(b21, b22), n);

//     int[][] c11 = summation(summation(p1, p4), subtraction(p7, p5));
//     int[][] c12 = summation(p3, p5);
//     int[][] c21 = summation(p2, p4);
//     int[][] c22 = summation(subtraction(p1, p2), summation(p3, p6));

//     return collectMatrix(c11, c12, c21, c22);
// }

void transpose_matrix_in_place(double* matr, int n){
    for(int i=0;i<n;i++){
        for(int j=i;j<n;j++){
            std::swap(matr[i * n + j], matr[j * n + i]);
        }
    }
}
double* get_transposed_matrix(double* matrix, int n){
    double* Tmatrix = new double[n*n];
    std::memcpy(Tmatrix, matrix, n*n*sizeof(double));
    transpose_matrix_in_place(Tmatrix, n);
    return Tmatrix;
}
void blockcpy(double* __restrict__ src, double* __restrict__ dest, int n, int block_n){

    const int block_line = n / block_n; // block number in src matrix line
    const int block_size = block_n * block_n; // element number in block
    const int block_line_size = block_n * block_n * block_line; // element number in block line

    for(int i=0;i<block_line;i++){ // i-th block row
        for(int j=0;j<block_line;j++){ // j-th block column
            for(int k=0;k<block_n;k++){ // k-th row in block
                std::memcpy(dest+i*block_line_size+j*block_size+k*block_n, src+i*n*block_n+k*n+j*block_n, block_n*sizeof(double));
            }
        }
    }
}

int increase_matrices(double* __restrict__ a, double* __restrict__ inc_a, double* __restrict__ b, double* __restrict__ inc_b, double* __restrict__ c, double* __restrict__  inc_c, int n){
    int inc_dim = 1;
    int offset = log2(n);
    if(offset<log2(n)){
        offset++;
    }
    inc_dim <<= offset;
    std::cout<<n<<" "<<inc_dim<<"\n";

    inc_a = new double[inc_dim*inc_dim];
    inc_b = new double[inc_dim*inc_dim];
    inc_c = new double[inc_dim*inc_dim];
    
    generate_zero_matrix(inc_a, inc_dim);
    generate_zero_matrix(inc_b, inc_dim);
    generate_zero_matrix(inc_c, inc_dim);

    for(int i=0;i<n;i++){
        std::memcpy(inc_a+i*inc_dim, a+i*n, n*sizeof(double));
        std::memcpy(inc_b+i*inc_dim, b+i*n, n*sizeof(double));
        // there is no memcpy for inc_c because c is zero matrix
    }
    return inc_dim;
}

void decrease_matrix(double* __restrict__ a, double* __restrict__ inc_a, int n, int inc_n){
    for(int i=0;i<n;i++){
        std::memcpy(a+i*n, inc_a+i*inc_n, n*sizeof(double));
    }
}

void split_matrices(double* __restrict__ a, double* __restrict__ a11, double* __restrict__ a12,double* __restrict__ a21, double* __restrict__ a22, int n){
    const int half_n = n>>1;
    for(int i=0;i<half_n;i++){
        std::memcpy(a11+i*half_n, a+i*n, half_n);
        std::memcpy(a12+i*half_n, a+i*(n+half_n), half_n);
    }
    for(int i=half_n;i<n;i++){
        std::memcpy(a21+i*half_n, a+i*n, half_n);
        std::memcpy(a22+i*half_n, a+i*(n+half_n), half_n);
    }
}
void collect_matrices(double* __restrict__ a, double* __restrict__ a11, double* __restrict__ a12,double* __restrict__ a21, double* __restrict__ a22, int n){
    const int half_n = n>>1;
    for(int i=0;i<half_n;i++){
        std::memcpy( a+i*n, a11+i*half_n, half_n);
        std::memcpy(a+i*(n+half_n), a12+i*half_n, half_n);
    }
    for(int i=half_n;i<n;i++){
        std::memcpy(a+i*n, a21+i*half_n, half_n);
        std::memcpy(a+i*(n+half_n), a22+i*half_n, half_n);
    }
}
void matrix_add(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            c[i*n+j]=a[i*n+j]+b[i*n+j];
        }
    }
}
void matrix_sub(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            c[i*n+j]=a[i*n+j]-b[i*n+j];
        }
    }
}

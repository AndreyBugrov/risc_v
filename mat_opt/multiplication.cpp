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

///
/// ADDING SIMD DIRECORY ALWAYS CAUSES WORSE PERFORMANCE!!!///
///

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
    omp_set_num_threads(4);
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
void row_matrix_mult_simd(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            #pragma omp simd
            for(int k=0;k<n;k++){
                c[i*n+k]+=a[i*n+j]*b[j*n+k];
            }
        }
    }
}
void row_matrix_mult_opt(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    for(int i=0;i<n;i++){
        double* c_i = c+i*n; // do not delete pointers because the memory can be used in the other functions
        for(int j=0;j<n;j++){
            double* b_j = b+j*n;
            double a_ij = a[i*n+j];
            for(int k=0;k<n;k++){
                c_i[k]+=a_ij*b_j[k];
            }
        }
    }
}
void row_matrix_mult_opt_simd(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){  
    for(int i=0;i<n;i++){
        double* c_i = c+i*n; // do not delete pointers because the memory can be used in the other functions
        for(int j=0;j<n;j++){
            double* b_j = b+j*n;
            double a_ij = a[i*n+j];
            #pragma omp for simd
                for(int k=0;k<n;k++){
                    c_i[k]+=a_ij*b_j[k];
                }
        }
    }
}

void row_matrix_mult_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    omp_set_num_threads(4);
    #pragma omp parallel for shared(a, b, c, n)
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                for(int k=0;k<n;k++){
                    c[i*n+k]+=a[i*n+j]*b[j*n+k];
                }
            }
        }
}
void row_matrix_mult_opt_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    omp_set_num_threads(4);
    #pragma omp parallel for shared(a, b, c, n)
        for(int i=0;i<n;i++){
            double* c_i = c+i*n; // do not delete pointers because the memory can be used in the other functions
            for(int j=0;j<n;j++){
                double* b_j = b+j*n;
                double a_ij = a[i*n+j];
                for(int k=0;k<n;k++){
                    c_i[k]+=a_ij*b_j[k];
                }
            }
        }
}
void row_matrix_mult_omp_simd(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    #pragma omp parallel for simd collapse(3)
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            for(int k=0;k<n;k++){
                c[i*n+k]+=a[i*n+j]*b[j*n+k];
            }
        }
    }
}
void row_matrix_mult_opt_omp_simd(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    omp_set_num_threads(4);
    #pragma omp parallel
    {
        for(int i=0;i<n;i++){
            double* c_i = c+i*n; // do not delete pointers because the memory can be used in the other functions
            for(int j=0;j<n;j++){
                double* b_j = b+j*n;
                double a_ij = a[i*n+j];
                #pragma omp for simd
                    for(int k=0;k<n;k++){
                        c_i[k]+=a_ij*b_j[k];
                    }
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
    omp_set_num_threads(4);
    #pragma omp parallel for shared(a, b, c, n)
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                for(int k=0;k<n;k++){
                    c[i*n+j] += a[i*n+k]*b[j*n+k];
                }
            }
        }
}
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
void base_strassen_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n, int opt_type){
    double* inc_a, *inc_b, *inc_c;

    int inc_n = 1;
    int offset = log2(n);
    if(offset<log2(n)){
        offset++;
    }
    inc_n <<= offset;

    inc_a = new double[inc_n*inc_n];
    inc_b = new double[inc_n*inc_n];
    inc_c = new double[inc_n*inc_n];
    
    generate_zero_matrix(inc_a, inc_n);
    generate_zero_matrix(inc_b, inc_n);
    generate_zero_matrix(inc_c, inc_n);

    for(int i=0;i<n;i++){
        std::memcpy(inc_a+i*inc_n, a+i*n, n*sizeof(double));
        std::memcpy(inc_b+i*inc_n, b+i*n, n*sizeof(double));
        // there is no memcpy for inc_c because c is zero matrix
    }
    switch (opt_type)
    {
    case 0:
        recursive_strassen_part(inc_a, inc_b, inc_c, inc_n);
        break;
    case 1:
        recursive_strassen_part_omp(inc_a, inc_b, inc_c, inc_n);
        break;
    case 2:
        recursive_strassen_part_rec_omp(inc_a, inc_b, inc_c, inc_n);
        break;
    
    default:
        break;
    }
    decrease_matrix(c, inc_c, n, inc_n);
}
void strassen_matrix_mult(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){  
    base_strassen_matrix_mult(a, b, c, n, 0);
}
void strassen_matrix_mult_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){  
    base_strassen_matrix_mult(a, b, c, n, 1);
}
void strassen_matrix_mult_rec_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    base_strassen_matrix_mult(a, b, c, n, 2);
}
void recursive_strassen_part(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    if(n<=kRecursiveStrassenMultLimit){
        row_matrix_mult(a, b, c, n);
        return;
    }
    n >>= 1;
    double *a11, *a12, *a21,*a22;
    double *b11, *b12, *b21, *b22;
    a11 = new double[n*n];
    a12 = new double[n*n];
    a21 = new double[n*n];
    a22 = new double[n*n];
    b11 = new double[n*n];
    b12 = new double[n*n];
    b21 = new double[n*n];
    b22 = new double[n*n];
    split_matrices(a, a11, a12, a21, a22, 2*n);
    split_matrices(b, b11, b12, b21, b22, 2*n);
    double** p = new double*[7];
    for(int i=0; i<7; i++){
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

    matrix_add(a11, a22, a11_add_a22, n);
    matrix_add(b11, b22, b11_add_b22, n);
    matrix_add(a21, a22, a21_add_a22, n);
    matrix_sub(b12, b22, b12_sub_b22, n);
    matrix_sub(b21, b11, b21_sub_b11, n);
    matrix_add(a11, a12, a11_add_a12, n);
    matrix_sub(a21, a11, a21_sub_a11, n);
    matrix_add(b11, b12, b11_add_b12, n);
    matrix_sub(a12, a22, a12_sub_a22, n);
    matrix_add(b21, b22, b21_add_b22, n);

    recursive_strassen_part(a11_add_a22, b11_add_b22, p[0], n); // + D
    recursive_strassen_part(a21_add_a22, b11, p[1], n); // + H2
    recursive_strassen_part(a11, b12_sub_b22, p[2], n); // + V2
    recursive_strassen_part(a22, b21_sub_b11, p[3], n); // + V1
    recursive_strassen_part(a11_add_a12, b22, p[4], n); // + H1
    recursive_strassen_part(a21_sub_a11, b11_add_b12, p[5], n); // + D2
    recursive_strassen_part(a12_sub_a22, b21_add_b22, p[6], n); // + D1

    delete[] a11_add_a22;
    delete[] b11_add_b22;
    delete[] a21_add_a22;
    delete[] b12_sub_b22;
    delete[] b21_sub_b11;
    delete[] a11_add_a12;
    delete[] a21_sub_a11;
    delete[] b11_add_b12;
    delete[] a12_sub_a22;
    delete[] b21_add_b22;

    double *c11 = new double[n*n], *c12 = new double[n*n], *c21 = new double[n*n], *c22 = new double[n*n];
    double* p1_add_p4 = new double[n*n], *p7_sub_p5 = new double[n*n], *p1_sub_p2 = new double[n*n], *p3_sum_p6 = new double[n*n]; // can be reduced to reusing c_ij

    matrix_add(p[0], p[3], p1_add_p4, n);
    matrix_sub(p[6], p[4], p7_sub_p5, n);
    matrix_sub(p[0], p[1], p1_sub_p2, n);
    matrix_add(p[2], p[5], p3_sum_p6, n);

    matrix_add(p1_add_p4, p7_sub_p5, c11, n);
    matrix_add(p[2], p[4], c12, n);
    matrix_add(p[1], p[3], c21, n);
    matrix_add(p1_sub_p2, p3_sum_p6, c22, n);

    for(int i=0;i<7;i++){
        delete[] p[i];
    }
    delete p;

    collect_matrices(c, c11, c12, c21, c22, 2*n); // may be delete c_ij here
    delete[] c11;
    delete[] c12;
    delete[] c21;
    delete[] c22;
}

void recursive_strassen_part_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    if(n<=kRecursiveStrassenMultLimit){
        row_matrix_mult_omp(a, b, c, n);
        return;
    }
    n >>= 1;
    double *a11, *a12, *a21,*a22;
    double *b11, *b12, *b21, *b22;
    a11 = new double[n*n];
    a12 = new double[n*n];
    a21 = new double[n*n];
    a22 = new double[n*n];
    b11 = new double[n*n];
    b12 = new double[n*n];
    b21 = new double[n*n];
    b22 = new double[n*n];
    split_matrices(a, a11, a12, a21, a22, 2*n);
    split_matrices(b, b11, b12, b21, b22, 2*n);
    double** p = new double*[7];
    for(int i=0; i<7; i++){
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

    matrix_add_omp(a11, a22, a11_add_a22, n);
    matrix_add_omp(b11, b22, b11_add_b22, n);
    matrix_add_omp(a21, a22, a21_add_a22, n);
    matrix_sub_omp(b12, b22, b12_sub_b22, n);
    matrix_sub_omp(b21, b11, b21_sub_b11, n);
    matrix_add_omp(a11, a12, a11_add_a12, n);
    matrix_sub_omp(a21, a11, a21_sub_a11, n);
    matrix_add_omp(b11, b12, b11_add_b12, n);
    matrix_sub_omp(a12, a22, a12_sub_a22, n);
    matrix_add_omp(b21, b22, b21_add_b22, n);

    recursive_strassen_part_omp(a11_add_a22, b11_add_b22, p[0], n); // + D
    recursive_strassen_part_omp(a21_add_a22, b11, p[1], n); // + H2
    recursive_strassen_part_omp(a11, b12_sub_b22, p[2], n); // + V2
    recursive_strassen_part_omp(a22, b21_sub_b11, p[3], n); // + V1
    recursive_strassen_part_omp(a11_add_a12, b22, p[4], n); // + H1
    recursive_strassen_part_omp(a21_sub_a11, b11_add_b12, p[5], n); // + D2
    recursive_strassen_part_omp(a12_sub_a22, b21_add_b22, p[6], n); // + D1

    delete[] a11_add_a22;
    delete[] b11_add_b22;
    delete[] a21_add_a22;
    delete[] b12_sub_b22;
    delete[] b21_sub_b11;
    delete[] a11_add_a12;
    delete[] a21_sub_a11;
    delete[] b11_add_b12;
    delete[] a12_sub_a22;
    delete[] b21_add_b22;

    double *c11 = new double[n*n], *c12 = new double[n*n], *c21 = new double[n*n], *c22 = new double[n*n];
    double* p1_add_p4 = new double[n*n], *p7_sub_p5 = new double[n*n], *p1_sub_p2 = new double[n*n], *p3_sum_p6 = new double[n*n]; // can be reduced to reusing c_ij

    matrix_add_omp(p[0], p[3], p1_add_p4, n);
    matrix_sub_omp(p[6], p[4], p7_sub_p5, n);
    matrix_sub_omp(p[0], p[1], p1_sub_p2, n);
    matrix_add_omp(p[2], p[5], p3_sum_p6, n);

    matrix_add_omp(p1_add_p4, p7_sub_p5, c11, n);
    matrix_add_omp(p[2], p[4], c12, n);
    matrix_add_omp(p[1], p[3], c21, n);
    matrix_add_omp(p1_sub_p2, p3_sum_p6, c22, n);

    for(int i=0;i<7;i++){
        delete[] p[i];
    }
    delete p;

    collect_matrices(c, c11, c12, c21, c22, 2*n); // may be delete c_ij here
    delete[] c11;
    delete[] c12;
    delete[] c21;
    delete[] c22;
}
void recursive_strassen_part_rec_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
        if(n<=kRecursiveStrassenMultLimit){
        row_matrix_mult(a, b, c, n);
        return;
    }
    n >>= 1;
    double *a11, *a12, *a21,*a22;
    double *b11, *b12, *b21, *b22;
    a11 = new double[n*n];
    a12 = new double[n*n];
    a21 = new double[n*n];
    a22 = new double[n*n];
    b11 = new double[n*n];
    b12 = new double[n*n];
    b21 = new double[n*n];
    b22 = new double[n*n];
    split_matrices(a, a11, a12, a21, a22, 2*n);
    split_matrices(b, b11, b12, b21, b22, 2*n);
    double** p = new double*[7];
    for(int i=0; i<7; i++){
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

    matrix_add(a11, a22, a11_add_a22, n);
    matrix_add(b11, b22, b11_add_b22, n);
    matrix_add(a21, a22, a21_add_a22, n);
    matrix_sub(b12, b22, b12_sub_b22, n);
    matrix_sub(b21, b11, b21_sub_b11, n);
    matrix_add(a11, a12, a11_add_a12, n);
    matrix_sub(a21, a11, a21_sub_a11, n);
    matrix_add(b11, b12, b11_add_b12, n);
    matrix_sub(a12, a22, a12_sub_a22, n);
    matrix_add(b21, b22, b21_add_b22, n);

    omp_set_num_threads(4);
    #pragma omp parallel sections
    {
        #pragma omp section
        {
            recursive_strassen_part(a11_add_a22, b11_add_b22, p[0], n); // + D
        }
        #pragma omp section
        {
            recursive_strassen_part(a21_add_a22, b11, p[1], n); // + H2
        }
        #pragma omp section 
        {
            recursive_strassen_part(a11, b12_sub_b22, p[2], n); // + V2
        }
        #pragma omp section
        {
            recursive_strassen_part(a22, b21_sub_b11, p[3], n); // + V1
        }
        #pragma omp section
        {
            recursive_strassen_part(a11_add_a12, b22, p[4], n); // + H1
        }
        #pragma omp section
        {
            recursive_strassen_part(a21_sub_a11, b11_add_b12, p[5], n); // + D2
        }
        #pragma omp section
        {
            recursive_strassen_part(a12_sub_a22, b21_add_b22, p[6], n); // + D1
        }
    }
    delete[] a11_add_a22;
    delete[] b11_add_b22;
    delete[] a21_add_a22;
    delete[] b12_sub_b22;
    delete[] b21_sub_b11;
    delete[] a11_add_a12;
    delete[] a21_sub_a11;
    delete[] b11_add_b12;
    delete[] a12_sub_a22;
    delete[] b21_add_b22;

    double *c11 = new double[n*n], *c12 = new double[n*n], *c21 = new double[n*n], *c22 = new double[n*n];
    double* p1_add_p4 = new double[n*n], *p7_sub_p5 = new double[n*n], *p1_sub_p2 = new double[n*n], *p3_sum_p6 = new double[n*n]; // can be reduced to reusing c_ij

    matrix_add(p[0], p[3], p1_add_p4, n);
    matrix_sub(p[6], p[4], p7_sub_p5, n);
    matrix_sub(p[0], p[1], p1_sub_p2, n);
    matrix_add(p[2], p[5], p3_sum_p6, n);

    matrix_add(p1_add_p4, p7_sub_p5, c11, n);
    matrix_add(p[2], p[4], c12, n);
    matrix_add(p[1], p[3], c21, n);
    matrix_add(p1_sub_p2, p3_sum_p6, c22, n);

    for(int i=0;i<7;i++){
        delete[] p[i];
    }
    delete p;

    collect_matrices(c, c11, c12, c21, c22, 2*n); // may be delete c_ij here
    delete[] c11;
    delete[] c12;
    delete[] c21;
    delete[] c22;
}

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

void decrease_matrix(double* __restrict__ a, double* __restrict__ inc_a, int n, int inc_n){
    for(int i=0;i<n;i++){
        std::memcpy(a+i*n, inc_a+i*inc_n, n*sizeof(double));
    }
}

void split_matrices(double* __restrict__ a, double* __restrict__ a11, double* __restrict__ a12,double* __restrict__ a21, double* __restrict__ a22, int n_a){
    const int half_n = n_a>>1;
    for(int i = 0; i < half_n; i++){
        std::memcpy(a11+i*half_n, a+i*n_a, half_n*sizeof(double));
    }
    for(int i = 0; i < half_n; i++){
        std::memcpy(a12+i*half_n, a+i*n_a+half_n, half_n*sizeof(double));
    }
    for(int i = 0; i < half_n; i++){
        std::memcpy(a21+i*half_n, a+(i+half_n)*n_a, half_n*sizeof(double));
    }
    for(int i = 0; i < half_n; i++){
        std::memcpy(a22+i*half_n, a+(i+half_n)*n_a+half_n, half_n*sizeof(double));
    }
}
void collect_matrices(double* __restrict__ a, double* __restrict__ a11, double* __restrict__ a12,double* __restrict__ a21, double* __restrict__ a22, int n_a){
    const int half_n = n_a>>1;
    for(int i = 0; i < half_n; i++){
        std::memcpy(a+i*n_a, a11+i*half_n, half_n*sizeof(double));
    }
    for(int i = 0; i < half_n; i++){
        std::memcpy(a+i*n_a+half_n, a12+i*half_n, half_n*sizeof(double));
    }
    for(int i = 0; i < half_n; i++){
        std::memcpy(a+(i+half_n)*n_a, a21+i*half_n, half_n*sizeof(double));
    }
    for(int i = 0; i < half_n; i++){
        std::memcpy(a+(i+half_n)*n_a+half_n, a22+i*half_n, half_n*sizeof(double));
    }
}
void matrix_add(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            c[i*n+j]=a[i*n+j]+b[i*n+j];
        }
    }
}
void matrix_add_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    omp_set_num_threads(4);
    #pragma omp parallel for shared(a, b, c, n)
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
void matrix_sub_omp(double* __restrict__ a, double* __restrict__ b, double* __restrict__ c, int n){
    omp_set_num_threads(4);
    #pragma omp parallel for shared(a, b, c, n)
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                c[i*n+j]=a[i*n+j]-b[i*n+j];
            }
        }
}

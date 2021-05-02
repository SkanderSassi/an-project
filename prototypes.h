#ifndef PROTOTYPES_H
#define PROTOTYPES_H



// Generators
double* generate_matrix_double(int dim, int low, int high, int verbose);
float* generate_matrix_float(int dim, int low, int high, int verbose);
// Printers
void print_matrix_float(int dim, float* matrix);
void print_matrix_double(int dim, double* matrix);
// Sanity check
int check_equal_matrix_double(int dim, double* matA, double* matB);
int check_equal_matrix_float(int dim, float* matA, float* matB);
// Multiplication
double* ijk_double(int dim, double* matA, double* matB, int verbose);
float* ijk_float(int dim, float* matA, float* matB, int verbose);

double* jik_double(int dim, double* matA, double* matB, int verbose);
float* jik_float(int dim, float* matA, float* matB, int verbose);

//Tests

void test_double_matrices(int dim, double* matA, double* matB);
void test_float_matrices(int dim, float* matA, float* matB);

#endif
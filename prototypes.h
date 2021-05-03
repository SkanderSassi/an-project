#ifndef PROTOTYPES_H
#define PROTOTYPES_H


struct FloatResult{
    float *matrix;
    int trial_id;
    double time_spent;
    char alg[3];
};

struct DoubleResult{
    double *matrix;
    int trial_id;
    double time_spent;
    char alg[3];
};

typedef struct FloatResult FloatResult;
typedef struct DoubleResult DoubleResult;


// Generators
double* generate_matrix_double(int dim, int low, int high, int verbose);
float* generate_matrix_float(int dim, int low, int high, int verbose);
// Printers
void print_matrix_float(int dim, float* matrix);
void print_matrix_double(int dim, double* matrix);

void print_float_result(int dim, FloatResult res_struct, int show_matrix, int show_time, int show_trial_id, int show_alg);
void print_double_result(int dim, DoubleResult res_struct, int show_matrix, int show_time, int show_trial_id, int show_alg);

// Sanity check
int check_equal_matrix_double(int dim, double* matA, double* matB);
int check_equal_matrix_float(int dim, float* matA, float* matB);

// Multiplication
DoubleResult ijk_double(int dim, double* matA, double* matB, int verbose);
FloatResult ijk_float(int dim, float* matA, float* matB, int verbose);

DoubleResult jik_double(int dim, double* matA, double* matB, int verbose);
FloatResult jik_float(int dim, float* matA, float* matB, int verbose);

DoubleResult kij_double(int dim, double* matA, double* matB, int verbose);
FloatResult kij_float(int dim, float* matA, float* matB, int verbose);

DoubleResult ikj_double(int dim, double* matA, double* matB, int verbose);
FloatResult ikj_float(int dim, float* matA, float* matB, int verbose);

DoubleResult jki_double(int dim, double* matA, double* matB, int verbose);
FloatResult jki_float(int dim, float* matA, float* matB, int verbose);

DoubleResult kji_double(int dim, double* matA, double* matB, int verbose);
FloatResult kji_float(int dim, float* matA, float* matB, int verbose);

//Tests

DoubleResult *test_double_matrices(int dim, double* matA, double* matB, int trial_id, DoubleResult *res_arr, int mult_algs);
FloatResult *test_float_matrices(int dim, float* matA, float* matB, int trial_id, FloatResult *res_arr, int mult_algs);

int write_double_to_csv(int dim, DoubleResult res, FILE* file);
int write_float_to_csv(int dim, FloatResult res, FILE* file);
#endif
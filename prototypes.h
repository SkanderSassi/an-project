#ifndef PROTOTYPES_H
#define PROTOTYPES_H
#include <stdio.h>
#include "datatypes.h"

// Generators
double generate_double(int low, int high);
float generate_float(int low, int high);

    double *generate_matrix_double(int dim, int low, int high, int verbose);
float* generate_matrix_float(int dim, int low, int high, int verbose);
// Printers
void print_matrix_float(int dim, float* matrix);
void print_matrix_double(int dim, double* matrix);

void print_float_result(int dim, FloatResult res_struct, int show_matrix, int show_time, int show_trial_id, int show_alg);
void print_double_result(int dim, DoubleResult res_struct, int show_matrix, int show_time, int show_trial_id, int show_alg);

// Sanity check
int check_equal_matrix_double(int dim, double* matA, double* matB);
int check_equal_matrix_float(int dim, float* matA, float* matB);

//Addition

DoubleResult ij_add_double(int dim, double* matA, double* matB, int verbose);
FloatResult ij_add_float(int dim, float* matA, float* matB, int verbose);

DoubleResult ji_add_double(int dim, double* matA, double* matB, int verbose);
FloatResult ji_add_float(int dim, float* matA, float* matB, int verbose);



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



int write_double_to_csv(int dim, DoubleResult res, FILE* file);
int write_float_to_csv(int dim, FloatResult res, FILE* file);


#endif
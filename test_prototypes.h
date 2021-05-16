#ifndef TEST_PROTOTYPES_H
#define TEST_PROTOTYPES_H
#include <stdio.h>
#include "datatypes.h"

DoubleResult *test_double_matrices(int dim, double* matA, double* matB, int trial_id, DoubleResult *res_arr, int mult_algs);
FloatResult *test_float_matrices(int dim, float* matA, float* matB, int trial_id, FloatResult *res_arr, int mult_algs);

void test_single_number(FILE *file, int low, int high,int trial_id);

#endif
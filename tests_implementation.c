#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "prototypes.h"
#include "test_prototypes.h"



void test_single_number(FILE *file, int low, int high, int trial_id)
{

    double dn1, dn2, double_result, time_spent;
    float fn1, fn2, float_result;
    clock_t begin, end;

    dn1 = generate_double(low, high);
    dn2 = generate_double(low, high);

    fn1 = generate_float(low, high);
    fn2 = generate_float(low, high);

    // Adding double time test
    begin = clock();
    double_result = dn1 + dn2;
    end = clock();
    // data_type, operation, time_spent, trial_id
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    fprintf(file, "double, ADD, %f, %d\n", time_spent, trial_id);

    // Multiplying double test

    begin = clock();
    double_result = dn1 * dn2;
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    fprintf(file, "double, MULT, %f, %d\n", time_spent, trial_id);

    // Adding float time test

    begin = clock();
    double_result = dn1 + dn2;
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    fprintf(file, "float, ADD, %f, %d\n", time_spent, trial_id);

    // Adding float time test

    begin = clock();
    double_result = dn1 * dn2;
    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    fprintf(file, "float, MULT, %f, %d\n", time_spent, trial_id);
}

DoubleResult *test_double_matrices(int dim, double *matA, double *matB, int trial_id, DoubleResult *res_arr, int mult_algs)
{

    int alg_id = 0;

    //Addition

    DoubleResult matCijk = ijk_double(dim, matA, matB, 0);
    DoubleResult matCjik = jik_double(dim, matA, matB, 0);
    DoubleResult matCkij = kij_double(dim, matA, matB, 0);
    DoubleResult matCikj = ikj_double(dim, matA, matB, 0);
    DoubleResult matCjki = jki_double(dim, matA, matB, 0);
    DoubleResult matCkji = kji_double(dim, matA, matB, 0);

    //Multiplication

    DoubleResult matCij = ij_add_double(dim, matA, matB, 0);
    DoubleResult matCji = ji_add_double(dim, matA, matB, 0);

    matCijk.trial_id = trial_id;
    matCjik.trial_id = trial_id;
    matCkij.trial_id = trial_id;
    matCikj.trial_id = trial_id;
    matCjki.trial_id = trial_id;
    matCkji.trial_id = trial_id;

    matCij.trial_id = trial_id;
    matCji.trial_id = trial_id;

    *(res_arr) = matCijk;
    *(res_arr + 1) = matCjik;
    *(res_arr + 2) = matCkij;
    *(res_arr + 3) = matCikj;
    *(res_arr + 4) = matCjki;
    *(res_arr + 5) = matCkji;

    // Addition

    *(res_arr + 6) = matCij;
    *(res_arr + 7) = matCji;

    return res_arr;
}
FloatResult *test_float_matrices(int dim, float *matA, float *matB, int trial_id, FloatResult *res_arr, int mult_algs)
{

    int alg_id = 0;

    //Multiplication

    FloatResult matCijk = ijk_float(dim, matA, matB, 0);
    FloatResult matCjik = jik_float(dim, matA, matB, 0);
    FloatResult matCkij = kij_float(dim, matA, matB, 0);
    FloatResult matCikj = ikj_float(dim, matA, matB, 0);
    FloatResult matCjki = jki_float(dim, matA, matB, 0);
    FloatResult matCkji = kji_float(dim, matA, matB, 0);

    // Addition

    FloatResult matCij = ij_add_float(dim, matA, matB, 0);
    FloatResult matCji = ji_add_float(dim, matA, matB, 0);

    matCijk.trial_id = trial_id;
    matCjik.trial_id = trial_id;
    matCkij.trial_id = trial_id;
    matCikj.trial_id = trial_id;
    matCjki.trial_id = trial_id;
    matCkji.trial_id = trial_id;

    matCij.trial_id = trial_id;
    matCji.trial_id = trial_id;

    //Multiplication

    *(res_arr) = matCijk;
    *(res_arr + 1) = matCjik;
    *(res_arr + 2) = matCkij;
    *(res_arr + 3) = matCikj;
    *(res_arr + 4) = matCjki;
    *(res_arr + 5) = matCkji;

    // Addition

    *(res_arr + 6) = matCij;
    *(res_arr + 7) = matCji;

    return res_arr;
}

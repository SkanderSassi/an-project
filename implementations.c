#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "prototypes.h"

void print_matrix_float(int dim, float *matrix)
{

    int i, j;
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            printf("%f |", matrix[i * dim + j]);
        }
        printf("\n");
    }
}

void print_matrix_double(int dim, double *matrix)
{

    int i, j;
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            printf("%f |", matrix[i * dim + j]);
        }
        printf("\n");
    }
}

void print_float_result(int dim, FloatResult res_struct, int show_matrix, int show_time, int show_trial_id, int show_alg)
{

    printf(" Data type: float ");
    if (show_matrix)
        print_matrix_float(dim, res_struct.matrix);
    if (show_alg)
        printf(" Algorithm: %s", res_struct.alg);
    if (show_time)
        printf(" Time spent : %f ", res_struct.time_spent);
    if (show_trial_id)
        printf(" Trial id : %d ", res_struct.trial_id);
}

void print_double_result(int dim, DoubleResult res_struct, int show_matrix, int show_time, int show_trial_id, int show_alg)
{

    printf(" Data type: double ");
    if (show_matrix)
        print_matrix_double(dim, res_struct.matrix);
    if (show_alg)
        printf(" Algorithm: %s ", res_struct.alg);
    if (show_time)
        printf("Time spent : %f ", res_struct.time_spent);
    if (show_trial_id)
        printf(" Trial id : %d ", res_struct.trial_id);
}

int check_equal_matrix_double(int dim, double *matA, double *matB)
{
    int isEqual = 1;
    int i, j;
    for (i = 0; i < dim && isEqual; i++)
    {
        for (j = 0; j < dim && isEqual; j++)
        {
            if (!(matA[i * dim + j] == matB[i * dim + j]))
                isEqual = 0;
        }
    }
    return isEqual;
}

int check_equal_matrix_float(int dim, float *matA, float *matB)
{
    int isEqual = 1;
    int i, j;
    for (i = 0; i < dim && isEqual; i++)
    {
        for (j = 0; j < dim && isEqual; j++)
        {
            if (!(matA[i * dim + j] == matB[i * dim + j]))
                isEqual = 0;
        }
    }
    return isEqual;
}

int write_double_to_csv(int dim, DoubleResult res, FILE *file)
{
    //dimension, data_type, algorithm, time_spent, trial_id
    int return_val;

    return_val = fprintf(file, "%d, %s, %s, %f, %d\n", dim, "double", res.alg, res.time_spent, res.trial_id);

    return return_val;
}
int write_float_to_csv(int dim, FloatResult res, FILE *file)
{
    //dimension, data_type, algorithm, time_spent, trial_id
    int return_val;

    return_val = fprintf(file, "%d, %s, %s, %f, %d\n", dim, "float", res.alg, res.time_spent, res.trial_id);

    return return_val;
}

double generate_double(int low, int high){
    double d_num, dn;
    
    d_num = (double)rand() / ((double)RAND_MAX + 1);
    dn = (low + d_num * (high - low));

    return dn;
}

float generate_float(int low, int high){
    float f_num, fn;

    f_num = (double)rand() / ((float)RAND_MAX + 1);
    fn = (low + f_num * (high - low));

    return fn;
}

double *generate_matrix_double(int dim, int low, int high, int verbose)
{

    double *matrix = (double *)malloc(dim * dim * sizeof(double));
    int i, j;
    double r_num;
    clock_t begin, end;
    if (verbose) //1 for show time 0 for no
        begin = clock();
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            
            matrix[i * dim + j] = generate_double(low, high);
        }
    }
    if (verbose)
    {
        end = clock();
        double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("Gen-- type: double dim: %d -- %f seconds\n", dim, time_spent);
    }
    return matrix;
}

float *generate_matrix_float(int dim, int low, int high, int verbose)
{
    float *matrix = (float *)malloc(dim * dim * sizeof(float));
    int i, j;
    clock_t begin, end;
    if (verbose) //1 for show time 0 for no
        begin = clock();
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            matrix[i * dim + j] = generate_float(low, high);
        }
    }
    if (verbose)
    {
        end = clock();
        double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        printf("Gen-- type: float -- dim: %d -- %f seconds\n", dim, time_spent);
    }

    return matrix;
}

DoubleResult ij_add_double(int dim, double *matA, double *matB, int verbose){
    DoubleResult res_struct;
    double *result = (double *)malloc(sizeof(double) * dim * dim);
    int i, j;
    clock_t begin, end;
    
    begin = clock();
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            result[i * dim + j] = matA[i * dim + j] + matB[i * dim + j];
            
        }
    }
    end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("IJ -- dim: %d -- type: double -- %f seconds\n", dim, time_spent);

    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "IJ");
    return res_struct;
}

FloatResult ij_add_float(int dim, float *matA, float *matB, int verbose){
    FloatResult res_struct;
    float *result = (float *)malloc(sizeof(float) * dim * dim);
    int i, j, k;
    clock_t begin, end;

    begin = clock();
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            result[i * dim + j] = matA[i * dim + j] + matB[i * dim + j];
        }
    }
    end = clock();
    float time_spent = (float)(end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("IJ -- dim: %d -- type: float -- %f seconds\n", dim, time_spent);
    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "IJ");
    return res_struct;
}

DoubleResult ji_add_double(int dim, double *matA, double *matB, int verbose){
    DoubleResult res_struct;
    double *result = (double *)malloc(sizeof(double) * dim * dim);
    int i, j;
    clock_t begin, end;

    begin = clock();
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            result[j * dim + i] = matA[j * dim + i] + matB[j * dim + i];
        }
    }
    end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("JI -- dim: %d -- type: double -- %f seconds\n", dim, time_spent);

    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "JI");
    return res_struct;
}

FloatResult ji_add_float(int dim, float *matA, float *matB, int verbose){
    FloatResult res_struct;
    float *result = (float *)malloc(sizeof(float) * dim * dim);
    int i, j, k;
    clock_t begin, end;

    begin = clock();
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            result[j * dim + i] = matA[j * dim + i] + matB[j * dim + j];
        }
    }
    end = clock();
    float time_spent = (float)(end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("JI -- dim: %d -- type: float -- %f seconds\n", dim, time_spent);
    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "JI");
    return res_struct;
}

DoubleResult ijk_double(int dim, double *matA, double *matB, int verbose)
{

    DoubleResult res_struct;
    double *result = (double *)malloc(sizeof(double) * dim * dim);
    int i, j, k;
    clock_t begin, end;
    double sum;

    begin = clock();
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            sum = 0;
            for (k = 0; k < dim; k++)
            {
                sum += matA[i * dim + k] * matB[k * dim + j];
            }
            result[i * dim + j] = sum;
        }
    }
    end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("IJK -- dim: %d -- type: double -- %f seconds\n", dim, time_spent);

    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "IJK");
    return res_struct;
}

FloatResult ijk_float(int dim, float *matA, float *matB, int verbose)
{
    FloatResult res_struct;
    float *result = (float *)malloc(sizeof(float) * dim * dim);
    int i, j, k;
    clock_t begin, end;
    float sum;

    begin = clock();
    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            sum = 0;
            for (k = 0; k < dim; k++)
            {
                sum += matA[i * dim + k] * matB[k * dim + j];
            }
            result[i * dim + j] = sum;
        }
    }
    end = clock();
    float time_spent = (float)(end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("IJK -- dim: %d -- type: float -- %f seconds\n", dim, time_spent);
    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "IJK");
    return res_struct;
}

DoubleResult jik_double(int dim, double *matA, double *matB, int verbose)
{
    DoubleResult res_struct;
    double *result = (double *)malloc(sizeof(double) * dim * dim);
    int i, j, k;
    clock_t begin, end;
    double sum;

    begin = clock();
    for (j = 0; j < dim; j++)
    {
        for (i = 0; i < dim; i++)
        {
            sum = 0;
            for (k = 0; k < dim; k++)
            {
                sum += matA[i * dim + k] * matB[k * dim + j];
            }
            result[i * dim + j] = sum;
        }
    }
    end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("JIK -- dim: %d -- type: double -- %f seconds\n", dim, time_spent);

    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "JIK");
    return res_struct;
}

FloatResult jik_float(int dim, float *matA, float *matB, int verbose)
{
    FloatResult res_struct;
    float *result = (float *)malloc(sizeof(float) * dim * dim);
    int i, j, k;
    clock_t begin, end;
    float sum;

    begin = clock();
    for (j = 0; j < dim; j++)
    {
        for (i = 0; i < dim; i++)
        {
            sum = 0;
            for (k = 0; k < dim; k++)
            {
                sum += matA[i * dim + k] * matB[k * dim + j];
            }
            result[i * dim + j] = sum;
        }
    }
    end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("JIK -- dim: %d -- type: float -- %f seconds\n", dim, time_spent);

    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "JIK");

    return res_struct;
}

DoubleResult kij_double(int dim, double *matA, double *matB, int verbose)
{
    DoubleResult res_struct;
    double *result = (double *)malloc(sizeof(double) * dim * dim);
    int i, j, k;
    clock_t begin, end;
    double r;

    begin = clock();

    for (k = 0; k < dim; k++)
    {
        for (i = 0; i < dim; i++)
        {
            r = matA[i * dim + k];

            for (j = 0; j < dim; j++)
            {
                result[i * dim + j] += r * matB[k * dim + j];
            }
        }
    }

    end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("KIJ -- dim: %d -- type: double -- %f seconds\n", dim, time_spent);
    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "KIJ");

    return res_struct;
}

FloatResult kij_float(int dim, float *matA, float *matB, int verbose)
{

    FloatResult res_struct;
    float *result = (float *)malloc(sizeof(float) * dim * dim);
    int i, j, k;
    clock_t begin, end;

    float r;

    begin = clock();

    for (k = 0; k < dim; k++)
    {
        for (i = 0; i < dim; i++)
        {
            r = matA[i * dim + k];

            for (j = 0; j < dim; j++)
            {
                result[i * dim + j] += r * matB[k * dim + j];
            }
        }
    }

    end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    if (verbose)
        printf("KIJ -- dim: %d -- type: float -- %f seconds\n", dim, time_spent);
    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "KIJ");

    return res_struct;
}

DoubleResult ikj_double(int dim, double *matA, double *matB, int verbose)
{
    DoubleResult res_struct;
    double *result = (double *)malloc(sizeof(double) * dim * dim);
    int i, j, k;
    clock_t begin, end;
    double r;

    begin = clock();
    for (i = 0; i < dim; i++)
    {
        for (k = 0; k < dim; k++)
        {
            r = matA[i * dim + k];
            for (j = 0; j < dim; j++)
            {
                result[i * dim + j] += r * matB[k * dim + j];
            }
        }
    }
    end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("IKJ -- dim: %d -- type: double -- %f seconds\n", dim, time_spent);

    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "IKJ");

    return res_struct;
}

FloatResult ikj_float(int dim, float *matA, float *matB, int verbose)
{
    FloatResult res_struct;
    float *result = (float *)malloc(sizeof(float) * dim * dim);
    int i, j, k;
    clock_t begin, end;

    float r;

    begin = clock();

    for (i = 0; i < dim; i++)
    {
        for (k = 0; k < dim; k++)
        {
            r = matA[i * dim + k];
            for (j = 0; j < dim; j++)
            {
                result[i * dim + j] += r * matB[k * dim + j];
            }
        }
    }
    end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    if (verbose)
        printf("IKJ -- dim: %d -- type: float -- %f seconds\n", dim, time_spent);

    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "KIJ");

    return res_struct;
}

DoubleResult jki_double(int dim, double *matA, double *matB, int verbose)
{
    DoubleResult res_struct;
    double *result = (double *)malloc(sizeof(double) * dim * dim);
    int i, j, k;
    clock_t begin, end;
    double r;

    begin = clock();
    for (j = 0; j < dim; j++)
    {

        for (k = 0; k < dim; k++)
        {
            r = matB[k * dim + j];
            for (i = 0; i < dim; i++)
            {
                result[i * dim + j] += r * matA[i * dim + k];
            }
        }
    }
    end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("JKI -- dim: %d -- type: double -- %f seconds\n", dim, time_spent);

    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "JKI");

    return res_struct;
}

FloatResult jki_float(int dim, float *matA, float *matB, int verbose)
{
    FloatResult res_struct;
    float *result = (float *)malloc(sizeof(float) * dim * dim);
    int i, j, k;
    clock_t begin, end;

    float r;

    begin = clock();

    for (j = 0; j < dim; j++)
    {
        for (k = 0; k < dim; k++)
        {
            r = matB[k * dim + j];
            for (i = 0; i < dim; i++)
            {
                result[i * dim + j] += r * matA[i * dim + k];
            }
        }
    }
    end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("JKI -- dim: %d -- type: float -- %f seconds\n", dim, time_spent);

    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "JKI");

    return res_struct;
}

DoubleResult kji_double(int dim, double *matA, double *matB, int verbose)
{

    DoubleResult res_struct;

    double *result = (double *)malloc(sizeof(double) * dim * dim);
    int i, j, k;
    clock_t begin, end;
    double r;

    begin = clock();
    for (k = 0; k < dim; k++)
    {

        for (j = 0; j < dim; j++)
        {
            r = matB[k * dim + j];
            for (i = 0; i < dim; i++)
            {
                result[i * dim + j] += r * matA[i * dim + k];
            }
        }
    }

    end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("KJI -- dim: %d -- type: double -- %f seconds\n", dim, time_spent);

    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "KJI");

    return res_struct;
}

FloatResult kji_float(int dim, float *matA, float *matB, int verbose)
{
    FloatResult res_struct;
    float *result = (float *)malloc(sizeof(float) * dim * dim);
    int i, j, k;
    clock_t begin, end;

    float r;

    begin = clock();

    for (k = 0; k < dim; k++)
    {
        for (j = 0; j < dim; j++)
        {
            r = matB[k * dim + j];
            for (i = 0; i < dim; i++)
            {
                result[i * dim + j] += r * matA[i * dim + k];
            }
        }
    }

    end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("KJI -- dim: %d -- type: float -- %f seconds\n", dim, time_spent);

    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "KJI");

    return res_struct;
}



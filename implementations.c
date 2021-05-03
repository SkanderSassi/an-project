#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "prototypes.h"

void print_matrix_float(int dim, float* matrix){

    int i,j;
    for(i=0; i<dim; i++){
        for(j=0; j<dim; j++){
            printf("%f |",matrix[i*dim+j]);
        }
        printf("\n");
    }
}

void print_matrix_double(int dim, double* matrix){

    int i,j;
    for(i=0; i<dim; i++){
        for(j=0; j<dim; j++){
            printf("%f |",matrix[i*dim+j]);
        }
        printf("\n");
    }
}

void print_float_result(int dim, FloatResult res_struct, int show_matrix, int show_time, int show_trial_id, int show_alg){
    
    printf(" Data type: float ");
    if(show_matrix)
        print_matrix_float(dim, res_struct.matrix);
    if(show_alg)
        printf(" Algorithm: %s", res_struct.alg);
    if(show_time)
        printf(" Time spent : %f ", res_struct.time_spent);
    if(show_trial_id)
        printf(" Trial id : %d ", res_struct.trial_id);
    

}

void print_double_result(int dim, DoubleResult res_struct, int show_matrix, int show_time, int show_trial_id, int show_alg){
    
    printf(" Data type: double ");
    if(show_matrix)
        print_matrix_double(dim, res_struct.matrix);
    if(show_alg)
        printf(" Algorithm: %s ", res_struct.alg);
    if(show_time)
        printf( "Time spent : %f ", res_struct.time_spent);
    if(show_trial_id)
        printf(" Trial id : %d ", res_struct.trial_id);
}

int check_equal_matrix_double(int dim, double* matA, double* matB){
    int isEqual = 1;
    int i,j;
    for(i=0; i<dim && isEqual; i++){
        for(j=0; j<dim && isEqual; j++){
            if(!(matA[i*dim+j] == matB[i*dim+j]))
                isEqual = 0;
        }
    }
    return isEqual;
}


int check_equal_matrix_float(int dim, float* matA, float* matB){
    int isEqual = 1;
    int i,j;
    for(i=0; i<dim && isEqual; i++){
        for(j=0; j<dim && isEqual; j++){
            if(!(matA[i*dim+j] == matB[i*dim+j]))
                isEqual = 0;
        }
    }
    return isEqual;
}

double* generate_matrix_double(int dim, int low, int high, int verbose){
    
    double* matrix = (double *) malloc(dim * dim * sizeof(double));
    int i,j;
    double r_num;
    clock_t begin, end;
    if(verbose) //1 for show time 0 for no
        begin = clock();
    for (i=0; i<dim; i++){
        for(j=0; j<dim; j++){
            r_num = (double) rand() / ((double) RAND_MAX + 1);
            matrix[i*dim+j] = (low + r_num * (high - low));
        }
    }
    if (verbose){
        end = clock();
        double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
        printf("Gen-- type: double dim: %d -- %f seconds\n", dim, time_spent);
    }
    return matrix;
}


float* generate_matrix_float(int dim, int low, int high, int verbose){
    float* matrix = (float *) malloc(dim * dim * sizeof(float));
    int i,j;
    float r_num;
    clock_t begin, end;
    if(verbose) //1 for show time 0 for no
        begin = clock();
    for (i=0; i<dim; i++){
        for(j=0; j<dim; j++){
            r_num = (float) rand() / ((float) RAND_MAX + 1);
            matrix[i*dim+j] = (low + r_num * (high - low));
        }
    }
    if (verbose){
        end = clock();
        double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
        printf("Gen-- type: float -- dim: %d -- %f seconds\n", dim, time_spent);
    }
    
    return matrix;
}


DoubleResult ijk_double(int dim, double* matA, double* matB, int verbose){

    DoubleResult res_struct;
    double* result = (double*) malloc(sizeof(double) * dim * dim);
    int i,j,k;
    clock_t begin, end;
    double sum;
    if(verbose)
        begin = clock();
    for(i=0; i<dim; i++){
        for(j=0; j<dim; j++){
            sum =0;
            for(k=0; k< dim; k++){
                sum+= matA[i*dim + k] * matB[k*dim + j]; 
            }
            result[i*dim + j] = sum;
        }
    }
    end = clock();
    double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("IJK -- dim: %d -- type: double -- %f seconds\n", dim, time_spent);

    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "IJK");
    return res_struct;
}

FloatResult ijk_float(int dim, float* matA, float* matB, int verbose){
    FloatResult res_struct;
    float* result = (float*) malloc(sizeof(float) * dim * dim);
    int i,j,k;
    clock_t begin, end;
    float sum;
    if(verbose)
        begin = clock();
    for(i=0; i<dim; i++){
        for(j=0; j<dim; j++){
            sum =0;
            for(k=0; k< dim; k++){
                sum+= matA[i*dim + k] * matB[k*dim + j]; 
            }
            result[i*dim + j] = sum;
        }
    }
    end = clock();
    float time_spent = (float) (end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("IJK -- dim: %d -- type: float -- %f seconds\n", dim, time_spent);
    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "IJK");
    return res_struct;
}

DoubleResult jik_double(int dim, double* matA, double* matB, int verbose){
    DoubleResult res_struct;
    double* result = (double*) malloc(sizeof(double) * dim * dim);
    int i,j,k;
    clock_t begin, end;
    double sum;
    if(verbose)
        begin = clock();
    for(j=0; j<dim; j++){
        for(i=0; i<dim; i++){
            sum =0;
            for(k=0; k< dim; k++){
                sum+= matA[i*dim + k] * matB[k*dim + j]; 
            }
            result[i*dim + j] = sum;
        }
    }
    end = clock();
    double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("JIK -- dim: %d -- type: double -- %f seconds\n", dim, time_spent);

    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "JIK");
    return res_struct;
}

FloatResult jik_float(int dim, float* matA, float* matB, int verbose){
    FloatResult res_struct;
    float* result = (float*) malloc(sizeof(float) * dim * dim);
    int i,j,k;
    clock_t begin, end;
    float sum;
    if(verbose)
        begin = clock();
    for(j=0; j<dim; j++){
        for(i=0; i<dim; i++){
            sum =0;
            for(k=0; k< dim; k++){
                sum+= matA[i*dim + k] * matB[k*dim + j]; 
            }
            result[i*dim + j] = sum;
        }
    }
    end = clock();
    double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("JIK -- dim: %d -- type: float -- %f seconds\n", dim, time_spent);

    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "JIK");

    return res_struct;
}

DoubleResult kij_double(int dim, double* matA, double* matB, int verbose){
    DoubleResult res_struct;
    double* result = (double*) malloc(sizeof(double) * dim * dim);
    int i,j,k;
    clock_t begin, end;
    double r;

    if(verbose)
        begin = clock();

    for (k=0; k<dim; k++){
        for(i=0; i<dim; i++){
            r = matA[i*dim + k];
            
            for(j=0; j<dim; j++ ){
                result[i*dim + j] += r * matB[k*dim + j];
            }
        }
    }

    end = clock();
    double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("KIJ -- dim: %d -- type: double -- %f seconds\n", dim, time_spent);
    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "KIJ");

    return res_struct;
}

FloatResult kij_float(int dim, float* matA, float* matB, int verbose){
    
    FloatResult res_struct;
    float* result = (float*) malloc(sizeof(float) * dim * dim);
    int i,j,k;
    clock_t begin, end;

    float r;

    if(verbose)
        begin = clock();

    for (k=0; k<dim; k++){
        for(i=0; i<dim; i++){
            r = matA[i*dim + k];
            
            for(j=0; j<dim; j++ ){
                result[i*dim + j] += r * matB[k*dim + j];
            }
        }
    }

    end = clock();
    double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;

    if (verbose)
        printf("KIJ -- dim: %d -- type: float -- %f seconds\n", dim, time_spent);
    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "KIJ");


    return res_struct;
}

DoubleResult ikj_double(int dim, double* matA, double* matB, int verbose){
    DoubleResult res_struct;
    double* result = (double*) malloc(sizeof(double) * dim * dim);
    int i,j,k;
    clock_t begin, end;
    double r;

    if(verbose)
        begin = clock();
    for (i=0; i<dim;i++){
        for(k=0; k<dim; k++){
            r = matA[i*dim + k];
            for(j=0; j<dim; j++ ){
                result[i*dim + j] += r * matB[k*dim + j];
            }
        }
    }
    end = clock();
    double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
    if (verbose)    
        printf("IKJ -- dim: %d -- type: double -- %f seconds\n", dim, time_spent);
    
    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "IKJ");

    return res_struct;
}

FloatResult ikj_float(int dim, float* matA, float* matB, int verbose){
    FloatResult res_struct;
    float* result = (float*) malloc(sizeof(float) * dim * dim);
    int i,j,k;
    clock_t begin, end;

    float r;

    if(verbose)
        begin = clock();

    for (i=0; i<dim;i++){
        for(k=0; k<dim; k++){
            r = matA[i*dim + k];
            for(j=0; j<dim; j++ ){
                result[i*dim + j] += r * matB[k*dim + j];
            }
        }
    }
    end = clock();
    double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
    
    if (verbose)
        printf("IKJ -- dim: %d -- type: float -- %f seconds\n", dim, time_spent);

    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "KIJ");

    
    return res_struct;
}

DoubleResult jki_double(int dim, double* matA, double* matB, int verbose){
    DoubleResult res_struct;
    double* result = (double*) malloc(sizeof(double) * dim * dim);
    int i,j,k;
    clock_t begin, end;
    double r;

    if(verbose)
        begin = clock();
    for (j=0; j<dim;j++){
        
        for(k=0; k<dim; k++){
            r = matB[k*dim + j];
            for(i=0; i<dim; i++){
                result[i*dim + j] += r * matA[i*dim + k];
            }
        }
    }
    end = clock();
    double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("JKI -- dim: %d -- type: double -- %f seconds\n", dim, time_spent);
    
    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "JKI");

    return res_struct;
}

FloatResult jki_float(int dim, float* matA, float* matB, int verbose){
    FloatResult res_struct;
    float* result = (float*) malloc(sizeof(float) * dim * dim);
    int i,j,k;
    clock_t begin, end;

    float r;

    
    begin = clock();

    for (j=0; j<dim;j++){
        for(k=0; k<dim; k++){
            r = matB[k*dim + j];
            for(i=0; i<dim; i++ ){
                result[i*dim + j] += r * matA[i*dim + k];
            }
        }
    }
    end = clock();
    double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("JKI -- dim: %d -- type: float -- %f seconds\n", dim, time_spent);
    
        
    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "JKI");

    

    return res_struct;
}

DoubleResult kji_double(int dim, double* matA, double* matB, int verbose){
    
    DoubleResult res_struct;

    double* result = (double*) malloc(sizeof(double) * dim * dim);
    int i,j,k;
    clock_t begin, end;
    double r;

    if(verbose)
        begin = clock();
    for (k=0; k<dim;k++){
        
        for(j=0; j<dim; j++){
            r = matB[k*dim + j];
            for(i=0; i<dim; i++){
                result[i*dim + j] += r * matA[i*dim + k];
            }
        }
    }

    end = clock();
    double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("KJI -- dim: %d -- type: double -- %f seconds\n", dim, time_spent);
    
    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "KJI");

    

    return res_struct;
}

FloatResult kji_float(int dim, float* matA, float* matB, int verbose){
    FloatResult res_struct;
    float* result = (float*) malloc(sizeof(float) * dim * dim);
    int i,j,k;
    clock_t begin, end;

    float r;

    
    begin = clock();

    for (k=0; k<dim;k++){
        for(j=0; j<dim; j++){
            r = matB[k*dim + j];
            for(i=0; i<dim; i++ ){
                result[i*dim + j] += r * matA[i*dim + k];
            }
        }
    }

    end = clock();
    double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
    if (verbose)
        printf("KJI -- dim: %d -- type: float -- %f seconds\n", dim, time_spent);
    
    res_struct.matrix = result;
    res_struct.time_spent = time_spent;
    strcpy(res_struct.alg, "KJI");

    

    return res_struct;
}

DoubleResult *test_double_matrices(int dim, double* matA, double* matB,int trial_id,  DoubleResult *res_arr, int mult_algs){
    
    int alg_id;

    DoubleResult matCijk = ijk_double(dim, matA, matB, 0);
    DoubleResult matCjik = jik_double(dim, matA, matB, 0);
    DoubleResult matCkij = kij_double(dim, matA, matB, 0);
    DoubleResult matCikj = ikj_double(dim, matA, matB, 0);
    DoubleResult matCjki = jki_double(dim, matA, matB, 0);
    DoubleResult matCkji = kji_double(dim, matA, matB, 0);

    *(res_arr) = matCijk;
    *(res_arr+1) = matCjik;
    *(res_arr+2) = matCkij;
    *(res_arr+3) = matCikj;
    *(res_arr+4) = matCjki;
    *(res_arr+5) = matCkji;

    for(alg_id=0; alg_id < mult_algs; alg_id++)
        (res_arr+alg_id)->trial_id;


    return res_arr;
}
FloatResult *test_float_matrices(int dim, float* matA, float* matB, int trial_id , FloatResult *res_arr, int mult_algs){
    
    int alg_id;

    FloatResult matCijk = ijk_float(dim, matA, matB, 0);
    FloatResult matCjik = jik_float(dim, matA, matB, 0);
    FloatResult matCkij = kij_float(dim, matA, matB, 0);
    FloatResult matCikj = ikj_float(dim, matA, matB, 0);
    FloatResult matCjki = jki_float(dim, matA, matB, 0);
    FloatResult matCkji = kji_float(dim, matA, matB, 0);

    
    *(res_arr) = matCijk;
    *(res_arr+1) = matCjik;
    *(res_arr+2) = matCkij;
    *(res_arr+3) = matCikj;
    *(res_arr+4) = matCjki;
    *(res_arr+5) = matCkji;

    for(alg_id=0; alg_id < mult_algs; alg_id++)
        (res_arr+alg_id)->trial_id;


    return res_arr;
}



#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "prototypes.h"

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

double* ijk_double(int dim, double* matA, double* matB, int verbose){
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
    if (verbose){
        end = clock();
        double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
        printf("IJK -- dim: %d -- type: double -- %f seconds\n", dim, time_spent);
    }
    return result;
}

float* ijk_float(int dim, float* matA, float* matB, int verbose){
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
    if (verbose){
        end = clock();
        float time_spent = (float) (end - begin) / CLOCKS_PER_SEC;
        printf("IJK -- dim: %d -- type: float -- %f seconds\n", dim, time_spent);
    }
    return result;
}

double* jik_double(int dim, double* matA, double* matB, int verbose){
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
    if (verbose){
        end = clock();
        double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
        printf("JIK -- dim: %d -- type: double -- %f seconds\n", dim, time_spent);
    }
    return result;
}

float* jik_float(int dim, float* matA, float* matB, int verbose){
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
    if (verbose){
        end = clock();
        double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
        printf("JIK -- dim: %d -- type: float -- %f seconds\n", dim, time_spent);
    }
    return result;
}

void test_double_matrices(int dim, double* matA, double* matB){
    
    double* matCijk = (double*) malloc(dim*dim*sizeof(double));
    double* matCjik = (double*) malloc(dim*dim*sizeof(double));
    
    matCijk = ijk_double(dim, matA, matB, 1);
    matCjik = jik_double(dim, matA, matB, 1);


}
void test_float_matrices(int dim, float* matA, float* matB){
    
    float* matCijk = (float*) malloc(dim*dim*sizeof(float));
    float* matCjik = (float*) malloc(dim*dim*sizeof(float));
    
    matCijk = ijk_float(dim, matA, matB, 1);
    matCjik = jik_float(dim, matA, matB, 1);
}



#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#define HIGH 20
#define LOW 10
#define TRIES 5
#define MULT_ALGS 6 



int main(int argc, char* argv[]){

    if (argc < 2){
        printf("Please input matrix dimension\n");
        return EXIT_FAILURE;
    }

    FloatResult *fl_res_arr = malloc(MULT_ALGS*sizeof(FloatResult));
    DoubleResult *dbl_res_arr = malloc(MULT_ALGS*sizeof(DoubleResult));

    int dim, mult_id;
    double *matAd, *matBd, *matCd;
    float *matAf, *matBf, *matCf;
    srand(1);
    int index;
    int trial_id;
    for(index=1; index<argc; index++){ // For every dimensions to test on 
        
        dim = strtol(argv[index], NULL, 10);

        matAd = generate_matrix_double(dim, LOW, HIGH, 0);
        matBd = generate_matrix_double(dim, LOW, HIGH, 0);
        
        matAf = generate_matrix_float(dim, LOW, HIGH, 0);
        matBf = generate_matrix_float(dim, LOW, HIGH, 0);
        printf("Dimension : %d\n", dim);
        for(trial_id=0; trial_id<TRIES; trial_id++){ 

            dbl_res_arr = test_double_matrices(dim, matAd, matBd, trial_id, dbl_res_arr, MULT_ALGS);
            fl_res_arr = test_float_matrices(dim, matAf, matBf, trial_id, fl_res_arr, MULT_ALGS);

            for(mult_id=0; mult_id<MULT_ALGS; mult_id++){
                print_double_result(dim, *(dbl_res_arr + mult_id), 0, 1, 1, 1);
                printf("\n");
                print_float_result(dim, *(fl_res_arr + mult_id), 0, 1, 1, 1);
                printf("\n");
            }


        }
    }

    
    
    
    
    

    return EXIT_SUCCESS;
}
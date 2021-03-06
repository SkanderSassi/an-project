#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#include "test_prototypes.h"

#define HIGH 20
#define LOW 10
#define TRIES 5
#define ALGS 8 



int main(int argc, char* argv[]){

    if (argc < 2){
        printf("Please input matrix dimension\n");
        return EXIT_FAILURE;
    }
    FloatResult *fl_res_arr = malloc(ALGS*sizeof(FloatResult));
    DoubleResult *dbl_res_arr = malloc(ALGS*sizeof(DoubleResult));

    int dim, mult_id, index, trial_id, printed_vals;
    double *matAd, *matBd, *matCd;
    float *matAf, *matBf, *matCf;
    FILE *csv_file, *csv_file2;

    // Matrices
    csv_file = fopen("data.csv", "w+");
    // Single numbers

    // csv_file2 = fopen("single_num.csv", "w+");

    if (csv_file == NULL){
        puts("Unable to open file");
        return EXIT_FAILURE;
    }
    //Set up columns

    fprintf(csv_file, "dimension,data_type,algorithm,time_spent,trial_id\n");
    // fprintf(csv_file2, "data_type, operation, time_spent, trial_id\n");
    // Random Seed
    srand(1);
    
    for(index=1; index<argc; index++){ // For every dimensions to test on 
        
        dim = strtol(argv[index], NULL, 10);

        matAd = generate_matrix_double(dim, LOW, HIGH, 0);
        matBd = generate_matrix_double(dim, LOW, HIGH, 0);

        matAf = generate_matrix_float(dim, LOW, HIGH, 0);
        matBf = generate_matrix_float(dim, LOW, HIGH, 0);


        printf("----------------------------- Dimension : %d ----------------------------\n", dim);
        for(trial_id=0; trial_id<TRIES; trial_id++){

            // test_single_number(csv_file2, LOW, HIGH, trial_id);
            dbl_res_arr = test_double_matrices(dim, matAd, matBd, trial_id, dbl_res_arr, ALGS);
            fl_res_arr = test_float_matrices(dim, matAf, matBf, trial_id, fl_res_arr, ALGS);

            for(mult_id=0; mult_id<ALGS; mult_id++){
                print_double_result(dim, *(dbl_res_arr + mult_id), 0, 1, 1, 1); 
                
                printed_vals = write_double_to_csv(dim, *(dbl_res_arr + mult_id), csv_file);
                printf("\n");
            }
            for(mult_id=0; mult_id<ALGS; mult_id++){
                print_float_result(dim, *(fl_res_arr + mult_id), 0, 1, 1, 1);
                printed_vals = write_float_to_csv(dim, *(fl_res_arr + mult_id), csv_file);
                printf("\n");
            }
        }
    }

    fclose(csv_file);
    // fclose(csv_file2);
    
    
    

    return EXIT_SUCCESS;
}
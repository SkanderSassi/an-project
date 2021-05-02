#include <stdio.h>
#include <stdlib.h>
#include "prototypes.h"
#define HIGH 20
#define LOW 10

int main(int argc, char* argv[]){

    if (argc < 2){
        printf("Please input matrix dimension\n");
        return EXIT_FAILURE;
    }

    int dim = strtol(argv[1], NULL, 10);
    double *matAd, *matBd, *matCd;
    float *matAf, *matBf, *matCf;
    srand(1);
    
    matAd = generate_matrix_double(dim, LOW, HIGH, 0);
    matBd = generate_matrix_double(dim, LOW, HIGH, 0);
    
    matAf = generate_matrix_float(dim, LOW, HIGH, 0);
    matBf = generate_matrix_float(dim, LOW, HIGH, 0);

    test_double_matrices(dim, matAd, matBd);
    test_float_matrices(dim, matAf, matBf);


    
    
    
    
    

    return EXIT_SUCCESS;
}
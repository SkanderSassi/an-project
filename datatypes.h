#ifndef DATATYPES_H
#define DATATYPES_H

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

#endif
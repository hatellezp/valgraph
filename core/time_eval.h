//
// Created by horacio on 29/06/2020.
//

#include <time.h>
#include <stdlib.h>
#include <math.h>

#include "matrix.h"
#include "algrank.h"

#ifndef VALGRAPHCORE_TIME_EVAL_H
#define VALGRAPHCORE_TIME_EVAL_H

double * gen_rand_array(int rows, int cols, double limit);

double single_execution_time(int matrix_size, double limit, int method);
double * execution_time_stats(int matrix_size, double limit, int reps, int method);
double * execution_time_array(int beg, int end, double limit, int reps, int method);

#endif //VALGRAPHCORE_TIME_EVAL_H

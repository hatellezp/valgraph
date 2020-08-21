//
// Created by horacio on 29/06/2020.
//

#include "time_eval.h"

double * gen_rand_array(int rows, int cols, double limit) {
    // force a new random seed each time
    srand(time(NULL));

    int n = rows * cols;
    double * m = malloc(sizeof(double )*n);

    for(int i=0; i<n; i++) {
        double x = ((double)rand()/(double)(RAND_MAX)) * limit;
        m[i] = x;
    }
    return m;
}

// measuring time functions
double single_execution_time(int matrix_size, double limit, int method) {
    clock_t t;
    t = clock();

    int N = matrix_size;
    double * array = gen_rand_array(N, N, limit);

    gsl_matrix_complex * m = from_array_dynamic(N, N, array);
    double tolerance = 0.0001; // it's not important
    bound_only_computation(m, tolerance, method); // I don't need to retrieve the value

    t = clock() - t;
    double time_taken = ((double) t) / CLOCKS_PER_SEC;

    return time_taken;
}

double * execution_time_stats(int matrix_size, double limit, int reps, int method) {
    double * times = malloc(sizeof(double)*3);

    double current_min = RAND_MAX;
    double current_max = -1;
    double average = 0;
    double current_time;

    // loop 'reps' times
    for(int i=0; i<reps; i++) {
        current_time = single_execution_time(matrix_size, limit, method);
        current_max = fmax(current_max, current_time);
        current_min = fmin(current_min, current_time);

        if (i == 0) {
            average = current_time;
        } else {
            double id = (double)i;
            average = (average*id + current_time) / (id + 1.);
        }
    }

    times[0] = current_min;
    times[1] = average;
    times[2] = current_max;

    return times;
}

double * execution_time_array(int beg, int end, double limit, int reps, int method) {
    int n = end - beg + 1;
    double * stats = malloc(sizeof(double)*n*3);
    double * in_stats = malloc(sizeof(double)*3);

    for(int i=0; i<n; i++) {
        in_stats = execution_time_stats(i+beg, limit, reps, method);

        stats[i*3] = in_stats[0];
        stats[i*3 + 1] = in_stats[1];
        stats[i*3 + 2] = in_stats[2];
    }

    free(in_stats);
    return stats;
}
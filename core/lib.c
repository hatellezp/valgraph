//
// Created by horacio on 20/08/2020.
//

#include "lib.h"

// the matrix must be square
gsl_matrix * create_matrix_dyn(size_t n, double * dyn_array) {
    gsl_matrix  * m = gsl_matrix_alloc(n , n);
    size_t i,j;
    for(i=0; i<n; i++) {
        for(j=0; i<n; i++) {
            gsl_matrix_set(m, i, j, dyn_array[i*n + j]);
        }
    }

    return m;
}

gsl_matrix * create_matrix_stat(size_t n, double stat_array[n][n]) {
    gsl_matrix * m = gsl_matrix_alloc(n, n);
    size_t i,j;

    for(i=0; i<n; i++) {
        for(j=0; j<n; j++) {
            gsl_matrix_set(m, i, j, stat_array[i][j]);
        }
    }

    return m;
}

gsl_matrix_complex * create_matrix_complex_dyn(double * dyn_array, size_t n) {
    return from_array_dynamic(n, n, dyn_array);
}

gsl_matrix_complex * create_matrix_complex_stat(double stat_array[], size_t n) {
    return from_array_static(n, n, stat_array);
}

vgraph_result solve(gsl_matrix * m, double a, double b, double c);
vgraph_result is_equivalent(gsl_matrix * m, double a, double b, double a2, double b2);
vgraph_result find_bound_A(gsl_matrix * m, double b);
vgraph_result find_bound_B(gsl_matrix * m, double a);







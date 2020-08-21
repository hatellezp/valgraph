//
// Created by horacio on 20/08/2020.
//



#ifndef VALGRAPHCORE_LIB_H
#define VALGRAPHCORE_LIB_H

// external includes
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

// internal includes
#include "algrank.h"
#include "c_errno.h"
#include "matrix.h"
#include "polynomial.h"
#include "polynomial_vector.h"
#include "polynomial_matrix.h"
#include "polynomial_matrix.h"

enum State {
    Success = 0,
    Uninvertible = 1,
    Unkown = 2,
};

typedef struct {
    enum State state;
    double number;
}vgraph_result;

gsl_matrix * create_matrix_dyn(size_t n, double * dyn_array);
gsl_matrix * create_matrix_stat(size_t n, double stat_array[n][n]);

gsl_matrix_complex * create_matrix_complex_dyn(size_t n, double * dyn_array);
gsl_matrix_complex * create_matrix_complex_stat(size_t n, double stat_array[n][n]);

vgraph_result solve(gsl_matrix * m, double a, double b, double c);
vgraph_result is_equivalent(gsl_matrix * m, double a, double b, double a2, double b2);
vgraph_result find_bound_A(gsl_matrix * m, double b);
vgraph_result find_bound_B(gsl_matrix * m, double a);


#endif //VALGRAPHCORE_LIB_H

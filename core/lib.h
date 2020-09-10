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

#include "test_lib.h"

#define VALGRAPH_TOLERANCE 0.00001

enum State {
    SUCCESS = 0,
    UNINVERTIBLE = 1,
    UNKOWN = 2,
    EQUIVALENT = 3,
    UNEQUIVALENT = 4,
};

typedef struct {
    enum State state;
    double number;
}vgraph_result;

vgraph_result BAD_RESULT = { UNKOWN, -1. };

gsl_matrix * create_matrix_dyn(size_t n, double * dyn_array);
gsl_matrix * create_matrix_stat(size_t n, double stat_array[n][n]);

gsl_matrix_complex * create_matrix_complex_dyn(size_t n, double * dyn_array);
gsl_matrix_complex * create_matrix_complex_stat(size_t n, double stat_array[n][n]);

enum State solve(gsl_matrix * m, gsl_vector * x, double a, double b, double c);
enum State is_equivalent(gsl_matrix * m, double a, double b, double a2, double b2);
vgraph_result find_bound_A(gsl_matrix * m, double b);
vgraph_result find_bound_B(gsl_matrix * m, double a);


#endif //VALGRAPHCORE_LIB_H

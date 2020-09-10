//
// Created by horacio on 29/06/2020.
//

#include <math.h>

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_matrix_complex_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

#include <fftw3.h>

#include "matrix.h"

// the method for computing the bound
#define FULL_V 0
#define SINGLE_V 1

#define UR_FORWARD 0
#define UR_BACKWARDS 1


#ifndef VALGRAPHCORE_ALGRANK_H
#define VALGRAPHCORE_ALGRANK_H

// functions for the bound and ranking calculation
gsl_complex compute_v(gsl_matrix_complex * m, gsl_complex mmod, gsl_complex vmod, size_t ind);
void compute_full_v(gsl_matrix_complex * m, gsl_vector_complex * x, gsl_complex identity_mod, gsl_complex vector_mod);
void create_unity_roots(gsl_vector_complex * units, size_t n, int inverse);
double bound_only_computation(gsl_matrix_complex * m, double tolerance, int method);

// this functions are more front api that
void compute_ranking(gsl_matrix_complex * m, double a, double b, double c, int * error, double * ranking);
void compute_ranking_stable(gsl_matrix_complex * m, double tolerance, int method, int * error, double * ranking);

// some interpolation and utilities algorithms
void polynomial_coefficient_computation(double * x, double * y, double * cof, size_t n);

#endif //VALGRAPHCORE_ALGRANK_H

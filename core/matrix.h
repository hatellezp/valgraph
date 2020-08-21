//
// Created by horacio on 27/06/2020.
//

#ifndef VALGRAPHCORE_MATRIX_H
#define VALGRAPHCORE_MATRIX_H

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_matrix_complex_double.h>

// printers of vectors and matrix, real and complex ones
void print_gsl_matrix(gsl_matrix * m);
void print_gls_vector(gsl_vector * v);

void print_gls_vector_complex(gsl_vector_complex * v);
void print_gsl_matrix_complex(gsl_matrix_complex * m);

void print_array_dyn(int rows, int cols, double * array);

// matrix and vector manipulation, creation and transformation
gsl_matrix * gsl_sq_matrix(size_t sz);
gsl_vector * gsl_vector_ones(int n);
void gsl_matrix_chrow(gsl_matrix * m, gsl_vector * v, int col_index);
void gsl_matrix_chcol(gsl_matrix * m, gsl_vector * v, int row_index);

gsl_matrix_complex * gsl_sq_matrix_complex(size_t sz);
gsl_vector_complex * gsl_vector_complex_ones(int n);
void gsl_matrix_complex_chrow(gsl_matrix_complex * m, gsl_vector_complex * v, int row_index);
void gsl_matrix_complex_chcol(gsl_matrix_complex * m, gsl_vector_complex * v, int col_index);


// create matrix from different types of array
gsl_matrix_complex * from_gsl_matrix(gsl_matrix * m);
gsl_matrix_complex * from_array_dynamic(int rows, int cols, double * array);
gsl_matrix_complex * from_array_static(int rows, int cols, double array[rows][cols]);

#endif //VALGRAPHCORE_MATRIX_H

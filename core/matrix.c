//
// Created by horacio on 27/06/2020.
//

#include "matrix.h"

/*---------------------------------------------------------------------------*/
// to print matrix and vectors, complex or full real ones
void print_gls_vector(gsl_vector * v) {
    int vsize = v -> size;
    for(int i=0; i<vsize; i++) {
        printf("%f ", gsl_vector_get(v, i));
    }
}

void print_gls_vector_complex(gsl_vector_complex * v) {
    int vsize = v -> size;
    for(int i=0; i<vsize; i++) {
        gsl_complex z = gsl_vector_complex_get(v, i);
        double real, imag;
        real = GSL_REAL(z);
        imag = GSL_IMAG(z);
        printf("(%f + %fi) ", real, imag);
    }
}

void print_gsl_matrix(gsl_matrix * m) {

    int r = m -> size1;
    int c = m -> size2;

    for(int i=0; i<r; i++){
        for(int j=0; j<c; j++) {
            printf("%f ", gsl_matrix_get(m, i, j));
        }
        printf("\n");
    }
}

void print_gsl_matrix_complex(gsl_matrix_complex * m) {

    int r = m -> size1;
    int c = m -> size2;

    for(int i=0; i<r; i++){
        for(int j=0; j<c; j++) {
            gsl_complex z = gsl_matrix_complex_get(m, i, j);
            double real, imag;
            real = GSL_REAL(z);
            imag = GSL_IMAG(z);
            printf("(%f + %fi) ", real, imag);
        }
        printf("\n");
    }
}

void print_array_dyn(int rows, int cols, double * array);
/*----------------------------------------------------------------------------*/

// matrix and vector manipulation, creation and transformation
gsl_matrix * gsl_sq_matrix(size_t sz) {
    gsl_matrix * m = gsl_matrix_alloc(sz, sz);
    return m;
}

gsl_vector * gsl_vector_ones(int n) {
    gsl_vector * v = gsl_vector_alloc(n);
    for(int i=0; i<n; i++) {
        gsl_vector_set(v, i, 1.);
    }
    return v;
}

void gsl_matrix_chrow(gsl_matrix * m, gsl_vector * v, int col_index);
void gsl_matrix_chcol(gsl_matrix * m, gsl_vector * v, int row_index);

gsl_matrix_complex * gsl_sq_matrix_complex(size_t sz);

gsl_vector_complex * gsl_vector_complex_ones(int n) {
    gsl_vector_complex * v = gsl_vector_complex_alloc(n);
    for(int i=0; i<n; i++) {
        gsl_vector_complex_set(v, i, gsl_complex_rect(1., 0.));
    }
    return v;
}

void gsl_matrix_complex_chrow(gsl_matrix_complex * m, gsl_vector_complex * v, int row_index);

void gsl_matrix_complex_chcol(gsl_matrix_complex * m, gsl_vector_complex * v, int col_index) {
    int i, j, rows, cols;
    rows = m -> size1;
    cols = m -> size2;
    for(i=0; i<rows; i++) {
        gsl_complex z = gsl_vector_complex_get(v, i);
        gsl_matrix_complex_set(m, i, col_index, z);
    }
}

/*---------------------------------------------------------------------------*/

// create matrix from different types of array
gsl_matrix_complex * from_gsl_matrix(gsl_matrix * m);

gsl_matrix_complex * from_array_dynamic(int rows, int cols, double * array) {
    gsl_matrix_complex * m = gsl_matrix_complex_alloc(rows, cols);
    for(int i=0; i<rows; i++) {
        for(int j=0; j<cols; j++) {
            gsl_complex z = gsl_complex_rect(array[i*rows + j], 0.);
            gsl_matrix_complex_set(m, i, j, z);
        }
    }
    return m;
}

gsl_matrix_complex * from_array_static(int rows, int cols, double array[rows][cols]);


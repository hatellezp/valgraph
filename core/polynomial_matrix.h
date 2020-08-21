//
// Created by horacio on 08/07/2020.
//

#include <math.h>
#include <gsl/gsl_matrix.h>

#include "polynomial.h"
#include "polynomial_vector.h"
#include "c_errno.h"

#ifndef VALGRAPHCORE_POLYNOMIAL_MATRIX_H
#define VALGRAPHCORE_POLYNOMIAL_MATRIX_H

typedef struct {
    size_t rows;
    size_t cols;
    vgraph_poly ** values;
}vgraph_pmatrix;

vgraph_pmatrix * vgraph_pmatrix_alloc(const size_t rs, const size_t cs);
void vgraph_pmatrix_free(vgraph_pmatrix * pm);

vgraph_pmatrix *
vgraph_pmatrix_from_array_double(double * ar, const size_t rows, const size_t cols, const size_t degree);

void vgraph_pmatrix_set(vgraph_pmatrix * pm, const size_t r, const size_t c, vgraph_poly * p);
vgraph_poly * vgraph_pmatrix_get(vgraph_pmatrix * pm, const size_t r, const size_t c);

void vgraph_pmatrix_copy(vgraph_pmatrix * dest, vgraph_pmatrix * source);

void vgraph_pmatrix_set_zero(vgraph_pmatrix * pm);
void vgraph_pmatrix_set_identity(vgraph_pmatrix * pm);

int vgraph_pmatrix_equal(vgraph_pmatrix * pm1, vgraph_pmatrix * pm2);
void vgraph_pmatrix_add(vgraph_pmatrix * add1, vgraph_pmatrix * add2);
void vgraph_pmatrix_sub(vgraph_pmatrix * sub1, vgraph_pmatrix * sub2);

// multiplication generate a new matrix
vgraph_pmatrix * vgraph_pmatrix_mul(vgraph_pmatrix * mul1, vgraph_pmatrix * mul2);
void vgraph_pmatrix_mul_sq(vgraph_pmatrix * mul1, vgraph_pmatrix * mul2);
void vgraph_pmatrix_scale(vgraph_pmatrix * pm, vgraph_poly * p);

void vgraph_pmatrix_swap_rows(vgraph_pmatrix * pm, const size_t i, const size_t j);
void vgraph_pmatrix_swap_columns(vgraph_pmatrix * pm, const size_t i, const size_t j);
void vgraph_pmatrix_swap_rowcol(vgraph_pmatrix * pm, const size_t i, const size_t j);
void vgraph_pmatrix_transpose (vgraph_pmatrix * pm);

#endif //VALGRAPHCORE_POLYNOMIAL_MATRIX_H

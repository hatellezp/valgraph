//
// Created by horacio on 13/07/2020.
//

#ifndef VALGRAPHCORE_POLYNOMIAL_LINALG_H
#define VALGRAPHCORE_POLYNOMIAL_LINALG_H

#include "polynomial.h"
#include "polynomial_vector.h"
#include "polynomial_matrix.h"

/*
 * My idea to adapt the Gauss-Jordan procedure to a matrix of polynomials is to evaluate
 * polynomials to a value and works as if real numbers.
 */
void vgraph_linalg_GJdecomp(vgraph_pmatrix * A, vgraph_pvector * b, double xval);

#endif //VALGRAPHCORE_POLYNOMIAL_LINALG_H

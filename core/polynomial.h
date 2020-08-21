//
// Created by horacio on 05/07/2020.
//

#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_complex.h>

#include "c_errno.h"

#ifndef VALGRAPHCORE_POLYNOMIAL_H
#define VALGRAPHCORE_POLYNOMIAL_H

/*
 * a polynomial type with a maximum degree and
 * a real effective degree
 */
typedef struct {
    size_t max_degree;
    size_t real_degree;
    gsl_vector * coeffs;
}vgraph_poly;

// struct poly;

vgraph_poly * vgraph_poly_alloc(size_t degree);
void vgraph_poly_free(vgraph_poly * p);

void vgraph_poly_swap(vgraph_poly ** p1, vgraph_poly ** p2);

vgraph_poly * vgraph_poly_from_array(double * ar, size_t size);
vgraph_poly * vgraph_poly_from_double(size_t , double d);
vgraph_poly * vgraph_poly_monome(size_t degree, size_t rdegree, double d);

void vgraph_poly_print(vgraph_poly * p);
void vgraph_poly_detail_print(vgraph_poly * p);

void vgraph_poly_set(vgraph_poly * p, size_t index, double val);
double vgraph_poly_get(vgraph_poly * p, size_t index);

void vgraph_poly_copy(vgraph_poly * dest, vgraph_poly * source);

int vgraph_poly_is_null(vgraph_poly * p);
int vgraph_poly_is_number(vgraph_poly * p);

void vgraph_poly_set_zero(vgraph_poly * p);
void vgraph_poly_set_identity(vgraph_poly * p);

int vgraph_poly_equal(vgraph_poly * p1, vgraph_poly * p2);
void vgraph_poly_add(vgraph_poly * add1, vgraph_poly * add2);
void vgraph_poly_scale(vgraph_poly * p, double d);
void vgraph_poly_sub(vgraph_poly * sub1, vgraph_poly * sub2);
void vgraph_poly_mul(vgraph_poly * mul1, vgraph_poly * mul2);

void vgraph_poly_div(vgraph_poly * dividend, vgraph_poly * divisor, vgraph_poly * quotient, vgraph_poly * rest); // ? implement or not?

/*
 * returns: 1 -> p1(val) > p2(val)
 *          0 -> p1(val) = p2(val)
 *         -1 -> p1(val) < p2(val)
 */
int vgraph_poly_compare(vgraph_poly * p1, vgraph_poly * p2, double val);

double vgraph_poly_eval(vgraph_poly * p, double val, int met);
double vgraph_poly_Cauchy_rbound(vgraph_poly * p);

#endif //VALGRAPHCORE_POLYNOMIAL_H

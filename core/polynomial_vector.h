//
// Created by horacio on 09/07/2020.
//

#include "polynomial.h"
#include "c_errno.h"

#ifndef VALGRAPHCORE_POLYNOMIAL_VECTOR_H
#define VALGRAPHCORE_POLYNOMIAL_VECTOR_H

typedef struct {
    size_t size;
    vgraph_poly ** values;
}vgraph_pvector;

vgraph_pvector * vgraph_pvector_alloc(const size_t s);
void vgraph_pvector_free(vgraph_pvector * pv);

void vgraph_pvector_print(vgraph_pvector * pv);

void vgraph_pvector_set(vgraph_pvector * pv, const size_t index, vgraph_poly * p);
vgraph_poly * vgraph_pvector_get(vgraph_pvector * pv, const size_t index);

void vgraph_pvector_set_zero(vgraph_pvector * pv);
void vgraph_pvector_set_all(vgraph_pvector * pv, vgraph_poly * p);

int vgraph_pvector_is_null(vgraph_pvector * pv);

int vgraph_pvector_equal(vgraph_pvector * pv1, vgraph_pvector * pv2);
void vgraph_pvector_add(vgraph_pvector * add1, vgraph_pvector * add2);
void vgraph_pvector_sub(vgraph_pvector * sub1, vgraph_pvector * sub2);
void vgraph_pvector_scale(vgraph_pvector * pv, vgraph_poly * p);

vgraph_poly * vgraph_pvector_prod(vgraph_pvector * prod1, vgraph_pvector * prod2);

#endif //VALGRAPHCORE_POLYNOMIAL_VECTOR_H

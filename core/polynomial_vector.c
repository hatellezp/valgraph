//
// Created by horacio on 09/07/2020.
//

#include "polynomial_vector.h"

vgraph_pvector * vgraph_pvector_alloc(const size_t s) {
    vgraph_pvector * pv;
    pv = (vgraph_pvector *) malloc(sizeof(vgraph_pvector));

    if (pv == 0) {
        VALGRAPH_ERROR_VAL("failed to allocate for the vgraph_pvector struct", VALGRAPH_NOMEM, 0);
    }

    vgraph_poly ** v;
    v = (vgraph_poly **) malloc(sizeof(vgraph_poly *) * s);

    if (v == 0) {
        VALGRAPH_ERROR_VAL("failed to allocate for the vgraph_pvector struct", VALGRAPH_NOMEM, 0);
    }

    pv -> size = s;
    pv -> values = v;
    return pv;
}

void vgraph_pvector_free(vgraph_pvector * pv) {
    free(pv -> values);
    free(pv);
}

void vgraph_pvector_print(vgraph_pvector * pv) {
    printf("[");
    for(size_t i=0; i<(pv -> size)-1; i++) {
        vgraph_poly_print((pv -> values)[i]);
        printf(", ");
    }
    vgraph_poly_print((pv -> values)[(pv -> size)-1]);
    printf("]");
}

void vgraph_pvector_set(vgraph_pvector * pv, const size_t index, vgraph_poly * p) {
    if (index < (pv -> size)) {
        (pv -> values)[index] = p;
    } else {
        VALGRAPH_ERROR_VOID("index out of bounds", VALGRAPH_OUTOFB);
    }
}

vgraph_poly * vgraph_pvector_get(vgraph_pvector * pv, const size_t index) {
    if (index < (pv -> size)) {
        return (pv -> values)[index];
    } else {
        VALGRAPH_ERROR_VAL("index out of bounds", VALGRAPH_OUTOFB, 0);
    }
}

void vgraph_pvector_set_zero(vgraph_pvector * pv) {
    for(size_t i=0; i<(pv -> size); i++) {
        vgraph_poly_set_zero((pv -> values)[i]);
    }
}

void vgraph_pvector_set_all(vgraph_pvector * pv, vgraph_poly * p) {
    for(size_t i=0; i<(pv -> size); i++) {
       vgraph_poly_copy((pv -> values)[i], p);
    }
}

int vgraph_pvector_is_null(vgraph_pvector * pv) {
    for(size_t i=0; i<(pv -> size); i++) {
        if (vgraph_poly_is_null((pv -> values)[i])) {
            continue;
        } else {
            return 0;
        }
    }
    return 1;
}

int vgraph_pvector_equal(vgraph_pvector * pv1, vgraph_pvector * pv2) {
    size_t s1 = pv1 -> size;
    size_t s2 = pv2 -> size;

    if (s1 == s2) {
        size_t i;
       for(i=0; i<=s1; i++) {
           if (vgraph_poly_equal(vgraph_pvector_get(pv1, i), vgraph_pvector_get(pv2, i))) {
               continue;
           } else {
               return 0;
           }
       }
       return 1;
    } else {
        return 0;
    }
}

void vgraph_pvector_add(vgraph_pvector * add1, vgraph_pvector * add2) {
    if (add1 -> size == add2 -> size) {
        size_t i;
        for(i=0; i < (add1 -> size); i++) {
            vgraph_poly_add((add1 -> values)[i], (add2 -> values)[i]); // use the set and get functions only when needed
        }
    } else {
        VALGRAPH_ERROR_VOID("vectors of different length", VALGRAPH_VBADLEN);
    }
}

void vgraph_pvector_sub(vgraph_pvector * sub1, vgraph_pvector * sub2) {
    if (sub1 -> size == sub2 -> size) {
        size_t i;
        for(i=0; i < (sub1 -> size); i++) {
            vgraph_poly_sub((sub1 -> values)[i], (sub2 -> values)[i]);
        }
    } else {
        VALGRAPH_ERROR_VOID("vectors of different length", VALGRAPH_VBADLEN);
    }
}

void vgraph_pvector_scale(vgraph_pvector * pv, vgraph_poly * p) {
    size_t i;
    for(i=0; i < (pv -> size); i++) {
        vgraph_poly_mul((pv -> values)[i], p);
    }
}

vgraph_poly * vgraph_pvector_prod(vgraph_pvector * prod1, vgraph_pvector * prod2) {

    if (prod1 -> size == prod2 -> size) {

        size_t md, s, i;

        md = (prod1 -> values)[0] -> max_degree;
         s = prod1 -> size;

        vgraph_poly * p = vgraph_poly_alloc(md);
        vgraph_poly * xtemp = vgraph_poly_alloc(md);

        for(i=0; i<s; i++) {
            vgraph_poly_copy(xtemp, (prod1 -> values)[i]); // copy to xtemp
            vgraph_poly_mul(xtemp, (prod2 -> values)[i]); // multiply both components
            vgraph_poly_add(p, xtemp);
        }

        vgraph_poly_free(xtemp);

        return p;
    } else {
        VALGRAPH_ERROR_VAL("vectors of different length", VALGRAPH_VBADLEN, 0);
    }
}



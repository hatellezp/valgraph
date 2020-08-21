//
// Created by horacio on 08/07/2020.
//

#include "polynomial_matrix.h"

vgraph_pmatrix * vgraph_pmatrix_alloc(const size_t rs, const size_t cs) {
    vgraph_pmatrix * pm;
    pm = (vgraph_pmatrix *) malloc(sizeof(vgraph_pmatrix));

    if (pm == 0) {
        VALGRAPH_ERROR_VAL("failed to allocate for vgraph_pmatrix struct", VALGRAPH_NOMEM, 0);
    }

    vgraph_poly ** vs;
    vs = (vgraph_poly **) malloc(sizeof(vgraph_poly * ) * cs * rs);

    if (vs == 0) {
        VALGRAPH_ERROR_VAL("failed to allocate for vgraph_pmatrix struct", VALGRAPH_NOMEM, 0);
    }

    pm -> rows = rs;
    pm -> cols = cs;
    pm -> values = vs;
    return pm;
}

void vgraph_pmatrix_free(vgraph_pmatrix * pm) {
    free(pm -> values);
    free(pm);
}

vgraph_pmatrix *
vgraph_pmatrix_from_array_double(double * ar, const size_t rows, const size_t cols, const size_t degree) {
   vgraph_pmatrix * pm = vgraph_pmatrix_alloc(rows, cols);
   size_t i,j;

   for(i=0; i<rows; i++) {
       for(j=0; j<cols; j++) {
           (pm -> values)[i*rows + j] = vgraph_poly_from_double(degree, ar[i*rows + j]);
       }
   }

   return pm;
}

void vgraph_pmatrix_set(vgraph_pmatrix * pm, const size_t r, const size_t c, vgraph_poly * p) {
    size_t rs = pm -> rows;
    size_t cs = pm -> cols;

    if ((r >= rs) || (c >= cs)) {
        VALGRAPH_ERROR_VOID("index out of bounds", VALGRAPH_OUTOFB);
    }

    (pm -> values)[r*(pm -> rows) + c] = p;
}

vgraph_poly * vgraph_pmatrix_get(vgraph_pmatrix * pm, const size_t r, const size_t  c) {
    size_t rs = pm -> rows;
    size_t cs = pm -> cols;

    if ((r >= rs) || (c >= cs)) {
        VALGRAPH_ERROR_VAL("index out of bounds", VALGRAPH_OUTOFB, 0);
    }

    return (pm -> values)[r*(pm -> rows) + c];
}

void vgraph_pmatrix_copy(vgraph_pmatrix * dest, vgraph_pmatrix * source) {
    size_t rsd, csd, rss, css, i, j;
    rsd = dest -> rows;
    csd = dest -> cols;
    rss = source -> rows;
    css = source -> cols;

    if ((rsd != rss) || (csd != css)) {
        VALGRAPH_ERROR_VOID("matrices of different size", VALGRAPH_MBADSIZE);
    }

    for(i=0; i<rsd; i++) {
        for(j=0; j<csd; j++) {
            (dest -> values)[rsd*i + j] = (source -> values)[rsd*i + j];
        }
    }

}

void vgraph_pmatrix_set_zero(vgraph_pmatrix * pm) {
    size_t i, j;
    for(i=0; i < (pm -> rows); i++) {
        for(j=0; j < (pm -> cols); j++) {
            vgraph_poly_set_zero((pm -> values)[i*(pm -> rows) + j]);
        }
    }
}

void vgraph_pmatrix_set_identity(vgraph_pmatrix * pm) {
    size_t rs, cs, i, j;
    rs = pm -> rows;
    cs = pm -> cols;

    if (rs != cs) {
        VALGRAPH_ERROR_VOID("the pmatrix is not square", VALGRAPH_MNOTSQR);
    }

    /*
    vgraph_pmatrix_set_zero(pm);
    for(int i=0; i < (pm -> rows); i++) {
        vgraph_poly_set_identity(vgraph_pmatrix_get(pm, i, i));
    }
     */

    // faster (even if it's a little faster, faster is faster)
    for(i=0; i<rs; i++) {
        for(j=0; j<cs; j++) {
            if (i == j) {
                vgraph_poly_set_identity((pm -> values)[rs*i + j]);
            } else {
                vgraph_poly_set_zero((pm -> values)[rs*i + j]);
            }
        }
    }

}

int vgraph_pmatrix_equal(vgraph_pmatrix * pm1, vgraph_pmatrix * pm2) {
    size_t rs1, rs2;
    rs1 = pm1 -> rows;
    rs2 = pm2 -> rows;

    if (rs1 == rs2) {
       size_t cs1, cs2;
       cs1 = pm1 -> cols;
       cs2 = pm2 -> cols;

       if (cs1 == cs2) {
          size_t i,j;
          for(i=0; i<rs1; i++) { // thanks CLion for the 'unreachable code' warning
              for(j=0; j<cs1; j++) {
                 if (vgraph_poly_equal((pm1 -> values)[rs1*i + j], (pm2 -> values)[rs1*i + j])) {
                     continue;
                 } else {
                     return 0;
                 }
              }
          }
          return 1;
       } else {
           return 0;
       }
    } else {
        return 0;
    }
}

// all operations are done in place with respect to the first argument
void vgraph_pmatrix_add(vgraph_pmatrix * add1, vgraph_pmatrix * add2) {
    size_t rs1, cs1, rs2, cs2;
    rs1 = add1 -> rows;
    cs1 = add1 -> cols;
    rs2 = add2 -> rows;
    cs2 = add2 -> cols;

    if ((rs1 != rs2) || (cs1 != cs2)) {
        VALGRAPH_ERROR_VOID("matrices of different size", VALGRAPH_MBADSIZE);
    }

    size_t i, j;
    // TODO: does this work !?
    for(i=0; i < rs1; i++) {
        for(j=0; j < cs1; j++) {
            vgraph_poly_add((add1 -> values)[i*rs1 + j], (add2 -> values)[i*rs1 + j]);
        }
    }
}

void vgraph_pmatrix_sub(vgraph_pmatrix * sub1, vgraph_pmatrix * sub2) {
    size_t rs1, cs1, rs2, cs2;
    rs1 = sub1 -> rows;
    cs1 = sub1 -> cols;
    rs2 = sub2 -> rows;
    cs2 = sub2 -> cols;

    if ((rs1 != rs2) || (cs1 != cs2)) {
        VALGRAPH_ERROR_VOID("matrices of different size", VALGRAPH_MBADSIZE);
    }

    size_t i, j;
    for(i=0; i < rs1; i++) {
        for(j=0; j < cs1; j++) {
            vgraph_poly_sub((sub1 -> values)[i*rs1 + j], (sub2 -> values)[i*rs1 + j]);
        }
    }
}

void vgraph_pmatrix_scale(vgraph_pmatrix * pm, vgraph_poly * p) {
    size_t rs, cs, i, j;
    rs = pm -> rows;
    cs = pm -> cols;
    for(i=0; i < rs; i++) {
        for(j=0; j < cs; j++) {
            vgraph_poly_mul((pm -> values)[i*rs + j], p);
        }
    }
}

void vgraph_pmatrix_mul_elements(vgraph_pmatrix * mul1, vgraph_pmatrix * mul2,
        const size_t r1, const size_t c2, vgraph_poly * p) {

    size_t cs1, rs2;
    cs1 = mul1 -> cols;
    rs2 = mul2 -> rows;

    if (cs1 != rs2) {
        VALGRAPH_ERROR_VOID("bad dimension for pmatrix multiplication", VALGRAPH_MBADMUL);
    }

    // put the target polynomial to zero
    vgraph_poly_set_zero(p);
    vgraph_poly * xtemp = vgraph_poly_alloc(p -> max_degree);

    size_t rs1, i;
    rs1 = mul1 -> rows;

    for(i=0; i<cs1; i++) {
       vgraph_poly_copy(xtemp, (mul1 -> values)[rs1*r1 + i]);
       vgraph_poly_mul(xtemp, (mul2 -> values)[rs2*i + c2]);
       vgraph_poly_add(p, xtemp);
    }

    vgraph_poly_free(xtemp);
}

// TODO: verify this method !!!
//     : haven't done the copy and swap between pm and mul1 yet
vgraph_pmatrix * vgraph_pmatrix_mul(vgraph_pmatrix * mul1, vgraph_pmatrix * mul2) {
    size_t rs1, cs1, rs2, cs2;
    rs1 = mul1 -> rows;
    cs1 = mul1 -> cols;
    rs2 = mul2 -> rows;
    cs2 = mul2 -> cols;

    if (cs1 != rs2) {
        VALGRAPH_ERROR_VAL("bad dimension for pmatrix multiplication", VALGRAPH_MBADMUL, 0);
    }

    size_t r1, c2;
    vgraph_pmatrix * pm = vgraph_pmatrix_alloc(rs1, cs2);

    for(r1=0; r1<rs1; r1++) {
        for(c2=0; c2<cs2; c2++) {
            vgraph_pmatrix_mul_elements(mul1, mul2, r1, c2, (pm -> values)[r1*rs1 + c2]);
        }
    }

    return pm;
}

void vgraph_pmatrix_mul_sq(vgraph_pmatrix * mul1, vgraph_pmatrix * mul2) {
    size_t rs1, cs1, rs2, cs2;
    rs1 = mul1 -> rows;
    cs1 = mul1 -> cols;
    rs2 = mul2 -> rows;
    cs2 = mul2 -> cols;

    if (rs1 != rs2 || cs1 != cs2) {
        VALGRAPH_ERROR_VOID("matrix haven't the same dimension", VALGRAPH_MBADMUL);
    }

    vgraph_pmatrix * pm = vgraph_pmatrix_alloc(rs1, cs1);
    pm = vgraph_pmatrix_mul(mul1, mul2);

    vgraph_pmatrix_copy(mul1, pm);
    vgraph_pmatrix_free(pm);
}

// I'm going to use the swap method of polynomials to avoid multiplications

// TODO: verify the swap works as intended ... I don't think it does ...  now I do
void vgraph_pmatrix_swap_rows(vgraph_pmatrix * pm, const size_t i, const size_t j) {
    size_t k, rs, cs;
    rs = pm -> rows;

    if (i >= rs || j >= rs) {
        VALGRAPH_ERROR_VOID("index out of bounds", VALGRAPH_OUTOFB);
    }

    if (i != j) {
        cs = pm -> cols;
        // vgraph_poly * p = vgraph_poly_alloc((pm -> values)[0] -> max_degree);

        for(k=0; k<cs; k++) {
            vgraph_poly_swap(&((pm -> values)[rs*i + k]), &((pm -> values) [rs*j + k]));

            /*
            vgraph_poly_copy(p, (pm -> values)[rs*i + k]);
            vgraph_poly_copy((pm -> values)[rs*i + k], (pm -> values)[rs*j + k]);
            vgraph_poly_copy((pm -> values)[rs*j + k], p);
             */
        }

        // vgraph_poly_free(p);
    }
}

void vgraph_pmatrix_swap_columns(vgraph_pmatrix * pm, const size_t i, const size_t j) {
    size_t k, rs, cs;
    cs = pm -> cols;

    if (i >= cs || j >= cs) {
        VALGRAPH_ERROR_VOID("index out of bounds", VALGRAPH_OUTOFB);
    }

    if (i != j) {
        rs = pm -> rows;
        // vgraph_poly * p = vgraph_poly_alloc((pm -> values)[0] -> max_degree);

        for(k=0; k<cs; k++) {
            vgraph_poly_swap(&((pm -> values)[rs*k + i]), &((pm -> values)[rs*k + j]));

            /*
            vgraph_poly_copy(p, (pm -> values)[rs*k + i]);
            vgraph_poly_copy((pm -> values)[rs*k + i], (pm -> values)[rs*k + j]);
            vgraph_poly_copy((pm -> values)[rs*k + j], p);
             */
        }

        // vgraph_poly_free(p);
    }

}

void vgraph_pmatrix_swap_rowcol(vgraph_pmatrix * pm, const size_t i, const size_t j) {
    size_t rs, cs;
    rs = pm -> rows;
    cs = pm -> cols;

    if (rs != cs)  {
        VALGRAPH_ERROR_VOID("matrix must be square to swap row and column", VALGRAPH_MNOTSQR);
    }

    if (i >= rs || j >= cs) {
        VALGRAPH_ERROR_VOID("index out of bounds", VALGRAPH_OUTOFB);
    }

    // vgraph_poly * p = vgraph_poly_alloc((pm -> values)[0] -> max_degree);
    size_t k, r, c;

    for(k=0; k<rs; k++) {
       r = rs*i + k;
       c = rs*k + j;

       if (r == c) {
           continue;
       } else {
           vgraph_poly_swap(&((pm -> values)[r]), &((pm -> values)[c]));

           /*
           vgraph_poly_copy(p, (pm -> values)[r]);
           vgraph_poly_copy((pm -> values)[r], (pm -> values)[c]);
           vgraph_poly_copy((pm -> values)[c], p);
            */
       }
    }

    // vgraph_poly_free(p);
}

// for the moment I going to follow gsl and transpose only square matrices
// TODO: maybe implement for all matrices
void vgraph_pmatrix_transpose (vgraph_pmatrix * pm) {
   size_t rs, cs;
   rs = pm -> rows;
   cs = pm -> cols;

   if (rs != cs) {
       VALGRAPH_ERROR_VOID("matrix must be square to transpose", VALGRAPH_MNOTSQR);
   }

   // vgraph_poly * p = vgraph_poly_alloc((pm -> values)[0] -> max_degree);
   size_t i,j;

   for(i=0; i<rs; i++) {
       for(j=i+1; j<cs; j++) {
           vgraph_poly_swap(&((pm -> values)[rs*i + j]), &((pm -> values)[rs*j + i]));

           /*
           vgraph_poly_copy(p, (pm -> values)[rs*i + j]);
           vgraph_poly_copy((pm -> values)[rs*i + j], (pm -> values)[rs*j + i]);
           vgraph_poly_copy((pm -> values)[rs*j + i], p);
            */
       }
   }

   // vgraph_poly_free(p);
}



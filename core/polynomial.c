//
// Created by horacio on 05/07/2020.
//

#include "polynomial.h"

// I want these for me here
size_t maxszt(size_t a, size_t b) {
    if (a > b) {
        return a;
    } else {
        return b;
    }
}

size_t minszt(size_t a, size_t b) {
    if (a < b) {
        return a;
    } else {
        return b;
    }
}


double maxdb(double a, double b) {
    if (a > b) {
        return a;
    } else {
        return b;
    }
}

double mindb(double a, double b) {
    if (a < b) {
        return a;
    } else {
        return b;
    }
}


/*----------------------------------------------------------------------------*/

vgraph_poly * vgraph_poly_alloc(const size_t degree) {
    vgraph_poly * p;
    p = (vgraph_poly *) malloc(sizeof(vgraph_poly));

    if (p == 0) {
        VALGRAPH_ERROR_VAL("failed to allocate space for polynomial struct", VALGRAPH_NOMEM, 0);
    }

    p -> real_degree = 0;
    p -> max_degree = degree;
    p -> coeffs = gsl_vector_calloc(degree + 1); // I want all values to zero from the start

    return p;
}

void vgraph_poly_free(vgraph_poly * p) {
    gsl_vector_free(p -> coeffs);
    free(p);
}

void vgraph_poly_swap(vgraph_poly ** p1, vgraph_poly ** p2) {
    vgraph_poly *temp = *p1;
    *p1 = *p2;
    *p2 = temp;
}
void vgraph_poly_recompute_real_degree(vgraph_poly * p, size_t from) {
    size_t f,i;
    if (from == 0) {
        f = p -> max_degree;
    } else {
        f = from;
    }
    for(i=f; (0 <= i && i <= f); i--) {
        if (gsl_vector_get(p -> coeffs, i) != 0.) {
            p -> real_degree = i;
            break;
        } else {
            continue;
        }
    }

    // outside without findind the new real degree means the new real degree is zero
    // but size_t is unsigned, so it overflows, this solution is for the moment
    if (i <= 0 || i > (p -> max_degree)) { // TODO: I don't like at all this solution
        p -> real_degree = 0;
    }
}

vgraph_poly * vgraph_poly_from_array(double * ar, const size_t size) {
    vgraph_poly * p = vgraph_poly_alloc(size - 1); // the degree of an array of size s is s-1
    int ddone = 0;
    size_t i = 0;
    for(i=size-1; (i>= 0 && i <= (size-1)); i--) {
        gsl_vector_set(p -> coeffs, i, ar[i]);

        // not need of doing a second pass
        if (!ddone && ar[i] != 0.) {
            p -> real_degree = i;
            ddone = 1;
        }
    }

    // vgraph_poly_recompute_real_degree(p, -1);
    return p;
}

void vgraph_poly_print(vgraph_poly * p) {
    if (vgraph_poly_is_null(p)) {
        printf("0.");
    } else if (vgraph_poly_is_number(p)) {
       printf("%lf", gsl_vector_get(p -> coeffs, 0));
    }
    else {

        double val;
        size_t rd = p -> real_degree;
        val = gsl_vector_get(p -> coeffs, rd);

        if (rd == 1) {
            printf("%lf*x", val);
        } else {
            printf("%lf*x^%zu", val, rd);
        }

        size_t i;
        for(i=(rd-1); (i >= 0 && i <= (rd-1)); i--) {
            val = gsl_vector_get(p -> coeffs, i);
            if (val == 0.) {
                continue;
            } else {
                if (i == 1) {
                    printf(" + %lf*x", val);
                } else if (i == 0) {
                    printf(" + %lf", val);
                } else {
                    printf(" + %lf*x^%zu", val, i); // TODO: see what this does
                }
            }
        }
    }
}

// not implemented yet
void vgraph_poly_detail_print(vgraph_poly * p) {
    VALGRAPH_ERROR_VOID("method not implemented", VALGRAPH_NOIMPL);
}

vgraph_poly * vgraph_poly_from_double(const size_t degree, double d) {
    vgraph_poly * p = vgraph_poly_alloc(degree);
    gsl_vector_set(p -> coeffs, 0, d);
    return p;
}

vgraph_poly * vgraph_poly_monome(const size_t degree, const size_t rdegree, double d) {
    if (rdegree > degree) {
        VALGRAPH_ERROR_VAL("the current degree must be lower or equal than the maximal degree", VALGRAPH_PBADDEG, 0);
    }

    vgraph_poly * p = vgraph_poly_alloc(degree);
    gsl_vector_set(p -> coeffs, rdegree, d);
    p -> real_degree = rdegree;
    return p;
}

size_t vgraph_poly_get_max_real_degree(vgraph_poly * p1, vgraph_poly * p2) {
    size_t max_real_degree = maxszt(p1 -> real_degree, p2 -> real_degree);
    return max_real_degree;
}

void vgraph_poly_set(vgraph_poly * p, const size_t index, double val) {
    if (index <= p -> max_degree) {

        size_t rd = p -> real_degree;

        if (index == rd && val == 0.) { // TODO: verify this
            gsl_vector_set(p -> coeffs, index, val);
            vgraph_poly_recompute_real_degree(p, rd);
        } else if (index >= rd) {
            gsl_vector_set(p -> coeffs, index, val);
            p -> real_degree = index;
        } else {
            gsl_vector_set(p -> coeffs, index, val);
        }

    } else {
        VALGRAPH_ERROR_VOID("index out of bounds", VALGRAPH_OUTOFB);
    }
}

// here the function I think will output an error using the gsl_vector interface
double vgraph_poly_get(vgraph_poly * p, const size_t index) {
    if (index > (p -> max_degree)) {
        VALGRAPH_ERROR_VAL("index out of bounds", VALGRAPH_OUTOFB, 0);
    }
    double val = gsl_vector_get(p -> coeffs, index);
    return val;
}

// all operations are done in place, the first argument is changed because of the operation
void vgraph_poly_copy(vgraph_poly * dest, vgraph_poly * source) {
    if (dest -> max_degree == source -> max_degree) {
        size_t max_real_degree = vgraph_poly_get_max_real_degree(dest, source);
        size_t i;
        for(i=0; i<= max_real_degree; i++) { // here there is no problem, we begin at 0 and go up
            gsl_vector_set(dest -> coeffs, i, gsl_vector_get(source -> coeffs, i));
        }
        dest -> real_degree = source -> real_degree; // TODO: be sure of this
    } else {
        VALGRAPH_ERROR_VOID("polynomials of different length", VALGRAPH_PBADDEG);
    }
}

int vgraph_poly_is_null(vgraph_poly * p) {
    return (p -> real_degree) == 0 && gsl_vector_get(p -> coeffs, 0) == 0.;
}

int vgraph_poly_is_number(vgraph_poly * p) {
    return (p -> real_degree) == 0;
}

void vgraph_poly_set_zero(vgraph_poly * p) {
    size_t i;
    for(i=0; i<= (p -> real_degree); i++) {
        gsl_vector_set(p -> coeffs, i, 0.);
    }
    p -> real_degree = 0;
}

void vgraph_poly_set_identity(vgraph_poly * p) {
    gsl_vector_set(p -> coeffs, 0, 1.);
    size_t i;
    for(i=1; i<= (p -> real_degree); i++) {
        gsl_vector_set(p -> coeffs, i, 0.);
    }
    p -> real_degree = 0;
}

int vgraph_poly_equal(vgraph_poly * p1, vgraph_poly * p2) {
    size_t md1 = p1 -> max_degree;
    size_t md2 = p2 -> max_degree;

    if (md1 == md2) {
        size_t rd1 = p1 -> real_degree;
        size_t rd2 = p2 -> real_degree;

        if (rd1 == rd2) {
            size_t i;
           for(i=0; i <= rd1; i++) {
               if (gsl_vector_get(p1 -> coeffs, i) != gsl_vector_get(p2 -> coeffs, i)) {
                   return 0;
               } else {
                   continue;
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

// TODO: repair this and access with gsl_vector_(set|get)
//     : super expensive calculation each set induces real_degree recomputation
void vgraph_poly_add(vgraph_poly * add1, vgraph_poly * add2) {
    if (add1 -> max_degree == add2 -> max_degree) {

        size_t max_real_degree = vgraph_poly_get_max_real_degree(add1, add2);
        size_t i;
        for(i=0; i<= max_real_degree; i++) {
            gsl_vector_set(add1 -> coeffs, i, gsl_vector_get(add1 -> coeffs, i) + gsl_vector_get(add2 -> coeffs, i));
        }

        // the sum made the leading coefficient vanish
        if (add1 -> real_degree == add2 -> real_degree && gsl_vector_get(add1 -> coeffs, add1 -> real_degree) == 0.) {
           vgraph_poly_recompute_real_degree(add1, add1 -> real_degree);
        } else if (add1 -> real_degree != add2 -> real_degree) { // the coeff does not vanish, take the leading one
            add1 -> real_degree = max_real_degree;
        }

    } else {
        VALGRAPH_ERROR_VOID("polynomials of different length", VALGRAPH_PBADDEG);
    }
}

void vgraph_poly_scale(vgraph_poly * p, double d) {
    size_t i;
    for(i=0; i <= (p -> real_degree); i++) {
        gsl_vector_set(p -> coeffs, i, d * gsl_vector_get(p -> coeffs, i));
    }
}

void vgraph_poly_sub2(vgraph_poly * sub1, vgraph_poly * sub2) {
    vgraph_poly * p = vgraph_poly_alloc(sub2 -> max_degree);
    vgraph_poly_copy(p, sub2);
    vgraph_poly_scale(p, -1.);
    vgraph_poly_add(sub1, p);
    vgraph_poly_free(p);
}

// this version is less beautiful, the first one uses the already written functions, but I think this one
// is the optimised one
void vgraph_poly_sub(vgraph_poly * sub1, vgraph_poly * sub2) {
    if (sub1 -> max_degree == sub2 -> max_degree) {

        size_t max_real_degree = vgraph_poly_get_max_real_degree(sub1, sub2);
        size_t i;
        for(i=0; i <= max_real_degree; i++) {
            gsl_vector_set(sub1 -> coeffs, i, gsl_vector_get(sub1 -> coeffs, i) - gsl_vector_get(sub2 -> coeffs, i));
        }

        // verify which coeff is the leading coeff
        if (sub1 -> real_degree == sub2 -> real_degree && gsl_vector_get(sub1 -> coeffs, sub1 -> real_degree) == 0.) {
            vgraph_poly_recompute_real_degree(sub1, sub1 -> real_degree);
        } else if (sub1 -> real_degree != sub2 -> real_degree) {
            sub1 -> real_degree = max_real_degree;
        }

    } else {
        VALGRAPH_ERROR_VOID("polynomials of different length", VALGRAPH_PBADDEG);
    }
}

void vgraph_poly_mul(vgraph_poly * mul1, vgraph_poly * mul2) {
    if (mul1 -> max_degree == mul2 -> max_degree) {
        size_t imp_degree = minszt(mul1 -> max_degree, (mul1 -> real_degree) + (mul2 -> real_degree));

        if ((mul1 -> max_degree) < (mul1 -> real_degree) + (mul2 -> real_degree)) {
            printf("\nWARNING: multiplication will result in higher degree polynomial, cast will occur.\n");
        }

        size_t i, j;
        double xtemp;
        vgraph_poly * p = vgraph_poly_alloc(mul1 -> max_degree);

        for(i=0; i <= imp_degree; i++) { // for each i compute the i-th coefficient
            xtemp = 0;
            for(j=0; j <= i; j++) {
                xtemp += gsl_vector_get(mul1 -> coeffs, j) * gsl_vector_get(mul2 -> coeffs, i-j);
            }
            gsl_vector_set(p -> coeffs, i, xtemp);
        }

        vgraph_poly_copy(mul1, p);
        vgraph_poly_free(p);

    } else {
        VALGRAPH_ERROR_VOID("polynomials of different length", VALGRAPH_PBADDEG);
    }
}

// TODO: verify method
void vgraph_poly_div(vgraph_poly * dividend, vgraph_poly * divisor, vgraph_poly * quotient, vgraph_poly * rest) {
    size_t amd, bmd, qmd, rmd;
    amd = dividend -> max_degree;
    bmd = divisor -> max_degree;
    qmd = quotient -> max_degree;
    rmd = rest -> max_degree;

    if (amd != bmd || amd != qmd || amd != rmd) {
        VALGRAPH_ERROR_VOID("polynomials of different max degree", VALGRAPH_PBADDEG);
    }

    size_t ard, brd;
    ard = dividend -> real_degree;
    brd = divisor -> real_degree;

    vgraph_poly_set_zero(quotient); // make the quotient zero, this works even if deg(divisor) > deg(dividend)

    if (ard < brd) {  // if the divisor is of a higher degree the copy to the rest and the quotient remains zero
        vgraph_poly_copy(rest, divisor);
    } else {
        vgraph_poly_set_zero(rest);

        vgraph_poly * xtemp1 = vgraph_poly_alloc(amd);
        vgraph_poly * xtemp2 = vgraph_poly_alloc(amd);

        vgraph_poly_copy(xtemp1, dividend);
        size_t rd1 = xtemp1 -> real_degree;
        size_t currentd;
        double rcoeff;

        while(rd1 >= brd) {
            vgraph_poly_set_zero(xtemp2);

            currentd = (xtemp1 -> real_degree) - (divisor -> real_degree);
            rcoeff = gsl_vector_get((xtemp1 -> coeffs), rd1) /
                        gsl_vector_get((divisor -> coeffs), divisor -> real_degree);

            // create new element of the sum
            // we do this manually
            gsl_vector_set((xtemp2 -> coeffs), currentd, rcoeff);
            xtemp2 -> real_degree = currentd;

            // add to the quotient
            vgraph_poly_add(quotient, xtemp2);

            // transform xtemp2 in the next substractor
            vgraph_poly_mul(xtemp2, divisor);

            // substract from xtemp1
            vgraph_poly_sub(xtemp1, xtemp2);
            rd1 = xtemp1 -> real_degree;
        }

        vgraph_poly_copy(rest, xtemp1);

        vgraph_poly_free(xtemp1);
        vgraph_poly_free(xtemp2);
    }
}

// TODO: verify this method
double vgraph_poly_eval(vgraph_poly * p, double val, int met) {
    if (met == -1) {
        return 1.;
    } else if (vgraph_poly_is_number(p)) {
        return gsl_vector_get(p -> coeffs, 0);
    } else { // other methods maybe ?
        double res = gsl_vector_get(p -> coeffs, 0);
        for(size_t i=1; i<=(p -> real_degree); i++) {
            res += gsl_vector_get(p -> coeffs, i) * pow(val, (double)i);
        }
        return res;
    }
}

int vgraph_poly_compare(vgraph_poly * p1, vgraph_poly * p2, const double val) {
    vgraph_poly * p = vgraph_poly_alloc(p1 -> max_degree);

    vgraph_poly_copy(p, p1); // put p1 in p
    vgraph_poly_sub(p, p2); // put p1-p2 in p

    double res = vgraph_poly_eval(p, val, 1);

    // don't forget to free p, you don't need it anymore
    vgraph_poly_free(p);

    if (res > 0.) { // p1(val) > p2(val)
       return 1;
    } else if (res == 0.) { // p1(val) = p2(val)
        return  0;
    } else {
        return -1;
    }
}

double vgraph_poly_Cauchy__rbound(vgraph_poly * p) {
    if (vgraph_poly_is_number(p)) {
        return 0.;
    } else {
        double den = fabs(gsl_vector_get(p -> coeffs, p -> real_degree));

        if (den == 0.) {
            VALGRAPH_ERROR_VAL("the leading coefficient can't be zero", VALGRAPH_PCORRP, 0);
        }

        double maxv = GSL_FLT_MIN;
        size_t i;
        double num;

        for(i=0; i<=(p -> real_degree -1); i++) {
            num = fabs(gsl_vector_get(p -> coeffs, i));
            maxv = maxdb(maxv, num / den);
        }

        return 1+maxv;
    }
}




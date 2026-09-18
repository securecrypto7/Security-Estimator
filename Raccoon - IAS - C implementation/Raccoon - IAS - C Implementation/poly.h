#ifndef IAS_POLY_H
#define IAS_POLY_H

#include <stdint.h>
#include "params.h"
#include "fips202.h"

typedef unsigned __int128 zq_t;   /* values in [0, q)            */
typedef __int128          zs_t;   /* signed / centered values    */

extern zq_t   IAS_Q;              /* the prime modulus            */
extern int    IAS_LOGQ_CEIL;     /* exact ceil(log2 q) = bitlen  */
extern int    IAS_LOGQ_SIZE;     
extern zq_t   IAS_NINV;          /* n^{-1} mod q                 */

typedef struct { zq_t c[IAS_N]; } poly;    /* canonical-unsigned coeffs   */
typedef struct { zs_t c[IAS_N]; } spoly;   /* signed (short) coeffs       */

/* ---- initialisation: compute q, zetas, n^{-1} ---- */
void ias_poly_init(void);

/* ---- modular scalar ops (mod IAS_Q) ---- */
zq_t ias_mulmod(zq_t a, zq_t b);
zq_t ias_addmod(zq_t a, zq_t b);
zq_t ias_submod(zq_t a, zq_t b);
zq_t ias_powmod(zq_t a, zq_t e);      /* a^e mod q, e up to 128-bit  */

/* ---- polynomial ops ---- */
void poly_zero(poly *r);
void poly_copy(poly *r, const poly *a);
void poly_add(poly *r, const poly *a, const poly *b);
void poly_sub(poly *r, const poly *a, const poly *b);

void poly_ntt(poly *a);                        /* in place, natural -> NTT   */
void poly_invntt(poly *a);                     /* in place, NTT -> natural   */
void poly_pointwise(poly *r, const poly *a, const poly *b);  /* NTT-domain  */
void poly_mul(poly *r, const poly *a, const poly *b);        /* natural dom */
void poly_mul_school(poly *r, const poly *a, const poly *b); /* reference   */

/* signed <-> unsigned */
void poly_from_signed(poly *r, const spoly *a);   /* reduce mod q          */
void poly_center(spoly *r, const poly *a);        /* lift to (-q/2, q/2]   */

zq_t ias_q_nu(int nu);
void poly_round(poly *r, const poly *a, int nu);        /* r in [0,q_nu)    */
void poly_lift_shift(poly *r, const poly *a, int nu);

long double poly_sqnorm_shift(const poly *a, int nu);

long double spoly_sqnorm(const spoly *a);

void poly_uniform(poly *r, keccak_state *st);

int u128_to_str(char *buf, zq_t x);

#endif

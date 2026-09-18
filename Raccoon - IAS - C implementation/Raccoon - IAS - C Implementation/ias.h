#ifndef IAS_IAS_H
#define IAS_IAS_H

#include <stddef.h>
#include "poly.h"
#include "sample.h"

/* ---- flat-array layout helpers ---- */
#define IAS_COMMIT_POLYS   (IAS_REP * IAS_K)   /* per signer: w_{i,b}, b in [rep] */
#define IAS_STATE_POLYS    (IAS_REP * IAS_L)   /* per signer: r_{i,b}             */

static inline int ias_cidx(int i, int b, int kk)
{ return ((i * IAS_REP) + b) * IAS_K + kk; }
static inline int ias_sidx(int b, int ll)
{ return b * IAS_L + ll; }

/* ---- public parameters (matrix A, natural + NTT domain) ---- */
typedef struct {
    uint8_t seed[IAS_SEEDBYTES];
    poly    A[IAS_K][IAS_L];        
    poly    Antt[IAS_K][IAS_L];    
} ias_pp;

/* ---- keys ---- */
typedef struct { spoly s[IAS_L]; } ias_sk;   /* short secret               */
typedef struct {
    poly tT[IAS_K];        
    poly t_full[IAS_K];    
} ias_pk;

/* ---- aggregate key ---- */
typedef struct {
    poly ttilde[IAS_K];            /* t~                                    */
    poly a[IAS_NSIGNERS];          /* a_i in T (signed monomials)           */
} ias_agg;

/* ---- coordinator context (round 1 output) ---- */
typedef struct {
    poly beta[IAS_REP];                    /* (beta_b), beta_1 = 1          */
    poly ttilde[IAS_K];                    /* t~                            */
    poly w[IAS_K];                         /* aggregate rounded commitment  */
    poly a[IAS_NSIGNERS];                  /* a_i                           */
    poly wi[IAS_NSIGNERS][IAS_K];          /* per-signer w_i (R_q, unrounded)*/
} ias_coord_ctx;

/* ---- per-signer response ---- */
typedef struct { poly c; spoly z[IAS_L]; } ias_resp;

/* ---- aggregate signature ---- */
typedef struct { poly c; spoly z[IAS_L]; poly h[IAS_K]; } ias_sig;

/* ---- scheme bounds (computed from estimator formulas at runtime) ---- */
typedef struct {
    long double B2ind;   /* per-signer L2 bound                             */
    long double B2;      /* aggregate L2 bound                              */
    long double beta;    /* forgery / SIS bound (informational)            */
    int         z_bits;  /* ceil(log2(2*B2+1)) -> bits per z coefficient    */
    int         t_bits;  /* logq_ceil - nu_t  -> bits per t_i^T coefficient */
} ias_bounds;

void ias_compute_bounds(ias_bounds *b);

/* ---- algorithm ---- */
void ias_setup(ias_pp *pp, const uint8_t seed[IAS_SEEDBYTES]);

/* KeyGen: draws s_i,e_i from the *currently seeded* secret RNG            */
void ias_keygen(const ias_pp *pp, ias_sk *sk, ias_pk *pk);

void ias_keyagg(const ias_pp *pp, const ias_pk *pks, int nsigners,
                ias_agg *agg);

/* Signer round 1: commit[] must have IAS_COMMIT_POLYS polys,
 * st[]     must have IAS_STATE_POLYS  spolys.                             */
void ias_sign_pre(const ias_pp *pp, poly *commit, spoly *st);

/* Coordinator round 1.  commits: flat nsigners*rep*k polys.              */
void ias_coord(const ias_pp *pp, const ias_pk *pks, const poly *commits,
               int nsigners, const uint8_t *msgs, size_t mlen_each,
               ias_coord_ctx *cc);

/* Signer round 2 (index idx).  st: this signer's rep*l spolys.           */
void ias_sign_resp(const ias_pp *pp, const ias_sk *sk, int idx,
                   const ias_agg *agg, const ias_coord_ctx *cc,
                   const spoly *st, const uint8_t *msgs, size_t mlen_each,
                   ias_resp *resp);

/* Coordinator combine.  Returns 0 on success, -1 on norm/consistency abort. */
int ias_coord_combine(const ias_pp *pp, const ias_pk *pks,
                      const ias_coord_ctx *cc, const ias_resp *resps,
                      int nsigners, const ias_bounds *bnd, ias_sig *sig);

/* Verify.  Returns 1 = accept, 0 = reject.                               */
int ias_verify(const ias_pp *pp, const ias_agg *agg,
               const uint8_t *msgs, size_t mlen_each, int nsigners,
               const ias_bounds *bnd, const ias_sig *sig);

#endif

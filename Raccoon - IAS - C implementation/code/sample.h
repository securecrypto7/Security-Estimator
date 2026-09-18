#ifndef IAS_SAMPLE_H
#define IAS_SAMPLE_H

#include "poly.h"
#include "fips202.h"

/* ---- secret randomness stream ---- */
void ias_rng_seed(const uint8_t seed[IAS_SEEDBYTES]);
void ias_rng_bytes(uint8_t *out, size_t len);

void sample_gauss(spoly *r, double sigma);

void absorb_poly(keccak_state *st, const poly *p);
void absorb_polyvec(keccak_state *st, const poly *v, int len);
void absorb_u32(keccak_state *st, uint32_t x);

void monomial_poly(poly *p, int pos, int sign);

void H_agg(poly *a, const poly *PK, int nsigners, int i, const poly *ti_k);
void H_beta(poly *beta, int rep,
            const poly *t_all, int nsigners,       
            const uint8_t *msgs, size_t msglen_each,
            const poly *wbar, int wbar_len);
void H_c(poly *c, const poly *ttilde, int tk,
         const uint8_t *msgs, size_t msglen_each, int nsigners,
         const poly *w, int wk);

#endif

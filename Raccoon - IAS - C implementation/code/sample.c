#include "sample.h"
#include <math.h>
#include <string.h>

/* ---------------- secret randomness stream ---------------- */
static keccak_state grng;

void ias_rng_seed(const uint8_t seed[IAS_SEEDBYTES])
{
    shake256_init(&grng);
    const uint8_t dom = 0x00;
    shake256_absorb(&grng, &dom, 1);
    shake256_absorb(&grng, seed, IAS_SEEDBYTES);
    shake256_finalize(&grng);
}
void ias_rng_bytes(uint8_t *out, size_t len) { shake256_squeeze(out, len, &grng); }

/* uniform double in (0,1) from 52 random bits */
static double rng_unit(void)
{
    uint8_t b[7];
    ias_rng_bytes(b, 7);
    uint64_t x = 0;
    for (int i = 0; i < 7; i++) x |= (uint64_t)b[i] << (8 * i);
    x &= (((uint64_t)1 << 52) - 1);
    double u = (double)x / (double)((uint64_t)1 << 52);
    if (u <= 0.0) u = 1e-15;
    return u;
}

void sample_gauss(spoly *r, double sigma)
{
    for (int i = 0; i < IAS_N; i += 2) {
        double u1 = rng_unit(), u2 = rng_unit();
        double mag = sigma * sqrt(-2.0 * log(u1));
        double z0 = mag * cos(2.0 * M_PI * u2);
        double z1 = mag * sin(2.0 * M_PI * u2);
        r->c[i]   = (zs_t)llroundl((long double)z0);
        if (i + 1 < IAS_N) r->c[i + 1] = (zs_t)llroundl((long double)z1);
    }
}

/* ---------------- absorb helpers ---------------- */
void absorb_u32(keccak_state *st, uint32_t x)
{
    uint8_t b[4] = { (uint8_t)x, (uint8_t)(x >> 8), (uint8_t)(x >> 16), (uint8_t)(x >> 24) };
    shake256_absorb(st, b, 4);
}
void absorb_poly(keccak_state *st, const poly *p)
{
    int nb = (IAS_LOGQ_CEIL + 7) / 8;
    uint8_t buf[16];
    for (int i = 0; i < IAS_N; i++) {
        zq_t v = p->c[i];
        for (int b = 0; b < nb; b++) { buf[b] = (uint8_t)(v & 0xff); v >>= 8; }
        shake256_absorb(st, buf, nb);
    }
}
void absorb_polyvec(keccak_state *st, const poly *v, int len)
{ for (int i = 0; i < len; i++) absorb_poly(st, &v[i]); }

/* ---------------- signed monomial ---------------- */
void monomial_poly(poly *p, int pos, int sign)
{
    poly_zero(p);
    p->c[pos] = (sign >= 0) ? 1 : (IAS_Q - 1);
}

/* derive one signed monomial from 4 squeezed bytes */
static void mono_from_stream(poly *p, keccak_state *st)
{
    uint8_t b[4];
    shake256_squeeze(b, 4, st);
    uint32_t x = (uint32_t)b[0] | ((uint32_t)b[1] << 8) |
                 ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
    int pos = (int)(x % IAS_N);
    int sign = (x & (1u << 31)) ? -1 : 1;
    monomial_poly(p, pos, sign);
}

/* ---------------- H_agg ---------------- */
void H_agg(poly *a, const poly *PK, int nsigners, int i, const poly *ti_k)
{
    keccak_state st;
    shake256_init(&st);
    const uint8_t dom = 0x01;
    shake256_absorb(&st, &dom, 1);
    absorb_polyvec(&st, PK, nsigners * IAS_K);   
    absorb_u32(&st, (uint32_t)i);
    absorb_polyvec(&st, ti_k, IAS_K);
    shake256_finalize(&st);
    mono_from_stream(a, &st);
}

/* ---------------- H_beta ---------------- */
void H_beta(poly *beta, int rep,
            const poly *t_all, int nsigners,
            const uint8_t *msgs, size_t msglen_each,
            const poly *wbar, int wbar_len)
{
    keccak_state st;
    shake256_init(&st);
    const uint8_t dom = 0x02;
    shake256_absorb(&st, &dom, 1);
    absorb_polyvec(&st, t_all, nsigners * IAS_K);
    if (msgs && msglen_each)
        shake256_absorb(&st, msgs, msglen_each * (size_t)nsigners);
    absorb_polyvec(&st, wbar, wbar_len);
    shake256_finalize(&st);
    monomial_poly(&beta[0], 0, +1);
    for (int b = 1; b < rep; b++) mono_from_stream(&beta[b], &st);
}

/* ---------------- H_c : weight-omega ternary challenge ---------------- */
void H_c(poly *c, const poly *ttilde, int tk,
         const uint8_t *msgs, size_t msglen_each, int nsigners,
         const poly *w, int wk)
{
    keccak_state st;
    shake256_init(&st);
    const uint8_t dom = 0x03;
    shake256_absorb(&st, &dom, 1);
    absorb_polyvec(&st, ttilde, tk);
    if (msgs && msglen_each)
        shake256_absorb(&st, msgs, msglen_each * (size_t)nsigners);
    absorb_polyvec(&st, w, wk);
    shake256_finalize(&st);

    poly_zero(c);
    signed char sgn[IAS_N];
    for (int i = 0; i < IAS_N; i++) sgn[i] = 0;
    int placed = 0;
    while (placed < IAS_OMEGA) {
        uint8_t b[2];
        shake256_squeeze(b, 2, &st);
        int j = ((int)b[0] | ((int)(b[1] & 1) << 8)) % IAS_N; /* 9-bit index */
        if (sgn[j] != 0) continue;
        sgn[j] = (b[1] & 2) ? -1 : 1;
        placed++;
    }
    for (int i = 0; i < IAS_N; i++)
        if (sgn[i]) c->c[i] = (sgn[i] > 0) ? 1 : (IAS_Q - 1);
}

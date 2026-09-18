#include "poly.h"
#include <string.h>
#include <math.h>

zq_t IAS_Q;
int  IAS_LOGQ_CEIL;
int  IAS_LOGQ_SIZE;
zq_t IAS_NINV;

static zq_t zetas[IAS_N];        
static int  LOGN;

static inline zq_t mulmod_m(zq_t a, zq_t b, zq_t m)
{
    zq_t blo = b & (((zq_t)1 << 35) - 1);
    zq_t bhi = b >> 35;
    zq_t r = (((a * bhi) % m) << 35) % m;
    r = (r + (a * blo) % m) % m;
    return r;
}

zq_t ias_mulmod(zq_t a, zq_t b) { return mulmod_m(a, b, IAS_Q); }
zq_t ias_addmod(zq_t a, zq_t b) { zq_t r = a + b; return r >= IAS_Q ? r - IAS_Q : r; }
zq_t ias_submod(zq_t a, zq_t b) { return a >= b ? a - b : a + IAS_Q - b; }

static zq_t powmod_m(zq_t a, zq_t e, zq_t m)
{
    zq_t r = 1 % m;
    a %= m;
    while (e) {
        if (e & 1) r = mulmod_m(r, a, m);
        a = mulmod_m(a, a, m);
        e >>= 1;
    }
    return r;
}
zq_t ias_powmod(zq_t a, zq_t e) { return powmod_m(a, e, IAS_Q); }

/* ---------- deterministic Miller-Rabin for m < 2^70 ---------- */
static int is_prime_u128(zq_t n)
{
    if (n < 2) return 0;
    static const uint32_t small[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (int i = 0; i < 12; i++) {
        if (n == small[i]) return 1;
        if (n % small[i] == 0) return 0;
    }
    zq_t d = n - 1;
    int r = 0;
    while ((d & 1) == 0) { d >>= 1; r++; }
    for (int i = 0; i < 12; i++) {
        zq_t a = small[i] % n;
        if (a == 0) continue;
        zq_t x = powmod_m(a, d, n);
        if (x == 1 || x == n - 1) continue;
        int composite = 1;
        for (int j = 0; j < r - 1; j++) {
            x = mulmod_m(x, x, n);
            if (x == n - 1) { composite = 0; break; }
        }
        if (composite) return 0;
    }
    return 1;
}

/* ntt_prime(logq): smallest prime >= 2^logq with q = 1 (mod 2n).*/
static zq_t ntt_prime(int logq)
{
    zq_t step = 2 * (zq_t)IAS_N;
    zq_t q = ((((zq_t)1 << logq) / step) * step) + 1;
    if (q < ((zq_t)1 << logq)) q += step;
    while (!is_prime_u128(q)) q += step;
    return q;
}

static unsigned brv(unsigned x, int bits)
{
    unsigned r = 0;
    for (int i = 0; i < bits; i++) { r = (r << 1) | (x & 1); x >>= 1; }
    return r;
}

void ias_poly_init(void)
{
    IAS_Q = ntt_prime(IAS_LOGQ);

    IAS_LOGQ_CEIL = 0;
    { zq_t t = IAS_Q - 1; while (t) { IAS_LOGQ_CEIL++; t >>= 1; } }

    IAS_LOGQ_SIZE = (int)ceil(log2((double)IAS_Q));

    LOGN = 0; { int t = IAS_N; while (t > 1) { t >>= 1; LOGN++; } }

    zq_t psi = 0;
    zq_t exp = (IAS_Q - 1) / (2 * (zq_t)IAS_N);
    for (zq_t g = 2;; g++) {
        zq_t cand = powmod_m(g, exp, IAS_Q);
        if (powmod_m(cand, IAS_N, IAS_Q) == IAS_Q - 1) { psi = cand; break; }
    }
    for (int i = 0; i < IAS_N; i++)
        zetas[i] = powmod_m(psi, brv(i, LOGN), IAS_Q);

    zq_t half = (IAS_Q + 1) / 2;
    IAS_NINV = 1 % IAS_Q;
    for (int i = 0; i < LOGN; i++) IAS_NINV = mulmod_m(IAS_NINV, half, IAS_Q);
}

/* ---------- polynomial basics ---------- */
void poly_zero(poly *r) { for (int i = 0; i < IAS_N; i++) r->c[i] = 0; }
void poly_copy(poly *r, const poly *a) { memcpy(r, a, sizeof(poly)); }
void poly_add(poly *r, const poly *a, const poly *b)
{ for (int i = 0; i < IAS_N; i++) r->c[i] = ias_addmod(a->c[i], b->c[i]); }
void poly_sub(poly *r, const poly *a, const poly *b)
{ for (int i = 0; i < IAS_N; i++) r->c[i] = ias_submod(a->c[i], b->c[i]); }

/* forward NTT, Cooley-Tukey, full splitting (down to len=1). */
void poly_ntt(poly *a)
{
    unsigned len, start, j, k = 1;
    for (len = IAS_N / 2; len >= 1; len >>= 1) {
        for (start = 0; start < IAS_N; start = j + len) {
            zq_t zeta = zetas[k++];
            for (j = start; j < start + len; j++) {
                zq_t t = ias_mulmod(zeta, a->c[j + len]);
                a->c[j + len] = ias_submod(a->c[j], t);
                a->c[j]       = ias_addmod(a->c[j], t);
            }
        }
        if (len == 1) break;
    }
}

/* inverse NTT, Gentleman-Sande, using -zetas in reverse order. */
void poly_invntt(poly *a)
{
    unsigned len, start, j, k = IAS_N;
    for (len = 1; len < IAS_N; len <<= 1) {
        for (start = 0; start < IAS_N; start = j + len) {
            zq_t zeta = IAS_Q - zetas[--k];        /* -psi^{...} */
            for (j = start; j < start + len; j++) {
                zq_t t = a->c[j];
                zq_t u = a->c[j + len];
                a->c[j]       = ias_addmod(t, u);
                a->c[j + len] = ias_mulmod(zeta, ias_submod(t, u));
            }
        }
    }
    for (j = 0; j < IAS_N; j++) a->c[j] = ias_mulmod(a->c[j], IAS_NINV);
}

void poly_pointwise(poly *r, const poly *a, const poly *b)
{ for (int i = 0; i < IAS_N; i++) r->c[i] = ias_mulmod(a->c[i], b->c[i]); }

void poly_mul(poly *r, const poly *a, const poly *b)
{
    poly ta, tb;
    poly_copy(&ta, a); poly_copy(&tb, b);
    poly_ntt(&ta); poly_ntt(&tb);
    poly_pointwise(r, &ta, &tb);
    poly_invntt(r);
}


void poly_mul_school(poly *r, const poly *a, const poly *b)
{
    zq_t acc[2 * IAS_N];
    for (int i = 0; i < 2 * IAS_N; i++) acc[i] = 0;
    for (int i = 0; i < IAS_N; i++)
        for (int j = 0; j < IAS_N; j++)
            acc[i + j] = ias_addmod(acc[i + j], ias_mulmod(a->c[i], b->c[j]));
    for (int i = 0; i < IAS_N; i++)
        r->c[i] = ias_submod(acc[i], acc[i + IAS_N]);   /* X^n = -1 */
}

/* ---------- signed <-> unsigned ---------- */
void poly_from_signed(poly *r, const spoly *a)
{
    for (int i = 0; i < IAS_N; i++) {
        zs_t v = a->c[i];
        zs_t q = (zs_t)IAS_Q;
        v %= q; if (v < 0) v += q;
        r->c[i] = (zq_t)v;
    }
}
void poly_center(spoly *r, const poly *a)
{
    zs_t q = (zs_t)IAS_Q, h = q / 2;
    for (int i = 0; i < IAS_N; i++) {
        zs_t v = (zs_t)a->c[i];
        if (v > h) v -= q;
        r->c[i] = v;
    }
}

/* ---------- modulus rounding ---------- */
zq_t ias_q_nu(int nu) { return IAS_Q >> nu; }       

void poly_round(poly *r, const poly *a, int nu)
{
    zq_t qnu = ias_q_nu(nu);
    zq_t add = (nu > 0) ? ((zq_t)1 << (nu - 1)) : 0;
    for (int i = 0; i < IAS_N; i++) {
        zq_t v = (a->c[i] + add) >> nu;   
        if (v >= qnu) v -= qnu;         
        r->c[i] = v;
    }
}

void poly_lift_shift(poly *r, const poly *a, int nu)
{
    for (int i = 0; i < IAS_N; i++) {
        zq_t v = a->c[i] << nu;           
        r->c[i] = v % IAS_Q;
    }
}

long double poly_sqnorm_shift(const poly *a, int nu)
{
    long double s = 0;
    zs_t q = (zs_t)IAS_Q, h = q / 2;
    for (int i = 0; i < IAS_N; i++) {
        zq_t v = (a->c[i] << nu) % IAS_Q;
        zs_t sv = (zs_t)v;
        if (sv > h) sv -= q;
        long double d = (long double)sv;
        s += d * d;
    }
    return s;
}

long double spoly_sqnorm(const spoly *a)
{
    long double s = 0;
    for (int i = 0; i < IAS_N; i++) {
        long double d = (long double)a->c[i];
        s += d * d;
    }
    return s;
}

/* ---------- uniform sampling mod q ---------- */
void poly_uniform(poly *r, keccak_state *st)
{
    int nbytes = (IAS_LOGQ_CEIL + 7) / 8;      /* bytes per candidate */
    zq_t mask = ((zq_t)1 << IAS_LOGQ_CEIL) - 1;
    int i = 0;
    uint8_t buf[16];
    while (i < IAS_N) {
        shake256_squeeze(buf, nbytes, st);
        zq_t v = 0;
        for (int b = 0; b < nbytes; b++) v |= (zq_t)buf[b] << (8 * b);
        v &= mask;
        if (v < IAS_Q) r->c[i++] = v;          /* rejection sampling */
    }
}

/* ---------- 128-bit decimal printer ---------- */
int u128_to_str(char *buf, zq_t x)
{
    char tmp[40];
    int n = 0;
    if (x == 0) { buf[0] = '0'; buf[1] = 0; return 1; }
    while (x > 0) { tmp[n++] = (char)('0' + (int)(x % 10)); x /= 10; }
    for (int i = 0; i < n; i++) buf[i] = tmp[n - 1 - i];
    buf[n] = 0;
    return n;
}

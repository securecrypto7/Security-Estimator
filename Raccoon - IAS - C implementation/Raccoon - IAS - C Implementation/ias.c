#include "ias.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

/* ============================ bounds ==================================== */
void ias_compute_bounds(ias_bounds *b)
{
    const long double n = IAS_N, k = IAS_K, l = IAS_L;
    const long double ss = ldexpl(1.0L, IAS_U_S);
    const long double sr = IAS_U_R_HALF ? ldexpl(1.0L, IAS_U_R) * 0.70710678118654752440L
                                        : ldexpl(1.0L, IAS_U_R);
    const long double omega = IAS_OMEGA, rep = IAS_REP;
    const int nu_t = IAS_NU_T, nu_w = IAS_NU_W, nu_wp = IAS_NU_WP;
    const long double ALPHA_raw = 1.0L + sqrtl(4.0L * IAS_KAPPA * logl(2.0L) / n);
    const long double ALPHA = ALPHA_raw > 2.0L ? ALPHA_raw : 2.0L;
    const long double rt_full = sqrtl(n * (k + l));
    const long double rt_k    = sqrtl(n * k);
    const long double AGG_POW = 3.0L;

    /* q_perp = |q - 2^nu_w * round(q/2^nu_w)| :  */
    zq_t rem = IAS_Q & ((( zq_t)1 << nu_w) - 1);        /* q mod 2^nu_w */
    zq_t halfstep = (zq_t)1 << (nu_w - 1);
    zq_t qperp_i  = (rem <= halfstep) ? rem : (((zq_t)1 << nu_w) - rem);
    const long double q_perp = (long double)qperp_i;

    long double B2ind = 2.0L*omega*ALPHA*ss*sqrtl(n*(k+l))
                      + ALPHA*sr*sqrtl(n*(k+l)*rep)
                      + rep*ldexpl(1.0L, nu_wp-1)*rt_k;

    long double B2 = (long double)IAS_NSIGNERS * B2ind
                   + (3.0L*ldexpl(1.0L, nu_w-1) + q_perp
                      + omega*ldexpl(1.0L, nu_t-1)) * rt_k;

    long double beta = 8.0L*omega*B2
                     + (ldexpl(1.0L, nu_w+2)*omega
                        + ldexpl(1.0L, nu_t+2)*powl(omega, AGG_POW)) * rt_k
                     + 16.0L*powl(omega, AGG_POW)*ss*rt_full;

    b->B2ind = B2ind;
    b->B2    = B2;
    b->beta  = beta;
    b->z_bits = (int)ceill(log2l(2.0L*B2 + 1.0L));
    b->t_bits = IAS_LOGQ_SIZE - nu_t;
}

/* ============================ matrix A ================================= */
static void expandA(ias_pp *pp)
{
    for (int row = 0; row < IAS_K; row++)
        for (int col = 0; col < IAS_L; col++) {
            keccak_state st;
            shake256_init(&st);
            uint8_t dom = 0x10;
            shake256_absorb(&st, &dom, 1);
            shake256_absorb(&st, pp->seed, IAS_SEEDBYTES);
            absorb_u32(&st, (uint32_t)row);
            absorb_u32(&st, (uint32_t)col);
            shake256_finalize(&st);
            poly_uniform(&pp->A[row][col], &st);
            poly_copy(&pp->Antt[row][col], &pp->A[row][col]);
            poly_ntt(&pp->Antt[row][col]);
        }
}

void ias_setup(ias_pp *pp, const uint8_t seed[IAS_SEEDBYTES])
{
    memcpy(pp->seed, seed, IAS_SEEDBYTES);
    expandA(pp);
}

/* out[k] = A * v  (v short); result natural domain, mod q. */
static void matvec_short(const ias_pp *pp, const spoly v[IAS_L], poly out[IAS_K])
{
    poly vn[IAS_L];
    for (int col = 0; col < IAS_L; col++) {
        poly_from_signed(&vn[col], &v[col]);
        poly_ntt(&vn[col]);
    }
    for (int row = 0; row < IAS_K; row++) {
        poly acc; poly_zero(&acc);
        for (int col = 0; col < IAS_L; col++) {
            poly t;
            poly_pointwise(&t, &pp->Antt[row][col], &vn[col]);
            poly_add(&acc, &acc, &t);
        }
        poly_invntt(&acc);
        poly_copy(&out[row], &acc);
    }
}

/* ============================ KeyGen =================================== */
void ias_keygen(const ias_pp *pp, ias_sk *sk, ias_pk *pk)
{
    spoly e[IAS_K];
    for (int ll = 0; ll < IAS_L; ll++) sample_gauss(&sk->s[ll], ias_sigma_s());
    for (int kk = 0; kk < IAS_K; kk++) sample_gauss(&e[kk], ias_sigma_s());

    poly As[IAS_K];
    matvec_short(pp, sk->s, As);
    for (int kk = 0; kk < IAS_K; kk++) {
        poly ek, t;
        poly_from_signed(&ek, &e[kk]);
        poly_add(&t, &As[kk], &ek);              /* A s + e            */
        poly_add(&t, &t, &t);                    /* 2 (A s + e) = t    */
        poly_copy(&pk->t_full[kk], &t);          /* unrounded key      */
        poly_round(&pk->tT[kk], &t, IAS_NU_T);   /* t^T = round_{nu_t} */
    }
}

/* ============================ KeyAgg ================================== */
static void build_PKflat(const ias_pk *pks, int nsigners, poly *out)
{
    for (int i = 0; i < nsigners; i++)
        for (int kk = 0; kk < IAS_K; kk++)
            poly_copy(&out[i*IAS_K + kk], &pks[i].tT[kk]);
}

static void keyagg_core(const ias_pk *pks, int nsigners, poly *PKflat,
                        poly *a, poly ttilde[IAS_K])
{
    for (int kk = 0; kk < IAS_K; kk++) poly_zero(&ttilde[kk]);
    for (int i = 0; i < nsigners; i++) {
        H_agg(&a[i], PKflat, nsigners, i, pks[i].tT);
        for (int kk = 0; kk < IAS_K; kk++) {
            poly lift, term;
            poly_lift_shift(&lift, &pks[i].tT[kk], IAS_NU_T); /* 2^{nu_t} t^T */
            poly_mul(&term, &a[i], &lift);                    /* a_i * (..)  */
            poly_add(&ttilde[kk], &ttilde[kk], &term);
        }
    }
}

void ias_keyagg(const ias_pp *pp, const ias_pk *pks, int nsigners, ias_agg *agg)
{
    (void)pp;
    poly *PKflat = malloc(sizeof(poly) * (size_t)nsigners * IAS_K);
    build_PKflat(pks, nsigners, PKflat);
    keyagg_core(pks, nsigners, PKflat, agg->a, agg->ttilde);
    free(PKflat);
}

/* ========================= Sign round 1 ============================== */
void ias_sign_pre(const ias_pp *pp, poly *commit, spoly *st)
{
    for (int b = 0; b < IAS_REP; b++) {
        spoly r[IAS_L], ep[IAS_K];
        for (int ll = 0; ll < IAS_L; ll++) {
            sample_gauss(&r[ll], ias_sigma_r());
            st[ias_sidx(b, ll)] = r[ll];         
        }
        for (int kk = 0; kk < IAS_K; kk++) sample_gauss(&ep[kk], ias_sigma_r());

        poly Ar[IAS_K];
        matvec_short(pp, r, Ar);
        for (int kk = 0; kk < IAS_K; kk++) {
            poly ek, t;
            poly_from_signed(&ek, &ep[kk]);
            poly_add(&t, &Ar[kk], &ek);              /* A r + e'          */
            poly_round(&commit[ias_cidx(0,b,kk)], &t, IAS_NU_WP);
        }
    }
}

/* ========================= Coordinator round 1 ======================= */
void ias_coord(const ias_pp *pp, const ias_pk *pks, const poly *commits,
               int nsigners, const uint8_t *msgs, size_t mlen_each,
               ias_coord_ctx *cc)
{
    (void)pp;
    poly *PKflat = malloc(sizeof(poly) * (size_t)nsigners * IAS_K);
    build_PKflat(pks, nsigners, PKflat);

    /* a_i and t~ */
    keyagg_core(pks, nsigners, PKflat, cc->a, cc->ttilde);

    /* (beta_b) = H_beta(L, (w_{i,b})) ; beta_1 = 1 */
    H_beta(cc->beta, IAS_REP, PKflat, nsigners, msgs, mlen_each,
           commits, nsigners * IAS_REP * IAS_K);

    /* per-signer w_i = sum_b 2^{nu'_w} beta_b w_{i,b} (R_q, unrounded) */
    poly sumw[IAS_K];
    for (int kk = 0; kk < IAS_K; kk++) poly_zero(&sumw[kk]);
    for (int i = 0; i < nsigners; i++) {
        for (int kk = 0; kk < IAS_K; kk++) poly_zero(&cc->wi[i][kk]);
        for (int b = 0; b < IAS_REP; b++)
            for (int kk = 0; kk < IAS_K; kk++) {
                poly lift, term;
                poly_lift_shift(&lift, &commits[ias_cidx(i,b,kk)], IAS_NU_WP);
                poly_mul(&term, &cc->beta[b], &lift);
                poly_add(&cc->wi[i][kk], &cc->wi[i][kk], &term);
            }
        for (int kk = 0; kk < IAS_K; kk++)
            poly_add(&sumw[kk], &sumw[kk], &cc->wi[i][kk]);
    }
    /* aggregate rounded commitment w = round_{nu_w}(sum_j w_j) */
    for (int kk = 0; kk < IAS_K; kk++)
        poly_round(&cc->w[kk], &sumw[kk], IAS_NU_W);

    free(PKflat);
}

/* ========================= Sign round 2 ============================== */
void ias_sign_resp(const ias_pp *pp, const ias_sk *sk, int idx,
                   const ias_agg *agg, const ias_coord_ctx *cc,
                   const spoly *st, const uint8_t *msgs, size_t mlen_each,
                   ias_resp *resp)
{
    (void)pp;
    /* c = H_c(t~, (m_i), w) */
    H_c(&resp->c, cc->ttilde, IAS_K, msgs, mlen_each, IAS_NSIGNERS,
        cc->w, IAS_K);

    /* Precompute NTT(a_idx * c) and NTT(beta_b). */
    poly acc_poly, ac_ntt;
    poly_mul(&acc_poly, &agg->a[idx], &resp->c);   /* a_i * c (natural) */
    poly_copy(&ac_ntt, &acc_poly);
    poly_ntt(&ac_ntt);

    poly beta_ntt[IAS_REP];
    for (int b = 0; b < IAS_REP; b++) {
        poly_copy(&beta_ntt[b], &cc->beta[b]);
        poly_ntt(&beta_ntt[b]);
    }

    for (int ll = 0; ll < IAS_L; ll++) {
        /* term1 = 2 * (a_i c) * s_ll */
        poly s_ntt, z_ntt, t;
        poly_from_signed(&s_ntt, &sk->s[ll]);
        poly_ntt(&s_ntt);
        poly_pointwise(&z_ntt, &ac_ntt, &s_ntt);
        poly_add(&z_ntt, &z_ntt, &z_ntt);          /* times 2           */
        /* term2 = sum_b beta_b * r_{i,b,ll} */
        for (int b = 0; b < IAS_REP; b++) {
            poly r_ntt;
            poly_from_signed(&r_ntt, &st[ias_sidx(b, ll)]);
            poly_ntt(&r_ntt);
            poly_pointwise(&t, &beta_ntt[b], &r_ntt);
            poly_add(&z_ntt, &z_ntt, &t);
        }
        poly_invntt(&z_ntt);
        poly_center(&resp->z[ll], &z_ntt);         /* short signed z_ll  */
    }
}

/* ===================== Coordinator combine =========================== */
int ias_coord_combine(const ias_pp *pp, const ias_pk *pks,
                      const ias_coord_ctx *cc, const ias_resp *resps,
                      int nsigners, const ias_bounds *bnd, ias_sig *sig)
{
    /* all c_i equal? */
    for (int i = 1; i < nsigners; i++)
        for (int j = 0; j < IAS_N; j++)
            if (resps[i].c.c[j] != resps[0].c.c[j]) return -1;
    poly_copy(&sig->c, &resps[0].c);

    const long double B2ind2 = bnd->B2ind * bnd->B2ind;

    /* per-signer check: || (z_i, y_i) ||_2 <= B2ind,
     * y_i = w_i - (A z_i - c a_i 2^{nu_t} t_i^T)                          */
    for (int i = 0; i < nsigners; i++) {
        poly Azi[IAS_K];
        matvec_short(pp, resps[i].z, Azi);

        poly cai;
        poly_mul(&cai, &sig->c, &cc->a[i]);        /* c * a_i            */

        long double sq = spoly_sqnorm(&resps[i].z[0]);
        for (int ll = 1; ll < IAS_L; ll++) sq += spoly_sqnorm(&resps[i].z[ll]);

        for (int kk = 0; kk < IAS_K; kk++) {
            poly term, tmp, yy;
            spoly yc;
            poly_mul(&term, &cai, &pks[i].t_full[kk]);         /* c a_i t_i     */
            poly_sub(&tmp, &Azi[kk], &term);                   /* A z_i - term  */
            poly_sub(&yy, &cc->wi[i][kk], &tmp);               /* y_i            */
            poly_center(&yc, &yy);
            sq += spoly_sqnorm(&yc);
        }
        if (sq > B2ind2) return -1;                /* abort: signer i    */
    }

    for (int ll = 0; ll < IAS_L; ll++) {
        for (int j = 0; j < IAS_N; j++) {
            zs_t acc = 0;
            for (int i = 0; i < nsigners; i++) acc += resps[i].z[ll].c[j];
            sig->z[ll].c[j] = acc;
        }
    }

    poly Az[IAS_K];
    matvec_short(pp, sig->z, Az);
    zq_t qnw = ias_q_nu(IAS_NU_W);
    for (int kk = 0; kk < IAS_K; kk++) {
        poly ct, tmp, u;
        poly_mul(&ct, &sig->c, &cc->ttilde[kk]);
        poly_sub(&tmp, &Az[kk], &ct);
        poly_round(&u, &tmp, IAS_NU_W);
        for (int j = 0; j < IAS_N; j++)
            sig->h[kk].c[j] = (cc->w[kk].c[j] + qnw - u.c[j]) % qnw;  /* w - u */
    }
    return 0;
}

/* ============================ Verify ================================== */
int ias_verify(const ias_pp *pp, const ias_agg *agg,
               const uint8_t *msgs, size_t mlen_each, int nsigners,
               const ias_bounds *bnd, const ias_sig *sig)
{
    (void)nsigners;
    /* w' = round_{nu_w}(A z - c t~) + h  mod q_{nu_w} */
    poly Az[IAS_K], wp[IAS_K];
    matvec_short(pp, sig->z, Az);
    zq_t qnw = ias_q_nu(IAS_NU_W);
    for (int kk = 0; kk < IAS_K; kk++) {
        poly ct, tmp, u;
        poly_mul(&ct, &sig->c, &agg->ttilde[kk]);
        poly_sub(&tmp, &Az[kk], &ct);
        poly_round(&u, &tmp, IAS_NU_W);
        for (int j = 0; j < IAS_N; j++)
            wp[kk].c[j] = (u.c[j] + sig->h[kk].c[j]) % qnw;
    }

    /* c' = H_c(t~, (m_i), w') */
    poly cprime;
    H_c(&cprime, agg->ttilde, IAS_K, msgs, mlen_each, IAS_NSIGNERS, wp, IAS_K);
    for (int j = 0; j < IAS_N; j++)
        if (cprime.c[j] != sig->c.c[j]) return 0;

    /* norm: || (z, 2^{nu_w} h mod q) ||_2 <= B2 */
    long double sq = 0;
    for (int ll = 0; ll < IAS_L; ll++) sq += spoly_sqnorm(&sig->z[ll]);
    for (int kk = 0; kk < IAS_K; kk++) sq += poly_sqnorm_shift(&sig->h[kk], IAS_NU_W);
    if (sq > bnd->B2 * bnd->B2) return 0;

    return 1;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "ias.h"
#include "serialize.h"

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static void seedbuf(uint8_t s[IAS_SEEDBYTES], uint32_t x)
{ for (int i = 0; i < IAS_SEEDBYTES; i++) s[i] = (uint8_t)(x*2654435761u + i*40503u); }

/* adaptive timing: run body until >= MIN_ITERS and >= MIN_SECS elapsed */
#define MIN_ITERS 5
#define MAX_ITERS 2000
#define MIN_SECS  0.40

/* globals reused across timed loops  */
static ias_pp   *pp;
static ias_sk   *sks;
static ias_pk   *pks;
static ias_agg  *agg;
static ias_coord_ctx *cc;
static ias_resp *resps;
static ias_sig  *sig;
static poly     *commits;
static spoly    *states;
static uint8_t  *msgs;
static ias_bounds bnd;
static int NS;
static size_t MLEN = 32;

int main(void)
{
    ias_poly_init();
    ias_compute_bounds(&bnd);
    NS = IAS_NSIGNERS;

    char qb[48]; u128_to_str(qb, IAS_Q);
    printf("========================================================\n");
    printf(" %s\n", IAS_NAME);
    printf("========================================================\n");
    printf("q = %s   (ceil(log2 q) = %d)\n", qb, IAS_LOGQ_CEIL);
    printf("n=%d  k=%d  l=%d  nsigners=%d  omega=%d  rep=%d\n",
           IAS_N, IAS_K, IAS_L, NS, IAS_OMEGA, IAS_REP);
    printf("nu_t=%d  nu_w=%d  nu_w'=%d\n", IAS_NU_T, IAS_NU_W, IAS_NU_WP);
    printf("B2ind=2^%.3Lf  B2=2^%.3Lf  beta(SIS)=2^%.3Lf\n",
           log2l(bnd.B2ind), log2l(bnd.B2), log2l(bnd.beta));
    printf("z_bits=%d  t_bits=%d\n\n", bnd.z_bits, bnd.t_bits);

    pp      = malloc(sizeof *pp);
    sks     = malloc(sizeof(ias_sk) * NS);
    pks     = malloc(sizeof(ias_pk) * NS);
    agg     = malloc(sizeof *agg);
    cc      = malloc(sizeof *cc);
    resps   = malloc(sizeof(ias_resp) * NS);
    sig     = malloc(sizeof *sig);
    commits = malloc(sizeof(poly)  * (size_t)NS * IAS_COMMIT_POLYS);
    states  = malloc(sizeof(spoly) * (size_t)NS * IAS_STATE_POLYS);
    msgs    = malloc(MLEN * NS);
    if (!pp||!sks||!pks||!agg||!cc||!resps||!sig||!commits||!states||!msgs) {
        printf("alloc failed\n"); return 2;
    }
    for (size_t i = 0; i < MLEN*(size_t)NS; i++) msgs[i] = (uint8_t)(i*131 + 7);

    uint8_t seed[IAS_SEEDBYTES]; seedbuf(seed, 0xA5A5A5A5u);

   
    ias_setup(pp, seed);
    for (int i = 0; i < NS; i++) {
        uint8_t sd[IAS_SEEDBYTES]; seedbuf(sd, 1000u + i);
        ias_rng_seed(sd); ias_keygen(pp, &sks[i], &pks[i]);
    }
    ias_keyagg(pp, pks, NS, agg);
    for (int i = 0; i < NS; i++) {
        uint8_t sd[IAS_SEEDBYTES]; seedbuf(sd, 7000u + i);
        ias_rng_seed(sd);
        ias_sign_pre(pp, &commits[(size_t)i*IAS_COMMIT_POLYS],
                         &states[(size_t)i*IAS_STATE_POLYS]);
    }
    ias_coord(pp, pks, commits, NS, msgs, MLEN, cc);
    for (int i = 0; i < NS; i++)
        ias_sign_resp(pp, &sks[i], i, agg, cc,
                      &states[(size_t)i*IAS_STATE_POLYS], msgs, MLEN, &resps[i]);
    int rc = ias_coord_combine(pp, pks, cc, resps, NS, &bnd, sig);
    int ok = ias_verify(pp, agg, msgs, MLEN, NS, &bnd, sig);
    if (rc != 0 || !ok) { printf("FATAL: honest run failed (rc=%d ok=%d)\n", rc, ok); return 3; }

    /* ================= timing ================= */
    printf("----------- timings  -------------------\n");
    double t0, el; long it;

    /* Setup */
    t0 = now_s(); it = 0;
    do { ias_setup(pp, seed); it++; el = now_s()-t0; }
    while (it < MIN_ITERS || (el < MIN_SECS && it < MAX_ITERS));
    double t_setup = el/it*1e3;

    /* KeyGen (one signer) */
    { uint8_t sd[IAS_SEEDBYTES]; seedbuf(sd, 1000u);
      t0 = now_s(); it = 0;
      do { ias_rng_seed(sd); ias_keygen(pp, &sks[0], &pks[0]); it++; el = now_s()-t0; }
      while (it < MIN_ITERS || (el < MIN_SECS && it < MAX_ITERS));
    }
    double t_kg = el/it*1e3;
    for (int i = 0; i < NS; i++) {
        uint8_t sd[IAS_SEEDBYTES]; seedbuf(sd, 1000u + i);
        ias_rng_seed(sd); ias_keygen(pp, &sks[i], &pks[i]);
    }
    ias_keyagg(pp, pks, NS, agg);

    /* KeyAgg (all signers) */
    t0 = now_s(); it = 0;
    do { ias_keyagg(pp, pks, NS, agg); it++; el = now_s()-t0; }
    while (it < MIN_ITERS || (el < MIN_SECS && it < MAX_ITERS));
    double t_ka = el/it*1e3;

    /* Sign_pre (one signer) */
    t0 = now_s(); it = 0;
    do { uint8_t sd[IAS_SEEDBYTES]; seedbuf(sd, 7000u); ias_rng_seed(sd);
         ias_sign_pre(pp, &commits[0], &states[0]); it++; el = now_s()-t0; }
    while (it < MIN_ITERS || (el < MIN_SECS && it < MAX_ITERS));
    double t_sp = el/it*1e3;
    /* restore all commitments/states */
    for (int i = 0; i < NS; i++) {
        uint8_t sd[IAS_SEEDBYTES]; seedbuf(sd, 7000u + i);
        ias_rng_seed(sd);
        ias_sign_pre(pp, &commits[(size_t)i*IAS_COMMIT_POLYS],
                         &states[(size_t)i*IAS_STATE_POLYS]);
    }

    /* Coord (all signers) */
    t0 = now_s(); it = 0;
    do { ias_coord(pp, pks, commits, NS, msgs, MLEN, cc); it++; el = now_s()-t0; }
    while (it < MIN_ITERS || (el < MIN_SECS && it < MAX_ITERS));
    double t_co = el/it*1e3;

    /* Sign_resp (one signer) */
    t0 = now_s(); it = 0;
    do { ias_sign_resp(pp, &sks[0], 0, agg, cc, &states[0], msgs, MLEN, &resps[0]);
         it++; el = now_s()-t0; }
    while (it < MIN_ITERS || (el < MIN_SECS && it < MAX_ITERS));
    double t_sr = el/it*1e3;
    for (int i = 0; i < NS; i++)
        ias_sign_resp(pp, &sks[i], i, agg, cc,
                      &states[(size_t)i*IAS_STATE_POLYS], msgs, MLEN, &resps[i]);

    /* Coord_combine (all signers) */
    t0 = now_s(); it = 0;
    do { rc = ias_coord_combine(pp, pks, cc, resps, NS, &bnd, sig); it++; el = now_s()-t0; }
    while (it < MIN_ITERS || (el < MIN_SECS && it < MAX_ITERS));
    double t_cc = el/it*1e3;

    /* Verify */
    t0 = now_s(); it = 0;
    do { ok = ias_verify(pp, agg, msgs, MLEN, NS, &bnd, sig); it++; el = now_s()-t0; }
    while (it < MIN_ITERS || (el < MIN_SECS && it < MAX_ITERS));
    double t_vf = el/it*1e3;

    printf("Setup            : %9.3f ms\n", t_setup);
    printf("KeyGen  /signer  : %9.3f ms   (x%d = %.3f ms)\n", t_kg, NS, t_kg*NS);
    printf("KeyAgg  (all)    : %9.3f ms\n", t_ka);
    printf("Sign_pre/signer  : %9.3f ms   (x%d = %.3f ms)\n", t_sp, NS, t_sp*NS);
    printf("Coord   (all)    : %9.3f ms\n", t_co);
    printf("Sign_resp/signer : %9.3f ms   (x%d = %.3f ms)\n", t_sr, NS, t_sr*NS);
    printf("Coord_combine    : %9.3f ms\n", t_cc);
    printf("Verify           : %9.3f ms\n", t_vf);
    double signer_total = t_kg + t_sp + t_sr;
    double coord_total  = t_ka + t_co + t_cc;
    printf("  per-signer total (KeyGen+Sign_pre+Sign_resp): %.3f ms\n", signer_total);
    printf("  coordinator total (KeyAgg+Coord+Combine)    : %.3f ms\n", coord_total);

    /* ================= serialization + sizes ================= */
    size_t vkb = ias_vk_bytes(&bnd);
    size_t zb  = ias_sig_z_bytes(&bnd);
    size_t cbz = ias_sig_c_bytes();
    size_t hb  = ias_sig_h_bytes();
    size_t fb  = ias_sig_full_bytes(&bnd);

    uint8_t *vk = malloc(vkb);
    uint8_t *zp = malloc(zb);
    uint8_t *cp = malloc(cbz);
    uint8_t *hp = malloc(hb);
    ias_pack_vk(vk, pp, &pks[0], &bnd);
    ias_pack_z (zp, sig, &bnd);
    ias_pack_c (cp, &sig->c);
    ias_pack_h (hp, sig);

    /* unpack z and check exact equality with the signed source */
    spoly *z2 = malloc(sizeof(spoly)*IAS_L);
    ias_unpack_z(z2, zp, &bnd);
    int zmatch = 1;
    for (int ll = 0; ll < IAS_L && zmatch; ll++)
        for (int j = 0; j < IAS_N; j++)
            if (z2[ll].c[j] != sig->z[ll].c[j]) { zmatch = 0; break; }

    printf("\n--- sizes (bytes) ---------------------------------------\n");
    printf("VK  = seed(%d) + t^T packed = %zu\n", IAS_SEEDBYTES, vkb);
    printf("SIG (z only, estimator sig) = %zu\n", zb);
    printf("  + c = %zu, + h = %zu  => full on-wire sig = %zu\n", cbz, hb, fb);

    /* estimator closed forms */
    size_t est_vk  = (size_t)IAS_SEEDBYTES + (size_t)IAS_K*IAS_N*(size_t)bnd.t_bits/8;
    size_t est_sig = (size_t)IAS_L*IAS_N*(size_t)bnd.z_bits/8;
    int vk_ok  = (vkb == est_vk);
    int sig_ok = (zb  == est_sig);

    int pass = (rc==0) && ok && zmatch && vk_ok && sig_ok;
    printf("\nRESULT: %s\n", pass ? "PASS" : "FAIL");

    free(vk); free(zp); free(cp); free(hp); free(z2);
    free(pp); free(sks); free(pks); free(agg); free(cc);
    free(resps); free(sig); free(commits); free(states); free(msgs);
    return pass ? 0 : 1;
}

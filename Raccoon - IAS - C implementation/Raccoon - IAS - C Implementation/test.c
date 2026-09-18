#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ias.h"

static void seedbuf(uint8_t s[IAS_SEEDBYTES], uint32_t x)
{ for (int i = 0; i < IAS_SEEDBYTES; i++) s[i] = (uint8_t)(x*2654435761u + i*40503u); }

int main(void)
{
    ias_poly_init();
    ias_bounds bnd; ias_compute_bounds(&bnd);

    char qb[48]; u128_to_str(qb, IAS_Q);
    printf("=== %s ===\n", IAS_NAME);
    printf("q=%s  ceil(log2 q)=%d  n=%d  (k,l)=(%d,%d)  nsigners=%d\n",
           qb, IAS_LOGQ_CEIL, IAS_N, IAS_K, IAS_L, IAS_NSIGNERS);
    printf("B2ind=2^%.3Lf  B2=2^%.3Lf  beta=2^%.3Lf  z_bits=%d  t_bits=%d\n",
           log2l(bnd.B2ind), log2l(bnd.B2), log2l(bnd.beta), bnd.z_bits, bnd.t_bits);

    const int NS = IAS_NSIGNERS;
    const size_t MLEN = 32;

    ias_pp   *pp   = malloc(sizeof *pp);
    ias_sk   *sks  = malloc(sizeof(ias_sk)  * NS);
    ias_pk   *pks  = malloc(sizeof(ias_pk)  * NS);
    ias_agg  *agg  = malloc(sizeof *agg);
    ias_coord_ctx *cc = malloc(sizeof *cc);
    ias_resp *resps = malloc(sizeof(ias_resp) * NS);
    ias_sig  *sig  = malloc(sizeof *sig);
    poly  *commits = malloc(sizeof(poly)  * (size_t)NS * IAS_COMMIT_POLYS);
    spoly *states  = malloc(sizeof(spoly) * (size_t)NS * IAS_STATE_POLYS);
    uint8_t *msgs  = malloc(MLEN * NS);
    if (!pp||!sks||!pks||!agg||!cc||!resps||!sig||!commits||!states||!msgs) {
        printf("alloc failed\n"); return 2;
    }
    for (size_t i = 0; i < MLEN*NS; i++) msgs[i] = (uint8_t)(i*131 + 7);

    /* Setup */
    uint8_t seed[IAS_SEEDBYTES]; seedbuf(seed, 0xA5A5A5A5u);
    ias_setup(pp, seed);

    /* KeyGen */
    for (int i = 0; i < NS; i++) {
        uint8_t sd[IAS_SEEDBYTES]; seedbuf(sd, 1000u + i);
        ias_rng_seed(sd);
        ias_keygen(pp, &sks[i], &pks[i]);
    }

    /* KeyAgg */
    ias_keyagg(pp, pks, NS, agg);

    /* Sign round 1 (per signer) */
    for (int i = 0; i < NS; i++) {
        uint8_t sd[IAS_SEEDBYTES]; seedbuf(sd, 7000u + i);
        ias_rng_seed(sd);
        ias_sign_pre(pp, &commits[(size_t)i*IAS_COMMIT_POLYS],
                         &states[(size_t)i*IAS_STATE_POLYS]);
    }

    /* Coordinator round 1 */
    ias_coord(pp, pks, commits, NS, msgs, MLEN, cc);

    /* Sign round 2 (per signer) */
    for (int i = 0; i < NS; i++)
        ias_sign_resp(pp, &sks[i], i, agg, cc,
                      &states[(size_t)i*IAS_STATE_POLYS], msgs, MLEN, &resps[i]);

    /* Coordinator combine */
    int rc = ias_coord_combine(pp, pks, cc, resps, NS, &bnd, sig);
    printf("coord_combine: %s\n", rc == 0 ? "OK" : "ABORT (norm/consistency)");

    /* report actual aggregate norm vs B2 */
    long double sq = 0;
    for (int ll = 0; ll < IAS_L; ll++) sq += spoly_sqnorm(&sig->z[ll]);
    long double zsq = sq;
    for (int kk = 0; kk < IAS_K; kk++) sq += poly_sqnorm_shift(&sig->h[kk], IAS_NU_W);
    printf("actual ||z||=2^%.3Lf   ||(z,2^nw h)||=2^%.3Lf   (B2=2^%.3Lf)\n",
           0.5L*log2l(zsq), 0.5L*log2l(sq), log2l(bnd.B2));

    /* Verify */
    int ok = ias_verify(pp, agg, msgs, MLEN, NS, &bnd, sig);
    printf("verify: %s\n", ok ? "ACCEPT" : "REJECT");

    /* Negative test: flip one z coefficient -> must reject */
    ias_sig *bad = malloc(sizeof *bad); memcpy(bad, sig, sizeof *bad);
    bad->z[0].c[0] += 1;
    int ok2 = ias_verify(pp, agg, msgs, MLEN, NS, &bnd, bad);
    printf("verify(tampered z): %s\n", ok2 ? "ACCEPT (BUG!)" : "reject (good)");

    /* Negative test: wrong message -> must reject */
    uint8_t *msgs2 = malloc(MLEN*NS); memcpy(msgs2, msgs, MLEN*NS);
    msgs2[0] ^= 0xff;
    int ok3 = ias_verify(pp, agg, msgs2, MLEN, NS, &bnd, sig);
    printf("verify(wrong msg): %s\n", ok3 ? "ACCEPT (BUG!)" : "reject (good)");

    int pass = (rc==0) && ok && !ok2 && !ok3;
    printf("\nRESULT: %s\n", pass ? "PASS" : "FAIL");

    free(bad); free(msgs2);
    free(pp); free(sks); free(pks); free(agg); free(cc);
    free(resps); free(sig); free(commits); free(states); free(msgs);
    return pass ? 0 : 1;
}

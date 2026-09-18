#include <string.h>
#include <assert.h>
#include "serialize.h"

/* ---- little-endian bit writer / reader over a byte buffer ---- */
static void bitpack(uint8_t *out, size_t *bitpos, uint64_t val, int nbits)
{
    for (int i = 0; i < nbits; i++) {
        if (val & (1ULL << i)) {
            size_t bp = *bitpos + (size_t)i;
            out[bp >> 3] |= (uint8_t)(1u << (bp & 7));
        }
    }
    *bitpos += (size_t)nbits;
}
static uint64_t bitunpack(const uint8_t *in, size_t *bitpos, int nbits)
{
    uint64_t val = 0;
    for (int i = 0; i < nbits; i++) {
        size_t bp = *bitpos + (size_t)i;
        if (in[bp >> 3] & (uint8_t)(1u << (bp & 7)))
            val |= (1ULL << i);
    }
    *bitpos += (size_t)nbits;
    return val;
}

/* ---- sizes ---- */
size_t ias_vk_bytes(const ias_bounds *b)
{
    return (size_t)IAS_SEEDBYTES
         + (size_t)IAS_K * ((size_t)IAS_N * (size_t)b->t_bits / 8);
}
size_t ias_sig_z_bytes(const ias_bounds *b)
{
    return (size_t)IAS_L * ((size_t)IAS_N * (size_t)b->z_bits / 8);
}
size_t ias_sig_c_bytes(void)
{
    /* omega nonzero coeffs: 9-bit position (ceil(log2 512)) + 1 sign bit  */
    size_t bits = (size_t)IAS_OMEGA * (9 + 1);
    return (bits + 7) / 8;
}
size_t ias_sig_h_bytes(void)
{
    int h_bits = IAS_LOGQ_SIZE - IAS_NU_W;
    return (size_t)IAS_K * ((size_t)IAS_N * (size_t)h_bits / 8);
}
size_t ias_sig_full_bytes(const ias_bounds *b)
{
    return ias_sig_z_bytes(b) + ias_sig_c_bytes() + ias_sig_h_bytes();
}

/* ---- VK = seed || packed t^T (rounded key, t_bits per coeff) ---- */
void ias_pack_vk(uint8_t *out, const ias_pp *pp, const ias_pk *pk,
                 const ias_bounds *b)
{
    size_t vkb = ias_vk_bytes(b);
    memset(out, 0, vkb);
    memcpy(out, pp->seed, IAS_SEEDBYTES);
    size_t bitpos = (size_t)IAS_SEEDBYTES * 8;
    for (int kk = 0; kk < IAS_K; kk++)
        for (int j = 0; j < IAS_N; j++) {
            uint64_t v = (uint64_t)pk->tT[kk].c[j];   /* in [0, q_{nu_t})  */
            assert(b->t_bits >= 64 || v < (1ULL << b->t_bits));
            bitpack(out, &bitpos, v, b->t_bits);
        }
}

/* ---- z vector : signed coeffs, biased into [0,2^z_bits) then packed ---- */
void ias_pack_z(uint8_t *out, const ias_sig *sig, const ias_bounds *b)
{
    size_t zb = ias_sig_z_bytes(b);
    memset(out, 0, zb);
    size_t bitpos = 0;
    uint64_t bias = 1ULL << (b->z_bits - 1);
    for (int ll = 0; ll < IAS_L; ll++)
        for (int j = 0; j < IAS_N; j++) {
            zs_t zc = sig->z[ll].c[j];
            /* every |coeff| is far below the L2 bound B2, so it fits */
            assert(zc > -(zs_t)bias && zc < (zs_t)bias);
            uint64_t v = (uint64_t)((zs_t)zc + (zs_t)bias);
            bitpack(out, &bitpos, v, b->z_bits);
        }
}
int ias_unpack_z(spoly z[IAS_L], const uint8_t *in, const ias_bounds *b)
{
    size_t bitpos = 0;
    zs_t bias = (zs_t)1 << (b->z_bits - 1);
    for (int ll = 0; ll < IAS_L; ll++)
        for (int j = 0; j < IAS_N; j++) {
            uint64_t v = bitunpack(in, &bitpos, b->z_bits);
            z[ll].c[j] = (zs_t)v - bias;
        }
    return 0;
}

/* ---- challenge c : omega signed unit coeffs (position + sign) ---- */
void ias_pack_c(uint8_t *out, const poly *c)
{
    size_t cb = ias_sig_c_bytes();
    memset(out, 0, cb);
    size_t bitpos = 0;
    zq_t qm1 = IAS_Q - 1;
    int written = 0;
    for (int j = 0; j < IAS_N && written < IAS_OMEGA; j++) {
        if (c->c[j] == 0) continue;
        int sign = (c->c[j] == qm1) ? 1 : 0;   /* -1 == q-1 -> sign bit    */
        bitpack(out, &bitpos, (uint64_t)j, 9);
        bitpack(out, &bitpos, (uint64_t)sign, 1);
        written++;
    }
}

/* ---- hint h : k polys, rounded to [0,q_{nu_w}), h_bits per coeff ---- */
void ias_pack_h(uint8_t *out, const ias_sig *sig)
{
    int h_bits = IAS_LOGQ_SIZE - IAS_NU_W;
    size_t hb = ias_sig_h_bytes();
    memset(out, 0, hb);
    size_t bitpos = 0;
    for (int kk = 0; kk < IAS_K; kk++)
        for (int j = 0; j < IAS_N; j++) {
            uint64_t v = (uint64_t)sig->h[kk].c[j];
            assert(h_bits >= 64 || v < (1ULL << h_bits));
            bitpack(out, &bitpos, v, h_bits);
        }
}

#include "fips202.h"
#include <string.h>

#define ROL(x, s) (((x) << (s)) | ((x) >> (64 - (s))))

static const uint64_t KECCAK_RC[24] = {
    0x0000000000000001ULL, 0x0000000000008082ULL, 0x800000000000808aULL,
    0x8000000080008000ULL, 0x000000000000808bULL, 0x0000000080000001ULL,
    0x8000000080008081ULL, 0x8000000000008009ULL, 0x000000000000008aULL,
    0x0000000000000088ULL, 0x0000000080008009ULL, 0x000000008000000aULL,
    0x000000008000808bULL, 0x800000000000008bULL, 0x8000000000008089ULL,
    0x8000000000008003ULL, 0x8000000000008002ULL, 0x8000000000000080ULL,
    0x000000000000800aULL, 0x800000008000000aULL, 0x8000000080008081ULL,
    0x8000000000008080ULL, 0x0000000080000001ULL, 0x8000000080008008ULL
};

static void keccak_f1600(uint64_t st[25])
{
    int round;
    uint64_t bc[5], t;
    for (round = 0; round < 24; round++) {
        /* Theta */
        for (int i = 0; i < 5; i++)
            bc[i] = st[i] ^ st[i + 5] ^ st[i + 10] ^ st[i + 15] ^ st[i + 20];
        for (int i = 0; i < 5; i++) {
            t = bc[(i + 4) % 5] ^ ROL(bc[(i + 1) % 5], 1);
            for (int j = 0; j < 25; j += 5)
                st[j + i] ^= t;
        }
        /* Rho + Pi */
        t = st[1];
        {
            static const int piln[24] = {10,7,11,17,18,3,5,16,8,21,24,4,
                                         15,23,19,13,12,2,20,14,22,9,6,1};
            static const int rotc[24] = {1,3,6,10,15,21,28,36,45,55,2,14,
                                         27,41,56,8,25,43,62,18,39,61,20,44};
            for (int i = 0; i < 24; i++) {
                int j = piln[i];
                bc[0] = st[j];
                st[j] = ROL(t, rotc[i]);
                t = bc[0];
            }
        }
        /* Chi */
        for (int j = 0; j < 25; j += 5) {
            for (int i = 0; i < 5; i++)
                bc[i] = st[j + i];
            for (int i = 0; i < 5; i++)
                st[j + i] ^= (~bc[(i + 1) % 5]) & bc[(i + 2) % 5];
        }
        /* Iota */
        st[0] ^= KECCAK_RC[round];
    }
}

void shake256_init(keccak_state *st)
{
    memset(st->s, 0, sizeof(st->s));
    st->pos = 0;
}

void shake256_absorb(keccak_state *st, const uint8_t *in, size_t inlen)
{
    unsigned int pos = st->pos;
    while (inlen > 0) {
        size_t take = SHAKE256_RATE - pos;
        if (take > inlen) take = inlen;
        for (size_t i = 0; i < take; i++)
            ((uint8_t *)st->s)[pos + i] ^= in[i];
        pos += (unsigned int)take;
        in += take;
        inlen -= take;
        if (pos == SHAKE256_RATE) {
            keccak_f1600(st->s);
            pos = 0;
        }
    }
    st->pos = pos;
}

void shake256_finalize(keccak_state *st)
{
    /* SHAKE domain separation 0x1F, pad10*1 with final 0x80 */
    ((uint8_t *)st->s)[st->pos] ^= 0x1F;
    ((uint8_t *)st->s)[SHAKE256_RATE - 1] ^= 0x80;
    keccak_f1600(st->s);
    st->pos = 0;
}

void shake256_squeeze(uint8_t *out, size_t outlen, keccak_state *st)
{
    unsigned int pos = st->pos;
    while (outlen > 0) {
        if (pos == SHAKE256_RATE) {
            keccak_f1600(st->s);
            pos = 0;
        }
        size_t give = SHAKE256_RATE - pos;
        if (give > outlen) give = outlen;
        memcpy(out, ((uint8_t *)st->s) + pos, give);
        out += give;
        outlen -= give;
        pos += (unsigned int)give;
    }
    st->pos = pos;
}

void shake256(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen)
{
    keccak_state st;
    shake256_init(&st);
    shake256_absorb(&st, in, inlen);
    shake256_finalize(&st);
    shake256_squeeze(out, outlen, &st);
}

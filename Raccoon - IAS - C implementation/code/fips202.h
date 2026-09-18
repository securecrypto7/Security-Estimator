#ifndef IAS_FIPS202_H
#define IAS_FIPS202_H

#include <stddef.h>
#include <stdint.h>

#define SHAKE256_RATE 136

typedef struct {
    uint64_t s[25];
    unsigned int pos;  
} keccak_state;

/* Incremental SHAKE256 (XOF) */
void shake256_init(keccak_state *st);
void shake256_absorb(keccak_state *st, const uint8_t *in, size_t inlen);
void shake256_finalize(keccak_state *st);
void shake256_squeeze(uint8_t *out, size_t outlen, keccak_state *st);

/* One-shot SHAKE256 */
void shake256(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen);

#endif
